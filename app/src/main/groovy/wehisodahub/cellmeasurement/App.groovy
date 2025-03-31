package wehisodahub.cellmeasurement

import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Option

import java.awt.geom.Point2D

import ij.IJ
import ij.process.ColorProcessor

import qupath.imagej.processing.RoiLabeling
import qupath.imagej.tools.IJTools

import qupath.lib.roi.interfaces.ROI
import qupath.lib.roi.ROIs
import qupath.lib.scripting.QP
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane
import qupath.lib.io.PathIO.GeoJsonExportOptions
import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.images.servers.PixelCalibration
import qupath.lib.images.servers.bioformats.BioFormatsServerBuilder

@Command(name = 'cellmeasurement',
         mixinStandardHelpOptions = true,
         version = '0.1',
         description = 'Extract cell measurements from nuclear and whole-cell segmentation masks.')
class App implements Runnable {

    @Option(names = ['-n', '--nuclear-mask'],
            description = 'Nuclear segmentation mask file in TIFF format',
            required = true)
    String nuclearMaskFilePath

    @Option(names = ['-w', '--whole-cell-mask'],
            description = 'Whole-cell segmentation mask file in TIFF format',
            required = true)
    String wholeCellMaskFilePath

    @Option(names = ['-t', '--tiff-file'],
            description = 'TIFF file containing multi-channel image data',
            required = true)
    String tiffFilePath

    @Option(names = ['-o', '--output-file'],
            description = 'Output path for GeoJSON file',
            required = true)
    String outputFilePath

    @Option(names = ['-d', '--downsample-factor'],
            description = 'Downsample factor',
            required = false)
    double downsampleFactor = 1.0

    @Option(names = ['-p', '--pixel-size-microns'],
            description = 'Pixel size in microns (default: 0.5)',
            required = false)
    double pixelSizeMicrons = 0.5

    @Option(names = ['--skip-measurements'],
            description = 'Skip adding measurements',
            required = false)
    boolean skipMeasurements = false


    @Option(names = ['-i', '--dist-threshold'],
            description = 'Distance threshold (in pixels) for matching ROIs',
            required = false)
    double distThreshold = 10.0

    // Extract ROIs from a given image
    static List<ROI> extractROIs(image, downsampleFactor) {
        def ip = image.getProcessor()
        if (ip instanceof ColorProcessor) {
            throw new IllegalArgumentException('RGB images are not supported!')
        }

        int n = image.getStatistics().max as int
        if (n == 0) {
            println 'No objects found in mask!'
            return []
        }

        // Convert mask to ROIs
        def roisIJ = RoiLabeling.labelsToConnectedROIs(ip, n)
        println 'Number of ROIs found: ' + roisIJ.size()

        return roisIJ.collect {
            if (it == null) return
            return IJTools.convertToROI(it, 0, 0, downsampleFactor, ImagePlane.getDefaultPlane())
        }.findAll { it != null }
    }

    static List<PathObject> createCellObjects(List<ROI> wholeCellROIs, List<ROI> nuclearROIs, double distThreshold) {
        def matchedPairs = matchROIs(nuclearROIs, wholeCellROIs, distThreshold)
        return matchedPairs.collect { nucleus, cell ->
            if (cell != null) {
                return PathObjects.createCellObject(cell, nucleus)
            }
        }.findAll { it != null }
    }

    static List<List<ROI>> matchROIs(List<ROI> nuclearROIs, List<ROI> wholeCellROIs, double distThreshold) {
        def matchedPairs = []

        nuclearROIs.each { nuclearROI ->
            def nuclearCentroid = new Point2D.Double(nuclearROI.getCentroidX(), nuclearROI.getCentroidY())
            def nearestCell = findNearestROI(nuclearCentroid, wholeCellROIs, distThreshold)
            matchedPairs << [nuclearROI, nearestCell]
        }

        return matchedPairs
    }

    static ROI findNearestROI(Point2D centroid, List<ROI> rois, double distThreshold) {
        ROI nearestROI = null
        double minDistance = Double.MAX_VALUE

        rois.each { roi ->
            def roiCentroid = new Point2D.Double(roi.getCentroidX(), roi.getCentroidY())
            double distance = centroid.distance(roiCentroid)
            if (distance < minDistance && distance < distThreshold) {
                minDistance = distance
                nearestROI = roi
            }
        }

        return nearestROI
    }

    @Override
    void run() {
        // Load whole cell mask image
        def wholeCellImp = IJ.openImage(wholeCellMaskFilePath)
        println 'Loaded whole cell mask width: ' + wholeCellImp.getWidth()

        // Load nuclear mask image
        def nuclearImp = IJ.openImage(nuclearMaskFilePath)
        println 'Loaded nuclear mask width: ' + nuclearImp.getWidth()

        // Build a server with supplied TIFF file
        def uri = new File(tiffFilePath).toURI()
        def builder = new BioFormatsServerBuilder()
        def server = builder.buildServer(uri)

        // Extract ROIs from whole cell and nuclear masks
        def wholeCellROIs = extractROIs(wholeCellImp, downsampleFactor)
        def nuclearROIs = extractROIs(nuclearImp, downsampleFactor)

        //[wholeCellROIs,nuclearROIs].transpose().collect { a, b -> println a ; println b }
        println 'Total whole cell ROIs: ' + wholeCellROIs.size()
        println 'Total nuclear ROIs: ' + nuclearROIs.size()

        // Convert QuPath ROIs to objects and add them to the hierarchy
        def pathObjects = createCellObjects(wholeCellROIs, nuclearROIs, distThreshold)
        println 'Total path objects: ' + pathObjects.size()

        if (!skipMeasurements) {
            // Set the pixel calibration
            PixelCalibration cal = new PixelCalibration.Builder()
                .pixelSizeMicrons(pixelSizeMicrons, pixelSizeMicrons)
                .build()
            println 'Set pixel calibration: ' + cal

            println 'Adding cell measurements...'
            ObjectMeasurements.addShapeMeasurements(
                pathObjects,
                cal,
                ObjectMeasurements.ShapeFeatures.AREA,
                ObjectMeasurements.ShapeFeatures.CIRCULARITY,
                ObjectMeasurements.ShapeFeatures.LENGTH,
                ObjectMeasurements.ShapeFeatures.MAX_DIAMETER,
                ObjectMeasurements.ShapeFeatures.MIN_DIAMETER,
                ObjectMeasurements.ShapeFeatures.NUCLEUS_CELL_RATIO,
                ObjectMeasurements.ShapeFeatures.SOLIDITY
            )

            // Define measurements
            def measurements = [
                ObjectMeasurements.Measurements.MEAN,
                ObjectMeasurements.Measurements.MEDIAN,
                ObjectMeasurements.Measurements.MIN,
                ObjectMeasurements.Measurements.MAX,
                ObjectMeasurements.Measurements.STD_DEV
            ]

            // Define compartments
            def compartments = [
                ObjectMeasurements.Compartments.CELL,
                ObjectMeasurements.Compartments.CYTOPLASM,
                ObjectMeasurements.Compartments.MEMBRANE,
                ObjectMeasurements.Compartments.NUCLEUS
            ]

            println 'Adding intensity measurements...'
            // Add intensity measurements
            pathObjects.each { pathObject ->
                ObjectMeasurements.addIntensityMeasurements(
                    server,
                    pathObject,
                    downsampleFactor,
                    measurements,
                    compartments
                )
            }
        }

        // Create a top-level annotation object for the whole image
        def width = server.getWidth()
        def height = server.getHeight()
        def roi = ROIs.createRectangleROI(0, 0, width, height, null)
        def annotation = PathObjects.createAnnotationObject(roi)

        // Add the annotation object to the start of the pathObjects list
        pathObjects.add(0, annotation)

        // Export the objects to GeoJSON
        QP.exportObjectsToGeoJson(
            pathObjects,
            outputFilePath,
            GeoJsonExportOptions.PRETTY_JSON,
            GeoJsonExportOptions.FEATURE_COLLECTION
        )
    }

    static void main(String[] args) {
        int exitCode = new CommandLine(new App()).execute(args)
        System.exit(exitCode)
    }

}

