package wehisodahub.cellmeasurement

import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Option

import java.awt.geom.Point2D
import java.nio.file.Paths

import ij.IJ
import ij.process.ColorProcessor

import qupath.imagej.processing.RoiLabeling
import qupath.imagej.tools.IJTools

import qupath.lib.roi.interfaces.ROI
import qupath.lib.roi.ROIs
import qupath.lib.roi.GeometryTools
import qupath.lib.scripting.QP
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.objects.CellTools
import qupath.lib.regions.ImagePlane
import qupath.lib.io.PathIO.GeoJsonExportOptions
import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.images.servers.PixelCalibration
import qupath.lib.images.servers.bioformats.BioFormatsServerBuilder

/**
* Entry point for the cell measurement application.
* This application takes nuclear and whole-cell segmentation masks, matches the nuclei
* to cells and uses the QuPath API to add cell shape and intensity measurements.
*/
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
    BigDecimal downsampleFactor = 1.0

    @Option(names = ['-p', '--pixel-size-microns'],
            description = 'Pixel size in microns (default: 0.5)',
            required = false)
    BigDecimal pixelSizeMicrons = 0.5

    @Option(names = ['--skip-measurements'],
            description = 'Skip adding measurements',
            required = false)
    boolean skipMeasurements = false

    @Option(names = ['-i', '--dist-threshold'],
            description = 'Distance threshold (in pixels) for matching ROIs',
            required = false)
    BigDecimal distThreshold = 10.0

    @Option(names = ['-e', '--cell-expansion'],
            description = 'Expansion factor for cell boundary estimation in pixels (default = 3.0)',
            required = false)
    BigDecimal cellExpansion = 3.0

    /**
    * Extract ROIs from a binary mask image.
    */
    static List<ROI> extractROIs(image, downsampleFactor) {
        def ip = image.getProcessor()
        if (ColorProcessor.class.isAssignableFrom(ip.getClass())) {
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
            if (it == null) { return }
            return IJTools.convertToROI(it, 0, 0, downsampleFactor, ImagePlane.getDefaultPlane())
        }.findAll { it != null }
    }

    /**
    * Create cell objects from matched nuclear and whole cell ROIs.
    */
    static List<PathObject> makeCellObjects(List<ROI> wholeCellROIs, List<ROI> nuclearROIs,
                                              BigDecimal distThreshold, BigDecimal cellExpansion) {
        def matchedPairs = matchROIs(nuclearROIs, wholeCellROIs, distThreshold, cellExpansion)
        def pathObjects = matchedPairs.collect { nucleus, cell ->
            if (cell != null) {
                return PathObjects.createCellObject(cell, nucleus)
            }
        }.findAll { it != null }
        return CellTools.constrainCellOverlaps(pathObjects)
        //return pathObjects
    }

    /**
    * Match nuclear ROIs to whole cell ROIs based on distance between centroids.
    * Estimate cell boundaries for unmatched nuclear ROIs with cell expansion
    */
    static List<List<ROI>> matchROIs(List<ROI> nuclearROIs, List<ROI> wholeCellROIs,
                                     BigDecimal distThreshold, BigDecimal cellExpansion) {
        def matchedPairs = []

        nuclearROIs.each { nuclearROI ->
            def nuclearCentroid = new Point2D.Double(nuclearROI.getCentroidX(), nuclearROI.getCentroidY())
            def nearestCell = findNearestROI(nuclearCentroid, wholeCellROIs, distThreshold)
            if (nearestCell == null) {
                def geom = CellTools.estimateCellBoundary(nuclearROI.getGeometry(), cellExpansion, 1.0)
                nearestCell = GeometryTools.geometryToROI(geom, nuclearROI.getImagePlane())
            }
            matchedPairs << [nuclearROI, nearestCell]
        }

        return matchedPairs
    }

    /**
    * Find the nearest ROI to a given centroid within a list of ROIs.
    */
    static ROI findNearestROI(Point2D centroid, List<ROI> rois, BigDecimal distThreshold) {
        ROI nearestROI = null
        BigDecimal minDistance = Double.MAX_VALUE

        rois.each { roi ->
            def roiCentroid = new Point2D.Double(roi.getCentroidX(), roi.getCentroidY())
            BigDecimal distance = centroid.distance(roiCentroid)
            if (distance < minDistance && distance < distThreshold) {
                minDistance = distance
                nearestROI = roi
            }
        }

        return nearestROI
    }

    /**
    * Filter out cells that have a membrane outside of the image bounds.
    */
    static List<PathObject> removeOutOfBoundsCells(List<PathObject> pathObjects, int imageWidth, int imageHeight) {
        return pathObjects.findAll { cell ->
            def roi = cell.getROI()
            def boundsX = roi.getBoundsX()
            def boundsY = roi.getBoundsY()
            def boundsWidth = roi.getBoundsWidth()
            def boundsHeight = roi.getBoundsHeight()
            return boundsX >= 0 &&
                   boundsY >= 0 &&
                   boundsX + boundsWidth <= imageWidth &&
                   boundsY + boundsHeight <= imageHeight
        }
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
        def uri = Paths.get(tiffFilePath).toUri()
        def builder = new BioFormatsServerBuilder()
        def server = builder.buildServer(uri)

        // Extract ROIs from whole cell and nuclear masks
        def wholeCellROIs = extractROIs(wholeCellImp, downsampleFactor)
        def nuclearROIs = extractROIs(nuclearImp, downsampleFactor)

        //[wholeCellROIs,nuclearROIs].transpose().collect { a, b -> println a ; println b }
        println 'Total whole cell ROIs: ' + wholeCellROIs.size()
        println 'Total nuclear ROIs: ' + nuclearROIs.size()

        // Convert QuPath ROIs to objects and add them to the hierarchy
        def pathObjects = makeCellObjects(wholeCellROIs, nuclearROIs, distThreshold, cellExpansion)
        println 'Total path objects: ' + pathObjects.size()

        // Filter out any cells that have a membrane outside of the image bounds
        def imageWidth = server.getWidth()
        def imageHeight = server.getHeight()
        pathObjects = removeOutOfBoundsCells(pathObjects, imageWidth, imageHeight)

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
        def roi = ROIs.createRectangleROI(0, 0, imageWidth, imageHeight, null)
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
        if (exitCode != 0) {
            println "Application exited with code: $exitCode"
        }
    }

}

