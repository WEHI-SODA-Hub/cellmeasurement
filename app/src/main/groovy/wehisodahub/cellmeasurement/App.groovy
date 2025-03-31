package wehisodahub.cellmeasurement

import ij.IJ
import ij.process.ColorProcessor
import qupath.lib.roi.interfaces.ROI
import qupath.lib.roi.ROIs
import qupath.lib.scripting.QP
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.objects.hierarchy.PathObjectHierarchy
import qupath.lib.regions.ImagePlane
import qupath.lib.io.PathIO.GeoJsonExportOptions
import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.images.servers.PixelCalibration
import qupath.imagej.processing.RoiLabeling
import qupath.imagej.tools.IJTools
import qupath.lib.images.servers.bioformats.BioFormatsServerBuilder
import java.awt.geom.Point2D

class App {

    // Extract ROIs from a given image
    static List<ROI> extractROIs(image, downsampleFactor) {
        def ip = image.getProcessor()
        if (ip instanceof ColorProcessor) {
            throw new IllegalArgumentException("RGB images are not supported!")
        }

        int n = image.getStatistics().max as int
        if (n == 0) {
            println 'No objects found in mask!'
            return []
        }

        // Convert mask to ROIs
        def roisIJ = RoiLabeling.labelsToConnectedROIs(ip, n)
        println "Number of ROIs found: " + roisIJ.size()

        return roisIJ.collect {
            if (it == null) return
            return IJTools.convertToROI(it, 0, 0, downsampleFactor, ImagePlane.getDefaultPlane())
        }.findAll { it != null }
    }

    static List<PathObject> createCellObjects(List<ROI> wholeCellROIs, List<ROI> nuclearROIs) {
        // TODO: we need to handle the case where the whole cell and nuclear masks
        // have different numbers of ROIs
        def matchedPairs = matchROIs(wholeCellROIs, nuclearROIs)
        return matchedPairs.collect { cell, nucleus ->
            if (cell != null) {
                return PathObjects.createCellObject(cell, nucleus)
            }
        }.findAll { it != null }
    }

    static List<List<ROI>> matchROIs(List<ROI> wholeCellROIs, List<ROI> nuclearROIs) {
        def matchedPairs = []

        wholeCellROIs.each { cellROI ->
            def cellCentroid = new Point2D.Double(cellROI.getCentroidX(), cellROI.getCentroidY())
            def nearestNucleus = findNearestROI(cellCentroid, nuclearROIs)
            matchedPairs << [cellROI, nearestNucleus]
        }

        return matchedPairs
    }

    static ROI findNearestROI(Point2D centroid, List<ROI> rois) {
        ROI nearestROI = null
        double minDistance = Double.MAX_VALUE

        rois.each { roi ->
            def roiCentroid = new Point2D.Double(roi.getCentroidX(), roi.getCentroidY())
            double distance = centroid.distance(roiCentroid)
            if (distance < minDistance) {
                minDistance = distance
                nearestROI = roi
            }
        }

        return nearestROI
    }

    static void main(String[] args) {
        // TODO: better arg parsing
        if (args.length < 6) {
            println "Usage: gradlew run --args=\"<wholeCellMaskFilePath> <nuclearMaskFilePath> <tiffFilePath> <outputFilePath> <downsampleFactor> <pixelSizeMicrons>\""
            return
        }

        // TODO: input validation
        // Input file paths
        def wholeCellMaskFilePath = args[0]
        def nuclearMaskFilePath = args[1]
        def tiffFilePath = args[2]
        def outputFilePath = args[3]
        def downsampleFactor = args[4].toDouble()
        def pixelSizeMicrons = args[5].toDouble()

        // Load whole cell mask image
        def wholeCellImp = IJ.openImage(wholeCellMaskFilePath)
        println "Loaded whole cell mask width: " + wholeCellImp.getWidth()

        // Load nuclear mask image
        def nuclearImp = IJ.openImage(nuclearMaskFilePath)
        println "Loaded nuclear mask width: " + nuclearImp.getWidth()

        // Extract ROIs from whole cell and nuclear masks
        def wholeCellROIs = extractROIs(wholeCellImp, downsampleFactor)
        def nuclearROIs = extractROIs(nuclearImp, downsampleFactor)

        //[wholeCellROIs,nuclearROIs].transpose().collect { a, b -> println a ; println b }
        println "Total whole cell ROIs: " + wholeCellROIs.size()
        println "Total nuclear ROIs: " + nuclearROIs.size()

        // Convert QuPath ROIs to objects and add them to the hierarchy
        def pathObjects = createCellObjects(wholeCellROIs, nuclearROIs)
        println "Total path objects: " + pathObjects.size()

        // Set the pixel calibration
        PixelCalibration cal = new PixelCalibration.Builder()
            .pixelSizeMicrons(pixelSizeMicrons, pixelSizeMicrons)
            .build()
        println "Set pixel calibration: " + cal

        // Add cell measurements
        println "Adding cell measurements..."
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

        // Build a server with supplied TIFF file
        def uri = new File(tiffFilePath).toURI()
        def builder = new BioFormatsServerBuilder()
        def server = builder.buildServer(uri)

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

        println "Adding intensity measurements..."
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

}

