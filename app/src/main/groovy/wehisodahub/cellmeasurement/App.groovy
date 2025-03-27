package wehisodahub.cellmeasurement

import ij.IJ
import ij.process.ColorProcessor
import qupath.lib.roi.interfaces.ROI
import qupath.lib.scripting.QP
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.objects.hierarchy.PathObjectHierarchy
import qupath.lib.regions.ImagePlane
import qupath.imagej.processing.RoiLabeling
import qupath.imagej.tools.IJTools
import qupath.lib.io.PathIO.GeoJsonExportOptions
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
        if (args.length < 3) {
            println "Usage: gradlew run --args=\"<wholeCellMaskFilePath> <nuclearMaskFilePath> <outputFilePath>\""
            return
        }

        // TODO: input validation
        // Input file paths
        def wholeCellMaskFilePath = args[0]
        def nuclearMaskFilePath = args[1]
        def outputFilePath = args[2]

        // Load whole cell mask image
        def wholeCellImp = IJ.openImage(wholeCellMaskFilePath)
        println "Loaded whole cell mask width: " + wholeCellImp.getWidth()

        // Load nuclear mask image
        def nuclearImp = IJ.openImage(nuclearMaskFilePath)
        println "Loaded nuclear mask width: " + nuclearImp.getWidth()

        // Extract ROIs from whole cell and nuclear masks
        def downsampleFactor = 1.0
        def wholeCellROIs = extractROIs(wholeCellImp, downsampleFactor)
        def nuclearROIs = extractROIs(nuclearImp, downsampleFactor)

        //[wholeCellROIs,nuclearROIs].transpose().collect { a, b -> println a ; println b }
        println "Total whole cell ROIs: " + wholeCellROIs.size()
        println "Total nuclear ROIs: " + nuclearROIs.size()

        // Convert QuPath ROIs to objects and add them to the hierarchy
        def pathObjects = createCellObjects(wholeCellROIs, nuclearROIs)
        println "Total path objects: " + pathObjects.size()
         
        // Export the objects to GeoJSON
        QP.exportObjectsToGeoJson(
            pathObjects,
            outputFilePath,
            GeoJsonExportOptions.PRETTY_JSON,
            GeoJsonExportOptions.FEATURE_COLLECTION
        )
    }

}

