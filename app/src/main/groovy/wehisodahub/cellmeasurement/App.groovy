package wehisodahub.cellmeasurement

import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Option

import groovyx.gpars.GParsPool

import java.awt.geom.Point2D
import java.nio.file.Paths
import java.awt.Rectangle
import java.awt.image.BufferedImage

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

import ij.IJ
import ij.process.ColorProcessor
import ij.process.ByteProcessor
import ij.process.ImageProcessor

import qupath.imagej.processing.RoiLabeling
import qupath.imagej.tools.IJTools

import qupath.lib.roi.RectangleROI
import qupath.lib.roi.interfaces.ROI
import qupath.lib.roi.ROIs
import qupath.lib.roi.GeometryTools
import qupath.lib.scripting.QP
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.objects.CellTools
import qupath.lib.regions.ImagePlane
import qupath.lib.regions.RegionRequest
import qupath.lib.io.PathIO.GeoJsonExportOptions
import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.measurements.MeasurementList
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

    @Option(names = ['-f', '--tiff-file'],
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

    @Option(names = ['--percentiles'],
            description = 'Calculate intensity percentiles. Only works if not skipping measurements.',
            required = false)
    boolean percentiles = false

    @Option(names = ['-i', '--dist-threshold'],
            description = 'Distance threshold (in pixels) for matching ROIs',
            required = false)
    BigDecimal distThreshold = 10.0

    @Option(names = ['-e', '--estimate-cell-boundary-dist'],
            description = 'Where no matching membrane ROI exists, expand the nucleus by this many pixels (default = 3.0)',
            required = false)
    BigDecimal estimateCellBoundaryDist = 3.0

    @Option(names = ['-t', '--threads'],
            description = 'Number of threads to use for parallel processing (default: 1)',
            required = false)
    int threads = 1

    /**
    * Extract ROIs from a binary mask image.
    */
    static List<ROI> extractROIs(image, downsampleFactor, threads = 1) {
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

        GParsPool.withPool(threads) {
            return roisIJ.collectParallel {
                if (it == null) { return }
                return IJTools.convertToROI(it, 0, 0, downsampleFactor, ImagePlane.getDefaultPlane())
            }.findAll { it != null }
        }
    }

    /**
    * Create cell objects from matched nuclear and whole cell ROIs.
    */
    static List<PathObject> makeCellObjects(List<ROI> wholeCellROIs, List<ROI> nuclearROIs,
                                              BigDecimal distThreshold, BigDecimal estimateCellBoundaryDist) {
        def matchedPairs = matchROIs(nuclearROIs, wholeCellROIs, distThreshold, estimateCellBoundaryDist)
        def pathObjects = matchedPairs.collect { nucleus, cell ->
            if (cell != null) {
                return PathObjects.createCellObject(cell, nucleus)
            }
        }.findAll { it != null }
        return CellTools.constrainCellOverlaps(pathObjects)
    }

    /**
    * Match nuclear ROIs to whole cell ROIs based on distance between centroids.
    * Estimate cell boundaries for unmatched nuclear ROIs with cell expansion
    */
    static List<List<ROI>> matchROIs(List<ROI> nuclearROIs, List<ROI> wholeCellROIs,
                                     BigDecimal distThreshold, BigDecimal estimateCellBoundaryDist,
                                     threads = 1) {
        GParsPool.withPool(threads) {
            nuclearROIs.collectParallel { nuclearROI ->
                def nuclearCentroid = new Point2D.Double(nuclearROI.getCentroidX(), nuclearROI.getCentroidY())
                def nearestCell = findNearestROI(nuclearCentroid, wholeCellROIs, distThreshold)
                if (nearestCell == null) {
                    def geom = CellTools.estimateCellBoundary(nuclearROI.getGeometry(), estimateCellBoundaryDist, 1.0)
                    nearestCell = GeometryTools.geometryToROI(geom, nuclearROI.getImagePlane())
                }
                return [nuclearROI, nearestCell]
            }
        }
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

    /**
     * Add percentile measurements for cell objects by compartment
     * @param server ImageServer containing the pixel data
     * @param pathObject PathObject to measure (MeasurementList will be updated)
     * @param downsampleFactor Resolution at which to request pixels
     * @param percentiles List of percentiles to calculate (default: [70, 80, 90, 95, 96, 97, 98, 99])
     * @param compartments Set of compartments to measure ('NUCLEUS', 'CYTOPLASM', 'MEMBRANE', 'CELL')
     */
    def addPercentileMeasurements(server, PathObject pathObject, double downsampleFactor = 1.0,
                                  List<Double> percentiles = [70, 80, 90, 95, 96, 97, 98, 99],
                                  Set<String> compartments = ['NUCLEUS', 'CYTOPLASM', 'MEMBRANE', 'CELL']) {

        // Only process cells
        if (!pathObject.isCell()) {
            return
        }

        def cell = pathObject
        def cellROI = cell.getROI()
        def nucleusROI = cell.getNucleusROI()

        if (cellROI == null) {
            return
        }

        // Get bounding box for the cell
        def bounds = cellROI.getBoundsX() != Double.POSITIVE_INFINITY ?
                     new RectangleROI((cellROI.getBoundsX() / downsampleFactor),
                                      (cellROI.getBoundsY() / downsampleFactor),
                                      (cellROI.getBoundsWidth() / downsampleFactor + 1),
                                      (cellROI.getBoundsHeight() / downsampleFactor + 1)) :
                     new RectangleROI(0, 0, server.getWidth(), server.getHeight())

        // Create region request
        def request = RegionRequest.createInstance(server.getPath(), downsampleFactor, bounds)

        try {
            // Get the image data
            def img = server.readRegion(request)
            if (img == null) return

            int width = img.getWidth()
            int height = img.getHeight()
            int nChannels = server.nChannels()

            // Create compartment masks
            def pathImage = IJTools.convertToImagePlus(server, request);
            def masks = createCompartmentMasks(cellROI, nucleusROI, width, height, pathImage, compartments)

            // Extract pixel values for each channel and compartment
            def measurements = cell.getMeasurementList()

            for (int c = 0; c < nChannels; c++) {
                def channelName = server.getChannel(c).getName()
                def pixelValues = extractChannelPixels(img, c, width, height)

                masks.each { compartment, mask ->
                    if (compartments.contains(compartment)) {
                        def compartmentPixels = getCompartmentPixels(pixelValues, mask)
                        if (compartmentPixels.size() > 0) {
                            addPercentileMeasurementsForCompartment(measurements, compartmentPixels,
                                                                    channelName, compartment, percentiles)
                        }
                    }
                }
            }

        } catch (Exception e) {
            println("Error processing ${pathObject}: ${e.getMessage()}")
        }
    }

    /**
     * Create binary masks for different cell compartments
     */
    def createCompartmentMasks(cellROI, nucleusROI, int width, int height, pathImage, compartments) {
        def masks = [:]

        // Create cell mask
        if (compartments.contains('CELL')) {
            masks['CELL'] = createROIMask(cellROI, width, height, pathImage)
        }

        // Create nucleus mask
        def nucleusMask = null
        if (nucleusROI != null && (compartments.contains('NUCLEUS') || compartments.contains('CYTOPLASM'))) {
            nucleusMask = createROIMask(nucleusROI, width, height, pathImage)
            if (compartments.contains('NUCLEUS')) {
                masks['NUCLEUS'] = nucleusMask
            }
        }

        // Create cytoplasm mask (cell - nucleus)
        if (compartments.contains('CYTOPLASM') && masks.containsKey('CELL') && nucleusMask != null) {
            masks['CYTOPLASM'] = subtractMasks(masks['CELL'], nucleusMask)
        }

        // Create membrane mask (cell boundary)
        if (compartments.contains('MEMBRANE') && masks.containsKey('CELL')) {
            masks['MEMBRANE'] = createMembraneMask(masks['CELL'])
        }

        return masks
    }

    /**
     * Create binary mask from ROI
     */
    def createROIMask(roi, int width, int height, pathImage) {
        def mask = new ByteProcessor(width, height)

        def roiIJ = IJTools.convertToIJRoi(roi, pathImage)

        mask.setColor(255)
        mask.fill(roiIJ)

        return mask
    }

    /**
     * Subtract one mask from another
     */
    def subtractMasks(mask1, mask2) {
        def result = mask1.duplicate()
        def pixels1 = result.getPixels() as byte[]
        def pixels2 = mask2.getPixels() as byte[]

        for (int i = 0; i < pixels1.length; i++) {
            if (pixels2[i] != 0) {
                pixels1[i] = 0
            }
        }

        return result
    }

    /**
     * Create membrane mask by finding boundary pixels
     */
    def createMembraneMask(cellMask) {
        def membrane = cellMask.duplicate()
        membrane.findEdges()
        return membrane
    }

    /**
     * Extract pixel values for a specific channel
     */
    def extractChannelPixels(BufferedImage img, int channel, int width, int height) {
        def pixels = new float[width * height]
        def raster = img.getRaster()

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                pixels[y * width + x] = raster.getSampleFloat(x, y, channel)
            }
        }

        return pixels
    }

    /**
     * Get pixel values within a compartment mask
     */
    def getCompartmentPixels(float[] allPixels, ByteProcessor mask) {
        def maskPixels = mask.getPixels() as byte[]
        def compartmentPixels = []

        for (int i = 0; i < maskPixels.length; i++) {
            if (maskPixels[i] != 0) {
                compartmentPixels.add(allPixels[i] as double)
            }
        }

        return compartmentPixels
    }

    /**
     * Calculate and add percentile measurements for a compartment
     */
    def addPercentileMeasurementsForCompartment(MeasurementList measurements, List<Double> pixels,
                                               String channelName, String compartment, List<Double> percentiles) {
        def stats = new DescriptiveStatistics()
        pixels.each { stats.addValue(it) }

        percentiles.each { percentile ->
            def value = stats.getPercentile(percentile)
            def capitalisedCompartment = compartment.toLowerCase().capitalize()
            def measurementName = "${channelName}: ${capitalisedCompartment}: Percentile: ${percentile}"
            measurements.putMeasurement(measurementName, value)
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
        println 'Extracting ROIs...'
        def wholeCellROIs = extractROIs(wholeCellImp, downsampleFactor, threads)
        def nuclearROIs = extractROIs(nuclearImp, downsampleFactor, threads)

        //[wholeCellROIs,nuclearROIs].transpose().collect { a, b -> println a ; println b }
        println 'Total whole cell ROIs: ' + wholeCellROIs.size()
        println 'Total nuclear ROIs: ' + nuclearROIs.size()

        // Convert QuPath ROIs to objects and add them to the hierarchy
        def pathObjects = makeCellObjects(wholeCellROIs, nuclearROIs, distThreshold, estimateCellBoundaryDist)
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
            GParsPool.withPool(threads) {
                // Add cell shape measurements
                pathObjects.eachParallel { pathObject ->
                    ObjectMeasurements.addShapeMeasurements(
                        pathObject,
                        cal,
                        ObjectMeasurements.ShapeFeatures.AREA,
                        ObjectMeasurements.ShapeFeatures.CIRCULARITY,
                        ObjectMeasurements.ShapeFeatures.LENGTH,
                        ObjectMeasurements.ShapeFeatures.MAX_DIAMETER,
                        ObjectMeasurements.ShapeFeatures.MIN_DIAMETER,
                        ObjectMeasurements.ShapeFeatures.NUCLEUS_CELL_RATIO,
                        ObjectMeasurements.ShapeFeatures.SOLIDITY
                    )
                }
            }

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
            GParsPool.withPool(threads) {
                // Add intensity measurements
                pathObjects.eachParallel { pathObject ->
                    ObjectMeasurements.addIntensityMeasurements(
                        server,
                        pathObject,
                        downsampleFactor,
                        measurements,
                        compartments
                    )
                }
            }

            if (percentiles) {
                println 'Adding intensity percentiles...'
                GParsPool.withPool(threads) {
                    pathObjects.eachParallel { pathObject ->
                        addPercentileMeasurements(
                            server,
                            pathObject,
                            downsampleFactor
                        )
                    }
                }
            }
        }

        // Create a top-level annotation object for the whole image
        def roi = ROIs.createRectangleROI(0, 0, imageWidth, imageHeight, null)
        def annotation = PathObjects.createAnnotationObject(roi)

        // Add the annotation object to the start of the pathObjects list
        pathObjects.add(0, annotation)

        println 'Exporting to GeoJSON...'
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

