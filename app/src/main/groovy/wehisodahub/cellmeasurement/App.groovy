package wehisodahub.cellmeasurement

import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Option

import groovyx.gpars.GParsPool

import java.awt.geom.Point2D
import java.nio.file.Paths
import java.awt.Rectangle
import java.awt.image.BufferedImage
import java.util.concurrent.atomic.AtomicInteger

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

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
import qupath.lib.objects.PathObjectTools
import qupath.lib.regions.ImagePlane
import qupath.lib.regions.RegionRequest
import qupath.lib.io.PathIO.GeoJsonExportOptions
import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.measurements.MeasurementList
import qupath.lib.images.servers.PixelCalibration
import qupath.lib.images.servers.bioformats.BioFormatsServerBuilder

import org.locationtech.jts.geom.Envelope
import org.locationtech.jts.index.strtree.STRtree
import org.locationtech.jts.simplify.DouglasPeuckerSimplifier

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

    @Option(names = ['--simplify-rois'],
            description = 'Simplify ROIs, can be used to speed up processing of complex masks',
            required = false)
    boolean simplifyROIs = false

    @Option(names = ['--tolerance'],
            description = 'Use this tolerance value when simplifying ROIs, must be greater than 0. Default: 0.5',
            required = false)
    BigDecimal simplifyTolerance = 0.5

    @Option(names = ['--percentiles'],
            description = 'Calculate specified comma-separated intensity percentiles. Only works if not skipping measurements. E.g. "70,80,90,95,96,97,98,99"',
            required = false)
    String percentiles = ''

    @Option(names = ['--erosion-steps'],
            description = 'Comma-separated erosion steps in pixels. E.g. "4,7,11,14,18". Set to empty to disable.',
            required = false)
    String erosionSteps = ''

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

    @Option(names = ['--tile-size'],
            description = 'Tile size in pixels for batch image I/O during measurement (default: 2048)',
            required = false)
    int tileSize = 2048

    @Option(names = ['--tile-overlap'],
            description = 'Tile overlap in pixels to reduce edge overflows (default: 200)',
            required = false)
    int tileOverlap = 200

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
                                              BigDecimal distThreshold, BigDecimal estimateCellBoundaryDist,
                                              threads = 1) {
        def matchedPairs = matchROIs(
            nuclearROIs, wholeCellROIs, distThreshold, estimateCellBoundaryDist, threads
        )
        def pathObjects = matchedPairs.collect { nucleus, cell ->
            if (cell != null) {
                return PathObjects.createCellObject(cell, nucleus)
            }
        }.findAll { it != null }
        def constrained = CellTools.constrainCellOverlaps(pathObjects)

        // constrainCellOverlaps uses Voronoi clipping which can produce GeometryCollections
        // with nested MultiPolygons. These cause infinite recursion in QuPath's GeoJSON writer
        // (ROITypeAdapters.writeCoordinates), which dispatches on top-level type but then
        // recurses back through the generic dispatcher for each child — so a MultiPolygon
        // whose children are themselves MultiPolygons (not flat Polygons) loops forever.
        // Applying buffer(0) re-builds each geometry from scratch via JTS's sweep-line
        // algorithm, guaranteeing a flat, valid Polygon or MultiPolygon(Polygon...) structure.
        // Cells that are still degenerate (zero area) after fixing are dropped.
        int fixed = 0
        int dropped = 0
        def valid = constrained.collect { cell ->
            def roi = cell.getROI()
            if (roi == null) { dropped++; return null }
            def geom = roi.getGeometry()
            if (geom == null || geom.isEmpty()) { dropped++; return null }

            if (!isWritableGeometry(geom)) {
                // Flatten complex/nested geometry via buffer(0)
                def fixedGeom = geom.buffer(0)
                if (fixedGeom.isEmpty() || fixedGeom.getArea() <= 0) { dropped++; return null }
                def fixedRoi = GeometryTools.geometryToROI(fixedGeom, roi.getImagePlane())
                def nucleusRoi = cell.getNucleusROI()
                fixed++
                return PathObjects.createCellObject(fixedRoi, nucleusRoi)
            }

            if (geom.getArea() <= 0) { dropped++; return null }
            return cell
        }.findAll { it != null }

        if (fixed > 0)
            println "Fixed ${fixed} cell(s) with nested/complex geometry using buffer(0)"
        if (dropped > 0)
            println "Dropped ${dropped} degenerate cell(s) produced by constrainCellOverlaps (zero area or non-polygon geometry)"
        return valid
    }

    /**
    * Returns true if the geometry is safe to pass to QuPath's GeoJSON writer.
    * The writer's writeCoordinates dispatches on type but recurses through a generic
    * dispatcher for children, so a MultiPolygon whose children are not flat Polygons
    * (e.g. nested MultiPolygons from GeometryCollection coercion) causes infinite recursion.
    */
    static boolean isWritableGeometry(org.locationtech.jts.geom.Geometry geom) {
        def type = geom.getGeometryType()
        if (type == 'Polygon') return true
        if (type == 'MultiPolygon') {
            for (int i = 0; i < geom.getNumGeometries(); i++) {
                if (geom.getGeometryN(i).getGeometryType() != 'Polygon') return false
            }
            return true
        }
        return false
    }

    /**
    * Match nuclear ROIs to whole cell ROIs based on distance between centroids.
    * Estimate cell boundaries for unmatched nuclear ROIs with cell expansion.
    * Uses an STRtree spatial index over whole-cell centroids to avoid an O(N×M)
    * linear scan — each nucleus queries only candidates within distThreshold.
    */
    static List<List<ROI>> matchROIs(List<ROI> nuclearROIs, List<ROI> wholeCellROIs,
                                     BigDecimal distThreshold, BigDecimal estimateCellBoundaryDist,
                                     threads = 1) {
        // Build spatial index once over whole-cell centroids — O(M log M), read-only so thread-safe
        def index = new STRtree()
        wholeCellROIs.each { roi ->
            double cx = roi.getCentroidX()
            double cy = roi.getCentroidY()
            index.insert(new Envelope(cx, cx, cy, cy), roi)
        }
        index.build()

        GParsPool.withPool(threads) {
            nuclearROIs.collectParallel { nuclearROI ->
                def nearestCell = findNearestROI(nuclearROI, index, distThreshold as double)
                if (nearestCell == null) {
                    def geom = CellTools.estimateCellBoundary(
                        nuclearROI.getGeometry(), estimateCellBoundaryDist, 1.0
                    )
                    nearestCell = GeometryTools.geometryToROI(geom, nuclearROI.getImagePlane())
                }
                return [nuclearROI, nearestCell]
            }
        }
    }

    /**
    * Find the nearest ROI to a given nuclear ROI's centroid using a pre-built STRtree spatial index.
    * Queries only candidates within distThreshold bounding box, then applies exact distance check.
    */
    static ROI findNearestROI(ROI nuclearROI, STRtree index, double distThreshold) {
        double cx = nuclearROI.getCentroidX()
        double cy = nuclearROI.getCentroidY()
        def searchEnv = new Envelope(cx - distThreshold, cx + distThreshold,
                                     cy - distThreshold, cy + distThreshold)
        def candidates = index.query(searchEnv)

        ROI nearestROI = null
        double minDistance = Double.MAX_VALUE
        candidates.each { roi ->
            double distance = Math.sqrt(
                Math.pow(cx - roi.getCentroidX(), 2) + Math.pow(cy - roi.getCentroidY(), 2)
            )
            if (distance < minDistance && distance < distThreshold) {
                minDistance = distance
                nearestROI = roi
            }
        }
        return nearestROI
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
     * Create a bounding box RectangleROI for a cell, accounting for downsampling
     * @param cellROI The cell's region of interest
     * @param downsampleFactor Resolution factor to apply
     * @param server ImageServer to get full image dimensions if needed
     * @return RectangleROI representing the bounding box
     */
    static RectangleROI getCellBoundingBox(ROI cellROI, double downsampleFactor, server) {
        return cellROI.getBoundsX() != Double.POSITIVE_INFINITY ?
               new RectangleROI((cellROI.getBoundsX() / downsampleFactor),
                                (cellROI.getBoundsY() / downsampleFactor),
                                (cellROI.getBoundsWidth() / downsampleFactor + 1),
                                (cellROI.getBoundsHeight() / downsampleFactor + 1)) :
               new RectangleROI(0, 0, server.getWidth(), server.getHeight())
    }

    /**
     * Load image data for a cell once, to be shared across multiple measurement types.
     * Reads the image region, creates compartment masks, and pre-extracts all channel pixel arrays.
     * @param server ImageServer containing the pixel data
     * @param pathObject PathObject (must be a cell)
     * @param downsampleFactor Resolution at which to request pixels
     * @param compartments Set of compartment masks to create ('CELL', 'NUCLEUS', 'CYTOPLASM', 'MEMBRANE')
     * @return Map with keys: masks (Map), allChannelPixels (List<float[]>); or null if loading fails
     */
    def loadCellData(server, PathObject pathObject, double downsampleFactor,
                     Set<String> compartments = ['CELL', 'NUCLEUS', 'CYTOPLASM', 'MEMBRANE']) {
        if (!pathObject.isCell()) return null

        def cellROI = pathObject.getROI()
        def nucleusROI = pathObject.getNucleusROI()
        if (cellROI == null) return null

        def bounds = getCellBoundingBox(cellROI, downsampleFactor, server)
        def request = RegionRequest.createInstance(server.getPath(), downsampleFactor, bounds)
        def img = server.readRegion(request)
        if (img == null) return null

        int width = img.getWidth()
        int height = img.getHeight()
        def pathImage = IJTools.convertToImagePlus(server, request)
        def masks = createCompartmentMasks(cellROI, nucleusROI, width, height, pathImage, compartments)
        def allChannelPixels = (0..<server.nChannels()).collect {
            c -> extractChannelPixels(img, c, width, height)
        }

        return [masks: masks, allChannelPixels: allChannelPixels]
    }

    /**
     * Group path objects by the tile their centroid falls in.
     * Tiles are defined in full image coordinates; tileSize is the nominal tile edge length.
     * @return Map from [tileCol, tileRow] to list of path objects
     */
    static Map<List<Integer>, List<PathObject>> groupCellsByTile(
            List<PathObject> pathObjects, int tileSize, int imageWidth, int imageHeight) {
        def tileMap = [:].withDefault { [] }
        pathObjects.each { pathObject ->
            def roi = pathObject.getROI()
            if (roi == null) return
            int tileCol = (int)(roi.getCentroidX() / tileSize)
            int tileRow = (int)(roi.getCentroidY() / tileSize)
            tileMap[[tileCol, tileRow]] << pathObject
        }
        return tileMap
    }

    /**
     * Extract a rectangular sub-region from a flat row-major pixel array.
     * @param tilePixels Full tile pixel array (row-major, tileWidth × tileHeight elements)
     * @param tileWidth  Width of the full tile in pixels
     * @param offX       X offset of sub-region within the tile (pixels)
     * @param offY       Y offset of sub-region within the tile (pixels)
     * @param subW       Width of sub-region (pixels)
     * @param subH       Height of sub-region (pixels)
     * @return New float[] with sub-region pixels in row-major order
     */
    static float[] extractSubRegionPixels(float[] tilePixels, int tileWidth,
                                           int offX, int offY, int subW, int subH) {
        float[] result = new float[subW * subH]
        for (int row = 0; row < subH; row++) {
            System.arraycopy(tilePixels, (offY + row) * tileWidth + offX,
                             result,     row * subW,
                             subW)
        }
        return result
    }

    /**
     * Create a binary mask for an ROI without requiring a PathImage, using explicit
     * origin coordinates. xOrigin and yOrigin are in the same coordinate space as
     * the ROI (full image coordinates); the resulting mask is in downsampled pixels.
     */
    def createROIMaskFromOrigin(roi, int width, int height,
                                 double xOrigin, double yOrigin, double downsampleFactor) {
        def mask = new ByteProcessor(width, height)
        try {
            def roiIJ = IJTools.convertToIJRoi(roi, xOrigin, yOrigin, downsampleFactor)
            mask.setColor(255)
            mask.fill(roiIJ)
        } catch (UnsupportedOperationException e) {
            println "Warning: could not convert ROI to ImageJ ROI, returning empty mask: ${e.getMessage()}"
        }
        return mask
    }

    /**
     * Create compartment masks without requiring a PathImage.
     * xOrigin/yOrigin are the full-image coordinates of the cell bounding box top-left corner.
     */
    def createCompartmentMasksFromOrigin(cellROI, nucleusROI, int width, int height,
                                          double xOrigin, double yOrigin, double downsampleFactor,
                                          Set<String> compartments) {
        def masks = [:]

        if (compartments.contains('CELL')) {
            masks['CELL'] = createROIMaskFromOrigin(cellROI, width, height, xOrigin, yOrigin, downsampleFactor)
        }

        def nucleusMask = null
        if (nucleusROI != null && (compartments.contains('NUCLEUS') || compartments.contains('CYTOPLASM'))) {
            nucleusMask = createROIMaskFromOrigin(nucleusROI, width, height, xOrigin, yOrigin, downsampleFactor)
            if (compartments.contains('NUCLEUS')) {
                masks['NUCLEUS'] = nucleusMask
            }
        }

        if (compartments.contains('CYTOPLASM') && masks.containsKey('CELL') && nucleusMask != null) {
            masks['CYTOPLASM'] = subtractMasks(masks['CELL'], nucleusMask)
        }

        if (compartments.contains('MEMBRANE') && masks.containsKey('CELL')) {
            masks['MEMBRANE'] = createMembraneMask(masks['CELL'])
        }

        return masks
    }

    /**
     * Load cell data from pre-fetched tile pixel arrays, avoiding per-cell image I/O.
     * Returns null if the cell's bounding box does not lie fully within the tile.
     * @param pathObject        Cell to load data for (must pass isCell())
     * @param downsampleFactor  Resolution factor used when reading the tile
     * @param compartments      Set of compartment masks to create
     * @param tileChannelPixels Pre-extracted pixel arrays per channel (tile pixel space, row-major)
     * @param tilePxX           Tile x origin in downsampled pixel coordinates
     * @param tilePxY           Tile y origin in downsampled pixel coordinates
     * @param tilePxW           Tile width in downsampled pixels
     * @param tilePxH           Tile height in downsampled pixels
     * @return Map with 'masks' and 'allChannelPixels', or null
     */
    def loadCellDataFromTile(PathObject pathObject, double downsampleFactor,
                              Set<String> compartments, List<float[]> tileChannelPixels,
                              int tilePxX, int tilePxY, int tilePxW, int tilePxH) {
        if (!pathObject.isCell()) return null

        def cellROI = pathObject.getROI()
        def nucleusROI = pathObject.getNucleusROI()
        if (cellROI == null) return null

        // Cell bounding box in downsampled pixel space
        int cPxX = (int)(cellROI.getBoundsX() / downsampleFactor)
        int cPxY = (int)(cellROI.getBoundsY() / downsampleFactor)
        int cPxW = (int)(cellROI.getBoundsWidth()  / downsampleFactor) + 1
        int cPxH = (int)(cellROI.getBoundsHeight() / downsampleFactor) + 1

        // Offset within the tile
        int offX = cPxX - tilePxX
        int offY = cPxY - tilePxY

        // Cell must lie fully within the tile
        if (offX < 0 || offY < 0 || offX + cPxW > tilePxW || offY + cPxH > tilePxH) return null

        // Sub-region pixels for each channel — no server read
        def allChannelPixels = tileChannelPixels.collect { tilePixels ->
            extractSubRegionPixels(tilePixels, tilePxW, offX, offY, cPxW, cPxH)
        }

        // xOrigin = -cellBoundsX/ds so that IJTools maps cell top-left → (0,0)
        double xOrigin = -cellROI.getBoundsX() / downsampleFactor
        double yOrigin = -cellROI.getBoundsY() / downsampleFactor
        def masks = createCompartmentMasksFromOrigin(
            cellROI, nucleusROI, cPxW, cPxH, xOrigin, yOrigin, downsampleFactor, compartments
        )

        return [masks: masks, allChannelPixels: allChannelPixels]
    }

    /**
     * Add standard intensity measurements (mean, median, min, max, std dev) for all
     * compartments from pre-loaded cell data. Replaces ObjectMeasurements.addIntensityMeasurements,
     * avoiding a redundant per-cell image read.
     * Measurement name format: "{Compartment}: {channel}: {stat}" — e.g. "Cell: DAPI: Mean"
     */
    def addStandardIntensityMeasurements(server, PathObject pathObject, Map cellData) {
        if (cellData == null) return
        def ml = pathObject.getMeasurementList()
        int nChannels = server.nChannels()
        def compartmentLabels = [CELL: 'Cell', NUCLEUS: 'Nucleus', CYTOPLASM: 'Cytoplasm', MEMBRANE: 'Membrane']

        for (int c = 0; c < nChannels; c++) {
            def channelName = server.getChannel(c).getName()
            def pixelValues = cellData.allChannelPixels[c]

            cellData.masks.each { compartment, mask ->
                def label = compartmentLabels.get(compartment, compartment)
                def pixels = getCompartmentPixels(pixelValues, mask)
                if (pixels.length == 0) return

                def stats = new DescriptiveStatistics(pixels)
                ml.put("${label}: ${channelName}: Mean",     stats.getMean())
                ml.put("${label}: ${channelName}: Median",   stats.getPercentile(50))
                ml.put("${label}: ${channelName}: Min",      stats.getMin())
                ml.put("${label}: ${channelName}: Max",      stats.getMax())
                ml.put("${label}: ${channelName}: Std.Dev.", stats.getStandardDeviation())
            }
        }
    }

    /**
     * Add percentile measurements for cell objects by compartment
     * @param server ImageServer containing the pixel data
     * @param pathObject PathObject to measure (MeasurementList will be updated)
     * @param downsampleFactor Resolution at which to request pixels
     * @param percentiles List of percentiles to calculate (default: [70, 80, 90, 95, 96, 97, 98, 99])
     * @param compartments Set of compartments to measure ('NUCLEUS', 'CYTOPLASM', 'MEMBRANE', 'CELL')
     * @param cellData Optional pre-loaded cell data from loadCellData(); if null, data is loaded internally
     */
    def addPercentileMeasurements(server, PathObject pathObject, double downsampleFactor = 1.0,
                                  List<Double> percentiles = [70, 80, 90, 95, 96, 97, 98, 99],
                                  Set<String> compartments = ['NUCLEUS', 'CYTOPLASM', 'MEMBRANE', 'CELL'],
                                  Map cellData = null) {

        // Only process cells
        if (!pathObject.isCell()) {
            return
        }

        if (pathObject.getROI() == null) {
            return
        }

        try {
            def data = cellData ?: loadCellData(server, pathObject, downsampleFactor, compartments)
            if (data == null) return

            def measurements = pathObject.getMeasurementList()
            int nChannels = server.nChannels()

            for (int c = 0; c < nChannels; c++) {
                def channelName = server.getChannel(c).getName()
                def pixelValues = data.allChannelPixels[c]

                data.masks.each { compartment, mask ->
                    if (compartments.contains(compartment)) {
                        def compartmentPixels = getCompartmentPixels(pixelValues, mask)
                        if (compartmentPixels.length > 0) {
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
     * Simplify a ROI using Douglas-Peucker algorithm
     * @param roi ROI to simplify
     * @param tolerance Douglas-Peucker tolerance in pixels
     * @return Simplified ROI
     */
    static ROI simplifyROI(ROI roi, double tolerance = 0.5) {
        if (roi == null) return null
        try {
            def geometry = roi.getGeometry()
            def simplified = DouglasPeuckerSimplifier.simplify(geometry, tolerance)
            return GeometryTools.geometryToROI(simplified, roi.getImagePlane())
        } catch (Exception e) {
            println "Warning: Could not simplify ROI, using original: ${e.message}"
            return roi
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
        try {
            def roiIJ = IJTools.convertToIJRoi(roi, pathImage)
            mask.setColor(255)
            mask.fill(roiIJ)
        } catch (UnsupportedOperationException e) {
            println "Warning: could not convert ROI to ImageJ ROI, returning empty mask: ${e.getMessage()}"
        }
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
     * Extract pixel values for a specific channel using a bulk raster read.
     */
    def extractChannelPixels(BufferedImage img, int channel, int width, int height) {
        return img.getRaster().getSamples(0, 0, width, height, channel, (float[]) null)
    }

    /**
     * Get pixel values within a compartment mask as a primitive double array, avoiding boxing overhead.
     */
    def getCompartmentPixels(float[] allPixels, ByteProcessor mask) {
        def maskPixels = mask.getPixels() as byte[]

        // Count non-zero mask pixels first to allocate exact array size
        int count = 0
        for (int i = 0; i < maskPixels.length; i++) {
            if (maskPixels[i] != 0) count++
        }

        double[] compartmentPixels = new double[count]
        int j = 0
        for (int i = 0; i < maskPixels.length; i++) {
            if (maskPixels[i] != 0) compartmentPixels[j++] = allPixels[i]
        }

        return compartmentPixels
    }

    /**
     * Calculate and add percentile measurements for a compartment
     */
    def addPercentileMeasurementsForCompartment(MeasurementList measurements, double[] pixels,
                                               String channelName, String compartment, List<Double> percentiles) {
        def stats = new DescriptiveStatistics(pixels)

        percentiles.each { percentile ->
            def value = stats.getPercentile(percentile)
            def compartmentName = compartment.toLowerCase().capitalize()
            def measurementName = "${channelName}: ${compartmentName}: Percentile: ${percentile}"
            measurements.put(measurementName, value)
        }
    }

    /**
     * Add erosion-based measurements
     * @param server ImageServer containing the pixel data
     * @param pathObject PathObject to measure (must be a cell)
     * @param downsampleFactor Resolution at which to request pixels
     * @param erosionSteps List of erosion distances in pixels (must be positive integers)
     * @param compartments Set of compartments to measure ('CELL', 'NUCLEUS')
     * @param cellData Optional pre-loaded cell data from loadCellData(); if null, data is loaded internally
     */
    def addErosionMeasurements(server, PathObject pathObject, double downsampleFactor,
                               List<Integer> erosionSteps, Set<String> compartments = ['CELL', 'NUCLEUS'],
                               Map cellData = null) {
        try {
            if (!pathObject.isCell()) return
            if (erosionSteps == null || erosionSteps.isEmpty()) return
            
            def cellROI = pathObject.getROI()
            def nucleusROI = pathObject.getNucleusROI()
            if (cellROI == null) return

            def data = cellData ?: loadCellData(server, pathObject, downsampleFactor, compartments)
            if (data == null) return

            def measurements = pathObject.getMeasurementList()

            // Process each compartment
            compartments.each { compartment ->
                def baseMask = (compartment == 'NUCLEUS' && nucleusROI != null)
                               ? data.masks['NUCLEUS']
                               : data.masks['CELL']
                if (baseMask == null) return
                
                // Use the first channel's pixel values to count mask area
                def basePixels = getCompartmentPixels(data.allChannelPixels[0], baseMask)
                int baseArea = basePixels.length
                
                if (baseArea == 0) return
                
                // Erode incrementally through sorted steps, reusing each prior erosion result.
                // e.g. for steps [4,7,11]: apply 4 dilations, then 3 more, then 4 more (17 total)
                // rather than 4+7+11=22 dilations if computed from baseMask each time.
                def sortedSteps = erosionSteps.toSorted()
                def currentMask = baseMask
                int prevSteps = 0

                sortedSteps.each { steps ->
                    int additionalSteps = steps - prevSteps
                    currentMask = erodeMask(currentMask, additionalSteps)
                    prevSteps = steps

                    def erodedPixels = getCompartmentPixels(data.allChannelPixels[0], currentMask)
                    int erodedArea = erodedPixels.length
                    
                    def compartmentName = compartment.toLowerCase().capitalize()
                    
                    // Add area fraction (only once per erosion step, not per channel)
                    def areaFraction = erodedArea / (double)baseArea
                    measurements.put(
                        "${compartmentName}: Eroded_${steps}px: Area_Fraction",
                        areaFraction
                    )
                    
                    // Add intensity measurements if pixels remain
                    if (erodedArea > 0) {
                        for (int c = 0; c < server.nChannels(); c++) {
                            def channelName = server.getChannel(c).getName()
                            def channelErodedPixels = getCompartmentPixels(data.allChannelPixels[c], currentMask)
                            
                            if (channelErodedPixels.length > 0) {
                                def stats = new DescriptiveStatistics(channelErodedPixels)
                                
                                measurements.put(
                                    "${channelName}: ${compartmentName}: Eroded_${steps}px: Mean",
                                    stats.getMean()
                                )
                                measurements.put(
                                    "${channelName}: ${compartmentName}: Eroded_${steps}px: Median",
                                    stats.getPercentile(50)
                                )
                            }
                        }
                    }
                }
            }
        } catch (Exception e) {
            println("Error in erosion measurements for ${pathObject}: ${e.getMessage()}")
        }
    }

    /**
     * Erode a binary mask by specified number of pixels
     * Note: ImageJ's erode() operates on foreground (255) pixels. In our masks, ROI pixels are 255.
     * However, erosion in ImageJ removes boundary pixels by treating them as background,
     * which actually expands the 255 region inward. We use dilate() instead which shrinks
     * the 255 region, giving us the desired erosion effect for our ROI measurements.
     */
    def erodeMask(ByteProcessor mask, int erosionSteps) {
        def eroded = mask.duplicate()
        for (int i = 0; i < erosionSteps; i++) {
            eroded.dilate()
        }
        return eroded
    }

    /**
     * Convert a Bio-Formats ImageServer to an ImageJ ImagePlus object to ensure compatibility with ROI conversion and feature generation code
     * @param server ImageServer to convert
     * @return ImagePlus containing the image data
     */
    static convertServerToImagePlus(server) {
        def request = RegionRequest.createInstance(server)
        return IJTools.convertToImagePlus(server, request)
    }

    /**
     * Load a mask image file using Bio-Formats and convert to ImagePlus
     * @param filePath Path to the mask image file
     * @return ImagePlus containing the mask image data
     */
    static loadMaskAsImagePlus(String filePath) {
        def uri = Paths.get(filePath).toUri()
        def builder = new BioFormatsServerBuilder()
        def server = builder.buildServer(uri)
        return convertServerToImagePlus(server).getImage()
    }

    @Override
    void run() {
        if (tileSize <= 0) {
            throw new IllegalArgumentException('tileSize must be > 0')
        }
        if (tileOverlap < 0) {
            throw new IllegalArgumentException('tileOverlap must be >= 0')
        }

        // Load whole cell mask using Bio-Formats via QuPath, more robust than ImageJ for very large images
        def wholeCellImp = loadMaskAsImagePlus(wholeCellMaskFilePath)
        println 'Loaded whole cell mask width: ' + wholeCellImp.getWidth()

        // Load nuclear mask using Bio-Formats via QuPath, more robust than ImageJ for very large images
        def nuclearImp = loadMaskAsImagePlus(nuclearMaskFilePath)
        println 'Loaded nuclear mask width: ' + nuclearImp.getWidth()

        // Build a server with supplied TIFF file using bioformats
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

        if (simplifyROIs && simplifyTolerance > 0) {
            println 'Simplifying ROIs with tolerance: ' + simplifyTolerance
            GParsPool.withPool(threads) {
                wholeCellROIs = wholeCellROIs.collectParallel { pathObject ->
                    simplifyROI(pathObject, simplifyTolerance)
                }
                nuclearROIs = nuclearROIs.collectParallel { pathObject ->
                    simplifyROI(pathObject, simplifyTolerance)
                }
            }
            println "Simplified ${wholeCellROIs.size()} whole cell ROIs"
            println "Simplified ${nuclearROIs.size()} nuclear ROIs"
        }

        // Convert QuPath ROIs to objects and add them to the hierarchy
        def pathObjects = makeCellObjects(
            wholeCellROIs, nuclearROIs, distThreshold, estimateCellBoundaryDist, threads
        )
        println 'Total path objects: ' + pathObjects.size()

        // Filter out any cells that have a membrane outside of the image bounds
        def imageWidth = server.getWidth()
        def imageHeight = server.getHeight()
        pathObjects = removeOutOfBoundsCells(pathObjects, imageWidth, imageHeight)

        if (!skipMeasurements) {
            // Set the pixel calibration
            PixelCalibration pixelCalibration = new PixelCalibration.Builder()
                .pixelSizeMicrons(pixelSizeMicrons, pixelSizeMicrons)
                .build()
            println 'Set pixel calibration: ' + pixelCalibration

            // Pre-parse optional measurement parameters before entering the parallel block
            def percentileList = percentiles ? percentiles.split(',').collect { it as Double } : null
            def erosionList = null
            if (erosionSteps && erosionSteps.trim()) {
                def parsedSteps = erosionSteps.split(',').collect { it.trim() as Integer }.findAll { it > 0 }
                if (parsedSteps.isEmpty()) {
                    println 'Warning: No valid positive erosion steps provided, skipping erosion measurements'
                } else {
                    erosionList = parsedSteps
                }
            }
            if (percentileList) 
                println 'Will add intensity percentiles: ' + percentileList
            if (erosionList)
                println 'Will add erosion measurements at steps: ' + erosionList

            def allCompartments      = ['CELL', 'NUCLEUS', 'CYTOPLASM', 'MEMBRANE'] as Set
            def erosionCompartments  = ['CELL', 'NUCLEUS'] as Set
            def percentileCompartments = ['NUCLEUS', 'CYTOPLASM', 'MEMBRANE', 'CELL'] as Set
            def ds = downsampleFactor as double

            // Group cells by tile (centroid determines tile assignment)
            def tileGroups = groupCellsByTile(pathObjects, tileSize, imageWidth, imageHeight)
            println "Adding cell measurements in ${tileGroups.size()} tiles (tile size: ${tileSize}px, overlap: ${tileOverlap}px)..."
            // Overlap ensures cells near tile edges have their full bounding box within the read region

            int processedCells = 0
            int totalCells = pathObjects.size()
            def tileReadFallbackTiles = new AtomicInteger(0)
            def tileFallbackCells = new AtomicInteger(0)
            def skippedCells = new AtomicInteger(0)

            tileGroups.each { tileKey, tileCells ->
                int tileCol = tileKey[0]
                int tileRow = tileKey[1]

                // Tile region in full image coordinates, padded by tileOverlap on all sides
                int txFull   = Math.max(0,           tileCol * tileSize - tileOverlap)
                int tyFull   = Math.max(0,           tileRow * tileSize - tileOverlap)
                int txEnd    = Math.min(imageWidth,  (tileCol + 1) * tileSize + tileOverlap)
                int tyEnd    = Math.min(imageHeight, (tileRow + 1) * tileSize + tileOverlap)
                int twFull   = txEnd - txFull
                int thFull   = tyEnd - tyFull

                // Read this tile once from the server
                def tileRequest = RegionRequest.createInstance(
                    server.getPath(), ds, txFull, tyFull, twFull, thFull
                )
                def tileImg = server.readRegion(tileRequest)
                if (tileImg == null) {
                    println "Warning: tile read returned null for tile (${tileCol}, ${tileRow}); falling back to per-cell reads for ${tileCells.size()} cells"
                    tileReadFallbackTiles.incrementAndGet()
                    tileFallbackCells.addAndGet(tileCells.size())

                    GParsPool.withPool(threads) {
                        tileCells.eachParallel { pathObject ->
                            ObjectMeasurements.addShapeMeasurements(
                                pathObject,
                                pixelCalibration,
                                ObjectMeasurements.ShapeFeatures.AREA,
                                ObjectMeasurements.ShapeFeatures.CIRCULARITY,
                                ObjectMeasurements.ShapeFeatures.LENGTH,
                                ObjectMeasurements.ShapeFeatures.MAX_DIAMETER,
                                ObjectMeasurements.ShapeFeatures.MIN_DIAMETER,
                                ObjectMeasurements.ShapeFeatures.NUCLEUS_CELL_RATIO,
                                ObjectMeasurements.ShapeFeatures.SOLIDITY
                            )

                            def cellData = loadCellData(server, pathObject, ds, allCompartments)
                            if (cellData == null) {
                                skippedCells.incrementAndGet()
                                return
                            }

                            addStandardIntensityMeasurements(server, pathObject, cellData)
                            if (percentileList) {
                                addPercentileMeasurements(
                                    server, pathObject, ds, percentileList,
                                    percentileCompartments, cellData
                                )
                            }
                            if (erosionList) {
                                addErosionMeasurements(
                                    server, pathObject, ds, erosionList,
                                    erosionCompartments, cellData
                                )
                            }
                        }
                    }

                    processedCells += tileCells.size()
                    println "Progress: ${processedCells}/${totalCells} cells (${(processedCells * 100.0 / totalCells).toInteger()}%)"
                    return
                }

                int tilePxW = tileImg.getWidth()
                int tilePxH = tileImg.getHeight()
                int tilePxX = (int)(txFull / ds)
                int tilePxY = (int)(tyFull / ds)

                // Extract all channel pixel arrays from the tile once, shared across all cells
                def nChannels = server.nChannels()
                def tileChannelPixels = (0..<nChannels).collect { c ->
                    extractChannelPixels(tileImg, c, tilePxW, tilePxH)
                }

                GParsPool.withPool(threads) {
                    tileCells.eachParallel { pathObject ->
                        // Shape measurements (geometry only, no image I/O)
                        ObjectMeasurements.addShapeMeasurements(
                            pathObject,
                            pixelCalibration,
                            ObjectMeasurements.ShapeFeatures.AREA,
                            ObjectMeasurements.ShapeFeatures.CIRCULARITY,
                            ObjectMeasurements.ShapeFeatures.LENGTH,
                            ObjectMeasurements.ShapeFeatures.MAX_DIAMETER,
                            ObjectMeasurements.ShapeFeatures.MIN_DIAMETER,
                            ObjectMeasurements.ShapeFeatures.NUCLEUS_CELL_RATIO,
                            ObjectMeasurements.ShapeFeatures.SOLIDITY
                        )

                        // Load all cell data from the shared tile pixels — no per-cell I/O
                        def cellData = loadCellDataFromTile(
                            pathObject, ds, allCompartments,
                            tileChannelPixels, tilePxX, tilePxY, tilePxW, tilePxH
                        )
                        if (cellData == null) {
                            tileFallbackCells.incrementAndGet()
                            cellData = loadCellData(server, pathObject, ds, allCompartments)
                        }

                        if (cellData == null) {
                            skippedCells.incrementAndGet()
                            return
                        }

                        addStandardIntensityMeasurements(server, pathObject, cellData)
                        if (percentileList) {
                            addPercentileMeasurements(
                                server, pathObject, ds, percentileList,
                                percentileCompartments, cellData
                            )
                        }
                        if (erosionList) {
                            addErosionMeasurements(
                                server, pathObject, ds, erosionList,
                                erosionCompartments, cellData
                            )
                        }
                    }
                }

                processedCells += tileCells.size()
                println "Progress: ${processedCells}/${totalCells} cells (${(processedCells * 100.0 / totalCells).toInteger()}%)"
            }

            if (tileReadFallbackTiles.get() > 0 || tileFallbackCells.get() > 0 || skippedCells.get() > 0) {
                println "Tile fallback summary: null-tile-fallback-tiles=${tileReadFallbackTiles.get()}, per-cell-fallbacks=${tileFallbackCells.get()}, skipped-cells=${skippedCells.get()}"
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
