package wehisodahub.cellmeasurement

import spock.lang.Specification
import spock.lang.TempDir
import spock.lang.Shared
import spock.lang.Subject
import spock.lang.Unroll

import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.Files
import java.awt.geom.Point2D
import java.net.URI
import java.awt.image.BufferedImage

import ij.ImagePlus
import ij.process.ByteProcessor

import qupath.imagej.tools.IJTools

import qupath.lib.roi.interfaces.ROI
import qupath.lib.roi.ROIs
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane
import qupath.lib.regions.RegionRequest
import qupath.lib.images.PathImage
import qupath.lib.images.servers.ImageChannel
import qupath.lib.images.servers.ImageServer
import qupath.lib.images.servers.bioformats.BioFormatsServerBuilder
import qupath.lib.measurements.MeasurementList

class AppSpec extends Specification {

    @TempDir
    Path tempDir

    @Shared
    @Subject
    App app = new App()

    def "should parse command line arguments correctly"() {
        given:
        def nuclearMask = createTempFile("nuclear.tiff")
        def wholeCellMask = createTempFile("wholecell.tiff")
        def tiffFile = createTempFile("image.tiff")
        def outputFile = tempDir.resolve("output.geojson").toString()

        when:
        app.nuclearMaskFilePath = nuclearMask
        app.wholeCellMaskFilePath = wholeCellMask
        app.tiffFilePath = tiffFile
        app.outputFilePath = outputFile
        app.downsampleFactor = 2.0
        app.pixelSizeMicrons = 1.0
        app.distThreshold = 15.0

        then:
        app.nuclearMaskFilePath == nuclearMask
        app.wholeCellMaskFilePath == wholeCellMask
        app.tiffFilePath == tiffFile
        app.outputFilePath == outputFile
        app.downsampleFactor == 2.0
        app.pixelSizeMicrons == 1.0
        app.distThreshold == 15.0
    }

    def "should have correct default values"() {
        given:
        def freshApp = new App()

        expect:
        freshApp.downsampleFactor == 1.0
        freshApp.pixelSizeMicrons == 0.5
        freshApp.skipMeasurements == false
        freshApp.distThreshold == 10.0
        freshApp.estimateCellBoundaryDist == 3.0
        freshApp.threads == 1
        freshApp.tileSize == 2048
        freshApp.tileOverlap == 200
    }

    def "extractROIs should return empty list for image with no objects"() {
        given:
        def processor = new ByteProcessor(100, 100)
        def image = new ImagePlus("test", processor)

        when:
        def result = App.extractROIs(image, 1.0, 1)

        then:
        result.isEmpty()
    }

    def "should return null when ROI is null"() {
        when:
        def result = App.simplifyROI(null, 0.5)

        then:
        result == null
    }

    def "should reduce point count with higher tolerance"() {
        given: "a complex polygon with many points"
        def numPoints = 100
        def x = (0..numPoints).collect { it * 10 + Math.sin(it * 0.1) * 2 } as double[]
        def y = (0..numPoints).collect { 50 + Math.cos(it * 0.1) * 2 } as double[]
        def roi = ROIs.createPolygonROI(x, y, ImagePlane.getDefaultPlane())
        def originalPoints = roi.getNumPoints()

        when:
        def simplifiedLow = App.simplifyROI(roi, 0.5)
        def simplifiedMid = App.simplifyROI(roi, 2.0)
        def simplifiedHigh = App.simplifyROI(roi, 5.0)

        then: "higher tolerance produces fewer points"
        simplifiedLow.getNumPoints() >= simplifiedMid.getNumPoints()
        simplifiedMid.getNumPoints() >= simplifiedHigh.getNumPoints()
        simplifiedHigh.getNumPoints() < originalPoints
        simplifiedHigh.getNumPoints() > 2 // Should still retain shape
    }

    def "should preserve ImagePlane information"() {
        given: "a polygon on a specific z and t plane"
        def plane = ImagePlane.getPlane(2, 5)
        def x = [10, 100, 100, 10] as double[]
        def y = [10, 10, 100, 100] as double[]
        def roi = ROIs.createPolygonROI(x, y, plane)

        when:
        def simplified = App.simplifyROI(roi, 0.5)

        then:
        simplified.getImagePlane() == plane
        simplified.getZ() == 2
        simplified.getT() == 5
    }

    def "should preserve area within tolerance for #roiType"() {
        given:
        def roi = createROI(roiType)
        def originalArea = roi.getArea()

        when:
        def simplified = App.simplifyROI(roi, tolerance)
        def simplifiedArea = simplified.getArea()
        def areaChange = Math.abs(simplifiedArea - originalArea) / originalArea * 100

        then:
        simplified != null
        areaChange < maxAreaChange

        where:
        roiType     | tolerance | maxAreaChange
        "polygon"   | 0.5       | 1.0
        "ellipse"   | 1.0       | 5.0
        "rectangle" | 0.5       | 0.1
    }

    def "should handle different ROI types without error"() {
        when:
        def simplified = App.simplifyROI(roi, 1.0)

        then:
        simplified != null
        noExceptionThrown()

        where:
        roi << [
            ROIs.createPolygonROI([10, 100, 100, 10] as double[], [10, 10, 100, 100] as double[], ImagePlane.getDefaultPlane()),
            ROIs.createEllipseROI(100, 100, 50, 30, ImagePlane.getDefaultPlane()),
            ROIs.createRectangleROI(10, 10, 100, 50, ImagePlane.getDefaultPlane()),
            ROIs.createPolylineROI([10, 50, 100] as double[], [10, 30, 25] as double[], ImagePlane.getDefaultPlane()),
            ROIs.createPointsROI([50] as double[], [50] as double[], ImagePlane.getDefaultPlane())
        ]
    }

    def "findNearestROI should return null when no ROIs within threshold"() {
        given:
        def centroid = new Point2D.Double(50, 50)
        def roi1 = createMockROI(10, 10)
        def roi2 = createMockROI(90, 90)
        def rois = [roi1, roi2]
        def threshold = 5.0

        when:
        def result = App.findNearestROI(centroid, rois, threshold)

        then:
        result == null
    }

    def "findNearestROI should return nearest ROI within threshold"() {
        given:
        def centroid = new Point2D.Double(50, 50)
        def roi1 = createMockROI(52, 48)  // distance ~2.8
        def roi2 = createMockROI(45, 55)  // distance ~7.1
        def rois = [roi1, roi2]
        def threshold = 10.0

        when:
        def result = App.findNearestROI(centroid, rois, threshold)

        then:
        result == roi1
    }

    def "removeOutOfBoundsCells should filter cells outside image bounds"() {
        given:
        def inBoundsCell = createMockPathObject(10, 10, 20, 20)
        def outOfBoundsCell = createMockPathObject(95, 95, 20, 20)  // extends beyond 100x100
        def pathObjects = [inBoundsCell, outOfBoundsCell]

        when:
        def result = App.removeOutOfBoundsCells(pathObjects, 100, 100)

        then:
        result.size() == 1
        result[0] == inBoundsCell
    }

    def "removeOutOfBoundsCells should keep all cells when within bounds"() {
        given:
        def cell1 = createMockPathObject(10, 10, 20, 20)
        def cell2 = createMockPathObject(50, 50, 30, 30)
        def pathObjects = [cell1, cell2]

        when:
        def result = App.removeOutOfBoundsCells(pathObjects, 100, 100)

        then:
        result.size() == 2
        result.containsAll([cell1, cell2])
    }

    def "groupCellsByTile should group path objects by centroid tile"() {
        given:
        def cellA = createMockCellPathObjectWithCentroid(100, 100)
        def cellB = createMockCellPathObjectWithCentroid(1999, 1999)
        def cellC = createMockCellPathObjectWithCentroid(2100, 100)
        def tileSize = 2048

        when:
        def grouped = App.groupCellsByTile([cellA, cellB, cellC], tileSize, 5000, 5000)

        then:
        grouped.size() == 2
        grouped[[0, 0]].containsAll([cellA, cellB])
        grouped[[1, 0]].contains(cellC)
    }

    def "extractSubRegionPixels should copy rectangular subregion in row-major order"() {
        given:
        // 4x4 tile:
        //  1  2  3  4
        //  5  6  7  8
        //  9 10 11 12
        // 13 14 15 16
        def tilePixels = [
            1f, 2f, 3f, 4f,
            5f, 6f, 7f, 8f,
            9f, 10f, 11f, 12f,
            13f, 14f, 15f, 16f
        ] as float[]

        when:
        // 2x2 starting at (1,1) => [6,7,10,11]
        def sub = App.extractSubRegionPixels(tilePixels, 4, 1, 1, 2, 2)

        then:
        sub as List == [6f, 7f, 10f, 11f]
    }

    def "loadCellDataFromTile should return masks and cropped pixels when cell fits tile"() {
        given:
        def pathObject = createRealCellPathObject(2, 3, 2, 2, 2, 3, 1, 1)
        def tilePixels = new float[10 * 10]
        for (int i = 0; i < tilePixels.length; i++) tilePixels[i] = i as float

        when:
        def data = app.loadCellDataFromTile(
            pathObject, 1.0, ['CELL', 'NUCLEUS'] as Set,
            [tilePixels], 0, 0, 10, 10
        )

        then:
        data != null
        data.masks.keySet().containsAll(['CELL', 'NUCLEUS'])
        data.allChannelPixels.size() == 1
        // width/height include +1 in App.getCellBoundingBox logic: (2x2 ROI -> 3x3 pixels)
        data.allChannelPixels[0].length == 9
    }

    def "loadCellDataFromTile should return null when cell bounds fall outside tile"() {
        given:
        def pathObject = createRealCellPathObject(2, 3, 2, 2, 2, 3, 1, 1)
        def tilePixels = new float[10 * 10]

        when:
        def data = app.loadCellDataFromTile(
            pathObject, 1.0, ['CELL'] as Set,
            [tilePixels], 5, 5, 10, 10
        )

        then:
        data == null
    }

    def "makeCellObjects should handle empty ROI lists"() {
        when:
        def result = App.makeCellObjects([], [], 10.0, 3.0)

        then:
        result.isEmpty()
    }

    def "makeCellObjects should not filter valid non-overlapping cells"() {
        given:
        def plane = ImagePlane.getDefaultPlane()
        // Two well-separated cells — constrainCellOverlaps won't clip either, so both should survive the filter
        def cell1 = ROIs.createRectangleROI(0, 0, 20, 20, plane)
        def nuc1  = ROIs.createRectangleROI(5, 5, 10, 10, plane)
        def cell2 = ROIs.createRectangleROI(200, 200, 20, 20, plane)
        def nuc2  = ROIs.createRectangleROI(205, 205, 10, 10, plane)

        when:
        def result = App.makeCellObjects([cell1, cell2], [nuc1, nuc2], 20.0, 3.0)

        then:
        result.size() == 2
        result.every { cell ->
            def geom = cell.getROI().getGeometry()
            geom.getArea() > 0 && (geom.getGeometryType() == 'Polygon' || geom.getGeometryType() == 'MultiPolygon')
        }
    }

    def "matchROIs should return pairs for each nuclear ROI"() {
        given:
        def nuclearROI = createMockROI(50, 50)
        def wholeCellROI = createMockROI(52, 48)
        def nuclearROIs = [nuclearROI]
        def wholeCellROIs = [wholeCellROI]

        when:
        def result = App.matchROIs(nuclearROIs, wholeCellROIs, 10.0, 3.0, 1)

        then:
        result.size() == 1
        result[0].size() == 2
        result[0][0] == nuclearROI
    }
 
    def "addPercentileMeasurements should return early when pathObject is not a cell"() {
        given:
        def server = createMockImageServer(100, 100)
        def pathObject = createMockPathObject(10, 10, 20, 20)
        pathObject.isCell() >> false

        when:
        app.addPercentileMeasurements(server, pathObject)

        then:
        0 * _
    }

    def "addPercentileMeasurements should return early when cellROI is null"() {
        given:
        def server = createMockImageServer(100, 100)
        def pathObject = createMockPathObject(10, 10, 20, 20)
        pathObject.isCell() >> true
        pathObject.getROI() >> null

        when:
        app.addPercentileMeasurements(server, pathObject)

        then:
        0 * server._
    }

    @Unroll
    def "addPercentileMeasurements should process cell with #compartments compartments"() {
        given:
        def server = createMockImageServer(100, 100)
        def pathObject = createMockPathObject(10, 10, 50, 50)
        pathObject.isCell() >> true
        
        when:
        app.addPercentileMeasurements(server, pathObject, 1.0, [95.0], compartments as Set)

        then:
        noExceptionThrown()

        where:
        compartments << [['CELL'], ['NUCLEUS'], ['CYTOPLASM'], ['MEMBRANE']]
    }

    def "createCompartmentMasks should create only requested compartment masks"() {
        given:
        def cellROI = ROIs.createRectangleROI(0, 0, 50, 50, ImagePlane.getDefaultPlane())
        def nucleusROI = ROIs.createRectangleROI(10, 10, 20, 20, ImagePlane.getDefaultPlane())
        def tiffFilePath = getClass().getResource('/synthetic_test.ome.tif').toURI()
        def pathImage = getPathImage(tiffFilePath)
        def compartments = ['CELL', 'NUCLEUS'] as Set

        when:
        def masks = app.createCompartmentMasks(cellROI, nucleusROI, 100, 100, pathImage, compartments)

        then:
        masks.keySet() == ['CELL', 'NUCLEUS'] as Set
        masks['CELL'] instanceof ByteProcessor
        masks['NUCLEUS'] instanceof ByteProcessor
    }

    def "createCompartmentMasks should create cytoplasm mask when both cell and nucleus are requested"() {
        given:
        def cellROI = ROIs.createRectangleROI(0, 0, 50, 50, ImagePlane.getDefaultPlane())
        def nucleusROI = ROIs.createRectangleROI(10, 10, 20, 20, ImagePlane.getDefaultPlane())
        def tiffFilePath = getClass().getResource('/synthetic_test.ome.tif').toURI()
        def pathImage = getPathImage(tiffFilePath)
        def compartments = ['CELL', 'NUCLEUS', 'CYTOPLASM'] as Set

        when:
        def masks = app.createCompartmentMasks(cellROI, nucleusROI, 100, 100, pathImage, compartments)

        then:
        masks.keySet().containsAll(['CELL', 'NUCLEUS', 'CYTOPLASM'])
        masks['CYTOPLASM'] instanceof ByteProcessor
    }

    def "createROIMask should create ByteProcessor with correct dimensions"() {
        given:
        def roi = ROIs.createRectangleROI(10, 10, 20, 20, ImagePlane.getDefaultPlane())
        def tiffFilePath = getClass().getResource('/synthetic_test.ome.tif').toURI()
        def pathImage = getPathImage(tiffFilePath)

        when:
        def mask = app.createROIMask(roi, 50, 30, pathImage)

        then:
        mask instanceof ByteProcessor
        mask.getWidth() == 50
        mask.getHeight() == 30
    }

    def "createROIMask should return empty mask for unsupported ROI geometry"() {
        given:
        // PointsROI triggers UnsupportedOperationException in the pathImage overload of
        // IJTools.convertToIJRoi. The fix should catch this and return an empty mask.
        def roi = ROIs.createPointsROI(50.0, 50.0, ImagePlane.getDefaultPlane())
        def tiffFilePath = getClass().getResource('/synthetic_test.ome.tif').toURI()
        def pathImage = getPathImage(tiffFilePath)

        when:
        def mask = app.createROIMask(roi, 50, 30, pathImage)

        then:
        noExceptionThrown()
        mask instanceof ByteProcessor
        mask.getWidth() == 50
        mask.getHeight() == 30
        // mask should be all zeros (no pixels filled) because the ROI could not be converted
        countNonZeroPixels(mask) == 0
    }

    def "createROIMaskFromOrigin should return a valid mask for any ROI geometry type"() {
        given:
        // A LineROI has a non-polygon geometry — verify the method is resilient and always
        // returns a valid ByteProcessor. The UnsupportedOperationException catch is a safety
        // net for degenerate geometries produced by constrainCellOverlaps in production;
        // those geometries are difficult to reproduce without real Voronoi-clipped cells.
        def roi = ROIs.createLineROI(0.0, 0.0, 10.0, 0.0, ImagePlane.getDefaultPlane())

        when:
        def mask = app.createROIMaskFromOrigin(roi, 100, 100, 0.0, 0.0, 1.0)

        then:
        noExceptionThrown()
        mask instanceof ByteProcessor
        mask.getWidth() == 100
        mask.getHeight() == 100
    }

    def "getCompartmentPixels should extract only masked pixels"() {
        given:
        def allPixels = [10.0f, 20.0f, 30.0f, 40.0f] as float[]
        def mask = new ByteProcessor(2, 2)
        mask.set(0, 0, 255)  // Include first pixel
        mask.set(1, 1, 255)  // Include fourth pixel

        when:
        def compartmentPixels = app.getCompartmentPixels(allPixels, mask)

        then:
        compartmentPixels.length == 2
        10.0d in compartmentPixels
        40.0d in compartmentPixels
    }

    @Unroll
    def "extractChannelPixels should extract pixels for channel #channel"() {
        given:
        def img = createMockBufferedImage(2, 2, channel, expectedValue)

        when:
        def pixels = app.extractChannelPixels(img, channel, 2, 2)

        then:
        pixels.length == 4
        pixels.every { it == expectedValue }

        where:
        channel | expectedValue
        0       | 100.0f
        1       | 200.0f
        2       | 50.0f
    }

    @Unroll
    def "addPercentileMeasurementsForCompartment should add measurements for percentiles #percentiles"() {
        given:
        def measurements = Mock(MeasurementList)
        def pixels = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0] as double[]
        def channelName = "TestChannel"
        def compartment = "NUCLEUS"

        when:
        app.addPercentileMeasurementsForCompartment(measurements, pixels, channelName, compartment, percentiles)

        then:
        percentiles.size() * measurements.put(_, _)

        where:
        percentiles << [[50.0], [70.0, 90.0], [95.0, 99.0]]
    }

    def "addErosionMeasurements should return early when pathObject is not a cell"() {
        given:
        def server = createMockImageServer(100, 100)
        def pathObject = createMockPathObject(10, 10, 20, 20)
        pathObject.isCell() >> false

        when:
        app.addErosionMeasurements(server, pathObject, 1.0, [4, 7])

        then:
        0 * _
    }

    def "addErosionMeasurements should return early when erosionSteps is null"() {
        given:
        def server = createMockImageServer(100, 100)
        def pathObject = createMockPathObject(10, 10, 20, 20)
        pathObject.isCell() >> true

        when:
        app.addErosionMeasurements(server, pathObject, 1.0, null)

        then:
        0 * server._
    }

    def "addErosionMeasurements should return early when erosionSteps is empty"() {
        given:
        def server = createMockImageServer(100, 100)
        def pathObject = createMockPathObject(10, 10, 20, 20)
        pathObject.isCell() >> true

        when:
        app.addErosionMeasurements(server, pathObject, 1.0, [])

        then:
        0 * server._
    }

    def "addErosionMeasurements should return early when cellROI is null"() {
        given:
        def server = createMockImageServer(100, 100)
        def pathObject = createMockPathObject(10, 10, 20, 20)
        pathObject.isCell() >> true
        pathObject.getROI() >> null

        when:
        app.addErosionMeasurements(server, pathObject, 1.0, [4])

        then:
        0 * server._
    }

    @Unroll
    def "addErosionMeasurements should write area fraction measurements for #compartments compartments"() {
        given:
        def server = createMockImageServerWithChannels(100, 100, 1, ['DAPI'])
        def pathObject = createRealCellPathObject(0, 0, 50, 50, 10, 10, 20, 20)
        def cellData = createTestCellData(50, 50, 1)

        when:
        app.addErosionMeasurements(server, pathObject, 1.0, [4], compartments as Set, cellData)

        then:
        pathObject.getMeasurementList().containsKey(expectedKey)

        where:
        compartments           | expectedKey
        ['CELL']               | 'Cell: Eroded_4px: Area_Fraction'
        ['NUCLEUS']            | 'Nucleus: Eroded_4px: Area_Fraction'
        ['CELL', 'NUCLEUS']    | 'Cell: Eroded_4px: Area_Fraction'
    }

    @Unroll
    def "addErosionMeasurements should write channel intensity measurements"() {
        given:
        def server = createMockImageServerWithChannels(100, 100, 1, ['DAPI'])
        def pathObject = createRealCellPathObject(0, 0, 50, 50, 10, 10, 20, 20)
        def cellData = createTestCellData(50, 50, 1)

        when:
        app.addErosionMeasurements(server, pathObject, 1.0, [4], ['CELL'] as Set, cellData)

        then:
        def ml = pathObject.getMeasurementList()
        ml.containsKey('DAPI: Cell: Eroded_4px: Mean')
        ml.containsKey('DAPI: Cell: Eroded_4px: Median')
        !Double.isNaN(ml.get('DAPI: Cell: Eroded_4px: Mean') as double)
    }

    @Unroll
    def "addErosionMeasurements should write measurements for all erosion steps #steps"() {
        given:
        def server = createMockImageServerWithChannels(100, 100, 1, ['DAPI'])
        def pathObject = createRealCellPathObject(0, 0, 50, 50, 10, 10, 20, 20)
        def cellData = createTestCellData(50, 50, 1)

        when:
        app.addErosionMeasurements(server, pathObject, 1.0, steps, ['CELL'] as Set, cellData)

        then:
        def ml = pathObject.getMeasurementList()
        steps.every { step -> ml.containsKey("Cell: Eroded_${step}px: Area_Fraction") }

        where:
        steps << [[4], [4, 7], [4, 7, 11, 14]]
    }

    def "erodeMask should process mask without errors"() {
        given:
        def mask = new ByteProcessor(20, 20)
        // Fill center area with white
        for (int y = 5; y < 15; y++) {
            for (int x = 5; x < 15; x++) {
                mask.set(x, y, 255)
            }
        }
        def originalCount = countNonZeroPixels(mask)

        when:
        def erodedMask = app.erodeMask(mask, 2)
        def erodedCount = countNonZeroPixels(erodedMask)

        then:
        erodedMask != null
        erodedMask instanceof ByteProcessor
        erodedMask.getWidth() == 20
        erodedMask.getHeight() == 20
        // Dilate produces shrinking effect (erosion) in our masks
        erodedCount < originalCount
        erodedCount > 0
    }

    def "erodeMask with zero steps should return identical mask"() {
        given:
        def mask = new ByteProcessor(10, 10)
        mask.set(5, 5, 255)
        
        def originalCount = countNonZeroPixels(mask)

        when:
        def erodedMask = app.erodeMask(mask, 0)
        def erodedCount = countNonZeroPixels(erodedMask)

        then:
        erodedCount == originalCount
    }

    def "erodeMask with large steps should eliminate small masks"() {
        given:
        def mask = new ByteProcessor(10, 10)
        // Small 3x3 area
        for (int y = 4; y < 7; y++) {
            for (int x = 4; x < 7; x++) {
                mask.set(x, y, 255)
            }
        }

        when:
        def erodedMask = app.erodeMask(mask, 10)
        def erodedCount = countNonZeroPixels(erodedMask)

        then:
        // After troubleshooting, large erosion steps should eliminate small regions
        erodedCount == 0
    }

    def "getCellBoundingBox should return scaled bounds for valid ROI"() {
        given:
        def roi = ROIs.createRectangleROI(100, 200, 50, 60, ImagePlane.getDefaultPlane())
        def server = createMockImageServer(1000, 1000)
        def downsampleFactor = 2.0

        when:
        def bounds = App.getCellBoundingBox(roi, downsampleFactor, server)

        then:
        bounds.getBoundsX() == 50.0  // 100 / 2
        bounds.getBoundsY() == 100.0  // 200 / 2
        bounds.getBoundsWidth() == 26.0  // 50 / 2 + 1
        bounds.getBoundsHeight() == 31.0  // 60 / 2 + 1
    }

    def "getCellBoundingBox should return full image bounds for infinite ROI"() {
        given:
        def roi = Mock(ROI)
        roi.getBoundsX() >> Double.POSITIVE_INFINITY
        def server = createMockImageServer(500, 400)
        def downsampleFactor = 1.0

        when:
        def bounds = App.getCellBoundingBox(roi, downsampleFactor, server)

        then:
        bounds.getBoundsX() == 0.0
        bounds.getBoundsY() == 0.0
        bounds.getBoundsWidth() == 500.0
        bounds.getBoundsHeight() == 400.0
    }

    private PathObject createRealCellPathObject(double x, double y, double w, double h,
                                                double nx, double ny, double nw, double nh) {
        def plane = ImagePlane.getDefaultPlane()
        def cellROI = ROIs.createRectangleROI(x, y, w, h, plane)
        def nucleusROI = ROIs.createRectangleROI(nx, ny, nw, nh, plane)
        return PathObjects.createCellObject(cellROI, nucleusROI, null, null)
    }

    private Map createTestCellData(int width, int height, int nChannels) {
        def cellMask = new ByteProcessor(width, height)
        for (int y = 5; y < height - 5; y++) {
            for (int x = 5; x < width - 5; x++) {
                cellMask.set(x, y, 255)
            }
        }
        def nuclearMask = new ByteProcessor(width, height)
        for (int y = 10; y < height - 10; y++) {
            for (int x = 10; x < width - 10; x++) {
                nuclearMask.set(x, y, 255)
            }
        }
        def allChannelPixels = (0..<nChannels).collect { c ->
            def pixels = new float[width * height]
            for (int i = 0; i < pixels.length; i++) pixels[i] = (c + 1) * 100.0f
            pixels
        }
        return [masks: ['CELL': cellMask, 'NUCLEUS': nuclearMask], allChannelPixels: allChannelPixels]
    }

    private ImageServer createMockImageServerWithChannels(int width, int height, int nChannels,
                                                          List<String> channelNames) {
        def server = Mock(ImageServer)
        server.getWidth() >> width
        server.getHeight() >> height
        server.nChannels() >> nChannels
        channelNames.eachWithIndex { name, i ->
            def ch = ImageChannel.getInstance(name, null)
            server.getChannel(i) >> ch
        }
        return server
    }

    // Helper method for erosion tests
    private int countNonZeroPixels(ByteProcessor mask) {
        def pixels = mask.getPixels() as byte[]
        int count = 0
        for (int i = 0; i < pixels.length; i++) {
            if (pixels[i] != 0) count++
        }
        return count
    }

    def "should process real image files"() {
        given:
        app.downsampleFactor = 1.0
        app.pixelSizeMicrons = 0.5
        app.skipMeasurements = false
        app.distThreshold = 10.0
        app.estimateCellBoundaryDist = 3.0
        app.percentiles = '90,95'
        app.threads = 1

        def nuclearMaskPath = getClass().getResource('/synthetic_test_nuclear.tiff').toURI()
        def wholeCellMaskPath = getClass().getResource('/synthetic_test_whole-cell.tiff').toURI()
        def tiffFilePath = getClass().getResource('/synthetic_test.ome.tif').toURI()

        def nuclearMask = Paths.get(nuclearMaskPath).toString()
        def wholeCellMask = Paths.get(wholeCellMaskPath).toString()
        def tiffFile = Paths.get(tiffFilePath).toString()

        app.nuclearMaskFilePath = nuclearMask
        app.wholeCellMaskFilePath = wholeCellMask
        app.tiffFilePath = tiffFile
        app.outputFilePath = tempDir.resolve("synthetic_test_segmentation.geojson").toString()

        when:
        app.run()

        then:
        Files.exists(Path.of(app.outputFilePath))
    }

    // Helper methods
    private String createTempFile(String filename) {
        def file = tempDir.resolve(filename)
        Files.createFile(file)
        return file.toString()
    }

    private ROI createROI(String type) {
        switch(type) {
            case "polygon":
                def x = [0, 100, 100, 0] as double[]
                def y = [0, 0, 100, 100] as double[]
                return ROIs.createPolygonROI(x, y, ImagePlane.getDefaultPlane())
            case "ellipse":
                return ROIs.createEllipseROI(100, 100, 50, 30, ImagePlane.getDefaultPlane())
            case "rectangle":
                return ROIs.createRectangleROI(10, 10, 100, 50, ImagePlane.getDefaultPlane())
            default:
                throw new IllegalArgumentException("Unknown ROI type: ${type}")
        }
    }

    private ROI createMockROI(double centroidX, double centroidY) {
        def roi = Mock(ROI)
        roi.getCentroidX() >> centroidX
        roi.getCentroidY() >> centroidY
        roi.getImagePlane() >> ImagePlane.getDefaultPlane()
        return roi
    }

    private PathObject createMockCellPathObjectWithCentroid(double centroidX, double centroidY) {
        // Use a real PathCellObject to avoid mocking limitations on PathObject
        def plane = ImagePlane.getDefaultPlane()
        def cellROI = ROIs.createRectangleROI(centroidX - 1, centroidY - 1, 2, 2, plane)
        def nucleusROI = ROIs.createRectangleROI(centroidX - 0.5, centroidY - 0.5, 1, 1, plane)
        return PathObjects.createCellObject(cellROI, nucleusROI, null, null)
    }

    private PathObject createMockPathObject(double x, double y, double width, double height) {
        def roi = ROIs.createRectangleROI(x, y, width, height, ImagePlane.getDefaultPlane())
        def pathObject =  PathObjects.createDetectionObject(roi)
        pathObject.getROI() >> roi
        return pathObject
    }

    private ImageServer createMockImageServer(int width, int height) {
        def server = Mock(ImageServer)
        server.getWidth() >> width
        server.getHeight() >> height
        return server
    }

    def createMockBufferedImage(int width, int height, int targetChannel, float expectedValue) {
        def image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB)
        def raster = image.getRaster()

        // Fill the target channel with expected values
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                raster.setSample(x, y, targetChannel, expectedValue)
            }
        }

        return image
    }

    private PathImage getPathImage(URI imageUri) {
        def builder = new BioFormatsServerBuilder()
        def server = builder.buildServer(imageUri)

        def request = RegionRequest.createInstance(
            server.getPath(),
            1.0, 0, 0, server.getWidth(), server.getHeight()
        )

        return IJTools.convertToImagePlus(server, request);
    }
}
