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

    def "makeCellObjects should handle empty ROI lists"() {
        when:
        def result = App.makeCellObjects([], [], 10.0, 3.0)

        then:
        result.isEmpty()
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

    def "getCompartmentPixels should extract only masked pixels"() {
        given:
        def allPixels = [10.0f, 20.0f, 30.0f, 40.0f] as float[]
        def mask = new ByteProcessor(2, 2)
        mask.set(0, 0, 255)  // Include first pixel
        mask.set(1, 1, 255)  // Include fourth pixel

        when:
        def compartmentPixels = app.getCompartmentPixels(allPixels, mask)

        then:
        compartmentPixels.size() == 2
        compartmentPixels.contains(10.0d)
        compartmentPixels.contains(40.0d)
    }

    @Unroll
    def "addPercentileMeasurementsForCompartment should add measurements for percentiles #percentiles"() {
        given:
        def measurements = Mock(MeasurementList)
        def pixels = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        def channelName = "TestChannel"
        def compartment = "NUCLEUS"

        when:
        app.addPercentileMeasurementsForCompartment(measurements, pixels, channelName, compartment, percentiles)

        then:
        percentiles.size() * measurements.putMeasurement(_, _)

        where:
        percentiles << [[50.0], [70.0, 90.0], [95.0, 99.0]]
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

    private ROI createMockROI(double centroidX, double centroidY) {
        def roi = Mock(ROI)
        roi.getCentroidX() >> centroidX
        roi.getCentroidY() >> centroidY
        roi.getImagePlane() >> ImagePlane.getDefaultPlane()
        return roi
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
