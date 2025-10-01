package wehisodahub.cellmeasurement

import spock.lang.Specification
import spock.lang.TempDir
import spock.lang.Shared
import spock.lang.Subject

import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.Files
import java.awt.geom.Point2D

import ij.ImagePlus
import ij.process.ByteProcessor

import qupath.lib.roi.interfaces.ROI
import qupath.lib.roi.ROIs
import qupath.lib.objects.PathObject
import qupath.lib.objects.PathObjects
import qupath.lib.regions.ImagePlane


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

    def "should process real image files"() {
        given:
        app.downsampleFactor = 1.0
        app.pixelSizeMicrons = 0.5
        app.skipMeasurements = false
        app.distThreshold = 10.0
        app.estimateCellBoundaryDist = 3.0
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
}
