> [!IMPORTANT]
> This codebase is now deprecated. Please see [cellmeasurement-py](https://github.com/WEHI-SODA-Hub/cellmeasurement-py) for the Python rewrite

# Cellmeasurement

This is a Groovy application that extracts Regions of Interest (ROIs) from whole cell and nuclear mask images, matches them based on their centroids, calculates measurements and exports the cell objects to a GeoJSON file for importing into QuPath or other image viewer.

## Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Requirements

- Java Development Kit (JDK) 8 or higher
- Groovy
- QuPath
- ImageJ

As long as you have Java installed, gradle should be able to handle the rest.

## Installation

1. Clone the repository:
   ```sh
   git clone https://github.com/WEHI-SODA-Hub/cellmeasurement
   cd cellmeasurement
   ```

2. Build the project:
   ```sh
   ./gradlew build
   ```

3. Test the project
   ```sh
   ./gradlew test
   ```

## Usage

Here is an example of running the app:

```sh
./gradlew run \
    --args="--nuclear-mask=$PWD/app/src/test/resources/synthetic_test_nuclear.tiff \
            --whole-cell-mask=$PWD/app/src/test/resources/synthetic_test_whole-cell.tiff \
            --tiff-file=$PWD/app/src/test/resources/synthetic_test.ome.tif \
            --output-file=$PWD/segmentation.geojson \
            --skip-measurements=false \
            --percentiles=70,80,90,95,96,97,98,99 \
            --erosion-steps=4,7,11,14,18"
```

Make sure to use absolute paths.

### Erosion Measurements

The `--erosion-steps` option enables spatial analysis of marker distribution from cell edges toward the interior. This is useful for:

- **Membrane vs cytoplasm gradients** - detecting membrane-localized receptors vs cytoplasmic proteins
- **Nuclear envelope analysis** - identifying proteins enriched at the nuclear boundary
- **Quality control** - detecting edge staining artifacts or segmentation issues
- **Cell morphology** - characterizing marker distribution patterns

**Default values:** `4,7,11,14,18` pixels (optimized for Lunaphore COMET imaging at 0.28 μm/pixel, corresponding to ~1.1, 2.0, 2.8, 3.9, 5.0 μm)

**Generated measurements:**
- `{Channel}: {Compartment}: Eroded_{N}px: Mean` - mean intensity after N pixels of erosion
- `{Channel}: {Compartment}: Eroded_{N}px: Median` - median intensity after N pixels of erosion  
- `{Compartment}: Eroded_{N}px: Area_Fraction` - fraction of compartment area remaining (0.0 = fully eroded, 1.0 = no erosion)

**Note:** Small cells may have `Area_Fraction = 0.0` at larger erosion steps, which is expected. These measurements can be filtered during analysis. Set `--erosion-steps=""` to disable erosion measurements.

Full arguments:

```
Usage: cellmeasurement [-hV] [--skip-measurements] [-d=<downsampleFactor>]
                       [-e=<estimateCellBoundaryDist>]
                       [--erosion-steps=<erosionSteps>] -f=<tiffFilePath>
                       [-i=<distThreshold>] -n=<nuclearMaskFilePath>
                       -o=<outputFilePath> [-p=<pixelSizeMicrons>]
                       [--percentiles=<percentiles>] [-t=<threads>]
                       -w=<wholeCellMaskFilePath>
Extract cell measurements from nuclear and whole-cell segmentation masks.
  -d, --downsample-factor=<downsampleFactor>
                            Downsample factor
  -e, --estimate-cell-boundary-dist=<estimateCellBoundaryDist>
                            Where no matching membrane ROI exists, expand the
                              nucleus by this many pixels (default = 3.0)
      --erosion-steps=<erosionSteps>
                            Comma-separated erosion steps in pixels. Default:
                              4,7,11,14,18. Set to empty to disable.
  -f, --tiff-file=<tiffFilePath>
                            TIFF file containing multi-channel image data
  -h, --help                Show this help message and exit.
  -i, --dist-threshold=<distThreshold>
                            Distance threshold (in pixels) for matching ROIs
  -n, --nuclear-mask=<nuclearMaskFilePath>
                            Nuclear segmentation mask file in TIFF format
  -o, --output-file=<outputFilePath>
                            Output path for GeoJSON file
  -p, --pixel-size-microns=<pixelSizeMicrons>
                            Pixel size in microns (default: 0.5)
      --percentiles=<percentiles>
                            Calculate specified comma-separated intensity
                              percentiles. Only works if not skipping
                              measurements. E.g. "70,80,90,95,96,97,98,99"
      --skip-measurements   Skip adding measurements
  -t, --threads=<threads>   Number of threads to use for parallel processing
                              (default: 1)
  -V, --version             Print version information and exit.
  -w, --whole-cell-mask=<wholeCellMaskFilePath>
                            Whole-cell segmentation mask file in TIFF format
```

Calculating cell measurements is the most time-consuming step. If you only want to check the
segmentations first, it is recommended to run with `--skip-measurements=true`.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any changes.

## Test Data

Test data was derived from [nf-core/test-datasets](https://github.com/nf-core/test-datasets)
under the MIT License.

### Source Files
All source files were derived from [nuclear_image.tif](https://github.com/nf-core/test-datasets/blob/modules/data/imaging/segmentation/nuclear_image.tif):

- `app/src/test/resources/synthetic_test.ome.tif` -- derived from [script in spatialproteomics pipeline](https://github.com/WEHI-SODA-Hub/spatialproteomics/blob/main/tests/data/comet/make_comet_test_data.py)
  and run through background subtraction step.
- `app/src/test/resources/synthetic_test_nuclear.tiff` -- `parquettotiff` output from [spatialproteomics pipeline](https://github.com/WEHI-SODA-Hub/spatialproteomics).
- `app/src/test/resources/test_data_whole-cell.tiff` -- as above.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
