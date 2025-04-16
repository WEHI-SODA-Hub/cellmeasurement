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

## Usage

Here is an example of running the app:

```sh
./gradlew run \
    --args="--nuclear-mask=nuclear_mask.tiff \
            --whole-cell-mask=whole_cell_mask.tiff \
            --tiff-file=image.tiff \
            --output-file=annotations.geojson"
```

Full arguments:

```
Usage: cellmeasurement [-hV] [--skip-measurements] [-d=<downsampleFactor>]
                       [-e=<cellExpansion>] [-i=<distThreshold>]
                       -n=<nuclearMaskFilePath> -o=<outputFilePath>
                       [-p=<pixelSizeMicrons>] -t=<tiffFilePath>
                       -w=<wholeCellMaskFilePath>
Extract cell measurements from nuclear and whole-cell segmentation masks.
  -d, --downsample-factor=<downsampleFactor>
                            Downsample factor
  -e, --cell-expansion=<cellExpansion>
                            Expansion factor for cell boundary estimation in
                              pixels (default = 3.0)
  -h, --help                Show this help message and exit.
  -i, --dist-threshold=<distThreshold>
                            Distance threshold (in pixels) for matching ROIs
  -n, --nuclear-mask=<nuclearMaskFilePath>
                            Nuclear segmentation mask file in TIFF format
  -o, --output-file=<outputFilePath>
                            Output path for GeoJSON file
  -p, --pixel-size-microns=<pixelSizeMicrons>
                            Pixel size in microns (default: 0.5)
      --skip-measurements   Skip adding measurements
  -t, --tiff-file=<tiffFilePath>
                            TIFF file containing multi-channel image data
  -V, --version             Print version information and exit.
  -w, --whole-cell-mask=<wholeCellMaskFilePath>
                            Whole-cell segmentation mask file in TIFF format
```

Calculating cell measurements is the most time-consuming step. If you only want to check the
segmentations first, it is recommended to run with `--skip-measurements=true`.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any changes.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
