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

To run the application, use the following command:

```sh
./gradlew run --args="<wholeCellMaskFilePath> <nuclearMaskFilePath> <outputFilePath>"
```

- `wholeCellMaskFilePath`: Path to the whole cell mask image file.
- `nuclearMaskFilePath`: Path to the nuclear mask image file.
- `outputFilePath`: Path to the output GeoJSON file.

Example:

```sh
./gradlew run --args="data/whole_cell_mask.tif data/nuclear_mask.tif output/cell_objects.geojson"
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any changes.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
