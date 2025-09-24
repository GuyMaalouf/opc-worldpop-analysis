# Population Density Analysis for Ol Pejeta Conservancy (WorldPop 2020)

This repository contains a Python script and instructions for processing gridded population data within the Ol Pejeta Conservancy in Kenya, using the 2020 WorldPop raster dataset.

## 📂 Files Required

- `kenya_opc_boundary.kml`: Polygon boundary of Ol Pejeta Conservancy.
    - You can export this from Google Earth or GIS tools.
- `ken_ppp_2020.tif`: 2020 WorldPop raster for Kenya. Download from:
    - [https://hub.worldpop.org/geodata/summary?id=6530]([url](https://hub.worldpop.org/geodata/summary?id=6530))

Place both files in the root directory before running the script.

## 🧪 Output

- Average and maximum population density within the conservancy at:
  - 100m x 100m resolution
  - 1km x 1km resolution (aggregated)
- Side-by-side plots for visual analysis

## 🐍 Dependencies

We recommend running inside a clean virtual environment. Required Python packages:

```bash
pip install numpy==1.24.4 scipy==1.8.0 geopandas rasterio matplotlib
