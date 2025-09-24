# OPC WorldPop Analysis

This repository contains a Python script and necessary data instructions to reproduce population density estimates for Ol Pejeta Conservancy (OPC), using the 2020 WorldPop dataset.

## Overview

The script performs the following tasks:
- Clips the WorldPop raster dataset to the OPC boundary
- Converts cell values to population density (people/km²)
- Aggregates data from 100m to 1km resolution for MAUP analysis
- Outputs summary statistics for each resolution
- Plots density heatmaps

This supports the population exposure assessment in accordance with the JARUS SORA 2.5 framework (Step 2: Ground Risk).

## Repository Contents

- `opc_worldpop_analysis.py` — Main Python script for processing and plotting
- `kenya_opc_boundary.kml` — KML polygon defining the Ol Pejeta Conservancy boundary
- `ken_ppp_2020.tif` — WorldPop Kenya raster dataset (not included - can be downloaded [here](https://hub.worldpop.org/geodata/summary?id=6530))

## Prerequisites

Install Python dependencies using pip:

```bash
pip install geopandas rasterio numpy matplotlib scipy
```

Recommended Python version: **3.10+**  
Tested with:
- `numpy==1.24.4`
- `scipy==1.8.0`
- `matplotlib==3.7.2`
- `rasterio==1.3.9`
- `geopandas==0.14.3`

## Input Files

Download the required files manually:

- WorldPop Kenya 2020 raster (`ken_ppp_2020.tif`)  
  → [WorldPop Dataset Link](https://hub.worldpop.org/geodata/summary?id=6530)

- Ol Pejeta Conservancy boundary KML (`kenya_opc_boundary.kml`)  
  → Use official GIS boundaries or digitise from public maps

## Running the Script

```bash
python opc_worldpop_analysis.py
```

## Outputs

- Summary statistics printed to terminal for both 100m and 1km resolution
- Side-by-side density heatmaps shown in a matplotlib window

## License

This project is released under the MIT License. See `LICENSE` for details.

## Citation

If you use this code in your research, please cite the corresponding academic paper (to be linked upon publication).
