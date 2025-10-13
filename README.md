# OPC WorldPop Analysis

This repository contains Python scripts and data instructions for reproducing the population density estimates used in the Ol Pejeta Conservancy (OPC) operational and adjacent-area analyses, based on the 2020 WorldPop dataset.

## Overview

The scripts support population exposure assessment in accordance with the **JARUS SORA 2.5** framework (Steps 2 and 8) by:
- Clipping WorldPop rasters to operational and adjacent-area boundaries,
- Converting raw cell counts to population density (people/km²),
- Aggregating from 100 m to 1 km resolution for **MAUP** (Modifiable Areal Unit Problem) analysis,
- Generating summary statistics for each region and resolution, and
- Plotting logarithmic-scale population density maps.

---

## Repository Contents

| File | Description |
|------|-------------|
| `opc_worldpop_analysis.py` | Processes and plots WorldPop data for the **operational area**. |
| `OPC_buffers_from_contingency.kml` | Defines the **1 km** and **5 km** spatial buffers (adjacent areas) derived from the contingency volume. |
| `adjacent_area_stats.py` | Computes statistics and generates log-scale population density maps for the **adjacent area** (5 km ring). |
| `kenya_opc_boundary.kml` | KML boundary for the Ol Pejeta Conservancy. |
| `ken_ppp_2020.tif` | WorldPop Kenya 2020 population raster (not included — download [here](https://hub.worldpop.org/geodata/summary?id=6530)). |

---

## Prerequisites

Install Python dependencies with:

```bash
pip install numpy scipy matplotlib rasterio geopandas shapely
```

Recommended Python version: **3.10+**  
Tested with:
- `numpy==1.24.4`
- `scipy==1.8.0`
- `matplotlib==3.7.2`
- `rasterio==1.3.9`
- `geopandas==0.14.3`
- `shapely==2.0.2`

## Input Files

Download the required files manually:

- WorldPop Kenya 2020 raster (`ken_ppp_2020.tif`)  
  → [WorldPop Dataset Link](https://hub.worldpop.org/geodata/summary?id=6530)

- Ol Pejeta Conservancy boundary KML (`kenya_opc_boundary.kml`)  
  → Use official GIS boundaries or digitise from public maps
  
- Ol Pejeta Conservancy adjacent area boundary KML (`OPC_buffers_from_contingency.kml`)  
  → `OPC_buffers_from_contingency.kml` — defines contingency, 1 km, and 5 km rings.

## Running the Script

1️⃣ Operational Area Analysis
```bash
python opc_worldpop_analysis.py
```

2️⃣ Adjacent Area Analysis
```bash
python adjacent_area_stats.py
```

## Outputs

- **Operational Area (Step 2):** summary statistics and heatmaps
- **Adjacent Area (Step 8):** 100 m and 1 km log-scale heatmaps

## Example Plots

**Operational Area Population Density**
<img width="800" height="600" alt="WorldPop_OPC_v2" src="https://github.com/user-attachments/assets/196d0503-9255-4e17-973e-2e83efec4e9c" />

**Adjacent Area Population Density**

<img width="345" height="350" alt="adjacent_100m" src="https://github.com/user-attachments/assets/80eec103-26be-493b-b186-bf5446a11ec9" /> <img width="345" height="350" alt="adjacent_1km" src="https://github.com/user-attachments/assets/286c1ccf-a22b-4c7e-a3eb-e26d180a754d" />


## License

This project is released under the MIT License. See `LICENSE` for details.

## Citation

If you use this code in your research, please cite the corresponding academic paper (to be linked upon publication).
