# OPC WorldPop Analysis

Population-exposure analysis for Ol Pejeta Conservancy (OPC) using WorldPop gridded data.  
This repository reproduces the **top-down (raster)** components of the study and outputs:
- Multi-resolution (100/200/500 m) **population density maps** for the operational area and the adjacent ring
- A **density distribution** (histogram) for the operational area with mean/median/99th-percentile overlays

> **SORA Annex F alignment:** Throughout, we follow Annex F guidance to match the **population grid size** to the **aircraft’s dispersion area**. For flights below 500 ft AGL, a **200 m** grid is the optimal scale; **500 m** is an acceptable coarser alternative; 1 km is outside the recommended maximum for this operational height. We refer to this as **grid-size sensitivity** (not “MAUP”).

---

## Project structure

```
data/
  Safety.kml
  OPC_buffers_from_contingency.kml
  ken_ppp_2020.tif   # or ken_ppp_2018.tif (see note below)
figs/
  # output images are written here
scripts/
  opc_multigrid_density_maps.py
  adjacent_ring_multigrid_density.py
  opc_density_distribution.py
```

**Note on WorldPop vintage:**  
The paper references WorldPop **2018**; the scripts default to `ken_ppp_2020.tif`. You can use either (or both). Just point `--raster` to the file you have. WorldPop (unconstrained, ppp) is available at: https://hub.worldpop.org/

---

## Installation

We recommend **conda/mamba** to avoid GDAL/PROJ headaches.

```bash
mamba create -n opc-pop python=3.10 -y
mamba activate opc-pop

mamba install -y geopandas rasterio shapely matplotlib numpy
# Optional (nice to have): scipy
# mamba install -y scipy
```

Alternatively with `pip` (works, but ensure system GDAL/PROJ are present):

```bash
python -m venv .venv && source .venv/bin/activate
pip install geopandas rasterio shapely matplotlib numpy
# pip install scipy
```

---

## Usage

All scripts accept simple path arguments; defaults are sensible for the repo layout.

### 1) Multi-grid maps for **operational area**
Generates 100/200/500 m density maps and prints summary stats (area, total pop, avg, max).

```bash
python scripts/opc_multigrid_density_maps.py   --kml data/Safety.kml   --raster data/ken_ppp_2020.tif   --outdir figs/   --vmin 0 --vmax 100
```

**Outputs (examples):**
- `figs/opc_density_100m.png`
- `figs/opc_density_200m.png`
- `figs/opc_density_500m.png`

Each panel is annotated with **Avg** and **Max** ppl/km².

### 2) Multi-grid maps for **adjacent ring (5 km minus contingency)**
```bash
python scripts/adjacent_ring_multigrid_density.py   --kml data/OPC_buffers_from_contingency.kml   --raster data/ken_ppp_2020.tif   --name-outer "Adjacent Area 5 km (from contingency)"   --name-inner "Contingency Volume (from KMZ)"   --outdir figs/   --vmin 0 --vmax 900
```

### 3) **Density distribution** for operational area (100 m grid)
Shows histogram (% of area per density bin), with mean/median/99th-percentile and optional cumulative curve.

```bash
python scripts/opc_density_distribution.py   --kml data/Safety.kml   --raster data/ken_ppp_2020.tif   --bin-width 1   --show-cumulative
```

**Why this matters:** The operational distribution reveals that **99%** of cells lie below ~**36.7 ppl/km²**, while a few isolated pixels can exceed **100 ppl/km²**. Because Step 2 of SORA uses the **single maximum cell** to set iGRC, such **outliers** can disproportionately affect classification. A justified percentile threshold (e.g., 99th) gives a **more representative** “typical maximum exposure”. We don’t prescribe a cutoff; we highlight the need to **identify and account for outliers**.

---

## Interpreting the results (Annex F grid-size sensitivity)

- **100 m** grid captures **localized peaks** near settlements and can **inflate the maximum**.
- **200 m** grid (optimal for < 500 ft AGL) smooths small clusters while preserving meaningful spatial structure.
- **500 m** grid further smooths peaks and **lowers the maximum**, still within recommended scale.

Example (operational area, WorldPop):
- Max at **100 m**: **104.84** ppl/km²  
- Max at **200 m**: **58.07** ppl/km²  
- Max at **500 m**: **53.35** ppl/km²

---

## Figures (suggested for README)

- **Multi-grid maps (100/200/500 m)** for OPC: highlight grid-size sensitivity  
  *(insert your updated PNG here)*
- **Operational density distribution** (100 m): mean/median/99th percentile with cumulative line  
  *(insert `figs/operational_density_distribution.png`)*

---

## Data & reproducibility

- **Source**: WorldPop Unconstrained Population (people per pixel, ~100 m at equator).  
  Each 100 m cell = **0.01 km²**, so ppl/km² = ppp × **100**.
- **CRS**: inputs are converted to `EPSG:4326` for masking; raster pixel size is in geographic degrees but the **ppp product** is already per ~100 m cell (WorldPop’s native gridding).  
- **Outputs**: saved under `figs/` by default.  
- **Large files**: keep GeoTIFFs in `data/` and consider adding them to `.gitignore`.

---

## Terminology

We avoid calling this “MAUP” in the code/docs. The effect we evaluate is **grid-size sensitivity** per **Annex F**: the choice of **population-grid size** relative to **dispersion area** materially affects **maximum-cell density** and hence iGRC.

---

## License

MIT (or your preferred license). Please cite WorldPop when using the datasets.
