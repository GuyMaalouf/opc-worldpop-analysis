# adjacent_area_stats.py  — v3 (adds 200 m & 500 m; separate plots)
# Requirements: geopandas, rasterio, shapely, numpy, matplotlib

import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from shapely.ops import unary_union
from pathlib import Path
import matplotlib


# ========= USER INPUTS =========
KML_PATH = "data/adjacent_area_5km.kml"
RASTER_PATH = "data/ken_ppp_2020.tif"  # WorldPop people per 100 m x 100 m cell
NAME_5KM = "Adjacent Area 5 km (from contingency)"
NAME_CONTING = "Contingency Volume (from KMZ)"

# Output images
OUT_100M = "figs/adjacent_area_100m.png"
OUT_200M = "figs/adjacent_area_200m.png"
OUT_500M = "figs/adjacent_area_500m.png"

# Plotting style
CMAP_NAME = "plasma"
USE_LOG   = False  # set True to use logarithmic color scale


# ========= HELPERS =========
def get_geom_by_name(gdf: gpd.GeoDataFrame, name: str):
    sel = gdf[gdf["Name"] == name]
    if sel.empty:
        return None
    return unary_union(sel.geometry.values)

def aggregate_density_from_people(people_array: np.ndarray, factor: int):
    """
    Aggregate a 100 m people-per-cell raster to a coarser grid (factor x factor).
    Returns:
      density_m  - masked array (ppl/km²) at coarse resolution
      people_m   - masked array (people per coarse cell)
      cell_area_km2 - coarse cell area in km²
    """
    rows, cols = people_array.shape
    R = (rows // factor) * factor
    C = (cols // factor) * factor
    block = people_array[:R, :C]

    new_rows = R // factor
    new_cols = C // factor
    # reshape to (new_rows, factor, new_cols, factor)
    b = block.reshape(new_rows, factor, new_cols, factor)

    # valid-count per coarse block (ignore NaNs)
    valid_count = np.sum(~np.isnan(b), axis=(1, 3))

    # sum people per coarse block, ignoring NaNs
    people_sum = np.nansum(b, axis=(1, 3))

    # set blocks with 0 valid fine cells to NaN
    people_sum[valid_count == 0] = np.nan

    cell_area_km2 = (0.1 * factor) ** 2  # (0.1 km = 100 m)
    density = people_sum / cell_area_km2

    # mask invalid
    density_m = np.ma.masked_invalid(density)
    people_m  = np.ma.masked_invalid(people_sum)
    return density_m, people_m, cell_area_km2


def positive_min(arrmasked):
    """Smallest positive value for LogNorm vmin."""
    data = np.asarray(arrmasked).astype(float)
    data = data[data > 0]
    return float(data.min()) if data.size else 1e-3

def summarize_density(density_array: np.ma.MaskedArray, cell_area_km2: float):
    n = density_array.count()
    if n == 0:
        return 0.0, 0.0, np.nan, np.nan
    area_km2   = n * cell_area_km2
    total_pop  = float(density_array.sum() * cell_area_km2)
    avg_density = float(density_array.mean())
    max_density = float(density_array.max())
    return area_km2, total_pop, avg_density, max_density

def save_density_plot(arr_m, title, outfile,
                      cmap_name=CMAP_NAME, use_log=USE_LOG,
                      vmin=None, vmax=None):
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 6))
    cmap = plt.colormaps.get_cmap(cmap_name).copy()
    cmap.set_bad(color='none')

    if use_log:
        # keep your log logic if needed
        from matplotlib.colors import LogNorm
        if vmin is None: vmin = positive_min(arr_m)
        if vmax is None: vmax = float(arr_m.max()) if arr_m.count() else vmin * 10
        norm = LogNorm(vmin=vmin, vmax=vmax)
        im = plt.imshow(arr_m, cmap=cmap, norm=norm)
    else:
        im = plt.imshow(arr_m, cmap=cmap, vmin=vmin, vmax=vmax)

    cbar = plt.colorbar(im)
    cbar.set_label("People per km²")
    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()
    print(f"🖼️ Saved: {outfile}")



# ========= LOAD KML & PICK GEOMETRIES =========
kml = gpd.read_file(KML_PATH, driver="KML").to_crs("EPSG:4326")

geom_5km = get_geom_by_name(kml, NAME_5KM)
if geom_5km is None:
    raise RuntimeError(f"Could not find placemark named: {NAME_5KM}")

inner_geom = get_geom_by_name(kml, NAME_CONTING)
if inner_geom is None:
    raise RuntimeError(f"Could not find placemark named: {NAME_CONTING}")

adjacent_area = geom_5km.difference(inner_geom)
if adjacent_area.is_empty:
    raise RuntimeError("Adjacent Area ring is empty after difference; check inputs.")

adjacent_gdf = gpd.GeoDataFrame({"Name": ["Adjacent Area Ring"]},
                                geometry=[adjacent_area],
                                crs="EPSG:4326")


# ========= RASTER MASK & STATS =========
with rasterio.open(RASTER_PATH) as src:
    out_image, out_transform = mask(
        src, adjacent_gdf.geometry, crop=True, filled=True, nodata=0
    )

# WorldPop: people per 100 m cell; 100 m cell area = 0.01 km²
people_100m = out_image[0].astype(float)
people_100m = np.where(people_100m > 0, people_100m, np.nan)

# 100 m density (ppl/km²)
dens_100m = people_100m / 0.01  # = *100
dens_100m_m = np.ma.masked_invalid(dens_100m)

# ---- Stats @ 100 m ----
area_100, pop_100, avg_100, max_100 = summarize_density(dens_100m_m, 0.01)
print("\n🟪 100 m × 100 m (native) — Adjacent Area (5 km ring)")
print(f"✅ Area (approx.): {area_100:.2f} km²")
print(f"✅ Total Population: {pop_100:.0f} people")
print(f"✅ Avg Density: {avg_100:.2f} ppl/km²")
print(f"✅ Max Density: {max_100:.2f} ppl/km²")

# ---- Aggregate to 200 m ----
dens_200m, ppl_200m, cell_area_200 = aggregate_density_from_people(people_100m, factor=2)
area_200 = dens_200m.count() * cell_area_200
total_200 = ppl_200m.sum()
avg_200   = dens_200m.mean()
max_200   = dens_200m.max()

print("\n🟦 200 m × 200 m aggregated — Adjacent Area (5 km ring)")
print(f"✅ Area (approx.): {area_200:.2f} km²")
print(f"✅ Total Population: {total_200:.0f} people")
print(f"✅ Avg Density: {avg_200:.2f} ppl/km²")
print(f"✅ Max Density: {max_200:.2f} ppl/km²")

# ---- Aggregate to 500 m ----
dens_500m, ppl_500m, cell_area_500 = aggregate_density_from_people(people_100m, factor=5)
area_500 = dens_500m.count() * cell_area_500
total_500 = ppl_500m.sum()
avg_500   = dens_500m.mean()
max_500   = dens_500m.max()

print("\n🟦 500 m × 500 m aggregated — Adjacent Area (5 km ring)")
print(f"✅ Area (approx.): {area_500:.2f} km²")
print(f"✅ Total Population: {total_500:.0f} people")
print(f"✅ Avg Density: {avg_500:.2f} ppl/km²")
print(f"✅ Max Density: {max_500:.2f} ppl/km²")

# ========= PLOTS (separate figures) =========
# --- Compute global color scale (shared across all plots) ---
all_values = np.concatenate([
    dens_100m_m.compressed(),
    dens_200m.compressed(),
    dens_500m.compressed()
])
vmin_global = np.nanpercentile(all_values, 0)   # or 0 for strict min
vmax_global = np.nanpercentile(all_values, 100)  # clip top 0.5% outliers

print(f"🎨 Shared color scale: vmin={vmin_global:.2f}, vmax={vmax_global:.2f}")

save_density_plot(dens_100m_m, "Adjacent Area Population Density (100 m grid size)", OUT_100M,
                  vmin=vmin_global, vmax=vmax_global)
save_density_plot(dens_200m,   "Adjacent Area Population Density (200 m grid size)", OUT_200M,
                  vmin=vmin_global, vmax=vmax_global)
save_density_plot(dens_500m,   "Adjacent Area Population Density (500 m grid size)", OUT_500M,
                  vmin=vmin_global, vmax=vmax_global)


