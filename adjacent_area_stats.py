# adjacent_area_stats.py  — v2 (separate plots + log scale)
# Requirements: geopandas, rasterio, shapely, numpy, matplotlib

import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from shapely.ops import unary_union
from pathlib import Path

# ========= USER INPUTS =========
KML_PATH = "OPC_buffers_from_contingency.kml"
RASTER_PATH = "ken_ppp_2020.tif"  # people per 100m x 100m cell
NAME_5KM = "Adjacent Area 5 km (from contingency)"
NAME_CONTING = "Contingency Volume (from KMZ)"

# Output images
OUT_100M = "adjacent_area_100m_log.png"
OUT_1KM  = "adjacent_area_1km_log.png"

# ========= LOAD KML & PICK GEOMETRIES =========
kml = gpd.read_file(KML_PATH, driver="KML").to_crs("EPSG:4326")

def get_geom_by_name(gdf: gpd.GeoDataFrame, name: str):
    sel = gdf[gdf["Name"] == name]
    if sel.empty:
        return None
    return unary_union(sel.geometry.values)

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

    people_per_cell = out_image[0]  # WorldPop people per 100m cell

    # Convert to density (people/km²): each 100m cell area = 0.01 km²
    density_100m = people_per_cell * 100.0
    masked_100m = np.ma.masked_less_equal(density_100m, 0)

    # ----- Stats @ 100 m -----
    cell_area_km2 = 0.01
    area_km2 = masked_100m.count() * cell_area_km2
    total_pop = np.ma.masked_less_equal(people_per_cell, 0).sum()
    avg_density_100m = masked_100m.mean()
    max_density_100m = masked_100m.max()

    print("\n🟪 100 m x 100 m (native) — Adjacent Area (5 km ring)")
    print(f"✅ Area (approx.): {area_km2:.2f} km²")
    print(f"✅ Total Population: {total_pop:.0f} people")
    print(f"✅ Avg Density: {avg_density_100m:.2f} ppl/km²")
    print(f"✅ Max Density: {max_density_100m:.2f} ppl/km²")

    # ----- Aggregate to 1 km -----
    factor = 10
    rows, cols = people_per_cell.shape
    new_rows = rows // factor
    new_cols = cols // factor
    cropped = people_per_cell[:new_rows*factor, :new_cols*factor]
    reshaped = cropped.reshape((new_rows, factor, new_cols, factor))
    people_1km = reshaped.sum(axis=(1, 3))    # people per 1 km² pixel
    density_1km = np.ma.masked_less_equal(people_1km, 0)  # ppl/km²

    area_1km2 = density_1km.count()  # ~km²
    total_pop_1km = density_1km.sum()
    avg_density_1km = density_1km.mean()
    max_density_1km = density_1km.max()

    print("\n🟦 1 km x 1 km aggregated — Adjacent Area (5 km ring)")
    print(f"✅ Area (approx.): {area_1km2:.0f} km²")
    print(f"✅ Total Population: {total_pop_1km:.0f} people")
    print(f"✅ Avg Density: {avg_density_1km:.2f} ppl/km²")
    print(f"✅ Max Density: {max_density_1km:.2f} ppl/km²")

    # ========= PLOTS (separate figures with logarithmic color scale) =========
    # Helper: get a positive vmin for LogNorm
    def positive_min(arrmasked):
        data = np.asarray(arrmasked).astype(float)
        data = data[data > 0]
        return float(data.min()) if data.size else 1e-3

    # ---- Figure 1: 100 m density ----
    vmin_100 = positive_min(masked_100m)
    vmax_100 = float(masked_100m.max()) if masked_100m.count() else 1.0

    plt.figure(figsize=(10, 6))
    cmap = plt.cm.plasma
    cmap.set_bad(color='none')  # Transparent for no-data
    im = plt.imshow(masked_100m, cmap=cmap, vmin=0)
    cbar = plt.colorbar(im)
    cbar.set_label("People per km²")
    plt.title("Population Density (100 m), Adjacent Area")
    plt.axis("off")
    Path(OUT_100M).parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(OUT_100M, dpi=200)
    plt.close()

    # ---- Figure 2: 1 km density ----
    vmin_1k = positive_min(density_1km)
    vmax_1k = float(density_1km.max()) if density_1km.count() else 1.0

    plt.figure(figsize=(10, 6))
    cmap = plt.cm.plasma
    cmap.set_bad(color='none')  # Transparent for no-data
    im2 = plt.imshow(density_1km, cmap=cmap, vmin=0)
    cbar2 = plt.colorbar(im2)
    cbar2.set_label("People per km²")
    plt.title("Population Density (1 km), Adjacent Area")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_1KM, dpi=200)
    plt.close()

    print(f"\n🖼️ Saved: {OUT_100M}")
    print(f"🖼️ Saved: {OUT_1KM}")
