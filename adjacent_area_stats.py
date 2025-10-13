# adjacent_area_stats.py
#
# Computes average & maximum population density within the 5 km Adjacent Area
# using the KML generated from the contingency volume.
#
# Requirements: geopandas, rasterio, shapely, numpy, matplotlib (optional for plots)

import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
import matplotlib.pyplot as plt
from shapely.ops import unary_union

# ========= USER INPUTS =========
# Use the KML that was created alongside the KMZ (same geometries, easier to read programmatically)
KML_PATH = "OPC_buffers_from_contingency.kml"

# WorldPop raster (people per 100m x 100m pixel). Adjust to your file:
# e.g., "ken_ppp_2020.tif" or your latest 2025/2020 dataset.
RASTER_PATH = "ken_ppp_2020.tif"

# Layer/placemark names as written in the KML I generated:
NAME_5KM = "Adjacent Area 5 km (from contingency)"
NAME_CONTING = "Contingency Volume (from KMZ)"

# If you have a Ground Risk Buffer layer available in another file (or the original KMZ),
# set NAME_INNER = "Ground Risk Buffer" and provide that geometry instead of the contingency.
USE_GROUND_RISK_BUFFER_IF_AVAILABLE = False  # set True if you load GRB instead of contingency

# ========= LOAD KML & PICK GEOMETRIES =========
# GeoPandas reads KML into a single layer with a "Name" column for placemark names
kml = gpd.read_file(KML_PATH, driver="KML").to_crs("EPSG:4326")

def get_geom_by_name(gdf: gpd.GeoDataFrame, name: str):
    sel = gdf[gdf["Name"] == name]
    if sel.empty:
        return None
    # Merge in case the placemark has multiple parts
    return unary_union(sel.geometry.values)

geom_5km = get_geom_by_name(kml, NAME_5KM)
if geom_5km is None:
    raise RuntimeError(f"Could not find placemark named: {NAME_5KM}")

# Inner geometry: by default use the contingency polygon from the KML we created
inner_geom = get_geom_by_name(kml, NAME_CONTING)
if inner_geom is None:
    raise RuntimeError(f"Could not find placemark named: {NAME_CONTING}")

# If you want to subtract the Ground Risk Buffer instead (recommended, if available):
# if USE_GROUND_RISK_BUFFER_IF_AVAILABLE:
#     # Example: read another KML/KMZ that contains the GRB, then:
#     # grb_kml = gpd.read_file("your_original_file.kml", driver="KML").to_crs("EPSG:4326")
#     # inner_geom = get_geom_by_name(grb_kml, "Ground Risk Buffer")
#     # if inner_geom is None: raise RuntimeError("GRB layer not found")
#     pass

# Build the Adjacent Area ring: 5 km polygon minus inner polygon
adjacent_area = geom_5km.difference(inner_geom)
if adjacent_area.is_empty:
    raise RuntimeError("Adjacent Area ring is empty after difference; check inputs.")

adjacent_gdf = gpd.GeoDataFrame({"Name": ["Adjacent Area Ring"]},
                                geometry=[adjacent_area],
                                crs="EPSG:4326")

# ========= RASTER MASK & STATS =========
with rasterio.open(RASTER_PATH) as src:
    # Mask (clip) the raster to the adjacent area
    out_image, out_transform = mask(
        src,
        adjacent_gdf.geometry,  # list-like is accepted
        crop=True,
        filled=True,
        nodata=0
    )

    # WorldPop typical convention: band 1 = people per 100m x 100m cell (0.01 km^2)
    data_people_per_cell = out_image[0]

    # Convert to density (people per km^2)
    # Each 100m x 100m pixel = 0.01 km^2 → density = (people_per_cell / 0.01) = *100
    data_density_km2_100m = data_people_per_cell * 100.0

    # Mask zero/negative as "no data"
    masked_100m = np.ma.masked_less_equal(data_density_km2_100m, 0)

    # --- 100 m results (density already in ppl/km^2) ---
    cell_area_km2 = 0.01
    valid_cells = masked_100m.count()
    area_km2 = valid_cells * cell_area_km2

    # Total population (sum people per cell)
    total_pop = (np.ma.masked_less_equal(data_people_per_cell, 0)).sum()

    avg_density_100m = masked_100m.mean()
    max_density_100m = masked_100m.max()

    print("\n🟪 100 m x 100 m (native) resolution — Adjacent Area (5 km ring)")
    print(f"✅ Area (approx.): {area_km2:.2f} km²")
    print(f"✅ Total Population: {total_pop:.0f} people")
    print(f"✅ Avg Density: {avg_density_100m:.2f} ppl/km²")
    print(f"✅ Max Density: {max_density_100m:.2f} ppl/km²")

    # --- Aggregate to 1 km x 1 km ---
    # We aggregate "people per 100m cell" by summing 10x10 blocks → people per 1 km² pixel,
    # then treat that as density at 1 km resolution.
    factor = 10  # 10 * 100 m = 1 km

    rows, cols = data_people_per_cell.shape
    new_rows = rows // factor
    new_cols = cols // factor

    # Crop to multiples of 10 to form full 1 km blocks
    cropped = data_people_per_cell[:new_rows * factor, :new_cols * factor]

    # Reshape to (new_rows, factor, new_cols, factor) and sum over the small dims
    reshaped = cropped.reshape((new_rows, factor, new_cols, factor))
    summed_blocks = reshaped.sum(axis=(1, 3))  # each pixel = people per 1 km²

    # Mask zero/negative as no data
    density_1km = np.ma.masked_less_equal(summed_blocks, 0)

    area_1km2 = density_1km.count()  # each valid pixel ~1 km²
    total_pop_1km = density_1km.sum()
    avg_density_1km = density_1km.mean()  # ppl/km²
    max_density_1km = density_1km.max()   # ppl/km²

    print("\n🟦 1 km x 1 km aggregated — Adjacent Area (5 km ring)")
    print(f"✅ Area (approx.): {area_1km2:.0f} km²")
    print(f"✅ Total Population: {total_pop_1km:.0f} people")
    print(f"✅ Avg Density: {avg_density_1km:.2f} ppl/km²")
    print(f"✅ Max Density: {max_density_1km:.2f} ppl/km²")

    # Optional quicklook plots (comment out if running headless)
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))

    axs[0].imshow(masked_100m, cmap='plasma', vmin=0)
    axs[0].set_title("100 m density (ppl/km²) — Adjacent Area")
    axs[0].axis('off')

    axs[1].imshow(density_1km, cmap='plasma', vmin=0)
    axs[1].set_title("1 km aggregated density (ppl/km²) — Adjacent Area")
    axs[1].axis('off')

    plt.tight_layout()
    plt.show()

