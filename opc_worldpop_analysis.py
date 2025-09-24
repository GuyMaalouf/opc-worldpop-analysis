import geopandas as gpd
import rasterio
from rasterio.mask import mask
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import uniform_filter

# === Load the KML boundary ===
kml_path = "kenya_opc_boundary.kml"
boundary_gdf = gpd.read_file(kml_path, driver='KML').to_crs("EPSG:4326")

# === Load and mask the WorldPop raster ===
raster_path = "ken_ppp_2020.tif"

with rasterio.open(raster_path) as src:
    out_image, out_transform = mask(src, boundary_gdf.geometry, crop=True)
    data_raw = out_image[0]  # people per 0.01 km²
    data_km2 = data_raw * 100  # convert to people/km²

    masked = np.ma.masked_less_equal(data_km2, 0)

    # === Original stats (100m resolution)
    cell_area_km2 = 0.01
    valid_cells = masked.count()
    area_km2 = valid_cells * cell_area_km2
    total_pop = (masked * cell_area_km2).sum()
    avg_density = masked.mean()
    max_density = masked.max()

    print("\n🟪 100m x 100m Results")
    print(f"✅ Area: {area_km2:.2f} km²")
    print(f"✅ Total Population: {total_pop:.0f}")
    print(f"✅ Avg Density: {avg_density:.2f} ppl/km²")
    print(f"✅ Max Density: {max_density:.2f} ppl/km²")

    # === Aggregate to 1km x 1km ===
    factor = 10
    rows, cols = data_km2.shape
    new_rows = rows // factor
    new_cols = cols // factor

    cropped = data_raw[:new_rows * factor, :new_cols * factor]
    reshaped = cropped.reshape((new_rows, factor, new_cols, factor))
    summed_blocks = reshaped.sum(axis=(1, 3))
    density_1km = summed_blocks

    masked_1km = np.ma.masked_less_equal(density_1km, 0)
    area_1km2 = masked_1km.count()
    total_pop_1km = masked_1km.sum()
    avg_density_1km = masked_1km.mean()
    max_density_1km = masked_1km.max()

    print("\n🟦 1km x 1km Aggregated Results")
    print(f"✅ Area: {area_1km2:.0f} km²")
    print(f"✅ Total Population: {total_pop_1km:.0f}")
    print(f"✅ Avg Density: {avg_density_1km:.2f} ppl/km²")
    print(f"✅ Max Density: {max_density_1km:.2f} ppl/km²")

    # === Plotting
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    axs[0].imshow(masked, cmap='plasma', vmin=0, vmax=100)
    axs[0].set_title("100m x 100m Population Density (ppl/km²)")
    axs[0].axis('off')

    axs[1].imshow(masked_1km, cmap='plasma', vmin=0, vmax=100)
    axs[1].set_title("1km x 1km Population Density (ppl/km²)")
    axs[1].axis('off')

    plt.tight_layout()
    plt.show()