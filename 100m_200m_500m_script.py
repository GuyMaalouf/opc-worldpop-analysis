import geopandas as gpd
import rasterio
from rasterio.mask import mask
import matplotlib.pyplot as plt
import numpy as np

# --- Inputs ---
kml_path = "Safety.kml"                 # Boundary KML
raster_path = "ken_ppp_2020.tif"        # WorldPop 'ppp' (people per pixel) at ~100 m native grid

# --- Load boundary ---
boundary_gdf = gpd.read_file(kml_path, driver='KML').to_crs("EPSG:4326")

with rasterio.open(raster_path) as src:
    # Clip and read raster (ppp: people per 100 m pixel, not density)
    out_image, out_transform = mask(src, boundary_gdf.geometry, crop=True)
    ppp = out_image[0].astype(np.float64)  # people per pixel
    nodata = src.nodata

# Build a mask: invalid if nodata or <=0 population (treat zeros as empty)
mask_invalid = np.zeros_like(ppp, dtype=bool)
if nodata is not None:
    mask_invalid |= (out_image[0] == nodata)
mask_invalid |= ~(ppp > 0)

ppp = np.where(mask_invalid, np.nan, ppp)

# Constants for 100 m grid (area = 0.1 km * 0.1 km = 0.01 km²)
cell_area_km2 = 0.01

# Convert to density at 100 m: people/km²
dens_100m = ppp / cell_area_km2  # = ppp / 0.01 = ppp * 100

def summarize_density(density_array, cell_area_km2):
    """Return (area_km2, total_pop, avg_density, max_density)."""
    valid = ~np.isnan(density_array)
    n = np.count_nonzero(valid)
    if n == 0:
        return 0.0, 0.0, np.nan, np.nan
    area = n * cell_area_km2
    # total population = sum(density * area_per_cell)
    total_pop = np.nansum(density_array[valid]) * cell_area_km2
    avg_density = np.nanmean(density_array[valid])
    max_density = np.nanmax(density_array[valid])
    return area, total_pop, avg_density, max_density

def aggregate_mean(density_array, factor):
    """
    Aggregate a density raster by averaging over (factor x factor) non-overlapping blocks.
    Assumes all fine cells have equal area. Ignores blocks that are entirely NaN.
    """
    # Trim to a multiple of factor
    rows, cols = density_array.shape
    R = (rows // factor) * factor
    C = (cols // factor) * factor
    d = density_array[:R, :C]

    # Reshape into (new_rows, factor, new_cols, factor)
    new_rows = R // factor
    new_cols = C // factor
    d4 = d.reshape(new_rows, factor, new_cols, factor)

    # Compute nanmean per block
    # If a whole block is NaN, result is NaN
    block_means = np.nanmean(np.nanmean(d4, axis=3), axis=1)
    return block_means

# --- Summaries ---
area_100, pop_100, avg_100, max_100 = summarize_density(dens_100m, cell_area_km2)
print("\n🟪 100 m × 100 m (native)")
print(f"Area: {area_100:.2f} km²  |  Total Pop: {pop_100:.0f}  |  Avg: {avg_100:.2f} ppl/km²  |  Max: {max_100:.2f} ppl/km²")

# Aggregate to 200 m (factor 2), 500 m (factor 5), 1 km (factor 10)
for label, factor in [("200 m × 200 m", 2), ("500 m × 500 m", 5), ("1 km × 1 km", 10)]:
    dens_coarse = aggregate_mean(dens_100m, factor)
    # Coarse cell area in km² = (0.1 * factor)^2
    cell_area_coarse = (0.1 * factor) ** 2
    area, pop, avg, maxv = summarize_density(dens_coarse, cell_area_coarse)
    print(f"\n🟦 {label}")
    print(f"Area: {area:.2f} km²  |  Total Pop: {pop:.0f}  |  Avg: {avg:.2f} ppl/km²  |  Max: {maxv:.2f} ppl/km²")

# --- Optional quick-look plots (comment out if running headless) ---
vmin, vmax = 0, 95  # common scale for all 3

fig, axs = plt.subplots(1, 3, figsize=(16, 5))
plt.subplots_adjust(wspace=0.001)

im0 = axs[0].imshow(dens_100m, cmap='plasma', vmin=vmin, vmax=vmax)
axs[0].set_title("100 m grid resolution", fontsize=15, pad=10); axs[0].axis('off')

dens_200 = aggregate_mean(dens_100m, 2)
axs[1].imshow(dens_200, cmap='plasma', vmin=vmin, vmax=vmax)
axs[1].set_title("200 m grid resolution", fontsize=15, pad=10); axs[1].axis('off')

dens_500 = aggregate_mean(dens_100m, 5)
axs[2].imshow(dens_500, cmap='plasma', vmin=vmin, vmax=vmax)
axs[2].set_title("500 m grid resolution", fontsize=15, pad=10); axs[2].axis('off')

# --- shared colorbar in its own fixed axis (no effect on subplot spacing) ---
import matplotlib as mpl

vmin, vmax = 0, 100  # or whatever range you’re using
sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), cmap='plasma')
sm.set_array([])  # required for older Matplotlib versions

fig.subplots_adjust(right=0.88)  # leave room on the right for the colorbar
cax = fig.add_axes([0.94, 0.1, 0.015, 0.80])  # [left, bottom, width, height] in figure coords
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label("People per km²", fontsize=15)

# --- Annotate average and max population densities under each plot --- 
stats = [ summarize_density(dens_100m, 0.01), # 100 m grid 
summarize_density(dens_200, (0.1 * 2) ** 2), # 200 m grid 
summarize_density(dens_500, (0.1 * 5) ** 2), # 500 m grid 
] 
for ax, (area, pop, avg, mx) in zip(axs, stats): 
 ax.text( 
  0.5, -0.05, f"Avg: {avg:.2f} ppl/km² - Max: {mx:.2f} ppl/km²", 
  transform=ax.transAxes, 
  ha='center', va='top', fontsize=15 
  )
  
# (keep your annotation code here)

plt.tight_layout(rect=[0, 0, 0.95, 1])  # leave space for the colorbar
plt.show()

