# operational_density_distribution.py
# Create a distribution (histogram) of population density for the operational area
# Requirements: geopandas, rasterio, numpy, matplotlib

import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
import matplotlib.pyplot as plt
from shapely.ops import unary_union
from pathlib import Path

# ========= USER INPUTS =========
KML_PATH = "Safety.kml"                 # your operational KML
RASTER_PATH = "ken_ppp_2020.tif"        # WorldPop 'ppp' (people per 100 m pixel)
NAME_OPERATIONAL = None                 # e.g. "Operational", or None to use all polygons
OUT_PLOT = "operational_density_distribution.png"

# Histogram settings
BIN_WIDTH = 2.0     # ppl/km² per bin
SHOW_CUMULATIVE = True   # also show cumulative % as a line on a twin y-axis

# ========= HELPERS =========
def get_geom(gdf: gpd.GeoDataFrame, name_substring: str | None):
    """Return a single polygon: unary_union of matched KML geometries.
    If name_substring is None, union all; otherwise do case-insensitive substring match."""
    if name_substring is None:
        sel = gdf
    else:
        sel = gdf[gdf["Name"].astype(str).str.contains(name_substring, case=False, na=False)]
        if sel.empty:
            # If not found, list names to help user
            print("\nAvailable placemark names in KML:")
            for n in sorted(gdf["Name"].astype(str).unique()):
                print(" -", n)
            raise ValueError(f'Could not find placemark containing: "{name_substring}"')
    return unary_union(sel.geometry.values)

# ========= LOAD OPERATIONAL AREA =========
kml = gpd.read_file(KML_PATH, driver="KML").to_crs("EPSG:4326")
operational_geom = get_geom(kml, NAME_OPERATIONAL)
operational_gdf = gpd.GeoDataFrame({"Name": ["Operational Area"]},
                                   geometry=[operational_geom],
                                   crs="EPSG:4326")

# ========= RASTER MASK =========
with rasterio.open(RASTER_PATH) as src:
    out_image, out_transform = mask(src, operational_gdf.geometry, crop=True, filled=True, nodata=0)
    ppp = out_image[0].astype(np.float64)  # WorldPop people per 100 m cell
    nodata = src.nodata

# Mask invalid (nodata or <=0)
invalid = np.zeros_like(ppp, dtype=bool)
if nodata is not None:
    invalid |= (out_image[0] == nodata)
invalid |= ~(ppp > 0)  # treat zeros/non-positive as empty
ppp = np.where(invalid, np.nan, ppp)

# Convert to density: ppl/km² (each 100 m cell is 0.01 km²)
density = ppp * 100.0

# Flatten and drop NaNs
vals = density.ravel()
vals = vals[~np.isnan(vals)]
if vals.size == 0:
    raise RuntimeError("No valid density pixels found inside the operational area.")

# ========= BUILD HISTOGRAM (10 ppl/km² bins) =========
vmax = np.nanmax(vals)
# start bins at 0, end at ceil to include max
bins = np.arange(0, np.ceil(vmax / BIN_WIDTH) * BIN_WIDTH + BIN_WIDTH, BIN_WIDTH)
counts, edges = np.histogram(vals, bins=bins)
percentages = (counts / counts.sum()) * 100.0
centers = edges[:-1]  # left-aligned bars (we'll align='edge')

# Quick stats
mean_val = float(np.nanmean(vals))
median_val = float(np.nanmedian(vals))
max_val = float(vmax)
n_at_max_bin = int(counts[-2] if max_val < edges[-1] else counts[-1])  # conservative
# Alternatively: exact-match count (rarely all pixels will equal exact max)
n_exact_max = int(np.sum(np.isclose(vals, max_val, rtol=0, atol=1e-9)))

print("\n📊 Operational Area Density Distribution")
print(f"Pixels (valid): {vals.size}")
print(f"Mean:   {mean_val:.2f} ppl/km²")
print(f"Median: {median_val:.2f} ppl/km²")
print(f"Max:    {max_val:.2f} ppl/km²")
print(f'Pixels exactly at max value: {n_exact_max}')
print(f'Pixels in top bin (>= {edges[-2]:.0f}): {int(counts[-1])}')

# ========= PLOT (robust legend) =========
fig, ax = plt.subplots(figsize=(15, 6))
ax.tick_params(axis='both', which='major', labelsize=14)

# histogram
bars = ax.bar(centers, percentages, width=BIN_WIDTH * 0.9, align='edge',
              edgecolor='k', linewidth=0.5, color='skyblue', label="Density distribution")

ax.set_xlabel("Population density (people per km²)",fontsize=14)
ax.set_ylabel("Percentage of area (%)",fontsize=14)
#ax.set_title("Population Density Distribution — Operational Area")
ax.grid(axis='y', linestyle='--', alpha=0.4)

# reference lines (capture handles!)
p99_val = np.nanpercentile(vals, 99)
line_mean   = ax.axvline(mean_val,   color='tab:red',    linestyle='--', linewidth=2, label=f"Mean = {mean_val:.1f}")
line_median = ax.axvline(median_val, color='tab:orange', linestyle='-.', linewidth=2, label=f"Median = {median_val:.1f}")
line_p99    = ax.axvline(p99_val,    color='tab:green',      linestyle=':',  linewidth=2, label=f"99th percentile = {p99_val:.1f}")
line_max = ax.axvline(max_val, color='gold', linestyle='--', linewidth=2, label=f"Max = {max_val:.1f}")


handles = [bars, line_mean, line_median, line_max, line_p99]
labels  = [h.get_label() for h in handles]

# optional cumulative on twin axis (capture its handle explicitly)
if SHOW_CUMULATIVE:
    ax2 = ax.twinx()
    ax2.tick_params(axis='both', which='major', labelsize=14)
    cum = np.cumsum(percentages)
    (line_cum,) = ax2.plot(edges[:-1] + BIN_WIDTH * 0.45, cum, linewidth=2, color='navy', label="Cumulative (%)")
    ax2.set_ylabel("Cumulative (%)",fontsize=14)
    ax2.set_ylim(0, 100)

    handles.append(line_cum)
    labels.append(line_cum.get_label())

# -- after creating fig, ax, bars, lines, and optional ax2 cumulative:
ax.set_xlim(0, 110)  # start at zero, end at last bin edge

legend = ax.legend(
    handles, labels,
    loc='upper left',
    bbox_to_anchor=(1.05, 1.0),
    fontsize=14,
    frameon=False,
    borderaxespad=0.0
)

# leave a little room on right; bbox_inches='tight' will trim to legend edge
plt.tight_layout(rect=[0, 0, 0.88, 1])
plt.savefig(OUT_PLOT, dpi=200, bbox_inches='tight', pad_inches=0.1)
plt.close()


print(f"✅ Saved: {OUT_PLOT}")
print(f"📈 99th percentile value: {p99_val:.2f} ppl/km²")

