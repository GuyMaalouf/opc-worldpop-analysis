#!/usr/bin/env python3
"""Recalculate 2026 adjacent-area WorldPop exposure using the 2025 method.

The 0-1 km outdoor-assembly band and the 1-5 km remaining adjacent-area band
are unioned before raster masking. The native WorldPop product and the
100 m -> 1 km aggregation logic reproduce the 2025 processing. A separate
area-preserving mean is also reported because aggregation should not change
the mean when the same population and geographic area are retained.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from exactextract import exact_extract
from pyproj import Geod
import rasterio
from rasterio.mask import mask
from shapely import wkt
from shapely.geometry import mapping


def geodesic_area_km2(geom) -> float:
    geod = Geod(ellps="WGS84")
    area_m2, _ = geod.geometry_area_perimeter(geom)
    return abs(area_m2) / 1_000_000.0


def aggregate_density_from_people(people_array: np.ndarray, factor: int):
    """Replicate the 2025 block aggregation from people per native cell."""
    rows, cols = people_array.shape
    r = (rows // factor) * factor
    c = (cols // factor) * factor
    block = people_array[:r, :c]
    reshaped = block.reshape(r // factor, factor, c // factor, factor)
    valid_count = np.sum(~np.isnan(reshaped), axis=(1, 3))
    people_sum = np.nansum(reshaped, axis=(1, 3))
    people_sum[valid_count == 0] = np.nan
    cell_area_km2 = (0.1 * factor) ** 2
    density = people_sum / cell_area_km2
    return np.ma.masked_invalid(density), np.ma.masked_invalid(people_sum), cell_area_km2


def summarize_density(density: np.ma.MaskedArray, cell_area_km2: float) -> dict:
    count = int(density.count())
    if count == 0:
        raise RuntimeError("No populated raster cells intersect the assessment area")
    return {
        "effective_area_km2": count * cell_area_km2,
        "total_population": float(density.sum() * cell_area_km2),
        "mean_density_people_per_km2": float(density.mean()),
        "max_density_people_per_km2": float(density.max()),
        "included_cells": count,
    }


def exact_population_sum(raster_path: Path, geom) -> float:
    gdf = gpd.GeoDataFrame({"id": [1]}, geometry=[geom], crs="EPSG:4326")
    result = exact_extract(str(raster_path), gdf, ["sum"], output="pandas")
    return float(result.loc[0, "sum"])


def save_plot(array, title: str, output: Path, mean_density: float, max_density: float, norm) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    cmap = plt.colormaps.get_cmap("plasma").copy()
    cmap.set_bad(color="white", alpha=0)
    image = ax.imshow(array, cmap=cmap, norm=norm)
    cbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.025)
    cbar.set_label("Population density (people/km²)")
    ax.set_title(title, pad=10)
    ax.axis("off")
    ax.text(
        0.02,
        0.02,
        f"Area-wide mean: {mean_density:.1f} people/km²\nMaximum: {max_density:.1f} people/km²",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.9, "edgecolor": "0.5"},
    )
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_combined_plot(array_100m, array_1km, output: Path, mean_density: float, max_100m: float, max_1km: float, norm) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 6.1))
    cmap = plt.colormaps.get_cmap("plasma").copy()
    cmap.set_bad(color="white", alpha=0)
    images = []
    for ax, array, title, maximum in [
        (axes[0], array_100m, "100 m grid", max_100m),
        (axes[1], array_1km, "1 km grid", max_1km),
    ]:
        image = ax.imshow(array, cmap=cmap, norm=norm)
        images.append(image)
        ax.set_title(title, pad=8)
        ax.axis("off")
        ax.text(
            0.02,
            0.02,
            f"Mean: {mean_density:.1f} people/km²\nMaximum: {maximum:.1f} people/km²",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.9, "edgecolor": "0.5"},
        )
    cbar = fig.colorbar(images[0], ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("Population density (people/km²)")
    fig.suptitle("2026 Adjacent-Area Population Density", y=0.98)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wkt", type=Path, required=True)
    parser.add_argument("--raster", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    wkt_text = args.wkt.read_text(encoding="utf-8").strip()
    geometry = wkt.loads(wkt_text)
    if geometry.is_empty or not geometry.is_valid:
        raise RuntimeError("Combined 2026 adjacent-area geometry is invalid")

    polygon_area_km2 = geodesic_area_km2(geometry)

    with rasterio.open(args.raster) as src:
        out_image, _ = mask(src, [mapping(geometry)], crop=True, filled=True, nodata=0)

    people_raw = out_image[0].astype(float)
    people_100m = np.where(people_raw > 0, people_raw, np.nan)
    density_100m = np.ma.masked_invalid(people_100m / 0.01)
    stats_100m_compatible = summarize_density(density_100m, 0.01)

    density_1km, people_1km, cell_area_1km = aggregate_density_from_people(people_100m, factor=10)
    stats_1km_compatible = summarize_density(density_1km, cell_area_1km)

    population_exact = exact_population_sum(args.raster, geometry)
    area_wide_mean = population_exact / polygon_area_km2

    positive_values = np.concatenate([density_100m.compressed(), density_1km.compressed()])
    positive_values = positive_values[positive_values > 0]
    norm = LogNorm(vmin=float(positive_values.min()), vmax=float(positive_values.max()))

    args.outdir.mkdir(parents=True, exist_ok=True)
    save_plot(
        density_100m,
        "2026 Adjacent-Area Population Density (100 m grid)",
        args.outdir / "adjacent_area_2026_100m.png",
        area_wide_mean,
        stats_100m_compatible["max_density_people_per_km2"],
        norm,
    )
    save_plot(
        density_1km,
        "2026 Adjacent-Area Population Density (1 km grid)",
        args.outdir / "adjacent_area_2026_1km.png",
        area_wide_mean,
        stats_1km_compatible["max_density_people_per_km2"],
        norm,
    )
    save_combined_plot(
        density_100m,
        density_1km,
        args.outdir / "adjacent_area_2026_100m_1km.png",
        area_wide_mean,
        stats_100m_compatible["max_density_people_per_km2"],
        stats_1km_compatible["max_density_people_per_km2"],
        norm,
    )

    output = {
        "input": {
            "worldpop_file": args.raster.name,
            "worldpop_product": "WorldPop unconstrained Kenya 2020 population, people per native cell",
            "geometry": "Union of the 0-1 km outdoor-assemblies band and 1-5 km remaining adjacent-area band",
            "geometry_sha256": hashlib.sha256(wkt_text.encode("utf-8")).hexdigest(),
        },
        "combined_geometry": {
            "geodesic_area_km2": polygon_area_km2,
        },
        "recommended_area_preserving_result": {
            "population_fractional_cell_coverage": population_exact,
            "mean_density_people_per_km2": area_wide_mean,
            "max_density_100m_people_per_km2": stats_100m_compatible["max_density_people_per_km2"],
            "max_density_1km_people_per_km2": stats_1km_compatible["max_density_people_per_km2"],
        },
        "2025_compatible_processing_diagnostics": {
            "native_100m": stats_100m_compatible,
            "aggregated_1km": stats_1km_compatible,
            "population_lost_by_1km_edge_cropping": stats_100m_compatible["total_population"] - stats_1km_compatible["total_population"],
            "note": "These means reproduce the 2025 positive-cell/occupied-block denominators and are retained only as diagnostics. The area-preserving mean is the appropriate common mean across resolutions.",
        },
        "outputs": [
            "adjacent_area_2026_100m.png",
            "adjacent_area_2026_1km.png",
            "adjacent_area_2026_100m_1km.png",
        ],
    }

    result_path = args.outdir / "adjacent_area_2026_results.json"
    result_path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
