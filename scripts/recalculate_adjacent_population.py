#!/usr/bin/env python3
"""Recalculate adjacent-area WorldPop exposure consistently.

This script:
1. reproduces the published 100 m and 1 km statistics using the original logic;
2. calculates the adjacent-area population using fractional raster-cell coverage;
3. divides that population by the exact polygon area;
4. reports centre-cell and all-touched boundary sensitivities.

WorldPop raster values are people per grid cell.
"""

from __future__ import annotations

import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from exactextract import exact_extract
from pyproj import Geod, Transformer
from rasterio.mask import mask
from shapely.geometry import Polygon, mapping
from shapely.ops import transform


KML_NS = {"kml": "http://www.opengis.net/kml/2.2"}


def parse_kml_polygon(path: Path, placemark_name: str) -> Polygon:
    root = ET.parse(path).getroot()
    for placemark in root.findall(".//kml:Placemark", KML_NS):
        name_node = placemark.find("kml:name", KML_NS)
        if name_node is None or name_node.text != placemark_name:
            continue
        coord_node = placemark.find(
            ".//kml:Polygon/kml:outerBoundaryIs/kml:LinearRing/kml:coordinates",
            KML_NS,
        )
        if coord_node is None or not coord_node.text:
            raise RuntimeError(f"No polygon coordinates found for {placemark_name!r}")
        coords = []
        for token in coord_node.text.split():
            lon, lat, *_ = token.split(",")
            coords.append((float(lon), float(lat)))
        geom = Polygon(coords)
        if not geom.is_valid:
            geom = geom.buffer(0)
        if geom.is_empty or not geom.is_valid:
            raise RuntimeError(f"Invalid polygon for {placemark_name!r}")
        return geom
    raise RuntimeError(f"Placemark not found: {placemark_name!r}")


def geodesic_area_km2(geom: Polygon) -> float:
    geod = Geod(ellps="WGS84")
    area_m2, _ = geod.geometry_area_perimeter(geom)
    return abs(area_m2) / 1_000_000.0


def projected_area_km2(geom: Polygon) -> float:
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32637", always_xy=True)
    geom_utm = transform(transformer.transform, geom)
    return geom_utm.area / 1_000_000.0


def original_processing(raster_path: Path, geom: Polygon) -> dict:
    """Reproduce the original script logic as closely as possible."""
    with rasterio.open(raster_path) as src:
        out_image, _ = mask(
            src,
            [mapping(geom)],
            crop=True,
            filled=True,
            nodata=0,
        )

    people_per_cell = out_image[0].astype(float)

    density_100m = people_per_cell * 100.0
    masked_100m = np.ma.masked_less_equal(density_100m, 0)

    area_100_km2 = masked_100m.count() * 0.01
    total_100 = float(np.ma.masked_less_equal(people_per_cell, 0).sum())
    mean_100 = float(masked_100m.mean())
    max_100 = float(masked_100m.max())

    factor = 10
    rows, cols = people_per_cell.shape
    new_rows = rows // factor
    new_cols = cols // factor
    cropped = people_per_cell[: new_rows * factor, : new_cols * factor]
    reshaped = cropped.reshape(new_rows, factor, new_cols, factor)
    people_1km = reshaped.sum(axis=(1, 3))
    density_1km = np.ma.masked_less_equal(people_1km, 0)

    area_1km_km2 = float(density_1km.count())
    total_1km = float(density_1km.sum())
    mean_1km = float(density_1km.mean())
    max_1km = float(density_1km.max())

    return {
        "native_100m": {
            "effective_area_km2_positive_cells_only": area_100_km2,
            "total_population": total_100,
            "mean_density_people_per_km2": mean_100,
            "max_density_people_per_km2": max_100,
        },
        "aggregated_1km": {
            "effective_area_km2_occupied_blocks": area_1km_km2,
            "total_population": total_1km,
            "mean_density_people_per_km2": mean_1km,
            "max_density_people_per_km2": max_1km,
        },
        "population_lost_by_edge_cropping": total_100 - total_1km,
    }


def masked_population_sum(
    raster_path: Path,
    geom: Polygon,
    *,
    all_touched: bool,
) -> float:
    with rasterio.open(raster_path) as src:
        values, _ = mask(
            src,
            [mapping(geom)],
            crop=True,
            filled=False,
            all_touched=all_touched,
        )
        band = values[0]
        if src.nodata is not None:
            band = np.ma.masked_equal(band, src.nodata)
        band = np.ma.masked_invalid(band)
        return float(band.sum())


def exact_population_sum(raster_path: Path, geom: Polygon) -> float:
    gdf = gpd.GeoDataFrame({"id": [1]}, geometry=[geom], crs="EPSG:4326")
    result = exact_extract(
        str(raster_path),
        gdf,
        ["sum"],
        output="pandas",
    )
    return float(result.loc[0, "sum"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kml", type=Path, required=True)
    parser.add_argument("--raster", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/adjacent_population_recalculation.json"),
    )
    args = parser.parse_args()

    outer = parse_kml_polygon(
        args.kml,
        "Adjacent Area 5 km (from contingency)",
    )
    inner = parse_kml_polygon(
        args.kml,
        "Contingency Volume (from KMZ)",
    )
    adjacent = outer.difference(inner)
    if adjacent.is_empty or not adjacent.is_valid:
        raise RuntimeError("Adjacent-area ring geometry is invalid")

    area_geodesic = geodesic_area_km2(adjacent)
    area_utm = projected_area_km2(adjacent)

    original = original_processing(args.raster, adjacent)

    pop_exact = exact_population_sum(args.raster, adjacent)
    pop_centre = masked_population_sum(args.raster, adjacent, all_touched=False)
    pop_all_touched = masked_population_sum(args.raster, adjacent, all_touched=True)

    corrected = {
        "polygon_area_geodesic_km2": area_geodesic,
        "polygon_area_utm37n_km2": area_utm,
        "fractional_coverage_population": pop_exact,
        "fractional_coverage_mean_density_people_per_km2": pop_exact / area_geodesic,
        "centre_cell_population": pop_centre,
        "centre_cell_mean_density_people_per_km2": pop_centre / area_geodesic,
        "all_touched_population": pop_all_touched,
        "all_touched_mean_density_people_per_km2": pop_all_touched / area_geodesic,
    }

    diagnostics = {
        "published_target_values": {
            "mean_100m": 52.3,
            "max_100m": 900.7,
            "mean_1km": 43.5,
            "max_1km": 460.2,
        },
        "reproduction_differences": {
            "mean_100m": original["native_100m"]["mean_density_people_per_km2"] - 52.3,
            "max_100m": original["native_100m"]["max_density_people_per_km2"] - 900.7,
            "mean_1km": original["aggregated_1km"]["mean_density_people_per_km2"] - 43.5,
            "max_1km": original["aggregated_1km"]["max_density_people_per_km2"] - 460.2,
        },
        "boundary_method_spread_people": pop_all_touched - pop_centre,
        "fractional_vs_centre_people": pop_exact - pop_centre,
        "fractional_vs_all_touched_people": pop_exact - pop_all_touched,
    }

    output = {
        "input": {
            "worldpop_file": args.raster.name,
            "worldpop_units": "people per native grid cell",
            "polygon": "Adjacent Area 5 km minus Contingency Volume",
        },
        "original_processing_reproduction": original,
        "corrected_area_preserving_calculation": corrected,
        "diagnostics": diagnostics,
        "interpretation": {
            "mean_should_be_resolution_invariant": True,
            "original_means_used_different_effective_denominators": True,
            "maximum_is_expected_to_decline_with_coarser_aggregation": True,
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
