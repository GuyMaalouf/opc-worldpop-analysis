#!/usr/bin/env python3
"""Calculate the adjacent-area mean from WorldPop's official 1 km density raster."""

from __future__ import annotations

import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import geopandas as gpd
from exactextract import exact_extract
from pyproj import Geod
from shapely.geometry import Polygon

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
            raise RuntimeError(f"No coordinates found for {placemark_name!r}")
        coords = []
        for token in coord_node.text.split():
            lon, lat, *_ = token.split(",")
            coords.append((float(lon), float(lat)))
        geom = Polygon(coords)
        if not geom.is_valid:
            geom = geom.buffer(0)
        return geom
    raise RuntimeError(f"Placemark not found: {placemark_name!r}")


def geodesic_area_km2(geom: Polygon) -> float:
    area_m2, _ = Geod(ellps="WGS84").geometry_area_perimeter(geom)
    return abs(area_m2) / 1_000_000.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kml", type=Path, required=True)
    parser.add_argument("--raster", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    outer = parse_kml_polygon(args.kml, "Adjacent Area 5 km (from contingency)")
    inner = parse_kml_polygon(args.kml, "Contingency Volume (from KMZ)")
    adjacent = outer.difference(inner)
    area_km2 = geodesic_area_km2(adjacent)

    gdf = gpd.GeoDataFrame({"id": [1]}, geometry=[adjacent], crs="EPSG:4326")
    result = exact_extract(
        str(args.raster),
        gdf,
        [
            "coverage_weighted_mean=mean",
            "population=sum(coverage_weight=area_spherical_km2)",
        ],
        output="pandas",
    )

    mean_density = float(result.loc[0, "coverage_weighted_mean"])
    total_population = float(result.loc[0, "population"])
    output = {
        "source": "WorldPop Kenya 2020 unconstrained population-density raster, official 1 km product",
        "polygon_area_geodesic_km2": area_km2,
        "coverage_weighted_mean_density_people_per_km2": mean_density,
        "area_weighted_population_people": total_population,
        "population_divided_by_polygon_area_people_per_km2": total_population / area_km2,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
