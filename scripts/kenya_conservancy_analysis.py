"""Reproduce the WorldPop comparison in dissertation Section 4.3."""

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from xml.etree import ElementTree

from affine import Affine
import numpy as np
from pyproj import Transformer
import rasterio
from rasterio.features import geometry_mask, geometry_window
from rasterio.windows import Window
from scipy.stats import spearmanr
from shapely.geometry import Polygon, shape
from shapely.ops import transform, unary_union
from shapely.validation import make_valid

ROOT = Path(__file__).resolve().parents[1]
KWCA_QUERY = (
    "https://services.arcgis.com/RgR03SugzTGUBkDp/arcgis/rest/services/"
    "KWCA_Member_Conservancies_2024/FeatureServer/0/query"
)
KWCA_PARAMS = {"where": "1=1", "f": "geojson", "outFields": "*", "outSR": 4326}
STATUSES = {"Operational", "Mature", "Emerging"}
NATIVE_AREA_KM2 = 0.01
OPC_NAME = "Ol Pejeta Conservancy"


def get_kwca(path):
    if path:
        content = path.read_bytes()
    else:
        with urlopen(f"{KWCA_QUERY}?{urlencode(KWCA_PARAMS)}", timeout=90) as response:
            content = response.read()
    data = json.loads(content)
    if data.get("type") != "FeatureCollection" or data.get("exceededTransferLimit"):
        raise ValueError("Expected a complete KWCA GeoJSON FeatureCollection.")
    if data.get("crs", {}).get("properties", {}).get("name") not in (
        None, "urn:ogc:def:crs:OGC:1.3:CRS84", "urn:ogc:def:crs:EPSG::4326",
        "EPSG:4326",
    ):
        raise ValueError("The KWCA polygons must be in WGS 84 longitude/latitude.")
    return data["features"], hashlib.sha256(content).hexdigest()


def read_opc(path):
    kml = ElementTree.parse(path)
    polygons = kml.findall(".//{*}Polygon")
    if len(polygons) != 1:
        raise ValueError("The thesis OPC KML must contain one polygon.")
    polygon = polygons[0]

    def ring(node):
        return [tuple(map(float, point.split(",")[:2])) for point in node.text.split()]

    outer = ring(polygon.find("./{*}outerBoundaryIs/{*}LinearRing/{*}coordinates"))
    holes = [ring(node) for node in polygon.findall(
        "./{*}innerBoundaryIs/{*}LinearRing/{*}coordinates"
    )]
    return Polygon(outer, holes)


def polygon_parts(geom):
    if geom.geom_type == "Polygon":
        return [geom]
    if hasattr(geom, "geoms"):
        return [part for child in geom.geoms for part in polygon_parts(child)]
    return []


def projected(geom):
    centre = geom.centroid
    epsg = (32600 if centre.y >= 0 else 32700) + int((centre.x + 180) // 6) + 1
    forward = Transformer.from_crs(4326, epsg, always_xy=True).transform
    backward = Transformer.from_crs(epsg, 4326, always_xy=True).transform
    return transform(forward, geom), backward


def select_sites(features, opc_path):
    selected = {}
    excluded = Counter()
    for feature in features:
        info = feature["properties"]
        if info["status"] not in STATUSES:
            excluded[info["status"]] += 1
            continue
        name = info["offic_name"].strip()
        geom = shape(feature["geometry"]) if feature["geometry"] else Polygon()
        repaired = not geom.is_valid
        if repaired:
            geom = unary_union(polygon_parts(make_valid(geom)))
        if geom.is_empty or geom.geom_type not in ("Polygon", "MultiPolygon"):
            excluded["no polygon"] += 1
            continue
        metric, backward = projected(geom)
        if metric.area / 1e6 < 50:
            excluded["under 50 km²"] += 1
            continue
        identifier = str(info["Kwca_No"])
        if identifier in selected:
            raise ValueError(f"Duplicate KWCA conservancy identifier {identifier}.")
        source = "KWCA"
        source_id = identifier
        if name == OPC_NAME:
            geom = read_opc(opc_path)
            if not geom.is_valid or geom.is_empty:
                raise ValueError("The thesis OPC polygon is invalid or empty.")
            metric, backward = projected(geom)
            source, source_id, repaired = "OPC thesis boundary", opc_path.name, False
        selected[identifier] = (
            geom, metric, backward,
            dict(name=name, category=info["consv_type"].lower(), source_status=info["status"],
                 area_km2=metric.area / 1e6, kwca_id=identifier, boundary_source=source,
                 source_id=source_id, boundary_repaired=repaired),
        )
    if sum(site[3]["name"] == OPC_NAME for site in selected.values()) != 1:
        raise ValueError("The selection must contain one Ol Pejeta Conservancy record.")
    return list(selected.values()), excluded


def mask(shape, grid, dimensions):
    if shape.is_empty:
        return np.zeros(dimensions, dtype=bool)
    return geometry_mask([shape], out_shape=dimensions, transform=grid, invert=True)


def block_sum(values, factor):
    height, width = values.shape
    return values.reshape(height // factor, factor, width // factor, factor).sum(
        axis=(1, 3), dtype=np.float64
    )


def calculate(src, site):
    polygon, projected_polygon, backward, info = site
    envelope = transform(backward, projected_polygon.buffer(5000, quad_segs=32))
    core = transform(backward, projected_polygon.buffer(-1000, quad_segs=32))

    # Align windows to the native raster origin for both 2 × 2 and 5 × 5 blocks.
    crop = geometry_window(src, [envelope], boundless=True)
    left, top = (int(offset) // 10 * 10 for offset in (crop.col_off, crop.row_off))
    right, bottom = ((int(end) + 9) // 10 * 10 for end in
                     (crop.col_off + crop.width, crop.row_off + crop.height))
    window = Window(left, top, right - left, bottom - top)
    grid = src.window_transform(window)
    raster = src.read(1, window=window, masked=True, boundless=True)
    values = raster.data
    valid = ~np.ma.getmaskarray(raster) & np.isfinite(values) & (values >= 0)
    inside = mask(polygon, grid, values.shape)
    outside = mask(envelope, grid, values.shape) & ~inside
    inner_valid = inside & valid
    outer_valid = outside & valid
    inner = np.where(inner_valid, values, 0)
    n_inside = int(inner_valid.sum())
    n_outside = int(outer_valid.sum())
    total = inner.sum(dtype=np.float64)
    record = {
        **info,
        "native_inside_valid_cells": n_inside,
        "native_outside_valid_cells": n_outside,
        "inside_raster_coverage_pct": 100 * n_inside / inside.sum() if inside.any() else np.nan,
        "outside_ring_raster_coverage_pct": 100 * n_outside / outside.sum() if outside.any() else np.nan,
        "total_raster_population": total if n_inside else np.nan,
    }
    mean_inside = total / (n_inside * NATIVE_AREA_KM2) if n_inside else np.nan
    mean_outside = (
        values[outer_valid].sum(dtype=np.float64) / (n_outside * NATIVE_AREA_KM2)
        if n_outside else np.nan
    )
    for factor in (1, 2, 5):
        population = block_sum(inner, factor)
        support = block_sum(inner_valid, factor)
        represented = support > 0
        density = np.divide(population * (1 / NATIVE_AREA_KM2), support,
                            out=np.full(support.shape, np.nan), where=represented)
        is_core = mask(core, grid * Affine.scale(factor), support.shape) & represented
        full = support == factor ** 2
        peak = density[represented].max() if n_inside else np.nan
        peaks = density == peak
        edge_peaks = peaks & ~is_core
        if not n_inside:
            location = ""
            median, p95, p99 = np.nan, np.nan, np.nan
        else:
            median, p95, p99 = np.percentile(density[represented], [50, 95, 99])
            if np.all(~peaks | ~is_core):
                location = "edge"
            elif not edge_peaks.any():
                location = "core"
            else:
                location = "tied across edge and core"
        stats = dict(mean_p_km2=mean_inside, median_p_km2=median, p95_p_km2=p95,
                     p99_p_km2=p99, max_p_km2=peak, max_location=location,
                     max_cells=int(peaks.sum()), max_cells_within_1km=int(edge_peaks.sum()),
                     max_core_p_km2=density[is_core].max() if is_core.any() else np.nan,
                     max_full_block_p_km2=density[full].max() if full.any() else np.nan,
                     mean_5km_outside_p_km2=mean_outside,
                     pop_within_1km_pct=100 * population[~is_core].sum() / total if total else np.nan,
                     represented_area_within_1km_pct=(
                         100 * support[~is_core].sum() / n_inside if n_inside else np.nan
                     ))
        record.update({f"{100 * factor}m_{key}": value for key, value in stats.items()})
    return record


def summarise(records):
    mean, p99, peak, outside, peak_100, peak_500 = (
        np.array([row[key] for row in records]) for key in
        ("200m_mean_p_km2", "200m_p99_p_km2", "200m_max_p_km2", "200m_mean_5km_outside_p_km2",
         "100m_max_p_km2", "500m_max_p_km2")
    )
    paired = np.isfinite(mean) & np.isfinite(outside)
    rho = spearmanr(mean[paired], outside[paired])[0] if paired.sum() > 1 else np.nan
    peak_ratio = peak[np.isfinite(p99) & (p99 > 0)] / p99[np.isfinite(p99) & (p99 > 0)]
    lowest_outside = paired & (outside <= np.quantile(outside[paired], .25)) if paired.any() else paired
    return dict(
        conservancies=len(records), paired_inside_outside=int(paired.sum()),
        stats_200m={label: dict(n=int(np.isfinite(values).sum()),
                                below_50=int((values < 50).sum()),
                                median_people_km2=finite_json(np.nanmedian(values)))
                    for label, values in (("mean", mean), ("p99", p99), ("maximum", peak))},
        median_peak_to_p99=finite_json(np.median(peak_ratio) if peak_ratio.size else np.nan),
        count_peak_above_twice_p99=int((peak_ratio > 2).sum()),
        count_max_drop_gt_one_third_100m_to_500m=int((peak_500 < peak_100 * (2 / 3)).sum()),
        count_200m_max_within_1km=int(sum(row["200m_max_location"] == "edge" for row in records)),
        median_mean_lowest_outside_quartile=finite_json(
            np.median(mean[lowest_outside]) if lowest_outside.any() else np.nan
        ),
        spearman_mean_inside_outside=finite_json(rho),
    )


def finite_json(value):
    return float(value) if np.isfinite(value) else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raster", type=Path, default=ROOT / "data/ken_ppp_2020.tif")
    parser.add_argument("--kwca", type=Path, help="Local KWCA GeoJSON snapshot. Default fetches the source layer.")
    parser.add_argument("--opc-kml", type=Path, default=ROOT / "data/operational_area.kml.kml")
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs/kenya_conservancies")
    args = parser.parse_args()
    features, boundary_hash = get_kwca(args.kwca)
    sites, excluded = select_sites(features, args.opc_kml)
    with rasterio.open(args.raster) as src:
        grid = src.transform
        if (src.crs is None or src.crs.to_epsg() != 4326 or grid.b or grid.d
                or not np.isclose(grid.a, 1 / 1200) or not np.isclose(grid.e, -1 / 1200)):
            raise ValueError("Expected the WGS 84, north-up, 3 arc-second WorldPop counts raster.")
        records = [calculate(src, site) for site in sites]
    records.sort(key=lambda row: (row["category"], row["name"].casefold()))
    report = dict(kwca_features=len(features), kwca_geojson_sha256=boundary_hash,
                  excluded=dict(excluded), **summarise(records))
    args.outdir.mkdir(parents=True, exist_ok=True)
    with (args.outdir / "kenya_conservancy_stats.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=records[0].keys())
        writer.writeheader()
        for record in records:
            writer.writerow({key: value if not isinstance(value, (float, np.floating)) or np.isfinite(value)
                             else "" for key, value in record.items()})
    summary = json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)
    (args.outdir / "kenya_conservancy_summary.json").write_text(summary + "\n", encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
