"""区画ポリゴン(GeoPackage、平面直角座標)を、アプリで読む GeoJSON(緯度経度)にする。

  python tools/make_area.py 区画.gpkg data/area.geojson [EPSG]

- 壊れたポリゴン(自己交差など)は直す。座標は小数7桁(約1cm)
必要なもの:pip install shapely pyproj
"""
import sys, sqlite3, json
from shapely import wkb, make_valid
from shapely.geometry import mapping
from shapely.ops import transform
from pyproj import Transformer
src, out = sys.argv[1], sys.argv[2]
db = sqlite3.connect(src)
t, col, srs = db.execute('select table_name, column_name, srs_id from gpkg_geometry_columns').fetchone()
epsg = int(sys.argv[3]) if len(sys.argv) > 3 else srs
tr = Transformer.from_crs(epsg, 6668, always_xy=True)
cols = [r[1] for r in db.execute(f'pragma table_info("{t}")') if r[1] != col]
feats = []; fixed = 0
for row in db.execute(f'select {", ".join(chr(34) + c + chr(34) for c in cols)}, "{col}" from "{t}"'):
    g = row[-1]; flags = g[3]; env = [0, 32, 48, 48, 64][(flags >> 1) & 7]; s = wkb.loads(bytes(g[8 + env:]))
    if not s.is_valid:
        s = make_valid(s); fixed += 1
        if s.geom_type == 'GeometryCollection':
            from shapely.geometry import MultiPolygon
            s = MultiPolygon([p for gg in s.geoms for p in (gg.geoms if gg.geom_type == 'MultiPolygon' else [gg]) if p.geom_type == 'Polygon'])
    p = {c: (round(v, 3) if isinstance(v, float) else v) for c, v in zip(cols, row[:-1])}
    p['区画面積_m2'] = round(s.area, 1)
    gl = transform(lambda x, y, z=None: tr.transform(x, y), s)
    def rnd(o): return [rnd(v) for v in o] if isinstance(o[0], (list, tuple)) else [round(o[0], 7), round(o[1], 7)]
    gm = mapping(gl); feats.append({'type': 'Feature', 'properties': p, 'geometry': {'type': gm['type'], 'coordinates': rnd(gm['coordinates'])}})
json.dump({'type': 'FeatureCollection', 'features': feats}, open(out, 'w', encoding='utf-8'), ensure_ascii=False, separators=(',', ':'))
print(len(feats), '区画', '直したポリゴン', fixed)
