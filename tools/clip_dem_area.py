"""1mDEM(GeoTIFF、何枚でも)を、区画ポリゴンの周り(既定20m)だけに切り詰め、GitHubに載せられる大きさのタイルにする。

  python tools/clip_dem_area.py 区画.gpkg data/dem1m 1mのtif ... [--buffer 20]

- 区画の外(バッファの外)はデータなしにするので、圧縮でとても小さくなる
- 出力は data/dem1m_01.tif …(1枚20MB未満、cm単位の整数、GDALスケール0.01付き)
必要なもの:pip install numpy tifffile imagecodecs shapely pillow
"""
import sys, os, sqlite3, numpy as np, tifffile
from shapely import wkb, make_valid
from shapely.ops import unary_union
from PIL import Image, ImageDraw
args = sys.argv[1:]; buf = 20.0
if '--buffer' in args: i = args.index('--buffer'); buf = float(args[i + 1]); del args[i:i + 2]
gpkg, prefix, tifs = args[0], args[1], args[2:]
db = sqlite3.connect(gpkg)
t, col = db.execute('select table_name, column_name from gpkg_geometry_columns').fetchone()
geoms = []
for (g,) in db.execute(f'select "{col}" from "{t}"'):
    flags = g[3]; env = [0, 32, 48, 48, 64][(flags >> 1) & 7]; s = wkb.loads(bytes(g[8 + env:]))
    geoms.append(make_valid(s) if not s.is_valid else s)
area = unary_union([g.buffer(buf, 4) for g in geoms]).simplify(1.0)
parts = [area] if area.geom_type == 'Polygon' else [p for p in area.geoms if p.geom_type == 'Polygon']
md = '<GDALMetadata><Item name="SCALE" sample="0" role="scale">0.01</Item><Item name="OFFSET" sample="0" role="offset">0</Item></GDALMetadata>'
os.makedirs(os.path.dirname(prefix) or '.', exist_ok=True)
n = 0; total = 0
def write(arr, x0, y1, gk, sx):
    global n, total
    ok = arr > -999999
    if not ok.any(): return
    rr = np.where(ok.any(1))[0]; cc = np.where(ok.any(0))[0]
    a = arr[rr[0]:rr[-1] + 1, cc[0]:cc[-1] + 1]; x0 += cc[0] * sx; y1 -= rr[0] * sx
    fn = None
    for _ in range(2):
        n += 1; fn = f'{prefix}_{n:02d}.tif'
        tifffile.imwrite(fn, a, tile=(256, 256), compression='zlib', compressionargs={'level': 9}, predictor=2,
                         extratags=[(33550, 'd', 3, (sx, sx, 0.0)), (33922, 'd', 6, (0, 0, 0, x0, y1, 0)), (34735, 'H', len(gk), gk),
                                    (42113, 's', 0, '-999999'), (42112, 's', 0, md)])
        mb = os.path.getsize(fn) / 1e6
        if mb < 20 or min(a.shape) < 512: break
        # 20MB以上なら、長い方の辺で2つに割る(2セル重ねる)
        os.remove(fn); n -= 1
        if a.shape[0] >= a.shape[1]:
            h = a.shape[0] // 2; write(a[:h + 2].copy(), x0, y1, gk, sx); write(a[h - 2:].copy(), x0, y1 - (h - 2) * sx, gk, sx)
        else:
            w = a.shape[1] // 2; write(a[:, :w + 2].copy(), x0, y1, gk, sx); write(a[:, w - 2:].copy(), x0 + (w - 2) * sx, y1, gk, sx)
        return
    total += mb; print(f'{fn}  {a.shape[1]}×{a.shape[0]}  {mb:.1f} MB', flush=True)
for tif in tifs:
    pg = tifffile.TiffFile(tif).pages[0]; tp = pg.tags[33922].value; sx = pg.tags[33550].value[0]
    gk = pg.tags[34735].value; x0, y1 = tp[3], tp[4]; H, W = pg.shape
    hit = [p for p in parts if p.bounds[2] > x0 and p.bounds[0] < x0 + W * sx and p.bounds[3] > y1 - H * sx and p.bounds[1] < y1]
    if not hit: continue
    z = pg.asarray()
    if z.dtype != np.int32:   # float(m)のDEMは cm の整数に
        nd = float(pg.tags[42113].value) if 42113 in pg.tags else -9999
        bad = ~np.isfinite(z) | (z == nd) | (z < -9000); z = np.where(bad, -999999, np.rint(z * 100)).astype(np.int32)
    m = Image.new('1', (W, H), 0); dr = ImageDraw.Draw(m)
    for p in hit:
        dr.polygon([((x - x0) / sx, (y1 - y) / sx) for x, y in p.exterior.coords], fill=1)
        for h in p.interiors: dr.polygon([((x - x0) / sx, (y1 - y) / sx) for x, y in h.coords], fill=0)
    z[~np.array(m)] = -999999
    write(z, x0, y1, gk, sx)
print(f'合計 {total:.0f} MB、{n} 枚')
