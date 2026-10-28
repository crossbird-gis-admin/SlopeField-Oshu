"""基盤地図情報の1mDEM(DEM1A)のダウンロードzipを、2次メッシュごとのGeoTIFFに変換する。

  python tools/fgd_dem1_tiles.py ダウンロードした.zip 出力フォルダ [EPSG] [2次メッシュ番号,...] [half]

- 最後に half を付けると、2次メッシュを南北に分けて oshu_1m_<番号>_s.tif / _n.tif にする(メモリや時間が足りないとき)

- ダウンロードzip(中に2次メッシュごとのzip)をそのまま読む。QGISのプラグインに通す必要はない
- 出力:出力フォルダ/oshu_1m_<2次メッシュ>.tif(平面直角座標・1m、標高はcm単位の整数、GDALスケール0.01付き)
- 既にある出力はとばすので、途中で止まっても同じコマンドで続きから
必要なもの:pip install numpy tifffile imagecodecs pyproj
"""
import zipfile, io, re, sys, os, time, numpy as np, tifffile
from pyproj import Transformer
SRC, OUT = sys.argv[1], sys.argv[2]
EPSG = int(sys.argv[3]) if len(sys.argv) > 3 else 6678
ONLY = set(sys.argv[4].split(',')) if len(sys.argv) > 4 and sys.argv[4] not in ('', 'all') else None
HALF = len(sys.argv) > 5 and sys.argv[5] == 'half'
R = 1.0
os.makedirs(OUT, exist_ok=True)
fw = Transformer.from_crs(6668, EPSG, always_xy=True); iv = Transformer.from_crs(EPSG, 6668, always_xy=True)
def parse(t):
    g = lambda p: re.search(p, t, re.S).group(1)
    lat0, lon0 = map(float, g(r'<gml:lowerCorner>(.*?)</gml:lowerCorner>').split())
    lat1, lon1 = map(float, g(r'<gml:upperCorner>(.*?)</gml:upperCorner>').split())
    hx, hy = map(int, g(r'<gml:high>(.*?)</gml:high>').split()); W, H = hx + 1, hy + 1
    sx, sy = map(int, g(r'<gml:startPoint>(.*?)</gml:startPoint>').split())
    body = g(r'<gml:tupleList>(.*?)</gml:tupleList>').strip()
    vals = np.array(re.findall(r',\s*(-?[0-9.]+)', body), np.float32)
    bad = np.array([k == 'データなし' for k in re.findall(r'^\s*([^,\s]+),', body, re.M)])
    vals[bad | (vals < -9000)] = np.nan
    v = np.full(W * H, np.nan, np.float32); s = sy * W + sx; v[s:s + len(vals)] = vals[:W * H - s]
    return lat0, lon0, lat1, lon1, W, H, v.reshape(H, W)
md = '<GDALMetadata><Item name="SCALE" sample="0" role="scale">0.01</Item><Item name="OFFSET" sample="0" role="offset">0</Item></GDALMetadata>'
gk = (1, 1, 0, 4, 1024, 0, 1, 1, 1025, 0, 1, 1, 3072, 0, 1, EPSG, 3076, 0, 1, 9001)
def convert(code, fn, tiles):
    t0 = time.time()
    W0 = min(t[1] for t in tiles); E0 = max(t[3] for t in tiles); S0 = min(t[0] for t in tiles); N0 = max(t[2] for t in tiles)
    CX = (tiles[0][3] - tiles[0][1]) / tiles[0][4]; CY = (tiles[0][2] - tiles[0][0]) / tiles[0][5]
    NX = round((E0 - W0) / CX); NY = round((N0 - S0) / CY)
    z = np.full((NY, NX), np.nan, np.float32)
    for lat0, lon0, lat1, lon1, W, H, v in tiles:
        c = round((lon0 - W0) / CX); r = round((N0 - lat1) / CY)
        sub = z[r:r + H, c:c + W]; k = ~np.isnan(v); sub[k] = v[k]
    del tiles
    ok = ~np.isnan(z); rr = np.where(ok.any(1))[0]; cc = np.where(ok.any(0))[0]
    lon = W0 + np.array([cc[0], cc[-1] + 1, cc[0], cc[-1] + 1]) * CX; lat = N0 - np.array([rr[0], rr[0], rr[-1] + 1, rr[-1] + 1]) * CY
    xs, ys = fw.transform(lon, lat)
    x0 = np.floor(min(xs)) - 2; x1 = np.ceil(max(xs)) + 2; y0 = np.floor(min(ys)) - 2; y1 = np.ceil(max(ys)) + 2
    W = int(x1 - x0); H = int(y1 - y0)
    out = np.full((H, W), -999999, np.int32); X = x0 + np.arange(W) + 0.5
    for r0 in range(0, H, 512):
        r1 = min(H, r0 + 512); Y = y1 - (np.arange(r0, r1) + 0.5)
        XX, YY = np.meshgrid(X, Y); lo, la = iv.transform(XX, YY)
        u = (lo - W0) / CX - 0.5; v = (N0 - la) / CY - 0.5
        c = np.floor(u).astype(np.int64); r = np.floor(v).astype(np.int64); fu = (u - c).astype(np.float32); fv = (v - r).astype(np.float32)
        inb = (c >= 0) & (r >= 0) & (c < NX - 1) & (r < NY - 1); c = np.clip(c, 0, NX - 2); r = np.clip(r, 0, NY - 2)
        val = (z[r, c] * (1 - fu) + z[r, c + 1] * fu) * (1 - fv) + (z[r + 1, c] * (1 - fu) + z[r + 1, c + 1] * fu) * fv
        nn = z[np.clip(np.rint(v).astype(np.int64), 0, NY - 1), np.clip(np.rint(u).astype(np.int64), 0, NX - 1)]
        val = np.where(np.isnan(val), nn, val); val[~inb] = np.nan
        out[r0:r1] = np.where(np.isnan(val), -999999, np.rint(val * 100)).astype(np.int32)
    del z
    tmp = fn + '.part'
    tifffile.imwrite(tmp, out, tile=(256, 256), compression='zlib', compressionargs={'level': 6}, predictor=2,
                     extratags=[(33550, 'd', 3, (R, R, 0.0)), (33922, 'd', 6, (0, 0, 0, x0, y1, 0)), (34735, 'H', len(gk), gk),
                                (42113, 's', 0, '-999999'), (42112, 's', 0, md)])
    os.replace(tmp, fn)
    print(code, f'{W}×{H}', f'データ {(out > -999999).mean():.0%}', f'{os.path.getsize(fn) / 1e6:.0f} MB', f'{time.time() - t0:.0f}秒', flush=True)

src = zipfile.ZipFile(SRC)
for info in src.infolist():
    m = re.search(r'FG-GML-(\d{6})-DEM1', info.filename)
    if not m: continue
    code = m.group(1)
    if ONLY and code not in ONLY: continue
    inner = zipfile.ZipFile(io.BytesIO(src.read(info)))
    xmls = [j for j in inner.infolist() if j.filename.lower().endswith('.xml')]
    # 3次メッシュの行番号(FG-GML-5841-53-XY の X、0=南)で南北に分ける
    row = lambda j: int(re.search(r'FG-GML-\d{4}-\d{2}-(\d)\d', j.filename).group(1))
    for part, sel in ([('_s', [j for j in xmls if row(j) < 5]), ('_n', [j for j in xmls if row(j) >= 5])] if HALF else [('', xmls)]):
      fn = os.path.join(OUT, f'oshu_1m_{code}{part}.tif')
      if not sel: continue
      if os.path.exists(fn) or (not HALF and os.path.exists(os.path.join(OUT, f'oshu_1m_{code}_s.tif'))): print(code + part, 'あり(とばす)'); continue
      convert(code + part, fn, [parse(inner.read(j).decode('utf-8')) for j in sel])
print('完了')
