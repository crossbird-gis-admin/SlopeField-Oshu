"""DEM(GeoTIFF)を、GitHubのWebアップロード(1ファイル25MBまで)に収まるタイルに分ける。

  python tools/dem_tiles.py 入力.tif data/名前 [列数] [行数]

- 標高を cm 単位の整数(int32)にして圧縮する。GDALのスケール(0.01)を書くので、QGISでもmで表示される
- タイルは2セルずつ重ねる(境目でも傾斜が計算できるように)
- データの無いタイルは作らない。各タイルはデータのある範囲に切り詰める
必要なもの:pip install numpy tifffile imagecodecs
"""
import sys, os, numpy as np, tifffile
src, prefix = sys.argv[1], sys.argv[2]
NC = int(sys.argv[3]) if len(sys.argv) > 3 else 3
NR = int(sys.argv[4]) if len(sys.argv) > 4 else 2
OV = 2
p = tifffile.TiffFile(src).pages[0]
z = p.asarray().astype(np.float64)
sx, sy, _ = p.tags[33550].value; tp = p.tags[33922].value; gk = p.tags[34735].value
nd = p.tags[42113].value if 42113 in p.tags else None
ok = np.isfinite(z) & (z > -9000)
if nd not in (None, ''): ok &= z != float(nd)
H, W = z.shape
cm = np.where(ok, np.rint(z * 100), -999999).astype(np.int32)
md = '<GDALMetadata><Item name="SCALE" sample="0" role="scale">0.01</Item><Item name="OFFSET" sample="0" role="offset">0</Item></GDALMetadata>'
os.makedirs(os.path.dirname(prefix) or '.', exist_ok=True)
out = []
for j in range(NR):
    for i in range(NC):
        r0, r1 = max(0, H * j // NR - OV), min(H, H * (j + 1) // NR + OV)
        c0, c1 = max(0, W * i // NC - OV), min(W, W * (i + 1) // NC + OV)
        m = ok[r0:r1, c0:c1]
        if not m.any(): continue
        rr = np.where(m.any(1))[0]; cc = np.where(m.any(0))[0]
        a0, a1, b0, b1 = r0 + rr[0], r0 + rr[-1] + 1, c0 + cc[0], c0 + cc[-1] + 1
        fn = f'{prefix}_r{j + 1}c{i + 1}.tif'
        tifffile.imwrite(fn, cm[a0:a1, b0:b1], tile=(256, 256), compression='zlib', compressionargs={'level': 9}, predictor=2,
                         extratags=[(33550, 'd', 3, (sx, sy, 0.0)), (33922, 'd', 6, (0, 0, 0, tp[3] + b0 * sx, tp[4] - a0 * sy, 0)),
                                    (34735, 'H', len(gk), gk), (42113, 's', 0, '-999999'), (42112, 's', 0, md)])
        mb = os.path.getsize(fn) / 1e6; out.append((fn, mb))
        print(f'{fn}  {b1 - b0}×{a1 - a0}  {mb:.1f} MB' + ('  ※25MBを超えています。列数・行数を増やしてください' if mb > 24.5 else ''))
print('CONFIG.dems に書く行:')
for fn, _ in out: print(f"    {{ url: '{fn.replace(os.sep, '/')}', name: '{os.path.basename(fn)[:-4]}' }},")
