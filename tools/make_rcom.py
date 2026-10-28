"""農業集落境界(GeoPackage の rcom レイヤ)から、DEMの範囲にかかる集落を取り出して GeoJSON(緯度経度)にする。

  python tools/make_rcom.py 集落界.gpkg DEM.tif data/rcom.geojson

DEM・集落界とも平面直角座標 X系(EPSG:6678)の前提(別の系なら 6678 を変更)。
必要なもの:pip install numpy tifffile imagecodecs shapely pyproj
"""
import sqlite3, json, sys, collections, tifffile, numpy as np
from shapely import wkb
from shapely.geometry import box, mapping
from shapely.ops import transform
import warnings; warnings.filterwarnings('ignore')
from pyproj import Transformer
# DEMのデータがある範囲(6678)
p=tifffile.TiffFile(sys.argv[2]).pages[0]; z=p.asarray(); tp=p.tags[33922].value
ok=z>-9000; rows=np.where(ok.any(1))[0]; cols=np.where(ok.any(0))[0]
x0=tp[3]+cols[0]*5; x1=tp[3]+(cols[-1]+1)*5; y1=tp[4]-rows[0]*5; y0=tp[4]-(rows[-1]+1)*5
# データのある所の外形(粗いマスク)で判定
k=40; m=ok[:ok.shape[0]//k*k,:ok.shape[1]//k*k].reshape(ok.shape[0]//k,k,ok.shape[1]//k,k).any((1,3))
from shapely.geometry import Polygon
from shapely.ops import unary_union
cells=[box(tp[3]+j*k*5, tp[4]-(i+1)*k*5, tp[3]+(j+1)*k*5, tp[4]-i*k*5) for i,j in zip(*np.where(m))]
cover=unary_union(cells)
db=sqlite3.connect(sys.argv[1]); cur=db.execute('select geom,KEY,PREF_NAME,CITY_NAME,KCITY_NAME,RCOM_NAME,RCOM_KANA,HININTEI from rcom')
tr=Transformer.from_crs(6678,6668,always_xy=True)
feats=[]; cnt=collections.Counter()
for g,key,pn,cn,kn,rn,rk,hn in cur:
    flags=g[3]; env=[0,32,48,48,64][(flags>>1)&7]; geom=wkb.loads(bytes(g[8+env:]))
    if not geom.intersects(cover): continue
    inside=geom.intersection(cover).area/geom.area
    gs=geom.simplify(0.5, preserve_topology=True)
    gl=transform(lambda x,y,z=None: tr.transform(x,y), gs)
    def rnd(o): return [rnd(v) for v in o] if isinstance(o[0],(list,tuple)) else [round(o[0],7),round(o[1],7)]
    gm=mapping(gl); gm={'type':gm['type'],'coordinates':rnd(gm['coordinates'])}
    feats.append({'type':'Feature','properties':{'集落名':rn,'よみ':rk,'旧市町村':kn,'市町村':cn,'集落コード':key,'面積_ha':round(geom.area/1e4,1)},'geometry':gm})
    cnt[cn]+=1
fc={'type':'FeatureCollection','features':feats}
s=json.dumps(fc,ensure_ascii=False,separators=(',',':'))
open(sys.argv[3],'w',encoding='utf-8').write(s)
print(len(feats), cnt, round(len(s)/1e6,2),'MB')

