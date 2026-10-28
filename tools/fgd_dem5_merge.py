"""基盤地図情報DEM(XML)→ モザイク(緯度経度)。5A優先、5Bで穴埋め"""

import zipfile, re, sys, numpy as np, tifffile
from pyproj import Transformer
CELL = 0.2/3600   # 5mDEM: 0.2秒
def parse(t):
    g=lambda p: re.search(p,t,re.S).group(1)
    lat0,lon0=map(float,g(r'<gml:lowerCorner>(.*?)</gml:lowerCorner>').split())
    lat1,lon1=map(float,g(r'<gml:upperCorner>(.*?)</gml:upperCorner>').split())
    hx,hy=map(int,g(r'<gml:high>(.*?)</gml:high>').split()); W,H=hx+1,hy+1
    sx,sy=map(int,g(r'<gml:startPoint>(.*?)</gml:startPoint>').split())
    rows=g(r'<gml:tupleList>(.*?)</gml:tupleList>').strip().split('\n')
    v=np.full(W*H,np.nan,np.float32)
    vals=np.array([float(r.split(',')[1]) for r in rows],np.float32)
    kind=[r.split(',')[0] for r in rows]
    vals[np.array([k=='データなし' for k in kind])]=np.nan
    vals[vals<-9000]=np.nan
    s=sy*W+sx; v[s:s+len(vals)]=vals[:W*H-s]
    return lat0,lon0,lat1,lon1,W,H,v.reshape(H,W)
def load(z):
    out=[]; zf=zipfile.ZipFile(z)
    for n in zf.namelist():
        if n.endswith('.xml'): out.append(parse(zf.read(n).decode('utf-8')))
    return out
OUT=sys.argv[3]
A=load(sys.argv[1]); B=load(sys.argv[2])
allt=A+B
W0=min(t[1] for t in allt); E0=max(t[3] for t in allt); S0=min(t[0] for t in allt); N0=max(t[2] for t in allt)
NX=round((E0-W0)/CELL); NY=round((N0-S0)/CELL)
print('grid',NX,NY,W0,S0,E0,N0)
z=np.full((NY,NX),np.nan,np.float32); src=np.zeros((NY,NX),np.uint8)
def put(ts,code):
    for lat0,lon0,lat1,lon1,W,H,v in ts:
        c=round((lon0-W0)/CELL); r=round((N0-lat1)/CELL)
        sub=z[r:r+H,c:c+W]; m=np.isnan(sub)&~np.isnan(v)
        sub[m]=v[m]; src[r:r+H,c:c+W][m]=code
put(A,1); put(B,2)

W0,N0,C=W0,N0,CELL
fw=Transformer.from_crs(6668,6678,always_xy=True); iv=Transformer.from_crs(6678,6668,always_xy=True)
ok=~np.isnan(z); rows=np.where(ok.any(1))[0]; cols=np.where(ok.any(0))[0]
lon=W0+np.array([cols[0],cols[-1]+1,cols[0],cols[-1]+1])*C; lat=N0-np.array([rows[0],rows[0],rows[-1]+1,rows[-1]+1])*C
xs,ys=fw.transform(lon,lat); R=5.0
x0=np.floor(min(xs)/R)*R; x1=np.ceil(max(xs)/R)*R; y0=np.floor(min(ys)/R)*R; y1=np.ceil(max(ys)/R)*R
W=int((x1-x0)/R); H=int((y1-y0)/R); print('out',W,H)
out=np.full((H,W),-9999,np.float32); osrc=np.zeros((H,W),np.uint8)
X=x0+(np.arange(W)+0.5)*R
for r0 in range(0,H,256):
    r1=min(H,r0+256); Y=y1-(np.arange(r0,r1)+0.5)*R
    XX,YY=np.meshgrid(X,Y); lo,la=iv.transform(XX,YY)
    u=(lo-W0)/C-0.5; v=(N0-la)/C-0.5
    c=np.floor(u).astype(int); r=np.floor(v).astype(int); fu=u-c; fv=v-r
    inb=(c>=0)&(r>=0)&(c<NX-1)&(r<NY-1); c=np.clip(c,0,NX-2); r=np.clip(r,0,NY-2)
    a=z[r,c]; b=z[r,c+1]; d=z[r+1,c]; e=z[r+1,c+1]
    val=(a*(1-fu)+b*fu)*(1-fv)+(d*(1-fu)+e*fu)*fv
    # 端で4点そろわないときは最近傍
    rn=np.clip(np.rint(v).astype(int),0,NY-1); cn=np.clip(np.rint(u).astype(int),0,NX-1)
    nn=z[rn,cn]; val=np.where(np.isnan(val),nn,val); val[~inb]=np.nan
    out[r0:r1]=np.where(np.isnan(val),-9999,val); osrc[r0:r1]=np.where(np.isnan(val),0,src[rn,cn])
gk=(1,1,0,4, 1024,0,1,1, 1025,0,1,1, 3072,0,1,6678, 3076,0,1,9001)
tags=[(33550,'d',3,(R,R,0.0)),(33922,'d',6,(0,0,0,x0,y1,0)),(34735,'H',len(gk),gk),(42113,'s',0,'-9999')]
tifffile.imwrite(OUT+'/奥州_5m_5A5B_EPSG6678.tif',out,tile=(256,256),compression='zlib',predictor=3,extratags=tags)
tifffile.imwrite(OUT+'/奥州_5m_出典_EPSG6678.tif',osrc,tile=(256,256),compression='zlib',extratags=[t for t in tags if t[0]!=42113]+[(42113,'s',0,'0')])
v=out[out>-9000]; print('valid',v.size/out.size,'z',v.min(),v.max(),'5A',(osrc==1).sum(),'5B',(osrc==2).sum())
