#Downscaling code for GCM, developed by Mamad Tamamadin, Department of Meteorology, Institut Teknologi Bandung. The method is adopted from Lima et al (2021) with paper title "A Bayesian Kriging model applied for spatial downscaling of daily rainfall from GCMs".  
#Guide:
# Install python 
# Install the libraries: netCDF4, scikit-gstat, pykrige, matplotlib, numpy, pandas, scipy, pyshp, pyproj, pathlib, sridentify
# Create the home directory, for example: "climate_change", and create some directories under this directory, as follows:
# Create the directory "script" and store the this code file in this directory. 
# Create the directory "csv" and store the observed rainfall data in each csv files in this directory.  
# Create the directory "subdas" and store the shapefile of subdas in this directory.
# Create the directory "output". You can see the output result in this directory.
# Create the directory "original". You can see the original csv file of GCM after you finished the downscaling.
# Create the directory "GCM" and store the GCM files in this directory.
# Command for executing this code: "python downscaling_GCM_bayesian_kriging.py subbasin_name gcm_file start_year end_year resolution_in_km"
# For example: "python downscaling_GCM_bayesian_kriging.py indonesia 01b-pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20400101-20641231 2020 2030 20"

from sridentify import Sridentify
import netCDF4
from os import listdir
from os.path import isfile, join
from skgstat import Variogram
from pykrige.ok import OrdinaryKriging
from pykrige.uk import UniversalKriging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import datetime
import scipy.stats as stats 
from scipy.stats import invgamma
import shapefile
import sys
from pyproj import Proj, transform
from pathlib import Path

da_das = pd.DataFrame(columns=['id','lonmin','lonmax','latmin','latmax','dx','dy'])

dasname=sys.argv[1]
gcm=str(sys.argv[2])
year1=sys.argv[3]
year2=sys.argv[4]
res=sys.argv[5]

home='/' #locate the home directory
pathshp=home+'subdas/'+dasname+".shp"
pathprj=home+'subdas/'+dasname+".prj"
pathgcm=home+'GCM/'+gcm+".nc"
pathout=home+'output/'+dasname+"-"+year1+"-"+year2+"-"+gcm+"/"
csvf = home+'csv/'+dasname+"/"

startrun=datetime.datetime.now()

sf=shapefile.Reader(pathshp)
shapes = sf.shapes()

ident=Sridentify()
ident.from_file(pathprj)
epsg=str(ident.get_epsg())

if (epsg == '32748'):
    inProj = Proj(init='epsg:32748')
elif (epsg == '4326'): 
    inProj = Proj(init='epsg:4326')
outProj = Proj(init='epsg:4326')

dom_all=sf.bbox
xmin,ymin = transform(inProj,outProj,dom_all[0],dom_all[1])
xmax,ymax = transform(inProj,outProj,dom_all[2],dom_all[3])
xmin=xmin-0.05
ymin=ymin-0.05
xmax=xmax+0.05
ymax=ymax+0.05
print(xmin,xmax)
print(ymin,ymax)

da=pd.DataFrame(columns=['station','lon','lat','alpha','beta'])

onlyfiles = [f for f in listdir(csvf) if isfile(join(csvf, f))]
nsheet=len(onlyfiles)

for i,sh in enumerate(onlyfiles):
    stasiun=sh.replace('.csv','')
    df = pd.read_csv(csvf+sh)
    lon = df['lon'].tolist()
    lon1 = float(lon[0])
    lat = df['lat'].tolist()
    lat1 = float(lat[0])
    df['ch']=df['ch'].astype(str)
    df['ch']=df['ch'].replace(np.nan, 0)
    df['ch']=df['ch'].replace('nan', 0)
    data = df['ch'].tolist()
    data=[0 if i==' ' else i for i in data]
    data=[0 if i=='' else i for i in data]
    data=[0 if i=="-" else i for i in data]
    data=[0 if i==" -" else i for i in data]
    data=[0 if i=="*-" else i for i in data]
    data=[0 if i==" " else i for i in data]
    data=[0 if i=="TH" else i for i in data]
    data=[0 if i=="TAD" else i for i in data]
    data=[i if i else 0 for i in data]    
    data=[x for x in data if str(x) != 'nan']
        

    dataa=[]
    for i,ch in enumerate(data):
        if ch=="":
            dataa.append(0)
        elif ch=='':
            dataa.append(0)
        elif ch==None:
            dataa.append(0)        
        else:
            try:
                dataa.append(float(ch))
            except:
                continue
    
    data=dataa

    data = [i for i in data if not math.isnan(i)]
    data2=[float(i)**2 for i in data]
    
    data2 = [i for i in data2 if not math.isnan(i)]  
    mean_of_distribution = np.mean(data)
    variance_of_distribution = np.var(data)
    g2=np.mean(data2)
    
    def gamma(g1, g2, size):
        g_alpha = g1**2/(g2-g1**2)
        g_beta = (g2-g1**2)/g1
        #print(g_alpha, g_beta) 
        da.loc[i]=[stasiun]+[lon1]+[lat1]+[g_alpha]+[g_beta]       
    gamma(mean_of_distribution,g2,len(data))
da.dropna(subset = ["alpha"], inplace=True)
da.dropna(subset = ["beta"], inplace=True)
print(da)   


#calculating the number of grid to get resolution of 10 km
dx=xmax-xmin
dy=ymax-ymin
g10km=float(float(res)/111)
ngridx=int(dx/g10km)
ngridy=int(dy/g10km)

grid_lon=np.linspace(xmin,xmax,ngridx) #grid size
grid_lat=np.linspace(ymin,ymax,ngridy) #grid size
model = 'exponential'

#----------krigging alpha--------------------
V2= Variogram(da[['lon','lat']].values,da.alpha.values, normalize=False)
V2.maxlag = 0.3
V2.n_lags = 25
V2.plot(show=False)
V2.model = model
ok2 = OrdinaryKriging(da["lon"], da["lat"],da["alpha"],variogram_model=model, verbose=True, nlags=30)
field2,ss1 = ok2.execute('grid',grid_lon,grid_lat)
xx1,yy1 = np.meshgrid(grid_lon,grid_lat)

xxs=xx1.tolist()
yys=yy1.tolist()
alphas=field2.tolist()
print(len(xxs),len(yys),len(alphas))
da_a_krg = pd.DataFrame(columns=['idgrid','lon','lat','alpha'])
i=0
for m in range(0, len(yys)):
    for n in range(0, len(xxs[0])):    
        da_a_krg.loc[i]=[i+1]+[xxs[m][n]]+[yys[m][n]]+[alphas[m][n]]
        i=i+1
print(da_a_krg)

#----------krigging beta--------------------
V2= Variogram(da[['lon','lat']].values,da.beta.values, normalize=False)
V2.maxlag = 0.3
V2.n_lags = 25
V2.plot(show=False)
V2.model = model
ok2 = OrdinaryKriging(da["lon"], da["lat"],da["beta"],variogram_model=model, verbose=True, nlags=30)
field2,ss1 = ok2.execute('grid',grid_lon,grid_lat)
xx1,yy1 = np.meshgrid(grid_lon,grid_lat)

xxs=xx1.tolist()
yys=yy1.tolist()
betas=field2.tolist()
da_b_krg = pd.DataFrame(columns=['idgrid','lon','lat','beta'])
i=0
for m in range(0, len(yys)):
    for n in range(0, len(xxs[0])):    
        da_b_krg.loc[i]=[i+1]+[xxs[m][n]]+[yys[m][n]]+[betas[m][n]]
        i=i+1

print(da_b_krg)
#------------------GCM parameter alpha & beta --------------------------
da_gcm=pd.DataFrame(columns=['lon','lat','alpha','beta'])
def gamma_process(k,lon,lat,prec_list):
    mean_of_distribution = np.mean(prec_list)
    variance_of_distribution = np.var(prec_list)
    data2=[float(i)**2 for i in prec_list]
    g2=np.mean(data2)
    def gamma(g1, g2, size):
        g_alpha = g1**2/(g2-g1**2)
        g_beta = (g2-g1**2)/g1
        da_gcm.loc[k]=[lon]+[lat]+[g_alpha]+[g_beta]
    gamma(mean_of_distribution,g2,len(prec_list))

df = pd.DataFrame(columns=['date','lon','lat','prec'])
nc = netCDF4.Dataset(pathgcm)
nc.variables.keys()

lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
time_var = nc.variables['time']
dtime = netCDF4.num2date(time_var[:],time_var.units)
def near(array,value):
        idx=(abs(array-value)).argmin()
        return idx

start = datetime.datetime(int(year1),1,1,0,0,0)
stop = datetime.datetime(int(year2),3,31,0,0,0)

istart = netCDF4.date2index(start,time_var,select='nearest')+1 
istop = netCDF4.date2index(stop,time_var,select='nearest') 

ixa = near(lon, xmin)
iya = near(lat, ymin)
ixb = near(lon, xmax) 
iyb = near(lat, ymax) 

nx = ixb-ixa
ny = iyb-iya

var = nc.variables['pr']
tim = dtime[istart:istop]

i=0; k=0
for m in range(0, nx+1):
    for n in range(0, ny+1):
        pr1 = var[istart:istop,iya+n,ixa+m]
        for t in range(0,len(pr1)):
            df.loc[i]=[str(tim[t]).split(' ')[0]]+[lon[ixa+m]]+[lat[iya+n]]+[float(pr1[t])*86400]
            i=i+1
        #print(df)
        gamma_process(k,lon[ixa+m],lat[iya+n],df["prec"])
        k=k+1
print(da_gcm)
df.to_csv(home+'original/'+dasname+"-"+year1+"-"+year2+"-"+gcm+'.csv')   

#----------- The Prediction of GCM----------------
print(' ')
print('Calculating the projection.......')
print(gcm)
print('Time range: ' + str(year1) + ' - ' + str(year2))
timstart = datetime.datetime(int(year1),1,1,0,0,0)
timstop = datetime.datetime(int(year2),12,31,0,0,0)
istart1 = netCDF4.date2index(timstart,time_var,select='nearest')
istop1 = netCDF4.date2index(timstop,time_var,select='nearest')+1
var = nc.variables['pr']
tim = dtime[istart1:istop1]
zl_list=[]
tima=[]
lona1=[]
lata1=[]
idgridl=[]
na = len(da_a_krg["alpha"])
m=0
for i in range(0,na):
    idgrid=da_a_krg["idgrid"][i]
    lona=da_a_krg["lon"][i]
    lata=da_a_krg["lat"][i]
    akrg=da_a_krg["alpha"][i]
    bkrg=da_b_krg["beta"][i]
    ix = near(da_gcm['lon'], lona)
    iy = near(da_gcm['lat'], lata)
    lon_gcm = da_gcm['lon'][ix]
    lat_gcm = da_gcm['lat'][iy]
    for latG,lonG,agcm,bgcm in zip(da_gcm['lat'],da_gcm['lon'],da_gcm['alpha'],da_gcm['beta']) :
        if (lonG==lon_gcm and latG==lat_gcm):
            ixa = near(lon, lonG)
            iya = near(lat, latG)
            zgcm = var[istart1:istop1,iya,ixa]
            def cal(zgcm,t,agcm,bgcm,akrg,bkrg):
                zr = stats.gamma.cdf((zgcm*86400), a=agcm, scale=(bgcm))
                zl = stats.gamma.ppf(zr, a=akrg, scale=(bkrg)) 
                print(zl)
                return zl
            zl_lista=[stats.gamma.ppf(stats.gamma.cdf((zgcm[t]*86400), a=agcm, scale=(bgcm)), a=akrg, scale=(bkrg)) for t in range(0,len(zgcm))]
            timaa=[str(tim[t]).split(' ')[0] for t in range(0,len(zgcm))]

            zl_list.extend(zl_lista)
            tima.extend(timaa)
            idgrida=[idgrid for i in range(0,len(zgcm))]
            idgridl.extend(idgrida)
            lonaa=[lona for i in range(0,len(zgcm))]
            lona1.extend(lonaa)
            lataa=[lata for i in range(0,len(zgcm))]
            lata1.extend(lataa)

            m=m+1
da_pred=pd.DataFrame({"idgrid":idgridl,'date':tima,'lon':lona1,'lat':lata1,'prec':zl_list})
print(da_pred)
print('Creating output files....')
da_pred.to_csv(home+'output/'+dasname+"-"+year1+"-"+year2+"-"+gcm+'-Res_'+res+'.csv')

print(' ')

nc.close()

print('Done.')
endrun=datetime.datetime.now()
runtime=endrun-startrun
print("runtime : "+str(runtime))
