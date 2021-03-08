from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from gal_data import gal_data
from astropy import units as u
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
from astropy.modeling.models import Sersic1D
import warnings
warnings.filterwarnings("ignore")

xnew = np.arange(0.25,2.0,.01)
radii = [.25,.5,.75,1.0,1.25,1.5,2.0]

def func(x,m,b):
    return(m*x+b)

def estimateHalfLight(df,galnames,b):
    if galnames.dtype == 'int':
        inds = df.index.values
    else:
        galnames = [galnames] if isinstance(galnames , str) else galnames
        pgcs =gal_data(galnames)['PGC']
        inds = [df[df.PGC==i].index[0] for i in pgcs]
    #bands = ['w1','w2','w3','w4','nuv','fuv']

    pgcs = []
    halflights = []
    nanmasks = []
    whys=[]
    bands = []
    
    for i in np.arange(len(inds)):
        galname=galnames[i]
        try:
            xnew,ypred,numnans,why = fitLogFunc(df,inds[i],b)
            if type(ypred) == float:
                halflights.append(np.nan)
                pgcs.append(galname)
                nanmasks.append(0)
                whys.append(why)
            else:
                
                nanmask = np.isnan(ypred)
                lennans = len(nanmask[nanmask==True])
                
                cs=np.cumsum(ypred[~nanmask])
                area = cs[-1]
                halflightloc = np.where(cs>(area/2.0))[0][0]
                
                halflight = xnew[halflightloc+lennans]*1.0
                halflights.append(halflight)
                pgcs.append(galname)
                nanmasks.append(numnans)
                whys.append(why)
        except:
            halflights.append(np.nan)
            pgcs.append(galname)
            nanmasks.append(np.nan)
            whys.append('Something Else')
            
    return(halflights,pgcs,nanmasks,whys)
        
def fitLogFuncAp(df,ind,b):
    try:
        dist =df.loc[ind].DIST_MPC
        r25 = df.loc[ind].R25_DEG
        radii_arcsec = np.asarray([r*r25*3600 for r in radii])
        
        suffixes=['0p25R25','0p5R25','0p75R25','1p0R25','1p25R25','1p5R25','2p0R25']
        cols = [b.upper()+'_'+s for s in suffixes]
        yval = df.loc[ind][cols].values.astype(float)

        areas = (np.pi*((radii_arcsec*u.arcsec.to(u.radian))*dist*1e6)**2)
        
        sb = yval/(areas)
        nanmask = np.isnan(yval)
        radii_arcsec = radii_arcsec[~nanmask]
        yval = yval[~nanmask]

        numnans = len(nanmask[nanmask==True])
        
        if len(yval)==0:
            return(np.nan,np.nan,numnans,'All NaN')
        else:
            try:
                popt, pcov = curve_fit(func, np.log10(radii_arcsec),np.log10(yval))
                xnew = np.logspace(1e-3,np.log10(radii_arcsec[-1]),500)
                ypred = func(xnew,*popt)
                return(xnew,ypred,numnans,'')
                
            except:
                return(np.nan,np.nan,numnans,'Fit wouldnt work')
    except:
        return(np.nan,np.nan,numnans,'Something weird')
        
def expfunc(x, a, b, c):
    return a*np.exp(-b*x) + c
    
def logfunc(x, a, b, c): # x-shifted log
    return a*np.log(x + b)+c
        
def fitExponentialFunc(rp):
    norm_x = rp.r_arcsec.min()
    norm_y = rp.I.max()
    fx2 = rp.r_arcsec - norm_x + 1
    fy2 = rp.I/norm_y
    
    popt, pcov = curve_fit(expfunc, fx2,fy2, p0=(1,1,1), maxfev=6000)
    
    mask = fy2>expfunc(fx2,*popt)
    
    rp = rp[~mask]
    
    return(rp,norm_x,norm_y,fx2,fy2)
    
def fitLogFunc(rp):
    norm_x = rp.r_arcsec.min()
    norm_y = np.cumsum(rp.I).max()
    fx2 = rp.r_arcsec - norm_x + 1
    fy2 = np.cumsum(rp.I)/norm_y
    
    popt, pcov = curve_fit(logfunc, fx2,fy2, p0=(1,1,1), maxfev=6000)
    
    return(rp,norm_x,norm_y,fx2,fy2,popt,pcov)
        
def findHalfLightSersicdf(pgcs,df1,df2):
    halflights = []
    amps = []
    ns = []
    datamasks = []
    nummasked = []
    #print('0')
    for i in np.arange(len(pgcs)):
        print(len(pgcs)-i)
        galmask = df2.gal.isin([pgcs[i]])
        rp = df2.loc[galmask]
        if np.isnan(rp.r_arcsec).all()==True:
            halflights.append('All NaN')
            amps.append(np.nan)
            ns.append(np.nan)
            datamasks.append(np.nan)
            nummasked.append(np.nan)
            #print('1')
        elif rp.r_arcsec.min()>200:
            halflights.append('Bad Decomp')
            amps.append(np.nan)
            ns.append(np.nan)
            datamasks.append(np.nan)
            nummasked.append(np.nan)
            #print('2')
        else:
            try:
                r25 = df1.loc[i].R25_DEG*3600
                mask = rp.r_arcsec<2*r25
                rp = rp[mask]
                #mask = rp.I>0
                #rp = rp[mask]
                
                sersic = models.Sersic1D()
                outlier_fit = fitting.FittingWithOutlierRemoval(fitting.LevMarLSQFitter(),sigma_clip, niter=3, sigma=2.0)
                fitted_model,filtered_data = outlier_fit(sersic,rp.r_arcsec,rp.I)
                
                    
                halflights.append(fitted_model.r_eff.value)
                amps.append(fitted_model.amplitude.value)
                ns.append(fitted_model.n.value)
                datamasks.append(filtered_data)
                nummasked.append(len(filtered_data==True))
            except:
                halflights.append('fit fail')
                amps.append('fit fail')
                ns.append('fit fail')
                datamasks.append('fit fail')
                nummasked.append('fit fail')

                
    return(pgcs,halflights,amps,ns,datamasks,nummasked)
        
def findHalfLightdf(pgcs,df1,df2):
    halflights = []
    #print('0')
    for i in np.arange(len(pgcs)):
        print(len(pgcs)-i)
        galmask = df2.gal.isin([pgcs[i]])
        rp = df2.loc[galmask]
        if np.isnan(rp.r_arcsec).all()==True:
            halflights.append('All NaN')
            #print('1')
        elif rp.r_arcsec.min()>200:
            halflights.append('Bad Decomp')
            #print('2')
        else:
            try:
                r25 = df1.loc[i].R25_DEG*3600
                mask = rp.r_arcsec<2*r25
                rp = rp[mask]
                
                rp_clean,norm_x,norm_y,fx2,fy2 = fitExponentialFunc(rp)
                #print('3')
                rp_cum,norm_x,norm_y,fx2,fy2,popt,pcov = fitLogFunc(rp_clean)
                #print('4')
                xr = np.linspace(1e-3,rp_cum.r_arcsec.max(),100)
                #print(popt)
                ys =logfunc(xr,*popt)
                ind = np.where(ys>0.5)[0][0]
                mask = ys>0
                halflight = np.abs(norm_x-xr[mask][ind])
                halflights.append(halflight)
            except:
                halflights.append('fit fail')
                
    return(pgcs,halflights)
