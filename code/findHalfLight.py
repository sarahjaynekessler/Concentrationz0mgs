from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from gal_data import gal_data
from astropy import units as u


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
        
def fitLogFunc(df,ind,b):
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
