from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
from gal_data import gal_data
from astropy import units as u
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
from astropy.modeling.models import Sersic1D
import warnings
from scipy.stats import chisquare
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
    #df1 is galbasedf with including r25- get pgc names from this file
    #df2 is pickle file
    
    halflights = []
    amps = []
    ns = []
    mses = []

    #print('0')
    for i in np.arange(len(pgcs)):
        print(len(pgcs)-i)
        galmask = df2.PGC.isin([pgcs[i]])
        rp = df2.loc[galmask]
        #print(rp)
        pgc=pgcs[i]
        
        if np.isnan(rp.r_arcsec).all()==True:
            halflights.append(np.nan)
            amps.append(np.nan)
            ns.append(np.nan)
            mses.append(np.nan)
            #print('1')
        elif rp.r_arcsec.min()>200:
            halflights.append(np.nan)
            amps.append(np.nan)
            ns.append(np.nan)
            mses.append(np.nan)

            #print('2')
        else:
            try:
                #print('3')
                rp.r_arcsec/=3600.
                r25 = df1.loc[i].R25_DEG
                #print(r25)`=-0
                #r25= r25.tolist()[0]
                #print(r25)
                mask = rp.r_arcsec<2*r25
                #print(mask)
                rp = rp[mask]
                
                #mask = rp.I>0
                #rp = rp[mask]
                ind = np.where(rp.r_arcsec<.5*r25)[0][-1]

                sersic = models.Sersic1D(bounds = {'n':(0,14)})
                outlier_fit = fitting.FittingWithOutlierRemoval(fitting.LevMarLSQFitter(),sigma_clip, niter=3, sigma=2.5)
                fitted_model,filtered_data = outlier_fit(sersic,rp.r_arcsec,rp.I)#,weights=0.1*rp.I)
                filtered_data[:ind]=False

                fit = fitting.LevMarLSQFitter()
                fitted_model = fit(sersic,rp.r_arcsec[~filtered_data],rp.I[~filtered_data])#,weights=(0.1*rp.I[~filtered_data]))
                mse = np.nanmean((rp.I[~filtered_data] - fitted_model(rp.I[~filtered_data]))**2)
                print('')
                print(pgc)
                #print(np.round(mse,decimals=2))
                
                """
                if mse>3.:
                    ns_mse = []
                    for n in nrange:
                        sersic = models.Sersic1D(bounds = {'n':(n,n+1)})
                        outlier_fit = fitting.FittingWithOutlierRemoval(fitting.LevMarLSQFitter(),sigma_clip, niter=3, sigma=3)
                        fitted_model,filtered_data_chisq = outlier_fit(sersic,rp.r_arcsec[~filtered_data],rp.I[~filtered_data])#,weights=(0.1*rp.I[~filtered_data]))
                        m = np.nanmean((rp.I[~filtered_data] - fitted_model(rp.I[~filtered_data]))**2)
                        if m<mse:
                            mse=m
                            ns_mse.append(n)
                            print(mse,n)
                        else:
                            pass
                    if len(ns_mse)==0:
                        print(ns_mse)
                        print('none better')
                        pass
                    else:
                        sersic = models.Sersic1D(bounds = {'n':(ns_mse[-1]-1,ns_mse[-1]+1)})
                        fitted_model = fit(sersic,rp.r_arcsec[~filtered_data],rp.I[~filtered_data])
                """
                re = np.round(fitted_model.r_eff.value*3600*0.9,decimals=3)
                #re*=0.9
                #re = np.round(re,decimals=3)
                n = np.round(fitted_model.n.value,decimals=3)
                #print(re,saloratio,mm15ratio)
                #print(n,mmres.n[i],np.round(mmres['T'][i],decimals=2))

                
                print(re)
                
                if re>250.:
                    halflights.append(np.nan)
                    amps.append('fit fail')
                    ns.append('fit fail')
                    mses.append('fit fail')
                    
                elif re<5.:
                    halflights.append(np.nan)
                    amps.append('fit fail')
                    ns.append('fit fail')
                    mses.append('fit fail')
                    
                else:
                    halflights.append(re)
                    amps.append(fitted_model.amplitude.value)
                    ns.append(n)
                    mses.append(mse)
            except:
                halflights.append(np.nan)
                amps.append('fit fail')
                ns.append('fit fail')
                mses.append('fit fail')


    redf = pd.DataFrame({'PGC':pgcs,'re':halflights,'amp': amps,'n':ns})
    #pgcs = np.asarray(pgcs)
    #halflights = np.asarray(halflights,dtype='float')
    #amps = np.asarray(amps)
    #ns = np.asarray(ns)
    #return(pgcs,halflights,amps,ns,datamasks,nummasked)
    return(redf)
    
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
