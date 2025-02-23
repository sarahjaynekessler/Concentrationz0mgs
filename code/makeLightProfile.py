from glob import glob
from astropy.io import fits
from UsefulAstroData import getPixelSize
from astropy import wcs
from UsefulAstroData import radec_maps,deproject
from reproject import reproject_interp
import numpy as np
import astropy.units as u
import astropy.constants as const
import pandas as pd

def runApertureLoop(path,df,bands,inputres,ind=None):
    pgcnames = ['PGC'+ str(i) for i in df.PGC.astype('int')]
    if ind==None:
        ind = len(pgcnames)
    else:
        ind=ind
    bands = [i.lower() for i in bands]
    wavesum = {'fuv':1540*1e-4,'nuv':2310*1e-4,'w1':3.4,'w2':4.6,'w3':12,'w4':22}

    for i in np.arange(len(pgcnames[:ind])):
        print(len(pgcnames[:ind]) - i)
        for band in bands:
            file = glob(path+pgcnames[i]+'_*'+band+'*'+inputres+'.fits')
            stars = glob(path+pgcnames[i]+'_*'+band+'*'+inputres+'_stars.fits')
            if len(file) == 0:
                pass
            else:
                if len(stars) == 0:
                    lightProfile(file[0],0,df,i,band,wavesum[band],pgcnames[i])
                else:
                    lightProfile(file[0],stars[0],df,i,band,wavesum[band],pgcnames[i])


def lightProfile(f,stars,df,i,band,wavlength,pgc):
    hdulist = fits.open(f)[0]
    data = hdulist.data
    w = wcs.WCS(hdulist.header)
    ster = getPixelSize(w)
    data[np.isnan(data)] = 0
    data=data*ster*1e6*1e-23*(const.c.to(u.micron/u.s)/(wavlength*u.micron)).value #convert from sr to pixel and then to erg/s/cm^2/Hz
    if stars != 0:
        starmask = fits.open(stars)[0].data
        data[stars==1] = np.nan
    else:
        pass
        
    racen,deccen = df.iloc[i].RA_DEG,df.iloc[i].DEC_DEG
    pa,incl = df.iloc[i].POSANG_DEG,df.iloc[i].INCL_DEG
    
    ra,dec = radec_maps(hdulist.header)
    rgrid,tgrid = deproject(ra,dec,racen,deccen,pa,incl)
    
    medback = np.nanmedian(data[np.logical_and(rgrid>1.5*df.iloc[i].R25_DEG,
                                rgrid<2*df.iloc[i].R25_DEG)])
    data-=medback
    
    ind = np.unravel_index(np.argsort(rgrid, axis=None), rgrid.shape)
    
    xs = rgrid[ind]*3600 #in arcsec
    sorteddata = hdulist.data[ind]
    ys = np.cumsum(sorteddata)/np.cumsum(sorteddata).max()
    mask = xs<400.0
    xs = xs[mask]
    sorteddata = sorteddata[mask]
    ys = ys[mask]
    name = [str(pgc) for i in len(xs)]
    galdf = pd.DataFrame({'r_arcsec':xs,'I':sorteddata,'Norm_I':ys,'PGC':name})
    galdf.to_csv('galaxycsvs/'+pgc+'_'+band+'_df.csv',index=None)
    

def makeLargePickle(path,df,bands,inputres,ind=None):
    pgcnames = ['PGC'+ str(i) for i in df.PGC.astype('int')]
    if ind==None:
        ind = len(pgcnames)
    else:
        ind=ind
    bands = [i.lower() for i in bands]

    wavesum = {'fuv':1540*1e-4,'nuv':2310*1e-4,'w1':3.4,'w2':4.6,'w3':12,'w4':22}

    for i in np.arange(len(pgcnames[:ind])):
        print(len(pgcnames[:ind]) - i)
        for band in bands:
            file = glob(path+pgcnames[i]+'_*'+band+'*'+inputres+'.fits')
            #print(path+pgcnames[i]+'_*'+band+'*'+inputres+'.fits')
            stars = glob(path+pgcnames[i]+'_*'+band+'*'+inputres+'_stars.fits')
            wavlength = wavesum[band]


                
            if len(stars) == 0:
                stars = 0
            else:
                stars=stars[0]
                

            if len(file) == 0:
                print('no file')
                pass
                
            else:
                f = file[0]
                pgc = pgcnames[i]
                hdulist = fits.open(f)[0]
                data = hdulist.data
                w = wcs.WCS(hdulist.header)
                ster = getPixelSize(w)
                data[np.isnan(data)] = 0
                data=data*ster*1e6*1e-23*(const.c.to(u.micron/u.s)/(wavlength*u.micron)).value #convert from sr to pixel and then to erg/s/cm^2/Hz
                if stars != 0:
                    starmask = fits.open(stars)[0].data
                    data[stars==1] = np.nan
                else:
                    pass
                    
                racen,deccen = df.iloc[i].RA_DEG,df.iloc[i].DEC_DEG
                pa,incl = df.iloc[i].POSANG_DEG,df.iloc[i].INCL_DEG
                #print('made it past ra stuff')
                ra,dec = radec_maps(hdulist.header)
                rgrid,tgrid = deproject(ra,dec,racen,deccen,pa,incl)
                
                medback = np.nanmedian(data[np.logical_and(rgrid>1.5*df.iloc[i].R25_DEG,
                                            rgrid<2*df.iloc[i].R25_DEG)])
                data-=medback
                
                inds = np.unravel_index(np.argsort(rgrid, axis=None), rgrid.shape)
                
                xs = rgrid[inds]*3600 #in arcsec
                sorteddata = hdulist.data[inds]
                ys = np.cumsum(sorteddata)/np.cumsum(sorteddata).max()
                mask = xs<400.0
                xs = xs[mask]
                if len(xs) == 0:
                    galdf = pd.DataFrame({'r_arcsec':[np.nan],'I':[np.nan],'Norm_I':[np.nan],'PGC':[pgc]})
                else:
                    sorteddata = sorteddata[mask]
                    ys = ys[mask]
                    name = [pgc for j in np.arange(len(xs))]
                    galdf = pd.DataFrame({'r_arcsec':xs,'I':sorteddata,'Norm_I':ys,'PGC':name})
                #print(i)
                if i==0:
                    fulldf = galdf.copy()
                    #print('galdf copy')
                    #print(fulldf)
                else:
                    fulldf = fulldf.append(galdf,ignore_index=True)
                    #print('other gals')
                    #print(fulldf)

    return(fulldf)
