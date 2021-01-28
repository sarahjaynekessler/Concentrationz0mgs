import pandas as pd
import numpy as np

def makeConcentrationVars(df,bands,suffix1,suffix2,prefix=None):
    lsun=3.846e+33
    if prefix!=None:

        varname1 = prefix+'_'+suffix1
        varname2 = prefix+'_'+suffix2

        totalConcentration = df[varname1]/df[varname2]

        sigma1 = df[varname1]/df['area_2kpc_pc']
        sigma2 = df[varname2]/df['area_R25_pc']

        sigmaConcentration = sigma1/sigma2

        df['c_'+prefix] = totalConcentration
        df['c_sigma_'+prefix]= sigmaConcentration

        df['sigma_'+prefix+'_'+suffix1]= sigma1
        df['sigma_'+prefix+'_'+suffix2]= sigma2


    else:    

        for b in bands:
            b=b.upper()

            varname1 = b.upper()+'_'+suffix1
            varname2 = b.upper()+'_'+suffix2

            totalConcentration = df[varname1]/df[varname2]

            sigma1 = df[varname1]/df['area_2kpc_pc']
            sigma2 = df[varname2]/df['area_R25_pc']

            sigmaConcentration = sigma1/sigma2
            
            Lsununits1 = df[varname1]/lsun
            Lsununits2 = df[varname2]/lsun

            LsunConcentration = Lsununits1/Lsununits2

            sigmaLsunUnits1 = Lsununits1/df['area_2kpc_pc']
            sigmaLsunUnits2 = Lsununits2/df['area_R25_pc']

            sigmaLsunConcentration = sigmaLsunUnits1/sigmaLsunUnits2

            df['c_'+b] = totalConcentration
            df['c_sigma_'+b]= sigmaConcentration
            df['c_sigma_'+b+'_Lsun']= sigmaLsunConcentration 
            df['c_'+b+'_Lsun']= LsunConcentration 

            df['sigma_'+b+'_'+suffix1]= sigma1
            df['sigma_'+b+'_'+suffix2]= sigma2
            df['sigma_'+b+'_Lsun_'+suffix1]= sigmaLsunUnits1 
            df['sigma_'+b+'_Lsun_'+suffix2]= sigmaLsunUnits2
            df[b+'_Lsun_'+suffix1]= Lsununits1
            df[b+'_Lsun_'+suffix2]= Lsununits2

    return(df)

        

