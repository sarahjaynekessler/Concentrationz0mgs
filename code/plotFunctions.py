import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import pandasFunctions
import warnings


sns.set_style('whitegrid')

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


warnings.filterwarnings('ignore')
bands = ['w1','w2','w3','w4','nuv','fuv']
cs = ['goldenrod','darkorange','red','darkred','dodgerblue','darkviolet']

def factorplot(df):

    g = sns.FacetGrid(df, col='GALTYPE_COURSE',row=2,height = 4,aspect = 0.8)

    for b in np.arange(len(bands)):
        g.map(sns.histplot,bands[b],stat='probability',binwidth=0.1,label = bands[b],color=cs[b],element='step',fill=True)

        g.map(sns.histplot,'c_'+bands[b],stat='probability',binwidth=0.1,label = bands[b],color=cs[b],element='step',fill=True)

    axes = g.axes.flatten()

    axes[0].set_title("Late Type")
    axes[1].set_title("Early Type")
    axes[1].set_xlabel('')
    axes[0].set_xlabel('')

    g.fig.text(x=0.5,y=0.05,horizontalalignment='center',s = r"$\Sigma$(Band)$_{2kpc}$ / $\Sigma$(Band)$_{R25}$")

    plt.xlim(-1.2,10)
    plt.legend()
    plt.show()

def fourPanelHist(df):

    fig, axes = plt.subplots(2,2, sharex=False, sharey=False,figsize = (8,8))

    for b in np.arange(len(bands)):
        sns.histplot(df[df['GALTYPE_COURSE']=='LT'][bands[b]],stat='density',binwidth=0.1,label = bands[b],color=cs[b],element='step',fill=True,ax=axes[0,0])
        axes[0,0].set_xlim(-1.2,2)
        axes[0,0].set_ylim(0,3)
        axes[0,0].set_xlabel('')
        axes[0,0].set_ylabel('')


        sns.histplot(df[df['GALTYPE_COURSE']=='ET'][bands[b]],stat='density',binwidth=0.1,color=cs[b],element='step',fill=True,ax=axes[0,1])
        axes[0,1].set_xlim(-1.2,2)
        axes[0,1].set_ylim(0,3)
        axes[0,1].set_xlabel('')
        axes[0,1].set_ylabel('')
        axes[0,1].tick_params(labelleft=False) 

        sns.histplot(df[df['GALTYPE_COURSE']=='LT']['c_'+bands[b]],stat='density',binwidth=0.1,color=cs[b],element='step',fill=True,ax=axes[1,0])
        axes[1,0].set_xlim(-1.2,10)
        axes[1,0].set_ylim(0,.5)
        axes[1,0].set_xlabel('')
        axes[1,0].set_ylabel('')

        sns.histplot(df[df['GALTYPE_COURSE']=='ET']['c_'+bands[b]],stat='density',binwidth=0.1,color=cs[b],element='step',fill=True,ax=axes[1,1])
        axes[1,1].set_xlim(-1.2,10)
        axes[1,1].set_ylim(0,.5)
        axes[1,1].set_xlabel('')
        axes[1,1].set_ylabel('')
        axes[1,1].tick_params(labelleft=False) 


    axes[0,0].set_title('Late Type')

    axes[0,1].set_title('Early Type')


    plt.subplots_adjust(wspace = .05)
    fig.legend(bbox_to_anchor=(1., .4))
    fig.subplots_adjust(right=.8)
    plt.show()


def galTypeLumSigma(df):

    fig, axes = plt.subplots(2,2, sharex=False, sharey=False,figsize = (8,8))

    for b in np.arange(len(bands)):
        b=bands[b].upper()
        ystr = 'sigma_'+str(bands[b])+'_1p0R25'
        ystr_L = bands[b].upper()+'_1p0R25'

        mask = df.GALTYPE_COURSE == 'LT'

        med = pandasFunctions.rollingmedainXY(df[mask],'sigma_Mstar_R25',np.arange(-1.5,2.5,.2))
        axes[1,0].plot(med.sigma_Mstar_R25,np.log10(med[ystr]),color=cs[b],lw=3)

        med = pandasFunctions.rollingmedainXY(df[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[0,0].plot(med.LOGMASS,np.log10(med[ystr_L]),color=cs[b],lw=3)


        mask = df.GALTYPE_COURSE == 'ET'
        med = pandasFunctions.rollingmedainXY(df[mask],'sigma_Mstar_R25',np.arange(-1.5,2.5,.2))
        axes[1,1].plot(med.sigma_Mstar_R25,np.log10(med[ystr]),color=cs[b],lw=3)

        med = pandasFunctions.rollingmedainXY(df[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[0,1].plot(med.LOGMASS,np.log10(med[ystr_L]),label=bands[b].upper(),color=cs[b],lw=3)



    axes[0,0].set_xlabel(r'log(M$_*$) [M$_{\odot}$]')
    axes[0,0].set_ylabel(r'log($\nu$L$_{\nu}$(Band))$_{R25}$ [erg/s]')


    axes[0,1].set_xlabel(r'log(M$_*$) [M$_{\odot}$]')
    axes[0,1].set_ylabel('')

    axes[1,0].set_xlabel(r'log($\Sigma$(M$_*$)) [M$_{\odot}$/pc$^{2}$]')
    axes[1,0].set_ylabel(r'log($\Sigma$(Band)$_{R25}$) [erg/s/kpc$^{2}$]')

    axes[1,1].set_xlabel(r'log($\Sigma$(M$_*$)) [M$_{\odot}$/pc$^{2}$]')
    axes[1,1].set_ylabel('')

    axes[0,0].set_title('Late Type',fontsize=16)
    axes[0,1].set_title('Early Type',fontsize=16)

    fig.legend(loc="center right",borderaxespad=0.1)
    plt.subplots_adjust(right=0.85)
    plt.savefig('/Users/kessler.363/Thesis/Concentrationz0mgs/plots/R25_L_sigma_galtype.png',bbox_inches='tight')



def galtypeConcentrationBySigma(df):

    fig, axes = plt.subplots(2,2, sharex=False, sharey=False,figsize = (8,8))

    for b in np.arange(len(bands)):
        b=bands[b].upper()
        ystr = 'c_sigma_'+str(bands[b])
        ystr_L ='c_'+str(bands[b])

        mask = df.GALTYPE_COURSE == 'LT'

        med = pandasFunctions.rollingmedainXY(df[mask],'sigma_Mstar_R25',np.arange(-1.5,2.5,.2))
        axes[1,0].plot(med.sigma_Mstar_R25,med[ystr],color=cs[b],lw=3)

        med = pandasFunctions.rollingmedainXY(df[mask],'LOGMASS',np.arange(8.5,12,.5))
        axes[0,0].plot(med.LOGMASS,med[ystr_L],color=cs[b],lw=3)


        mask = df.GALTYPE_COURSE == 'ET'
        med = pandasFunctions.rollingmedainXY(df[mask],'sigma_Mstar_R25',np.arange(-1.5,2.5,.2))
        axes[1,1].plot(med.sigma_Mstar_R25,med[ystr],color=cs[b],lw=3)

        med = pandasFunctions.rollingmedainXY(df[mask],'LOGMASS',np.arange(8.5,12,.5))
        axes[0,1].plot(med.LOGMASS,med[ystr_L],label=bands[b].upper(),color=cs[b],lw=3)


    axes[0,0].set_xlabel(r'log(M$_*$) [M$_{\odot}$]')
    axes[0,0].set_ylabel(r'$\nu$L$_{\nu}$(Band)$_{2kpc}$/$\nu$L$_{\nu}$(Band)$_{R25}$')


    axes[0,1].set_xlabel(r'log(M$_*$) [M$_{\odot}$]')
    axes[0,1].set_ylabel('')


    axes[1,0].set_xlabel(r'log($\Sigma$(M$_*$)) [M$_{\odot}$/pc$^{2}$]')
    axes[1,0].set_ylabel(r'$\Sigma$(Band)$_{2kpc}$ /$\Sigma$(Band)$_{R25}$')

    axes[1,1].set_xlabel(r'log($\Sigma$(M$_*$)) [M$_{\odot}$/pc$^{2}$]')
    axes[1,1].set_ylabel('')

    axes[1,0].set_xlim(xmin=-0.5)
    axes[1,1].set_xlim(xmin=0)

    axes[1,0].set_ylim(0,5)
    axes[1,1].set_ylim(0,11)


    axes[0,0].set_title('Late Type',fontsize=16)
    axes[0,1].set_title('Early Type',fontsize=16)

    fig.legend(loc="center right",borderaxespad=0.1)
    plt.subplots_adjust(right=0.85)
    plt.savefig('/Users/kessler.363/Thesis/Concentrationz0mgs/plots/concentration_sigma_galtype.png',bbox_inches='tight')


def galtypeConcentrationByMass(df):

    fig, axes = plt.subplots(2,2, sharex=False, sharey=False,figsize = (8,8))

    for b in np.arange(len(bands)):
        b=bands[b].upper()
        ystr = 'c_sigma_'+str(bands[b])
        ystr_L ='c_'+str(bands[b])

        mask = df.GALTYPE_COURSE == 'LT'


        med = pandasFunctions.rollingmedainXY(df[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[1,0].plot(med.LOGMASS,np.log10(med[ystr]),color=cs[b],lw=3)
        axes[0,0].plot(med.LOGMASS,np.log10(med[ystr_L]),color=cs[b],lw=3)


        mask = df.GALTYPE_COURSE == 'ET'

        med = pandasFunctions.rollingmedainXY(df[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[1,1].plot(med.LOGMASS,np.log10(med[ystr]),color=cs[b],lw=3)
        axes[0,1].plot(med.LOGMASS,np.log10(med[ystr_L]),label=bands[b].upper(),color=cs[b],lw=3)


    axes[0,0].set_ylabel(r'log($\nu$L$_{\nu}$(Band)$_{2kpc}$/$\nu$L$_{\nu}$(Band)$_{R25}$)')
    axes[0,1].set_ylabel('')
    axes[1,0].set_ylabel(r'log($\Sigma$(Band)$_{2kpc}$ /$\Sigma$(Band)$_{R25}$)')
    axes[1,1].set_ylabel('')

    fig.text(0.5, 0.04, r'log(M$_*$) [M$_{\odot}$]', ha='center',fontsize=14)


    axes[0,0].set_title('Late Type',fontsize=16)
    axes[0,1].set_title('Early Type',fontsize=16)

    fig.legend(loc="center right",borderaxespad=0.1)
    plt.subplots_adjust(right=0.85)
    plt.savefig('/Users/kessler.363/Thesis/Concentrationz0mgs/plots/concentration_sigma_byLMASS_galtype.png',bbox_inches='tight')


#now in solar lum

def galTypeLumSigmaSolLum(df):
    fig, axes = plt.subplots(2,2, sharex=False, sharey=False,figsize = (8,8))

    for b in np.arange(len(bands)):
        ystr = 'sigma_'+str(bands[b].upper())+'_Lsun_1p0R25'
        ystr_L = bands[b].upper()+'_Lsun_1p0R25'
        
        dfn = df.copy()
        dfn[ystr] = np.log10(dfn[ystr])
        dfn[ystr_L] = np.log10(dfn[ystr_L])
        
        mask = df.GALTYPE_COURSE == 'LT'

        med = pandasFunctions.rollingmedainXY(dfn.copy()[mask],'sigma_Mstar_R25',np.arange(-1.5,2.5,.2))
        axes[1,0].plot(med.sigma_Mstar_R25,med[ystr],color=cs[b],lw=3)

        med = pandasFunctions.rollingmedainXY(dfn.copy()[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[0,0].plot(med.LOGMASS,med[ystr_L],color=cs[b],lw=3)


        mask = df.GALTYPE_COURSE == 'ET'
        med = pandasFunctions.rollingmedainXY(dfn.copy()[mask],'sigma_Mstar_R25',np.arange(-1.5,2.5,.2))
        axes[1,1].plot(med.sigma_Mstar_R25,med[ystr],color=cs[b],lw=3)

        med = pandasFunctions.rollingmedainXY(dfn.copy()[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[0,1].plot(med.LOGMASS,med[ystr_L],label=bands[b].upper(),color=cs[b],lw=3)



    axes[0,0].set_xlabel(r'log(M$_*$) [M$_{\odot}$]')
    axes[0,0].set_ylabel(r'log($\nu$L$_{\nu}$(Band))$_{R25}$ [L$_{\odot}$]')


    axes[0,1].set_xlabel(r'log(M$_*$) [M$_{\odot}$]')
    axes[0,1].set_ylabel('')

    axes[1,0].set_xlabel(r'log($\Sigma$(M$_*$)) [M$_{\odot}$/pc$^{2}$]')
    axes[1,0].set_ylabel(r'log($\Sigma$(Band)$_{R25}$) [L$_{\odot}$ kpc$^{-2}$]')

    axes[1,1].set_xlabel(r'log($\Sigma$(M$_*$)) [M$_{\odot}$/pc$^{2}$]')
    axes[1,1].set_ylabel('')

    axes[0,0].set_title('Late Type',fontsize=16)
    axes[0,1].set_title('Early Type',fontsize=16)

    fig.legend(loc="center right",borderaxespad=0.1)
    plt.subplots_adjust(right=0.85)
    plt.savefig('/Users/kessler.363/Thesis/Concentrationz0mgs/plots/R25_Lsun_sigma_galtype.png',bbox_inches='tight')



def galtypeConcentrationByMassSolLum(df):

    fig, axes = plt.subplots(2,2, sharex=False, sharey=False,figsize = (8,8))

    for b in np.arange(len(bands)):
        ystr = 'c_sigma_'+str(bands[b].upper())+'_Lsun'
        ystr_L = 'c_'+bands[b].upper()+'_Lsun'
        
        dfn = df.copy()
        dfn[ystr] = np.log10(dfn[ystr])
        dfn[ystr_L] = np.log10(dfn[ystr_L])

        mask = df.GALTYPE_COURSE == 'LT'


        med = pandasFunctions.rollingmedainXY(dfn.copy()[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[1,0].plot(med.LOGMASS,med[ystr],color=cs[b],lw=3)
        axes[0,0].plot(med.LOGMASS,med[ystr_L],color=cs[b],lw=3)


        mask = df.GALTYPE_COURSE == 'ET'

        med = pandasFunctions.rollingmedainXY(dfn.copy()[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[1,1].plot(med.LOGMASS,med[ystr],color=cs[b],lw=3)
        axes[0,1].plot(med.LOGMASS,med[ystr_L],label=bands[b].upper(),color=cs[b],lw=3)


    axes[0,0].set_ylabel(r'log($\nu$L$_{\nu}$(Band)$_{2kpc}$/$\nu$L$_{\nu}$(Band)$_{R25}$)')
    axes[0,1].set_ylabel('')
    axes[1,0].set_ylabel(r'log($\Sigma$(Band)$_{2kpc}$ /$\Sigma$(Band)$_{R25}$)')
    axes[1,1].set_ylabel('')

    fig.text(0.5, 0.04, r'log(M$_*$) [M$_{\odot}$]', ha='center',fontsize=14)


    axes[0,0].set_title('Late Type',fontsize=16)
    axes[0,1].set_title('Early Type',fontsize=16)

    fig.legend(loc="center right",borderaxespad=0.1)
    plt.subplots_adjust(right=0.85)
    plt.savefig('/Users/kessler.363/Thesis/Concentrationz0mgs/plots/concentration_sigma_byLMASS_Lsun_galtype.png',bbox_inches='tight')


#Residual plots!

def ResidsGalTypeLumSigmaSolLum(df):
    fig, axes = plt.subplots(2,2, sharex='row', sharey='row',figsize = (8,8))

    sigmabins = np.arange(-1.5,2.5,.2)[:-1]
    massbins = np.arange(8.5,12,.2)[:-1]

    for b in np.arange(len(bands)):
        ystr = 'sigma_'+str(bands[b].upper())+'_Lsun_1p0R25'
        ystr_L = bands[b].upper()+'_Lsun_1p0R25'

        dfn = df.copy()
        dfn[ystr] = np.log10(dfn[ystr])
        dfn[ystr_L] = np.log10(dfn[ystr_L])
        mask = df.GALTYPE_COURSE == 'LT'

        med = pandasFunctions.rollingstdXY(dfn.copy()[mask],'sigma_Mstar_R25',np.arange(-1.5,2.5,.2))
        axes[1,0].bar(sigmabins+.1,med[ystr],color=cs[b],alpha=0.6,width=0.2)

        med = pandasFunctions.rollingstdXY(dfn.copy()[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[0,0].bar(massbins+.1,med[ystr_L],color=cs[b],alpha=0.6,width=0.2)

        mask = df.GALTYPE_COURSE == 'ET'
        med = pandasFunctions.rollingstdXY(dfn.copy()[mask],'sigma_Mstar_R25',np.arange(-1.5,2.5,.2))
        axes[1,1].bar(sigmabins+.1,med[ystr],color=cs[b],alpha=0.6,width=0.2)

        med = pandasFunctions.rollingstdXY(dfn.copy()[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[0,1].bar(massbins+.1,med[ystr_L],color=cs[b],alpha=0.6,width=0.2,label = bands[b].upper())



    axes[0,0].set_xlabel(r'log(M$_*$) [M$_{\odot}$]')
    axes[0,1].set_xlabel(r'log(M$_*$) [M$_{\odot}$]')
    axes[0,0].set_ylabel(r'$\sigma$ ($\nu$L$_{\nu}$(Band)$_{R25}$) [dex]')

    axes[1,0].set_xlabel(r'log($\Sigma$(M$_*$)) [M$_{\odot}$/pc$^{2}$]')
    axes[1,1].set_xlabel(r'log($\Sigma$(M$_*$)) [M$_{\odot}$/pc$^{2}$]')
    axes[1,0].set_ylabel(r'$\sigma$($\Sigma$(Band)$_{R25}$) [dex]')


    axes[0,0].set_title('Late Type',fontsize=16)
    axes[0,1].set_title('Early Type',fontsize=16)

    fig.legend(loc="center right",borderaxespad=0.1)
    plt.subplots_adjust(right=0.85)

    plt.savefig('/Users/kessler.363/Thesis/Concentrationz0mgs/plots/Resids_R25_Lsun_sigma_galtype.png',bbox_inches='tight')


def ResidsGaltypeConcentrationByMassSolLum(df):

    fig, axes = plt.subplots(2,2, sharex='row', sharey='row',figsize = (8,8))

    for b in np.arange(len(bands)):

        massbins = np.arange(8.5,12,.2)[:-1]
        
        ystr = 'c_sigma_'+str(bands[b].upper())+'_Lsun'
        ystr_L = 'c_'+bands[b].upper()+'_Lsun'
        
        dfn = df.copy()
        dfn[ystr] = np.log10(dfn[ystr])
        dfn[ystr_L] = np.log10(dfn[ystr_L])

        mask = df.GALTYPE_COURSE == 'LT'


        med = pandasFunctions.rollingstdXY(dfn.copy()[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[0,0].bar(massbins+.1,med[ystr_L],color=cs[b],alpha=0.6,width=0.2)
        axes[1,0].bar(massbins+.1,med[ystr],color=cs[b],alpha=0.6,width=0.2)


        mask = df.GALTYPE_COURSE == 'ET'

        med = pandasFunctions.rollingstdXY(dfn.copy()[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[1,1].bar(massbins+.1,med[ystr],color=cs[b],alpha=0.6,width=0.2)
        axes[0,1].bar(massbins+.1,med[ystr_L],color=cs[b],alpha=0.6,width=0.2,label = bands[b].upper())


    axes[0,0].set_ylabel(r'$\sigma$ ($\nu$L$_{\nu}$(Band)$_{2kpc}$/$\nu$L$_{\nu}$(Band)$_{R25}$) [dex]')
    axes[0,1].set_ylabel('')
    axes[1,0].set_ylabel(r'$\sigma$ ($\Sigma$(Band)$_{2kpc}$ /$\Sigma$(Band)$_{R25}$) [dex]')
    axes[1,1].set_ylabel('')

    fig.text(0.5, 0.04, r'log(M$_*$) [M$_{\odot}$]', ha='center',fontsize=14)


    axes[0,0].set_title('Late Type',fontsize=16)
    axes[0,1].set_title('Early Type',fontsize=16)

    fig.legend(loc="center right",borderaxespad=0.1)
    plt.subplots_adjust(right=0.85)
    plt.savefig('/Users/kessler.363/Thesis/Concentrationz0mgs/plots/Resids_concentration_sigma_byLMASS_Lsun_galtype.png',bbox_inches='tight')


def ResidsGalTypeLumSigma(df):
    fig, axes = plt.subplots(2,2, sharex='row', sharey='row',figsize = (8,8))

    sigmabins = np.arange(-1.5,2.5,.2)[:-1]
    massbins = np.arange(8.5,12,.2)[:-1]

    for b in np.arange(len(bands)):
        ystr = 'c_sigma_'+str(bands[b].upper())+'_Lsun'
        ystr_L = 'c_'+bands[b].upper()+'_Lsun'

        dfn = df.copy()
        dfn[ystr] = np.log10(dfn[ystr])
        dfn[ystr_L] = np.log10(dfn[ystr_L])
        mask = df.GALTYPE_COURSE == 'LT'

        med = pandasFunctions.rollingstdXY(dfn.copy()[mask],'sigma_Mstar_R25',np.arange(-1.5,2.5,.2))
        axes[1,0].bar(sigmabins+.1,med[ystr],color=cs[b],alpha=0.6,width=0.2)

        med = pandasFunctions.rollingstdXY(dfn.copy()[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[0,0].bar(massbins+.1,med[ystr_L],color=cs[b],alpha=0.6,width=0.2)

        mask = df.GALTYPE_COURSE == 'ET'
        med = pandasFunctions.rollingstdXY(dfn.copy()[mask],'sigma_Mstar_R25',np.arange(-1.5,2.5,.2))
        axes[1,1].bar(sigmabins+.1,med[ystr],color=cs[b],alpha=0.6,width=0.2)

        med = pandasFunctions.rollingstdXY(dfn.copy()[mask],'LOGMASS',np.arange(8.5,12,.2))
        axes[0,1].bar(massbins+.1,med[ystr_L],color=cs[b],alpha=0.6,width=0.2,label = bands[b].upper())



    axes[0,0].set_xlabel(r'log(M$_*$) [M$_{\odot}$]')
    axes[0,1].set_xlabel(r'log(M$_*$) [M$_{\odot}$]')
    axes[0,0].set_ylabel(r'$\sigma$ ($\nu$L$_{\nu}$(Band)$_{2kpc}$/$\nu$L$_{\nu}$(Band)$_{R25}$) [dex]')

    axes[1,0].set_xlabel(r'log($\Sigma$(M$_*$)) [M$_{\odot}$/pc$^{2}$]')
    axes[1,1].set_xlabel(r'log($\Sigma$(M$_*$)) [M$_{\odot}$/pc$^{2}$]')
    axes[1,0].set_ylabel(r'$\sigma$ ($\Sigma$(Band)$_{2kpc}$ /$\Sigma$(Band)$_{R25}$) [dex]')


    axes[0,0].set_title('Late Type',fontsize=16)
    axes[0,1].set_title('Early Type',fontsize=16)

    fig.legend(loc="center right",borderaxespad=0.1)
    plt.subplots_adjust(right=0.85)

    plt.savefig('/Users/kessler.363/Thesis/Concentrationz0mgs/plots/Resids_R25_sigma_galtype.png',bbox_inches='tight')
    
    
def sixPanelSFMS(df,galtype,cols):

    bandsu = [b.upper() for b in bands]
    
    test = df[cols].copy()
    c = dict(zip(test.columns[:6],bandsu))
    test = test.rename(columns=c)
    
    lt = test[test.GALTYPE_COURSE == galtype]
    lt = lt.drop(columns = 'GALTYPE_COURSE')
    
    for b in bands:
        lt[b.upper()] = np.log10(lt[b.upper()])
        
    ltmelt = lt.melt(id_vars='LOGMASS')
    med = pandasFunctions.rollingmedainXY(lt,'LOGMASS',np.arange(8.5,12,.2))
    
    g = sns.FacetGrid(ltmelt,col='variable',col_wrap=2,height=2.75,aspect=1,sharex =True,sharey=True,hue = 'variable',palette=cs,despine=False)
    g.map_dataframe(sns.scatterplot,x='LOGMASS',y='value',alpha=0.6,marker='.')
    axes = g.fig.axes
    for b in np.arange(len(bandsu)):
        axes[b].plot(med.LOGMASS,med[bandsu[b]],color=cs[b],lw=3)
        axes[b].annotate('n='+str(len(lt[bandsu[b]][~np.isnan(lt[bandsu[b]])])),(10,38),fontsize=14)
    g.fig.subplots_adjust(wspace=0.05, hspace=0.05)
    g.add_legend()

    g._legend.set_title('')
    g.set_titles('')
    for lh in g._legend.legendHandles:
        lh.set_alpha(1)
        lh._sizes = [50]
    g.fig.text(x=0.5,y=0.05,horizontalalignment='center',s = r'log(M$_*$) [M$_{\odot}$]',fontsize=14)
    g.fig.text(x=0.05,y=0.5,verticalalignment='center',s = r'log($\nu$L$_{\nu}$(Band)$_{R25}$) [erg/s]',rotation=90,fontsize=14)
    g.fig.subplots_adjust(bottom=.1)
    g.fig.subplots_adjust(left=.15)
    plt.subplots_adjust(top=0.95)
    if galtype == 'LT':
        g.fig.suptitle('Late Type')
    else:
        g.fig.suptitle('Early Type')
        
        
    plt.savefig('/Users/kessler.363/Thesis/Concentrationz0mgs/plots/SFMS_scatter_'+galtype+'.png',bbox_inches='tight')
    plt.close('all')

