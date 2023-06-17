import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.units as u
from antares_client.search import search
from elasticsearch_dsl import Search
from .bandpassdict import *
import sncosmo
from astropy.table import Table
from astro_ghost.ghostHelperFunctions import getTransientHosts
# from yse import yse_object
from . import alerce_api,yse
import tempfile
from astropy.coordinates import SkyCoord

def query_antares(query=None):
    '''
    Queries ANTARES for new objects. See https://nsf-noirlab.gitlab.io/csdc/antares/client/tutorial/searching.html#advanced-searches for information on how to construct a query.

    Parameters
    ----------
    query : elasticsearch_dsl.search.Search, optional
        The query to be sent to ANTARES. The default is None, which will return all objects with at least 4 detections and at least one detection in the last 7 days.

    Returns
    -------
    s : Iterator
        An iterator over the result loci of the query.
    '''
    today = Time.now()
    if query ==None:
        query = (
        Search()
        .filter("range", **{"properties.num_mag_values": {"gte": 4, "lte": 100}})
        .filter('range',**{'properties.oldest_alert_observation_time': {"gte":(today-7*u.day).mjd}})
        # .filter("term", tags="extragalactic")
        # .filter("term", tags="high_amplitude_transient_candidate")
        .to_dict()
        )
    s = search(query)
    return s

class antares_object:

    def __init__(self,locus):
        '''
        Instantiates an antares object from an antares Locus object.
        See https://nsf-noirlab.gitlab.io/csdc/antares/client/index.html for information on the antares client.

        Parameters
        ----------
        locus : antares_client.models.Locus
            An antares Locus object.
        
        Attributes
        ----------
        locus : antares_client.models.Locus
            An antares Locus object.
        lc : DataFrame
            The lightcurve of the object.
        name : str
            The ZTF Object ID of the object.
        ra : float
            The right ascension of the object.
        dec : float
            The declination of the object.
        url : str
            The URL of the object in the antares database.
        '''
        self.locus =locus
        self.lc = self.locus.lightcurve
        self.name =self.locus.properties['ztf_object_id']
        self.ra,self.dec = locus.ra, locus.dec
        self.url = f"https://antares.noirlab.edu/loci/{self.locus.locus_id}"
    
    def check_yse(self,verbose=True):
        '''
        Cross-checks to see if object exists in YSE. Checks if the object is in the TNS public catalog, and if so, checks if it is in the YSE database.

        Parameters
        ----------
        verbose : bool, optional
            Whether to print out whether the object was found in YSE. The default is True.
        
        Returns
        -------
        yse_object : prep.source.yse.yse_object
            The YSE object if it exists, or None if it does not.
        '''
        if 'tns_public_objects' in self.locus.catalogs:
            name = self.locus.catalog_objects['tns_public_objects'][0]['name']
            try:
                yo=yse.yse_object(name)
                if verbose:
                    print('Found in YSE')
                
                return yse.yse_object(name)
            except:
                if verbose:
                    print('Not found in YSE')
                return
        else:
            print('Could not find TNS name')
            return
        
    def check_alerce(self,verbose=False):
        '''
        Cross-checks to see if object exists in ALeRCE. Uses the ZTF Object ID to query the ALeRCE database.

        Parameters
        ----------
        verbose : bool, optional
            Whether to print out whether the object was found in ALeRCE. The default is False.

        Returns
        -------
        alerce_object : prep.source.alerce_api.alerce_object
            The ALeRCE object if it exists, or None if it does not.
        '''
        q1 = alerce_api.query_alerce({'oid':self.name,'format':'pandas'})
        if q1.empty:
            if verbose:
                print('Not found in ALeRCE')
            return
        else:
            if verbose:
                print('Found in ALeRCE')
            return alerce_api.alerce_object(q1.iloc[0])
    
    def get_lc(self):
        '''
        Gets the lightcurve of the object. Redundant with the lc attribute. Made for consistency with other object classes.

        Returns
        -------
        lc : DataFrame
            The lightcurve of the object.
        '''
        return self.lc
    
    def salt3(self,plot=False):
        '''
        Fits a SALT3 model to the lightcurve of the object. Uses the sncosmo package.

        Parameters
        ----------
        plot : bool, optional
            Whether to plot the lightcurve and the fit. The default is False.
        
        Attributes Generated
        --------------------
        salt_params : dict
            The parameters of the fit and errors.
        '''
        self.lcm = self.lc[~self.lc[['ant_mjd','ant_mag','ant_magerr']].isna().any(axis=1)]
        # Need to figure out ztf passbands
        # ztf  - g and r
        # doing the same thing as for yse
        salt2band = []
        zpsys =[]
        mask = []
        for i in range(len(self.lcm.ant_mjd.values)):
            if self.lcm.ant_passband.values[i] == 'g':
                salt2band.append(bandpassdict['Band: ZTF-Cam - g-ZTF'])
                mask.append(True)
            elif self.lcm.ant_passband.values[i] == 'R':
                salt2band.append(bandpassdict['Band: ZTF-Cam - r-ZTF'])
                mask.append(True)
            else:
                salt2band.append('None')
                mask.append(False)
            zpsys.append('AB')
        
        model = sncosmo.Model(source='salt3')
        fitparams = ['z', 't0', 'x0', 'x1', 'c']
        mask = mask #& (self.pdata['FLUXCAL']<1e10) #& (pdata['MJD'].values>60052)
        # model.set(z=sn['Host Redshift'])
        salt2mjd = self.lcm['ant_mjd'][mask].values
        mag = self.lcm['ant_mag'][mask].values
        flux = 10**(-0.4*(mag-27.5))
        # flux = pdata['FLUXCAL'][mask].values
        mag_err=self.lcm['ant_magerr'][mask].values
        # mag_err[(mag_err < 0.01)|(mag_err ==None)]=0.01
        fluxerr = flux*mag_err*0.4*np.log(10)
        # fluxerr = pdata['FLUXCALERR'][mask].values
        zp = np.array([27.5]*len(salt2band))[mask]
        zpsys = np.array(zpsys)
        salt2band = np.array(salt2band)
        zpsys = zpsys[mask]
        salt2band = salt2band[mask]

        data = Table([salt2mjd,salt2band,flux,fluxerr,zp,zpsys],names=['mjd','band','flux','fluxerr','zp','zpsys'],meta={'t0':salt2mjd[flux == np.max(flux)]})
        result, fitted_model = sncosmo.fit_lc(
                    data, model, fitparams,
                    bounds={'t0':(salt2mjd[flux == np.max(flux)]-10, salt2mjd[flux == np.max(flux)]+10),
                            'z':(0.0,0.7),'x1':(-3,3),'c':(-0.3,0.3)})  # bounds on parameters (if any)
        print('points fitted =',len(salt2mjd))
        print('chisq =',result['chisq'])
        today = Time.now().mjd
        lcphase = today-result['parameters'][1]
        lcp =lcphase
        if plot:
            if lcphase > 0: lcphase = '+%.1f'%(lcphase)
            else: lcphase = '%.1f'%(lcphase)
            print("phase = %s days"%(lcphase))
            print("z = %.3f"%(result['parameters'][0]))
            print('t0 = %i'%(result['parameters'][1]))
            print('ms = %.2f'%(10.635-2.5*np.log10(result['parameters'][2])))
            print('x1 = %.2f'%(result['parameters'][3]))
            print('c = %.2f'%(result['parameters'][4]))
            # for i in range(len(result['param_names'])):
            #     if i==2:
            #         print('ms =',10.635-2.5*np.log10(result.parameters[i]))
            #     else:
            #         rp=result['parameters'][i]
            #         print(result['param_names'][i],'=',f'%i'%(rp))
        
        params = np.array([result['chisq'],lcp, result['parameters'][0],
                                      result['parameters'][1],
                                      10.635-2.5*np.log10(result['parameters'][2]),
                                      result['parameters'][3],
                                      result['parameters'][4],
                                      len(salt2mjd)]+list(result['errors'].values()))
        pnames = ['chisq','phase','z','t0','ms','x1','c','npoints']+list(key+'_err' for key in result['errors'].keys())

        self.salt_params = dict(zip(pnames,params))
        if plot:
            plt.figure(figsize=(11,8))
            plt.title(f"{self.name} (ANTARES)")
            plt.xlabel('MJD')
            plt.ylabel('Mag')
            p1 = self.lcm.query('ant_passband == "g"')
            plt.scatter(p1.ant_mjd,p1.ant_mag,label='ZTF - g',color='g')
            p2 = self.lcm.query('ant_passband == "R"')
            plt.scatter(p2.ant_mjd,p2.ant_mag,label='ZTF - r',color='r')
            plotmjd = np.arange(result['parameters'][1]-20,result['parameters'][1]+50,0.5)
            for band in np.unique(salt2band):
                salt2flux = fitted_model.bandflux(band, plotmjd, zp=27.5,zpsys=zpsys[band == salt2band][0])
                plt.plot(plotmjd,-2.5*np.log10(salt2flux)+27.5,label=band)
            plt.gca().invert_yaxis()
            plt.grid()
            plt.legend(loc=0)
        
        # return self.salt_params
    


    def plot_lc(self):
        '''
        Plots the light curve data.
        '''
        plt.title(self.name)
        plt.xlabel('MJD')
        plt.ylabel('Magnitude')
        plt.scatter(self.lc.ant_mjd,self.lc.ant_mag,c='k')
        plt.gca().invert_yaxis()

    def run_ghost(self):
        '''
        Runs GHOST for the object and returns the host galaxy information.

        Returns
        -------
        ghost_info : DataFrame
            DataFrame containing the host galaxy information.
        '''
        snCoord = SkyCoord(self.ra,self.dec,unit='deg',frame='icrs')
        with tempfile.TemporaryDirectory() as tmp:
            ghost_info = getTransientHosts(snCoord=[snCoord], verbose=0,starcut='normal',savepath=tmp,ascentMatch=True)
        return ghost_info
    