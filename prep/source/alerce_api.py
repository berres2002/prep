import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.units as u
from alerce.core import Alerce
from .bandpassdict import *
import sncosmo
from astropy.table import Table
from astro_ghost.ghostHelperFunctions import getTransientHosts
from antares_client.search import get_by_ztf_object_id
from . import antares
import tempfile
from astropy.coordinates import SkyCoord
alerce = Alerce()

def query_alerce(query=None):
    '''
    Queries ALeRCE for SNe Ia, Ib/c, and II in the last 7 days. Uses lc_classifier to find new objects.

    Parameters
    ----------
    query : dict, optional
        A dictionary of query parameters. The default is None. See https://alerce.readthedocs.io/en/latest/tutorials/ztf_api.html#query-objects for more information on query parameters.
    
    Returns
    -------
    objects : pd.DataFrame
        A DataFrame of the query results.
    '''
    sne =['SNIa','SNIbc','SNII']
    if query is None:
        query = {
        "classifier": "lc_classifier",
        "class_name": "SNIa",
        'firstmjd': [(Time.now()-7*u.day).mjd,Time.now().mjd], # and last 3 days
        'page_size':50
        }
    objects = alerce.query_objects(**query)
    return objects

class alerce_object:

    def __init__(self,aobject):
        '''
        Initializes a alerce_object from a row of the alerce query.

        Parameters
        ----------
        aobject : pd.Series
            A row from the alerce query DataFrame.
        
        Attributes
        ----------
        oid : str
            The object id.
        name : str
            The object name same as oid.
        ra : float
            The mean right ascension of the object.
        dec : float
            The mean declination of the object.
        url : str
            The url of the object in alerce.
        '''
        self.oid =aobject.oid
        self.name = self.oid
        self.ra, self.dec= aobject.meanra,aobject.meandec
        self.url = f"https://alerce.online/object/{self.oid}"


    def check_antares(self,verbose=False):
        '''
        Cross-checks if the object is in Antares. Uses ZTF object id to query ANTARES.

        Parameters
        ----------
        verbose : bool, optional
            If True, prints if the object is found in Antares. The default is False.

        Returns
        -------
        antares_object : prep.source.antares.antares_object or None
            Returns the antares_object if found in Antares, else None.
        '''
        cat = get_by_ztf_object_id(self.oid)
        if cat is None:
            return 
        else:
            if verbose:
                print('Found in Antares')
            return antares.antares_object(cat)
        
    def check_yse(self,verbose=False):
        '''
        Cross-checks if the object is in YSE. Uses Antares object to query YSE. Very hacky but alerce does not provide TNS names.

        Parameters
        ----------
        verbose : bool, optional
            If True, prints if the object is found in YSE. The default is False.

        Returns
        -------
        yse_object : prep.source.yse.yse_object or None
            Returns the yse_object if found in YSE, else None.
        '''
        return self.check_antares(verbose=verbose).check_yse(verbose=verbose)
    
    def get_probabilities(self):
        '''
        Get classifier probabilities from alerce for object.

        Returns
        -------
        probs : json dict
            The probabilities for each classifier from alerce.
        '''
        self.probs = alerce.query_probabilities(self.oid,format="json")
        return self.probs
    
    def get_top_class(self): 
        '''
        Get the top class from the classifier probabilities.

        Returns
        -------
        top_class : tuple
            The top class and probability.
        '''
        # From Patrick Aleo
        try:
            output = self.probs
        except:
            output=self.get_probabilities()

        top_class = None
        for item in output:
            if item['ranking'] == 1:
                if top_class is None or \
                    (item['classifier_name'] == 'lc_classifier' and \
                    item['classifier_version'] == 'hierarchical_rf_1.1.0') or \
                    (item['classifier_name'] == 'stamp_classifier' and \
                    item['classifier_version'] == 'stamp_classifier_1.0.4' and \
                    top_class['classifier_version'] != 'hierarchical_rf_1.1.0') or \
                    (item['classifier_name'] == 'stamp_classifier' and \
                    item['classifier_version'] == 'stamp_classifier_1.0.0' and \
                    top_class['classifier_version'] not in ['hierarchical_rf_1.1.0', 'stamp_classifier_1.0.4']):
                    top_class = item
        # print(top_class['class_name'],top_class['probability'])
        return (top_class['class_name'],top_class['probability']) if top_class is not None else ("Not_classified",None)
        
    def get_lc(self):
        '''
        Gets the light curve from alerce.

        Returns
        -------
        lc : pd.DataFrame
            The light curve from alerce. Real detections only.
        '''
        self.lc = alerce.query_detections(self.oid,format="pandas").query('has_stamp==True')
        return self.lc
    
    def plot_lc(self):
        '''
        Plot the light curve from alerce.
        '''
        try:
            self.lc
        except:
            raise ValueError('Need to run get_lc first')
        plt.errorbar(self.lc['mjd'],self.lc['magpsf'],yerr=self.lc['sigmapsf'],fmt='.')
        plt.gca().invert_yaxis()
        plt.xlabel('MJD')
        plt.ylabel('mag')
        plt.title(self.oid)
        plt.show()
    
    def salt3(self,plot=False):
        '''
        Run salt3 on the light curve from alerce. Need to run get_lc first.

        Parameters
        ----------
        plot : bool, optional
            If True, plots the light curve and the salt3 fit. The default is False.

        Returns
        -------
        salt_params : dict
            The salt3 parameters and errors.
        '''
        try:
            self.lc
        except:
            raise ValueError('Need to run get_lc first')
        self.lcm = self.lc[~self.lc[['mjd','magpsf','sigmapsf']].isna().any(axis=1)]
        # Need to figure out ztf passbands
        # ztf  - g and r
        # doing the same thing as for yse
        salt2band = []
        zpsys =[]
        mask = []
        for i in range(len(self.lcm.mjd.values)):
            if self.lcm.fid.values[i] == 1:
                salt2band.append(bandpassdict['Band: ZTF-Cam - g-ZTF'])
                mask.append(True)
            elif self.lcm.fid.values[i] == 2:
                salt2band.append(bandpassdict['Band: ZTF-Cam - r-ZTF'])
                mask.append(True)
            else:
                salt2band.append('None')
                mask.append(False)
            zpsys.append('AB')
        
        model = sncosmo.Model(source='salt3')
        fitparams = ['z', 't0', 'x0', 'x1', 'c']
        mask = mask
        salt2mjd = self.lcm['mjd'][mask].values
        mag = self.lcm['magpsf'][mask].values
        flux = 10**(-0.4*(mag-27.5))
        # flux = pdata['FLUXCAL'][mask].values
        mag_err=self.lcm['sigmapsf'][mask].values
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
            plt.title(f"{self.name} (ALeRCE)")
            plt.xlabel('MJD')
            plt.ylabel('Mag')
            p1= self.lcm.query('fid == 1')
            plt.scatter(p1.mjd,p1.magpsf,label='ZTF - g',color='g')
            p2 = self.lcm.query('fid == 2')
            plt.scatter(p2.mjd,p2.magpsf,label='ZTF - r',color='r')
            plotmjd = np.arange(result['parameters'][1]-20,result['parameters'][1]+50,0.5)
            for band in np.unique(salt2band):
                salt2flux = fitted_model.bandflux(band, plotmjd, zp=27.5,zpsys=zpsys[band == salt2band][0])
                plt.plot(plotmjd,-2.5*np.log10(salt2flux)+27.5,label=band)
            plt.gca().invert_yaxis()
            plt.grid()
            plt.legend(loc=0)


        return self.salt_params
    
    def run_ghost(self):
        '''
        Runs GHOST on the object to find potential host galaxies.

        Returns
        -------
        ghost_info : pd.DataFrame
            Dataframe with the information of the potential host galaxies.
        '''
        snCoord = SkyCoord(self.ra,self.dec,unit='deg',frame='icrs')
        with tempfile.TemporaryDirectory() as tmp:
            ghost_info = getTransientHosts(snCoord=snCoord, verbose=0,starcut='normal',savepath=tmp,ascentMatch=True)
        return ghost_info
