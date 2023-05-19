import requests as req
import pandas as pd
import os
from requests.auth import HTTPBasicAuth
from ..auth import *
from io import StringIO
import matplotlib.pyplot as plt
import sncosmo
from .bandpassdict import *
from astropy.table import Table
import numpy as np
from astropy.time import Time
import re
# from tns_search_download_csv import search_tns
from astro_ghost.ghostHelperFunctions import getTransientHosts
import tempfile
from astropy.coordinates import SkyCoord

# today = Time.now()

def young_and_fast(path=None):
    today = Time.now()
    q1=req.get('https://ziggy.ucolick.org/yse/explorer/254/download',auth=HTTPBasicAuth(login, password))
    q1.raise_for_status()
    qt1=q1.text
    time =today.to_value('datetime').strftime('%Y%m%d_%H%M%S')
    if path is None:
        path = os.getcwd()
    p1 = os.path.join(path,f'YSE_{time}.csv')
    with open(p1,'w') as f:
        f.write(qt1)
        f.close()
    return 0

class yse_object:

    def __init__(self,name,data=None):
        # Name format should be "20xxabc" but can have SN in front
        self.name =name
        if 'SN' in self.name:
                self.ns = self.name.split('SN')[-1].strip()
        else:
                self.ns = self.name
        self.data = data
        if self.check_pz()==False:
            print(f'{self.name} not found on YSE-PZ')
            raise ValueError('Not found on YSE-PZ')
            return
        else:
             self.ebool = True
             self.url = self._gen_YSE_PZ_url()

    def query_yse_object(self):
        r1 = req.get(f'https://ziggy.ucolick.org/yse/api/transients/?name={self.ns}',auth=HTTPBasicAuth(login, password)).json()['results']
        
        return r1


    def check_pz(self):
        r1 = req.get(f'https://ziggy.ucolick.org/yse/api/transients/?name={self.ns}',auth=HTTPBasicAuth(login, password)).json()['results']
        if r1 == []:
            return False
        else:
            self.ra = r1[0]['ra']
            self.dec = r1[0]['dec']
            return r1!=[]
    
    def _gen_YSE_PZ_url(self):
        # if self.ebool == None:
        #     self.check_PZ()
        # assert self.ebool
        url = f'https://ziggy.ucolick.org/yse/transient_detail/{self.ns}'
        return url
    
    def get_lc(self):
        
        r=req.get(f'https://ziggy.ucolick.org/yse/download_photometry/{self.ns}/',auth=HTTPBasicAuth(login, password))
        r.raise_for_status()
        li=r.text.split('\n')
        key1=[]
        val1=[]
        for i in range(len(li)):
            if li[i].startswith('#'):
                vs=li[i].removeprefix('#')
                key1.append(vs.split(':')[0].strip())
                val1.append(vs.split(re.findall('[A-Z]: ',vs)[0])[-1])
        header = dict(zip(key1,val1))
        self.header = header
        # DISC_DATE=Time(header['DISC_DATE']).mjd
        # self.dd = DISC_DATE
        # NON_DETECT_LIMIT=header['NON_DETECT_LIMIT']
        # NON_DETECT_BAND=header['NON_DETECT_BAND']
        # NON_DETECT_FLT=NON_DETECT_BAND.split(' - ')[-1].strip()
        # qstring = f'DQ!="Bad" & FLT != "Unknown" & MAG<{NON_DETECT_LIMIT} & MJD>{NON_DETECT_MJD} & FLT != "{NON_DETECT_FLT}"'
        snr_cut = 3
        qstring = f' FLT != "Unknown"  & MAGERR<{snr_cut}'
        self.lc_data=pd.read_csv(StringIO('\n'.join(li)),sep='\s+',comment='#')
        self.pdata=pd.read_csv(StringIO('\n'.join(li)),sep='\s+',comment='#').query(qstring).dropna()
    
    def salt3(self,plot=False):
        try:
            self.pdata
        except:
            raise ValueError('Need to run get_lc first')
        salt2band = np.array([])
        zpsys = np.array([])
        mask  = []
        # bs =[]
        for ins,flt in zip(self.pdata['INSTRUMENT'],self.pdata['FLT']):
            if ins == 'HKO':
                ins = 'ACAM1'
            # bs.append()
            if f"Band: {ins} - {flt}" not in bandpassdict.keys():
                
                mask.append(False)
                salt2band = np.append(salt2band,str(np.nan))
                zpsys = np.append(zpsys,str(np.nan))
            else:
                
                mask.append(True)
                salt2band = np.append(salt2band,bandpassdict[f"Band: {ins} - {flt}"])
                if 'bessell' in salt2band[-1]:
                    zpsys = np.append(zpsys,'Vega')
                else:
                    zpsys = np.append(zpsys,'AB')
        model = sncosmo.Model(source='salt3')

        # if sn.Redshift:
        #     model.set(z=sn.Redshift); fitparams = ['t0', 'x0', 'x1', 'c']
        try:
            if self.header['REDSHIFT']!='' or self.header['REDSHIFT']==None:
                model.set(z=float(self.header['REDSHIFT']))
                fitparams = [ 't0', 'x0', 'x1', 'c']
            else:
                pass
        except:
            fitparams = ['z', 't0', 'x0', 'x1', 'c']
        mask = mask & (self.pdata['FLUXCAL']<1e10) #& (pdata['MJD'].values>60052)
        # model.set(z=sn['Host Redshift'])
        salt2mjd = self.pdata['MJD'][mask].values
        mag = self.pdata['MAG'][mask].values
        flux = 10**(-0.4*(mag-27.5))
        # flux = pdata['FLUXCAL'][mask].values
        mag_err=self.pdata['MAGERR'][mask].values
        # mag_err[(mag_err < 0.01)|(mag_err ==None)]=0.01
        fluxerr = flux*mag_err*0.4*np.log(10)
        # fluxerr = pdata['FLUXCALERR'][mask].values
        zp = np.array([27.5]*len(salt2band))[mask]
        zpsys = zpsys[mask]
        salt2band = salt2band[mask]

        data = Table([salt2mjd,salt2band,flux,fluxerr,zp,zpsys],names=['mjd','band','flux','fluxerr','zp','zpsys'],meta={'t0':salt2mjd[flux == np.max(flux)]})

        # pkguess = np.atleast_1d(salt2mjd[flux/fluxerr > 3][flux[flux/fluxerr >3] == np.max(flux[flux/fluxerr >3])])
        # if len(pkguess):
        #     pkguess = pkguess[0]
        #     data = data[(salt2mjd > pkguess-20) & (salt2mjd < pkguess+40)]
        # print(salt2mjd)
        result, fitted_model = sncosmo.fit_lc(
                    data, model, fitparams,
                    bounds={'t0':(salt2mjd[flux == np.max(flux)]-10, salt2mjd[flux == np.max(flux)]+10),
                            'z':(0.0,0.7),'x1':(-3,3),'c':(-0.3,0.3)})  # bounds on parameters (if any)
        
        today = Time.now().mjd
        lcphase = today-result['parameters'][1]
        lcp =lcphase
        if plot:
            print('points fitted =',len(salt2mjd))
            print('chisq =',result['chisq'])
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
            plt.title(f"{self.ns} (YSE-PZ)")
            plt.xlabel('MJD')
            plt.ylabel('Mag')
            # plt.plot(pdata['MJD'],pdata['MAG'],marker='o',linestyle='none',label='not fitted points')
            j=0
            pdm =~self.pdata[['INSTRUMENT','FLT']].duplicated().values
            for i in range(len(self.pdata.INSTRUMENT.values[pdm])):
                wog = self.pdata
                ins = wog.INSTRUMENT.values[pdm][i]
                flt = wog.FLT.values[pdm][i]
                qs=f'INSTRUMENT == "{ins}" & FLT == "{flt}"'
                dq = self.pdata.query(qs)
                plt.plot(dq.MJD,dq.MAG,marker='o',linestyle='none',label=f"{ins} - {flt}")
                j+=1
            # plt.plot(result)
            plotmjd = np.arange(result['parameters'][1]-20,result['parameters'][1]+50,0.5)
            for band in np.unique(salt2band):
                salt2flux = fitted_model.bandflux(band, plotmjd, zp=27.5,zpsys=zpsys[band == salt2band][0])
                
                plt.plot(plotmjd,-2.5*np.log10(salt2flux)+27.5,label=band)
            plt.gca().invert_yaxis()
            plt.grid()
            plt.legend(loc=0)
    

    def run_ghost(self):
        snCoord = SkyCoord(self.ra,self.dec,unit='deg',frame='icrs')
        with tempfile.TemporaryDirectory() as tmp:
            ghost_info = getTransientHosts(snCoord=snCoord, verbose=0,starcut='normal',savepath=tmp,ascentMatch=True)
        return ghost_info

    def parsnip(self):
        # this will run parsnip on the light curve
        return
