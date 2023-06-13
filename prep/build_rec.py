import numpy as np 
import pandas as pd
import requests as req
from .auth import toku

class build_rec:
    def __init__(self, obj):
        self.obj =obj
        self.name = obj.name
        self.salt_params = obj.salt_params
        self.ra = obj.ra
        self.dec = obj.dec
        self.url = obj.url
        self.string = self.build_str()
        self.df = self.build_df()

    def build_str(self):
        z = self.salt_params['z']
        t = self.salt_params['phase']
        if t > 0:
            t = '+%i'%t
        else:
            t = '%i'%t
        x1 = self.salt_params['x1']
        c = self.salt_params['c']
        return f'<{self.obj.url}|{self.name}> z = {z:.3f}, t = {t}, x1 = {x1:.1f}, c = {c:.2f}, ra, dec = {self.ra:.3f}, {self.dec:.3f}; Found using automation'
    
    def build_df(self):
        df1=pd.DataFrame()
        df1[['name','ra','dec']] = [[self.name,self.ra,self.dec]]
        k1=self.salt_params.keys()
        k1 = list(k1)
        df1[k1]=list(self.salt_params.values())
        df1['url'] = self.url
        return df1
    
    def post(self,string=None,channel='D03BK3YKUQN'):
        if string is None:
            string = self.string
        p1=req.post('https://slack.com/api/chat.postMessage',
                 params={'channel':channel,
                         'text':string,
                         'mrkdwn':'true',
                         'parse':'none'},
                         headers={'Authorization': f'Bearer {toku}'})
        p1.raise_for_status()
        if p1.status_code == 200:
            print('Posted to Slack')

def post(string=None, channel='D03BK3YKUQN'):
        if string is None:
            raise ValueError('No string provided')
        p1=req.post('https://slack.com/api/chat.postMessage',
                 params={'channel':channel,
                         'text':string,
                         'mrkdwn':'true',
                         'parse':'none'},
                         headers={'Authorization': f'Bearer {toku}'})
        p1.raise_for_status()
        if p1.status_code == 200:
            print('Posted to Slack')