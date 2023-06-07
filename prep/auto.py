from prep.source import *
from astropy.time import Time
import tempfile
import os
import glob
from prep.build_rec import post as pst
from prep.build_rec import build_rec as bs
import time




def run(sources=['antares','alerce','yse'],post=True):
    ps= []
    if 'alerce' in sources:
        aq=alerce_api.query_alerce()
        for i in range(aq['oid'].values.size):
            try:
                ao=alerce_api.alerce_object(aq.iloc[i])
                ao.get_lc()
                ao.salt3()
                ps.append(bs(ao).string)
            except:
                print(f'failed on {aq.iloc[i].oid}')
                continue
    
    if 'antares' in sources:
        today = Time.now()
        query = (
        Search()
        .filter("range", **{"properties.num_mag_values": {"gte": 4, "lte": 100}})
        .filter('range',**{'properties.oldest_alert_observation_time': {"gte":(today-7*u.day).mjd}})
        #.filter("term", tags="extragalactic")
        # .filter("term", tags="high_amplitude_transient_candidate")
        .to_dict()
        )
        aq=antares.query_antares(query)
        j=0
        for locus in aq:
            if j>50:
                break
            try:
                ao=antares.antares_object(locus)
                print(ao.name)
                ao.get_lc()
                ao.salt3()
                ps.append(bs(ao).string)
            except:
                print(f'failed on {locus.properties["ztf_object_id"]}')
                
            j+=1
    
    if 'yse' in sources:
        # maybe introduce temp file to save csv
        with tempfile.TemporaryDirectory() as td:
            yq = yse.young_and_fast(td)
            gp = glob.glob(os.path.join(td,'YSE_*.csv'))[0]
            qd=pd.read_csv(gp,skiprows=1,names=['name','classification','first_detection','latest_detection','number_of_detection','group_name'])
        qd.sort_values(by='number_of_detection',ascending=False,inplace=True)
        f4=qd.query('number_of_detection>=5')
        for name in f4.name.values:
            try:
                y1=yse.yse_object(name)
                y1.get_lc()
                y1.salt3()
                ps.append(bs(y1).string)
            except:
                print(f'failed on {name}')
                continue

    if post:
        ps = '\n'.join(ps)
        pst(ps)

    return 0


def run_sched(sep= 86400): # 24 hours in seconds
    while True:
        t1 = time.monotonic()
        run()
        t2 = time.monotonic()
        td = t2-t1
        time.sleep(sep-td)
    return 0


if __name__ == '__main__':
    run_sched()