from prep.source import *
from astropy.time import Time
import tempfile
import os
import glob
import sqlite3
from prep.build_rec import post as pst
from prep.build_rec import build_rec as bs
import time

def check_post(name,cur):
    if cur == None:
        return False
    if cur.execute(f'SELECT * FROM posted_names WHERE name LIKE "{name}"').fetchone() != None:
        print(f'{name} already posted!')
        return True
    else:
        return False
    

    # if lines is None:
    #     return False
    
    # if name in lines:
    #     print(f'{name} already posted!')
    #     return True
    # else:
    #     return False


def run(sources=['antares','alerce','yse_yaf'],post=True,name_file='posted_names_test.db'):
    if os.path.exists(name_file):
        print(f'using existing posted names file: {name_file}') # checking if the names file exists, if it doesn't, will not check for duplicates and will not write new names to the file
        # with open(name_file,'r') as f:
        #     lines = f.read().split('\n')
        #     f.close()
        conn = sqlite3.connect(name_file)
        cur = conn.cursor()
    else:
        print(f'WARNING: posted names file not found, will not check for duplicates')
        cur = None
    ps= []
    names = []
    if 'alerce' in sources:
        aq=alerce_api.query_alerce()
        for i in range(aq['oid'].values.size):
            try:
                ao=alerce_api.alerce_object(aq.iloc[i])
                if check_post(ao.name,cur):
                    continue
                ao.get_lc()
                ao.salt3()
                bss =bs(ao)
                names.append(bss.name)
                ps.append(bss.string)
            except Exception as e:
                print(f'failed on {aq.iloc[i].oid}\n{e}')
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
                # print(ao.name)
                if check_post(ao.name,cur):
                    continue
                ao.get_lc()
                ao.salt3()
                bss = bs(ao)
                names.append(bss.name)
                ps.append(bss.string)
            except Exception as e:
                print(f'failed on {locus.properties["ztf_object_id"]}'+'\n'+str(e))
                
            j+=1
    
    if 'yse_yaf' in sources:
        # maybe introduce temp file to save csv
        with tempfile.TemporaryDirectory() as td:
            yq = yse.young_and_fast(td)
            gp = glob.glob(os.path.join(td,'YSE_*.csv'))[0]
            qd=pd.read_csv(gp,skiprows=1,names=['name','classification','first_detection','latest_detection','number_of_detection','group_name'])
        qd.sort_values(by='number_of_detection',ascending=False,inplace=True)
        f4=qd.query('number_of_detection>=4')
        for name in f4.name.values:
            try:
                y1=yse.yse_object(name)
                if y1.in_field:
                    note = f'Found in YSE fields: {", ".join(y1.fields)}'
                else:
                    note = None
                if check_post(y1.name,cur):
                    continue
                y1.get_lc()
                y1.salt3()
                bss = bs(y1,note=note)
                names.append(bss.name)
                ps.append(bss.string)
            except Exception as e:
                print(f'failed on {name}\n{e}')
                continue
    if 'yse_hst' in sources:
        with tempfile.TemporaryDirectory() as td:
            yq = yse.possible_hst(td)
            gp = glob.glob(os.path.join(td,'YSE_*.csv'))[0]
            qd=pd.read_csv(gp)
        qd.sort_values(by='disc_date',ascending=False,inplace=True)
        f4=qd
        for name in f4.name.values:
            try:
                y1=yse.yse_object(name)
                if y1.in_field:
                    note = f'Found in YSE fields: {", ".join(y1.fields)}'
                else:
                    note =None
                if check_post(y1.name,cur):
                    continue
                y1.get_lc()
                y1.salt3()
                bss = bs(y1,note=note)
                names.append(bss.name)
                ps.append(bss.string)
            except Exception as e:
                print(f'failed on {name}\n{e}')
                continue

    if post:
        if cur is None: # if the file doesn't exist, don't write to it
            pass
        else:
            cur.executemany('INSERT INTO posted_names VALUES (?)',[(n,) for n in names])
            # bnames = '\n'.join(names)
            # with open(name_file,'a') as f:
            #     f.write(bnames)
            #     f.close()
            conn.commit()
            conn.close()
        ps = '\n'.join(ps)
        pst(ps,
            # channel='C05E9AJ18HG'
            )

    return 0


def run_sched(sep= 86400): # 24 hours in seconds
    while True:
        t1 = time.monotonic()
        run()
        t2 = time.monotonic()
        td = t2-t1
        time.sleep(sep-td) # calculation offset so it runs at the same time every day
    return 0


if __name__ == '__main__':
    print('running')
    run_sched()

def main():
    print('running')
    run_sched()
    return 0