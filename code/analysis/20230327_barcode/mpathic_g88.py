import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/home/tom/mpathic')
from mpathic import learn_model
from multiprocessing import Pool
df = pd.read_csv("g88_barcodes.csv")

df.rename(columns={"promoter": "seq"}, inplace=True)

genes = ['rspAp',
         'araBp',
         'znuCp',
         'xylAp',
         'xylFp']

def run(gene):
    db = "../../../data/mpathic_footprints/20230327_barcode/{}_g88_dataset_db".format(gene)
    _df = df.loc[df.name == gene, :]
    _df = _df.loc[_df.ct < 1000, :]
    #_df = _df.loc[_df.ct_1/_df.ct > 0.02, :]
    mcmc_df = learn_model.main(
        df=_df,
        lm='IM',
        modeltype='MAT',
        LS_means_std=None,
        db=db,
        iteration=300000,
        burnin=30000,
        thin=30,
        runnum=0,
        initialize='rand',
        start=0,
        end=None,
        foreground=1,
        background=0,
        alpha=0,
        pseudocounts=1,
        test=False,
        drop_library=False,
        verbose=True,
    )
    mcmc_df.to_csv("footprints/g88_{}_mcmc_mpathic2.csv".format(gene))

pool = Pool(6)                         
pool.map(run, genes) 
