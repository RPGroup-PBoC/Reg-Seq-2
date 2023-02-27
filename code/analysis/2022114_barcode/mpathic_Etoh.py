import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/home/tom/mpathic')
from mpathic import learn_model
from multiprocessing import Pool
df = pd.read_csv("../../../data/Etoh_barcodes.csv")

df.rename(columns={"promoter": "seq"}, inplace=True)

genes = ['rspAp',
         'araBp',
         'znuCp',
         'xylAp',
         'xylFp',
         'dicCp',
         'relBp',
         'ftsKp1',
         'ftsKp2',
         'lacIp',
         'marRp',
         'dgoRp',
         'ompRp1',
         'ompRp2',
         'ompRp3',
         'ompRp4',
         'araCp']

def run(gene):
    db = "../../../data/mpathic_footprints/20221114_barcode/{}_Etoh_dataset_db".format(gene)
    _df = df.loc[df.name == gene, :]
    _df = _df.loc[_df.counts < 1000, :]
    _df = _df.loc[_df.ct_1/_df.ct > 0.02, :]
    mcmc_df = learn_model.main(
        df=_df,
        lm='IM',
        modeltype='MAT',
        LS_means_std=None,
        db=db,
        iteration=600000,
        burnin=60000,
        thin=60,
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
    mcmc_df.to_csv("footprints/Etoh_{}_mcmc_mpathic2.csv".format(gene))

pool = Pool(6)                         
pool.map(run, genes) 