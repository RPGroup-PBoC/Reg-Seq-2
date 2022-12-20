import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/home/tom/mpathic')
from mpathic import learn_model

df = pd.read_csv("../../../../RegSeq/data/sequencing_data/LB_data.csv", index_col=0)
df = df.loc[df['gene'] == 'rspA', :]
df.rename(
    columns={
        'counts_RNA': 'ct_1',
        'counts_DNA': 'ct_0',
    },
    inplace=True
)
df['ct'] = df['ct_0'] + df['ct_1']
df = df[['ct', 'ct_0', 'ct_1', 'seq']]
df.dropna(inplace=True)
df.reset_index(drop=True, inplace=True)



db = "../../../data/mpathic_footprints/20221114_barcode/rspAp_LB_dataset_db"

mcmc_df = learn_model.main(
    df=df,
    lm='IM',
    modeltype='MAT',
    LS_means_std=None,
    db=db,
    iteration=300000,
    burnin=100000,
    thin=1000,
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
mcmc_df.to_csv("rspA_mcmc_mpathic_regseq_data2.csv")