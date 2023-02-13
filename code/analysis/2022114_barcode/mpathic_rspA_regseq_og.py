import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/home/tom/mpathic')
from mpathic import learn_model

df = pd.read_csv("./data_files/rspA_LB_regseq", delim_whitespace=True)
df['seq'] = [x[0:160] for x in df['seq']]
db = "../../../data/mpathic_footprints/20221114_barcode/_rspAp_LB_dataset_db"

mcmc_df = learn_model.main(
    df=df,
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
mcmc_df.to_csv("rspA_mcmc_mpathic_regseq_og.csv")

