import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/home/tom/mpathic')
from mpathic import learn_model

# Custom analysis
df_new_LB_DNA = pd.read_csv(
    "../../../data/LB_heatshock_bc_by_gc/LB_DNA_113_identified.txt", 
    delim_whitespace=True, 
    names=["ct", "tag", "seq", "mapping_count", "name", "nmut"]
)

df_new_LB_RNA = pd.read_csv(
    "../../../data/LB_heatshock_bc_by_gc/LB_RNA_113_identified.txt", 
    delim_whitespace=True, 
    names=["ct", "tag", "seq", "mapping_count", "name", "nmut"]
)

df_new_LB = df_new_LB_DNA.merge(df_new_LB_RNA, on=['tag', 'seq', 'mapping_count', 'name', 'nmut'], how='outer', suffixes=('_0', '_1'))
df_new_LB.fillna(0, inplace=True)
df_new_LB = df_new_LB.loc[df_new_LB['ct_0'] < 1000, :]
df_new_LB = df_new_LB.loc[df_new_LB['ct_1'] < 1000, :]
df_new_LB = df_new_LB.loc[df_new_LB['nmut'] < 30, :]
df_new_LB = df_new_LB.loc[df_new_LB['mapping_count'] > 1, :]
df_new_LB = df_new_LB.loc[df_new_LB["name"] == "rspA", :]
df_new_LB.reset_index(inplace=True, drop=True)

# Fix type
df_new_LB["ct_1"] = [int(x) for x in df_new_LB.ct_1]
df_new_LB["ct_0"] = [int(x) for x in df_new_LB.ct_0]
df_new_LB["ct"] = df_new_LB.ct_1 + df_new_LB.ct_0
df_new_LB = df_new_LB[['ct', 'ct_0','ct_1', 'seq']]



db = "../../../data/mpathic_footprints/20221114_barcode/__rspAp_LB_dataset_db"

mcmc_df = learn_model.main(
    df=df_new_LB,
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
mcmc_df.to_csv("rspA_mcmc_mpathic_regseq_data2.csv")
