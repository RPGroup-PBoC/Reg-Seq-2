import pandas as pd
import numpy as np

import cmdstanpy
import polars as pl
import arviz as az

import glob as glob

import bebi103

from multiprocessing import Pool

df_gc = pd.read_csv("growth_conditions_short.csv", sep=';')
d = 2

sm_sin = cmdstanpy.CmdStanModel(stan_file='stan_models/single_gamma.stan')
sm_mix = cmdstanpy.CmdStanModel(stan_file='stan_models/mixture_gamma.stan')


df_results = pl.DataFrame()

dfwt = pl.read_csv("../../../data/twist_orders/2022-02-15_twist_order.csv")
proms = dfwt['promoter'].unique().to_list()
proms = [x for x in proms if x not in ["galEp", "ybeDp2"]]

def run_gc(gc):
    files = glob.glob(f"footprints/{gc}-*_footprints.csv")
    print(f'{df_gc.loc[df_gc.Index == gc, "Condition"].to_list()}')
    for file in files:
        df = pl.read_csv(file)
        proms = [ a for a in df['promoter'].unique().to_list() if a not in ["galEp", "ybeDp2"]]

        for prom in proms:
            print(prom)
            x = df.filter(pl.col("promoter") == prom)['mut_info'].to_list()
            x = np.convolve(x, np.ones(2*d+1)/(2*d+1), mode='valid') 
            data = dict(
                N=len(x), 
                n=x, 
                alpha_mu=-1,
                alpha_sigma=2,
                beta_mu=2,
                beta_sigma=3
            )
            rep = file.split('_')[0][-1]
            
    
            # Sample single gamma
            with bebi103.stan.disable_logging():
                samples_single = sm_sin.sample(
                    data=data,
                    chains=4,
                    iter_sampling=1000,
                )
            # Convert to ArviZ InferenceData instance
            samples_single = az.from_cmdstanpy(posterior=samples_single, posterior_predictive="n_ppc", log_likelihood="log_lik")
    
            # Sample mixture of gammas
            with bebi103.stan.disable_logging():
                samples_mix = sm_mix.sample(
                    data=data,
                    chains=1,
                    iter_sampling=5000,
                    iter_warmup=5000
                )
    
            # Convert to ArviZ InferenceData instance
            samples_mix = az.from_cmdstanpy(posterior=samples_mix, posterior_predictive="n_ppc", log_likelihood="log_lik")
    
            single_loo = az.loo(samples_single, scale="deviance")
            mix_loo = az.loo(samples_mix, scale="deviance")
            
            d_loo = mix_loo.elpd_loo - single_loo.elpd_loo
            w_single = np.exp(d_loo/2) / (1 + np.exp(d_loo/2))
            w_mix = 1 - w_single
    
            df_mcmc = bebi103.stan.arviz_to_dataframe(samples_mix)
            df_mcmc[['w', 'alpha[0]', 'alpha[1]', 'beta_[0]', 'beta_[1]', 'diverging__']].write_csv(f"inference_results/{gc}-{rep}-{prom}_inference.csv")
            
            with open("mixture_summary.txt", "a") as f:
                f.write(f"{prom}, {gc}, {rep}, {w_mix}")
if __name__ == '__main__': # Executed only on main process.
   
    with Pool(5) as p:
        result = sum(p.map(run_gc, df_gc.Index))
    
