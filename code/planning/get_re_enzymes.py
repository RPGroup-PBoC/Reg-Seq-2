from Bio.Restriction import *
import itertools
import pandas as pd
import git

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

enzyme_list = [(enzyme, enzyme.site) for enzyme in itertools.islice(CommOnly, len(CommOnly))]
pd.DataFrame(columns=["enzyme", "site"], data=enzyme_list).to_csv(f"/{homedir}/data/re_enzyme_list.csv", index=False)

