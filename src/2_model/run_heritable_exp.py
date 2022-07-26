#!/usr/bin/python

import subprocess
import glob
from pathlib import Path

savepath =  '../results/heritability_072522'

data_path = '../data/simulation'

data_path = '/scratch/users/aconard/sims/k_obs_pve/overlap_sim'
# data_path = '/scratch/users/aconard/sims/k_obs_pve/non_overlap_sim'
filelist = list(Path(data_path).rglob('Xsim_labeled.*.txt'))

# data_path = '../data/als_snps223095'
# data_path = '../data/diabetes_snps72820'
# filelist = list(Path(data_path).rglob('data_pd.*.pkl'))

assert len(filelist)>0, 'check path: {}'.format(data_path)
Path(savepath).mkdir(parents=False, exist_ok=True)

for path in filelist:
	simnum = path.name.split('.')[-2]
	pset = path.parent.name

	full_file = '{}/mod_{}-sim-{}.pkl'.format(savepath, pset, simnum)
	if not Path(full_file).exists():
		subprocess.run(["sbatch", "parallel_heritability.sbatch", "-a", pset, "-b", simnum, "-c", data_path, "-d", savepath])
	break