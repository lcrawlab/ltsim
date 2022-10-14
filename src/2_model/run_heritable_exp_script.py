#!/usr/bin/python

'''
Example Usage: python run_heritable_exp_script.py --data_path '../../data/simulation' --savepath '../../results/hscript_test'
Arguments:
    data_path: path to folder with simulated data (contains individual subfolders, each subfolder corresponding to one simulation parameter set)
    savepath: path to save trained model results
    use_slurm: optional flag to use slurm scheduler if running on cluster; if using slurm, may need to adjust SBATCH params depending on your cluster setup
'''
import subprocess
import argparse

import glob
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--savepath", required=True)
parser.add_argument("--data_path", required=True)
parser.add_argument("--use_slurm", type=bool)

if __name__ == '__main__':
	args = parser.parse_args()
	savepath = args.savepath
	data_path = args.data_path
	use_slurm = args.use_slurm
	filelist = list(Path(data_path).rglob('Xsim_labeled.*.txt'))

	assert len(filelist)>0, 'check path: {}'.format(data_path)
	Path(savepath).mkdir(parents=True, exist_ok=True)

	for path in filelist:
		simnum = path.name.split('.')[-2]
		pset = path.parent.name

		full_file = '{}/mod_{}-sim-{}.pkl'.format(savepath, pset, simnum)
		if use_slurm:
			call = 'sbatch'
		else:
			call = 'sh'

		if not Path(full_file).exists():
			subprocess.run([call, "parallel_heritability.sbatch", "-a", pset, "-b", simnum, "-c", data_path, "-d", savepath])
		break