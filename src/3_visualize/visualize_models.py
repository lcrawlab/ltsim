#!/usr/bin/python
'''
Example Usage: python visualize_models.py --trained_model_path '../../results/hscript_test' --save_plot_path '../../results/hvis'
Description: Create classification and interpretability line plots, across varied prevalences
Arguments:
    trained_model_path: path to folder with saved trained model pkl files
    save_plot_path: path to folder where plot images should be saved
'''

import joblib
import argparse
from pathlib import Path
import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt 


import sys
sys.path.append('../2_model/')
from heritability_model_interpretable import ModelSum


def load_trained_mods(res_path):
	print(res_path)
	res_df = pd.DataFrame()
	for file in list(glob.glob(res_path + '/*.pkl')):
		ms = joblib.load(file)
		ms.populate_pset_vars()
		cdf = pd.DataFrame({'pset':ms.pset, 'simnum':ms.simnum, 'data_type':ms.data_type,
							'rho':ms.rho, 'pve':ms.pve,'maf_frac':ms.maf_frac, 'obs':ms.obs, 'k':ms.k, 
							'overlap':ms.overlap, 'class_imbalance':ms.class_imbalance}, index={'{}.{}'.format(ms.pset, ms.simnum)})
		for mod in ms.real_class_auc.keys():
			cdf['{}_real_class_auc'.format(mod)]=np.mean(ms.real_class_auc[mod])
			cdf['{}_rand_class_auc'.format(mod)]=np.mean(ms.rand_class_auc[mod])
			if ms.real_interp_auc[mod][0] is not None:
				cdf['{}_real_interp_auc'.format(mod)]=np.mean(ms.real_interp_auc[mod])
				cdf['{}_rand_interp_auc'.format(mod)]=np.mean(ms.rand_interp_auc[mod]) 
		res_df = res_df.append(cdf)

	return res_df

def create_plot(res_df, save_plot_path, over=False, class_imbalance=0.5):
	"Create baseline lineplots for logistic regression and random forest, with varying k values, heritability, and number of observations"
	
	res_df = res_df.fillna(value=np.nan)
	rsubdf = res_df[(res_df['data_type']=='simulation')&(res_df['class_imbalance']==class_imbalance)].groupby(['pset'])

	# Calculate confidence interval for replicates
	mean_df = rsubdf.mean()
	ci_df = 1.96 * rsubdf.std()/np.sqrt(rsubdf.count())

	klist = mean_df['k'].unique()

	for k in klist:
		mlabel = {'log_reg': 'Logistic Regression', 'rand_forest':'Random Forest', 'class': 'Classification', 'interp':'Interpretability'}

		fig, axs = plt.subplots(2,2, figsize=(15,8), sharex=True, sharey=True)
		_ =fig.suptitle('Classification and Interpretability - Prevalence: {}'.format(k), fontsize=25)
		mean_df_k = mean_df[(mean_df['k']==k) & (mean_df['overlap']==over) & (mean_df['class_imbalance']==class_imbalance)].copy()
		mean_df_k = mean_df_k.sort_values(by=['obs', 'pve'])

		for i, task in enumerate(['class', 'interp']):
			for j, mod in enumerate(['log_reg', 'rand_forest']):
				for obs in mean_df_k['obs'].unique():
					mean_df_k_o = mean_df_k[(mean_df_k['obs']==obs)]
					ci_df_k_o = ci_df.loc[mean_df_k_o.index]

					if pd.notnull(ci_df_k_o['simnum']).all(): #remove if can't calculate CI

						_ = axs[i,j].plot(mean_df_k_o['pve'], mean_df_k_o['{}_rand_{}_auc'.format(mod, task)], '--', color='gray', label='Random Phenotype')
						_ = axs[i,j].plot(mean_df_k_o['pve'], mean_df_k_o['{}_real_{}_auc'.format(mod, task)], '-o', label='Obs:{}'.format(int(obs)))

						ci_up = mean_df_k_o['{}_real_{}_auc'.format(mod, task)] + ci_df_k_o['{}_real_{}_auc'.format(mod, task)] 
						ci_low = mean_df_k_o['{}_real_{}_auc'.format(mod, task)] - ci_df_k_o['{}_real_{}_auc'.format(mod, task)] 
						_ = axs[i,j].fill_between(mean_df_k_o['pve'], ci_up, ci_low, alpha=.05)


						ci_up = mean_df_k_o['{}_rand_{}_auc'.format(mod, task)] + ci_df_k_o['{}_rand_{}_auc'.format(mod, task)] 
						ci_low = mean_df_k_o['{}_rand_{}_auc'.format(mod, task)] - ci_df_k_o['{}_rand_{}_auc'.format(mod, task)] 
						_ = axs[i,j].fill_between(mean_df_k_o['pve'], ci_up, ci_low, color='gray', alpha=.02)


						_ = axs[i,j].set_title('{}: {}'.format(mlabel[task], mlabel[mod]), fontsize=20)
						if (i==0)&(j==0):
							format_axs(axs[i,j], legend=True)
						else:
							format_axs(axs[i,j], legend=False)
			fig.tight_layout()
			
			Path(save_plot_path).mkdir(parents=True, exist_ok=True)
			fspath = '{}/fig_class_interp_pve_{}.pdf'.format(save_plot_path, k)
			fig.savefig(fspath)

def format_axs(ax, xlabel='Heritability (pve)', ylabel='Cross. Val. AUC', legend=False, fontmin=10):
	if legend:
		handles, labels = ax.get_legend_handles_labels()
		by_label = dict(zip(labels, handles))
		ax.legend(by_label.values(), by_label.keys(), fontsize=fontmin, loc='upper left',frameon=False)

	ax.set_xlabel(xlabel, fontsize=fontmin+5)
	ax.set_ylabel(ylabel, fontsize=fontmin+5)

	ax.set_ylim([0.3, 1])
	ax.xaxis.set_tick_params(labelsize=fontmin)
	ax.yaxis.set_tick_params(labelsize=fontmin)


parser = argparse.ArgumentParser()
parser.add_argument("--trained_model_path", required=True)
parser.add_argument("--save_plot_path", required=True)

if __name__ == '__main__':
	args = parser.parse_args()

	rdf = load_trained_mods(args.trained_model_path)
	assert(rdf.shape[0]>0), 'check model path'
	create_plot(rdf, args.save_plot_path)