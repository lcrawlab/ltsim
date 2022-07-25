
import argparse
import joblib
from pathlib import Path
import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

from sklearn import metrics

import matplotlib.pyplot as plt 

from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

def format_data(respath, pset, simnum):
    dpath = Path('{}/{}/Xsim_labeled.{}.txt'.format(respath, pset, simnum)) 
    assert dpath.exists(), 'check X path: {}'.format(dpath)
    data =  pd.read_csv(dpath, sep=' ', index_col=0)
#     mask = pd.read_csv('{}/{}/masksim_3_labeled.{}.txt'.format(respath, pset, simnum), sep=' ')
    pheno = pd.read_csv('{}/{}/ysim_labeled.{}.txt'.format(respath, pset, simnum), sep=' ', index_col=0)
    causal_snps = pd.read_csv('{}/{}/causal_snps.{}.txt'.format(respath, pset, simnum), sep=' ', header=None)
    causal_snps = list(causal_snps[0])
    pheno['y']=pheno['y']-1
    pheno['label1']=pheno['y']
    pheno['label2']=1 - pheno['label1']


    data = data.join(pheno[['label1', 'label2']])
    return data, causal_snps

def rpheno(data, rand_state=1, rand_phenotype=False):
    # Randomize phenotype
    data_input = data.copy()
    if rand_phenotype:
        data_input['label1'] = np.random.RandomState(seed=rand_state).permutation(data_input['label1'])
        data_input['label2'] = 1-data_input['label1']
        
    X = data_input.drop(['label1','label2'], axis = 1)
    y = data_input[['label1']]

    return X, y

# def train_model(mod_type, X, y, rand_state=5, num_cv=5):
#     # Setup classifier
#     if mod_type == 'log_reg':
#         mlmod = LogisticRegression(penalty='l2', solver='liblinear', class_weight='balanced', random_state=rand_state)
        
#     elif mod_type == 'rand_forest':
#         mlmod = RandomForestClassifier(class_weight='balanced', random_state=rand_state)
    
#     # CV train
#     cv_results = cross_validate(mlmod, X, y.values.ravel(), scoring='roc_auc',cv=num_cv)
    
#     return cv_results

class ModelSum:
    """Handy class to store all information about each model"""
    def __init__(self, pset, simnum):
        self.pset = pset
        self.simnum = int(simnum)

        self.data_type = None
        # self.log_reg_real = None
        # self.rand_forest_real = None
        # self.log_reg_rand = None
        # self.rand_forest_rand = None
        self.rho=None
        self.pve=None
        self.maf_frac=None
        self.obs=None
        self.k=None
        self.overlap=None
        self.class_imbalance = 0.5

        self.real_interp_auc = {}
        self.rand_interp_auc = {}

        self.real_class_auc = {}
        self.rand_class_auc = {}

        self.classification_aucs = {}

    def populate_pset_vars(self):
        # Set up simlation run
        sp = self.pset.split('-')

        if 'prop_case' in sp:
            self.class_imbalance = float(sp[sp.index('prop_case')+1])
        if 'rho' in sp:
            self.rho = float(sp[sp.index('rho')+1])
        if 'pve' in sp:
            self.pve = float(sp[sp.index('pve')+1])
        if 'maf_frac' in sp:
            self.maf_frac = float(sp[sp.index('maf_frac')+1])
        if 'obs' in sp:
            self.obs = float(sp[sp.index('obs')+1])
        if 'k' in sp:
            self.k = float(sp[sp.index('k')+1])

        if 'non_overlap_degree' in sp:
            self.overlap = False
        else:
            self.overlap = True

# get interpretability and class ROC for one split (then average ROCs over 5 splits)

def split_data(X, y, test_size=0.1, rand_state=5):
    train_x, test_x, train_y, test_y = train_test_split(X, y, test_size=test_size, random_state=rand_state,
                                                    shuffle = True, stratify = y)
    snp_labels = list(X.columns)
    return train_x, test_x, train_y, test_y, snp_labels


def train_model(mod_type, train_x, test_x, train_y, test_y, var_labels=None, causal_vars=None, var_type='snps', rand_state=5, calc_interpretability=True):
    # Setup classifier
    implemented_list = ['log_reg', 'rand_forest']
    if not calc_interpretability:
        interp_auc = None

    assert mod_type in implemented_list, 'check mod type: {}'.format(mod_type)
    if mod_type == 'log_reg':
        mlmod = LogisticRegression(penalty='l2', solver='liblinear', class_weight='balanced', random_state=rand_state)
        mlmod.fit(train_x, np.ravel(train_y)) 
        test_class_auc = mlmod.score(test_x, test_y) #confirm that score is AUC

        if calc_interpretability:
            if var_type == 'snps':
                coefs = np.abs(mlmod.coef_) # check absolute value
#           elif causa_vars == 'pathways':
            vscores = pd.DataFrame(zip(var_labels, np.ravel(coefs)), columns=['var', 'score']).set_index('var')
            _, interp_auc = calc_ROC(vscores, causal_vars)
            
    elif mod_type == 'rand_forest':
        mlmod = RandomForestClassifier(class_weight='balanced', random_state=rand_state)
        mlmod.fit(train_x, np.ravel(train_y)) 
        test_class_auc = mlmod.score(test_x, test_y) #confirm that score is AUC

        if calc_interpretability:
            if var_type == 'snps':
                coefs = np.abs(mlmod.feature_importances_) # check absolute value
                vscores = pd.DataFrame(zip(var_labels, np.ravel(coefs)), columns=['var', 'score']).set_index('var')
            # elif var_type == 'patways':

            _, interp_auc = calc_ROC(vscores, causal_vars)
        
    return mlmod, test_class_auc, interp_auc

# right now, using snps; modify so can use either SNps or genens
def calc_ROC(df_scores, true_causal): # layer: 'snp' or 'gene', list_pos: 'list_snp' or 'list_genes'
    # from Ashley
    assert true_causal is not None, 'check true_causal list'

    npips = len(df_scores)
    all_vals = list(df_scores.index)
    Pos = true_causal
    Negs = set(all_vals) - set(Pos)

    df_res = pd.DataFrame(columns=["TPR", "FPR", "FDR", "PWR"], dtype=object) 
    df_scores.sort_values(by=['score'], ascending=False, inplace=True)
    
    for i in range(1,npips+1):
        v = df_scores[0:i].index
        z = set(all_vals) -  set(v)
#         TP = len(intersection(set(v), set(Pos)))
        TP = len(set(v) & set(Pos))
        FP = len(set(v)-set(Pos))
#         TN = len(intersection(set(z), set(Negs)))
        TN = len(set(z) & set(Negs))
        FN = len(set(z)-set(Negs))
#         print(TPs,FP,TN,FN)        
        TPR = TP/(TP+FN)
            
        FPR = FP/(FP+TN) 
        FDR = FP/(FP+TP)
        PWR = 1-(1-TPR)
        df_res.loc[len(df_res)] = [TPR, FPR, FDR, PWR]
        
    # add AUC calc here    
    auc = metrics.auc(np.array(df_res['FPR']), np.array(df_res['TPR']))
    return df_res, auc

def score_cv(model_type, X, y, true_causal=None, calc_interpretability=False, cv=5):
    mods = []
    class_aucs = []
    interp_aucs = []
    
    for i in range(cv):
        train_x, test_x, train_y, test_y, snp_labels = split_data(X, y, rand_state=5+i) #check random state set
        mlmod, test_class_auc, interp_auc = train_model(model_type, train_x, test_x, train_y, test_y, snp_labels, true_causal, calc_interpretability=calc_interpretability) 

        # print(i, test_class_auc, interp_auc)
        mods.append(mlmod)
        class_aucs.append(test_class_auc)
        interp_aucs.append(interp_auc)
        
    return mods, class_aucs, interp_aucs

parser = argparse.ArgumentParser()
parser.add_argument("--pset", required=True)
parser.add_argument("--simnum", required=True)
parser.add_argument("--respath", required=True)
parser.add_argument("--savepath", required=True)
#python heritability_model.py --pset "non_overlap_degree-4-ind-1e+06-tot_snp_sim-5000-frac_causal-0.1-k-0.5-obs-500-maf_frac-0.05-pve-1-rho-0" --simnum 1 --respath "../data/simulation"

if __name__ == '__main__':
    args = parser.parse_args()
    savepath = args.savepath

    ms = ModelSum(args.pset, args.simnum)
    ms.populate_pset_vars()

    # set up simulation run
    if (args.pset.split('-')[0]=='ind')|(args.pset.split('-')[0]=='non_overlap_degree'): #clean?
        ms.data_type='simulation'
        calc_interp = True
        dpkl, csnps = format_data(respath=args.respath, pset=args.pset, simnum=args.simnum)

    else:
        ms.data_type = 'real'
        calc_interp = False
        csnps = None
        rdpath = Path('{}/{}/data_pd.{}.pkl'.format(args.respath, args.pset, args.simnum)) 
        assert rdpath.exists(), 'check data_pd path: {}'.format(rdpath)
        dpkl = joblib.load(rdpath)

    # Fix this
    assert dpkl.dropna(how='any').shape == dpkl.shape, print('NaN: {}, {}'.format(args.pset, args.simnum))
    mod_list = ['log_reg', 'rand_forest']

    cX, cY = rpheno(dpkl, rand_phenotype=False)
    for model in mod_list:
        mmod, mcauc, miauc = score_cv(model, cX, cY, true_causal=csnps, calc_interpretability=calc_interp)
        ms.real_class_auc[model]=mcauc
        ms.real_interp_auc[model]=miauc

    # rand pheno
    cX, cY = rpheno(dpkl, rand_phenotype=True)
    for model in mod_list:
        mmod, mcauc, miauc = score_cv(model, cX, cY, true_causal=csnps, calc_interpretability=calc_interp)
        ms.rand_class_auc[model]=mcauc
        ms.rand_interp_auc[model]=miauc

    full_save = '{}/mod_{}-sim-{}.pkl'.format(savepath, args.pset, args.simnum)
    assert Path(savepath).exists(), 'check savepath: {}'.format(savepath)
    joblib.dump(ms, full_save)
