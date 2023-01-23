# LT-Sim

A liability threshold linear mixed model simulation tool for case-control studies.

### Toggles
For each dataset there are 100 replicates. For those replicates, the following are constant: $k$, number observations (i.e. samples), number SNPs, PVE, and MAF (0.001 - 0.5). 

### Replicates
For each replicate, change SNPs in groups 1-3.

### Simulation types 

We want to understand the influence of effect sizes on the ability to detect additive and epistatic effects. Thus we have these types of simulations:

(1) $c$: number of causal SNPs in groups 1, 2, and 3. Groups 1 and 2 have 10 each, rest go to group 3. Sparse (1% of SNP-sets are enriched for trait) vs. polygenic (10% of SNP-sets are enriched)

(2) $H^2$: broad-sense heritability that is additive.

(3) $\rho$: the contribution of epistatic effects (sampled from normal distribution).

(4) $f$: fraction of class imbalance.


### Steps
- Population size:
    - [x] ind = 1e6 

- Minor allel frequency:
    - [x] maf_frac = 0.05
    - [ ] maf_frac = 0.0001
    
- Sampling 3000 total observations (in real data als_cases = 2100; als_controls = 790):
    - [x] obs_cases = 1500, obs_controls = 1500 (i.e. balanced)
    - [ ] obs_cases = 1000, obs_controls = 2000
    - [ ] obs_cases = 2000, obs_controls = 1000
    - [ ] obs_cases = 2100, obs_controls = 900 (~0.7 skew as in real data)
    
- Simulation according to liability threshold model (5 per 100K people in U.S. pop.): 
	- [ ] $k$ = 0.01
	- [x] $k$ = 0.1
    - [ ] $k$ = 0.5

- Simulate fraction causal pathways, genes, and eQTLs (i.e. frac_causal):
    - [ ] $c$ = 0.01
    - [ ] $c$ = .65 (from ALS_EUR GWAS at p-val 0.05)https://www.projectmine.com/research/download-data/
    - [x] $c$ = .1

- Phenotypic variance explained (i.e. broad-sense heritability $H^2$)
    - [ ] $pve$ = 0.3 
    - [ ] $pve$ = 0.4 
    - [x] $pve$ = 0.6

- $H^2$ fraction from marignal (i.e. additive) effects
    - [ ] $\rho$ = 0
    - [x] $\rho$ = 0.5
        
### Data Generation
- Used liability threshold model to generate 1 million people.
- Sampled 1000 at 0.1% prevalence ($k$) with 6.5% causal SNPs.
- Simulated causal pairs of SNPs.
- Deemed genes associated with causal SNPs to be causal
- Used real mask for SNPs to genes (GWAS Catalog and UCSC table browser) and subsampled for simulations
- Used real mask (CP) for genes to pathways and subsampled for simulations

### To do:
X - done, O - in progress
- [x] Complete simulation set-up with causal SNP pairs from pathways
- [o] Set up 'Steps' for data generation for 100 replicates
- [x] Set up scripts to run for BANNs
- [x] Set up scripts to run sparsenn model
- [x] Discuss best comparison to make with sparsenn model
- [ ] Set up script to generate ROC curve input for sparsenn model
- [x] Set up ROC curve scripts
- [ ] Redo with more complex causal SNP generation scheme.

### To run (need to update *):
After you have started the `multio_case_ctrl` conda environment, 

`INPUT_DATA="/data/compbio/aconard/microsoft_research/mit_collab/data/"`
`SIM_DATADIR="/data/compbio/aconard/microsoft_research/mit_collab/sim_results/"`
`BANNS_RESULTS="/nbu/compbio/aconard/microsoft_research_collab/banns_results/01122022/"` 

0) Run `./sparsenn/simulations/lt_sims/run_sparsenn.ipynb` to generate mask1 (SNPs to genes) and mask2 (genes to pathways).

1) * `Rscript ./sparsenn/simulations/lt_sims/lt_simulations_updated_3.r $k $f $c $H2 $\rho $INPUT_DATA $SIM_DATADIR`

2) `python ./sparsenn/simulations/run_banns.py $SIM_DATADIR $BANNS_RESULTS`

3) `./sparsenn/src/run_sparsenn.sh $SIM_DATADIR`

4) * `python ./sparsenn/simulations/roc_calc.py --results-folder=$RES --two-layers=1 --tp-snps-genes=$TPs_regex`
