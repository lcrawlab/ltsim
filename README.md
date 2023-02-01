# LT-Sim

A liability threshold linear mixed model simulation tool for case-control studies.

### Overview
A major focus in human genetics is the development of computational models that can both accurately predict disease risk and identify mutations that influence phenotypic variation. Currently, it is a challenge for researchers to effectively assess the power of their methods due to the lack of scalable simulation tools that can emulate complex genetic interactions and realistic trait architectures. This is particularly difficult for simulating synthetic binary phenotypes which, in reality, suffer from ascertainment bias due to the unequal distribution of case and control populations in nature. In response, here we present LT-sim: a software package that uses a liability threshold linear mixed model to properly simulate ascertainment and other properties that naturally occur in case-control data. We apply LT-Sim in common and rare disease scenarios of type 1 diabetes and ALS, and demonstrate the software’s ability to enable *in-silico* hypothesis testing efficiently.

<img src="https://github.com/lcrawlab/ltsim/blob/main/figs/fig1_ltim.png" height="350">

The objective of LT-Sim is to generate synthetic genetic data that mimic biological relationships among individuals while maintaining population structure and accounting for selection biases frequently present in case-control studies. A. Such studies are commonly modeled using a liability threshold model, which provides a continuous liability distribution to generate case/control phenotypes.
On the latent liability model, the x-axis represents the continuous liability scale, while the y-axis represents the density of the population corresponding to a liability score. The liability threshold (red dotted line), separates individuals who are affected with the disease from those who are not. B. LT-Sim provides 13 parameters at the population (parameters 1 - 7), observation (parameters 8,9), and genotype levels (parameters 10 - 13) so that users can simulate the full spectrum of case-control studies.

### Installation
1) Clone the LT-Sim repo: 
    - `git clone https://github.com/lcrawlab/ltsim.git`
2) Install dependencies:
    - `R` (4.1.1)
    - `parallel` (4.1.1)
    - `truncnorm` (1.0.8)
    - `hash` (2.2.6.1)
    - `dplyr` (1.0.7)

### Running LT-Sim 

#### Overview
LT-Sim can be run in 2 ways, 

1. with single inputs for each input parameter using `lt-sim.r` 

2. with multiple inputs for one/multiple input parameters using `lt-sim.sh`.

#### To run

1) Given a parameter setting of interest, the user can immediately submit those values as input parameters to the LT-Sim R script. For example:

`Rscript lt-sim.r /path/for/outputs --ind 2000 --k 0.1 --obs 200 --tot_snp_sim 1000 --frac_causal 0.1`

Additional argument information and defaults for all other arguments can be seen with `Rscript lt-sim.r --help`:

```
usage: lt-sim.r [--] [--help] [--hierarchy] [--opts OPTS] [--seed SEED]
       [--mask_dir MASK_DIR] [--degree DEGREE] [--ind IND] [--k K]
       [--pve PVE] [--fst FST] [--obs OBS] [--prop_case PROP_CASE]
       [--tot_snp_sim TOT_SNP_SIM] [--frac_causal FRAC_CAUSAL]
       [--maf_lower MAF_LOWER] [--rho RHO] [--prop_pop_A PROP_POP_A]
       [--num_reps NUM_REPS] output_dir

LTSim

positional arguments:
  output_dir         Path for generated output

flags:
  -h, --help         show this help message and exit
  --hierarchy        Set this flag to sample causal SNPs according to
                     biological annotations (e.g., genes or pathways)

optional arguments:
  -x, --opts         RDS file containing argument values
  -s, --seed         Random seed to use in generating simulated data
                     [default: 0]
  -m, --mask_dir     Path to SNP/pathway hierarchy mask (default: NA)
  -d, --degree       Degree of non-overlapping pathways if mask used
                     [default: 0]
  -i, --ind          Number of individuals in population [default:
                     1000]
  -k, --k            Prevalence of disease/trait in population
                     [default: 0.15]
  -p, --pve          Broad sense heritability [default: 0.4]
  -f, --fst          Targeted fixation index as a measure of population
                     differentiation due to genetic structure [default:
                     0.05]
  -o, --obs          Total observed population [default: 200]
  --prop_case        Proportion of case observations [default: 0.5]
  -t, --tot_snp_sim  Total number of SNPs in population [default: 1000]
  --frac_causal      Proportion of causal SNPs in population [default:
                     0.1]
  --maf_lower        Minor allele frequency lower bound [default: 0.05]
  -r, --rho          Proportion of additive genetic variation [default:
                     0.5]
  --prop_pop_A       Proportion of population A relative to total
                     [default: 0.5]
  -n, --num_reps     Number of dataset replicates [default: 5]
```

2) Given a parameter sweep of interest, generate the file(s) for which the sweep is of interest, where each value of interest for a given input parameter is separated by a space, and pass those directly to this Bash script:

`./lt-sim.sh $seed $mask_dir $output_dir $non_overlap $ind $k $pve $fst $obs $prop_case_tot_snp_sim $frac_causal $hierarchy`
