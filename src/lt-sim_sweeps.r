# lt-sim_sweep.r

suppressPackageStartupMessages(library("argparser"))

# Clear console
cat("\014")

# Clear environment
rm(list = ls(all = TRUE))

# create parser object
parser <- arg_parser("LT-Sim Parameter Sweep Helper Script")

# specify our desired options 
parser <- add_argument(parser, "output_dir", type="character",
                       help="Path for generated output")
parser <- add_argument(parser, "--ind_file", type="character", default=NA,
                       help="Path to file with number of individuals in population to test")
parser <- add_argument(parser, "--k_file", type="character", default=NA,
                       help="Path to file with prevalences of disease/trait in population to test")
parser <- add_argument(parser, "--pve_file", type="character", default=NA,
                       help="Path to file with broad sense heritability values to test")
parser <- add_argument(parser, "--fst_file", type="character", default=NA,
                       help="Path to file with FST values to test")
parser <- add_argument(parser, "--obs_file", type="character", default=NA,
                       help="Path to file with total observed population values to test")
parser <- add_argument(parser, "--prop_case_file", type="character", default=NA,
                       help="Path to file with proportion of case observations to test")
parser <- add_argument(parser, "--tot_snp_sim_file", type="character", default=NA,
                       help="Path to file with total number of SNPs in population to test")
parser <- add_argument(parser, "--frac_causal_file", type="character", default=NA,
                       help="Path to file with proportion of causal SNPs in population to test")

# Other arguments to sweep over
parser <- add_argument(parser, "--maf_lower_file", type="character", default=NA,
                       help="Path to file with minor allele frequency lower bounds to test")
parser <- add_argument(parser, "--rho_file", type="character", default=NA,
                       help="Path to file with proportions of additive genetic variation to test")
parser <- add_argument(parser, "--prop_pop_A_file", type="character", default=NA,
                       help="Path to file with proportions of population A relative to total to test")

parser <- add_argument(parser, "--seed", type="integer", default=0,
                       help="Random seed to use in generating simulated data")
parser <- add_argument(parser, "--mask_dir", type="character", default=NA,
                       help="Path to SNP/pathway hierarchy mask (default: NA)")
parser <- add_argument(parser, "--degree", type="integer", default=0,
                       help="Degree of non-overlapping pathways if mask used")
parser <- add_argument(parser, "--hierarchy", flag=TRUE,
                       help="Set this flag to sample causal SNPs according to biological annotations (e.g., genes or pathways)")
parser <- add_argument(parser, "--num_reps", type="integer", default=5,
                       help="Number of dataset replicates")

parser <- add_argument(parser, "--dry_run", flag=TRUE,
                       help="Print commands needed to run parameter sweep to stdout without running lt-sim")


#### Argument assignments ####
argv <- parse_args(parser)

# Get files to read
file_names = names(argv)[grep('_file', names(argv))]
file_vars = sapply(file_names, function(x) {
  gsub('_file', '', x)
}, USE.NAMES = F)

# Read and store all non-default values from files and corresponding argument names
file_vals = list()
for (i in seq_along(file_names)) {
  var = file_vars[i]
  file = file_names[i]
  if (!is.na(argv[[file]])) {
    file_vals[[var]] = strsplit(readLines(argv[[file]]), " ")[[1]]
  }
}

# User error if no files were passed 
if (length(file_vals) == 0) {
  stop("No parameters to sweep; use lt-sim.r to simulate for single values instead.")
}

file_val_grid = expand.grid(file_vals)
file_val_grid[] <- lapply(file_val_grid, as.character)

# Call lt-sim.r
for (i in 1:nrow(file_val_grid)) {
  command = paste0('Rscript lt-sim.r ', argv$output_dir, ' --seed ', argv$seed,
                   ' --degree ', argv$degree, ' --num_reps ', argv$num_reps)
  # Some arguments handled specially
  if (!is.na(argv$mask_dir)) {
    command = paste0(command, ' --mask_dir ', argv$mask_dir)
  }
  if (isTRUE(argv$hierarchy)) {
    command = paste0(command, ' --hierarchy')
  }
  # add arguments to sweep over
  sweep_command = paste0(paste0(' --', names(file_vals), ' ', 
                                file_val_grid[i,]), 
                         collapse='')
  command = paste0(command, sweep_command)
  
  # System call
  if (argv$dry_run) {
    cat(command, '\n')
  } else {
    system(command)
  }
  
}
