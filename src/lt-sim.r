# lt-sim.r

suppressPackageStartupMessages(library("argparser"))

# Clear console
cat("\014")

# Clear environment
rm(list = ls(all = TRUE))

# create parser object
parser <- arg_parser("LT-Sim")

# specify our desired options 
parser <- add_argument(parser, "output_dir", type="character",
                       help="Path for generated output")
parser <- add_argument(parser, "--seed", type="integer", default=0,
                    help="Random seed to use in generating simulated data")
parser <- add_argument(parser, "--mask_dir", type="character", default=NA,
                    help="Path to SNP/pathway hierarchy mask (default: NA)")
parser <- add_argument(parser, "--degree", type="integer", default=0,
                       help="Degree of non-overlapping pathways if mask used")
parser <- add_argument(parser, "--ind", type="integer", default=1000,
                       help="Number of individuals in population")
parser <- add_argument(parser, "--k", type="double", default=0.15,
                       help="Prevalence of disease/trait in population")
parser <- add_argument(parser, "--pve", type="double", default=0.4,
                       help="Broad sense heritability")
parser <- add_argument(parser, "--fst", type="double", default=0.05,
                       help="Targeted fixation index as a measure of population \ differentiation due to genetic structure. \ Set to 0 to simulate single population.")
parser <- add_argument(parser, "--obs", type="integer", default=200,
                       help="Total observed population")
parser <- add_argument(parser, "--prop_case", type="double", default=0.5,
                       help="Proportion of case observations")
parser <- add_argument(parser, "--tot_snp_sim", type="integer", default=1000,
                       help="Total number of SNPs in population")
parser <- add_argument(parser, "--frac_causal", type="double", default=0.1,
                       help="Proportion of causal SNPs in population")
parser <- add_argument(parser, "--hierarchy", flag=TRUE,
                       help="Set this flag to sample causal SNPs according to biological annotations (e.g., genes or pathways)")

# Other optional arguments 
parser <- add_argument(parser, "--maf_lower", type="double", default=0.05,
                       help="Minor allele frequency lower bound")
parser <- add_argument(parser, "--rho", type="double", default=0.5,
                       help="Proportion of additive genetic variation")
parser <- add_argument(parser, "--prop_pop_A", type="double", default=0.5,
                       help="Proportion of population A relative to total")
parser <- add_argument(parser, "--num_reps", type="integer", default=5,
                       help="Number of dataset replicates")


# Load R libraries
library(truncnorm)
library(dplyr)
library(hash)
library(parallel)

#### Argument assignments ####
# Script settings
argv <- parse_args(parser)

SEED <- as.integer(argv$seed); set.seed(SEED)
MASK_DIR <- argv$mask_dir # NA or directory
SIM_OUTDIR <- argv$output_dir
NON_OVERLAP_USED <- as.integer(argv$degree)

# Population number params
ind <- as.integer(argv$ind) # Number of individuals in the population
k <- as.double(argv$k) # Prevalence of the disease/trait 
pve <- as.double(argv$pve) # Broad sense heritability             
fst <- as.double(argv$fst)# 0 if none, value above 0 otherwise (0.005)

# Sampling params
obs <- as.integer(argv$obs) # Number of cases and controls
prop_case <- as.double(argv$prop_case)# Fraction case samples

# Causal SNP params
tot_snp_sim <- as.integer(argv$tot_snp_sim) # Number of total SNPs in the data (e.g. 1e5)
frac_causal <- as.double(argv$frac_causal) # Fraction of SNPs that are causal (e.g. 0.005)
hierarchy <- as.integer(argv$hierarchy) # 0 if no, 1 if yes
#chunk.size <- as.integer(args[14])# use chunks to merge matrices e.g. 100000 (0 otherwise)

#### Other arguments ####
# Population stat params
propA = argv$prop_pop_A # proportion of population A (pick value 0 - 1)

# Simulate phenotypes params
maf_lower = argv$maf_lower # minor allele frequency lower bound
rho = argv$rho # proportion of additive genetic variation

# Set dataset replicate number to simulate 
ndatasets = argv$num_reps # number dataset replicates
#########################

cat("\nWelcome to LTSim. Lets get started.\n")
cat("\nParameters set: 
1) output_dir    ",SIM_OUTDIR,"
2) seed          ",SEED," 
3) mask_dir      ",MASK_DIR,"
4) non_overlap   ",NON_OVERLAP_USED,"
5) ind           ",ind,"
6) k             ",k,"
7) pve           ",pve,"
8) fst           ",fst,"
9) obs           ",obs,"
10) prop_case    ",prop_case,"
11) tot_snp_sim  ",tot_snp_sim,"
12) frac_causal  ",frac_causal,"
13) hierarchy    ",hierarchy,"
14) maf_lower     ",maf_lower," (MAF lower bound threshold)
15) rho          ",rho,"  (additive heritability proportion)
16) prop_pop_A   ",propA,"  (population A proportion in total population (if fst>0))
17) num_datasets ",ndatasets,"    (number dataset replicates)
")


# Check that sampling is possible
if(!(ind*k > obs*prop_case)) {
    cat("\n")
    stop("Please resolve calculation. The affected/diseased population (ind*k=",ind*k,") must be larger than the desired number of diseased individuals to sample (obs*prop_case=",obs*prop_case,").\n")
}

if(hierarchy){
    if(NON_OVERLAP_USED){
        cat("Non-overlap masks used.\n")
        SUB_OUTDIR_NAME = paste("non_overlap_degree",NON_OVERLAP_USED,"ind",ind,"tot_snp_sim",tot_snp_sim,"frac_causal",frac_causal,"k",k,"obs",obs,"maf_lower",maf_lower,"pve",pve,"rho", rho, "prop_case", prop_case,"fst", fst, sep="-")
        # Input biological masks
        MASK1_DIR <- paste(MASK_DIR,"masksim_1_non_overlap_labeled.txt",sep="/")
        MASK2_DIR <- paste(MASK_DIR,"masksim_2_non_overlap_labeled.txt",sep="/")
        MASK3_DIR <- paste(MASK_DIR,"masksim_3_non_overlap_labeled.txt",sep="/")

    }else{
        cat("Overlap masks used.\n")
        SUB_OUTDIR_NAME = paste("ind",ind,"tot_snp_sim",tot_snp_sim,"frac_causal",frac_causal,"k",k,"obs",obs,"maf_lower",maf_lower,"pve",pve,"rho", rho, "prop_case", prop_case, "fst", fst, sep="-")

        # Input biological masks
        MASK1_DIR <- paste(MASK_DIR,"mask1_pd.txt",sep="/")
        MASK2_DIR <- paste(MASK_DIR,"mask2_pd.txt",sep="/")
        MASK3_DIR <- paste(MASK_DIR,"mask3_pd.txt",sep="/")
    }
}else{
    cat("No hierarchy information used.\n")
    SUB_OUTDIR_NAME = paste("no_hier-ind",ind,"tot_snp_sim",tot_snp_sim,"frac_causal",frac_causal,"k",k,"obs",obs,"maf_lower",maf_lower,"pve",pve,"rho", rho, "prop_case", prop_case, "fst", fst, sep="-")
}
######################################################

# Create output dirs if needed
path = paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/")
dir.create(path, showWarnings = TRUE, recursive = TRUE, mode = "0777")

# Get fraction causal pathways
get_causal_pathways<- function(mask2_pd, tot_pathways, frac_causal){    
    causal_pathways = sample(colnames(mask2_pd), size=round(tot_pathways*frac_causal), replace=F)
    return(list(causal_pathways=causal_pathways, n_c_pathways=length(causal_pathways))) # causal pathways, number causal pathways
}

# Get fraction causal genes from each of those causal pathways
get_causal_genes<- function(mask2_pd, causal_pathways, frac_causal){
    tmp = mask2_pd[,causal_pathways]
    causal_genes = c()
    for(curr_cp in causal_pathways){ # for each current causal pathway (curr_cp)
        curr_cp_genes = rownames(tmp[curr_cp]%>%dplyr::filter_all(all_vars(.==1)))
        ncurr_cp_genes = length(curr_cp_genes)
        sample_size = round(ncurr_cp_genes*frac_causal)
        if(ncurr_cp_genes==1){
            tmp_causal_genes = curr_cp_genes
        }else if(sample_size == 0){
            tmp_causal_genes = sample(curr_cp_genes, size=1, replace=F)
        }else{
            tmp_causal_genes = sample(curr_cp_genes, size=sample_size, replace=F)
        }
        causal_genes = c(causal_genes, tmp_causal_genes)
    }
    
    # Keep set of genes (could have been sampled twice from same gene beloning to >1 pathway)
    return(list(causal_genes=unique(causal_genes), n_c_genes=length(causal_genes)))
}

# Get fraction causal eQTLs/SNPs from each of those causal genes
get_causal_snps <- function(mask1_pd, causal_genes, frac_causal){
    tmp = mask1_pd[,causal_genes]
    causal_snps = c()
    for(curr_cg in causal_genes){ # for each current causal gene (curr_cg)
        #print(curr_cg)
        curr_cg_snps = rownames(tmp[curr_cg]%>%dplyr::filter_all(all_vars(.==1)))
        #print(curr_cg_snps)
        ncurr_cg_snps = length(curr_cg_snps)
        sample_size = round(ncurr_cg_snps*frac_causal)
        if(ncurr_cg_snps==1){
            tmp_causal_snps = curr_cg_snps
        }else if(sample_size == 0){
            tmp_causal_snps = sample(curr_cg_snps, size=1, replace=F)
        }else{
            tmp_causal_snps = sample(curr_cg_snps, size=sample_size, replace=F)
        }
        causal_snps = c(causal_snps, tmp_causal_snps)
    }
    
    # Keep set of SNPs (could have been sampled twice from same eQTL SNP mapping to >1 genes)
    return(list(causal_snps=unique(causal_snps), n_c_snps=length(causal_snps)))
}


## MAIN
cores <- detectCores()-2

# Name all SNPs by 's<INT>' for total number of SNPs
all_snp_names = paste0("s",1:tot_snp_sim)

# Uniformly sample causal SNPs/eQTLs (via biological hierarchy or without)
if (hierarchy){
    # Import masks data
    mask1_pd = read.csv(MASK1_DIR, sep=" ", row.names=1)
    mask2_pd = read.csv(MASK2_DIR, sep=" ", row.names=1)
    mask3_pd = read.csv(MASK3_DIR, sep=" ", row.names=1)

    cat("Matrix Stats\n")
    tot_pathways = ncol(mask2_pd) #Number of pathways
    cat("total pathways in mask 2\n")
    cat(tot_pathways,"\n")
    tot_genes = nrow(mask2_pd)
    cat("total genes in mask 2\n")
    cat(tot_genes,"\n")
    tot_snp_sim = nrow(mask1_pd)
    cat("total SNPs in mask 1\n")
    cat(tot_snp_sim,"\n")

    cat("Mask 1 (SNPs to genes)\n")
    cat(dim(mask1_pd),"\n")
    cat("Sum of mask 1 (SNPs to genes)\n")
    cat(sum(mask1_pd),"\n")

    cat("Mask 2 (genes to pathways)\n")
    cat(dim(mask2_pd),"\n")
    cat("Sum of mask 2 (genes to pathways)\n")
    cat(sum(mask2_pd),"\n")

    cat("Mask 3 (SNPs to pathways)\n")
    cat(dim(mask3_pd),"\n")
    cat("Sum of mask 3 (SNPs to pathways)\n")
    cat(sum(mask3_pd),"\n")
    
    # Sample fraction causal sampled pathways' fraction causal sampled genes' fraction causal sampled snps
    print("Getting causal pathways...\n")
    causal_p = get_causal_pathways(mask2_pd, tot_pathways, frac_causal)
    print("Getting causal genes in causal pathways...\n")
    causal_g = get_causal_genes(mask2_pd, causal_p$causal_pathways, frac_causal)
    print("Getting causal eQTL/SNPs in causal genes...\n")
    causal_s = get_causal_snps(mask1_pd, causal_g$causal_genes, frac_causal)

    # Ensure no duplicate pathways
    stopifnot("Resolve: duplicated pathways." = length(unique(causal_p))==length(causal_p))

    # Print number of causal features
    cat("Number causal pathways:\n")
    cat(causal_p$n_c_pathways,"\n")
    cat("Number causal genes:\n")
    cat(causal_g$n_c_genes,"\n")
    cat("Number causal eQTLs/SNPs:\n")
    cat(causal_s$n_c_snps,"\n")

}else{
    # Uniformly sample SNPs/eQTLs
    cat("\nSample causal eQTL/SNPs uniformly at random.\n")
    causal_snps = unique(sample(all_snp_names, size=round(tot_snp_sim*frac_causal), replace=F))
    causal_s = list(causal_snps=causal_snps, n_c_snps=length(causal_snps))
}

# Generate causal SNP genome across individuals (ind)
cat("\nCreate genome:\n")
if (fst>0){
    # Create population structure
    maf <- maf_lower + 0.4*runif(tot_snp_sim)
    saveRDS(maf, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"maf.rds",sep="/"))
    cat("Creating genome with FST...\n")
    delta = sqrt(maf-maf^2-maf*(1-maf)*(1-fst))
    saveRDS(delta, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"delta.rds",sep="/"))       
    # Population A
    cat("\nCreating subsets pop. A ...\n")
    subsets <- splitIndices((ind*propA)*tot_snp_sim, cores) # Split indices
    #saveRDS(subsets, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"subsets_A.rds",sep="/"))
    set.seed(SEED)
    rand_numbers <- runif((ind*propA)*tot_snp_sim)
    rand_numbers2 <- runif((ind*propA)*tot_snp_sim)
    cat("\nCreating genome pop. A ...\n")
    geno_list <- mclapply(subsets, function(x) {
                 (rand_numbers[x] < (maf+delta)) + (rand_numbers2[x] < (maf+delta))
     }, mc.cores = cores, mc.set.seed = TRUE)
    rm(subsets);rm(rand_numbers);rm(rand_numbers2)
    saveRDS(geno_list, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno_list_A.rds",sep="/"))
    GenoA <- unlist(geno_list); rm(geno_list) 
    GenoA <- matrix(as.double(GenoA),(ind*propA),tot_snp_sim,byrow = TRUE)
    rownames(GenoA) <- paste0("popA",1:dim(GenoA)[1])
    colnames(GenoA) <- all_snp_names
    saveRDS(GenoA, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"GenoA.rds",sep="/"))
    cat("Genome population A dimension: ",dim(GenoA),"\n")
    
    # Population B
    cat("\nCreating subsets pop. B ...\n")
    subsets <- splitIndices((ind*(1-propA))*tot_snp_sim, cores) # Split indices
    saveRDS(subsets, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"subsets_B.rds",sep="/"))
    set.seed(SEED)
    rand_numbers <- runif((ind*(1-propA))*tot_snp_sim)
    rand_numbers2 <- runif((ind*(1-propA))*tot_snp_sim)
    cat("\nCreating genome pop. B ...\n")
    geno_list <- mclapply(subsets, function(x) {
                 (rand_numbers[x] < (maf-delta)) + (rand_numbers2[x] < (maf-delta))
     }, mc.cores = cores, mc.set.seed = TRUE)
    rm(subsets);rm(rand_numbers);rm(rand_numbers2)
    saveRDS(geno_list, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno_list_B.rds",sep="/"))
    GenoB <- unlist(geno_list); rm(geno_list)
    GenoB <- matrix(as.double(GenoB),(ind*(1-propA)),tot_snp_sim,byrow = TRUE)
    len_geno1 = dim(GenoA)[1]+1
    rownames(GenoB) <- paste0("popB",len_geno1:ind)
    colnames(GenoB) <- all_snp_names
    saveRDS(GenoB, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"GenoB.rds",sep="/"))
    cat("Genome population B dimension: ",dim(GenoB),"\n")
    
    cat("Combining populations A and B...\n")
    geno = rbind(GenoA,GenoB); rm(GenoA); rm(GenoB) # geno = big.matrix(cbind(GenoA, GenoB));rm(GenoA);rm(GenoB)
    colnames(geno) <- all_snp_names
    saveRDS(geno, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno.rds",sep="/")) #save(geno, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno.RData",sep="/"),compress = TRUE)
    
}else if(fst==0){
    maf <- maf_lower + 0.45*runif(tot_snp_sim)
    saveRDS(maf, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"maf.rds",sep="/"))
    # Creating genomes
    subsets <- splitIndices((ind/2)*tot_snp_sim, cores)
    set.seed(SEED)
    rand_numbers <- runif((ind/2)*tot_snp_sim)
    rand_numbers2 <- runif((ind/2)*tot_snp_sim)
    cat("\nCreating genome first half pop. ...\n")
    geno_list <- mclapply(subsets, function(x) {
                     (rand_numbers[x] < maf) + 
                     (rand_numbers2[x] < maf)
    }, mc.cores = cores, mc.set.seed = TRUE)
    rm(subsets);rm(rand_numbers);rm(rand_numbers2);rm(subsets)
    saveRDS(geno_list, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno_list_1st_half.rds",sep="/"))
    geno_first <- unlist(geno_list); rm(geno_list) 
    geno_first <- matrix(as.double(geno_first),(ind/2),tot_snp_sim,byrow = TRUE)
    rownames(geno_first) <- paste0("pop",1:dim(geno_first)[1])
    colnames(geno_first) <- all_snp_names
    saveRDS(geno_first, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno_first.rds",sep="/"))
    cat("Genome population first half dimension: ",dim(geno_first),"\n")
    
    subsets <- splitIndices((ind/2)*tot_snp_sim, cores)
    set.seed(SEED)
    rand_numbers <- runif((ind/2)*tot_snp_sim)
    rand_numbers2 <- runif((ind/2)*tot_snp_sim)
    cat("\nCreating genome first half pop. ...\n")
    geno_list <- mclapply(subsets, function(x) {
                     (rand_numbers[x] < maf) + 
                     (rand_numbers2[x] < maf)
    }, mc.cores = cores, mc.set.seed = TRUE)
    rm(subsets);rm(rand_numbers);rm(rand_numbers2)
    saveRDS(geno_list, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno_list_2nd_half.rds",sep="/"))
    geno_second <- unlist(geno_list); rm(geno_list) 
    geno_second <- matrix(as.double(geno_second),(ind/2),tot_snp_sim,byrow = TRUE)
    len_geno_first = dim(geno_first)[1]+1
    rownames(geno_second) <- paste0("pop",len_geno_first:ind)
    colnames(geno_second) <- all_snp_names
    saveRDS(geno_second, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno_second.rds",sep="/"))
    cat("Genome population second half dimension: ",dim(geno_second),"\n")

    cat("Combining populations (first and second half)...\n")
    geno = rbind(geno_first,geno_second); rm(geno_first); rm(geno_second) # geno = big.matrix(cbind(GenoA, GenoB));rm(GenoA);rm(GenoB)
    colnames(geno) <- all_snp_names
    saveRDS(geno, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno.rds",sep="/")) 
    ######
    
    ### PARALLEL VERSION ###
#     subsets <- splitIndices(ind*tot_snp_sim, cores)
#     set.seed(SEED)
#     rand_numbers <- runif(ind*tot_snp_sim)
#     rand_numbers2 <- runif(ind*tot_snp_sim)
#     geno_list <- mclapply(subsets, function(x) {
#                      (rand_numbers[x] < maf) + 
#                      (rand_numbers2[x] < maf)
#     }, mc.cores = cores, mc.set.seed = TRUE)
#     saveRDS(geno_list, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno_list.rds",sep="/"))
#     geno <- unlist(geno_list); rm(geno_list)
#     geno <- matrix(as.double(geno),ind,tot_snp_sim,byrow=TRUE) 
#     rownames(geno) <- paste0("pop", 1:ind)
#     colnames(geno) <- all_snp_names
#     cat("\ngeno parallel\n")
#     print(geno)
#     stop("smile")
#    saveRDS(geno, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno.rds",sep="/"))
    ######
    
}else{
    exit("Please enter an FST value >0, or type 0 to run LTSim without FST.\n")
}

cat("Population genomic complete and saved here:", paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"geno.rds",sep="/"), "\n")
cat("Genome dimensions: ",dim(geno),"\n")
cat("and it looks like: \n")
print(geno[1:4,1:4])

# Calculate mean and std. for X genotype matrix.
Xmean <- apply(geno, 2, mean); 
saveRDS(Xmean, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"Xmean.rds",sep="/"))
Xsd <- apply(geno, 2, sd);
saveRDS(Xsd, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"Xsd.rds",sep="/"))

# Generate replicates, where each has a new partition of causal SNPs into one of the interaction groups, or into the additive group. 
# NOTE: keeping causal genes and pathways indexed by dataset if we want to increase complexity
for (ndat in 1:ndatasets){
    cat("\nDataset replicate: ", ndat, "\n")
    cat("\nSample SNPs:\n")
    nsnp = round(tot_snp_sim*frac_causal)
    
    cat(paste("Split", tot_snp_sim*frac_causal,"causal SNPs into causal 1, causal 2, and causal 3"),"\n")
    nthird_causal_snps = floor(causal_s$n_c_snps/3)
    
    cat("Sample interacting (s1, s2) and additive (s3) SNPs.\n")
    s1=sample(causal_s$causal_snps, nthird_causal_snps, replace=F)
    ncausal1= length(s1)
    
    s2=sample(causal_s$causal_snps[causal_s$causal_snps%in%s1==FALSE], nthird_causal_snps, replace=F)
    ncausal2= length(s2)
    
    s3=sample(causal_s$causal_snps[causal_s$causal_snps%in%c(s1,s2)==FALSE], causal_s$n_c_snps-length(s1)-length(s2), replace=F)
    ncausal3= length(s3)
    
    cat("Center and scale causal SNP sets.\n")
    Xcausal1=t((t(geno[,s1])-Xmean[s1])/Xsd[s1]);
    Xcausal2=t((t(geno[,s2])-Xmean[s2])/Xsd[s2]);
    Xcausal3=t((t(geno[,s3])-Xmean[s3])/Xsd[s3]);
    
    cat("Number interacting causal SNPS s1:", ncausal1,"\n")
    cat("Number interacting causal SNPS s2:", ncausal2,"\n")
    cat("Number additive causal SNPS s3:", ncausal3,"\n")
    
    cat("\nCalculate epistatsis and additive effects on phenotype:\nCalculate marginal effects.\n")
    Xmarginal=cbind(Xcausal1,Xcausal2,Xcausal3)
    saveRDS(Xmarginal, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"Xmarginal.rds",sep="/"))
    beta=rnorm(dim(Xmarginal)[2])
    y_marginal=c(Xmarginal%*%beta)
    beta=beta*sqrt(pve*rho/var(y_marginal))
    y_marginal=Xmarginal%*%beta
    
    #Remove unneeded marginal parameters for storage
    rm(Xmarginal); rm(Xcausal3)
    
    cat("Create causal epistatic matrix between s1 and s2 SNPs with dimension.\n")
    Xepi=c()
    Xepi <- do.call(cbind,lapply(1:ncol(Xcausal1), function(x) Xcausal1[,x]*Xcausal2))
    saveRDS(Xepi, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"Xepi.rds",sep="/"))
   
    #Remove unneeded causal SNP parameters for storage"
    rm(Xcausal1); rm(Xcausal2)

    cat("Calculate epistatic effects.\n")
    beta=rnorm(dim(Xepi)[2])
    y_epi=c(Xepi%*%beta)
    beta=beta*sqrt(pve*(1-rho)/var(y_epi))
    y_epi=Xepi%*%beta

    #Remove unneeded parameters for storage.
    rm(Xepi)

    cat("Calculate liability threshold.\n")
    y_err=rnorm(ind)
    y_err=y_err*sqrt((1-pve)/var(y_err))

    cat("Continuous phenotype calculation.\n")
    y=y_marginal+y_epi+y_err

    #Remove unneeded parameters for storage.
    rm(y_marginal); rm(y_epi)

    # Set threshold                             
    thresh=qnorm(1-k,mean=0,sd=1)
    cat("Set threshold: ", thresh,"\n")

    # Creatae case control phenotype vector
    cat("Create case control phenotype vector.\n")
    pheno=(y>thresh)+1

    cat("\nGet number of unaffected and affected people:\n")
    pop_controls = sum(pheno==1)
    cat("Proportion of normal indivs: ", pop_controls/length(y),"\n")
    cat("Number unaffected: ", pop_controls, "\n") #todo: update variable names to be normal_indivs                             
    pop_cases = sum(pheno==2)
    cat("Proportion of diseased indivs: ", pop_cases/length(y),"\n")
    cat("Number affected: ", pop_cases, "\n") #todo: update variable names to be disease_indivs
    
    # Number diseased individuals to sample
    num_cases_to_sample = floor(obs*prop_case)
    cat("Number of affected/diseased individuals to sample: ", num_cases_to_sample,"\n")

    # Check that there are enough cases to sample for observations
    actual_cases = which(y>thresh)
    num_actual_cases = length(actual_cases)
    if(!(num_actual_cases)>floor(num_cases_to_sample)){
        cat("\nWARNING: There are barely enough affected/diseased observations from the population to sample under the threshold.\n")
        cat(num_actual_cases, "!>", num_cases_to_sample,"\n")
        
        # Recalculate possible observations to maintain k and proportion of case to control
        cases = sample(actual_cases, num_actual_cases, replace = FALSE)
        cat("cases", length(cases),"\n")
        obs_new = floor(num_actual_cases/num_cases_to_sample*obs)
        num_controls_new = obs_new - num_actual_cases
        controls = sample(which(y<=thresh),num_controls_new,replace = FALSE)
        cat("controls", length(controls),"\n")
    }else{
        cat("Sample cases and controls given calculated threshold.\n")
        cases = sample(which(y>thresh), num_cases_to_sample, replace = FALSE)
        cat("cases", length(cases),"\n")
        controls = sample(which(y<=thresh),floor(obs-num_cases_to_sample),replace = FALSE)
        cat("controls", length(controls),"\n")
    }
    
    y = pheno[c(cases,controls)] # vector of 2s (affected) and 1s (unaffected)
    X = geno[c(cases,controls),] 
    y_matrix = cbind(rownames(X), y)
    
    cat("\nFinalize results and save:\n")
    cat("Dataset",ndat,"dimensions for X:",dim(X),"matrix and y:",length(y),"vector.\n")           
                                 
    if (hierarchy){
        cat("Ensure no duplicate eQTL/SNP names.\n")
        stopifnot("Resolve: duplicate eQTL/SNP names."=length(unique(colnames(X)))==length(colnames(X)))
    }
                                
    cat("Center and scale final X (causal and non-causal SNPs/eQTLs).\n")
    Xmean=apply(X, 2, mean); 
    Xsd=apply(X, 2, sd); 
    Xcs=t((t(X)-Xmean)/Xsd)
    saveRDS(Xcs, file=paste(SIM_OUTDIR,SUB_OUTDIR_NAME,"Xcs.rds",sep="/"))
                                 
    if (hierarchy){
        # Filter mask 1 to include simulation eQTLs/SNPs mapping to genes")
        tmp = mask1_pd[colnames(X),]
        mask1_sim = tmp[,(which(colSums(tmp)>0))]

        print("mask1_sim")
        print(dim(mask1_sim))

        cat("\nCheck that X column names match mask1_sim")
        stopifnot("Resolve: X SNP order does not match mask1_sim."=identical(colnames(X), rownames(mask1_sim))==T)

        cat("\nFilter mask 2 to include simulation genes mapping to pathways")
        tmp = mask2_pd[colnames(mask1_sim),]
        mask2_sim = tmp[,(which(colSums(tmp)>0))]
        print("mask2_sim")
        print(dim(mask2_sim))

        cat("\nFilter mask 3 to include simulation genes mapping to pathways")
        tmp = mask3_pd[rownames(mask1_sim),]
        mask3_sim = tmp[,(which(colSums(tmp)>0))]
        print("mask3_sim")
        print(dim(mask3_sim))
    
        cat("Save masks, causal pathways, and causal genes.")
        #Save masks labeled row and columns
        write.table(file=paste(path,paste("masksim_1_labeled", as.character(ndat),"txt", sep="."),sep="/"), mask1_sim, sep=" ", quote=FALSE)
        write.table(file=paste(path,paste("masksim_2_labeled", as.character(ndat),"txt", sep="."),sep="/"), mask2_sim, sep=" ", quote=FALSE)
        write.table(file=paste(path,paste("masksim_3_labeled", as.character(ndat),"txt", sep="."),sep="/"), mask3_sim, sep=" ", quote=FALSE)

        # Save causal pathways, genes, and eQTLs/SNPs
        write.table(file=paste(path,paste("causal_pathways", as.character(ndat),"txt", sep="."),sep="/"), causal_p$causal_pathways, sep=" ", quote=FALSE,row.names=FALSE, col.names=FALSE)
        write.table(file=paste(path,paste("causal_genes", as.character(ndat),"txt", sep="."),sep="/"), causal_g$causal_genes, sep=" ", quote=FALSE,row.names=FALSE, col.names=FALSE)
    }
    
    cat("Save X, y, and causal SNPs.\n")
    write.table(file=paste(path,paste("Xsim_labeled", as.character(ndat),"txt", sep="."),sep="/"), X, sep=" ", quote=FALSE)
    write.table(file=paste(path,paste("Xsim_labeled_cent_scal", as.character(ndat),"txt", sep="."),sep="/"), Xcs, sep=" ", quote=FALSE)
    write.table(file=paste(path,paste("ysim_labeled", as.character(ndat),"txt", sep="."),sep="/"), y_matrix, sep=" ", quote=FALSE, row.names=FALSE)
    
    # Save causal SNPs/eQTLs
    write.table(file=paste(path,paste("causal_snps_epistatic", as.character(ndat),"txt", sep="."),sep="/"), c(s1,s2), sep=" ", quote=FALSE,row.names=FALSE, col.names=FALSE)
    write.table(file=paste(path,paste("causal_snps_additive", as.character(ndat),"txt", sep="."),sep="/"), s3, sep=" ", quote=FALSE,row.names=FALSE, col.names=FALSE)
    
    if (hierarchy){
        cat("\nMask 1 annotations:\n")
        #print(head(mask1_sim))
        cat(dim(mask1_sim))

        cat("\nMask 3 annotations:\n")
        #print(head(mask2_sim))
        cat(dim(mask3_sim))
        
        cat(paste("Number of SNP-sets:", dim(mask1_sim)[2], "\n"))
        cat(paste("Number of pathways:", dim(mask3_sim)[2], "\n"))
        cat(paste("Number of causal SNP-sets:", causal_g$n_c_genes, "\n"))
        cat(paste("Number of causal pathways:", causal_p$n_c_pathways, "\n"))
    }
}

cat("\nLT-Sim complete. Results are here:",path,"\n")