# lt_simulations_fst.r
# Last mod. 11/09/2022
# Ashley Mae Conard

# TO DO : MAKE PATHWAY->GENE->SNP HIERARCHY OPTIONAL

# Clear console
cat("\014")

# Clear environment
rm(list = ls(all = TRUE))

# Input
args = commandArgs(TRUE)
if (length(args)==0) {
	  stop("Input: SEED (int), MASK_DIR, SIM_OUTDIR, RES_OUTDIR, NON_OVERLAP_USED (0 or positive integer for degree), k, obs, pve, prop_case fst hierarchy, chunk.size", call.=FALSE)
} else if (length(args) < 12) {
	stop("Input: SEED (int), MASK_DIR, SIM_OUTDIR, RES_OUTDIR, NON_OVERLAP_USED (0 or positive integer for degree), k, obs, pve, prop_case fst hierarchy chunk.size", call.=FALSE)
}


# Load R libraries
library(truncnorm)
library(dplyr)
library(hash)
library(bigmemory)
library(parallel)

# Argument assignments
SEED <- as.integer(args[1]); set.seed(SEED)
MASK_DIR <- args[2] # NA or directory
SIM_OUTDIR <- args[3]
RES_OUTDIR <- args[4]
NON_OVERLAP_USED <-as.integer(args[5])
k <- as.double(args[6]) #Prevalence of the Disease/Trait 
obs <- as.integer(args[7]) #Number of Cases and Controls
pve <- as.double(args[8]) #Broad sense heritability             
prop_case <- as.double(args[9])# fraction case experiments (e.g. in our real ALS data: .883 case experiments)
fst <- as.double(args[10])# 0 if none, value above 0 otherwise (0.005)
hierarchy <- as.integer(args[11]) # 0 if no, 1 if yes
chunk.size <- as.integer(args[12])# use chunks to merge matrices e.g. 100000 (0 otherwise)

################## PARAMS ###########################
# Population Stats
ind = 10 #15e3 #1e6 # Number of Individuals in the Populations
propA = 0.5 # Proportion of Population A (pick value 0 - 1)
tot_snp_sim = 100#1e5 #5e3 #Number of Total SNPs in the Data
frac_causal = 0.5#0.005 # Fraction of SNPs that are causal

##### GOOD FOR TESTING Population params #######  
# ind = 1e3 # Number of Individuals in the Populations
# tot_snp_sim = 80 # Number of Total SNPs in the Data
# frac_causal = .1 # Fraction of SNPs that are causal
# k = 0.1 #Prevalence of the Disease/Trait 
# obs = 50 #Number of Cases and Controls
################################################

# Simulating phenotypes params
maf_frac = 0.05 # MAF fraction
rho=0.5 # proportion that is additive

# Set dataset number to simulate 
ndatasets = 5 #### 10

cat("params:ind",ind,"tot_snp_sim",tot_snp_sim,"frac_causal",frac_causal,"k",k,"obs",obs,"maf_frac",maf_frac,"pve",pve,"rho", rho, "prop_case", prop_case, "fst", fst, sep="-")

if(hierarchy){
    if(NON_OVERLAP_USED){
        print("non-overlap masks used")
        SUB_OUTDIR_NAME = paste("non_overlap_degree",NON_OVERLAP_USED,"ind",ind,"tot_snp_sim",tot_snp_sim,"frac_causal",frac_causal,"k",k,"obs",obs,"maf_frac",maf_frac,"pve",pve,"rho", rho, "prop_case", prop_case,"fst", fst, sep="-")
        # Input biological masks
        MASK1_DIR <- paste(MASK_DIR,"masksim_1_non_overlap_labeled.txt",sep="/")
        MASK2_DIR <- paste(MASK_DIR,"masksim_2_non_overlap_labeled.txt",sep="/")
        MASK3_DIR <- paste(MASK_DIR,"masksim_3_non_overlap_labeled.txt",sep="/")

    }else{
        print("overlap masks used")
        SUB_OUTDIR_NAME = paste("ind",ind,"tot_snp_sim",tot_snp_sim,"frac_causal",frac_causal,"k",k,"obs",obs,"maf_frac",maf_frac,"pve",pve,"rho", rho, "prop_case", prop_case, "fst", fst, sep="-")

        # Input biological masks
        MASK1_DIR <- paste(MASK_DIR,"mask1_pd.txt",sep="/")
        MASK2_DIR <- paste(MASK_DIR,"mask2_pd.txt",sep="/")
        MASK3_DIR <- paste(MASK_DIR,"mask3_pd.txt",sep="/")
    }
}else{
    print("No hierarchy information used")
    SUB_OUTDIR_NAME = paste("no_hier-ind",ind,"tot_snp_sim",tot_snp_sim,"frac_causal",frac_causal,"k",k,"obs",obs,"maf_frac",maf_frac,"pve",pve,"rho", rho, "prop_case", prop_case, "fst", fst, sep="-")
}
######################################################

# Create output dirs if needed
path = paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/")
dir.create(path, showWarnings = TRUE, recursive = TRUE, mode = "0777")
path = paste(RES_OUTDIR,SUB_OUTDIR_NAME,sep="/")
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

    cat("\nMatrix Stats")
    tot_pathways = ncol(mask2_pd) #Number of pathways
    cat("\ntotal pathways in mask 2\n")
    cat(tot_pathways)
    tot_genes = nrow(mask2_pd)
    cat("\ntotal genes in mask 2\n")
    cat(tot_genes)
    tot_snp_sim = nrow(mask1_pd)
    cat("\ntotal SNPs in mask 1\n")
    cat(tot_snp_sim)

    cat("\nMask 1 (SNPs to genes)\n")
    cat(dim(mask1_pd))
    cat("\nSum of mask 1 (SNPs to genes)\n")
    cat(sum(mask1_pd))

    cat("\nMask 2 (genes to pathways)\n")
    cat(dim(mask2_pd))
    cat("\nSum of mask 2 (genes to pathways)\n")
    cat(sum(mask2_pd))

    cat("\nMask 3 (SNPs to pathways)\n")
    cat(dim(mask3_pd))
    cat("\nSum of mask 3 (SNPs to pathways)\n")
    cat(sum(mask3_pd))
    
    # Sample fraction causal sampled pathways' fraction causal sampled genes' fraction causal sampled snps
    print("\nGetting causal pathways...\n")
    causal_p = get_causal_pathways(mask2_pd, tot_pathways, frac_causal)
    print("\nGetting causal genes in causal pathways...\n")
    causal_g = get_causal_genes(mask2_pd, causal_p$causal_pathways, frac_causal)
    print("\nGetting causal eQTL/SNPs in causal genes...\n")
    causal_s = get_causal_snps(mask1_pd, causal_g$causal_genes, frac_causal)

    # Ensure no duplicate pathways
    stopifnot("Resolve: duplicated pathways." = length(unique(causal_p))==length(causal_p))

    # Print number of causal features
    cat("\nNumber causal pathways:\n")
    cat(causal_p$n_c_pathways)
    cat("\nNumber causal genes:\n")
    cat(causal_g$n_c_genes)
    cat("\nNumber causal eQTLs/SNPs:\n")
    cat(causal_s$n_c_snps)

}else{
    # Uniformly sample SNPs/eQTLs
    print("\nSampling causal eQTL/SNPs uniformly at random...\n")
    causal_snps = unique(sample(all_snp_names, size=round(tot_snp_sim*frac_causal), replace=F))
    causal_s = list(causal_snps=causal_snps, n_c_snps=length(causal_snps))
    cat("Example:\n")
    cat(causal_s$causal_snps[1:6])
    cat("Complete.\n")
}

# Generate causal SNP genome across individuals (ind)
cat("Creating genome...")
maf <- maf_frac + 0.4*runif(tot_snp_sim)
if (fst>0){
    # Create population structure
    maf <- maf_frac + 0.4*runif(tot_snp_sim)
    delta = sqrt(maf-maf^2-maf*(1-maf)*(1-fst))
    
    # Split indices
    subsets <- splitIndices((ind*propA)*tot_snp_sim, cores)
    
    # Population A
    set.seed(SEED)
    geno_list <- mclapply(subsets, function(x) {
				    (runif(((ind*propA)*tot_snp_sim)[x]) < (maf+delta)) + (runif(((ind*propA)*tot_snp_sim)[x]) < (maf+delta))
     }, mc.cores = cores)
    Geno1 <- unlist(geno_list); rm(geno_list) 
    Geno1 <- matrix(as.double(Geno1),(ind*propA),tot_snp_sim,byrow = TRUE)
    rownames(Geno1) <- paste0("popA",1:dim(Geno1)[1])
    cat("\nGenome population A dimension:\n")
    cat(dim(Geno1),"\n")

    # Population B
    set.seed(SEED)
    geno_list <- mclapply(subsets, function(x) {
				    (runif(((ind*(1-propA))*tot_snp_sim)[x]) < (maf-delta)) + (runif(((ind*(1-propA))*tot_snp_sim)[x]) < (maf-delta))
     }, mc.cores = cores)
    Geno2 <- unlist(geno_list); rm(geno_list)
    Geno2 <- matrix(as.double(Geno2),(ind*(1-propA)),tot_snp_sim,byrow = TRUE)
    len_geno1 = dim(Geno1)[1]+1
    rownames(Geno2) <- paste0("popB",len_geno1:ind)
    cat("\nGenome population B dimension:\n")
    cat(dim(Geno2))
    
    # Combine populations
    if(chunk.size > 0){

	# Initialize empty matrix
	geno <- matrix(nrow = nrow(Geno1) + nrow(Geno2), ncol = ncol(Geno1))

	# read Geno1 in chunks and add to geno
	for(i in 1:ceiling(nrow(Geno1)/chunk.size)){
            start <- (i-1)*chunk.size + 1
	    end <- min(i*chunk.size, nrow(Geno1))
	    geno[start:end,] <- Geno1[start:end,]
	}

	# read Geno2 in chunks and add to geno
        for(i in 1:ceiling(nrow(Geno2)/chunk.size)){
	    start <- (i-1)*chunk.size + 1 + nrow(Geno1)
	    end <- min(i*chunk.size + nrow(Geno1), nrow(geno))
	    geno[start:end,] <- Geno2[start:end,]
	}
   
       	#remove Geno1 and Geno2
	rm(Geno1)
	rm(Geno2)
    
    }else{
        geno = big.matrix(cbind(Geno1, Geno2));rm(Geno1);rm(Geno2)
        #geno = rbind(Geno1,Geno2); rm(Geno1); rm(Geno2)
    }
}else if(fst==0){
    # Simulate genotypes without population structure
    subsets <- splitIndices(ind*tot_snp_sim, cores)
    set.seed(SEED)
    geno_list <- mclapply(subsets, function(x) {
				    (runif((ind*tot_snp_sim)[x]) < maf) + (runif((ind*tot_snp_sim)[x]) < maf)
     }, mc.cores = cores)
    geno <- unlist(geno_list)
    geno <- matrix(as.double(geno),ind,tot_snp_sim,byrow=TRUE) 
    rownames(geno) <- paste0("pop", 1:ind)
}else{
    exit("Please enter an FST value >0, or type 0 to run LTSim without FST.")
}

cat("\nGenome complete. 
Genome dimensions:\n")
cat(dim(geno),"\n")
colnames(geno) <- all_snp_names
print(geno[1:4,1:4])

# Calculate mean and std. for X genotype matrix.
Xmean <- apply(geno, 2, mean); 
Xsd <- apply(geno, 2, sd);

# Generate replicates, where each has a new partition of causal SNPs into one of the interaction groups, or into the additive group. 
# NOTE: keeping causal genes and pathways indexed by dataset if we want to increase complexity
for (ndat in 1:ndatasets){
    cat("Dataset replicate: ", ndat, "\n")
    cat("\nNumber of causal SNPs")
    nsnp = round(tot_snp_sim*frac_causal)
    
    cat(paste("\nSplit ", tot_snp_sim," into causal 1, causal 2, and causal 3"))
    nthird_causal_snps = floor(causal_s$n_c_snps/3)
    
    cat("\nSample interacting (s1, s2) and additive (s3) SNPs.")
    s1=sample(causal_s$causal_snps, nthird_causal_snps, replace=F)
    ncausal1= length(s1)
    
    s2=sample(causal_s$causal_snps[causal_s$causal_snps%in%s1==FALSE], nthird_causal_snps, replace=F)
    ncausal2= length(s2)
    
    s3=sample(causal_s$causal_snps[causal_s$causal_snps%in%c(s1,s2)==FALSE], causal_s$n_c_snps-length(s1)-length(s2), replace=F)
    ncausal3= length(s3)
    
    cat("\nCenter and scale causal SNP sets.")
    Xcausal1=t((t(geno[,s1])-Xmean[s1])/Xsd[s1]);
    Xcausal2=t((t(geno[,s2])-Xmean[s2])/Xsd[s2]);
    Xcausal3=t((t(geno[,s3])-Xmean[s3])/Xsd[s3]);
    
    cat("\nNumber interacting causal SNPS s1:\n")
    cat(dim(Xcausal1))
    cat("\nNumber interacting causal SNPS s2:\n")
    cat(dim(Xcausal2))
    cat("\nNumber additive causal SNPS s3:\n")
    cat(dim(Xcausal3))
    
    cat("\nCalculate marginal effects")
    Xmarginal=cbind(Xcausal1,Xcausal2,Xcausal3)
    beta=rnorm(dim(Xmarginal)[2])
    y_marginal=c(Xmarginal%*%beta)
    beta=beta*sqrt(pve*rho/var(y_marginal))
    y_marginal=Xmarginal%*%beta
    
    cat("\nRemove unneeded marginal parameters for storage")
    rm(Xmarginal); rm(Xcausal3)
    
    cat("\nCreate causal epistatic matrix between s1 and s2 SNPs with dimension:\n")
    Xepi=c()
    for(i in 1:ncausal1){
      Xepi=cbind(Xepi,Xcausal1[,i]*Xcausal2)
    }
   
    cat("\nRemove unneeded causal SNP parameters for storage")
    rm(Xcausal1); rm(Xcausal2)

    cat("\nCalculate epistatic effects")
    beta=rnorm(dim(Xepi)[2])
    y_epi=c(Xepi%*%beta)
    beta=beta*sqrt(pve*(1-rho)/var(y_epi))
    y_epi=Xepi%*%beta

    cat("\nRemove unneeded parameters for storage")
    rm(Xepi)

    cat("\nCalculate error (i.e. liability threshold)")
    y_err=rnorm(ind)
    y_err=y_err*sqrt((1-pve)/var(y_err))

    cat("\nContinuous phenotype calculation")
    y=y_marginal+y_epi+y_err

    cat("\nRemove unneeded parameters for storage")
    rm(y_marginal); rm(y_epi)

    cat("\nSet threshold:\n")
    thresh=qnorm(1-k,mean=0,sd=1)
    cat(thresh)

    cat("\nSet phenotype")
    pheno=(y>thresh)+1

    cat("\nGet number of affected and unaffected people.\n")
    ncases = sum(pheno==2)
    cat("Prop diseased indivs: ", ncases/length(y),"\n")
    cat("Number of affected:\n")
    cat(ncases) #todo: update variable names to be disease_indivs
    ncontrols = sum(pheno==1)
    cat("Prop normal indivs: ", ncontrols/length(y),"\n")
    cat("Number of unaffected:\n")
    cat(ncontrols) #todo: update variable names to be normal_indivs
    
    # Number diseased individuals to sample
    num_cases = floor(obs*prop_case)
    cat("Number of diseased individuals to sample:\n")
    cat(num_cases)

    # Check that there are enough cases to sample for observations
    if(!(length(which(y>thresh))>floor(num_cases))){
        cat("\nError: There should be enough obs. to sample under threshold\n")
        cat(length(which(y>thresh)))
        cat(">",num_cases,"\n")
        next
    }

    cat("Sample cases and controls given calculated threshold.\n")
    cases = sample(which(y>thresh),num_cases,replace = FALSE)
    controls = sample(which(y<=thresh),floor(obs-num_cases),replace = FALSE)
    y = pheno[c(cases,controls)] # vector of 2s (affected) and 1s (unaffected)
    X = geno[c(cases,controls),] # 
    
    if (hierarchy){
        cat("\nEnsure no duplicate eQTL/SNP names")
        stopifnot("Resolve: duplicate eQTL/SNP names."=length(unique(colnames(X)))==length(colnames(X)))
    }

    cat("\nSampled genotypes (X) across two populations look looks like:\n")
    print(X[1:4,1:4])
    cat("Dimensions: ",dim(X), "\n")

    cat("\nCenter and scale final X (causal and non-causal SNPs/eQTLs)")
    Xmean=apply(X, 2, mean); 
    Xsd=apply(X, 2, sd); 
    Xcs=t((t(X)-Xmean)/Xsd)
    
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
    
        #Save masks unlabeled row and columns
        write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("masksim_1", as.character(ndat),"txt", sep="."),sep="/"), mask1_sim, sep=" ", row.names=FALSE, , col.names=FALSE)
        write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("masksim_2", as.character(ndat),"txt", sep="."),sep="/"), mask2_sim, sep=" ", row.names=FALSE, , col.names=FALSE)
        write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("masksim_3", as.character(ndat),"txt", sep="."),sep="/"), mask3_sim, sep=" ", row.names=FALSE, , col.names=FALSE)

        #Save masks labeled row and columns
        write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("masksim_1_labeled", as.character(ndat),"txt", sep="."),sep="/"), mask1_sim, sep=" ", quote=FALSE)
        write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("masksim_2_labeled", as.character(ndat),"txt", sep="."),sep="/"), mask2_sim, sep=" ", quote=FALSE)
        write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("masksim_3_labeled", as.character(ndat),"txt", sep="."),sep="/"), mask3_sim, sep=" ", quote=FALSE)

        # Save causal pathways, genes, and eQTLs/SNPs
        write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("causal_pathways", as.character(ndat),"txt", sep="."),sep="/"), causal_p$causal_pathways, sep=" ", quote=FALSE,row.names=FALSE, col.names=FALSE)
        write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("causal_genes", as.character(ndat),"txt", sep="."),sep="/"), causal_g$causal_genes, sep=" ", quote=FALSE,row.names=FALSE, col.names=FALSE)
    }
    
    cat("\nSave files")
    write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"), paste("Xsim", as.character(ndat),"txt", sep="."),sep="/"), X, sep=" ", row.names=FALSE, col.names=FALSE)
    write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"), paste("Xsim_cent_scal", as.character(ndat),"txt", sep="."),sep="/"), Xcs, sep=" ", row.names=FALSE, col.names=FALSE)
    write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("ysim", as.character(ndat),"txt", sep="."),sep="/"), y, sep=" ", row.names=FALSE, col.names=FALSE)
    
    cat("\nSave files again but with the names of the SNPs/SNP-sets")
    #colnames(X) = paste("SNP",1:ncol(X),sep="")
    #rownames(X) = paste("sample",1:nrow(X),sep="")
    #rownames(Xcs) = paste("sample",1:nrow(Xcs),sep="")
    #colnames(mask) = paste("SNP_set",1:ncol(mask),sep="")
    #rownames(mask) = colnames(X)

    y_matrix = cbind(rownames(X), y)

    write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("Xsim_labeled", as.character(ndat),"txt", sep="."),sep="/"), X, sep=" ", quote=FALSE)
    write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("Xsim_labeled_cent_scal", as.character(ndat),"txt", sep="."),sep="/"), Xcs, sep=" ", quote=FALSE)
    write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("ysim_labeled", as.character(ndat),"txt", sep="."),sep="/"), y_matrix, sep=" ", quote=FALSE, row.names=FALSE)
    
    # Save causal SNPs/eQTLs
    write.table(file=paste(paste(SIM_OUTDIR,SUB_OUTDIR_NAME,sep="/"),paste("causal_snps", as.character(ndat),"txt", sep="."),sep="/"), c(s1,s2), sep=" ", quote=FALSE,row.names=FALSE, col.names=FALSE)

    # Output number of cases above threshold
    cat("\nNumber of cases above threshold:\n")
    cat(length(which(y>thresh)),"\n")

    ### Check dimensions and add SNP names ###
    cat("\nDimensions of X: ", dim(X), "\n")
    cat("Dimensions of y: ", length(y), "\n")

    cat(paste("\nNumber of SNPs:", tot_snp_sim, "\n"))
    cat(paste("\nNumber of causal SNPs:", causal_s$n_c_snps, "\n"))
    
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

    cat("LTSim complete.")
}
