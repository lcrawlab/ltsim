# create_hierarchy.r
# Last mod. 10/9/2022
# Ashley Mae Conard
# Purpose: Create non-overlapping masks that link SNP to gene, depending given degree and total number SNPS.

# Clear console
cat("\014")

# Clear environment
rm(list = ls(all = TRUE))

# Input
args = commandArgs(TRUE)
if (length(args)==0) {
	  stop("Input: OUTDIR, DEGREE", call.=FALSE)
} else if (length(args) < 2) {
	stop("Input: OUTDIR, DEGREE", call.=FALSE)
}

# Argument assignments
OUTDIR <- args[1]
TOT_NUM_SNP <- as.integer(args[2]) #Number of Total SNPs in the Data (5e3)
NDEGREE <- as.integer(args[2]) # the degree of connectivity between SNPs to each gene (e.g. NDEGREE=4, then 4 SNPs map to 1 gene)

# Simulating phenotypes params
SUB_OUTDIR_NAME = paste("masks_non_overlap_degree",NDEGREE,"tot_snp_sim",TOT_NUM_SNP, sep="-")
#####################################################

# Create output dirs if needed
path = paste(OUTDIR,SUB_OUTDIR_NAME,sep="/")
dir.create(path, showWarnings = TRUE, recursive = TRUE, mode = "0777")

# 1) Create mask 1
nsnp_sets= round(TOT_NUM_SNP/NDEGREE) #Number of SNP-sets
mask1_pd=matrix(0, TOT_NUM_SNP, nsnp_sets) #Initialize annotation mask as a matrix of zeros of size nSNPs by nSets

j_cnt = 1
i_cnt = 1
for(j in 0:nsnp_sets){
    for(i in (j_cnt:(j+NDEGREE))){
        mask1_pd[i_cnt,j_cnt]=1
        i_cnt = i_cnt+1
        if(i_cnt > TOT_NUM_SNP){
            break
        }
    }
    j_cnt=j_cnt+1
    if(j_cnt > nsnp_sets){
        break
    }
    
}
snps_enum = 1:TOT_NUM_SNP
genes_enum = 1:nsnp_sets

nsnp_sets = dim(mask1_pd)[2]

s_names <- paste("s", snps_enum, sep="")
g_names <- paste("g", genes_enum, sep="")

rownames(mask1_pd) <- s_names
colnames(mask1_pd) <- g_names

# 2) Create mask 2
npathways=round(nsnp_sets/NDEGREE) #Number of SNP-sets
mask2_pd=matrix(0, nsnp_sets, npathways) #Initialize annotation mask as a matrix of zeros of size nSNPs by nSets

j_cnt = 1
i_cnt = 1
for(j in 0:nsnp_sets){
    for(i in (j_cnt:(j+NDEGREE))){
        mask2_pd[i_cnt,j_cnt]=1
        i_cnt = i_cnt+1
        if(i_cnt > nsnp_sets){
            break
        }
    }
    j_cnt=j_cnt+1
    if(j_cnt > npathways){
        break
    }
}

pathways_enum = 1:npathways
p_names <- paste("p", pathways_enum, sep="")

rownames(mask2_pd) <- g_names
colnames(mask2_pd) <- p_names

# 3) Create mask3 of SNP to pathway. Map SNP to gene mask with gene to pathway mask
mask3_pd=matrix(0, TOT_NUM_SNP, npathways) #Initialize annotation mask as a matrix of zeros of size nSNPs by nSets
rownames(mask3_pd) = s_names
colnames(mask3_pd) = p_names

for(i in 1:TOT_NUM_SNP) {
    idx = which(mask1_pd[i,]==1)
    gene_name = colnames(mask1_pd)[idx][1]

    # Find pathway_name mapping to gene_name
    idx = which(mask2_pd[gene_name,]==1)
    pathway_name = colnames(mask2_pd)[idx]

    # Assign 1 to SNP assignment to pathway
    mask3_pd[i, pathway_name] = 1
}

# Save masks
write.table(file=paste(OUTDIR,SUB_OUTDIR_NAME,"masksim_1_non_overlap.txt",sep="/"), mask1_pd, sep=" ", row.names=FALSE, col.names=FALSE)
write.table(file=paste(OUTDIR,SUB_OUTDIR_NAME,"masksim_2_non_overlap.txt",sep="/"), mask2_pd, sep=" ", row.names=FALSE, col.names=FALSE)
write.table(file=paste(OUTDIR,SUB_OUTDIR_NAME,"masksim_3_non_overlap.txt",sep="/"), mask3_pd, sep=" ", row.names=FALSE, col.names=FALSE)

write.table(file=paste(OUTDIR,SUB_OUTDIR_NAME,"masksim_1_non_overlap_labeled.txt",sep="/"), mask1_pd, sep=" ", quote=FALSE)
write.table(file=paste(OUTDIR,SUB_OUTDIR_NAME,"masksim_2_non_overlap_labeled.txt",sep="/"), mask2_pd, sep=" ", quote=FALSE)
write.table(file=paste(OUTDIR,SUB_OUTDIR_NAME,"masksim_3_non_overlap_labeled.txt",sep="/"), mask3_pd, sep=" ", quote=FALSE)