# lt-sim_sweeps.sh

# Check to make sure input is correct
if [ $# -le 12 ]; then
	echo 'Usage: ./lt-sim_sweeps.sh 
            1) seed                  (integer) 
            2) /PATH/TO/output_dir   (path for generated output)
            3) /PATH/TO/ind          (FILE, number of individuals in the population)
            4) /PATH/TO/k            (FILE, prevalence of the disease/trait in population)
            5) /PATH/TO/pve          (FILE, broad sense heritability)
            6) /PATH/TO/fst          (FILE, targeted fixation index as a measure of population \\ differentiation due to genetic structure (suggest 0.05))
            7) /PATH/TO/obs          (FILE, total observed population)
            8) /PATH/TO/prop_case    (FILE, proportion of case observations) 
            9) /PATH/TO/tot_snp_sim  (FILE, total number SNPs in population)
            10) /PATH/TO/frac_causal (FILE, proportion of causal SNPs in population)
            11) hierarchy            (0:no or 1:yes, sample causal SNPs according to biological annotations (e.g., genes or pathways)
            12) /PATH/TO/mask        (NA:no, or path to mask) 
            13) degree               (0:no or positive integer for degree of each gene connections to SNP(s))
         '
	exit 1
fi


echo "
Generating simulations for 10 replicates per parameter set.
"

# Primary parameters
SEED=$1
SIM_OUTDIR=$2

IND_FILE=$3
K_FILE=$4
PVE_FILE=$5
FST_FILE=$6
OBS_FILE=$7
PROP_CASE=$8
TOT_SNP_SIM=$9
FRAC_CAUSAL=$10

HIER=$11
MASK_DIR=$12
DEG=$13

echo "Number processors: $NUM_PROCESSORS"
echo "K, observations, PVE: $K_FILE, $OBS_FILE, $PVE_FILE"

# Secondary parameters
NUM_SIM_RUNS=10 # Number replicates

# Create output directories as needed
mkdir -p $SIM_OUTDIR
mkdir -p $RES_OUTDIR

# Creating command script for HTSeq
COMMAND_SCRIPT=${SIM_OUTDIR}"/run_sims.txt"

# Removing command script if already made
if [ -f ${COMMAND_SCRIPT} ]; then
    echo -e "REPLACING ${COMMAND_SCRIPT}\n"
    rm -rf ${COMMAND_SCRIPT}
fi

# Setting the counter for the number of processes to run in a batch
START=0
i=$START

# Create grid of parameter search
for ind in $(awk '{print $ind}' $IND_FILE)
do  
    for k in $(awk '{print $k}' $K_FILE)
    do  
        for pve in $(awk '{print $pve}' $PVE_FILE)
        do  
            for fst in $(awk '{print $fst}' $FST_FILE);  
            do
                for obs in $(awk '{print $obs}' $OBS_FILE);  
                do
                    for prop_case in $(awk '{print $prop_case}' $PROP_CASE_FILE);  
                    do
                        for tot_snp_sim in $(awk '{print $tot_snp_sim}' $TOT_SNP_SIM_FILE);  
                        do
                            for frac_causal in $(awk '{print $frac_causal}' $FRAC_CAUSAL_FILE);  
                            do
                                echo "Combination: seed $SEED, k $k, samples $obs, and PVE (H^2) $pve"
                                ((i = i + 1)) # Incrementing counter
                                
                                # If number of processors reached, add wait to form a batch that will finish, and then process the next batch
                                if (( ${i}%${NUM_PROCESSORS}==0 )); then # was 
                                    echo "wait" >> $COMMAND_SCRIPT
                                fi 

                                # conda activate multio_case_ctrl
                                echo "( /usr/bin/time Rscript lt-sim.r $SEED $MASK_DIR $SIM_OUTDIR $NON_OVERLAP_USED $k $obs $pve > $RES_OUTDIR/run_out.seed$SEED-k$k-obs$obs-pve$pve.txt 2>&1) &">> $COMMAND_SCRIPT
                            done
                        done
                    done
                done
            done
        done  
    done
done 

# Run command
echo "LT-Sim prepared your script here: $COMMAND_SCRIPT for your review."
echo "LT-Sim script ready to run with ./$COMMAND_SCRIPT"