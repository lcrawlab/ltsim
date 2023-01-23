# lt-sim_sweeps.sh
# Ashley Conard
# Last Mod. 11/17/21

# Check to make sure input is correct
if [ $# -le 12 ]; then
	echo 'Usage: ./lt-sim_sweeps.sh 
            1) seed                  (integer) 
            2) /PATH/TO/mask         (path to mask or type NA) 
            3) /PATH/TO/output_dir   (path for generated output)
            4) degree                (0 or positive integer for degree)
            5) ind                   (number of individuals in the population)
            6) k                     (prevalence of the disease/trait in population)
            7) pve                   (broad sense heritability)
            8) fst                   (targeted fixation index as a measure of population \\ differentiation due to genetic structure (suggest 0.05))
            9) obs                   (total observed population)
            10) prop_case            (proportion of case observations) 
            11) tot_snp_sim          (total number SNPs in population)
            12) frac_causal          (proportion of causal SNPs in population)
            13) hierarchy            (sample causal SNPs according to biological annotations (e.g., genes or pathways)
         '
	exit 1
fi


echo "
Generating simulations for 10 replicates per parameter set.
"

# Primary parameters
SEED=$1
MASK_DIR=$2
SIM_OUTDIR=$3
RES_OUTDIR=$4
NUM_PROCESSORS=$5
NON_OVERLAP_USED=$6
K_FILE=$7
OBS_FILE=$8
PVE_FILE=$9


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
for k in $(awk '{print $k}' $K_FILE)
do  
    for obs in $(awk '{print $obs}' $OBS_FILE)
    do  
        for pve in $(awk '{print $pve}' $PVE_FILE);  
        do
            echo "Combination: seed $SEED, k $k, samples $obs, and PVE (H^2) $pve"
            ((i = i + 1)) # Incrementing counter
            #((SEED = SEED + 1)) # Incrementing seed
            # If number of processors reached, add wait to form a batch that will finish, and then process the next batch
            if (( ${i}%${NUM_PROCESSORS}==0 )); then # was 
                echo "wait" >> $COMMAND_SCRIPT
            fi 

            # conda activate multio_case_ctrl
            echo "( /usr/bin/time Rscript lt-sim.r $SEED $MASK_DIR $SIM_OUTDIR $RES_OUTDIR $NON_OVERLAP_USED $k $obs $pve > $RES_OUTDIR/run_out.seed$SEED-k$k-obs$obs-pve$pve.txt 2>&1) &">> $COMMAND_SCRIPT
        done
    done  
done 

# Run command
echo "LT-Sim prepared your script here: $COMMAND_SCRIPT for your review."
echo "LT-Sim script ready to run with ./$COMMAND_SCRIPT"