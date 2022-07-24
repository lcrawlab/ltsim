# run_sims.sh
# Ashley Conard
# Last Mod. 11/17/21

# Check to make sure input is correct
if [ $# -le 8 ]; then
	echo 'Usage: ./run_sims.sh 
                1) SEED
                2) FULL/PATH/TO/MASK_DIR/
                3) FULL/PATH/TO/SIM_OUTDIR/
                4) FULL/PATH/TO/RES_OUTDIR/
                5) NUM_PROCESSORS
                6) NON_OVERLAP_USED (0 or positive integer for degree)
                7) FULL/PATH/TO/K_FILE 
                8) FULL/PATH/TO/OBS_FILE
                9) FULL/PATH/TO/PVE_FILE
                NOTE: Parameter files can be separated by a space or left empty.
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
            echo "( /usr/bin/time Rscript /data/compbio/aconard/microsoft_research/mit_collab/sparsenn/simulations/lt_sims/lt_simulations.r $SEED $MASK_DIR $SIM_OUTDIR $RES_OUTDIR $NON_OVERLAP_USED $k $obs $pve > $RES_OUTDIR/run_out.seed$SEED-k$k-obs$obs-pve$pve.txt 2>&1) &">> $COMMAND_SCRIPT
        done
    done  
done 

# Run command
# bash $COMMAND_SCRIPT
echo "Script ready to run with ./$COMMAND_SCRIPT"



#################################
# echo "python ../../../BANNs/BANN_numpy/run_banns.py $SIMDIR $BANNS_RES" >> COMMAND_SCRIPT # TODO: get SIMDIR # TODO: make relative to another folder
# echo "python SCHD.py $SIMDIR $SCHD_RES" >> COMMAND_SCRIPT
# echo "python ../roc_curve.py $SIMDIR $OUTDIR" >> COMMAND_SCRIPT



# for ((i=1;i<=NUM_SIM_RUNS;i++))
#     do
#         echo "Combination: seed $SEED, k $k, samples $obs, and PVE (H^2) $pve"
#         #((i = i + 1)) # Incrementing counter
#         ((SEED = SEED + 1)) # Incrementing seed
#         # If number of processors reached, add wait to form a batch that will finish, and then process the next batch
#         if (( ${i}%${NUM_PROCESSORS}==0 )); then # was 
#             echo "wait" >> $COMMAND_SCRIPT
#         fi 

#         # conda activate multio_case_ctrl
#         echo "( /usr/bin/time Rscript /data/compbio/aconard/microsoft_research/mit_collab/sparsenn/simulations/lt_sims/lt_simulations.r $SEED $MASK_DIR $OUTDIR $METHOD_OUTDIR $NON_OVERLAP_USED $k $obs $pve 2> $OUTDIR/run_out.$SEED.txt) &">> $COMMAND_SCRIPT
# done