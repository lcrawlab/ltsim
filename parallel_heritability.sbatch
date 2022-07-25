#!/bin/bash
#SBATCH -n 2 #Request 8 tasks (cores)
#SBATCH -N 1 #Request 1 node
#SBATCH -p sched_mem1TB_centos7  #Run on sched_engaging_default partition
#SBATCH --mem=1GB #Request 1G of memory per CPU
#SBATCH -o ../results/slurm/output_%j.txt #redirect output to output_JOBID.txt
#SBATCH -e ../results/slurm/error_%j.txt #redirect errors to error_JOBID.txt
#SBATCH --exclude=node009,node011
#SBATCH --time 0-1:00

while getopts a:b:c:d: flag
do
    case "${flag}" in

        a) pset=${OPTARG};;
        b) simnum=${OPTARG};;
        c) respath=${OPTARG};;
		d) savepath=${OPTARG};;
    esac
done

python heritability_model_interpretable.py --pset $pset --simnum $simnum --respath $respath --savepath $savepath