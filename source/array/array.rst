.. _backbone-label:

https://hpc-unibe-ch.github.io/slurm/array-jobs.html
Array jobs
==============================
Job arrays are used for running the same job a large number of times with only slight differences between the jobs. For instance, let's say that you need to run 100 jobs, each with a different seed value for the random number generator. Or maybe you want to run the same analysis script on data for each of the 50 states in the USA. Job arrays are the best choice for such cases.

Below is an example Slurm script for running Python where there are 5 jobs in the array::

	#!/bin/bash
	#SBATCH --job-name=array-job     # create a short name for your job
	#SBATCH --output=slurm-%A.%a.out # stdout file
	#SBATCH --error=slurm-%A.%a.err  # stderr file
	#SBATCH --nodes=1                # node count
	#SBATCH --ntasks=1               # total number of tasks across all nodes
	#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
	#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
	#SBATCH --time=00:01:00          # total run time limit (HH:MM:SS)
	#SBATCH --array=0-4              # job array with index values 0, 1, 2, 3, 4
	#SBATCH --mail-type=all          # send email on job start, end and fault
	#SBATCH --mail-user=<YourNetID>@princeton.edu

	echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
	echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
	echo "Executing on the machine:" $(hostname)

	module purge
	module load anaconda3/2022.5
	conda activate myenv

	python myscript.py


The key line in the Slurm script above is::
	#SBATCH --array=0-4




blast array
~~~~~~~~~~~~~
This is previously used blast arry script::

	#!/bin/bash
	#SBATCH --account=plantpath
	#SBATCH --qos=plantpath
	#SBATCH --job-name=blast_array
	#SBATCH --mail-type=FAIL,END
	#SBATCH --mail-user=YOUREMAIL
	#SBATCH --output <blastp_%j.log>
	#SBATCH --ntasks=1
	#SBATCH --cpus-per-task=4
	#SBATCH --mem=8gb
	#SBATCH --time=4:00:00
	#SBATCH --array=1-100
	date;hostname;pwd
	 
	module load ncbi_blast
	 
	export INPUT_DIR="input"
	export OUTPUT_DIR="output"
	export LOG_DIR="logs"
	mkdir -p ${OUTPUT_DIR} ${LOG_DIR}
	 
	RUN_ID=$(( $SLURM_ARRAY_TASK_ID + 1 ))
	 
	QUERY_FILE=$( ls ${INPUT_DIR} | sed -n ${RUN_ID}p )
	QUERY_NAME="${QUERY_FILE%.*}"
	 
	QUERY="${INPUT_DIR}/${QUERY_FILE}"
	OUTPUT="${OUTPUT_DIR}/${QUERY_NAME}.out"
	 
	echo -e "Command:\nblastn –query ${QUERY} –db nt –out ${OUTPUT} –evalue 0.001 –outfmt 6 –num_threads 8"
	 
	blastn -query ${QUERY} -db nt -out ${OUTPUT} -evalue 0.001 -outfmt 6 -num_threads 8
	 
	date

