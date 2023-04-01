.. _backbone-label:

GNU Parallel
==============================
GNU Parallel is a shell tool for executing jobs in parallel on one or multiple computers. It's a helpful tool for automating the parallelization of multiple (often serial) jobs, in particular allowing one to group jobs into a single SLURM submission to take advantage of the multiple cores on a given Savio node.

A job can be a single core serial task, multi-core job, or MPI application. A job can also be a command that reads from a pipe. The typical input is a list of input parameters needed as input for all the jobs collectively. GNU parallel can then split the input and pipe it into commands in parallel. GNU parallel makes sure output from the commands is the same output as you would get had you run the commands sequentially, and output names can be easily tied to input file names for simple post-processing. This makes it possible to use output from GNU parallel as input for other programs.

Basic usage
~~~~~~~~~~~~~~~~~~~~~~
Reading command arguments on the command line::

	parallel    [-j N] [OPTIONS]    COMMAND {} ::: TASKLIST

The ``-j <N>`` option permits to define the jobs per machine. 

We can use multiple inputs per task, distinguishing the inputs by ``{1}, {2}``, etc.::

	parallel --link -j 2 cp file{1}.in file{2}.out ::: 1 2 3 ::: 4 5 6
	ls file*out
	# file4.out  file5.out  file6.out

Note that ``--link`` is needed so that 1 is paired with 4, 2 with 5, etc., instead of doing all possible pairs amongst the two sets.


Combinatorials
~~~~~~~~~~~~~~~~
You can keep adding ``:::`` and ``::::`` to add additional arguments, and these will be combined to generate all possible combinations. This is extremely useful for testing commands with different combinations of input parameters::

	parallel --dry-run -k -j 4 Rscript run_analysis.R {1} {2} ::: `seq 1 2` ::: A B C

The output will look like this::

	Rscript run_analysis.R 1 A
	Rscript run_analysis.R 1 B
	Rscript run_analysis.R 1 C
	Rscript run_analysis.R 2 A
	Rscript run_analysis.R 2 B
	Rscript run_analysis.R 2 C

Blast example
~~~~~~~~~~~~~~~~~
Here's our example task list, task.lst::

	../blast/data/protein1.faa
	../blast/data/protein2.faa

Here's the script we want to run BLAST on a single input file, run-blast.sh::

	#!/bin/bash
	blastp -query $1 -db ../blast/db/img_v400_PROT.00 -out $2  -outfmt 7 -max_target_seqs 10 -num_threads $3

Now let's use GNU parallel in the context of a SLURM job script::

	#!/bin/bash
	#SBATCH --job-name=job-name
	#SBATCH --account=account_name
	#SBATCH --partition=partition_name
	#SBATCH --nodes=2
	#SBATCH --cpus-per-task=2
	#SBATCH --time=2:00:00

	## Command(s) to run (example):
	module load bio/blast/2.6.0
	module load gnu-parallel/2019.03.22

	export WDIR=/your/desired/path
	cd $WDIR

	# set number of jobs based on number of cores available and number of threads per job
	export JOBS_PER_NODE=$(( $SLURM_CPUS_ON_NODE / $SLURM_CPUS_PER_TASK ))

	echo $SLURM_JOB_NODELIST |sed s/\,/\\n/g > hostfile

	parallel --jobs $JOBS_PER_NODE --slf hostfile --wd $WDIR --joblog task.log --resume --progress -a task.lst sh run-blast.sh {} output/{/.}.blst $SLURM_CPUS_PER_TASK


Some things to notice:
- Here BLAST will use multiple threads for each job, based on the SLURM_CPUS_PER_TASK variable that is set based on the -c (or --cpus-per-task) SLURM flag.
- We programmatically determine how many jobs to run on each node, accounting for the threading.
- Setting the working directory with --wd is optional; without that your home directory will be used (if using multiple nodes via --slf) or the current working directory will be used (if using one node).
- The --resume and --joblog flags allow you to easily restart interrupted work without redoing already completed tasks.
- The --progress flag causes a progress bar to be displayed.
- In this case, only one of the three inputs to run-blast.sh is provided in the task list. The second argument is determined from the first, after discarding the path and file extension, and the third is constant across tasks.


Bowtie2
~~~~~~~~~~~~~
We can design a bash script with the bowtie2 command for each pair separately . Obviously if dealing with hundreds of files this could become cumbersome::

	module load bowtie2
	bowtie2 --threads 4 -x tair -k1 -q -1 SRR4420293_1.fastq.gz -2 SRR4420293_2.fastq.gz -S first_R1.sam >& first.log
	..........................
	bowtie2 --threads 4 -x tair -k1 -q -1 SRR4420295_1.fastq.gz -2 SRR4420295_2.fastq.gz -S fifth_R1.sam >& third.log

GNU parallel let’s us automate this task by using a combination of substitution and separators notably ``:::`` and ``:::+``. We can also make optimum use of the available threads::

	module load bowtie2
	time parallel -j2 "bowtie2 --threads 4 -x tair -k1 -q -1 {1} -2 {2} -S {1/.}.sam >& {1/.}.log" ::: fastqfiles/*_1.fastq.gz :::+ fastqfiles/*_2.fastq.gz


You may have also noticed some new syntax where we are using {1}, {2}, {1/.} and {2/.}. The {1} will give us the first file and {2} the second from the list taking two at a time. The {1/.} and {2/.} will take the prefix of the file name before the “.”

