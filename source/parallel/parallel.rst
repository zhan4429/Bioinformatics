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


Blast example
~~~~~~~~~~~~~~~~~``
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

