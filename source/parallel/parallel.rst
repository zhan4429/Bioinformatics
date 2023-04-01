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


RNAseq by HISAT2
~~~~~~~~~~~~~~~~~~
Index the reference genome
++++++++++++++++++++++++++++
When aligning against a reference genome, we first need to prepare (index) that genome so that our software can read it efficiently. This indexing needs to be done only once; after indexing it once, the index can be reused when aligning the same organism. We even recommend placing the index outside of the project folder; for clarity, we will keep it local::

	# The reference genome.
	IDX=refs/genome.fa

	# Build the genome index.
	hisat2-build $IDX $IDX

	# Index the reference genome with samtools.
	samtools faidx $IDX

Generate the alignments
++++++++++++++++++++++++++++
A single alignment with hisat2 would be generated with::

	hisat2 -x refs/genome.fa -1 reads/BORED_1_R1.fq  -2 reads/BORED_1_R2.fq | head

we can see that it is a SAM file that gets produced::
	
	HD     VN:1.0  SO:unsorted
	@SQ     SN:Golden       LN:128756
	@PG     ID:hisat2       PN:hisat2       VN:2.2.1        CL:"/home/ialbert/miniconda3/envs/bioinfo/bin/hisat2-align-s --wrapper basic-0 -x refs/genome.fa --read-lengths 100 -1 reads/BORED_1_R1.fq -2 reads/BORED_1_R2.fq"
	HWI-ST718_146963544:7:1101:1307:81391   83      Golden  33031   60      100M    =       32827   -304    ACCGCCGCATACGGCCGATTGTCGCAGCCCGGGTCGATTATAACAACGGTGCAATCTCAGCTAAACCGACGCAGTTTTGCTCCTTGGATTCTGAGCCCGG  <BB><5CA<>5B@@B@C>@<@<5?<<>?@BB@>B=C@AADCBB?56HHDGGHIHEGFGEHIGGHFHAGDIGGEIIIIIHCC;AIIIGFHHHHFFFDFCCC                                         AS:i:0   XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100        YS:i:0  YT:Z:CP NH:i:1
	HWI-ST718_146963544:7:1101:1307:81391   163     Golden  32827   60      100M    =       33031   304     CATGCACAGCGGTCAAACAGTATGTCCCAAGGGGACTTAAGCGCGGTGGCCTCCCCTATCCCCTACGAGGCTACCCGGATCGATGACGCGAATTGGGGAC  @C?DFFFFHHHH?FHIGIIIHGHIIIIIGHIGBD@FHIIIIIGGED@EBC?;BDFDEDEEDDDCDD@BDDB?CCCBB59@DBBBAC>BDB9952>CBD05                                         AS:i:0   XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100        YS:i:0  YT:Z:CP NH:i:1
	HWI-ST718_146963544:7:1101:1328:22995   83      Golden  33653   60      100M    =       33534   -219    ATCGCGCACGTAGTAGCCTAGAGCGCCAGGGGCGGAAATTCGCCTGAAAAGTTTTGCCGGCGCACAAGCACGATCGGCTCCTAATAGGAGGTGAATTAGA  BBBBBDDBDDDCCCCCCCCA<-BBBABDDDDDDEDDDDDDDDEECEFFFDBEEHJJJIHEGGGFIGIHGFJIJJJJJJJIJJJJJIJHHHHHFFFFFCCC                                         AS:i:0   XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:100        YS:i:0  YT:Z:CP NH:i:1


To store this output in alignment files for each sample for BORED states, we’ll need to run::

	hisat2 -x refs/genome.fa -1 reads/BORED_1_R1.fq -2 reads/BORED_1_R2.fq | samtools sort > BORED_1.bam
	hisat2 -x refs/genome.fa -1 reads/BORED_2_R1.fq -2 reads/BORED_2_R2.fq | samtools sort > BORED_2.bam
	hisat2 -x refs/genome.fa -1 reads/BORED_3_R1.fq -2 reads/BORED_3_R2.fq | samtools sort > BORED_3.bam

and for EXCITED states, we need to run::

	hisat2 -x refs/genome.fa -1 reads/EXCITED_1_R1.fq -2 reads/EXCITED_1_R2.fq | samtools sort > EXCITED_1.bam
	hisat2 -x refs/genome.fa -1 reads/EXCITED_2_R1.fq -2 reads/EXCITED_2_R2.fq | samtools sort 	> EXCITED_2.bam
	hisat2 -x refs/genome.fa -1 reads/EXCITED_3_R1.fq -2 reads/EXCITED_3_R2.fq | samtools sort > EXCITED_3.bam


How to automate the alignments
++++++++++++++++++++++++++++++++
We could write it out our commands by hand like above, but you can appreciate how error-prone and non-reusable that process is. We’d rather automate the process. Recall that our design.csv file contains::

	sample,condition
	BORED_1, bored
	BORED_2, bored
	BORED_3, bored
	EXCITED_1, excited
	EXCITED_2, excited
	EXCITED_3, excited

from which we created an ids.txt file that contains the first column with no header::

	BORED_1
	BORED_2
	BORED_3
	EXCITED_1
	EXCITED_2
	EXCITED_3

now our commands become::

	# The index name.
	IDX=refs/genome.fa

	# Create the BAM folder.
	mkdir -p bam

	# Align the FASTQ files to the reference genome.
	cat ids.txt | parallel "hisat2 -x $IDX -1 reads/{}_R1.fq -2 reads/{}_R2.fq | samtools sort > bam/{}.bam"

	# Index each BAM file.
	cat ids.txt | parallel  "samtools index bam/{}.bam"



Making a bigWig coverage file
++++++++++++++++++++++++++++++++

Since what we primarily care about is the “coverage” rather than individual alignments, we can turn the BAM file into a so-called BigWig file.

You can’t just make a bigWig; that would be too easy. First, you turn your BAM files into bedGraph, then you take each bedGraph and turn it into a bigWig. Sigh… Thankfully we have it all automated. The following steps will be necessary::

	# Turn each BAM file into bedGraph coverage. The files will have the .bg extension.
	cat ids.txt | parallel "bedtools genomecov -ibam  bam/{}.bam -split -bg  > bam/{}.bg"

	# Convert each bedGraph coverage into bigWig coverage. The files will have the .bw extension.
	cat ids.txt | parallel "bedGraphToBigWig bam/{}.bg  ${IDX}.fai bam/{}.bw"


