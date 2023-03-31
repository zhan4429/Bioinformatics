.. _backbone-label:

Snakemake
==============================

rule all
~~~~~~~~~~~~~~
The simplest possible Snakemake pipeline is one that produces a single file. To begin with, you need to tell Snakemake the name of the file you want to produce::

	rule all:
		input: "hello_world.txt"

This tells Snakemake that the pipeline should create a file called ``hello_world.txt``. 

Next we can save the content into a file called ``snakefile``. Now run the workflow with::

	snakemake --snakefile snakefile -j 1 

However, you will get such result::
	
	MissingInputException in rule all in file snakefile, line 1:
	Missing input files for rule all:
    	affected files:
      	  hello_world.txt


The problem is that you have only told Snakemake what file to produce, not how to produce it. Snakemake has looked through your snakefile, but not found any instructions for how to produce the file. That is what the error message means.

Now we will tell Snakemake how to produce "hello/world.txt". To do that you need to create a new rule. You can call your rules almost anything, but it is best to choose a descriptive name.

Let's call the rule hello_world. Write it in your snakefile after your all rule::

	rule all:
   		input: "hello_world.txt"
	rule hello_world:
    	output: "hello_world.txt"
    	shell: "echo Hello World > hello_world.txt"


Now when you run the pipeline, it should work::

	Building DAG of jobs...
	Using shell: /usr/local/bin/bash
	Provided cores: 1 (use --cores to define parallelism)
	Rules claiming more threads will be scaled down.
	Job stats:
	job            count    min threads    max threads
	-----------  -------  -------------  -------------
	hello_world        1              1              1
	total              1              1              1

	Select jobs to execute...

	[Fri Mar 31 09:48:22 2023]
	rule hello_world:
	    output: hello_world.txt
	    jobid: 0
	    reason: Missing output files: hello_world.txt
	    resources: tmpdir=/tmp

	[Fri Mar 31 09:48:22 2023]
	Finished job 0.
	1 of 1 steps (100%) done
	Complete log: .snakemake/log/2023-03-31T094822.182463.snakemake.log
