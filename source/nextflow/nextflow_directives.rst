.. _backbone-label:

Nextflow-Directives
==============================
Directive declarations allow the definition of optional settings, like the number of cpus and amount of memory, that affect the execution of the current process without affecting the task itself.
They must be entered at the top of the process body, before any other declaration blocks (i.e. input, output, etc).

Note: You do not use ``=`` when assigning a value to a directive.

For example, the process below uses the three directives, ``tag``, ``cpus`` and ``label``. ::


	//process_index.nf
	nextflow.enable.dsl=2

	process BWA_INDEX {

	  tag {"BWA_INDEX $genome"}
	  label 'process_low'
	  cpus 1

	  input:
	  path genome

	  output:
	  path("*")

	  script:
	  """
	  bwa index ${genome}
	  """
	}

	ref_ch = Channel.fromPath("/workspace/nextflow_tutorial/data/ref_genome/ecoli_rel606.fasta")  

	workflow {
	  
	  BWA_INDEX( ref_ch )
	  
	}


The ``tag`` directive to allow you to give a custom tag to each process execution. This tag makes it easier to identify a particular task (executed instance of a process) in a log file or in the execution report.

The second directive ``label`` allows the annotation of processes with mnemonic identifier of your choice. Labels are useful to organise workflow processes in separate groups which can be referenced in the configuration file to select and configure subset of processes having similar computing requirements.

The third directive ``cpus`` allows you to define the number of CPUs required for each task.



