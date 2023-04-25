.. _backbone-label:

Nextflow-inputs
==============================
You can only define one input block at a time and it must contain one or more inputs declarations. The input block follows the syntax shown below::

	input:
	 <input qualifier> <input name>

The input qualifier declares the type of data to be received. Types of input qualifiers:

- val: Lets you access the received input value by its name in the process script.
- env: Lets you use the received value to set an environment variable named as > the specified input name.
- path: Lets you handle the received value as a path, staging the file properly in the execution context.
- stdin: Lets you forward the received value to the process stdin special file.
- tuple: Lets you handle a group of input values having one of the above qualifiers.
- each: Lets you execute the process for each entry in the input collection.


Input values
~~~~~~~~~~~~~~
The ``val`` qualifier allows you to receive value data as input. It can be accessed in the process script by using the specified input name, as shown in the following example::

	//process_input_value.nf
	nextflow.enable.dsl=2

	process PRINTCHR {

	  input:
	  val chr

	  script:
	  """
	  echo processing chromosome ${chr}
	  """
	}

	chr_ch = Channel.of( 1..22,'X','Y' )

	workflow {

	  PRINTCHR( chr_ch )
	}

Input files
~~~~~~~~~~~~~

When you need to handle files as input you need the path qualifier.
Using the ``path`` qualifier means that Nextflow will stage it in the process execution directory, and it can be accessed in the script by using the name specified in the input declaration.
The input file name can be defined dynamically by defining the input name as a Nextflow variable and referenced in the script using the ``$variable_name`` syntax.
For example in the script below we assign the variable name ``genome`` to the input files using the ``path`` qualifier. The file is referenced using the variable substitution syntax ``${genome}`` in the script block::

	//process_input_file.nf
	nextflow.enable.dsl=2

	/*
	 * Index the reference genome for use by bwa and samtools.
	 */
	process BWA_INDEX {

	  input:
	  path genome

	  script:
	  """
	  bwa index ${genome}
	  """
	}

	ref_ch = Channel.fromPath("/workspace/nextflow_tutorial/data/ref_genome/ecoli_rel606.fasta")  

	workflow {
	  BWA_INDEX( ref_ch )
	}

