.. _backbone-label:

Nextflow-outputs
==============================

The output declaration block allows us to define the channels used by the process to send out the files and values produced.
An output block is not required, but if it is present it can contain one or more outputs declarations.

The output block follows the syntax shown below::

	output:
	 <output qualifier> <output name>

Output values
~~~~~~~~~~~~~~
The type of output data is defined using output qualifiers.
The ``val`` qualifier allows us to output a value defined in the script.
If we want to capture a file instead of a value we can use the ``path`` qualifier that can capture one or more files produced by the process, over the specified channel.
Because Nextflow processes can only communicate through channels if we want to share a value input into one process as input to another process we would need to define that value in the output declaration block.
When an output file name contains a ``*`` or ``?`` character it is interpreted as a pattern match. This allows to capture multiple files into a list and output them as a one item channel.
Since all the files produced by the process are captured using ``"*"`` in the output block, when the task is completed all the output files are sent over the output channel. A downstream ``operator`` or ``process`` declaring the same channel as ``input`` will be able to receive it.

Grouped inputs and outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~
So far we have seen how to declare multiple input and output channels, but each channel was handling only one value at time. However Nextflow can handle groups of values using the ``tuple`` qualifiers.
In tuples the first item is the grouping key and the second item is the list of files::

	[group_key,[file1,file2,...]]

When using channel containing a ``tuple``, such a one created with ``.filesFromPairs`` factory method, the corresponding input declaration must be declared with a ``tuple`` qualifier, followed by definition of each item in the tuple.

Create a new file ``process_tuple_input.nf``; add the following and ``nextflow run process_tuple_input.nf -process.echo``::

	//process_tuple_input.nf
	nextflow.enable.dsl=2

	process TUPLEINPUT{
	  input:
	  tuple val(sample_id), path(reads)
	  
	  script:
	  """
	  echo ${sample_id}
	  echo ${reads}
	  """
	}

	reads_ch = Channel.fromFilePairs("/workspace/nextflow_tutorial/data/trimmed_fastq/SRR2584863_{1,2}.trim.fastq.gz")

	workflow {
	  
	  TUPLEINPUT(reads_ch)
	  
	}