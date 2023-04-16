.. _backbone-label:

Running Nextflow in HPC
==============================

Disadvantages of Nextflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Nextflow makes a lot of duplicate and intermediate file copies as it operates; processing large amounts of data can easily exhaust storage quotas
- Nextflow uses parallel file systems to synchronize tasks, which is a frequent source of undesired behavior

How to use Nextflow at Purdue RCAC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wrapped Nextflow submission to Slurm
+++++++++++++++++++++++++++++++++++++

The easiest method is placing the ``nextflow run ..`` command into a batch script and submitting it to Slurm with ``sbatch``. The manager process will run on the allocated compute node, and all tasks are configured to use the local executor.

The major benefit of this method, besides simplicity, is only the initial submission waits in a Slurm queue; it is a good pattern for a workflow which includes a very large number of small tasks. One should not combine this method with the Nextflow ``Slurm`` executor because the job running the Nextflow manager is likely to end before the requested task is finished waiting in a queue.

There are two significant caveats to running the Nextflow workflow process directly inside a Slurm job allocation:

- The Nextflow working directory must be placed on the $SCRATCH filesystem. Nextflow uses a file locking feature not available on any of the other filesystems.
- The workflow cannot run longer than the maximum wall time available to a single job in the Slurm QOS being used. This can be partially mitigated by using multiple Slurm submissions in series and passing the -resume flag to Nextflow, but only progress for completely finished tasks will be preserved from one submission to the next.

Nextflow Submits Tasks as Slurm Jobs
+++++++++++++++++++++++++++++++++++++


Basic concepts
~~~~~~~~~~~~~~~
Nextflow is a reactive workflow framework and a programming DSL that eases the writing of data-intensive computational pipelines.

It is designed around the idea that the Linux platform is the lingua franca of data science. Linux provides many simple but powerful command-line and scripting tools that, when chained together, facilitate complex data manipulations.

Nextflow extends this approach, adding the ability to define complex program interactions and a high-level parallel computational environment based on the dataflow programming model.

DSL2 syntax
~~~~~~~~~~~~~~~~
Nextflow (version > 20.07.1) provides a revised syntax to the original DSL, known as DSL2. The DSL2 syntax introduces several improvements such as modularity (separating components to provide flexibility and enable reuse), and improved data flow manipulation. This further simplifies the writing of complex data analysis pipelines, and enhances workflow readability, and reusability.

This feature is enabled by the following directive at the beginning a workflow script::

	nextflow.enable.dsl=2

Scripts that contain the directive ``nextflow.preview.dsl=2`` use an early version of the DSL2 syntax, which may include experimental features that have been changed or removed in the formal DSL2 syntax. Scripts without these directives use the first version of the Nextflow syntax which we refer to as DSL1. DSL1 workflows use many of the same concepts presented in this lesson, but some aspects such as the flow of data are written differently. DSL1 workflows are also written in a single script, unlike DSL2 workflows which can be spread across many files.



Processes, channels, and workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Nextflow workflows have three main parts: processes, channels, and workflows. Processes describe a task to be run. A process script can be written in any scripting language that can be executed by the Linux platform (Bash, Perl, Ruby, Python, etc.). Processes spawn a task for each complete input set. Each task is executed independently, and cannot interact with another task. The only way data can be passed between process tasks is via asynchronous queues, called channels.

Processes define inputs and outputs for a task. Channels are then used to manipulate the flow of data from one process to the next. The interaction between processes, and ultimately the pipeline execution flow itself, is then explicitly defined in a workflow section.


Process
+++++++++++++++++

In Nextflow, a process is the basic processing primitive to execute a user script.

The process definition starts with the keyword ``process``, followed by process name and finally the process body delimited by curly brackets. The process body must contain a string which represents the command or, more generally, a script that is executed by it. 

Syntax::

	process < name > {

	   [ directives ]

	   input:
	    < process inputs >

	   output:
	    < process outputs >

	   when:
	    < condition >

	   [script|shell|exec]:
	   < user script to be executed >

	}


Workflow execution
++++++++++++++++++
While a process defines what command or script has to be executed, the executor determines how that script is actually run in the target system.

If not otherwise specified, processes are executed on the local computer. The local executor is very useful for pipeline development, testing, and small scale workflows, but for large scale computational pipelines, a High Performance Cluster (HPC) or Cloud platform is often required.

Nextflow provides a separation between the pipeline’s functional logic and the underlying execution platform. This makes it possible to write a pipeline once, and then run it on your computer, compute cluster, or the cloud, without modifying the workflow, by defining the target execution platform in a configuration file.

Nextflow provides out-of-the-box support for major batch schedulers and cloud platforms such as Sun Grid Engine, SLURM job scheduler, AWS Batch service and Kubernetes. 




Files
~~~~~~~~~~~~~~~~~
There are several ways to pass file inputs into a process. Nextflow will create a symlink to the original file in the relevant task execution directory when a process is handling an input file through the methods described below. The file can then be accessed by the script using the name specified in the input declaration::

	proteins = Channel.fromPath('data/*.txt')
	process catThemAll {
	    input:
	    file query_file from proteins
	    shell:
	    "cat ${query_file}"
	}

Each file ending with .txt in the data/ directory will be processed by a separate task generated by the catThemAll process. Each task can be identified by a unique task name (``workflow:catThemAll (2)``) or task_id hash (eg. ``7a/7b3084``). Nextflow will stage a symlink to the specific file that the task will process in the tasks’ execution directory, eg. ``/work/tasks/7a/7b3084/.file-2.txt``.

Another way to process files is using the path qualifier. Both file and path qualifiers are similar, except that the former expects file objects, whereas the latter can also interpret strings as the path of the input file. Note that when using raw strings, the path qualifier does not interpret special characters (eg. wildcards), so this syntax works best if you know the absolute string path of your file. Here is the path qualifier in action::

	process catThemAll {
	    input:
	    path x1 from file('data/example-*.txt')
	    path x2 from 'file:///absolute/path/to/working-dir/data/ids.txt'
	    shell:
	    """
	    cat ${x1}
	    cat ${x2}
	    """
	}


Assume you had two files under your data/ directory: example-1.txt and example-2.txt. Will the above process spawn one task or two? In this case, the path-string qualifier for x2 behaves more like a Value Channel, allowing the process to consume it infinitely many times. As a result the process will spawn tasks until x1 runs out of files, resulting in two tasks.

Outputs
~~~~~~~~~~~~
Similar to Process Inputs, the Output block of a Process defines to which Channels the Process should send out the results. For example::

	customer_ids = Channel.from(1, 2, 3, 4)
	process get_data_for_ids {
	    input:
	    val id from customer_ids
	    output:
	    file data_for_id_*' into data_for_ids
	    shell:
	    '''
	    echo !{id} > data_for_id_!{id}.txt
	    '''
	}
	data_for_ids.view()
	>> /path/to/dir/work/36/1ecd790e4eeb3a786c2e5e288b/data_for_id_3.txt
	>> /path/to/dir/work/aa/19b1ac052387cd05ab04021f5f/data_for_id_2.txt
	>> /path/to/dir/work/6e/a5ce292656c2c8802108293c97/data_for_id_4.txt
	>> /path/to/dir/work/89/0efb6979d54a412274f0ad685d/data_for_id_1.txt


In the above example, the get_data_for_ids process sends the files generated by the shell command into the data_for_ids channel, which downstream processes can then consume. As with Inputs, anything from values to files to stdout can be output to the channel. 




Channel
~~~~~~~~~~
In Nextflow there are two kinds of channels: queue channels and value channels.
Queue channel
++++++++++++++
A queue channel is a non-blocking unidirectional FIFO queue which connects two processes, channel factories, or operators.

A queue channel is usually created using a factory method (``_channel-of``, ``_channel-path``, etc) or chaining it with a channel operator (map, flatMap, etc). Queue channels are also created by process output declarations.

Value channel
~~~~~~~~~~~~~~~
A value channel a.k.a. singleton channel is bound to a single value and can be read an unlimited number of times without consuming its content.  

A value channel is created using the value factory method or by operators returning a single value, such as ``first``, ``last``, ``collect``, ``count``, ``min``, ``max``, ``reduce``, ``sum``, etc.


A value channel is implicitly created by a process when it is invoked with a simple value. Furthermore, a value channel is also implicitly created as output for a process whose inputs are all value channels.

For example::
	
	process foo {
	  input:
	  val x

	  output:
	  path 'x.txt'

	  """
	  echo $x > x.txt
	  """
	}

	workflow {
	  result = foo(1)
	  result.view { "Result: ${it}" }
	}


In the above example, since the ``foo`` process is invoked with a simple value instead of a channel, the input is implicitly converted to a value channel, and the output is also provided as a value channel.


flat
+++++++

When true the matching files are produced as sole elements in the emitted tuples (default: false)::

	ch_reads = Channel
	    .fromFilePairs(params.read_path + '/**{1,2}.f*q*', flat: true)

	ch_reads.view()

The output will be like this::

	$ ./main.nf --read_path /data/reads/
	N E X T F L O W  ~  version 19.09.0-edge
	Launching `./main.nf` [elegant_volta] - revision: bb88634790
	WARN: DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE
	[SRR1950773, /data/reads/SRR1950773_1.fastq.gz, /data/reads/SRR1950773_2.fastq.gz]
	[SRR1950772, /data/reads/SRR1950772_1.fastq.gz, /data/reads/SRR1950772_2.fastq.gz]


If ``flat: false``, the output will be like this::

	[SRR1950773, [/data/reads/SRR1950773_1.fastq.gz, /data/reads/SRR1950773_2.fastq.gz]]
	[SRR1950772, [/data/reads/SRR1950772_1.fastq.gz, /data/reads/SRR1950772_2.fastq.gz]]


Your first script
~~~~~~~~~~~~~~~~~~~~
This is a Nextflow script. It contains:

- An optional interpreter directive (“Shebang”) line, specifying the location of the Nextflow interpreter.
- ``nextflow.enable.dsl=2`` to enable DSL2 syntax.
- A multi-line Nextflow comment, written using C style block comments, followed by a single line comment.
- A pipeline parameter ``params.input`` which is given a default value, of the relative path to the location of a compressed fastq file, as a string.
- An unnamed ``workflow`` execution block, which is the default workflow to run.
- A Nextflow channel used to read in data to the workflow.
- A call to the process ``NUM_LINES``.
- An operation on the process output, using the channel operator ``view()``.
- A Nextflow process block named ``NUM_LINES``, which defines what the process does.
- An ``input`` definition block that assigns the input to the variable read, and declares that it should be interpreted as a file path.
- An ``output`` definition block that uses the Linux/Unix standard output stream stdout from the script block.
- A script block that contains the bash commands  ``printf '${read}'`` to print the name of the read file, and ``gunzip -c ${read}	wc -l`` to count the number of lines in the gzipped read file.

/scratch/negishi/zhan4429/biocontainers/SRR23043636_1.fastq.gz

The contents of ``wc.nf``::


	#!/usr/bin/env nextflow

	nextflow.enable.dsl=2

	/*  Comments are uninterpreted text included with the script.
	    They are useful for describing complex parts of the workflow
	    or providing useful information such as workflow usage.

	    Usage:
	       nextflow run wc.nf --input <input_file>

	    Multi-line comments start with a slash asterisk /* and finish with an asterisk slash. */
	//  Single line comments start with a double slash // and finish on the same line

	/*  Workflow parameters are written as params.<parameter>
	    and can be initialised using the `=` operator. */
	params.input = "data/yeast/reads/ref1_1.fq.gz"

	//  The default workflow
	workflow {

	    //  Input data is received through channels
	    input_ch = Channel.fromPath(params.input)

	    /*  The script to execute is called by its process name,
	        and input is provided between brackets. */
	    NUM_LINES(input_ch)

	    /*  Process output is accessed using the `out` channel.
	        The channel operator view() is used to print
	        process output to the terminal. */
	    NUM_LINES.out.view()
	}

	/*  A Nextflow process block
	    Process names are written, by convention, in uppercase.
	    This convention is used to enhance workflow readability. */
	process NUM_LINES {

	    input:
	    path read

	    output:
	    stdout

	    script:
	    /* Triple quote syntax """, Triple-single-quoted strings may span multiple lines. The content of the string can cross line boundaries without the need to split the string in several pieces and without concatenation or newline escape characters. */
	    """
	    printf '${read} '
	    gunzip -c ${read} | wc -l
	    """
	}

You should see output similar to this::

	zhan4429@login01.negishi:[nextflow] $ nextflow run wc.nf 
	N E X T F L O W  ~  version 23.04.0
	Launching `wc.nf` [sad_neumann] DSL2 - revision: 9fba5dcc47
	executor >  local (1)
	[5c/065cb0] process > NUM_LINES (1) [100%] 1 of 1 ✔
	SRR23043636_1.fastq.gz 112903040

Nextflow scripting
~~~~~~~~~~~~~~~~~~~
Nextflow is a Domain Specific Language (DSL) implemented on top of the Groovy programming language, which in turn is a super-set of the Java programming language. This means that Nextflow can run any Groovy and Java code. It is not necessary to learn Groovy to use Nextflow DSL but it can be useful in edge cases where you need more functionality than the DSL provides.

Comments
++++++++++
When we write any code it is useful to document it using comments. In Nextflow comments use the same syntax as in the C-family programming languages ::

	// This is a single line comment. Everything after the // is ignored.

	/*
	   Comments can also
	   span multiple
	   lines.
	 */

Multi-line strings
+++++++++++++++++++
A block of text that span multiple lines can be defined by delimiting it with triple single ''' or double quotes """::

	text = """
	    This is a multi-line string
	    using triple quotes.
	    """

String interpolation
++++++++++++++++++++++
To use a variable inside a single or multi-line double quoted string "" prefix the variable name with a $ to show it should be interpolated::

	

Lists
+++++++++++++++++=
To store multiple values in a variable we can use a List. A List (also known as array) object can be defined by placing the list items in square brackets and separating items by commas ,::

	kmers = [11,21,27,31]

You can access a given item in the list with square-bracket notation []. These positions are numbered starting at ``0``, so the first element has an index of ``0``::

	kmers = [11,21,27,31]
	println(kmers[0])


We can use negative numbers as indices in Groovy. They count from the end of the list rather than the front: the index -1 gives us the last element in the list, -2 the second to last, and so on. Because of this, kmers[3] and kmers[-1] point to the same element in our example list::

	kmers = [11,21,27,31]
	//Lists can also be indexed with negative indexes
	println(kmers[3])
	println(kmers[-1])

The output::

	31
	31




FastQC
++++++++++++
The simple example for FastQC::

sample_ch=Challel.fromPath('data/sample.fastq')

process FASTQC {
	input: 
		file reads from sample_ch
	output:
		file 'fastqc_logs' into fastqc_ch

	scritp:
	"""
	mkdir fastqc_logs
	fastqc -o fastqc_logs -f fastq -q ${reads}
}






https://training.nextflow.io/basic_training/config/

process {
  executor="slurm"
  clusterOptions="--account=cancercenter-dept --qos=cancercenter-dept-b"
}

apptainer.autoMounts=true
apptainer.enabled=true
apptainer.libraryDir='/home/jobrant/.apptainer/cache/library'
apptainer.cacheDir='/home/jobrant/.apptainer/cache/'




process {
  executor='slurm'
  queueSize = 15
  pollInterval = '5 min'
  dumpInterval = '6 min'
  queueStatInterval = '5 min'
  exitReadTimeout = '13 min'
  killBatchSize = 30
  submitRateLimit = '20 min'
  clusterOptions = '-q debug -t 00:30:00 -C haswell'
}  


Configuration file
~~~~~~~~~~~~~~~~~~~~~
When a workflow script is launched, Nextflow looks for a file named nextflow.config in the current directory and in the script base directory (if it is not the same as the current directory). Finally, it checks for the file: ``$HOME/.nextflow/config``.

When more than one of the above files exists, they are merged, so that the settings in the first override the same settings that may appear in the second, and so on.

The default config file search mechanism can be extended by providing an extra configuration file by using the command line option: ``-c <config file>``.



publishDir 
+++++++++++++
The publishDir directive allows you to publish the process output files to a specified folder. For example::

	process foo {
	    publishDir '/data/chunks'

	    output:
	    path 'chunk_*'

	    '''
	    printf 'Hola' | split -b 1 - chunk_
	    '''
	}

The above example splits the string Hola into file chunks of a single byte. When complete the ``chunk_*`` output files are published into the ``/data/chunks`` folder.

By default files are published to the target folder creating a symbolic link for each process output that links the file produced into the process working directory. This behavior can be modified using the mode parameter.

Multiple glob patterns
+++++++++++++++++++++++++
Multiple glob patterns can be specified using a list::

	Channel.fromFilePairs( ['/some/data/SRR*_{1,2}.fastq', '/other/data/QFF*_{1,2}.fastq'] )

Operators
~~~~~~~~~~~~~~~~~
collect
++++++++++++
The collect operator ``collects`` all the items emitted by a channel to a List and return the resulting object as a sole emission. For example::

	Channel
	    .of( 1, 2, 3, 4 )
	    .collect()
	    .view()

	# outputs
	[1,2,3,4]

Files and I/O
~~~~~~~~~~~~~~~~~~~
Opening files
+++++++++++++++
To access and work with files, use the file method, which returns a file system object given a file path string::

	myFile = file('some/path/to/my_file.file')

The ``file`` method can reference either files or directories, depending on what the string path refers to in the file system.

When using the wildcard characters *, ?, [] and {}, the argument is interpreted as a glob path matcher and the file method returns a list object holding the paths of files whose names match the specified pattern, or an empty list if no match is found::

	listOfFiles = file('some/path/*.fa')


