.. _backbone-label:

Nextflow
==============================

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



