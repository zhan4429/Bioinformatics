.. _backbone-label:

find and xargs
==============================

find
~~~~~~~~~~
``find`` Finds files matching an expression. The ``find`` command is similar to ls but in many ways it is more powerful. It can be used to recursively search the directory tree for a specified path name, seeking files that match a given Boolean expression (a test which returns true or false.

``find . -name "*.embl"`` 
This command will return the files which name has the ``.embl`` suffix.

``find . -type d``. 
This command will return all the subdirectories contained in the current directory.

``-mtime`` search files by modifying date ``-atime`` search files by last access date ``-size`` search files by file size ``-user`` search files by user they belong to.


xargs
~~~~~~~~~~
The ``xargs`` command in UNIX is a command line utility for building an execution pipeline from standard input. Whilst tools like ``grep`` can accept standard input as a parameter, many other tools cannot. Using ``xargs`` allows tools like ``echo``, ``rm`` and ``mkdir`` to accept standard input as arguments.

How to use xargs
~~~~~~~~~~~~~~~~~
By default ``xargs`` reads items from standard input as separated by blanks and executes a command once for each argument. In the following example standard input is piped to ``xargs`` and the ``mkdir`` command is run for each argument, creating three folders::

    echo 'one two three' | xargs mkdir
    ls
    one two three

When filenames contains spaces you need to use ``-d`` option to change delimiter::

    find . -name '*.txt' | xargs -d '\n' rm


How to use xargs with find
~~~~~~~~~~~~~~~~~~~~~~~~~~
The most common usage of ``xargs`` is to use it with the ``find`` command. This uses ``find`` to search for files or directories and then uses xargs to operate on the results. Typical examples of this are changing the ownership of files or moving files.

``find`` and ``xargs`` can be used together to operate on files that match certain attributes. In the following example files older than two weeks in the temp folder are found and then piped to the xargs command which runs the rm command on each file and removes them::

    find /tmp -mtime +14 | xargs rm

How to print commands that are executed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``-t`` option prints each command that will be executed to the terminal. This can be helpful when debugging scripts::

    echo 'one two three' | xargs -t rm
        rm one two three

How to view the command and prompt for execution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``-p`` command will print the command to be executed and prompt the user to run it. This can be useful for destructive operations where you really want to be sure on the command to be run::

    echo 'one two three' | xargs -p touch
    touch one two three ?...