.. _backbone-label:

Nextflow-Metrics and reports
==============================
Nextflow is able to produce multiple reports and charts providing several runtime metrics and execution information.

The ``-with-report`` option enables the creation of the workflow execution report.

The ``-with-trace`` option enables the create of a tab separated file containing runtime information for each executed task, including: submission time, start time, completion time, cpu and memory used..

The ``-with-timeline`` option enables the creation of the workflow timeline report showing how processes where executed along time. This may be useful to identify most time consuming tasks and bottlenecks. See an example at this link.

The ``-with-dag`` option enables to rendering of the workflow execution direct acyclic graph representation. Note: this feature requires the installation of Graphviz, an open source graph visualization software, in your system.

