--
-- nextflow 21.10.6 modulefile
--
-- "URL: https://www.psc.edu/resources/software"
-- "Category: Other"
-- "Description: Nextflow enables scalable and reproducible scientific workflows using software containers."

whatis("Name: nextflow")
whatis("Version: 21.10.6")
whatis("Category: Other")
whatis("Description: Nextflow enables scalable and reproducible scientific workflows using software containers.")

help([[
Nextflow enables scalable and reproducible scientific workflows using software containers.

To load the module, type

> module load nextflow/21.10.6

To unload the module, type

> module unload nextflow/21.10.6

Documentation
-------------
For help, type

> nextflow -h

Tools included in this module are

* nextflow
]])

local package = "nextflow"
local version = "21.10.6"
local base    = pathJoin("/opt/packages",package,version)
prepend_path("PATH", base)