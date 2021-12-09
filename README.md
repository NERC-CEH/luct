# Tracking Land-Use Change

This is an open, reproducible, computational research project on land-use change in the UK.

The code is developed by:

* [Peter Levy][]

based on a combination of targets and workflowr packages in R.

The background to the method for the estimation of land-use change is described in [this paper][].

## Project dependencies

This is an open, shareable, reproducible, computational research
project.

-   All the computational work and document preparation is done with the
    [R](https://www.r-project.org/) statistical computing environment.

-   The research project is contained in a single directory, with the 
    exception that some data sets are too large to store on GitHub.

-   We use the [`renv`](https://rstudio.github.io/renv/) package to
    manage the R package versions used by the project

-   We are using the [`targets`](https://github.com/ropensci/targets)
    package to structure the project so that the work is computationally
    reproducible.

-   The project code and documents are shared publicly on GitHub at
    <https://github.com/NERC-CEH/luct>

-   The main report is produced using [`bookdown`](https://bookdown.org/) 
    and shared publicly on GitHub at <https://nerc-ceh.github.io/luct/>

-   We are exploring the
    [`workflowr`](https://github.com/jdblischak/workflowr) package to
    structure the project so that all the materials and outputs are
    available via an openly accessible, automatically generated website. 
    However, GitHub cannot currently show both the `bookdown` and `workflowr` 
    website documents simultaneously, so this is still under investigation.


## Workflow management

The project uses the R [`targets`](https://github.com/ropensci/targets)
package to structure and manage the workflow and to make it reproducible.
Central to this is the idea of the workflow as a "pipeline" - a defined
list of functions which transform data. 
Here, the core pipeline contains the computational steps that read, reformat 
and process the input data (time series and maps of land use and land-use 
change data), and run the data assimilation steps that estimate the matrices 
of land-use change, and produce the maps of past land use.
Potentially there can be multiple pipelines, which produce other analyses, 
reports, or publications, in addition to the core process. 
These are used to generate documentation in
the form of web pages with the `workflowr` package, but are not discussed further 
here.

The core pipeline is defined in the file `_targets.R` as a list of "targets".  
The targets represent the steps in the series of computations which make up the 
pipeline.
A target is defined with the syntax `tar_target(target_name, function_name(inputs))`.
The target is thus a named R data object which is the outcome of a named function
with specified inputs. The one exception to this is that the target may simply be 
a file for input or output.
In the current project, the core pipeline is a list of 87 targets which specify the 
input files, the reformatting and transformation of these data, and subsequent
calculations which make up the data assimilation algorithm.

The pipeline is managed using a ["Make"-like](https://en.wikipedia.org/wiki/Make_(software)) 
procedure, which analyses the dependencies between the different steps in the pipeline.
If there have been no changes to the code in the target functions or input data since 
the last time it was run, it identifies that everything is up-to-date, and no further 
computation is needed. If any the source code of target function or the content of any 
data file has changed, it identifies which parts of the pipeline are affected by this,
and all the dependencies are recomputed. This has several advantages: forcing the workflow 
to be declared at a higher level of abstraction; only running the necessary computation, 
so saving run-time for tasks that are already up to date; and most importantly, 
providing tangible evidence that the results match the underlying code and data, 
and confirm the computation is reproducible. So as to identify changes, each target is 
represented by its 
[hash value](https://en.wikipedia.org/wiki/Hash_function), stored in the 
`_targets` directory. 

## Project directory structure

### `_targets` directory

This directory is managed by the `targets` package. It contains the
metadata describing the status of the computational pipelines and the
cached results of those computations.

### `analysis` directory

[`workflowr`](https://github.com/jdblischak/workflowr) creates a set of
standard directories. See the package documentation for details on how
these directories are used. The `analysis` contains 
[`rmarkdown`](https://rmarkdown.rstudio.com/) notebooks which document
the workflow. These are still in development.

### `R` directory
This contains the bulk of the R source code for the functions used in the project.

### `data-raw` directory
This contains the raw data files for the project, in their original form 
as far as possible.
To avoid duplication, this is a symbolic link to an 
[earlier iteration](https://github.com/NERC-CEH/luc_track/tree/master/data-raw)
However, many of these are too large to share via GitHub, and would need to be 
shared by another mechanism (e.g. as binary assets). 

### `data` directory
This contains the processed data files resulting from transformations of the raw 
data. This typically involves reprojection, reclassification, filtering and 
unit conversions. Again, many of these are too large to share via GitHub.

### `docs` directory
This contains the html web pages generated by the Rmarkdown files in 
with `workflowr` or `bookdown`.

### `output` directory
This contains output files from the project, the results of the data assimlation.

### `slurm` directory
This contains files for the steps which require high-performance computing,
run via [slurm](https://en.wikipedia.org/wiki/Slurm_Workload_Manager), the 
widely used job scheduling system on HPC systems. These are generic enough 
to run on any HPC machine with slurm, and have been run on both JASMIN and 
POLAR, althouh the queue names, number of processors and memory limits
will be system-specific.

### `manuscripts` directory

The report is prepared and formatted using 
[`bookdown`](https://github.com/rstudio/bookdown) in a subdirectory of 
`manuscripts` that contains all
the necessary infrastructure files (templates, bibliographies, etc.).

### `renv` directory

The [`renv`](https://rstudio.github.io/renv/) package keeps track of the
R packages (and their versions) used by the project. It allows anyone to
reinstate the same packages and versions in their local copy of the
project.

The `renv` directory contains the information need by `renv` to
reinstate the local package environment

### `.gitignore`

`.gitignore` in the R project root directory is used for all manual
entries so that all the manual rules are in one place. Packages, such as
`renv`, may create their own `.gitignore` files in subdirectories that
they manage.


### Installation

Assuming you already have a current version of R 
installed, clone the project repository <https://github.com/NERC-CEH/luct>
    from GitHub.

When you open the project, you may get warning messages about packages
not being installed. This is because you need to use the
[`renv`](https://rstudio.github.io/renv/) package to reinstate the
packages that are used by the project.

1.  Install [`renv`](https://rstudio.github.io/renv/) in that project if
    it is not already installed

2.  Use `renv::restore()` to install all the needed packages in the
    project-specific library:

        renv::restore()

### Get data

Any files in `data`, `output` and `_targets` that are more than
trivially small are not shared via Git and GitHub. They will be shared
via a separate, yet to be determined, mechanism (e.g.
[Zenodo](https://about.zenodo.org/)).


### `renv` collaboration

The `renv` package is used to keep track of the installed packages and
their versions. See the `renv` [collaboration
guide](https://rstudio.github.io/renv/articles/collaborating.html) or
the workflow for synchronising package environments between
collaborators.

## Still to do

-   More detailed setup instructions and notes should go in this project-level
    `READ.md` file.
-   The `README.md` files in the subdirectories are currently generic, but should
    describe the purpose of each subdirectory and the files in that directory.
 
## Acknowledgements
The website is based on a template by [Ross W. Gayler][]
 
[Peter Levy]: https://github.com/peterlevy
[Ross W. Gayler]: https://www.rossgayler.com/

[this paper]: https://doi.org/10.5194/bg-15-1497-2018
