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
    and shared publicly on GitHub at <https://github.com/NERC-CEH/luct>

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
the form of web pages with the `workflowr` package, but not discussed further 
here.

## Project directory structure

### `_targets` directory

This directory is managed by the `targets` package. It contains the
metadata describing the status of the computational pipelines and the
cached results of those computations. You will normally only manipulate
these via functions from `targets`.

### `workflowr` directories

[`workflowr`](https://github.com/jdblischak/workflowr) creates a set of
standard directories. See the package documentation for details on how
these directories are used. The brief purposes are:

-   `analysis` - [`rmarkdown`](https://rmarkdown.rstudio.com/) analysis
    notebooks
-   `R` - R code not in analysis notebooks (changed from the `workflowr`
    default of `code`)
-   `data` - raw data and associated metadata
-   `docs` - automatically generated website
-   `output` - generated data and other objects

[`workflowr`](https://github.com/jdblischak/workflowr) only manages the
subset of files that it knows about, so you will need to manually stage
and commit any other files that need to be mirrored on GitHub.

If any files in `data` and `output` are more than trivially small, they
are not shared via Git and GitHub.

-   `.gitignore` is used to keep them out of Git.
-   There will be a separate mechanism (e.g.
    [Zenodo](https://about.zenodo.org/)) for sharing those large files.

### `manuscripts` directory

The `analysis` notebooks are for capturing all the analytical work that
was done, including exploratory work and abandoned directions. They
contain both the code and enough interpretation/explanation to make
sense of the results.

The notebooks will be too verbose, and inappropriately
structured/formatted for publication. Publishable documents are written
separately and kept in the `manuscripts` directory.

`manuscripts` contains a subdirectory for each
manuscript/document/presentation.

Each manuscript/document/presentation is prepared and formatted using a
package like [`rticles`](https://github.com/rstudio/rticles) or
[`bookdown`](https://github.com/rstudio/bookdown). Each document is
prepared in a separate subdirectory of `manuscripts` that contains all
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

## How to do things

-   All detailed setup instructions and notes go in this project-level
    `READ.md` file.
-   The `README.md` files in the subdirectories only state the purpose
    of each subdirectory and the files in that directory.

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

 
## Acknowledgements
The website is based on a template by [Ross W. Gayler][]
 
[Peter Levy]: https://github.com/peterlevy
[Ross W. Gayler]: https://www.rossgayler.com/

[this paper]: https://doi.org/10.5194/bg-15-1497-2018
