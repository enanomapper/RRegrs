RRegrs
======

eNanoMapper Developers |  [eNanoMapper Project] (http://www.enanomapper.net/)


The current tool is a collection of regression tools from R that could be used to search the best regression models for any dataset. The initial use of the script is aimed at finding QSAR models for chemoinformatics / nanotoxicology.

The full R script will contain: Loading dataset, Filter dataset, Scaling dataset, Feature selection, Regression models, Summary with top models, Statistics of the best model, etc.

The script will be modular in order to create flexible APIs.

The main authors are from the National Technical University of Athens (NTUA, Greece), Maastricht University (UM, Netherlands) and University of A Coruna (Spain).

Outputs:
- CSV files for statistics
- PDF files for plots

Releases
--------
* Version 0.0.5 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.32580.svg)](http://dx.doi.org/10.5281/zenodo.32580)
  * Bug in function ScalingDS: normalization
  * Incomplete option in function ScalingDS: standarization
  * Bug in plot 4, function SVMRFEreg
  * Bug in function RFRFEreg, removed abs() from resids calculation
  * Bug in function LASSOreg, removed abs() from resids calculation
* Version 0.0.4 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.21946.svg)](http://dx.doi.org/10.5281/zenodo.21946)
  * adds testing using the testthat package
  * removed the dependency on doMC so that it should install on Windows now (where doMC is not available)
  * added the tutorial PDF so that it is part of the distribution now
  * allows absolute paths now to point to the data
  * adds additional parameter value checking
* Version 0.0.3 [![DOI](https://zenodo.org/badge/6059/egonw/RRegrs.svg)](http://dx.doi.org/10.5281/zenodo.16446)
  * Completed the full author list
* Version 0.0.2: initial R package

Install
-------

    > install.packages(c("caret", "corrplot", "data.table")) # dependencies
    > install.packages("testthat") # if you want to test the package
    > install.packages("devtools") # to install from GitHub
    > library(devtools)
    > install_github("enanomapper/RRegrs", subdir="RRegrs")

Compile from source
-------------------

    $ R CMD build RRegrs
    $ R CMD check --as-cran RRegrs_0.0.4.tar.gz

Notes for Windows
-------------------

Before installing RRegrs:
* Update R to 3.2.2
* Install Rtools with the default configuration options from [https://cran.r-project.org/bin/windows/Rtools/] (https://cran.r-project.org/bin/windows/Rtools/). 

Caret version for RRegrs
----------------------------

The new version of caret is introducing some modifications that generates problems with some RRegrs methods.
Therefore, RRegrs is working well with caret version 6.0.52:

* please check caret version and report to this list

```> packageVersion('caret')```

* remove this version and the bug

```> remove.packages(c('caret'))```

* install Caret Package version 6.0.52, without the bug

```> library(devtools)```

```> install_version(package='caret',version='6.0-52')```


Help
-------------------

For questions, suggestions, problems or praise please use the Issues tab here or subscribe to [rregrs] (https://groups.google.com/forum/#!forum/rregrs) Google group.
