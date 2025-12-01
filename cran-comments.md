# spectre 1.0.5
Updated the DESCRIPTION and CITATION files to reflect new standards and removed broken links in README.md

# spectre 1.0.4
Removed C++ version requirement

# spectre 1.0.3

We added and updated the citation information of the package.

# spectre 1.0.2 

## Resubmission and update
This is a resubmission of the spectre package.

# LTO Issue:
We have manually rewritten RcppExports.cpp to account for the Catch2 framework.

# spectre 1.0.2 

## Resubmission and update
This is a resubmission of the spectre package in an updated form. Below is a list of changes made to address each issues from the previous submission:

Issue 1: "memtest"
The problem occurred because at one point in the code only the upper triangular matrix was calculated, not the complete matrix. Since the diagonal and the lower triangular matrix have values of "NA", the unsigned integer overflowed.
The problem was solved in two places
- the lower triangular matrix and diagonal are now set in the constructor of the C++ class "NA" instead of in the R function.
- the loops in the C++ code iterate only over the upper triangular matrix instead of iterating over the complete matrix and skipping "NA".

We also added a new unit test for this issue.

# Issue 2: LTO
The problem occurred because we use the Catch2 framework with testthat to test the C++ code (https://github.com/r-lib/testthat/issues/1230). This was fixed as described in the issue by manually adjusting the RcppExports.cpp.

Beyond addressing these issues we have also removed the `gdm` package dependency by creating a .RDS file of data to be used in tests and vignettes.

We have retested the package as follows:

## Test environments
* Local
  * macOS BigSur, R 4.0.5
  * Windows 10 Home, R 4.0.3
  * Ubuntu 20.04 LTS, R 4.0.5
* GitHub Actions 
  * windows-latest, R: 'release'
  * macOS-latest, R: 'release'
  * ubuntu-20.04, R: 'release'
* win-builder (release, devel)

## R CMD check results

0 errors | 0 warnings | 0 note

* This is an updated resubmission.

## Reverse dependencies
There are currently no reverse dependencies.
---

# spectre 1.0.1 

## Resubmission
This is a second resubmission. Below is a list of changes made to address each comment:

   Found the following (possibly) invalid URLs:
     URL: https://cran.rstudio.com/web/packages/spectre/index.html
       From: README.md
       Status: 200
       Message: OK
       CRAN URL not in canonical form

Please fix and resubmit.


* The link was updated, to the https canonical form 

---

# spectre 1.0.1 

## Resubmission
This is a resubmission. Below is a list of changes made to address each comment:

Found the following (possibly) invalid URLs:
URL: http://cran.rstudio.com/web/packages/spectre/index.html (moved to https://cran.rstudio.com/web/packages/spectre/index.html)
From: README.md
Status: 200
Message: OK
CRAN URL not in canonical form

Canonical: https://CRAN.R-project.org/package=spectre

Please fix and resubmit.


* The link was updated

---

Update Version: 1.0.1

Update addresses 2 notes and 1 error discovered by CRAN package check run 2021-05-30 23:49:52 CEST. Solutions to points brought up are listed below:

Check: installed package size
Result: NOTE
     installed size is 6.3Mb
     sub-directories of 1Mb or more:
     libs 5.6Mb
     
* We further cleaned the package to remove any garbage. Unfortunately, the package is still this note but we believe that any further removals from the package will hamper usability.

Check: dependencies in R code
Result: NOTE
    Namespace in Imports field not imported from: ‘RcppProgress’
     All declared Imports should be used.
     
* We moved ‘RcppProgress’ from "Imports" to "Suggests" within the DESCRIPTION.

Check: whether package can be installed
Result: ERROR
    Installation failed.
Flavor: r-patched-solaris-x86

* This error is produced due to problems between Solaris and testthat (r-lib/testthat#1257). To prevent this error we have exluded the unit tests via compiler macros.

---

# spectre v1.0.0

## Resubmission
This is a resubmission. Below is a list of changes made to address each comment:

Found the following (possibly) invalid URLs: URL: https://www.tidyverse.org/lifecycle/#stable (moved to https://lifecycle.r-lib.org/articles/stages.html)
From: README.md
Status: 200
Message: OK

Please change http --> https, add trailing slashes, or follow moved content as appropriate.

* The link was updated

Is there some reference about the method you can add in the Description field in the form Authors (year) <doi:.....>?

* We added a reference about the method in the Description field

## Test environments
* Local
  * macOS BigSur, R 4.0.5
  * Windows 10 Home, R 4.0.3
  * Ubuntu 20.04 LTS, R 4.0.5
* GitHub Actions 
  * windows-latest, R: 'release'
  * macOS-latest, R: 'release'
  * ubuntu-20.04, R: 'release'
* win-builder (release, devel)

## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release.

## Reverse dependencies
There are currently no reverse dependencies.