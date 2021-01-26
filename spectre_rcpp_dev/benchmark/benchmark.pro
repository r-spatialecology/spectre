TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

HEADERS += \
    ../../tests/benchmark/catch.hpp \
    ../../src/minconf.h

SOURCES += \
    ../../tests/benchmark/benchmark-minconf.cpp \
    ../../src/minconf.cpp

## Template from the example at http://dirk.eddelbuettel.com/blog/2011/03/25/#rinside_and_qt

## comment this out if you need a different version of R,
## and set set R_HOME accordingly as an environment variable
R_HOME = $$system(R RHOME)

## include headers and libraries for R
RCPPFLAGS = 		$$system($$R_HOME/bin/R CMD config --cppflags)
RLDFLAGS = 		$$system($$R_HOME/bin/R CMD config --ldflags)

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL = 		$$system($$R_HOME/bin/Rscript -e \"Rcpp:::CxxFlags\(\)\")
RCPPLIBS = 		$$system($$R_HOME/bin/Rscript -e \"Rcpp:::LdFlags\(\)\")

## include headers RcppProgress
RCPPPROGRESSINCL = $$system($$R_HOME/bin/Rscript -e \"RcppProgress:::CxxFlags\(\)\")

## compiler etc settings used in default make rules
QMAKE_CXXFLAGS += $$RCPPFLAGS $$RCPPINCL $$RCPPPROGRESSINCL
QMAKE_LIBS += $$RLDFLAGS $$RCPPLIBS
