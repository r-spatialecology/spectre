TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

HEADERS += \
        ../../src/minconf.h \
      #  ../../src/test-minconf.h \
        ../../src/optimizer.h

SOURCES += \
        ../../src/minconf.cpp \
      #  ../../src/test-minconf.cpp \
        ../../src/optimizer.cpp \
        main.cpp

# OpenMP support
QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

## Template from the example at http://dirk.eddelbuettel.com/blog/2011/03/25/#rinside_and_qt

## comment this out if you need a different version of R,
## and set set R_HOME accordingly as an environment variable
R_HOME = $$system(R RHOME)

## include headers and libraries for R
RCPPFLAGS = 		$$system($$R_HOME/bin/R CMD config --cppflags)
RLDFLAGS = 		$$system($$R_HOME/bin/R CMD config --ldflags)
RBLAS = 		$$system($$R_HOME/bin/R CMD config BLAS_LIBS)
RLAPACK = 		$$system($$R_HOME/bin/R CMD config LAPACK_LIBS)

## if you need to set an rpath to R itself, also uncomment
RRPATH =		-Wl,-rpath,$$R_HOME/lib

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL = 		$$system($$R_HOME/bin/Rscript -e \"Rcpp:::CxxFlags\(\)\")
RCPPLIBS = 		$$system($$R_HOME/bin/Rscript -e \"Rcpp:::LdFlags\(\)\")

## include headers and libraries for RInside embedding classes
RINSIDEINCL = 		$$system($$R_HOME/bin/Rscript -e \"RInside:::CxxFlags\(\)\")
RINSIDELIBS = 		$$system($$R_HOME/bin/Rscript -e \"RInside:::LdFlags\(\)\")

## include headers RcppProgress
RCPPPROGRESSINCL = $$system($$R_HOME/bin/Rscript -e \"RcppProgress:::CxxFlags\(\)\")

## compiler etc settings used in default make rules
QMAKE_CXXFLAGS +=	$$RCPPWARNING $$RCPPFLAGS $$RCPPINCL $$RINSIDEINCL $$RCPPPROGRESSINCL
QMAKE_LIBS += $$RLDFLAGS $$RBLAS $$RLAPACK $$RINSIDELIBS $$RCPPLIBS


