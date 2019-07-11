TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

HEADERS += \
        ../src/calculate_solution_commonness.h \
        ../src/mh_optimizer.h \
        ../src/rcpp_sample.h

SOURCES += \
        main.cpp \
        ../src/calculate_solution_commonness.cpp \
        ../src/mh_optimizer.cpp \
        ../src/rcpp_sample.cpp

# OpenMP support
QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

# Adding R, RInside & Rcpp
INCLUDEPATH += /usr/share/R/include
INCLUDEPATH += /usr/lib/R/site-library/Rcpp/include
INCLUDEPATH += /usr/lib/R/site-library/RInside/include/
LIBS += -L/usr/lib -lR
LIBS += -L/usr/lib/R/site-library/RInside/lib -lRInside


