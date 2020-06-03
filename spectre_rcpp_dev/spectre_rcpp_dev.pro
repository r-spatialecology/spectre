TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

HEADERS += \
        ../src/backtracking.h \
        ../src/calculate_solution_commonness.h \
        ../src/constraint_satisfaction_problem.h \
        ../src/mh_optimizer.h \
        ../src/minconf.h \
        ../src/optimizer.h \
        ../src/rcpp_sample.h

SOURCES += \
        ../src/backtracking.cpp \
        ../src/constraint_satisfaction_problem.cpp \
        ../src/mh_optimizer.cpp \
        ../src/minconf.cpp \
        ../src/optimizer.cpp \
        main.cpp \
        ../src/calculate_solution_commonness.cpp \
        ../src/rcpp_sample.cpp

# OpenMP support
QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

# Adding R, RInside & Rcpp
INCLUDEPATH += /usr/share/R/include
INCLUDEPATH += /usr/lib/R/site-library/Rcpp/include
INCLUDEPATH += /usr/lib/R/site-library/RInside/include
LIBS += -L/usr/lib -lR
LIBS += -L/usr/lib/R/site-library/RInside/lib -lRInside


