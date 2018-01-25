TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    simulation.cpp \
    particle.cpp \
    potential.cpp \
    trialfunction.cpp \
    system.cpp \
    Vec3/vec3.cpp \
    Potentials/simpleharmonicoscillator.cpp

HEADERS += \
    simulation.h \
    particle.h \
    potential.h \
    trialfunction.h \
    system.h \
    Vec3/vec3.h \
    Potentials/simpleharmonicoscillator.h
