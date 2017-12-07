
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include <GL/glut.h>
#include <CGAL/Timer.h>

#include "PPMC/configuration.h"
#include "PPMC/rangeCoder/rangecod.c"
#include "PPMC/rangeCoder/qsmodel.c"

#include "PPMC/frenetRotation.cpp"

#include "PPMC/mymeshComp.cpp"
#include "PPMC/mymeshCompTests.cpp"
#include "PPMC/mymeshDecomp.cpp"
#include "PPMC/mymeshAdaptiveQuantization.cpp"
#include "PPMC/mymeshLifting.cpp"
#include "PPMC/mymeshUtils.cpp"


#include "PPMC/mymeshBaseBuilder.h"
#include "PPMC/mymeshIO.cpp"
#include "PPMC/mymeshfunctions.hpp"
