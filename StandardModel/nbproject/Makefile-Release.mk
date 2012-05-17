#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/RGEvolutor.o \
	${OBJECTDIR}/src/EWSMThreeLoopEW2QCD.o \
	${OBJECTDIR}/src/StandardModel.o \
	${OBJECTDIR}/src/Particle.o \
	${OBJECTDIR}/src/EWSMTwoLoopQCD.o \
	${OBJECTDIR}/src/EWSM.o \
	${OBJECTDIR}/src/WilsonCoefficient.o \
	${OBJECTDIR}/src/EWSMTwoLoopEW.o \
	${OBJECTDIR}/src/EWSMcache.o \
	${OBJECTDIR}/src/QCD.o \
	${OBJECTDIR}/src/StandardModelMatching.o \
	${OBJECTDIR}/src/EWSMOneLoopEW.o \
	${OBJECTDIR}/src/CKM.o \
	${OBJECTDIR}/src/EWSMThreeLoopEW.o \
	${OBJECTDIR}/src/Meson.o \
	${OBJECTDIR}/src/EWSMThreeLoopQCD.o \
	${OBJECTDIR}/src/EWSMApproximateFormulae.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f2 \
	${TESTDIR}/TestFiles/f1

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libstandardmodel.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libstandardmodel.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libstandardmodel.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libstandardmodel.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libstandardmodel.a

${OBJECTDIR}/src/RGEvolutor.o: src/RGEvolutor.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/RGEvolutor.o src/RGEvolutor.cpp

${OBJECTDIR}/src/EWSMThreeLoopEW2QCD.o: src/EWSMThreeLoopEW2QCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMThreeLoopEW2QCD.o src/EWSMThreeLoopEW2QCD.cpp

${OBJECTDIR}/src/StandardModel.o: src/StandardModel.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/StandardModel.o src/StandardModel.cpp

${OBJECTDIR}/src/Particle.o: src/Particle.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Particle.o src/Particle.cpp

${OBJECTDIR}/src/EWSMTwoLoopQCD.o: src/EWSMTwoLoopQCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMTwoLoopQCD.o src/EWSMTwoLoopQCD.cpp

${OBJECTDIR}/src/EWSM.o: src/EWSM.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSM.o src/EWSM.cpp

${OBJECTDIR}/src/WilsonCoefficient.o: src/WilsonCoefficient.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/WilsonCoefficient.o src/WilsonCoefficient.cpp

${OBJECTDIR}/src/EWSMTwoLoopEW.o: src/EWSMTwoLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMTwoLoopEW.o src/EWSMTwoLoopEW.cpp

${OBJECTDIR}/src/EWSMcache.o: src/EWSMcache.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMcache.o src/EWSMcache.cpp

${OBJECTDIR}/src/QCD.o: src/QCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/QCD.o src/QCD.cpp

${OBJECTDIR}/src/StandardModelMatching.o: src/StandardModelMatching.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/StandardModelMatching.o src/StandardModelMatching.cpp

${OBJECTDIR}/src/EWSMOneLoopEW.o: src/EWSMOneLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMOneLoopEW.o src/EWSMOneLoopEW.cpp

${OBJECTDIR}/src/CKM.o: src/CKM.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/CKM.o src/CKM.cpp

${OBJECTDIR}/src/EWSMThreeLoopEW.o: src/EWSMThreeLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMThreeLoopEW.o src/EWSMThreeLoopEW.cpp

${OBJECTDIR}/src/Meson.o: src/Meson.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Meson.o src/Meson.cpp

${OBJECTDIR}/src/EWSMThreeLoopQCD.o: src/EWSMThreeLoopQCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMThreeLoopQCD.o src/EWSMThreeLoopQCD.cpp

${OBJECTDIR}/src/EWSMApproximateFormulae.o: src/EWSMApproximateFormulae.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMApproximateFormulae.o src/EWSMApproximateFormulae.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f2: ${TESTDIR}/tests/RunningMass.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f2 $^ ${LDLIBSOPTIONS} 

${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/StandardModelTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS} 


${TESTDIR}/tests/RunningMass.o: tests/RunningMass.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${TESTDIR}/tests/RunningMass.o tests/RunningMass.cpp


${TESTDIR}/tests/StandardModelTest.o: tests/StandardModelTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${TESTDIR}/tests/StandardModelTest.o tests/StandardModelTest.cpp


${OBJECTDIR}/src/RGEvolutor_nomain.o: ${OBJECTDIR}/src/RGEvolutor.o src/RGEvolutor.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/RGEvolutor.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/RGEvolutor_nomain.o src/RGEvolutor.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/RGEvolutor.o ${OBJECTDIR}/src/RGEvolutor_nomain.o;\
	fi

${OBJECTDIR}/src/EWSMThreeLoopEW2QCD_nomain.o: ${OBJECTDIR}/src/EWSMThreeLoopEW2QCD.o src/EWSMThreeLoopEW2QCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWSMThreeLoopEW2QCD.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMThreeLoopEW2QCD_nomain.o src/EWSMThreeLoopEW2QCD.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWSMThreeLoopEW2QCD.o ${OBJECTDIR}/src/EWSMThreeLoopEW2QCD_nomain.o;\
	fi

${OBJECTDIR}/src/StandardModel_nomain.o: ${OBJECTDIR}/src/StandardModel.o src/StandardModel.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/StandardModel.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/StandardModel_nomain.o src/StandardModel.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/StandardModel.o ${OBJECTDIR}/src/StandardModel_nomain.o;\
	fi

${OBJECTDIR}/src/Particle_nomain.o: ${OBJECTDIR}/src/Particle.o src/Particle.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/Particle.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Particle_nomain.o src/Particle.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/Particle.o ${OBJECTDIR}/src/Particle_nomain.o;\
	fi

${OBJECTDIR}/src/EWSMTwoLoopQCD_nomain.o: ${OBJECTDIR}/src/EWSMTwoLoopQCD.o src/EWSMTwoLoopQCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWSMTwoLoopQCD.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMTwoLoopQCD_nomain.o src/EWSMTwoLoopQCD.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWSMTwoLoopQCD.o ${OBJECTDIR}/src/EWSMTwoLoopQCD_nomain.o;\
	fi

${OBJECTDIR}/src/EWSM_nomain.o: ${OBJECTDIR}/src/EWSM.o src/EWSM.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWSM.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSM_nomain.o src/EWSM.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWSM.o ${OBJECTDIR}/src/EWSM_nomain.o;\
	fi

${OBJECTDIR}/src/WilsonCoefficient_nomain.o: ${OBJECTDIR}/src/WilsonCoefficient.o src/WilsonCoefficient.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/WilsonCoefficient.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/WilsonCoefficient_nomain.o src/WilsonCoefficient.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/WilsonCoefficient.o ${OBJECTDIR}/src/WilsonCoefficient_nomain.o;\
	fi

${OBJECTDIR}/src/EWSMTwoLoopEW_nomain.o: ${OBJECTDIR}/src/EWSMTwoLoopEW.o src/EWSMTwoLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWSMTwoLoopEW.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMTwoLoopEW_nomain.o src/EWSMTwoLoopEW.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWSMTwoLoopEW.o ${OBJECTDIR}/src/EWSMTwoLoopEW_nomain.o;\
	fi

${OBJECTDIR}/src/EWSMcache_nomain.o: ${OBJECTDIR}/src/EWSMcache.o src/EWSMcache.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWSMcache.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMcache_nomain.o src/EWSMcache.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWSMcache.o ${OBJECTDIR}/src/EWSMcache_nomain.o;\
	fi

${OBJECTDIR}/src/QCD_nomain.o: ${OBJECTDIR}/src/QCD.o src/QCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/QCD.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/QCD_nomain.o src/QCD.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/QCD.o ${OBJECTDIR}/src/QCD_nomain.o;\
	fi

${OBJECTDIR}/src/StandardModelMatching_nomain.o: ${OBJECTDIR}/src/StandardModelMatching.o src/StandardModelMatching.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/StandardModelMatching.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/StandardModelMatching_nomain.o src/StandardModelMatching.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/StandardModelMatching.o ${OBJECTDIR}/src/StandardModelMatching_nomain.o;\
	fi

${OBJECTDIR}/src/EWSMOneLoopEW_nomain.o: ${OBJECTDIR}/src/EWSMOneLoopEW.o src/EWSMOneLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWSMOneLoopEW.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMOneLoopEW_nomain.o src/EWSMOneLoopEW.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWSMOneLoopEW.o ${OBJECTDIR}/src/EWSMOneLoopEW_nomain.o;\
	fi

${OBJECTDIR}/src/CKM_nomain.o: ${OBJECTDIR}/src/CKM.o src/CKM.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/CKM.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/CKM_nomain.o src/CKM.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/CKM.o ${OBJECTDIR}/src/CKM_nomain.o;\
	fi

${OBJECTDIR}/src/EWSMThreeLoopEW_nomain.o: ${OBJECTDIR}/src/EWSMThreeLoopEW.o src/EWSMThreeLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWSMThreeLoopEW.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMThreeLoopEW_nomain.o src/EWSMThreeLoopEW.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWSMThreeLoopEW.o ${OBJECTDIR}/src/EWSMThreeLoopEW_nomain.o;\
	fi

${OBJECTDIR}/src/Meson_nomain.o: ${OBJECTDIR}/src/Meson.o src/Meson.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/Meson.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Meson_nomain.o src/Meson.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/Meson.o ${OBJECTDIR}/src/Meson_nomain.o;\
	fi

${OBJECTDIR}/src/EWSMThreeLoopQCD_nomain.o: ${OBJECTDIR}/src/EWSMThreeLoopQCD.o src/EWSMThreeLoopQCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWSMThreeLoopQCD.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMThreeLoopQCD_nomain.o src/EWSMThreeLoopQCD.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWSMThreeLoopQCD.o ${OBJECTDIR}/src/EWSMThreeLoopQCD_nomain.o;\
	fi

${OBJECTDIR}/src/EWSMApproximateFormulae_nomain.o: ${OBJECTDIR}/src/EWSMApproximateFormulae.o src/EWSMApproximateFormulae.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWSMApproximateFormulae.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -I. -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMApproximateFormulae_nomain.o src/EWSMApproximateFormulae.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWSMApproximateFormulae.o ${OBJECTDIR}/src/EWSMApproximateFormulae_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f2 || true; \
	    ${TESTDIR}/TestFiles/f1 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libstandardmodel.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc