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
CND_PLATFORM=MinGW-Windows
CND_DLIB_EXT=dll
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/104313830/Cspbridg.o \
	${OBJECTDIR}/_ext/104313830/Cspcard.o \
	${OBJECTDIR}/_ext/104313830/Cspmain.o \
	${OBJECTDIR}/_ext/104313830/Cspnet.o \
	${OBJECTDIR}/_ext/104313830/Jjsolver.o \
	${OBJECTDIR}/_ext/104313830/MT1RC.o \
	${OBJECTDIR}/_ext/104313830/My_time.o \
	${OBJECTDIR}/_ext/104313830/cspback.o \
	${OBJECTDIR}/_ext/104313830/cspbranc.o \
	${OBJECTDIR}/_ext/104313830/cspcapa.o \
	${OBJECTDIR}/_ext/104313830/cspcover.o \
	${OBJECTDIR}/_ext/104313830/cspdebug.o \
	${OBJECTDIR}/_ext/104313830/cspgomo.o \
	${OBJECTDIR}/_ext/104313830/cspheur.o \
	${OBJECTDIR}/_ext/104313830/cspprep.o \
	${OBJECTDIR}/_ext/104313830/cspprice.o \
	${OBJECTDIR}/_ext/104313830/cspsep.o \
	${OBJECTDIR}/_ext/104313830/cspsolve.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCSP.${CND_DLIB_EXT}

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCSP.${CND_DLIB_EXT}: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCSP.${CND_DLIB_EXT} ${OBJECTFILES} ${LDLIBSOPTIONS} -shared

${OBJECTDIR}/_ext/104313830/Cspbridg.o: ../CSP/src/Cspbridg.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/Cspbridg.o ../CSP/src/Cspbridg.c

${OBJECTDIR}/_ext/104313830/Cspcard.o: ../CSP/src/Cspcard.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/Cspcard.o ../CSP/src/Cspcard.c

${OBJECTDIR}/_ext/104313830/Cspmain.o: ../CSP/src/Cspmain.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/Cspmain.o ../CSP/src/Cspmain.c

${OBJECTDIR}/_ext/104313830/Cspnet.o: ../CSP/src/Cspnet.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/Cspnet.o ../CSP/src/Cspnet.c

${OBJECTDIR}/_ext/104313830/Jjsolver.o: ../CSP/src/Jjsolver.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/Jjsolver.o ../CSP/src/Jjsolver.c

${OBJECTDIR}/_ext/104313830/MT1RC.o: ../CSP/src/MT1RC.C 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.cc) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/MT1RC.o ../CSP/src/MT1RC.C

${OBJECTDIR}/_ext/104313830/My_time.o: ../CSP/src/My_time.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/My_time.o ../CSP/src/My_time.c

${OBJECTDIR}/_ext/104313830/cspback.o: ../CSP/src/cspback.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/cspback.o ../CSP/src/cspback.c

${OBJECTDIR}/_ext/104313830/cspbranc.o: ../CSP/src/cspbranc.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/cspbranc.o ../CSP/src/cspbranc.c

${OBJECTDIR}/_ext/104313830/cspcapa.o: ../CSP/src/cspcapa.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/cspcapa.o ../CSP/src/cspcapa.c

${OBJECTDIR}/_ext/104313830/cspcover.o: ../CSP/src/cspcover.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/cspcover.o ../CSP/src/cspcover.c

${OBJECTDIR}/_ext/104313830/cspdebug.o: ../CSP/src/cspdebug.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/cspdebug.o ../CSP/src/cspdebug.c

${OBJECTDIR}/_ext/104313830/cspgomo.o: ../CSP/src/cspgomo.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/cspgomo.o ../CSP/src/cspgomo.c

${OBJECTDIR}/_ext/104313830/cspheur.o: ../CSP/src/cspheur.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/cspheur.o ../CSP/src/cspheur.c

${OBJECTDIR}/_ext/104313830/cspprep.o: ../CSP/src/cspprep.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/cspprep.o ../CSP/src/cspprep.c

${OBJECTDIR}/_ext/104313830/cspprice.o: ../CSP/src/cspprice.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/cspprice.o ../CSP/src/cspprice.c

${OBJECTDIR}/_ext/104313830/cspsep.o: ../CSP/src/cspsep.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/cspsep.o ../CSP/src/cspsep.c

${OBJECTDIR}/_ext/104313830/cspsolve.o: ../CSP/src/cspsolve.c 
	${MKDIR} -p ${OBJECTDIR}/_ext/104313830
	${RM} $@.d
	$(COMPILE.c) -O2  -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/104313830/cspsolve.o ../CSP/src/cspsolve.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libCSP.${CND_DLIB_EXT}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
