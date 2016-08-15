######################################################################################
# Makefile for building: CSPlibraries
# Usage:
# "make CPLEX" to produce CSPlibCPLEX
# "make XPRESS" to produce CSPlibXPRESS
# "make SCIP" to produce CSPlibSCIP
# "make" to produce all three libraries in one go
######################################################################################

####### Compiler, tools and options
CC            = gcc
CXX           = g++

RM = rm -f
MKDIR=mkdir

# Macros
CND_PLATFORM=MinGW-Windows
CND_DLIB_EXT=dll
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}
INCTAUPATH = -I../tauhitas/src

####### Object Files
OBJECTS     = ${OBJECTDIR}/src/Cspbridg.o ${OBJECTDIR}/src/Cspcard.o ${OBJECTDIR}/src/Cspmain.o ${OBJECTDIR}/src/Cspnet.o ${OBJECTDIR}/src/Jjsolver.o ${OBJECTDIR}/src/MT1RC.o ${OBJECTDIR}/src/My_time.o ${OBJECTDIR}/src/cspback.o ${OBJECTDIR}/src/cspbranc.o ${OBJECTDIR}/src/cspcapa.o ${OBJECTDIR}/src/cspcover.o ${OBJECTDIR}/src/cspdebug.o ${OBJECTDIR}/src/cspgomo.o ${OBJECTDIR}/src/cspheur.o ${OBJECTDIR}/src/cspprep.o ${OBJECTDIR}/src/cspprice.o ${OBJECTDIR}/src/cspsep.o ${OBJECTDIR}/src/cspsolve.o

CSPCPX = CSPlibCPLEX
CPXFLAGS = -ggdb -Wall -DCPLEX7 -DMICROSOFT2
CPXDIR = ../Cplex
CPXINC = -I$(CPXDIR)/include/ilcplex
CPXLIBS = -L$(CPXDIR)/lib -lcplex75

CSPXPR = CSPlibXPRESS
XPRFLAGS = -ggdb -Wall -DXPRESS_13 -DMICROSOFT2
XPRDIR = ../XPress
XPRINC = -I$(XPRDIR)
XPRLIBS	= -L$(XPRDIR) -lxprl -lxprs

CSPSCIP = CSPlibSCIP
SCIPFLAGS = -ggdb -Wall -DVSCIP -DMICROSOFT2
DIRLPS	    = ../scip-3.1.1
DIRSOPLEX   = ../soplex-2.0.1
SOPLEXLIB   = soplex-2.0.1.mingw.x86.gnu.opt
NLPILIB     = nlpi.cppad-3.1.1.mingw.x86.gnu.opt
SCIPLIB     = scip-3.1.1.mingw.x86.gnu.opt
OBJSCIPLIB  = objscip-3.1.1.mingw.x86.gnu.opt
LPISPXLIB   = lpispx-3.1.1.mingw.x86.gnu.opt
SCIPINC     = -I$(DIRLPS)/src -I$(DIRSOPLEX)/src
SCIPLIBS    = -L$(DIRLPS)/lib -L$(DIRSOPLEX)/lib -L$(DIRLPS)/lib -l$(OBJSCIPLIB) -l$(SCIPLIB) -l$(NLPILIB) -l$(LPISPXLIB) -l$(SOPLEXLIB)


all : CPLEX XPRESS SCIP

CPLEX :
	$(RM) -r ${CND_BUILDDIR}/${CND_CONF}
	$(RM) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/${CSPCPX}.${CND_DLIB_EXT}
	${MKDIR} -p ${OBJECTDIR}/src
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	$(CXX) -c $(CPXFLAGS) $(CPXINC) -o ${OBJECTDIR}/src/Cspbridg.o src/Cspbridg.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) -o ${OBJECTDIR}/src/Cspcard.o src/Cspcard.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/Cspmain.o src/Cspmain.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/Cspnet.o src/Cspnet.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/Jjsolver.o src/Jjsolver.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) -o ${OBJECTDIR}/src/MT1RC.o src/MT1RC.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) -o ${OBJECTDIR}/src/My_time.o src/My_time.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspback.o src/cspback.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) -o ${OBJECTDIR}/src/cspbranc.o src/cspbranc.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) -o ${OBJECTDIR}/src/cspcapa.o src/cspcapa.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) -o ${OBJECTDIR}/src/cspcover.o src/cspcover.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) -o ${OBJECTDIR}/src/cspdebug.o src/cspdebug.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspgomo.o src/cspgomo.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspheur.o src/cspheur.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) -o ${OBJECTDIR}/src/cspprep.o src/cspprep.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) -o ${OBJECTDIR}/src/cspprice.o src/cspprice.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) -o ${OBJECTDIR}/src/cspsep.o src/cspsep.c
	$(CXX) -c $(CPXFLAGS) $(CPXINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspsolve.o src/cspsolve.c
	$(CXX) -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/${CSPCPX}.${CND_DLIB_EXT} ${OBJECTS} $(CPXFLAGS) ${CPXLIBS} -shared

XPRESS :
	$(RM) -r ${CND_BUILDDIR}/${CND_CONF}
	$(RM) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/${CSPXPR}.${CND_DLIB_EXT}
	${MKDIR} -p ${OBJECTDIR}/src
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o ${OBJECTDIR}/src/Cspbridg.o src/Cspbridg.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o ${OBJECTDIR}/src/Cspcard.o src/Cspcard.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/Cspmain.o src/Cspmain.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/Cspnet.o src/Cspnet.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/Jjsolver.o src/Jjsolver.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o ${OBJECTDIR}/src/MT1RC.o src/MT1RC.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o ${OBJECTDIR}/src/My_time.o src/My_time.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspback.o src/cspback.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o ${OBJECTDIR}/src/cspbranc.o src/cspbranc.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o ${OBJECTDIR}/src/cspcapa.o src/cspcapa.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o ${OBJECTDIR}/src/cspcover.o src/cspcover.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o ${OBJECTDIR}/src/cspdebug.o src/cspdebug.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspgomo.o src/cspgomo.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspheur.o src/cspheur.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o ${OBJECTDIR}/src/cspprep.o src/cspprep.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o ${OBJECTDIR}/src/cspprice.o src/cspprice.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) -o ${OBJECTDIR}/src/cspsep.o src/cspsep.c
	$(CXX) -c $(XPRFLAGS) $(XPRINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspsolve.o src/cspsolve.c
	$(CXX) -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/${CSPXPR}.${CND_DLIB_EXT} ${OBJECTS} $(XPRFLAGS) ${XPRLIBS} -shared

SCIP : 
	$(RM) -r ${CND_BUILDDIR}/${CND_CONF}
	$(RM) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/${CSPSCIP}.${CND_DLIB_EXT}
	${MKDIR} -p ${OBJECTDIR}/src
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}	
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o ${OBJECTDIR}/src/Cspbridg.o src/Cspbridg.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o ${OBJECTDIR}/src/Cspcard.o src/Cspcard.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/Cspmain.o src/Cspmain.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/Cspnet.o src/Cspnet.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/Jjsolver.o src/Jjsolver.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o ${OBJECTDIR}/src/MT1RC.o src/MT1RC.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o ${OBJECTDIR}/src/My_time.o src/My_time.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspback.o src/cspback.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o ${OBJECTDIR}/src/cspbranc.o src/cspbranc.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o ${OBJECTDIR}/src/cspcapa.o src/cspcapa.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o ${OBJECTDIR}/src/cspcover.o src/cspcover.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o ${OBJECTDIR}/src/cspdebug.o src/cspdebug.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspgomo.o src/cspgomo.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspheur.o src/cspheur.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o ${OBJECTDIR}/src/cspprep.o src/cspprep.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o ${OBJECTDIR}/src/cspprice.o src/cspprice.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) -o ${OBJECTDIR}/src/cspsep.o src/cspsep.c
	$(CXX) -c $(SCIPFLAGS) $(SCIPINC) $(INCTAUPATH) -o ${OBJECTDIR}/src/cspsolve.o src/cspsolve.c
	$(CXX) -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/${CSPSCIP}.${CND_DLIB_EXT} ${OBJECTS} $(SCIPFLAGS) ${SCIPLIBS} -shared

clean:
	$(RM) -r ${CND_BUILDDIR}/${CND_CONF}
	$(RM) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/*.${CND_DLIB_EXT}
