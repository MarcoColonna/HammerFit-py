#
# Makefile to compile the Project
#

SHELL = /bin/sh

# ROOT stuff
ROOTCFLAGS    = $(shell /opt/miniconda/envs/root_forge/bin/root-config --cflags)
ROOTLIBS      = $(shell /opt/miniconda/envs/root_forge/bin/root-config --libs)
ROOTGLIBS     = $(shell /opt/miniconda/envs/root_forge/bin/root-config --glibs)
ROOTLDFLAGS   = $(shell /opt/miniconda/envs/root_forge/bin/root-config --ldflags)

EVTGENLIBS    = -L/home/mcolonna/evtgen/install/lib64 -lEvtGen -lHepMC3search -lEvtGenExternal -lpythia8 -L/home/mcolonna/TAUOLA/lib -lTauolaCxxInterface -lTauolaFortran -lTauolaHEPEVT -lTauolaHepMC3 -lgfortran -lHepMC3
EVTGENHEADERS = -I/home/mcolonna/evtgen/install/include
EVTGENSHARED = -Wl,-rpath,/home/mcolonna/evtgen/install/lib64
HEPMCGENLIBS = -L/home/mcolonna/HepMC3/install/lib64 -lHepMC3
HEPMCGENHEADERS = -I/home/mcolonna/HepMC3/install/include
HEPMCGENSHARED = -Wl,-rpath,/home/mcolonna/HepMC3/install/lib64
TAUOLALIBS = -L/home/mcolonna/TAUOLA/lib 
TAUOLAHEADERS = -I/home/mcolonna/TAUOLA/include
TAUOLASHARED = -Wl,-rpath,/home/mcolonna/TAUOLA/lib


CXXFLAGS     += $(ROOTCFLAGS)

GCC = g++

#B0_tauola: B0_tauola.C
#	$(GCC) $(CXXFLAGS) $(ROOTLIBS) $(EVTGENLIBS) $(EVTGENHEADERS) $(EVTGENSHARED) $(HEPMCGENLIBS) $(HEPMCGENHEADERS) $(HEPMCGENSHARED) $(TAUOLALIBS) $(TAUOLAHEADERS) $(TAUOLASHARED) $@.C -o $@
B02DstTauNu_TauMu: B02DstTauNu_TauMu.C
	$(GCC) $(CXXFLAGS) $(ROOTLIBS) $(EVTGENLIBS) $(EVTGENHEADERS) $(EVTGENSHARED) $(HEPMCGENLIBS) $(HEPMCGENHEADERS) $(HEPMCGENSHARED) $(TAUOLALIBS) $(TAUOLAHEADERS) $(TAUOLASHARED) $@.C -o $@

B02DstMuNu: B02DstMuNu.C
	$(GCC) $(CXXFLAGS) $(ROOTLIBS) $(EVTGENLIBS) $(EVTGENHEADERS) $(EVTGENSHARED) $(HEPMCGENLIBS) $(HEPMCGENHEADERS) $(HEPMCGENSHARED) $(TAUOLALIBS) $(TAUOLAHEADERS) $(TAUOLASHARED) $@.C -o $@


# Clean up: remove executables and outdated files.
clean:
	rm -rf B02DstTauNu_TauMu
	rm -rf B02DstMuNu

