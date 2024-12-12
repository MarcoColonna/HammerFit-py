from math import floor, sqrt
from sys import path
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot
import numpy as np
from numpy import meshgrid, linspace, array
from io import FileIO
import matplotlib.pyplot as plt
import mplhep as hep
#plt.style.use([hep.style.LHCb2])
import ROOT
import json


import importlib.util
import sys
from pathlib import Path

# If you dispose of a conda environment having the Hammer-Version you want installed
# the next lines will be (luckily) outdated
# just replace with:

from hammer import hepmc, pdg
from hammer.hammerlib import (FourMomentum, Hammer, IOBuffer,
                              Particle, Process, RecordType)



# otherwise if you have hammer built in your personal workspace you can link the python interface as following

#hammer_path = '/home/mcolonna/Hammer_1.3/Hammer-install/lib64/python3.10/site-packages/hammer'

#spec = importlib.util.spec_from_file_location("hammer", f"{hammer_path}/__init__.py")
#hammer = importlib.util.module_from_spec(spec)
#sys.modules["hammer"] = hammer
#spec.loader.exec_module(hammer)

#hammerlib_path = f"{hammer_path}/hammerlib.so"
#spec_hammerlib = importlib.util.spec_from_file_location("hammer.hammerlib", hammerlib_path)
#hammerlib = importlib.util.module_from_spec(spec_hammerlib)
#sys.modules["hammer.hammerlib"] = hammerlib
#spec_hammerlib.loader.exec_module(hammerlib)

#Hammer = hammerlib.Hammer
#IOBuffer = hammerlib.IOBuffer
#RecordType = hammerlib.RecordType

#######

# the histo_info is a utility class mainly used for plotting
class histo_info:
    def __init__(self, axis_titles, binning):
        self._axis_titles = axis_titles
        self._binning = binning

# the hammer cacher class handles directly the hammer histogram
# it access it and it changes, if required, the FF and the WC d.o.f
# giving access to the histogram as it changes wrt them
class HammerCacher:
    def __init__(self, fileName, histoName, FFscheme, WilsonSet, FormFactors, WilsonCoefficients, scaleFactor, histo_infos, verbose=False):#, **kwargs):
        self._histoName = histoName 
        self._FFScheme = FFscheme
        self._WilsonSet = WilsonSet 
        self._scaleFactor = scaleFactor
        
        self._wcs = WilsonCoefficients
        self._FFs = FormFactors
        
        self._nobs = 1
        self._strides = [1]
        self._ham = Hammer()
        self._ham.set_units("GeV")

        self._histo_infos = histo_infos

        buf = IOBuffer(RecordType.UNDEFINED)
        if(verbose):
            print(f"fileName = {fileName}")
            print(f"histoName = {histoName}")

        with open(fileName, 'rb', buffering=0) as fin:
            if buf.load(fin) and self._ham.load_run_header(buf):
                self._ham.init_run()
                if buf.load(fin):
                    while buf.kind == RecordType.HISTOGRAM or buf.kind == RecordType.HISTOGRAM_DEFINITION:
                        if buf.kind == RecordType.HISTOGRAM_DEFINITION:
                            name = self._ham.load_histogram_definition(buf)
                        else:
                            info = self._ham.load_histogram(buf)
                        if not buf.load(fin):
                            break
    
        self._ham.set_ff_eigenvectors(self._FFScheme["Process"],self._FFScheme["SchemeVar"],self._FFs)
        self._ham.set_wilson_coefficients(self._WilsonSet, self._wcs)
        self._histo = self._ham.get_histogram(histoName, FFscheme["name"])
        dims = self._ham.get_histogram_shape(histoName)
        dims = dims[1:] + dims[:1]
        dims.pop()
        
        ndims = self._ham.get_histogram_shape(histoName)
        
        for ndim in ndims:
            self._nobs*=ndim
        for dim in dims:
            self._strides = [c * dim for c in self._strides]
            self._strides.append(1)
        
        self._normFactor = self.getHistoTotalSM()
        
    def calcPos(self, indices):
        sum_list = []  
        total = 0 
        sum_list = [i * j for i, j in zip(indices, self._strides)]  
        total = sum(sum_list)
        return total
    
    def checkWCCache(self, wcs):
        isCached = True
        for key in wcs:
            if key not in self._wcs.keys():
                self._wcs[key] = wcs[key]
                isCached = False
            elif not (self._wcs[key] - wcs[key])==0:
                self._wcs[key] = wcs[key]
                isCached = False
        return isCached

    def checkFFCache(self, FFs):
        isCached = True
        for key in FFs:
            if key not in self._FFs.keys():
                self._FFs[key] = FFs[key]
                isCached = False
            elif not (self._FFs[key] - FFs[key])==0:
                self._FFs[key] = FFs[key]
                isCached = False
        return isCached

    def getHistoElement(self ,indices, wcs, FFs):
        pos = calcPos(indices)
        if not self.checkFFCache(FFs):
            self._ham.set_ff_eigenvectors(self._FFScheme["Process"], self._FFScheme["SchemeVar"], FFs)
            self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        if not self.checkWCCache(wcs):
            self._ham.reset_wilson_coefficients(self._WilsonSet)
            self._ham.set_wilson_coefficients(self._WilsonSet, wcs)
            self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        return self._histo[pos].sum_wi * scaleFactor
    
    def getHistoTotalSM(self):
        total = 0
        wcs = {}
        for key, value in self._wcs.items():
            if key == 'SM':
                wcs[key] = 1.
            else:
                wcs[key] = 0.
        self._ham.reset_wilson_coefficients(self._WilsonSet)
        self._ham.set_wilson_coefficients(self._WilsonSet, wcs)
        self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        for ni in range(self._nobs):
            total += self._histo[ni].sum_wi
        return total
    
    def getHistoTotal(self, wcs, FFs):
        total = 0
        self._ham.set_ff_eigenvectors(self._FFScheme["Process"], self._FFScheme["SchemeVar"], FFs)
        self.ham.reset_wilson_coefficients(self._WilsonSet, wcs)
        self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        for ni in range(nobs):
            total += _histo[ni].sumWi
        return total
    
    def getHistoElementByPos(self, pos, wcs, FFs):
        if not self.checkFFCache(FFs):
            self._ham.reset_ff_eigenvectors(self._FFScheme["Process"], self._FFScheme["SchemeVar"])
            self._ham.set_ff_eigenvectors(self._FFScheme["Process"], self._FFScheme["SchemeVar"], FFs)
            self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        if not self.checkWCCache(wcs):
            self._ham.reset_wilson_coefficients(self._WilsonSet)
            self._ham.set_wilson_coefficients(self._WilsonSet, wcs)
            self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        return self._histo[pos].sum_wi * self._scaleFactor / self._normFactor
    
    def getHistoElementByPosNoScale(self, pos, wcs, FFs):
        if not self.checkFFCache(FFs):
            self._ham.reset_ff_eigenvectors(self._FFScheme["Process"], self._FFScheme["SchemeVar"])
            self._ham.set_ff_eigenvectors(self._FFScheme["Process"], self._FFScheme["SchemeVar"], FFs)
            self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        if not self.checkWCCache(wcs):
            self._ham.reset_wilson_coefficients(self._WilsonSet)
            self._ham.set_wilson_coefficients(self._WilsonSet, wcs)
            self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        return self._histo[pos].sum_wi

    def getHistoElementByPosSM(self, pos, wcs, FFs):
        for key in wcs.keys():
            if key != 'SM':
                wcs[key] = 0.
        if not self.checkFFCache(FFs):
            self._ham.set_ff_eigenvectors(self._FFScheme["Process"], self._FFScheme["SchemeVar"], FFs)
            self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        if not self.checkWCCache(wcs):
            self._ham.reset_wilson_coefficients(self._WilsonSet)
            self._ham.set_wilson_coefficients(self._WilsonSet, wcs)
            self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        return self._histo[pos].sum_wi * self._scaleFactor / self._normFactor

    def getHistoElementByPosNoScaleSM(self, pos, wcs, FFs):
        for key in wcs.keys():
            if key != 'SM':
                wcs[key] = 0.
        if not self.checkFFCache(FFs):
            self._ham.set_ff_eigenvectors(self._FFScheme["Process"], self._FFScheme["SchemeVar"], FFs)
            self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        if not self.checkWCCache(wcs):
            self._ham.reset_wilson_coefficients(self._WilsonSet)
            self._ham.set_wilson_coefficients(self._WilsonSet, wcs)
            self._histo = self._ham.get_histogram(self._histoName, self._FFScheme["name"])
        return self._histo[pos].sum_wi

# Multi hammer cacher allows you to store multiple histograms in multiple files
# and treat them like a single one for when you parallelize the hammer reweighting process
# that can be very time consuming
class MultiHammerCacher:
    def __init__(self, cacherList):
        self._cacherList = []
        self._normFactor = 0
        self._scaleFactor = cacherList[0]._scaleFactor
        self._nobs = cacherList[0]._nobs
        self._strides = cacherList[0]._strides
        self._wcs = cacherList[0]._wcs
        self._FFs = cacherList[0]._FFs
        self._histo_infos = cacherList[0]._histo_infos
        for cacher in cacherList:
            self._cacherList.append(cacher)
            self._normFactor += cacher.getHistoTotalSM()
    
    def getHistoElementByPos(self, pos, wcs, FFs):
        res = 0
        for i in range(len(self._cacherList)):
            res += self._cacherList[i].getHistoElementByPosNoScale(pos,wcs,FFs)
            self._wcs = wcs
            self._FFs = FFs
        return res * self._scaleFactor / self._normFactor
    
    def getHistoElementByPosSM(self, pos, wcs, FFs):
        res = 0
        for key in wcs.keys():
            if key != 'SM':
                wcs[key] = 0.
        for i in range(len(self._cacherList)):
            res += self._cacherList[i].getHistoElementByPosNoScale(pos,wcs,FFs)
            self._wcs = wcs
            self._FFs = FFs 
        return res * self._scaleFactor / self._normFactor

# the background cacher access not hammer reweighted histograms and gives us in a format
# similar to the HammerCacher (easier to handle them together later)
class BackgroundCacher:
    def __init__(self, fileName, histoName, strides):
        self._fileName = fileName
        self._histoName = histoName
        self._strides = strides

        file = ROOT.TFile.Open(self._fileName, "READ")
        if not file or file.IsZombie():
            print("Error: Could not open file.")
            return
        
        hist = file.Get(self._histoName)
        if not hist:
            print(f"Error: Histogram '{self._histoName}' not found in file '{self._fileName}'.")
            file.Close()
            return

        if not isinstance(hist, ROOT.TH1):
            print(f"Error: '{self._histoName}' is not a 1D histogram.")
            file.Close()
            return

        self._histo = np.array([hist.GetBinContent(bin_idx) for bin_idx in range(1, hist.GetNbinsX() + 1)])

        file.Close()
        self._nobs=len(self._histo)
        self._normFactor=self._histo.sum()

    def getHistoElementByPos(self, pos,wcs,FFs):
        return self._histo[pos] / self._normFactor

        
# here we define the multiplicative nuisance parameters to apply to the hammer
# reweighted histogram
# the wrapper can change the d.o.f of the contribution and returns the content of
# a given bin wrt to the current values of the d.o.f. with an evaluate function
class HammerNuisWrapper:
    def __init__(self, hac, **kwargs):
        self._hac = hac
        self._nobs = hac._nobs
        self._wcs = hac._wcs
        self._FFs = hac._FFs
        self._params = {}
        for key, value in kwargs.items():
            self._params[key] = value
        self._nbin = 0
        self._strides = hac._strides
        self._dim = len(hac._strides)
        self._histo_infos = hac._histo_infos
    
    def set_wcs(self,wcs):
        self._wcs = {"SM":wcs[list(wcs.keys())[0]],"S_qLlL": complex(wcs[list(wcs.keys())[1]], wcs[list(wcs.keys())[2]]),"S_qRlL": complex(wcs[list(wcs.keys())[3]], wcs[list(wcs.keys())[4]]),"V_qLlL": complex(wcs[list(wcs.keys())[5]], wcs[list(wcs.keys())[6]]),"V_qRlL": complex(wcs[list(wcs.keys())[7]], wcs[list(wcs.keys())[8]]),"T_qLlL": complex(wcs[list(wcs.keys())[9]], wcs[list(wcs.keys())[10]])}
    
    def set_FFs(self,FFs):
        FFs_temp = {}
        for key, value in FFs.items():
            if key in self._FFs.keys():
                FFs_temp[key] = float(value)
        self._FFs = FFs_temp
    
    def set_params(self,params):
        params_temp = {}
        for key, value in params.items():
            if key in self._params.keys():
                params_temp[key] = value
        self._params = params_temp
        
    def set_nbin(self,nbin):
        self._nbin = nbin
    
    def evaluate(self):
        val = self._hac.getHistoElementByPos(self._nbin, self._wcs, self._FFs)
        for key, value in self._params.items():
            val = val*value
        return val

# this is the same as HammerNuisWrapper but the evaluate method is
# ignoring the WCs
# if you for example want to inject new physics B2DTauNu and not in B2DMuNu
# you'll build a ordinary HammerNuisWrapper for B2DTauNu
# and a SM one for B2DMuNu
class HammerNuisWrapperSM:
    def __init__(self, hac, **kwargs):
        self._hac = hac
        self._nobs = hac._nobs
        self._wcs = hac._wcs
        self._FFs = hac._FFs
        self._params = {}
        for key, value in kwargs.items():
            self._params[key] = value
        self._nbin = 0
        self._strides = hac._strides
        self._dim = len(hac._strides)
        self._histo_infos = hac._histo_infos
    
    def set_wcs(self,wcs):
        #for key in self._wcs.keys():
        #    if key == 'SM':
        #        self._wcs[key] = wcs[key]
        #    else:
        #        self._wcs[key] = complex(wcs['Re_'+key],wcs['Im_'+key])
        self._wcs = {"SM":wcs[list(wcs.keys())[0]],"S_qLlL": complex(wcs[list(wcs.keys())[1]], wcs[list(wcs.keys())[2]]),"S_qRlL": complex(wcs[list(wcs.keys())[3]], wcs[list(wcs.keys())[4]]),"V_qLlL": complex(wcs[list(wcs.keys())[5]], wcs[list(wcs.keys())[6]]),"V_qRlL": complex(wcs[list(wcs.keys())[7]], wcs[list(wcs.keys())[8]]),"T_qLlL": complex(wcs[list(wcs.keys())[9]], wcs[list(wcs.keys())[10]])}

    def set_FFs(self,FFs):
        FFs_temp = {}
        for key, value in FFs.items():
            if key in self._FFs.keys():
                FFs_temp[key] = float(value)
        self._FFs = FFs_temp
    
    def set_params(self,params):
        params_temp = {}
        for key, value in params.items():
            if key in self._params.keys():
                params_temp[key] = value
        self._params = params_temp
        
    def set_nbin(self,nbin):
        self._nbin = nbin
    
    def evaluate(self):
        val = self._hac.getHistoElementByPosSM(self._nbin, self._wcs, self._FFs)
        for key, value in self._params.items():
            val = val*value
        return val

# this attaches Nuisance parameters to the BackgroundCacher
# notice that one of the Nuisance parameters here should always be the yield
class BackgroundNuisWrapper:
    def __init__(self, bkg, **kwargs):
        self._bkg = bkg
        self._nobs = bkg._nobs
        self._params = {}
        self._wcs={}
        self._FFs={}
        for key, value in kwargs.items():
            self._params[key] = value
        self._nbin = 0
        self._strides = bkg._strides
        self._dim = len(bkg._strides)


    def set_nbin(self,nbin):
        self._nbin = nbin

    def set_wcs(self,wcs):
        self._wcs = {}

    def set_FFs(self,FFs):
        self._FFs = {}

    def set_params(self,params):
        params_temp = {}
        for key, value in params.items():
            if key in self._params.keys():
                params_temp[key] = value
        self._params = params_temp
    
    def evaluate(self):
        val = self._bkg.getHistoElementByPos(self._nbin, self._wcs, self._FFs)
        for key, value in self._params.items():
            val = val*value
        return val
    
# the template class takes the wrapper and allows to generate templates, and toys
# wrt any set of d.o.f we want
class template:
    def __init__(self, name, wrap):
        self._name = name
        self._wrap = wrap
        self._nobs = wrap._nobs
        self._nwcs = len(self._wrap._wcs)
        self._nFFs = len(self._wrap._FFs)
        self._nparams = len(self._wrap._params)
        self._strides = wrap._strides
        self._histo_infos = wrap._histo_infos
    
    def generate_template(self, **kwargs):
        wcs = {}
        FFs = {}
        params = {}

        for i, (key, value) in enumerate(kwargs.items()):
            if i < self._nwcs*2-1:
                wcs[key] = value
            elif self._nwcs*2-1 <= i < self._nwcs*2-1+self._nFFs:
                FFs[key] = value
            else:
                params[key] = value
        self._wrap.set_wcs(wcs)
        self._wrap.set_FFs(FFs)
        self._wrap.set_params(params)

        bin_contents = np.zeros(self._nobs)

        for i in range(self._nobs):
            self._wrap.set_nbin(i)
            val=self._wrap.evaluate()
            bin_contents[i]+=val
        
        return bin_contents

    def generate_toy(self, **kwargs):
        wcs = {}
        FFs = {}
        params = {}

        for i, (key, value) in enumerate(kwargs.items()):
            if i < self._nwcs*2-1:
                wcs[key] = value
            elif self._nwcs*2-1 <= i < self._nwcs*2-1+self._nFFs:
                FFs[key] = value
            else:
                params[key] = value
        self._wrap.set_wcs(wcs)
        self._wrap.set_FFs(FFs)
        self._wrap.set_params(params)

        bin_contents = np.zeros(self._nobs)

        for i in range(self._nobs):
            self._wrap.set_nbin(i)
            val=self._wrap.evaluate()
            bin_contents[i]+=np.random.poisson(val)
        
        return bin_contents 

# the fitter contains a template list and data (toys in the examples)
# it contains the definition of a nul_pdf and an alternative_pdf to be injected in the definition of the modifier
# a small plotting interface is implemented to retireve the projected histograms (from the strides) and overlay data
class fitter:
    def __init__(self,template_list,nul_params):
        self._template_list = template_list
        self._data = np.array([])
        self._nul_params = nul_params
    
    def nul_pdf(self):
        def func():
            res = np.zeros(self._template_list[0]._nobs)
            for temp in self._template_list:
                res += temp.generate_template(**self._nul_params)
            return res
        return func        

    def alt_pdf(self):
        def func(**kwargs):
            res = np.zeros(self._template_list[0]._nobs)
            for temp in self._template_list:
                res += temp.generate_template(**kwargs)
            return res
        return func
    
    def upload_data(self,data):
        self._data = data

    def get_histos(self,input):
        v_out = []
        nobs = self._template_list[0]._nobs
        strides = self._template_list[0]._strides
        n_histos = len(strides)
        dim = np.zeros(n_histos)
        index = np.zeros(n_histos)

        for i in range(n_histos):
            if(i==0):
                dim[i]=int(nobs/strides[i])
            else:
                dim[i]=int(strides[i-1]/strides[i])

        for i in range(n_histos):
            histo = np.zeros(int(dim[i]))
            v_out.append(histo)    
        
        for i in range(nobs):
            for j in range(n_histos):
                index[j]=i
                for k in range(j):
                    index[j]-=(index[k]*strides[k])
                index[j] = int(index[j]/strides[j])
                v_out[j][int(index[j])]+=input[i]

        return v_out

    def plot(self, **kwargs):
        strides = self._template_list[0]._strides
        n_histos = len(strides)
        axis_titles = self._template_list[0]._histo_infos._axis_titles
        binning = self._template_list[0]._histo_infos._binning
        contributions = []
        contributions.append(self.get_histos(self._data))
        for k in range(len(self._template_list)):
            contributions.append(self.get_histos(self._template_list[k].generate_template(**kwargs)))
        n_histos = len(contributions[0])
        n_contributions = len(contributions)

        colors = plt.cm.viridis(np.linspace(0, 1, n_contributions))

        fig = plt.figure(figsize=(8, 6 * n_histos))
        gs = gridspec.GridSpec(n_histos * 2, 1, height_ratios=[4, 1] * n_histos) 

        axs = []
        for i in range(n_histos):
            ax_main = fig.add_subplot(gs[i * 2])
            ax_ratio = fig.add_subplot(gs[i * 2 + 1], sharex=ax_main)
            axs.append((ax_main, ax_ratio))

        for i in range(n_histos):
            total_contribution = np.zeros(len(contributions[0][i]))
            for j in range(n_contributions):
                bin_content = contributions[j][i]
                n_bins = len(bin_content)
                bin_edges = np.linspace(binning[i][0], binning[i][1], n_bins + 1)

                # Main plot
                if j == 0:
                    axs[i][0].errorbar(
                        bin_edges[:-1] + np.diff(bin_edges) / 2,
                        bin_content,                           
                        yerr=np.sqrt(bin_content),              
                        fmt='+',                                 
                        alpha=1.,                                 
                        label='Data'                             
                    )
                else:
                    axs[i][0].bar(
                        bin_edges[:-1], bin_content, width=np.diff(bin_edges),
                        align='edge', alpha=0.5, color=colors[j], label=self._template_list[j-1]._name
                    )
                    total_contribution += bin_content

            axs[i][0].step(
                bin_edges, np.append(total_contribution, total_contribution[-1]),
                where='post', color='black', linestyle='--', linewidth=1.5, label='Total'
            )

            ratio = contributions[0][i] / total_contribution
            ratio_err = np.sqrt(contributions[0][i]) / total_contribution
            axs[i][1].errorbar(
                bin_edges[:-1] + np.diff(bin_edges) / 2, ratio, yerr=ratio_err,
                fmt='+', alpha=1., label='Data'
            )
            axs[i][1].axhline(1, color='red', linestyle='--', linewidth=1)

            axs[i][0].set_xlim(binning[i][0], binning[i][1])
            axs[i][0].set_title('')
            axs[i][0].set_xlabel(axis_titles[i])
            axs[i][0].set_ylabel('Bin Content')
            axs[i][0].legend(loc='best')

            # Formatting ratio plot
            axs[i][1].set_xlim(binning[i][0], binning[i][1])
            axs[i][1].set_ylim(0.8, 1.2)
            axs[i][1].set_xlabel('')
            axs[i][1].set_ylabel('Data/Model')
            axs[i][1].tick_params(axis='y', which='both', right=False)

            plt.setp(axs[i][0].get_xticklabels(), visible=False)

        plt.tight_layout(h_pad=0.5)
        plt.show()

# The reader class aim is to make everything above not necessary to be fully undestood
# A config file is provided and the reader produces itself the necessary objects:
# Cachers -> Wrappers -> Templates -> Fitter (returned)
# giving access to a fitter with a toy stored inside as data (temporary)
class Reader:
    def __init__(self, filename):
        self.name = filename
        with open(filename, 'r') as f:
            self.config = json.load(f)
        
    def createFitter(self, verbose=False):
        template_list = []
        nul_params = {}

        for mode, mode_config in self.config.items():
            hac_list = []
            if verbose:
                print(f"Reading the mode: {mode}")
            fileNames = mode_config["fileNames"]
            histoname = mode_config["histoname"]
            ffscheme = mode_config["ffscheme"]
            wcscheme = mode_config["wcscheme"]
            formfactors = mode_config["formfactors"]
            wilsoncoefficients = mode_config["wilsoncoefficients"]
            scalefactor = mode_config["scalefactor"]
            nuisance = mode_config["nuisance"]
            is_hammer_weighted = mode_config["ishammerweighted"]
            injectNP = mode_config["injectNP"]
            histo_infos = histo_info(mode_config["axistitles"], mode_config["binning"])
            _wilsoncoefficients = {}
            for key, value in wilsoncoefficients.items():
                if key == 'SM':
                    if key not in nul_params:
                        nul_params[key] = value[0]
                else:
                    if 'Re_'+key not in nul_params:
                        nul_params['Re_'+key] = value[0]
                        nul_params['Im_'+key] = value[1]
            for key, value in wilsoncoefficients.items():
                _wilsoncoefficients[key] = complex(value[0],value[1])
            if is_hammer_weighted:
                for fileName in fileNames:
                    hac_list.append(HammerCacher(fileName, histoname, ffscheme, wcscheme, formfactors, _wilsoncoefficients, scalefactor, histo_infos))
                cacher = MultiHammerCacher(hac_list)
                if injectNP:
                    wrapper = HammerNuisWrapper(cacher, **nuisance)
                    temp = template(mode, wrapper)
                    template_list.append(temp)
                    parameters = formfactors | nuisance
                    for key, value in parameters.items():
                        if key not in nul_params:
                            nul_params[key] = value
                else:
                    wrapper = HammerNuisWrapperSM(cacher, **nuisance)
                    temp = template(mode, wrapper)
                    template_list.append(temp)
                    for key, value in parameters.items():
                        if key not in nul_params:
                            nul_params[key] = value
            else:
                for fileName in fileNames:
                    hac_list.append(BackgroundCacher(fileName, histoname))
                cacher = hac_list[0]
                wrapper = BackgroundNuisWrapper(cacher, **nuisance)
                temp = template(mode, wrapper)
                template_list.append(temp)
                parameters = nuisance
                for key, value in parameters.items():
                    if key not in nul_params:
                        nul_params[key] = value

        return fitter(template_list, nul_params)

