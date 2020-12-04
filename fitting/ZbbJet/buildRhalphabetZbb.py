#!/usr/bin/env python
import ROOT as r, sys, math, os
from ROOT import TFile, TTree, TChain, gPad, gDirectory
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array

r.gSystem.Load(os.getenv('CMSSW_BASE') + '/lib/' + os.getenv('SCRAM_ARCH') + '/libHiggsAnalysisCombinedLimit.so')
sys.path.insert(0, '../.')
from tools import *
from hist import *

MASS_BINS = 23
MASS_LO = 40       # mass range for RooVar
MASS_HI = 201
MASS_HIST_LO = 47   # mass range for histograms
#MASS_HIST_LO = 54
MASS_HIST_HI = 201
BLIND_LO = 110
BLIND_HI = 131
RHO_LO = -6
RHO_HI = -2.1

SF2018={
    'shift_SF'  : 0.970,           'shift_SF_ERR' : 0.012,
    'smear_SF'  : 0.9076,          'smear_SF_ERR' : 0.0146,
    'V_SF'      : 0.953,           'V_SF_ERR'     : 0.016,
    'BB_SF'     : 1.0,             'BB_SF_ERR' : 0.3,     ## M2 SF                                                                                                        
}
SF2017={
    'shift_SF'  : 0.978,           'shift_SF_ERR' : 0.012,
    'smear_SF'  : 0.9045,          'smear_SF_ERR' : 0.048,
    'V_SF'      : 0.924,           'V_SF_ERR'  : 0.018,
    'BB_SF'     : 1.0,             'BB_SF_ERR' : 0.3        , # prelim ddb SF                                                                                              
}
SF2016={
    'BB_SF'     : 1.0,             'BB_SF_ERR' : 0.23,     ## M2 SF
    'V_SF'      : 0.993,            'V_SF_ERR'  : 0.043,
    'shift_SF'  : 1.001,            'shift_SF_ERR' : 0.012   , # m_data/m_mc, sqrt((m_data_err/m_data)**2+(m_mc_err/m_mc)**2)                                           
    'smear_SF'  : 1.084,            'smear_SF_ERR' : 0.0905  ,
}

qcdTFpars_2018 ={'n_rho':2, 'n_pT':2,
                 'pars':[0.0139,-0.9680,2.3695 ,0.6775 ,1.0759 ,1.4427 ,0.1826 ,0.2077 ,1.8612 ,-0.8737],
                 'fitpath':'ddb2018_Jun24_v3/ddb_M2_full/TF22_MC_w2Fit/rhalphabase.root'}
qcdTFpars_2017={'n_rho':2, 'n_pT':2,
                'pars':[ 0.0151 , -1.0359, 2.3953 , 0.7093 , 1.0947 , 1.6930 , -0.1745, 0.1980 , 1.4567 , -0.0427],
                'fitpath':'ddb_Jun24_v2/ddb_M2_full/TF22_MC_w2Fit/rhalphabase.root'}
qcdTFpars_2016_dak8 ={#'n_rho':3,'n_pT':5, 
                      #'pars':[0.0201,
                      #        4.9506, -0.5974, 1.7446, 0.7833, -2.826, -1.3115, 6.0557, -1.7539,  2.2512, 
                      #        13.25453, -11.7266, 5.429, 7.832, -29.377, 27.9622, -15.8643, 0.27883, 29.89604, -11.7782, 14.5337],
                      #'fitpath':'dak8/dak8_M2_2016_May10_tightMatch_ttvetoM_withDDT/TF35_MC_w2Fitv2/rhalphabase.root',
                      #'n_rho':2, 'n_pT':2,
                      #'pars': [0.15477433346748115, 0.9122590385803235, 1.18977336688895, 0.9972772436722153, -0.830081470482444, 3.022917531298525, -0.5566103505875617, 2.615210712845645, -0.7017472329157322, 0.4469318465393677],
                      #'fitpath':'dak8/dak8_M2_2016_Sep10_MatejTemplates/TF22_MC_w2Fitv2/rhalphabase.root',
                      #Oct10
                      #'n_rho':4, 'n_pT':1,  
                      #'pars': [0.15417995041232063, -0.8207842585109368, 2.245589577757258, 0.5730538912760856, 1.2236886562588154, 1.0351736936041647, 1.1016206945092009, 1.9974507981591145, -2.763517461753434, 5.686545840579456, -4.016073689805843],
                      #'fitpath':'dak8/dak8_M2_2016_Oct8_MatejTemplates/TF41_MC_w2Fitv2_muonCR/rhalphabase.root',
                      'n_rho':4, 'n_pT':3,
                      'pars':[0.15554712421560596, 2.6165808326312714, -0.9454880352431942, 1.6443716645006674, -0.26245021260391965, 0.7971469920476046, -9.776790271572207, 15.318593552982087, -0.6101218164258739, 13.332759028448706, 3.2589267882801636, 17.125456560338428, -25.137397082525464, -12.144382161475114, -17.96603283237203, -14.488200097242665, -18.917520890017805, 29.999170853046905, 22.30966464096825, 24.437631270382226, 15.306685855146654],
                      'fitpath':'dak8/dak8_M2_2016_Oct8_MatejTemplates/TF43_MC_w2Fitv2_muonCR/rhalphabase.root',
                      #'n_rho':4, 'n_pT':3, 
                      #'pars': [-3.6605197900701256, 5.83200439853875, -1.8775811351152143, 2.4133811596545627, 0.5974912928687424, 20.068849557814502, -14.833823380891307, 9.401899963485967, 0.9336976299054527, -1.3486751443559122, -25.3906669701001, 16.63733561007404, -0.13284160536213463, -6.675899141425084, 3.3652333529318597, 29.976583306549053, -11.223702654217266, -4.0230139130106295, 17.89138973145517, -17.357952014597423],
                      #'fitpath':'dak8/tightTopVeto_2016/TF43_MC_w2Fitv2/rhalphabase.root',
                      #'n_rho':4, 'n_pT':3,
                      #'pars':[0.02012216832938594, -2.1245752517206924, 4.642350520743385, -1.0952863593351232, 2.125819045977998, 0.6400817635219482, 18.09067598406464, -13.524206910672952, 8.159227666670368, 1.5938158455076206, -1.5578556549905969, -21.721069806832972, 13.618575996538375, 2.7475908117196326, -8.303068594227735, 4.067853838371654, 26.164655421043356, -8.649124391553201, -6.678034076297433, 19.255508863004238, -17.797774928138846],
                      #'fitpath':'dak8/dak8_M2_2016_May10_tightMatch_ttvetoM_withDDT/TF43_MC_w2Fitv2/rhalphabase.root',
                      #'n_rho':3, 'n_pT':4,
                      #'pars':[0.02012216832938594, 4.882609080643348, -0.7825180500719675, 1.9452746069951488, 0.699752258725546, -6.188988214828072, 2.041255428456303, 4.2384682149208075, -1.2307220499218658, 15.699706656941188, 1.3643528905941835, -6.996883061214927, 2.946363828315384, -20.28729986279143, -7.770939022211714, 26.481927721827127, -13.404608598818896, 29.76224377560486, 14.174444785751547, -29.99996571889061, 5.807858832284083],
                      #'fitpath':'dak8/dak8_M2_2016_May10_tightMatch_ttvetoM_withDDT/TF34_MC_w2Fitv2/rhalphabase.root',
                      #'n_rho':2, 'n_pT':2,
                      #'pars':[0.020122168239502677, 0.38112978719457047, 1.8493237587961104, 0.7285856838418248, -0.4517927327288014, 3.779697562244934, -1.393937860245007\
                      #        , 3.0742854675930147, 0.15761645878317054, -1.6113936238170368],
                      #'fitpath':'dak8/dak8_M2_2016_May10_tightMatch_ttvetoM_withDDT/TF22_MC_w2Fitv2/rhalphabase.root'
                      }
qcdTFpars_2016_ddb = {'n_rho':2, 'n_pT':2,
                      'pars':[ 0.012287,-1.0841,2.0706, 0.82500,1.2746 ,1.8651 ,0.24751,0.43354 ,1.2328 ,-0.00224],
                      'fitpath':'ddb/ddb_M2_2016/TF22_MC_w2Fitv2/rhalphabase.root'}

def main(options, args):
    ifile = options.ifile
    odir = options.odir

    from rhalphabet_builder import RhalphabetBuilder, LoadHistograms, GetSF, TestSmoothing
    ifiles = {'qcd700': '/uscms_data/d3/cmantill/bacon/CMSSW_10_2_13/src/ZbbFit/fitting/ZbbJet/QCDHT700.root',
              'qcd1000': '/uscms_data/d3/cmantill/bacon/CMSSW_10_2_13/src/ZbbFit/fitting/ZbbJet/QCDHT1000.root',
              'qcd1500': '/uscms_data/d3/cmantill/bacon/CMSSW_10_2_13/src/ZbbFit/fitting/ZbbJet/QCDHT1500.root',
              'qcd2000': '/uscms_data/d3/cmantill/bacon/CMSSW_10_2_13/src/ZbbFit/fitting/ZbbJet/QCDHT2000.root',
          }
    TestSmoothing(ifiles,4,1)

    '''
    # Load the input histograms
    # 	- 2D histograms of pass and fail mass,pT distributions
    # 	- for each MC sample and the data
    f = r.TFile.Open(ifile)    
    fLoose = None
    qcdTFpars = {}
    if options.ifile_loose is not None:
        fLoose = r.TFile.Open(options.ifile_loose)
    if   options.year =='2018':
          sf=SF2018
          if not options.pseudo:    qcdTFpars = qcdTFpars_2018 
    elif options.year =='2017':
          sf=SF2017
          if not options.pseudo:    qcdTFpars = qcdTFpars_2017
    elif options.year =='2016':
          sf=SF2016
          if not options.pseudo:  
              if 'dak8' in ifile: # dak8
                  qcdTFpars = qcdTFpars_2016_dak8
              else: # ddb
                  qcdTFpars = qcdTFpars_2016_ddb

    (pass_hists,fail_hists) = LoadHistograms(f, options.pseudo, options.blind, options.useQCD, scale=options.scale, r_signal=options.r, mass_range=[MASS_HIST_LO, MASS_HIST_HI], blind_range=[BLIND_LO, BLIND_HI], rho_range=[RHO_LO,RHO_HI], fLoose=fLoose,sf_dict=sf,createPassFromFail=options.createPassFromFail,skipQCD=options.skipQCD)

    # Build the workspacees
    rhalphabuilder = RhalphabetBuilder(pass_hists, fail_hists, f, options.odir, nr=options.NR, np=options.NP, mass_nbins=MASS_BINS, mass_lo=MASS_LO, mass_hi=MASS_HI, blind_lo=BLIND_LO, blind_hi=BLIND_HI, rho_lo=RHO_LO, rho_hi=RHO_HI, blind=options.blind, mass_fit=options.massfit, freeze_poly=options.freeze, remove_unmatched=options.removeUnmatched, input_file_loose=fLoose,suffix=options.suffix,sf_dict=sf,mass_hist_lo=MASS_HIST_LO,mass_hist_hi=MASS_HIST_HI,qcdTFpars=qcdTFpars,exp=options.exp,multi=options.multi,pseudo=options.pseudo)
    rhalphabuilder.run()
    if options.prefit:
        rhalphabuilder.prefit()
    elif options.loadfit is not None:
        rhalphabuilder.loadfit(options.loadfit)
    if not options.pseudo:
        rhalphabuilder.createdeco()
    '''
##-------------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('-i', '--ifile', dest='ifile', default='hist_1DZbb.root', help='file with histogram inputs',
                      metavar='ifile')
    parser.add_option('--ifile-loose', dest='ifile_loose', default=None, help='second file with histogram inputs (looser b-tag cut to take W/Z/H templates)',
                      metavar='ifile_loose')
    parser.add_option('-o', '--odir', dest='odir', default='./', help='directory to write plots', metavar='odir')
    parser.add_option('--pseudo', action='store_true', dest='pseudo', default=False, help='use MC', metavar='pseudo')
    parser.add_option('--blind', action='store_true', dest='blind', default=False, help='blind signal region',
                      metavar='blind')
    parser.add_option('--use-qcd', type='int', dest='useQCD', default=1, help='use real QCD MC',
                      metavar='useQCD')
    parser.add_option('--massfit', action='store_true', dest='massfit', default=False, help='mass fit or rho',
                      metavar='massfit')
    parser.add_option('--exp', action='store_true', dest='exp', default=False, help='use exp(bernstein poly) transfer function',
                      metavar='exp')
    parser.add_option('--multi', action='store_true', dest='multi', default=False, help='define RooMultiPdf',
                      metavar='multi')
    parser.add_option('--freeze', action='store_true', dest='freeze', default=False, help='freeze pol values',
                      metavar='freeze')
    parser.add_option('--scale', dest='scale', default=1, type='float',
                      help='scale factor to scale MC (assuming only using a fraction of the data)')
    parser.add_option('--nr', dest='NR', default=2, type='int', help='order of rho (or mass) polynomial')
    parser.add_option('--np', dest='NP', default=1, type='int', help='order of pt polynomial')
    parser.add_option('-r', dest='r', default=1, type='float', help='signal strength for MC pseudodataset')
    parser.add_option('--remove-unmatched', action='store_true', dest='removeUnmatched', default =False,help='remove unmatched', metavar='removeUnmatched')
    parser.add_option('--prefit', action='store_true', dest='prefit', default =False,help='do prefit', metavar='prefit')
    parser.add_option('--loadfit', dest='loadfit', default=None, help='load qcd polynomial parameters from alternative rhalphabase.root',metavar='loadfit')
    parser.add_option('-y' ,'--year', type='choice', dest='year', default ='2016',choices=['2016','2017','2018'],help='switch to use different year ', metavar='year')
    parser.add_option('--suffix', dest='suffix', default='', help='suffix for conflict variables',metavar='suffix')
    parser.add_option('--createPassFromFail', action='store_true', dest='createPassFromFail', default=False, help='Creating data_obs pass from data_obs fail', metavar='createPassFromFail')
    parser.add_option('--skipQCD', action='store_true', dest='skipQCD', default=False, help='skipQCD MC template', metavar='skipQCD')

    (options, args) = parser.parse_args()

    import tdrstyle

    tdrstyle.setTDRStyle()
    r.gStyle.SetPadTopMargin(0.10)
    r.gStyle.SetPadLeftMargin(0.16)
    r.gStyle.SetPadRightMargin(0.10)
    r.gStyle.SetPalette(1)
    r.gStyle.SetPaintTextFormat("1.1f")
    r.gStyle.SetOptFit(0000)
    r.gROOT.SetBatch()

    main(options, args)
##-------------------------------------------------------------------------------------
