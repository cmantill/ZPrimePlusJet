import os
import glob
import math
import array
import sys
import time
from optparse import OptionParser
import ROOT
from ROOT import *
from math import sqrt

def getRatio(hist, reference):
    ratio = hist.Clone("%s_ratio"%hist.GetName())
    ratio.SetDirectory(0)
    ratio.SetLineColor(hist.GetLineColor())
    ratio.SetMarkerSize(hist.GetMarkerSize())
    ratio.GetXaxis().SetTitleFont(43)
    ratio.GetXaxis().SetLabelFont(43)
    ratio.GetYaxis().SetTitleFont(43)
    ratio.GetYaxis().SetLabelFont(43)
    ratio.GetYaxis().SetLabelSize(22)
    ratio.GetXaxis().SetLabelSize(22)
    ratio.GetXaxis().SetTitleSize(26)
    ratio.GetYaxis().SetTitleSize(26)
    ratio.GetYaxis().SetNdivisions(5)
    ratio.GetYaxis().SetNoExponent()
    for xbin in xrange(1,reference.GetNbinsX()+1):
        ref = reference.GetBinContent(xbin)
        val = hist.GetBinContent(xbin)
        refE = reference.GetBinError(xbin)
        valE = hist.GetBinError(xbin)
        try:
            ratio.SetBinContent(xbin, val/ref)
            ratio.SetBinError(xbin, math.sqrt( (val*refE/(ref**2))**2 + (valE/ref)**2 ))
        except ZeroDivisionError:
            ratio.SetBinError(xbin, 0.0)
    return ratio

def makeComparison(iHists,iLegs,iName,iTag):
    tag = ROOT.TLatex(0.3,0.5,iTag)
    tag.SetNDC();
    tag.SetTextFont(42)
    tag.SetTextSize(0.000001);

    iHists[0].SetLineColor(ROOT.kBlue)
    iHists[1].SetLineColor(ROOT.kOrange)
    for ic,h in enumerate(iHists):
        h.SetLineColor(ic+1)
        h.SetLineWidth(2)
        h.SetFillStyle(0)
        h.GetXaxis().SetLabelSize(0.04)
        h.GetXaxis().SetTitleOffset(1.1)
        h.GetXaxis().SetTitleSize(0.04)
        h.GetYaxis().SetLabelSize(0.04)
        h.GetYaxis().SetTitleOffset(1.2)
        h.GetYaxis().SetTitleSize(0.04)
        h.GetYaxis().SetTitle('Events (a.u)')
        #h.GetXaxis().SetTitle('
        #leg.AddEntry(h,iLegs[ic],"l")

    c = ROOT.TCanvas("c"+iName,"c"+iName,900,800)
    p2 = ROOT.TPad("pad2","pad2",0,0,1,0.31);
    p2.SetTopMargin(0);
    p2.SetBottomMargin(0.3);
    p2.SetLeftMargin(0.15)
    p2.SetRightMargin(0.03)
    p2.SetFillStyle(0);
    p2.Draw();
    p1 = ROOT.TPad("pad1","pad1",0,0.31,1,1);
    p1.SetBottomMargin(0);
    p1.SetLeftMargin(p2.GetLeftMargin())
    p1.SetRightMargin(p2.GetRightMargin())
    p1.Draw();
    p1.cd();

    mainframe = iHists[0].Clone('mainframe')
    mainframe.Reset('ICE')
    mainframe.GetXaxis().SetTitleFont(43)
    mainframe.GetXaxis().SetLabelFont(43)
    mainframe.GetYaxis().SetTitleFont(43)
    mainframe.GetYaxis().SetLabelFont(43)
    mainframe.GetYaxis().SetTitle('Events')
    mainframe.GetYaxis().SetLabelSize(22)
    mainframe.GetYaxis().SetTitleSize(26)
    mainframe.GetYaxis().SetTitleOffset(2.0)
    mainframe.GetXaxis().SetTitle('')
    mainframe.GetXaxis().SetLabelSize(0)
    mainframe.GetXaxis().SetTitleSize(0)
    mainframe.GetXaxis().SetTitleOffset(1.5)
    mainframe.GetYaxis().SetNoExponent()
    mainframe.Draw()

    for ic,s in enumerate(iHists):
        #s.Scale(1/s.Integral());
        if ic==0:
            s.SetMinimum(0.001)
            s.Draw("hist")
        else:
            s.Draw("histsame")

    tag.Draw()
    p1.BuildLegend();
    p2.cd()
    ratioframe = mainframe.Clone('ratioframe')
    ratioframe.Reset('ICE')
    ratioframe.GetYaxis().SetRangeUser(0.50,1.50)
    ratioframe.GetYaxis().SetTitle('Pass/Fail')
    ratioframe.GetXaxis().SetTitle('m_{SD} (GeV)')
    ratioframe.GetXaxis().SetLabelSize(22)
    ratioframe.GetXaxis().SetTitleSize(26)
    ratioframe.GetYaxis().SetNdivisions(5)
    ratioframe.GetYaxis().SetNoExponent()
    #ratioframe.GetYaxis().SetTitleOffset(mainframe.GetYaxis().GetTitleOffset())                                                                                                                            
    ratioframe.GetXaxis().SetTitleOffset(3.0)
    ratioframe.Draw()

    ratio = getRatio(iHists[0],iHists[1])
    ratio.SetMinimum(0.5)
    ratio.SetMaximum(1.5)
    ratioframe.GetYaxis().SetRangeUser(0.5,1.5)
    line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1.0,
                      ratio.GetXaxis().GetXmax(), 1.0)
    line.SetLineColor(ROOT.kGray)
    line.SetLineStyle(2)
    line.Draw()
    ratio.Draw()

    c.cd()
    c.Modified()
    c.Update()   
    c.SaveAs("plots/"+iName+".pdf")

def main(options,args):
    wp = options.wp
    proc = options.proc
    ofile = ROOT.TFile.Open(options.ifile,'read')
    ptbinlegs = [
        "450 < p_{T} < 500 GeV",
        "500 < p_{T} < 550 GeV",
        "550 < p_{T} < 600 GeV",
        "600 < p_{T} < 675 GeV",
        "675 < p_{T} < 800 GeV",
        "800 < p_{T} < 1200 GeV",
             ]
    for ptbin in range(0,6):
        for iq,qcdbin in enumerate(['qcd700','qcd1000','qcd1500','qcd2000']):
            print('pass%s_catp%i'%(qcdbin,ptbin))
            hp = ofile.Get('pass%s_catp%i'%(qcdbin,ptbin)).Clone()
            hp.SetDirectory(0)
            print('fail%s_catf%i'%(qcdbin,ptbin))
            hf = ofile.Get('fail%s_catf%i'%(qcdbin,ptbin)).Clone()
            hf.SetDirectory(0)
            hp_s = ofile.Get('pass%sout_pass_4_1__m_pt_catp%i'%(qcdbin,ptbin)).Clone()
            hp_s.SetDirectory(0)
            hf_s = ofile.Get('pass%sout_fail_4_1__m_pt_catf%i'%(qcdbin,ptbin)).Clone()
            hf_s.SetDirectory(0)
            makeComparison([hp,hp_s],['Pass','Pass Smooth'],'comppass_qcdbin%s_ptbin%i'%(qcdbin,ptbin),ptbinlegs[ptbin])
            makeComparison([hf,hf_s],['Fail','Fail Smooth'],'compfail_qcdbin%s_ptbin%i'%(qcdbin,ptbin),ptbinlegs[ptbin])
            if iq==0:
                hpass = hp
                hfail = hf
                hpass_s = hp_s
                hfail_s = hf_s
            else:
                hpass.Add(hp)
                hfail.Add(hf)
                hpass_s.Add(hp_s)
                hfail_s.Add(hf_s)

        makeComparison([hpass,hpass_s],['Pass','Pass Smooth'],'comppass_ptbin%i'%(ptbin),ptbinlegs[ptbin])
        print(hfail,hfail_s)
        makeComparison([hfail,hfail_s],['Fail','Fail Smooth'],'compfail_ptbin%i'%(ptbin),ptbinlegs[ptbin])

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('-i','--ifile', dest='ifile', default = 'input.root',help='directory with data', metavar='idir')
    parser.add_option('--wp',dest='wp', default = 'M', help='wp')
    parser.add_option('--proc',dest='proc',default = 'qcd', help='process')
    (options, args) = parser.parse_args()

    import tdrstyle
    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.10)
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetPaintTextFormat("1.1f")
    ROOT.gStyle.SetOptFit(0000)
    ROOT.gROOT.SetBatch()

    main(options,args)
