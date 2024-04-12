#!/usr/bin/env python

#PyRoot file created from template located at:

from optparse import OptionParser
import os.path
import os
import sys
import subprocess
import glob
from array import array
from math import sqrt, exp

#Root stuff
from ROOT import TFile, TTree, TBranch, TLorentzVector, TLorentzRotation, TVector3
from ROOT import TCanvas, TMath, TH2F, TH1F, TProfile, TRandom, TGraphErrors, TGraph, TLine, TLegend
from ROOT import gBenchmark, gDirectory, gROOT, gStyle, gPad, gSystem

#Other packages
# import uproot
# import numpy as np
# import pandas


#My Stuff
#Required: add to PYTHONPATH environment variable, e.g.
from jz_pyroot_helper import * # /u/home/jzarling/path/PythonPath on farm
#from jz_pyroot_FitMacros import *

# CUT_STRING = "((11<=PS_column_right&&PS_column_right<=24) || (28<=PS_column_right&&PS_column_right<=40) || (44<=PS_column_right&&PS_column_right<=57) || (125<=PS_column_right&&PS_column_right<=131))"
CUT_STRING = "28<=PS_column_right&&PS_column_right<=40"

# e.g. first one makes int_N0 -> TCanvas index 4
channel_to_4x4 = {
	0:4, 1:8, 2:12, 3:16, 
	4:3, 5:7, 6:11, 7:15, 
	8:2, 9:6, 10:10, 11:14, 
	12:1, 13:5, 14:9, 15:13, 
}


def main(argv):
	#Usage controls from OptionParser
	parser_usage = ""
	parser = OptionParser(usage = parser_usage)
	(options, args) = parser.parse_args(argv)
	if(len(args) != 2):
		parser.print_help()
		return
	
	fname=argv[0]
	strname_out=argv[1]
	

	f  = TFile.Open(fname)
	tr = f.Get("tr")
	gStyle.SetOptStat(0)

	# First create histogram
	c1 = TCanvas("c1","c1",1200,900)
	h_list = []
	for side in ["N","S"]:
	# for side in ["N"]:
		for i in range(16):
			var =  "ped_"+side+str(i)
			print "Creating hists for: " + var
			tr.Draw(var+":event>>h2(100,0,3E8,200,380,430)",CUT_STRING,"colz")
			h2_thischan = gPad.GetPrimitive("h2")
			h2_thischan.SetNameTitle("h_"+var,var)
			profile = h2_thischan.ProfileX()
			profile.SetErrorOption("s")
			h_list.append(profile.ProjectionX("_px","C=E"))
	# Then draw in 4x4 canvas
	del c1
	c1 = TCanvas("c1","c1",1200*2,900*2)
	c1.Divide(4,4)
	for side in ["N","S"]:
		for i in range(16):
			c1.cd( channel_to_4x4[i] )
			print "cd: " + str(channel_to_4x4[i])
			if(side=="N"): 
				h_list[i].GetYaxis().SetRangeUser(0,10)
				h_list[i].Draw()	
			if(side=="S"):
				h_list[i+16].GetYaxis().SetRangeUser(0,10)
				h_list[i+16].Draw()	
		c1.SaveAs("monitoring/"+strname_out+"_pedestalRMS_"+side+".png")
		
	print("Done ")



def SaveAllHists(fname,hist_list=[],delete_after_saving=True,open_opt="RECREATE",subdirname="",clearSubDir=False):
	f = TFile.Open(fname,open_opt)
	f.cd() # Need to cd before checking directories, otherwise f will be None
	# Make directory if it doesn't exist
	if(f.GetDirectory(subdirname)==None): 
		# print "Creating directory...."
		f.mkdir(subdirname)
	# Delete directory if it does exist
	elif clearSubDir:
		# print "Cleaning subdir..."
		curr_TDir = f.GetDirectory(subdirname)
		f.cd()
		f.rmdir(subdirname)
		f.mkdir(subdirname)
	f.cd()
	f.cd(subdirname)
	
	# If histograms to be saved were not specified, save ALL TH1F and TH2F objects currently open
	if(len(hist_list)==0):
		all_objects_list = gDirectory.GetList()
		for obj in all_objects_list:
			if("<class 'ROOT.TH1" or "<class 'ROOT.TH2" in str(type(obj))): obj.Write() # Save 1D and 2D histograms
	# Save histograms specified in list
	else:
		for h in hist_list: h.Write()
	f.Close()
	
	# Delete histograms, if desired
	if(delete_after_saving and len(hist_list)==0):
		for obj in all_objects_list: 
			if("<class 'ROOT.TH1" in str(type(obj)) or "<class 'ROOT.TH2" in str(type(obj))): 
				gDirectory.Delete(obj.GetName()) # Save 1D histograms			gDirectory.Delete(obj.GetName())
	if(delete_after_saving and len(hist_list)!=0): 
		for h in hist_list: del h
		
	return


if __name__ == "__main__":
   main(sys.argv[1:])

#####################################################
###### Things from my pyroot helper file ############
#####################################################


# shell_exec(command, stdin_str):
# # command is the string (space separated)

# returns a numpy array for whichever axis (or axis uncertainty) is specified
# str_which_array options are "X", "EX", "Y", or "EY", for whichever I want
# jz_tgraph2numpy(gr,str_which_array):
# jz_th1f2numpy(gr,str_which_array)   
# jzGetRandomColor()   
# jz_DressUpObject(tobj,NewNameStr,kColor=kBlack,kMarkerStyle=kFullCircle,title="",xtitle="",ytitle="",Opacity=1.0,MarkerSize=1.0):
# jz_get_hist_binnum(my_val,nbins,h_min,h_max) 
# jzPlotHistStack(h_list,legend_list,tag_name,rebin_factor=1,range_lo=-1000,range_hi=1000,legend_limits=[],SavePNG=True):


#####################################################
###### Simple PyRoot Examples            ############
#####################################################

# Some example code...

## Opening a TFile (read-only)
# f = TFile.Open("name.root")
## Creating a new TFile
# f = TFile.Open("name.root","RECREATE")
## Modifying an existing file 
# f = TFile.Open("name.root","UPDATE")

## Retrieving TObject from TFile
# h   = f.Get("hname")

## Arrays
# my_arr = array('d',[])

## TLegend
# legend = TLegend(xmin,ymin,xmax,ymax) #0.1 is lower limit of plot, 0.9 is upper limit (beyond on either side is labeling+whitespace)
# legend.AddEntry(h,"Label","pl")
# legend.Draw()

## Marking up a histogram for plotting
# def AddPlotCosmetics(fname,gr_name_str,kColor,kMarkerStyle):
## Need to return to this one...

## Saving a drawn histogram/graph/whatever
# c1.SaveAs("CompareSigmas.png") # .png, .pdf, .C, ...





