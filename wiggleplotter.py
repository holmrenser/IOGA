#!/usr/bin/python 
'''
Author: Rens Holmer
Plots wiggletracks with contigbreaks, mean coverage +/- stdev
'''
from __future__ import print_function
import sys, matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from pylab import *

y=[]
contigstart=[]
length=0
tracker='.'
with open(sys.argv[1],"rU") as infile:
	for line in infile:
		length=length+1
		y.append(float(line.split('\t')[2].strip()))
		if line.split('\t')[1] == '1':
			tracker+='.'
			#print(tracker,end='\r')
			contigstart.append(length)
		#if int(line.split('\t')[2].strip()) > 250:
			#print length, int(line.split('\t')[2].strip())


	cov_stdev=std(y)
	cov_mean=mean(y)

	ylim(0,1.5*max(y))

	for i in contigstart:
		pyplot.plot([i,i],[0,1.5*max(y)],color='0.50',linestyle='-',linewidth=0.5)

	pyplot.plot(y,'k-')
	pyplot.plot([0,length],[cov_mean,cov_mean],'r-') #mean
	pyplot.plot([0,length],[cov_mean+cov_stdev,cov_mean+cov_stdev],'r--') #+sdev
	pyplot.plot([0,length],[cov_mean-cov_stdev,cov_mean-cov_stdev],'r--') #-sdev


	suptitle('Mean coverage = ' + str(mean(y)),fontsize=18)
	figure=gcf()
	figure.set_size_inches(20,10)
	savefig(sys.argv[1]+'.png',dpi=100)
