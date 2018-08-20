#!/usr/bin/env python

"""
flag_bad_plcal_channels.py

Written by Craig Anderson 05/08/2018

This script will read in the leakage solutions calculated in the course of the XY phase and 
polariation calibration script process_polcal.py, and willflag channels where the leakage
amplitude exceeds some threshold. Typically, these channels are RFI afflicated and should be
flagged anyway. This step is important for obtaining good frequency-independent leakage estimates
subsequently, and greatly improves the quality of the data output from the polcal and 
cube imaging pipelines.
"""

# Import packages
import re
import numpy as np
import pyrap.tables as pt
from subprocess import call
import argparse
import glob

def main(args):

	# Init stuff

	# Set hard coded parameters

	####### Parse out cmdline args

	# Standardise inputs
	if args.baseDir[-1] != '/': args.baseDir += '/'
	baseDir = args.baseDir

	# Parse out args
	targetStr = args.targetStr
	if '-' in args.beams:
		beamsStrList = [str(i).zfill(2) for i in range(int(args.beams.split('-')[0]),int(args.beams.split('-')[1])+1)]
	else:
		beamsStrList = [str(args.beams).zfill(2)]

	flagLeakThreshLower = args.flagLeakThreshLower 
	flagLeakThreshUpper = args.flagLeakThreshUpper
	rotAntStr = str(args.rotAnt)
	anyAntBool = args.anyAnt


	#For each target ms:
	for targetBeam in beamsStrList:

		#Init stuff
		bad1MHzChansList = [] #Will be the list of 1 MHZ chans with bad leakage amps

		# Select target 
		targetMS = baseDir +'BPCAL/'+ '1934_SB'+targetStr+'_beam'+targetBeam+'_apply.ms'

		# Print info
		print '\n'*2
		print 'Flagging ' + targetMS
		print '\n'*2

		targetHandle = pt.table(targetMS,readonly=False) #table to write the data to
		targetFlagsTable = targetHandle.getcol('FLAG')


		#Determine the number of 1 MHz leakage solution channels that we have data for
		nChans = len(glob.glob(baseDir+'script_io/parset_leakages.b'+targetBeam+'_c*')) #get number of channels in ms

		#Check that the number of 1 MHZ channels we have data for matches the 1934 observations (which are fine-channelised)
		assert 54*nChans==targetFlagsTable.shape[1]


		###Pull in the leakages previously calculated by process_polcal.py for the XY phase correction

		#For each 1 MHz channel over which the leakages are calculated

		if anyAntBool: #If we are going to flag on bea leakage amplitudes on ALL antennas, then

			targetAntsTable = targetHandle.getcol('ANTENNA1')
			antsInArr = list(set(targetAntsTable))

			#Init some arrays that we can store solutions in
			L12LeakArr = np.zeros([nChans,len(antsInArr)],dtype='complex')
			L21LeakArr = np.zeros([nChans,len(antsInArr)],dtype='complex')

			for chan in range(nChans):
				for ant in antsInArr:

					try:
						#Open the parset with the leakage solutions
						with open(baseDir+'script_io/parset_leakages.b'+targetBeam+'_c'+str(chan)) as f:
							lines = f.read().splitlines()

						#Pull the leakages from the rotated PAF ant
						for line in lines:
							if re.findall('leakage.d12.'+str(ant)+'.*?\]',line):
								L12str = re.findall('\[.*?\]',line)[0]
								L12ComplexComps = L12str.replace('[','').replace(']','').split(',')
								L12Leak = float(L12ComplexComps[0])+float(L12ComplexComps[1])*1j
								L12LeakArr[chan,ant] = L12Leak

								#print L12LeakPhaseCorrectionAngle
							elif re.findall('leakage.d21.'+str(ant)+'.*?\]',line):
								L21str = re.findall('\[.*?\]',line)[0]
								L21ComplexComps = L21str.replace('[','').replace(']','').split(',')
								L21Leak = float(L21ComplexComps[0])+float(L21ComplexComps[1])*1j
								L21LeakArr[chan,ant] = L21Leak

					except Exception as e:
						print e

					#Once we have our leakages, determine if the channel needs flagging, and lag it if so
					maxLeakAmpL12 = np.nanmax(np.abs(L12LeakArr),axis=1)
					maxLeakAmpL21 = np.nanmax(np.abs(L21LeakArr),axis=1)
					minLeakAmpL12 = np.nanmin(np.abs(L12LeakArr),axis=1)
					minLeakAmpL21 = np.nanmin(np.abs(L21LeakArr),axis=1)

					#If the leakage amplitude in this channel fails the criterial set down, add it to the kill list
					for idx in range(len(maxLeakAmpL12)):

						if maxLeakAmpL12[idx] > flagLeakThreshUpper or maxLeakAmpL12[idx] > flagLeakThreshUpper or minLeakAmpL12[idx] < flagLeakThreshLower or minLeakAmpL21[idx] < flagLeakThreshLower:
							bad1MHzChansList.append(idx)
						else:
							pass


		else: #If we're just looking at the rotated ant, then only look at the rotated ant

			#Init some arrays that we can store solutions in
			L12LeakArr = np.zeros(nChans,dtype='complex')
			L21LeakArr = np.zeros(nChans,dtype='complex')

			for chan in range(nChans):

				try:
					#Open the parset with the leakage solutions
					with open(baseDir+'script_io/parset_leakages.b'+targetBeam+'_c'+str(chan)) as f:
						lines = f.read().splitlines()

					#Pull the leakages from the rotated PAF ant
					for line in lines:
						if re.findall('leakage.d12.'+rotAntStr+'.*?\]',line):
							L12str = re.findall('\[.*?\]',line)[0]
							L12ComplexComps = L12str.replace('[','').replace(']','').split(',')
							L12Leak = float(L12ComplexComps[0])+float(L12ComplexComps[1])*1j
							L12LeakArr[chan] = L12Leak

							#print L12LeakPhaseCorrectionAngle
						elif re.findall('leakage.d21.'+rotAntStr+'.*?\]',line):
							L21str = re.findall('\[.*?\]',line)[0]
							L21ComplexComps = L21str.replace('[','').replace(']','').split(',')
							L21Leak = float(L21ComplexComps[0])+float(L21ComplexComps[1])*1j
							L21LeakArr[chan] = L21Leak

					#Once we have our leakages, determine if the channel needs flagging, and lag it if so
					maxLeakAmp = np.nanmax([np.abs(L12Leak),np.abs(L21Leak)])
					minLeakAmp = np.nanmin([np.abs(L12Leak),np.abs(L21Leak)])

					#If the leakage amplitude in this channel fails the criterial set down, add it to the kill list
					if maxLeakAmp > flagLeakThreshUpper or minLeakAmp < flagLeakThreshLower:
						bad1MHzChansList.append(chan)
					else:
						pass

				except Exception as e:
					print e


		###Fine channel flagging fo the 1934 data on a per beam basis.

		#Get fine channel indexes corresponding to the bad !MHz channels

		# Print info
		print '\n'
		print 'Flagging the following 1 MHz channel indexes with bad leakage solutions for: ' + targetMS
		print bad1MHzChansList
		print '\n'

		fineChannelIndexesToFlag = [i for i in range(targetFlagsTable.shape[1]) if i/54 in bad1MHzChansList]
		targetFlagsTable[:,fineChannelIndexesToFlag,:]=True #Set all visibilities to True (flagged) for these channels

		#Write the flagging tables back to their measurement sets
		targetHandle.putcol('FLAG',targetFlagsTable)

		#Flush changes and close
		if targetHandle.datachanged():
			targetHandle.flush()
			targetHandle.close()
			print 'Flags have been changed and flushed to the target MS.'
		else:
			print 'Write to flags table in target MS failed. Leaving table open...'

ap = argparse.ArgumentParser(description="""flag_bad_plcal_channels.py

Written by Craig Anderson 06/08/2018

This script will read in the leakage solutions calculated in the course of the XY phase and 
polariation calibration script process_polcal.py, and will flag channels where the leakage
amplitude is outside user-defiend thresholds. Typically, these channels are RFI afflicated and should be
flagged anyway. This step is important for obtaining good frequency-independent leakage estimates
subsequently, and greatly improves the quality of the data output from the polcal and 
cube imaging pipelines.""")

ap.add_argument('baseDir',help='Directory that the BPCAL and script_io directories sit in.',type=str)
ap.add_argument('targetStr',help='Target SB (e.g. 5018)',type=str)
ap.add_argument('beams',help='The beams to which we will apply the flagging. Ranges are accepted, must be specified with a hypen (e.g. 2-7), and are inclusive. Single beams are specified with an integer.',type=str)
ap.add_argument('--flagLeakThreshUpper','-tu',help='The leakage amplitude above which a 1 MHz channel will be flagged. 8% is the expected leakage. [default=0.12]',default=0.12,type=float)
ap.add_argument('--flagLeakThreshLower','-tl',help='The leakage amplitude below which a 1 MHz channel will be flagged. 8% is the expected leakage. [default=0.0]',default=0.0,type=float)
ap.add_argument('--rotAnt','-r',help='Antenna index corresponding to rotated PAF [default 0]',default=0,type=int)
ap.add_argument('--anyAnt','-a',help='If the leakage values are bad on ANY antenna (again, typically indicating bad RFI), the channels will be flagged [default=false]',default=False,type=bool)
args = ap.parse_args()
main(args)

