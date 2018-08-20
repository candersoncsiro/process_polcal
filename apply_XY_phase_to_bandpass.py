#!/usr/bin/env python

"""
Script to apply the on-axis calibration to bandpass tables 
Original version by Craig Anderson
Updated by George Heald
"""

# Import packages
import re
import numpy as np
import pyrap.tables as pt
from subprocess import call
import argparse

def main(args):

	# Construct table names
	if args.basedir[-1] != '/': args.basedir += '/'
	while args.bptab[-1] == '/': args.bptab = args.bptab[:-1]
	inputTable = args.basedir + args.bptab
	outputTable = inputTable + args.extension

	# Clarify what's going on to the user
	print 'Input table is',inputTable
	print 'Output table will be',outputTable

	# Make a copy of the apply table before it is modified
	call(['rm', '-rvf', outputTable]);
	call(['cp', '-rv', inputTable, outputTable ]);

	# Get data
	t = pt.table(inputTable); #table to read the data from
	w = pt.table(outputTable,readonly=False); #table to write the data to
	d = t.getcol('BANDPASS');
	d_XYcorr = t.getcol('BANDPASS');

	#define the channel parameters.
	chans = args.bandwidth
	chanWidth = 54

	#Print particulars of the dataaset to the log file
	print('Specified bandwitdh is %d.'%chans)
	print('Hardcoded number of fine-channels per MHz is 54.')
	print('...so the bandpass table should have a length of %d'%(chans*chanWidth*2))
	print('The bandpass table actually has shape:')
	print(d.shape)

	assert d.shape[-1] == chans*chanWidth*2

	#specify the beam parameters
	beams = args.nbeams

	#calculate the corrected bandpass table
	for beam in range(beams):

		#create the beam string designator
		beamStr = "b%02d"%beam
		#beamStr = "b"+str(beam)  

		print('\n')
		print('Operating on beam '+beamStr+':')
		print('\n')

		# Collect XY phase across frequency to report stats later
		beamstat = []

		# init a stored (null) correction value in case that channel 0 leakage is bad
		L12LeakPhaseCorrectionMultiplierStor = 1+0j #multiply by length of 1, zero phase
		L21LeakPhaseCorrectionMultiplierStor = 1+0j
		L12LeakPhaseCorrectionAngleStor = np.angle(L12LeakPhaseCorrectionMultiplierStor)
		L21LeakPhaseCorrectionAngleStor = np.angle(L21LeakPhaseCorrectionMultiplierStor)

		#For each 1 MHz channel
		for i in range(chans):

			#create the channel string designator
			chanStr = "c%01d"%i

			#print 'Operating on beam '+beamStr+' and channel '+chanStr+'...'

			#Open the parset with the leakage solutions
			# NOTE: this string might later need to go in the arguments
			try:
				with open(args.basedir+'script_io/parset_leakages.'+beamStr+'_'+chanStr) as f:
					lines = f.read().splitlines()

					#Pull the leakages from the rotated PAF ant
					for line in lines: #go through each leakage solution entry in the parset...
						#If the calibration solution in question relates to the XY leakage of the rotated ant
						if re.findall('leakage.d12.%d.0 = \[.*?\]'%args.antenna,line):
							L12str = re.findall('\[.*?\]',line)[0]
							L12ComplexComps = re.findall(r"[-+]?\d*\.\d+|\d+",L12str)
							L12Leak = float(L12ComplexComps[0])+float(L12ComplexComps[1])*1j
							try:
								L12LeakPhaseCorrectionMultiplier = np.conj(L12Leak/np.abs(L12Leak)) #Signal in X feed is some fraction of the signal in the Y feed with a phase applied to it, which is calculated as / encapsulated in the complex leakages. We want to pre-apply that complex leakage phase to the X feed signal with calibration, so that the signal common to both feeds has a phase of zero.
								L12LeakPhaseCorrectionMultiplierStor = L12LeakPhaseCorrectionMultiplier
								L12LeakPhaseCorrectionAngle = np.angle(L12Leak/np.abs(L12Leak),deg=True)
								L12LeakPhaseCorrectionAngleStor = L12LeakPhaseCorrectionAngle
							 except:
								L12LeakPhaseCorrectionMultiplier = L12LeakPhaseCorrectionMultiplierStor
								L12LeakPhaseCorrectionAngle = L12LeakPhaseCorrectionAngleStor
						#If the calibration solution in question relates to the YX leakage of the rotated ant, do as above.
						elif re.findall('leakage.d21.%d.0 = \[.*?\]'%args.antenna,line): 
							L21str = re.findall('\[.*?\]',line)[0]
							L21ComplexComps = re.findall(r"[-+]?\d*\.\d+|\d+",L21str)
							L21Leak = float(L21ComplexComps[0])+float(L21ComplexComps[1])*1j
							try:
								L21LeakPhaseCorrectionMultiplier = np.conj(L21Leak/np.abs(L21Leak))
								L21LeakPhaseCorrectionMultiplierStor = L21LeakPhaseCorrectionMultiplier
								L21LeakPhaseCorrectionAngle = np.angle(L21Leak/np.abs(L21Leak),deg=True)
								L21LeakPhaseCorrectionAngleStor = L21LeakPhaseCorrectionAngle
							except:
								L21LeakPhaseCorrectionMultiplier = L21LeakPhaseCorrectionMultiplierStor
								L21LeakPhaseCorrectionAngle = L21LeakPhaseCorrectionAngleStor
						else:
							pass
			except:
				print('Exception for channel'+chanStr)
				L12LeakPhaseCorrectionMultiplier = L12LeakPhaseCorrectionMultiplierStor
				L21LeakPhaseCorrectionMultiplier = L21LeakPhaseCorrectionMultiplierStor

			#print the phase correction angle
			#print 'XY- and YX-leakage-derived phase corrections for channel %s and beam %s are %.2f and %.2f radians (mean %.2f radians).' %( chanStr, beamStr, np.angle(L12LeakPhaseCorrectionMultiplier), np.angle(np.conj(L21LeakPhaseCorrectionMultiplier)), np.angle(np.mean([L12LeakPhaseCorrectionMultiplier,np.conj(L21LeakPhaseCorrectionMultiplier)])) ) 
			beamstat.append(np.angle(np.mean([L12LeakPhaseCorrectionMultiplier,np.conj(L21LeakPhaseCorrectionMultiplier)])))

			#Calculate the fine filter-bank channel range to which the leakage corrections must be applied
			#Bandpass correction factors are formatted XYXYXY... (confirmed by comparison of table with 
			#pipeline-generated bandpass plots), so need twice the total number of channels to grab for a give freq chunk
			chLo = i*(2*chanWidth)+1 
			chHi = (i+1)*(2*chanWidth)

			#Apply the correction. For the dummy column, for the selected beam, for all antennas but the rotated one, 
			#for the selected channel range, apply the X-Y phase correction to the Y feed to offset it from the X feed.
			#Original Correction: d_XYcorr[0,beam,:,chLo:chHi:2] = d[0,beam,:,chLo:chHi:2]*L12LeakPhaseCorrectionMultiplier
			#Apply the leakage phase to the X pol on all ants

			#Correct for the XY phase in the current 1 MHz channel in the current beam. The correction will differ based on the sense of rotation of the PAF.
			if args.sense==1:
				d_XYcorr[0,beam,:,chLo:chHi:2] = d[0,beam,:,chLo:chHi:2]*L12LeakPhaseCorrectionMultiplier
				if np.abs(L12LeakPhaseCorrectionMultiplier)<1.01:
					print('Applying an XY phase correction gain of %.1f degrees in phase and %.1f in amplitude to channel range %d--%d, corresponding to 1 MHz channel index %d.'%(np.angle(L12LeakPhaseCorrectionMultiplier,deg=True),np.abs(L12LeakPhaseCorrectionMultiplier),chLo,chHi,i))
				else:
					print('BAD CHANNEL DETECTED! Using correction from last good channel...')
			elif args.sense==-1:  
				raise Exception                      
				#d_XYcorr[0,beam,:,chLo:chHi:2] = d[0,beam,:,chLo:chHi:2]*np.conj(L12LeakPhaseCorrectionMultiplier)
			else:
				raise Exception

		beamstat = np.array(beamstat)*180./np.pi
		print 'Mean and standard deviation of mean XY,YX phase corrections for beam %s are %.1f +/- %.1f degrees'%(beamStr,np.mean(beamstat),np.std(beamstat))

	#Write back to the calparameters table
	w.putcol('BANDPASS',d_XYcorr);

	#Flush changes and close
	if w.datachanged():
		w.flush();
		w.close();
		t.close();
		print 'Data has been changed and flushed to table.'
	else:
		print 'Write to table failed. Leaving table open...'

ap = argparse.ArgumentParser()
ap.add_argument('basedir',help='Directory within which the script_io subdirectory is found')
ap.add_argument('bptab',help='BP table on top of which to apply corrections, given relative to BASEDIR')
ap.add_argument('--antenna','-a',help='Antenna index corresponding to rotated PAF [default 0]',default=0,type=int)
ap.add_argument('--sense','-s',help='The sense of rotation of the rotated PAF [-1 == couterclockwise, 1 == clockwise; default -1]',default=-1,type=int)
ap.add_argument('--nbeams','-n',help='Number of beams in the observation [default 36]',default=36,type=int)
ap.add_argument('--bandwidth','-b',help='Bandwith of observation in MHz [default 192]',default=192,type=int)
ap.add_argument('--extension','-x',help='Extension to append to create output BP table [default .xy]',default='.xy')
args = ap.parse_args()
main(args)

