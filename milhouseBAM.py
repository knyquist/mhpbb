import matplotlib.pylab as pylab
import matplotlib.pyplot as pyplot
pylab.rcParams['figure.figsize'] = 20, 15
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = 20, 15
plt.rcParams.update({'font.size': 20})
import utiliBAM as ubm
import nyplot as ny
import numpy as np

def generateMilhousePlotsFromBAM( SMRTLinkInfo, num=1000 ):
	"""
	Generate Milhouse-style plots 
		<SMRTLinkIDs> contains SMRTLinkIDs for each chip. Can be a list of SMRTLink JobIDs or a single SMRTLink JobID
		<num>         is the number of alignments to subsample to. Default is 1000 alignments. 

	Kristofor Nyquist 1/26/2016
	"""

	SMRTLinkIDs    = [ chips[0] for chips in SMRTLinkInfo ]
	SMRTLinkLabels = [ chips[1] for chips in SMRTLinkInfo ]

	chip_data = ubm.getChipKPIs( SMRTLinkIDs, num )

	### generate the plots
	plotBar( chip_data, SMRTLinkIDs, SMRTLinkLabels, 'total nreads' )
	try:
		plotBar( chip_data, SMRTLinkIDs, SMRTLinkLabels, 'fraction mapped' )
	except:
		print 'was unable to generate fraction loaded barplot. Most likely was unable to find location of primary folder for opening up subreads.bam file'

	plotSurvival( chip_data, SMRTLinkIDs, SMRTLinkLabels, 'readlength' )
	plotBoxplot(  chip_data, SMRTLinkIDs, SMRTLinkLabels, 'readlength' )

	plotSurvival( chip_data, SMRTLinkIDs, SMRTLinkLabels, 'templatespan' )
	plotBoxplot(  chip_data, SMRTLinkIDs, SMRTLinkLabels, 'templatespan' )

	plotSurvival( chip_data, SMRTLinkIDs, SMRTLinkLabels, 'polreadlength' )
	plotBoxplot(  chip_data, SMRTLinkIDs, SMRTLinkLabels, 'polreadlength' )

	plotSurvival( chip_data, SMRTLinkIDs, SMRTLinkLabels, 'accuracy', legend_loc=3 )
	plotBoxplot(  chip_data, SMRTLinkIDs, SMRTLinkLabels, 'accuracy' )

	plotBoxplot(  chip_data, SMRTLinkIDs, SMRTLinkLabels, 'insertions', ylab='rate' )
	plotBoxplot(  chip_data, SMRTLinkIDs, SMRTLinkLabels, 'deletions',  ylab='rate' )
	plotBoxplot(  chip_data, SMRTLinkIDs, SMRTLinkLabels, 'mismatches', ylab='rate' )

	return chip_data

def plotBar( chip_data, SMRTLinkIDs, SMRTLinkLabels, metric, legend_loc=1, width=0.4 ):
	plt.figure()
	plt.hold( True )

	left = 1; left_record = [ left ]
	for chipID in SMRTLinkIDs:
		pyplot.bar( left, chip_data[ chipID ][ metric ], width )
		left += 1
		left_record.append( left )
	left_record = left_record[0:-1]

	tick_locs = [ rec+width/2. for rec in left_record ]
	plt.xticks( tick_locs, SMRTLinkLabels )
	plt.xlim( ( min(tick_locs)-width/2., max(tick_locs)+width/2. ) )
	plt.ylabel( metric )
	plt.title( metric )
	fig = plt.gcf()
	fig.set_size_inches( 10, 6 )

def plotSurvival( chip_data, SMRTLinkIDs, SMRTLinkLabels, metric, legend_loc=1 ):
	plt.figure()
	plt.hold( True )
	for chipID in SMRTLinkIDs:
		sorted_rl = np.sort( chip_data[ chipID ][ metric ] )
		x = sorted_rl[::-1]
		y = np.arange( sorted_rl.size ) / float( len( sorted_rl ) )
		trim = np.max( np.where( y <= 0.005 ) )
		plt.step( x[trim::], y[trim::] )

	plt.grid()
	plt.xlabel( metric )
	plt.ylabel( '1 - ecdf' )
	plt.title( metric + ' survival' )
	plt.legend( SMRTLinkLabels, loc=legend_loc )
	fig = plt.gcf()
	fig.set_size_inches( 10, 6 )

def plotBoxplot( chip_data, SMRTLinkIDs, SMRTLinkLabels, metric, ylab=None, setpositions=None, showmean=True ):
	data = []
	for chipID in SMRTLinkIDs:
		data.append( chip_data[ chipID ][ metric ] )

	ny.coloredBoxplot( data, displabels=SMRTLinkLabels, setpositions=setpositions, dispmeans=showmean )
	if not ylab:
		plt.ylabel( metric )
	else:
		plt.ylabel( ylab )
		plt.title(  metric )
	
	fig = plt.gcf()
	fig.set_size_inches( 10, 6 )
