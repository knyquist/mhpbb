from pbcore.io.dataset import DataSetIO as io
from pbcore.io.align import BamAlignment as ba
from os import listdir
import numpy as np
import nyplot as ny
import random


def openBAMALN( SMRTLinkID ):
	"""
	Open all bam alignments associated w/ SMRTLink job

	- Will probably need to change to search method once SMRTLink matures and the filesystem spreads out. 
	- primary_path may not remain a single fixed location

	Kristofor Nyquist 12/9/2015
	"""

	if type( SMRTLinkID ) is str and len( SMRTLinkID ) == 6:
		rootdir = SMRTLinkID[0:3];
		linkjob = SMRTLinkID;
		primary_path = '/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/jobs-root/' + rootdir + '/' + linkjob + '/tasks'
		alignment_dirs = getAlignmentDirs( primary_path )

		# generate list of all BAM alignment files
		dsets = []
		for alndir in alignment_dirs:
			contents = listdir( primary_path + '/' + alndir )
			
			flag = 0
			for item in contents:
				if '.alignmentset.xml' in item:
					alnfile = item
					fullpath = primary_path + '/' + alndir + '/' + alnfile
					dsets.append( fullpath )
					flag = 1

			if flag == 0:
				print 'alignmentset.xml not found in ' + primary_path + '/' + alndir

		# generate alignment dataset for entire read
		dset = io.openDataSet( *dsets )

	else:
		print 'SMRTLinkID is either non-existent, or it was not found in path'
		dsets = None

	return dset


def getAlignmentDirs( path ):

	contents = listdir( path )
		
	align_dirs = []
	for item in contents:
		if '.tasks.gather_alignmentset-' in item:
			align_dirs.append( item ) 

	return align_dirs

def convertHoleNumber2XY( holenumber ):
	"""
	Take PacBio holenumber and convert to X and Y coordinate

	Kristofor Nyquist 12/9/2015

	"""

	if holenumber.dtype == 'int32':
		
		mpv = np.iinfo( np.uint16 ).max + 1 # store the max # values [0, 65535] of an unsigned 16 bit integer (65536 values)
		X = holenumber // mpv
		Y = holenumber %  mpv
		return (X, Y)

	else:		
		print 'holenumber is not the accepted int32 type, cannot convert to XY coordinates'
		return None

def subsampleReads( alignmentSet, num ):
	"""
	Take full alignmentSet and randomly subsample to num # samples

	Returns a list of BamAlignments

	Kristofor Nyquist 12/10/2015

	"""

	if isinstance( alignmentSet, io.AlignmentSet ):
		ss = random.sample( alignmentSet, num )
		return ss

	else:
		print 'first input must be an instance of the AlignmentSet class from pbcore'
		return None

def getAlignmentIndicesOfHoleNumbers( alignmentSet, holenumbers, dict=False ):
	"""
	Take full alignmentSet or list of BamAlignments. Return (toggle) dict/list of indices

	corresponding to BamAlignments from specific holenumbers.

	Kristofor Nyquist 12/17/2015

	"""

	# allow AlignmentSet class or a python list of BamAlignments
	proceed = False
	if isinstance( alignmentSet, io.AlignmentSet ): 
		proceed = True

	if proceed or ( type( alignmentSet ) == list and isinstance( alignmentSet[0], ba ) ):
		indices = {}
		for hole in holenumbers:
			indices[ hole ] = list( np.where( alignmentSet.index.holeNumber == hole )[0] )

		# dict flag returns results as dictionary: keys are holenumbers, values are subread indices
		if dict:
			return indices
		
		# dict unflagged returns flattened list of subread indices
		else:
			flattened_indices = []
			for indexlist in indices.values():
				for index in indexlist:
					flattened_indices.append( index )

			return flattened_indices

	else:
		print 'first input must be an instance of the AlignmentSet class from pbcore or a list of BamAlignments from pbcore'
		return None

def subsampleReadsByIndices( alignmentSet, indices ):
	"""
	Take full alignmentSet and subsample by list of indices

	Returns a list of BamAlignments

	Kristofor Nyquist 12/17/2015

	"""

	proceed = False
	if isinstance( alignmentSet, io.AlignmentSet ) and type( indices ) == list:
		proceed = True

	elif proceed or ( type( alignmentSet ) == list and isinstance( alignmentSet[0], ba ) and type( indices ) == list ): 
		ss = []
		for index in indices:
			ss.append( alignmentSet[ index ] )
		return ss

	else:
		print 'first input must be an instance of the AlignmentSet class from pbcore or a list of BamAlignments from pbcore'

def getChipKPIs( SMRTLinkIDs, num=1000 ):
	"""
	Return key performance indicators from chip(s). 
		<SMRTLinkIDs> contains SMRTLinkIDs for each chip. Can be a list of SMRTLink JobIDs or a single SMRTLink JobID
		<num>         is the number of alignments to subsample to. Default is 1000 alignments. 

	Returns a dictionary of dictionary data structure. First layer of keys is each SMRTLink JobID, representing
	individual chips. Second layer of keys corresponds to each KPI for each chip.
		chips = { 'jobID_1': dict1,
				  'jobID_2': dict2  }

		chips[ 'jobID_1' ] = { 'readlength'  : list, 
							   'templatespan': list,
							   'insertions'  : list,
							   'deletions'   : list,
							   'mismatches'  : list,
							   'accuracy'    : list,
							   'IPD'         : list }

	Kristofor Nyquist 1/26/2016
	"""

	chips = {}								# chips dict will store KPIs for each chip
	for jobID in SMRTLinkIDs:
		a = openBAMALN( jobID )            	# open alignment reader
		if len( a ) > num:					# if # alignments > num,
			ss = subsampleReads( a, num )	# randomly subsample

		chips[ jobID ] = getAlnBamKPIs( a, ss )

	return chips

def getHoleNumberKPIs( SMRTLinkID, holenumbers ):
	"""
	Same concept as getChipKPIs except this time it grabs KPIs for specific set of holes
		<SMRTLinkID>  SMRTLinkID for chip
		<holenumbers> list of holenumbers

	Kristofor Nyquist 1/27/2016
	"""

	a = openBAMALN( SMRTLinkID )
	zmw_alnix_map = getAlignmentIndicesOfHoleNumbers( a, holenumbers, dict=True )
	
	zmws = {}
	for zmw in zmw_alnix_map:

		ss   = []
		for index in zmw_alnix_map[ zmw ]: 
			ss.append( a[ index ] )
			
		zmws[ index ] = getAlnBamKPIs( a, ss )

	return zmws

def getAlnBamKPIs( alignmentSet, subsampledSet ):
	"""
	Retrieve the KPIs for a single chip in a dictionary structure.
		<alignmentSet> is a pbcore alignmentSet reader object. Likely subsampled using routine subsampleReads(...).

	Kristofor Nyquist 1/26/2016

	"""
	chip                   = {}
	chip[ 'readlength' ]   = []
	chip[ 'templatespan' ] = []
	chip[ 'insertions' ]   = []
	chip[ 'deletions' ]    = []
	chip[ 'mismatches' ]   = []
	chip[ 'accuracy' ]     = []
	chip[ 'IPD' ]          = []

	for aln in subsampledSet:
		chip[ 'readlength'   ].append( float( aln.readEnd - aln.readStart ) )
		chip[ 'templatespan' ].append( float( aln.referenceEnd - aln.referenceStart ) )
		chip[ 'insertions'   ].append( float( aln.nIns ) / chip[ 'readlength' ][-1] )
		chip[ 'deletions'    ].append( float( aln.nDel ) / chip[ 'readlength' ][-1] )
		chip[ 'mismatches'   ].append( float( aln.nMM )  / chip[ 'readlength' ][-1] )
		error_rate =   ( aln.nIns + aln.nDel + aln.nMM ) / chip[ 'readlength' ][-1]
		chip[ 'accuracy'     ].append( 1 - error_rate )
		chip[ 'IPD'          ].append( aln.IPD() )

	chip[ 'polreadlength' ] = reconstructPolReadlength( alignmentSet, subsampledSet )			# add polymerase read length as KPI

	return chip

def reconstructPolReadlength( alignmentSet, subsampledSet ):
	"""
	Estimate the full polymerase readlength from the subsampled holenumber set.
		<alignmentSet>  is the full alignmentset reader object from pbcore
		<subsampledSet> is the subsampled dataset reader used for estimating the other KPIs

	Kristofor Nyquist 1/27/2016
	"""
	
	polreadlength = []
	for aln in subsampledSet:
		hn = aln.holeNumber
		polreadlength.append( max( alignmentSet.index[ alignmentSet.index.holeNumber == hn ].aEnd ) - min( alignmentSet.index[ alignmentSet.index.holeNumber == hn ].aStart ) )

	return polreadlength


	

