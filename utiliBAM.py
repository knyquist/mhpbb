from pbcore.io.dataset import DataSetIO as io
from pbcore.io.align import BamAlignment as ba
from pbcore.io import BamIO
from os import listdir
import xml.etree.ElementTree as ET
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

def getPrimaryFolder( SMRTLinkID ):

	rootdir  = SMRTLinkID[0:3];
	linkjob  = SMRTLinkID;
	filename1 = '/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/jobs-root/' + rootdir + '/' + linkjob + '/tasks/pbsmrtpipe.tasks.subreadset_align_scatter-1/chunk_subreadset_0.subreadset.xml'
	filename2 = '/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/jobs-root/' + rootdir + '/' + linkjob + '/tasks/pbcoretools.tasks.subreadset_align_scatter-1/chunk_subreadset_0.subreadset.xml'

	flag = 0
	try:
		xml_tree = ET.parse( filename1 )
		flag = 1
	except IOError:
		try:
			xml_tree = ET.parse( filename2 )
			flag = 1
		except:
			print 'Could not find primary location, failing gracefully w/o primary location' # did not find xml file that points to primary folder
			return None
	if flag == 1:
		root = xml_tree.getroot()
		for layer1 in root:
			if 'ExternalResources' in layer1.tag:
				for layer2 in layer1:
					if 'ExternalResource' in layer2.tag:
						for tup in layer2.items():
							if tup[0] == 'ResourceId':
								folder = tup[1].split('/')
								folder = folder[0:-1]
								folder = '/'.join( folder )
								return folder # found primary folder

		print 'Found xml pointer file, but was unable to find primary folder path'
		return None

def getTotalNumberOfAvailableZMWs( SMRTLinkID ):
	try:
		primary_path = getPrimaryFolder( SMRTLinkID )
	except:
		print 'could not find primary location, failing gracefully'
	else:
		if not primary_path:
			print 'Error: getPrimaryFolder() returned None'
			return None
		else:
			try: # to identify the subreads.bam file
				primary_files = listdir( primary_path )
				subread_files = [ f for f in primary_files if 'subreads.bam' in f ]
				subread_bam = min( subread_files, key=len )
				full_path = primary_path + '/' + subread_bam
			except:
				print 'could not find the subreads.bam file, failing gracefully'
				return None
			else:
				indexed_subread_reader = BamIO.IndexedBamReader( full_path )
				return len( indexed_subread_reader.index )

def getReferenceSequence( SMRTLinkID ):
	"""
	Find fasta file of reference sequence and return sequence as string

	Kristofor Nyquist 2/12/2016
	"""
	if type( SMRTLinkID ) is str and len( SMRTLinkID ) == 6:
		rootdir = SMRTLinkID[0:3];
		linkjob = SMRTLinkID;
		primary_path = '/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/jobs-root/' + rootdir + '/' + linkjob + '/tasks'
		xml_tree = ET.parse( primary_path + '/pbsmrtpipe.tasks.gather_alignmentset-1/file.alignmentset.xml' )
		root = xml_tree.getroot()

		# parse xml file to get path to reference fasta. This is hacky code, but gets the job done
		flag = 0
		for layer1 in root:
		    if 'ExternalResources' in layer1.tag and flag == 0:
		        for layer2 in layer1:
		            if 'ExternalResource' in layer2.tag and flag == 0:
		                for layer3 in layer2:
		                    if 'ExternalResources' in layer3.tag and flag == 0:
		                        for layer4 in layer3:
		                            if 'ExternalResource' in layer4.tag and flag == 0:
		                                for tup in layer4.items():
		                                    if tup[0] == 'ResourceId':
		                                        fastaPath = tup[1]
		                                        flag = 1
		
		# create list of the reference sequence
		ref_file = open( fastaPath, 'rb' )
		reference = ''
		ref_file.next() # skip the header
		for line in ref_file:
			reference = reference + line[0:-1] # skip \n symbol
		return reference

def calculateCoverageByGC( SMRTLinkID ):
	"""
	Use pbcore's _intervalContour to retrieve the number of reads covering each base in the reference sequence.
	Use utiliBAM's getReferenceSequence to retrieve the reference sequence.

	Calculate the coverage by GC percentage and return for convenient boxplotting

	Kristofor Nyquist 2/12/2016
	"""

	chip_data = openBAMALN( SMRTLinkID )
	print 'calculating coverage for ' + chip_data.refNames[0]
	coverage  = chip_data._intervalContour( chip_data.refNames[0] )
	coverage  = coverage[0:-1]
	reference = getReferenceSequence( SMRTLinkID )
	if len( coverage ) == len( reference ): # proceed
		start_ix = 12						# use a 25 bp window to calculate % GC
		end_ix   = len( reference ) - 13
		gc  = []
		cov = []
		for i in range( start_ix, end_ix+1 ): # +1 because right boundary is exclusive
			cov.append( coverage[i] )
			seq = reference[ i-12:i+13 ].upper()
			gc_cont = [sub for sub in seq if sub == 'G' or sub == 'C']
			prct_gc = float( len( gc_cont ) ) / len( seq )
			gc.append( prct_gc )

		# generate boxplot for each possible value
		coverage = []
		labels   = []
		for i in range( 0, 25+1 ): # +1 because right boundary is exclusive
		    i = float(i)
		    ixs = [ ix for ix,val in enumerate( gc ) if val == i/25 ]
		    
		    bin_coverage = []
		    for index in ixs:
		        bin_coverage.append( cov[index] )
		    coverage.append( bin_coverage )
		    labels.append( str(i/25)[1:] )
		labels[0] = 0
		labels[-1] = 1

		return ( coverage, labels )

	else:
		print 'Coverage map has a different length than the reference sequence. There is a problem'

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
		if len( alignmentSet ) > num:					# make sure there are enough alignments to subsample
			ss = random.sample( alignmentSet, num )
		else:
			ss = alignmentSet 							# if there are less alignments than desired subsample ss = a
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
		a  = openBAMALN( jobID )            	# open alignment reader
		ss = subsampleReads( a, num )			# randomly subsample
		chips[ jobID ] = getAlnBamKPIs( a, ss, jobID )

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

def getAlnBamKPIs( alignmentSet, subsampledSet, SMRTLinkID ):
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
		chip[ 'holenumber'   ] = aln.HoleNumber
		chip[ 'readlength'   ].append( float( aln.readEnd - aln.readStart ) )
		chip[ 'templatespan' ].append( float( aln.referenceEnd - aln.referenceStart ) )
		chip[ 'insertions'   ].append( float( aln.nIns ) / chip[ 'readlength' ][-1] )
		chip[ 'deletions'    ].append( float( aln.nDel ) / chip[ 'readlength' ][-1] )
		chip[ 'mismatches'   ].append( float( aln.nMM )  / chip[ 'readlength' ][-1] )
		error_rate =   ( aln.nIns + aln.nDel + aln.nMM ) / chip[ 'readlength' ][-1]
		chip[ 'accuracy'     ].append( 1 - error_rate )
		chip[ 'IPD'          ].append( aln.IPD() )

	chip[ 'polreadlength' ] = reconstructPolReadlength( alignmentSet, subsampledSet )			# add polymerase read length as KPI
	chip[ 'total nreads' ]  = len( alignmentSet.index )
	try: # getting the total number of ZMWs from primary
		 # I do a try, because if the file system changes, this is going to fail and I want it to fail gracefully
		total_ZMWs = getTotalNumberOfAvailableZMWs( SMRTLinkID )
	except:
		print 'getTotalNumberofAvailableZMWs failed, cannot calculate loading percentage'
	else:
		if not total_ZMWs:
			print 'getTotalNumberofAvailableZMWs returned None, cannot calculate loading percentage'
		else:
			chip[ 'fraction loaded' ] = float(chip[ 'total nreads' ])/total_ZMWs

	return chip

def reconstructPolReadlength( alignmentSet, subsampledSet ):
	"""
	Estimate the full polymerase readlength from the subsampled holenumber set.
		<alignmentSet>  is the full alignmentset reader object from pbcore
		<subsampledSet> is the subsampled dataset reader used for estimating the other KPIs

	Kristofor Nyquist 1/27/2016
		updated 3/14/2016 for aln_index exception catching
	"""
	
	polreadlength = []
	for aln in subsampledSet:
		hn = aln.holeNumber
		try:
			max_aln_index = max( alignmentSet.index[ alignmentSet.index.holeNumber == hn ].aEnd )
			min_aln_index = min( alignmentSet.index[ alignmentSet.index.holeNumber == hn ].aStart )
		except:
			print 'holenumber did not have alignment end or alignment start indices.'
		else:
			polreadlength.append( max_aln_index - min_aln_index )

	return polreadlength


	

