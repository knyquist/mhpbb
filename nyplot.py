import numpy as np
import matplotlib.pylab as pylab
pylab.rcParams['figure.figsize'] = 20, 15
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = 20, 15
plt.rcParams.update({'font.size': 20})

# TBD
#def groupedBoxplot(  ):

def coloredBoxplot( data, setpositions=None, displabels=None, dispmeans=False, dispoutliers=None ):
	"""
	Generate a boxplot figure with multiple boxplots, each a different color.

		<data>         is a list of lists. Each list forms the data for one of the boxes
		<displabels>   is a list of strings. Each string labels one of the boxes
		<dispmeans>    is a boolean that toggles whether to show mean of each distribution
		<dispoutliers> is a boolean that toggles whether to show the outliers of each distribution

	Kristofor Nyquist 1/25/2016
	
	FOR THE FUTURE: if you wanted to do multiple chips, with each chip having multiple 
	templates (like in milhouse), you would use this function to create a boxplot for 
	each chip, making sure to set hold=True for the figure and to set the positions 
	so there was a gap between chips. 

	see http://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
	"""
	
	plt.figure()
	bp = plt.boxplot( data, positions=setpositions, labels=displabels, showmeans=dispmeans, showfliers=dispoutliers )
	setBoxColors( bp, len(data) )
	plt.grid()
	return bp

def setBoxColors( h, nsets ):
	"""
	Secondary function that sets the boxplot figure colorscheme
		<h>     is the handle to the boxplot figure
		<nsets> is the number of datasets, i.e. the number of different colors required

	Kristofor Nyquist 1/25/2016
	"""

	colors     = [ 'k', 'y', 'm', 'c', 'r', 'g', 'b' ]
	color_dump = [] # will recycle colors if we need more than seven boxes in our figure

	iterable = 0
	for i in range( 0, nsets ):
		c = colors.pop()								# pop out a color to use it, 
		color_dump.append( c )							# but append it to the dump list
		if not colors:	        						# if all colors have been popped out, 
			colors     = list( reversed( color_dump ) )	# we reinitialize it from the color_dump
			color_dump = []

		plt.setp( h[ 'boxes'    ][ i        ], color = c )
		plt.setp( h[ 'caps'     ][ iterable ], color = c ); iterable += 1
		plt.setp( h[ 'caps'     ][ iterable ], color = c ); iterable -= 1
		plt.setp( h[ 'whiskers' ][ iterable ], color = c ); iterable += 1
		plt.setp( h[ 'whiskers' ][ iterable ], color = c ); iterable += 1
		plt.setp( h[ 'medians'  ][ i        ], color = c )


