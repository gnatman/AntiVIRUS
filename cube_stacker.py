#! /usr/bin/env python
#Stack up data cubes
#This code may or may not work, and has definitely not been proven to be scientifically rigorous, I started it and then had to abandon the project, so I intend to finish it at some point.

#numpy is required by pyfits
import numpy as np
#pyfits to manipulate fits files
import pyfits
#import glob to search for the right file names
import glob
#Use ds9 to align the cubes
from ds9 import *

#read in fits files, shift second image by integer, sum images, write to file
#output_dir = "/Users/jimmy/Astro/reduced/mcdonald/BCG_1261/"
#input_dir = "/Users/jimmy/Astro/reduced/mcdonald/BCG_1261/"

#offsets_file = "/Users/jimmy/Astro/reduced/mcdonald/BCG_1261/offsets.txt"


def stack(input_dir, target, offsets_file):
	input_fovs = glob.glob(input_dir+target+'*_fov.fits') #used the fov files because they're stacked
	input_cubes = [cube[0:-9]+'.fits' for cube in input_fovs]
	print("Stacking the following files into one cube: ")
	print(input_cubes)
	number_of_cubes = len(input_cubes)
	output_dir = input_dir
	
	
	if os.path.isfile(offsets_file):
		#read in the offsets file
		shift = np.genfromtxt(offsets_file,dtype=None)
	else:
		#make one using ds9
		d = ds9()
		for index, cube in enumerate(input_fovs):
			d.set("frame "+str(index+1))
			d.set("file "+cube)
			#d.set("datacube 29")
			d.set("mode crosshair")
	
		dummy = raw_input("Press enter once crosshairs are in the correct location");
	
		coords = np.zeros((number_of_cubes,2))
		shift = np.zeros((number_of_cubes,2))
	
		for frame in range(number_of_cubes):
			d.set("frame "+str(frame+1))
			coords_str = d.get("crosshair")
			coords[frame,:] = [int(n) for n in coords_str.split()]
		
		for frame in range(number_of_cubes):
			shift[frame] = coords[0] - coords[frame]
	
	min_y_shift = min(shift[:,1])
	print("min_y_shift: "+str(min_y_shift))
	min_x_shift = min(shift[:,0])
	print("min_x_shift: "+str(min_x_shift))
	if min_y_shift < 0:
		shift[:,1] = shift[:,1]-min_y_shift
	if min_x_shift < 0:
		shift[:,0] = shift[:,0]-min_x_shift
	print("Shifting by the following values: \n"+str(shift))
	np.savetxt(offsets_file, shift, fmt='%d')
	
	#print(input_cubes[0])
	input_data = pyfits.getdata(input_cubes[0])
	var_data, hdr = pyfits.getdata(input_cubes[0], 1, header=True) #read from 2nd header
	spectrum_size = int(input_data[:,0,0].size)
	y_size = int(input_data[0,:,0].size)
	x_size = int(input_data[0,0,:].size)
	max_y_shift = max(shift[:,1])
	max_x_shift = max(shift[:,0])
	output_data = np.zeros((spectrum_size,y_size+max_y_shift,x_size+max_x_shift))
	var_output_data = np.zeros((spectrum_size,y_size+max_y_shift,x_size+max_x_shift))
	stacker_tracker = np.zeros((y_size+max_y_shift,x_size+max_x_shift))
	
	for index, cube in enumerate(input_cubes):
		input_data = pyfits.getdata(cube)
		var_data = pyfits.getdata(cube, 1)
		data_header = pyfits.getheader(cube)
		shift_x = int(shift[index,0])
		shift_y = int(shift[index,1])
	
		for x in range(x_size):
			for y in range(y_size):
				output_data[:,y+shift_y,x+shift_x] = output_data[:,y+shift_y,x+shift_x] + input_data[:,y,x]
				var_output_data[:,y+shift_y,x+shift_x] = var_output_data[:,y+shift_y,x+shift_x] + var_data[:,y,x]
				if not sum(input_data[:,y,x]) == 0:
					stacker_tracker[y+shift_y,x+shift_x] = stacker_tracker[y+shift_y,x+shift_x] + 1
	
	pyfits.writeto(output_dir+'combined_cube.fits',output_data,data_header,clobber=True)
	pyfits.append(output_dir+'combined_cube.fits',var_output_data,data_header)
	pyfits.writeto(output_dir+'stacker_tracker.fits',stacker_tracker,data_header,clobber=True)