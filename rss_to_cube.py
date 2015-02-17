#! /usr/bin/env python
#Playing with making fits files that look like VIRUS FOV

#numpy is required by pyfits
import numpy as np
#pyfits to manipulate fits files
import pyfits
#To use things like pi
import math
#To mess with OS paths, use this.
import os
#import ds9 for the sky fiber selection
from ds9 import *
#To make a plot
import pylab as plt
#To fit a curve
from scipy.optimize import curve_fit
#To perform the sigma clipped mean
from astropy.stats import sigma_clip
#to try and find the mode if that's called for
from collections import Counter

plot='plot'

target_skyline = 5577 #Lots of routines will require this, so define it here

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))
def line(x, *p):
	m, b = p
	return m*x+b

def convert_to_pixels(input_angstroms, cdelt, crval):
	return (input_angstroms-crval)/cdelt
def convert_to_angstroms(input_pixels, cdelt, crval):
	return (input_pixels*cdelt)+crval

def sky_shift(input_filename,output_filename):
	#Implement a total shift vs. and individual shift
	individual_shift = 'n'
	if plot == 'plot':
		fig=plt.figure()
		plt.title('Dummy plot to start things off',fontsize=16)
		plt.plot(range(10),range(10))
		plt.ion()
		plt.show()
		plt.clf()
		sleep_time=3
	input_data = pyfits.getdata(input_filename)
	data_header = pyfits.getheader(input_filename, 'DATA')
	cdelt = float(data_header["CDELT1"])
	crval = float(data_header["CRVAL1"])
	crpix = float(data_header["CRPIX1"])
	spectral_resolution = cdelt #I think this is right, but am willing to be corrected
	number_of_fibers = input_data[:,0].size
	number_of_pixels = input_data[0,:].size
	original_peaks = np.zeros(number_of_fibers)
	difference = np.zeros(number_of_fibers)
	target_skyline_pix = convert_to_pixels(target_skyline, cdelt, crval)
	for index, spectrum in enumerate(input_data):
		#print('spectrum '+str(spectrum))
		original_peak = max(spectrum) #I am assuming that the sky line is the brightest pixel in the image
		#print('original_peak '+str(original_peak))
		original_peak_position = spectrum.argmax()
		#print('original_peak_position '+str(original_peak_position))
		p0 = [original_peak, original_peak_position, 1.]
		coeff, gauss_matrix = curve_fit(gauss, np.array(range(number_of_pixels)), spectrum, p0=p0)
		gauss_fit = gauss(np.array(range(number_of_pixels)), *coeff)
		original_peaks[index] = coeff[1]
		#print('coeff '+str(coeff))
		if coeff[1] > 0:
			difference[index] = round(coeff[1] - target_skyline_pix) #not sure which one of these is better, the smart fit
			#print('1st difference[index] '+str(difference[index]))
			difference[index] = round(spectrum.argmax() - target_skyline_pix) #or the dumb easy max find
			#print('2nd difference[index] '+str(difference[index]))
		else:
			difference[index] = 0
		#print('Fitted mean = ', coeff[1])
		if plot == 'plot':
			print("I SHOULD BE PLOTTING NOW")
			plt.title('Observed '+str(target_skyline)+' for each fiber',fontsize=16)
			plt.plot(range(number_of_pixels),spectrum)
			plt.plot(range(number_of_pixels),gauss_fit)
			plt.draw()
			time.sleep(sleep_time)
			plt.clf()
	#print(difference)
	if individual_shift == 'n':
		temp_list = Counter(difference)
		difference_mode = temp_list.most_common(1)
		#print(difference_mode[0][0])
		difference.fill(difference_mode[0][0])
	#print(difference)
	if plot == 'plot':
			plt.title('Observed '+str(target_skyline)+' for each fiber',fontsize=16)
			fiber = np.array(range(number_of_fibers))
			plt.plot(fiber,convert_to_angstroms(original_peaks, cdelt, crval))
			plt.axhline(y=(target_skyline - spectral_resolution))
			plt.axhline(y=(target_skyline + spectral_resolution))
			plt.axhline(y=(target_skyline))
			plt.ylim([target_skyline-(2*max(difference)), target_skyline+(2*max(difference))])
			plt.draw()
			time.sleep(sleep_time)
			plt.clf()
			
	#output_data = np.zeros((number_of_fibers,number_of_pixels))
	output_data = input_data
	#print(output_data.shape)
	for index, spectrum in enumerate(input_data):
		if difference[index] > 0:
			shifter = np.zeros(difference[index])
			temp_spectrum = np.concatenate((spectrum, shifter), axis=0)
			# = spectrum+shifter
			new_spectrum = temp_spectrum[difference[index]:number_of_pixels+difference[index]]
			#print(spectrum - new_spectrum)
			spectrum = new_spectrum
			output_data[index,:] = new_spectrum
		#if difference[index] == 0:
		#	output_data = input_data
			
	pyfits.writeto(output_filename,output_data,data_header,clobber=True)
	

def sky_subtraction(input_filename,output_filename,sky_fibers_file):
	if plot == 'plot':
		fig=plt.figure()
		plt.title('Dummy plot to start things off',fontsize=16)
		plt.plot(range(10),range(10))
		plt.ion()
		plt.show()
		sleep_time=3
	sky_fibers = np.zeros(2)
	if not os.path.isfile(sky_fibers_file):
		d = ds9()
		d.set("file "+input_filename)
		d.set("mode crosshair")
		d.set("scale mode 90")
		d.set("raise")
		d.set("raise")
		dummy = raw_input("Press enter once crosshairs are on the lower sky fiber");
		sky_fibers[0] = int(d.get("crosshair").split()[1])
		print("Lower Sky Fiber: "+str(sky_fibers[0]))
		d.set("raise")
		d.set("raise")
		dummy = raw_input("Press enter once crosshairs are on the upper sky fiber");
		sky_fibers[1] = int(d.get("crosshair").split()[1])
		print("Upper Sky Fiber: "+str(sky_fibers[1]))
		np.savetxt(sky_fibers_file, sky_fibers, fmt='%d')
	#read in the offsets file
	sky_fibers = np.genfromtxt(sky_fibers_file,dtype=None)
	print("Sky Fibers read in: "+str(sky_fibers))
	input_data = pyfits.getdata(input_filename)
	#print(input_filename)
	data_header = pyfits.getheader(input_filename, 'DATA',1)
	gain = float(data_header["GAIN1"]) #probably don't need this, but it was in my IDL, so I'm bringing it over
	cdelt = float(data_header["CDELT1"])
	crval = float(data_header["CRVAL1"])
	crpix = float(data_header["CRPIX1"])
	print("cdelt: "+str(cdelt))
	print("crval: "+str(crval))
	print("crpix: "+str(crpix))
	
	
	
	#input_data is the RSS image data
	
	wave_pix = np.zeros(2)
	number_of_pixels = input_data[0,:].size
	number_of_fibers = input_data[:,0].size
	print('Number of pixels in wavelength axis: '+str(number_of_pixels))
	print('Number of pixels in the fiber axis: '+str(number_of_fibers))
	#find first non-zero value of the array, use the skylines because we've already checked that they're good
	wave_pix[0] = next((i for i, x in enumerate(input_data[sky_fibers[0],:]) if x), None)
	wave_pix[1] = wave_pix[0]+int(0.1*number_of_pixels)
	print("Wave pix: "+str(wave_pix))
	Nrd = np.std(input_data[sky_fibers[0]:sky_fibers[1],wave_pix[0]:wave_pix[1]])
	
	
	
	if plot == 'plot_not':
		plt.title('Sky fibers (also where Nrd is being determined)',fontsize=16)
		fiber = np.array(range(number_of_fibers))
		fiber_intensity = np.array(range(number_of_fibers))
		for index in np.array(range(number_of_fibers)):
			fiber_intensity[index] = sum(input_data[index,:])
		plt.plot(fiber,fiber_intensity)
		plt.axvline(x=sky_fibers[0])
		plt.axvline(x=sky_fibers[1])
		plt.draw()
		time.sleep(sleep_time)
		plt.clf()


	var = np.zeros((number_of_fibers,number_of_pixels))
	var.fill(Nrd)
	
	
	#CALCULATE AND SUBTRACT THE NOISE FROM OUR SIGNAL
	#Calculate the vertical offset, which is just background noise, and subtract that away.
	target_skyline_window = int(40*cdelt) #based on cdelt, we pick a reasonable range  around the window to be fitting around. (~xx pixels for VIMOS HR ~80 pixels for VIMOS LR)
	
	#offset = linear fit to continuium
	#corrected_image = image-offset
	
	#;DETERMINE FLUX NORMALIZATION
	#;Perform a Gaussian fit on the data with the offset removed, and then come up with a normalization factor for each fiber.
	#Subtract the gaussian fit from the spectrum
	
	target_skyline_pixels = convert_to_pixels(target_skyline,cdelt,crval)
	#
	#
	#
	target_skyline_pixels = target_skyline_pixels+15 #!!!!!!!This is a cheat, need to implement the wavelength shifting code
	#
	#
	#
	original_peaks = np.zeros(number_of_fibers)
	integrated_skyline_flux = np.zeros(number_of_fibers)
	for index, spectrum in enumerate(input_data):
		original_peak = max(spectrum) #I am assuming that the sky line is the brightest pixel in the image
		original_peak_position = target_skyline_pixels
		p0 = [original_peak, original_peak_position, 1.]
		#print(original_peak_position)
		#if index == 65:
		#	print(spectrum)
		#	print(original_peak)
		#	
		try:
			coeff, gauss_matrix = curve_fit(gauss, np.array(range(number_of_pixels)), spectrum, p0=p0)
			gauss_fit = gauss(np.array(range(number_of_pixels)), *coeff)
			original_peaks[index] = coeff[1]
			peak = coeff[1] # peak
			sigma = coeff[2] #sigma
		except:
			print("Gaussian fit failed for some reason, defaulting to dumb detection of peak")
			original_peaks[index] = original_peak_position
			peak = original_peak_position
			sigma = cdelt #default to 4x the resolution
		
		p0 = [1.0, 1000]
		#from peak-5sigma to peak-3sigme and peak+3sigma to peak+5sigma
		if sigma > (0.1*number_of_pixels):
			print("Sigma is too high, resetting to 1/100 of wavelength range.")
			sigma = 0.01*number_of_pixels
		if sigma < cdelt/4:
			sigma = cdelt/3
		#print(sigma)
		before_peak = np.array(range(int(round((peak-10*sigma))),int(round((peak-5*sigma)))))
		after_peak = np.array(range(int(round((peak+5*sigma))),int(round((peak+10*sigma)))))
		linear_pixels = np.concatenate((before_peak, after_peak), axis=0)
		#stop()
		#if len(linear_pixels) <= 1:
		#	linear_pixels = 
		#print(linear_pixels)
		#print(len(linear_pixels))
		#coeff, linear_matrix = curve_fit(line, [0,1,2,3,4,5,6,7,8,9], [ 1100.2154541  ,1145.8223877 , 1194.49658203 , 1173.9029541 ,  1254.51281738,  1212.34179688 , 1264.14782715,  1208.92175293 , 1092.0267334  , 1016.29296875], p0=p0)
		if len(linear_pixels) > 1:
			#print(index)
			#print('Length of linear pixels:' + str(len(linear_pixels)))
			#try:
			coeff, linear_matrix = curve_fit(line, np.array(range(len(linear_pixels))), input_data[index,linear_pixels], p0=p0)
			line_fit = line(np.array(range(2*target_skyline_window)), *coeff)
			#except:
			#	coeff = [1.0, 1000]
			#	line_fit = line(range(2*target_skyline_window), *coeff)
			if plot == 'plot1':
				plt.title('Skyline)',fontsize=16)
				plt.plot(np.array(range(2*target_skyline_window)),input_data[index,target_skyline_pixels-target_skyline_window:target_skyline_pixels+target_skyline_window])
				plt.plot(np.array(range(2*target_skyline_window)),line_fit)
				plt.ylim(0,5000)
				plt.draw()
				#time.sleep(sleep_time)
				plt.clf()
			offset = line(np.array(range(number_of_pixels)), *coeff)
			offset_spectrum = spectrum-offset#/2 #not sure if this 2 is required...
			if plot == 'plot1':
				plt.title('Skyline)',fontsize=16)
				plt.plot(range(100), spectrum[790:890])
				plt.plot(range(100), offset_spectrum[790:890])
				plt.ylim(-1000,5000)
				plt.draw()
				plt.clf()
			p0 = [original_peak, original_peak_position, 1.]
			#print(index)
			#
			#Problem here with some fits not converging
			#
			#print(offset_spectrum)
			#fig=plt.figure()
			#plt.plot(range(number_of_pixels),offset_spectrum)
			#plt.show()
			try:
				coeff, gauss_matrix = curve_fit(gauss, np.array(range(number_of_pixels)), offset_spectrum, p0=p0)
				gauss_fit = gauss(np.array(range(number_of_pixels)), *coeff)
			except:
				print("Fit didn't work, defaulting to backup fit.")
				coeff, gauss_matrix = curve_fit(gauss, np.array(range(number_of_pixels)), spectrum, p0=p0)
				gauss_fit = gauss(np.array(range(number_of_pixels)), *coeff)
			if plot == 'plot1':
				plt.title('Skyline)',fontsize=16)
				plt.plot(range(100), spectrum[790:890])
				plt.plot(range(100), gauss_fit[790:890])
				plt.ylim(-1000,5000)
				plt.draw()
				plt.clf()
			#print(coeff[0])
			#print(coeff[2])
			integrated_skyline_flux[index]=coeff[0]*coeff[2]*math.sqrt(2*math.pi) #0 is height, 2 is width
			#print(integrated_skyline_flux)
		#else:
		#	line_fit = [0,0]
	median_skyline_flux = np.median(integrated_skyline_flux)
	#print(median_skyline_flux)
	scale = integrated_skyline_flux/median_skyline_flux
	if plot == 'plot1':
		plt.title('Skyline)',fontsize=16)
		plt.plot(range(number_of_fibers),scale)
		plt.draw()
		time.sleep(sleep_time)
		plt.clf()
	print('Scale is broken, do not use it yet.')
	#for index, spectrum in enumerate(input_data):
	#	#print('Scale: '+str(scale[index]))
	#	if scale[index] > 0.5:
	#		print('Before scale: ',str(input_data[index,10:20]))
	#		input_data[index,:] = input_data[index,:]/scale[index]
	#		print('After scale: ',str(input_data[index,10:20]))
	#		var[index,:] = var[index,:]/scale[index]
	#	else:
	#		#input_data[index,:] = input_data[index,:]*0
	#		var[index,:] = 99999999^2
	average_sky_vector = np.zeros((number_of_pixels))
	temp_storage = np.zeros((sky_fibers[1]-sky_fibers[0]))
	for spectrum_index in np.array(range(number_of_pixels)):
		for fiber_index in np.array(range(sky_fibers[0],sky_fibers[1])):
			#print('fiber_index: '+str(fiber_index))
			temp_storage[fiber_index-sky_fibers[0]] = input_data[fiber_index,spectrum_index]
			#print(input_data[fiber_index,spectrum_index])
		clipped_vector = sigma_clip(temp_storage, 5, 1) #5sigma clipped mean of temp_storage
		#print(temp_storage)
		average_sky_vector[spectrum_index] = np.mean(clipped_vector)
	#print(average_sky_vector)
	if plot == 'plot1':
		plt.title('Average Skyline to be subtract',fontsize=16)
		plt.plot(np.array(range(number_of_pixels)),average_sky_vector)
		plt.ylim(0,1000)
		plt.draw()
		time.sleep(sleep_time)
		plt.clf()
	for index, spectrum in enumerate(input_data):
		if plot == 'plot':
				plt.title('Skyline subtraction in action',fontsize=16)
				plt.plot(np.array(range(number_of_pixels)),input_data[index,:])
				#input_data[index,:] = input_data[index,:]-average_sky_vector
				plt.plot(np.array(range(number_of_pixels)),input_data[index,:]-average_sky_vector)
				plt.plot(np.array(range(number_of_pixels)),average_sky_vector)
				plt.ylim(0,5000)
				plt.draw()
				#time.sleep(sleep_time)
				plt.clf()
		input_data[index,:] = input_data[index,:]-average_sky_vector
	#print("Saving to file: "+output_filename)
	pyfits.writeto(output_filename,input_data,data_header,clobber=True)
	pyfits.append(output_filename, var, data_header)
		
					
	
	

def virus(input_filename, output_filename):
	#create a 400x400 pixel blank field
	#277x267
	radius = 5	#radius of the spaxels in pixels
	pixels_per_fiber = math.pi*radius*radius #this should calculate the "surface area" of the spaxels in pixels
	pixels_per_fiber = 81 #cheat for r=5 because the line above isn't working
	starting_y, starting_x = 5, 15
	increment_y, increment_x = 16, 19
	spaxels_in_x_direction = 15
	spaxels_in_y_direction = 17
	field_size = [(radius*2*spaxels_in_y_direction)+((increment_y-(2*radius))*(spaxels_in_y_direction-1))+1,(radius*2*spaxels_in_x_direction)+((increment_x-(2*radius))*(spaxels_in_x_direction-1))+1] #x and y dimension of field
	
	
	fiber_number = 0

	input_data = pyfits.getdata(input_filename)
	var_data, hdr = pyfits.getdata(input_filename, 1, header=True)
	input_header = pyfits.getheader(input_filename)

	#In order to perserve flux, divide by the number of pixels going in to make each fiber 
	input_data = input_data/pixels_per_fiber

	number_of_fibers = input_data[:,0].size
	number_of_pixels = input_data[0,:].size
	field = np.zeros((number_of_pixels,field_size[0],field_size[1]))
	var_field = np.zeros((number_of_pixels,field_size[0],field_size[1]))
	fov_field = np.zeros((field_size[0],field_size[1]))
	
	

	for y_index in range(17):
		for x_index in range(15):
			if y_index % 2:
				shift = 10
			else:
				shift = 0
			center_y, center_x = starting_y+y_index*increment_y, starting_x-shift+(x_index*increment_x)
			y,x = np.ogrid[-center_y:field_size[0]-center_y, -center_x:field_size[1]-center_x]
			mask = x*x + y*y <= radius*radius
			#spectrum = input_data[fiber_number,:]
			if y_index % 2:
				for spectrum_index in range(number_of_pixels):
					field[spectrum_index,mask] = input_data[fiber_number,spectrum_index]
					var_field[spectrum_index,mask] = var_data[fiber_number,spectrum_index]
					fov_field[mask] = sum(input_data[fiber_number,:])
				fiber_number = fiber_number+1
			else:
				if x_index <= 13:
					for spectrum_index in range(number_of_pixels):
						field[spectrum_index,mask] = input_data[fiber_number,spectrum_index]
						var_field[spectrum_index,mask] = var_data[fiber_number,spectrum_index]
						fov_field[mask] = sum(input_data[fiber_number,:])
					fiber_number = fiber_number+1
	cdelt = float(input_header["CDELT1"])
	crval = float(input_header["CRVAL1"])
	crpix = float(input_header["CRPIX1"])
	input_header["CDELT3"] = cdelt
	input_header["CRVAL3"] = crval
	input_header["CRPIX3"] = crpix
	pyfits.writeto(output_filename,field,input_header,clobber=True)
	pyfits.append(output_filename, var_field, input_header)
	pyfits.writeto(os.path.splitext(output_filename)[0]+"_fov.fits",fov_field,input_header,clobber=True)
plot="noplotting"
#input_filename = '/Users/jimmy/Astro/reduced/mcdonald/20150215/vp_0098_shifted.fits'
#input_filename = '/Users/jimmy/Astro/reduced/mcdonald/20140404/vp_0065_oextr1.fits'
#output_filename = '/Users/jimmy/Astro/reduced/mcdonald/20150215/vp_0098_subbed.fits'
#sky_shift(input_filename, output_filename)
#sky_fibers_file = '/Users/jimmy/Astro/reduced/mcdonald/L109092_2_fibers.txt'
#sky_subtraction(input_filename,output_filename,sky_fibers_file)
#input_filename = '/Users/jimmy/Astro/reduced/mcdonald/20150215/vp_0098_subbed.fits'
#output_filename = '/Users/jimmy/Astro/reduced/mcdonald/20150215/vp_0098_cubefits.fits'
#virus(input_filename, output_filename)
#sky_shift(input_filename)

#Shift the spaxels to line up with 5577
#	fit a Gaussian to the 5577 skyline
#	determine the spectral resolution
#	plot the original peak positions of the 5577 skyline
#	overplot a horizontal line at 5577 to show the target
#	overplot the lines of spectral resolution to see if we are within that limit
#	determine the number of pixels to shift by (round it off)
#	create a new image that's shifted by the appropriate amount in each spaxel
#	write that new image to the old image

#Calculate and subtract the linear offset from the signal
#	select a window where pixels are in the window around 5577
#	calculate the linear offset around the skyline
#	subtract the offset from the image data
#	plot the original image data and the offset to check that it's working properly
	
#Determine the flux normalization
#	fit a gaussian to the 5577 skyline
#	plot the data and the gaussian fit to make sure they line up
#	calculate the integrated flux in the skyline
#	find the median of these flux values
#	divide each sky flux value by the median
#	plot the scale values to check results
#	normalize the transmission for each fiber relative to the median
#	*shifting must occur by this point*
#	Create a sky vector for each fiber#, containing the intensities at that pixel in that region that is skyline, and subtract it from transmission scaled data.
#	stick the sky vectors from each fiber by using a 5 sigma clipped mean (each pixel gets stacked with the 5 sigma clipped mean)
#	store the 5 sigma clipped mean stacked specrum
#	plot the sky vector calculated, and show sky features
#	subtract the mean sky spectrum from all image data
#	plot the original, fit, and corrected data as it's being subtracted
#	reset totally negative fibers to 0, since they were 0 to begin with, and are now subracted, reset them to zero, or just never subtract 0 fibers to begin with
#	(track the sky image spaxels for diagnostic plotting)
#	output the spectra to file
#	(output the diagnostic file)
