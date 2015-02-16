#! /usr/bin/env python
#A tool to automatically reduce VIRUS-P data for a given night.

#pyfits to read in fits files
import pyfits
#numpy is required by pyfits
import numpy as np
#sys to read in the arguments
import sys
#Not sure why this is imported, maybe I can just use this to read the tildas and not use the next import line
import os
#To expand the tilda and define the home directory, also to call other python scripts from within this one
from os.path import expanduser
#Import shell utilities for file copy
import shutil
#Import time to slow things down as we're executing external calls
import time
#Import rss_to_cube to turn the row stacked spectra into a datacube
import rss_to_cube
#cube_stacker will combine cubes into one
import cube_stacker

#Set important directories and filenames
HOME = expanduser("~")
RAW_DATA_DIR = HOME+"/Astro/reduced/mcdonald/raw/"
OUTPUT_HEAD_DIR = HOME+"/Astro/reduced/mcdonald/"
database_file_name = "database.txt"
instrument_parfile = HOME+"/Astro/idl/p3d/data/instruments/virus/virus.prm"
linelist = HOME+"/Astro/idl/p3d/data/tables/linelists/necd.dat"
user_parfile = HOME+"/Astro/p3d/user_p3d.prm"
#arc_elements = "NeCd"
sleep_time = 10

#First argument is the night that the data was taken, in the form: YYYYMMDD
date = sys.argv[1]

OUTPUT_DIR = OUTPUT_HEAD_DIR+date
INPUT_DIR = RAW_DATA_DIR+date
if not os.path.isdir(OUTPUT_DIR):
	print("Directory "+OUTPUT_DIR+" doesn't exist, making it now.")
	os.makedirs(OUTPUT_DIR)

#Keywords I have to work with:
#UT, the time in UT, which may be good for solving the problem above

#If a file does not already exist with the table of information needed, create it.
#We make this tracker file so we don't have to go through and reclassify the whole list of 2000+ files every time
if not os.path.isfile(OUTPUT_DIR+"/"+database_file_name): 
	table_of_attributes = open(OUTPUT_DIR+"/"+database_file_name, "w")
	print("Database does not exist, creating it now. (This may take a bit.)")
	for fits_file in os.listdir(INPUT_DIR):
		if fits_file.endswith(".fits"):
			galaxy_hdu = pyfits.open(INPUT_DIR+"/"+fits_file)
			if str(galaxy_hdu[0].header["OBSERVAT"]) == "MCDONALD":
				table_of_attributes.write(INPUT_DIR+"/"+fits_file+"\t") #First Column, full file name
				table_of_attributes.write(fits_file.replace(" ", "")+"\t") #Second Column, just the file name, without directory
				table_of_attributes.write(str(galaxy_hdu[0].header["IMAGETYP"]).replace(" ", "_")+"\t") #Third column, Image Type
				table_of_attributes.write(galaxy_hdu[0].header["OBJECT"].replace(" ", "_")+"\t") #Fourth Column, object name
				table_of_attributes.write(str(galaxy_hdu[0].header["COMMENTS"]).replace(" ", "_")+"\t") #Fourth column, comments
				table_of_attributes.write(str(galaxy_hdu[0].header["DATE-OBS"])[0:4]+"\t") #Fifth column, year
				table_of_attributes.write(str(galaxy_hdu[0].header["DATE-OBS"])[5:7]+"\t") #Sixth column, month
				table_of_attributes.write(str(galaxy_hdu[0].header["DATE-OBS"])[8:10]+"\t") #Seventh column, day
				table_of_attributes.write(str(galaxy_hdu[0].header["UT"])+"\n") #Eight column, UT Time	 
	galaxy_hdu.close()
	table_of_attributes.close()

#Make the Master Bias Frame
#Read the table back in, just so we're consistent if it wasn't made in the previous step.
datatable=np.loadtxt(OUTPUT_DIR+"/"+database_file_name, dtype=str)
bias_indices = np.where(datatable[:,2]=="zero")[0]
table_of_biases = datatable[bias_indices]
master_bias_file = OUTPUT_DIR+"/"+os.path.splitext(table_of_biases[0,1])[0]+"_mbias1.fits"
if not os.path.isfile(master_bias_file):
	bias_files = " "
	for filename in table_of_biases[:,0]:
		bias_files = bias_files+","+filename
	#print("p3d_dispatch.py p3d_cmbias "+bias_files+" "+instrument_parfile+" outvar userparfile="+user_parfile+" opath="+OUTPUT_DIR)
	os.system("p3d_dispatch.py p3d_cmbias "+bias_files+" "+instrument_parfile+" outvar userparfile="+user_parfile+" opath="+OUTPUT_DIR)
	while not os.path.isfile(master_bias_file): 
		time.sleep(1)
		print("Waiting for Master Bias creation to finish...")
	print("Master Bias Created: "+master_bias_file)

#Make the Trace File
#Read the table back in, just so we're consistent if it wasn't made in the previous step.
trace_indices = np.where(datatable[:,2]=="flat")[0]
table_of_traces = datatable[trace_indices]
trace_mask_file = OUTPUT_DIR+"/"+os.path.splitext(table_of_traces[0,1])[0]+"_imcmb1_trace1.fits"
if not os.path.isfile(trace_mask_file):
	trace_files = " "
	for filename in table_of_traces[:,0]:
		trace_files = trace_files+","+filename
	os.system("p3d_dispatch.py p3d_ctrace "+trace_files+" "+instrument_parfile+" outvar masterbias="+master_bias_file+" userparfile="+user_parfile+" opath="+OUTPUT_DIR)
	while not os.path.isfile(trace_mask_file): 
		time.sleep(1)
		print("Waiting for Trace Mask creation to finish...")
	print("Trace Mask Created: "+trace_mask_file)

#Create Dispersion solution
#Read the table back in, just so we're consistent if it wasn't made in the previous step.
arc_indices = np.where(datatable[:,2]=="comp")[0]
table_of_arcs = datatable[arc_indices]
disp_mask_file = OUTPUT_DIR+"/"+os.path.splitext(table_of_arcs[0,1])[0]+"_imcmb1_dmask1.fits"
if not os.path.isfile(disp_mask_file):
	arc_files = " "
	for filename in table_of_arcs[:,0]:
		arc_files = arc_files+","+filename
	os.system("p3d_dispatch.py p3d_cdmask "+arc_files+" "+instrument_parfile+" savedvar outvar masterbias="+master_bias_file+" tracemask="+trace_mask_file+" userparfile="+user_parfile+" opath="+OUTPUT_DIR+" arclinefile="+linelist)
	while not os.path.isfile(disp_mask_file): 
		time.sleep(1)
		print("Waiting for Dispersion Mask creation to finish...")
	print("Dispersion Mask Created: "+disp_mask_file)
	
#Create a Fiber Flat
#Read the table back in, just so we're consistent if it wasn't made in the previous step.
flat_indices = np.where(datatable[:,2]=="flat")[0]
table_of_flats = datatable[flat_indices]
flat_mask_file = OUTPUT_DIR+"/"+os.path.splitext(table_of_flats[0,1])[0]+"_imcmb1_flatf1.fits"
if not os.path.isfile(flat_mask_file):
	flat_files = " "
	for filename in table_of_flats[:,0]:
		flat_files = flat_files+","+filename
	os.system("p3d_dispatch.py p3d_cflatf "+flat_files+" "+instrument_parfile+" outvar masterbias="+master_bias_file+" tracemask="+trace_mask_file+" dispmask="+disp_mask_file+" userparfile="+user_parfile+" opath="+OUTPUT_DIR)
	while not os.path.isfile(flat_mask_file): 
		time.sleep(1)
		print("Waiting for Fiber Flats creation to finish...")
	print("Fiber Flats Created: "+flat_mask_file)

#Reduce the science data
#Read the table back in, just so we're consistent if it wasn't made in the previous step.
object_indices = np.where(datatable[:,2]=="object")[0]
table_of_objects = datatable[object_indices]
for index, input_object in enumerate(table_of_objects[:,0]):
	object_file = OUTPUT_DIR+"/"+os.path.splitext(table_of_objects[index,1])[0]+"_oextr1.fits"
	CUBE_OUTPUT_DIR = OUTPUT_HEAD_DIR+table_of_objects[index,3][0:-9]+"/"
	if not os.path.isdir(CUBE_OUTPUT_DIR):
		print("Directory "+CUBE_OUTPUT_DIR+" doesn't exist, making it now.")
		os.makedirs(CUBE_OUTPUT_DIR)
	data_cube = CUBE_OUTPUT_DIR+table_of_objects[index,3]+"_"+str(index)+".fits"
	if not os.path.isfile(object_file):
		os.system("p3d_dispatch.py p3d_cobjex "+input_object+" "+instrument_parfile+" outvar bswvar masterbias="+master_bias_file+" tracemask="+trace_mask_file+" dispmask="+disp_mask_file+" flatfield="+flat_mask_file+" userparfile="+user_parfile+" opath="+OUTPUT_DIR)
		while not os.path.isfile(object_file): 
			time.sleep(1)
			print("Waiting for Object data reduction "+table_of_objects[index,1]+" to finish...")
		time.sleep(10)
		print("Science Object Reduced: "+object_file)
		print("Shifting Skylines.")
		shifted_file = OUTPUT_DIR+"/"+os.path.splitext(table_of_objects[index,1])[0]+"_shifted.fits"
		rss_to_cube.sky_shift(object_file, shifted_file)
		print("Subtracting the sky from file: "+shifted_file)
		sky_fibers_file = CUBE_OUTPUT_DIR+table_of_objects[index,3]+"_"+str(index)+"_fibers.txt"
		print("Using sky fibers list: "+sky_fibers_file)
		subbed_object_file = OUTPUT_DIR+"/"+os.path.splitext(table_of_objects[index,1])[0]+"_subbed.fits"
		print("And saving to file: "+subbed_object_file)
		rss_to_cube.sky_subtraction(shifted_file,subbed_object_file,sky_fibers_file)
		print("Making data cube: "+data_cube)
		rss_to_cube.virus(subbed_object_file,data_cube)

combine_question = raw_input("Combine cubes now?");
if combine_question == 'y':
	#once all the individual cubes are made, stack them 
	temp_list = table_of_objects[:,3]
	for index, input_object in enumerate(table_of_objects[:,3]):
		temp_list[index] = input_object[0:-9]
	object_list = list(set(temp_list))
	for index, target in enumerate(object_list):
		input_dir = OUTPUT_HEAD_DIR+target+"/"
		offsets_file = OUTPUT_HEAD_DIR+target+"/offsets.txt"
		if not os.path.isfile(input_dir+"combined_cube.fits"):
			print("Creating cube: "+input_dir+"combined_cube.fits")
			cube_stacker.stack(input_dir, target, offsets_file)
