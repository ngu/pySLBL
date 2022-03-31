# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# Raster2DAN3D.py
# Created on: 09.01.2018
# by: Pierrick.Nicolet@ngu.no
# ---------------------------------------------------------------------------
#
# This tool intends to create the input data required by DAN3D (se references
# below). It converts the DEM produced by the SLBL tool as a “Path topography”
# and the thickness file (optional in the SLBL tool) as the “Source depth”.
# Optionally, an “Erosion map” can be created. To date, it attributes the
# material 2 to all the cells with values above 0 in the input grid and the
# material 1 to all the cells where the value is 0. It is therefore intended to
# be used with the thickness raster of the SLBL tool (applied to loose material
# along the landslide path), but a more flexible solution should be implemented
# in a future release. All grids are saved in the defined folder as surfer
# grids, a format that can be read but not written by ArcMap. The header is
# modified to shift all the grids to (0,0), which is advised for DAN3D. An
# additional text file is saved with the original header that can be manually
# copy/pasted in the resulting grids (replaces the 4 first lines).
#
# References:
#
# McDougall, S. & Hungr, O.
# A model for the analysis of rapid landslide motion across three-dimensional terrain 
# Canadian Geotechnical Journal, 2004, 41, 1084-1097
#
# McDougall, S. & Hungr, O.
# Dynamic modelling of entrainment in rapid landslides 
# Canadian Geotechnical Journal, 2005, 42, 1437-1448 
#
# ---------------------------------------------------------------------------

import arcpy
import numpy as np
import os

def write_surfgrd(npgrid, filepath, cellSize):
	with open(filepath,'w') as file:
		file.write('DSAA\n')
		file.write('{0} {1}\n'.format(npgrid.shape[1], npgrid.shape[0])) #columns and rows
		file.write('0 {0}\n'.format(int(npgrid.shape[1]*cellSize))) #y coordinate shifted to 0
		file.write('0 {0}\n'.format(int(npgrid.shape[0]*cellSize))) #x coordinate shifted to 0
		file.write('{0} {1}\n'.format(np.amin(npgrid), np.amax(npgrid)))
		for i in np.flipud(range(0,npgrid.shape[0])):
			nb = 1
			for j in range(0,npgrid.shape[1]):
				if j == npgrid.shape[1]-1: #last line
					file.write("{0}\n\n".format(npgrid[i,j]))
				elif nb%10 == 0: # mutiple of ten
					file.write("{0}\n".format(npgrid[i,j]))
				else:
					file.write("{0} ".format(npgrid[i,j]))
				nb += 1

if __name__=="__main__":

	# Retrive the parameters from the GUI
	SlblFile = arcpy.GetParameterAsText(0)
	DepthFile = arcpy.GetParameterAsText(1)
	ErosionFile = arcpy.GetParameterAsText(2)
	outFolder = arcpy.GetParameterAsText(3)
	cellFactor = float(arcpy.GetParameterAsText(4))
	
	str_message = 'Saving DAN3D input files...'
	arcpy.AddMessage(str_message)
	
	SlblFile = arcpy.Raster(SlblFile)
	DepthFile = arcpy.Raster(DepthFile)
	if len(ErosionFile) > 0:
		Erosion = True
		ErosionFile = arcpy.Raster(ErosionFile)
	else:
		Erosion = False
	
	IncellSize = SlblFile.meanCellWidth
	if cellFactor != 1:
		cellSize = IncellSize * cellFactor
		str_message = "Inputs rasters are being aggreagated. Output cell size will be {}m".format(cellSize)
		arcpy.AddMessage(str_message)
		SlblFile = arcpy.sa.Aggregate(SlblFile, cellFactor, "MEDIAN" ,"EXPAND", "DATA")
		DepthFile = arcpy.sa.Aggregate(DepthFile, cellFactor, "MEDIAN" ,"EXPAND", "DATA")
		if Erosion == True:
			ErosionFile = arcpy.sa.Aggregate(ErosionFile, cellFactor, "MEDIAN" ,"EXPAND", "DATA")
	else:
		cellSize = IncellSize

	# The extent is moved by half a cell (cell center instead of border)
	xmin = SlblFile.extent.XMin + (cellSize/2.0)
	ymin = SlblFile.extent.YMin + (cellSize/2.0)
	xmax = SlblFile.extent.XMax - (cellSize/2.0)
	ymax = SlblFile.extent.YMax - (cellSize/2.0)
	
	# Imports the rasters in numpy for processing
	SlblArray = arcpy.RasterToNumPyArray (SlblFile,nodata_to_value=0)
	DepthArray = arcpy.RasterToNumPyArray (DepthFile,nodata_to_value=0)
	
	path_topo_file = os.path.join(outFolder,'Path_topography.grd')
	source_depth_file = os.path.join(outFolder,'Source_depth.grd')
	
	if Erosion == True:
		#The assumption is that material 2 is attributed in every cell where the value is higher than 0 and material 1 where it is 0
		ErosionArray = arcpy.RasterToNumPyArray (ErosionFile,nodata_to_value=0)
		erosion_map_file = os.path.join(outFolder,'Erosion_map.grd')
		ErosionArray[ErosionArray>0]=2
		ErosionArray[ErosionArray==0]=1
		write_surfgrd(ErosionArray,erosion_map_file,cellSize)
	
	write_surfgrd(SlblArray,path_topo_file,cellSize) #save outside database
	write_surfgrd(DepthArray,source_depth_file,cellSize)
	
	# Saves a header file with real coordinates
	with open(os.path.join(outFolder,'header.txt'),'w') as file:
		file.write("DSAA\n")
		file.write('{0} {1}\n'.format(SlblArray.shape[1], SlblArray.shape[0]))
		file.write('{0} {1}\n'.format(xmin, xmax)) #after resizing
		file.write('{0} {1}\n'.format(ymin, ymax)) #after resizing

