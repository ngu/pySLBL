# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# pySLBL.py
# Created on: 23.01.2017
# by: Pierrick.Nicolet@ngu.no
#  ---------------------------------------------------------------------------
#
# This tool is an implementation of the Sloping Local Base Level (SLBL) in ArcGIS.
# The SLBL method is described in the references below. It creates a 
# potential failure surface using a digital elevation model and a mask representing the
# unstable area (as a polygon shapefile). It may be use for other purposes as well,
# for example to reconstruct the initial topography under a deposit. The parameters that need to be defined 
# are the tolerance, that controls the curvature, the maximum thickness and a 
# stop criterion
#
# There is no publication direclty related to this implementation of the tool,
# but you may cite the Github repository:
# Nicolet, P. (2020): pySLBL, https://github.com/ngu/pySLBL
#
# References:
#
# Jaboyedoff, M.; Carrea, D.; Derron, M.-H.; Oppikofer, T.; Penna, I. M. & Rudaz, B.
# A review of methods used to estimate initial landslide failure surface depths and volumes 
# Engineering Geology, 2020, 267, 105478
#
# Jaboyedoff, M.; Couture, R. & Locat, P.
# Structural analysis of Turtle Mountain (Alberta) using digital elevation model: Toward a progressive failure 
# Geomorphology, 2009, 103, 5 - 16
#
# Jaboyedoff, M.; Baillifard, F.; Couture, R.; Locat, J.; Locat, P. & Rouiller, J.
# Toward preliminary hazard assessment using DEM topographic analysis and simple mechanic modeling 
# Landslides evaluation and stabilization. Balkema, Rotterdam, 2004, 191-197
#
# ---------------------------------------------------------------------------


import arcpy
import numpy as np
from copy import deepcopy
import os
import sys
from scipy.interpolate import griddata
import datetime
import inspect
import string, random


if __name__=="__main__":
	# Retrive the parameters from the GUI
	grid_dem_file = arcpy.GetParameterAsText(0)
	mask_file = arcpy.GetParameterAsText(1)
	grid_out_file = arcpy.GetParameterAsText(2)
	if grid_out_file.find('#verbose') != -1:
		verbose = True
		grid_out_file = grid_out_file.split('#')[0]
	else:
		verbose= False
	tol = float(arcpy.GetParameterAsText(3).replace(',','.'))
	maxt = arcpy.GetParameterAsText(4)
	if len(maxt)==0:
		maxt = np.inf
	else:
		maxt = float(maxt.replace(',','.'))
	maxv = arcpy.GetParameterAsText(5)
	if len(maxv)==0:
		maxv = np.inf
	else:
		maxv = float(maxv.replace(',','.'))
	stop = float(arcpy.GetParameterAsText(6).replace(',','.'))
	method = arcpy.GetParameterAsText(7)
	if method == '4 neighbours, average' or method == '4 neighbours, min/max':
		nb_neigh = 4
	elif method == '8 neighbours, average' or method == '8 neighbours, min/max':
		nb_neigh = 8
	else:
		arcpy.AddError('Unknown method')
	if method == '4 neighbours, min/max' or method == '8 neighbours, min/max':
		criteria = 'minmax'
	elif method == '4 neighbours, average' or method == '8 neighbours, average':
		criteria = 'average'
	not_deepen = arcpy.GetParameterAsText(8)
	inverse = arcpy.GetParameterAsText(9)
	grid_diff_out = arcpy.GetParameterAsText(10)
	grid_hill_out = arcpy.GetParameterAsText(11)
	
	# Check the necessary extensions
	arcpy.CheckOutExtension("3D")
	if arcpy.CheckOutExtension("3D")!="CheckedOut":
		arcpy.CheckOutExtension("Spatial")
		if arcpy.CheckOutExtension("Spatial")!="CheckedOut":
			arcpy.AddMessage("A 3D analyst or Spatial analyst license is needed to run this tool")

	# Check if the results are saved in a database
	desc = arcpy.Describe(os.path.dirname(grid_out_file))
	if desc.workspaceType == 'LocalDatabase':
		saveInGDB = True
	else:
		saveInGDB = False

	# Get the name and extension of the output
	if os.path.basename(grid_out_file).find('.') != -1:
		ext='.' + os.path.basename(grid_out_file).split('.')[-1]
		out_name = os.path.join(os.path.dirname(grid_out_file),'.'.join(os.path.basename(grid_out_file).split('.')[0:-1]))
	else:
		ext=''
		out_name = grid_out_file
		arcpy.env.workspace = os.path.dirname(grid_out_file)

	grid_dem_file = arcpy.Raster(grid_dem_file)
	

# Convert the polygon features to a raster mask
	try:
		mask_desc = arcpy.Describe(mask_file)
		try:
			# Retrieve the selected polygons (works if the input is a layer)
			Set = mask_desc.FIDSet
		except:
			#makes a layer first if the input is a file
			if arcpy.GetInstallInfo()['ProductName'] == 'Desktop':
				mask_file = arcpy.mapping.Layer(mask_file)
			else:
				mask_file = arcpy.mp.Layer(mask_file)
			mask_desc = arcpy.Describe(mask_file)
			# Retrieve the selected polygons
			Set = mask_desc.FIDSet

		if mask_desc.shapeType != "Polygon":
			arcpy.AddError('The mask layer must be a polygon layer')
		else:
			# Retrive key parameters from the DEM that will be used for the raster mask
			arcpy.env.extent = grid_dem_file.extent
			arcpy.env.outputCoordinateSystem = arcpy.Describe(grid_dem_file).spatialReference
			cellSize = grid_dem_file.meanCellWidth
			
			# Generate a random name for intermediate data (to avoid conflict if previous intermediate data weren't correctly deleted)
			lower_alphabet = string.ascii_lowercase
			random_string =''.join(random.choice(lower_alphabet) for i in range(3))
			grid_mask_file_temp = 'mask_temp_' + random_string
			if verbose:
				arcpy.AddMessage('grid_mask_file_temp initial = {}'.format(grid_mask_file_temp))
			grid_mask_file_temp = arcpy.ValidateTableName(grid_mask_file_temp,os.path.dirname(out_name))
			grid_mask_file_temp = os.path.join(os.path.dirname(out_name), grid_mask_file_temp + ext)
			if verbose:
				arcpy.AddMessage('grid_mask_file_temp validated = {}'.format(grid_mask_file_temp.encode('utf-8','replace')))
			
			# Generate a raster using the OID field as value
			oid_fieldname = mask_desc.OIDFieldName
			if verbose:
				arcpy.AddMessage('Processing extent: x:[{}:{}], y:[{}:{}]'.format(arcpy.env.extent.XMin,arcpy.env.extent.XMax,arcpy.env.extent.YMin,arcpy.env.extent.YMax))
				arcpy.AddMessage('Cellsize: {}'.format(cellSize))
				arcpy.AddMessage('Coordinate system: {}'.format(arcpy.env.outputCoordinateSystem.name))
				arcpy.AddMessage('OID field name = {}'.format(oid_fieldname))
			arcpy.PolygonToRaster_conversion (mask_file, oid_fieldname, grid_mask_file_temp, "CELL_CENTER", "#", cellSize)

			# Reclassify the raster to keep selected polygons (or all polygons if none is selected) as the mask
			listValue = []
			fields = ['OID@']
			rows = arcpy.da.SearchCursor(mask_file, fields)
			if len(Set) > 0:
				for row in rows:
					feat = row[0]
					if str(feat) in Set:
						listValue.append([feat,1])
					else:
						listValue.append([feat,0])
			else:
				#no selected polygons --> takes all
				for row in rows:
					feat = row[0]
					listValue.append([feat,1])
			listValue.append(['NoData',0]) #must be at the end of the list
			if verbose:
				msg = '[%s]' % ', '.join(map(str, listValue))
				arcpy.AddMessage(msg)
				arcpy.AddMessage('path to Reclassify: {}'.format(os.path.abspath(inspect.getfile(arcpy.sa.Reclassify))))
				arcpy.AddMessage('Workspace: {}'.format(arcpy.env.workspace))
				arcpy.AddMessage('Scratch workspace: {}'.format(arcpy.env.scratchWorkspace))
			myRemapVal = arcpy.sa.RemapValue(listValue)
			grid_mask_file = arcpy.sa.Reclassify(grid_mask_file_temp, "Value", myRemapVal, "NODATA")
			if verbose:
				arcpy.AddMessage('grid_mask_file reclassified')
				desc = arcpy.Describe(grid_mask_file)
				layersource = os.path.join(str(desc.path), str(desc.name))
				arcpy.AddMessage('reclassified mask (intermediate) saved in: {}'.format(layersource))
			arcpy.Delete_management(grid_mask_file_temp)
			if verbose:
				arcpy.AddMessage('grid_mask_file_temp deleted')
			grid_mask = arcpy.RasterToNumPyArray (grid_mask_file,nodata_to_value=0)
			if verbose:
				arcpy.AddMessage('grid_mask transformed to numpy array')
	except:
		arcpy.AddError("Unexpected error when converting polygon to mask")
		e = sys.exc_info()[1]
		arcpy.AddError(e.args[0])
		try:
			arcpy.Delete_management(grid_mask_file)
		except:
			pass
		sys.exit(0)

	# Reduces the processing extent to accelerate the computation (keeps 2 cells around the masked area)
	s = grid_mask.shape
	i_min = s[0]
	i_max = 0
	j_min = s[1]
	j_max = 0

	for i in range(0,s[0]-1):
		for j in range(0,s[1]-2):
			if grid_mask[i,j] > 0:
				if i < i_min:
					i_min = deepcopy(i)
				if i > i_max:
					i_max = deepcopy(i)
				if j < j_min:
					j_min = deepcopy(j)
				if j > j_max:
					j_max = deepcopy(j)

	i_min = i_min - 2
	i_max = i_max + 2
	j_min = j_min - 2
	j_max = j_max + 2

	i_extent=i_max-i_min
	j_extent=j_max-j_min

	x_min = grid_mask_file.extent.XMin + (j_min*cellSize)
	y_min = grid_mask_file.extent.YMin + ((s[0]-i_max)*cellSize)
	lowerLeftArea= arcpy.Point(x_min,y_min)
	if verbose:
		arcpy.AddMessage('x_min = {}, y_min = {}, i_extent = {}, j_extent = {}'.format(x_min,y_min,i_extent,j_extent))
		
	# Convert the grid in numpy for computation
	grid_dem = arcpy.RasterToNumPyArray (grid_dem_file,lowerLeftArea,j_extent,i_extent,nodata_to_value=-9999)
	grid_dem = grid_dem.astype(np.float32)
	s2=grid_dem.shape

	# interpolate to fill missing data in the DEM (inspired from:
	# https://stackoverflow.com/questions/12923593/interpolation-of-sparse-grid-using-python-preferably-scipy)
	grid_dem[grid_dem==-9999]=np.nan
	mask = np.isfinite(grid_dem)
	if not mask.all(): #there are cell with no data
		str_message = 'Filling missing data in the DEM using the nearest neighbours'
		arcpy.AddMessage(str_message)
		index = np.where(mask==True)
		points = np.dstack((index[0],index[1]))
		values = grid_dem[mask].flatten()
		values = np.array(values)
		points = points.reshape(len(values),2)
		values = values.reshape(len(values),1)
		grid_x, grid_y = np.mgrid[0:s2[0]:1,0:s2[1]:1]
		grid_dem = griddata(points, values, (grid_x, grid_y), method='nearest')
		grid_dem = np.squeeze(grid_dem, axis=2)

	grid_mask = arcpy.RasterToNumPyArray (grid_mask_file,lowerLeftArea,j_extent,i_extent,nodata_to_value=0)
	grid_mask = (grid_mask[:] > 0)

	# initilizes thickness and difference grids
	grid_thickn = np.zeros(s2)
	grid_diff = np.ones(s2)

	# Creates a matrice to store the values of the neighbouring cells in the previous iteration
	mat_neigh = np.zeros(s2)
	mat_neigh = np.expand_dims(mat_neigh,axis=2)
	if nb_neigh ==4:
		mat_neigh = np.tile(mat_neigh,(1,1,4))
	else: 
		mat_neigh = np.tile(mat_neigh,(1,1,8))

	# Creates a matrice where the proposed value and previous value are stored for comparison
	mat_comp = np.zeros(s2)
	mat_comp = np.expand_dims(mat_comp,axis=2)
	mat_comp = np.tile(mat_comp,(1,1,2))

	# Initiate the slbl grid (altitude)
	grid_slbl = deepcopy(grid_dem)

	# Retrieves the minimum altitude around the landslide if this condition (not deepening) is set
	if not_deepen == 'true':
		grid_mask_border = arcpy.sa.Expand(grid_mask_file, 1, [1])
		grid_mask_border_np = arcpy.RasterToNumPyArray (grid_mask_border,lowerLeftArea,j_extent,i_extent,nodata_to_value=0)
		grid_mask_border_np = (grid_mask_border_np[:] > 0)
		z_min = np.nanmin(grid_dem[grid_mask_border_np])
		str_message = 'Minimum altitude in or around the landslide area {} m'.format(str(z_min))
		arcpy.AddMessage(str_message)

	iter = 0
	volume = 0.

	# The SLBL strarts here
	while np.amax(grid_diff)>stop and np.amax(grid_thickn)<maxt and volume<maxv:
		iter=iter+1
		grid_thickn_prev = deepcopy(grid_thickn)
		grid_slbl_prev = deepcopy(grid_slbl)
		
		# writes the values of the neighbourings cells in the 3rd dimension of the matrix
		mat_neigh[:-1,:,0]=grid_slbl_prev[1:,:]
		mat_neigh[1:,:,1]=grid_slbl_prev[:-1,:]
		mat_neigh[:,:-1,2]=grid_slbl_prev[:,1:]
		mat_neigh[:,1:,3]=grid_slbl_prev[:,:-1]
		
		# diagonals
		if nb_neigh ==8:
			mat_neigh[:-1,:-1,4]=grid_slbl_prev[1:,1:]
			mat_neigh[:-1,1:,5]=grid_slbl_prev[1:,:-1]
			mat_neigh[1:,1:,6]=grid_slbl_prev[:-1,:-1]
			mat_neigh[1:,:-1,7]=grid_slbl_prev[:-1,1:]
		
		if criteria == 'minmax':
			mat_max=np.amax(mat_neigh,axis=2)
			mat_min=np.amin(mat_neigh,axis=2)
			mat_mean=(mat_max+mat_min)/2
		elif criteria == 'average':
			mat_mean=np.mean(mat_neigh,axis=2)
		mat_mean=mat_mean+tol
		if not_deepen == 'true':
			mat_mean=np.maximum(mat_mean,z_min)
		
		mat_comp[:,:,0]=mat_mean
		mat_comp[:,:,1]=grid_slbl
		
		# Check if the new value should be kept
		if inverse == 'true':
			grid_slbl=np.amax(mat_comp,axis=2)
		else:
			grid_slbl=np.amin(mat_comp,axis=2)
		
		# Replaces the values of the SLBL by the original values outside the masked area
		grid_slbl[~grid_mask]=grid_dem[~grid_mask]
		
		grid_thickn = np.absolute(grid_dem - grid_slbl)
		grid_diff = np.absolute(grid_thickn - grid_thickn_prev)
		
		volume = (np.sum(grid_thickn)*cellSize*cellSize)
		
		if iter%100==0:
			str_message = '{0} iterations. Max diff is {1}, max thickness is {2} and volume is {3}'.format(str(iter),str(np.amax(grid_diff)),str(np.amax(grid_thickn)),str(volume))
			arcpy.AddMessage(str_message)
	# The SLBL is finished

	# Output the main results to the console
	str_message = 'SLBL computed in {} iterations'.format(str(iter))
	arcpy.AddMessage(str_message)
	str_message = 'Maximum thickness: {} m'.format(str(np.amax(grid_thickn)))
	arcpy.AddMessage(str_message)
	str_message = 'Average thickness: {} m'.format(str(np.mean(grid_thickn[grid_mask])))
	arcpy.AddMessage(str_message)
	volume = (np.sum(grid_thickn)*cellSize*cellSize)
	str_message = 'Total volume: {} million m3'.format(str(volume/1000000))
	arcpy.AddMessage(str_message)
	if volume < 251464.767769637:
		scheidegger = 30.9637565320735
	else:
		scheidegger = np.rad2deg(np.arctan(np.power(volume,-0.15666)*np.power(10,0.62419)))
	str_message = 'Angle of reach: {} degrees'.format(str(scheidegger))
	arcpy.AddMessage(str_message)

	# prints the results if the command line is used
	print('SLBL done in {} iterations'.format(str(iter)))
	print('Volume is: {} million m3'.format(str((np.sum(grid_thickn)*cellSize*cellSize)/1000000)))
	print('Max height is: {} m'.format(str(np.amax(grid_thickn))))

	str_message = 'Saving elevation file...'
	arcpy.AddMessage(str_message)

	# Saves the results in the original grid extent
	lowerLeftExtended= arcpy.Point(grid_mask_file.extent.XMin,grid_mask_file.extent.YMin)
	grid_dem_after = arcpy.RasterToNumPyArray (grid_dem_file,lowerLeftExtended,s[1],s[0],nodata_to_value=-9999)
	grid_dem_after = grid_dem_after.astype(np.float32)

	grid_dem_after[i_min:i_max,j_min:j_max]=grid_slbl

	grid_dem_out = arcpy.NumPyArrayToRaster(grid_dem_after,lowerLeftExtended,grid_dem_file.meanCellWidth, value_to_nodata=-9999)
	grid_dem_out.save(grid_out_file)
	
	# Saves the thickness file (if selected) and tries to add it to the project
	if grid_diff_out== 'true':
		str_message = 'Saving thickness file...'
		arcpy.AddMessage(str_message)
		grid_thick_after=np.zeros(grid_dem_after.shape)
		grid_thick_after[i_min:i_max,j_min:j_max]=grid_thickn
		grid_out = arcpy.NumPyArrayToRaster(grid_thick_after,lowerLeftExtended,grid_dem_file.meanCellWidth)
		
		grid_out_thick_file = out_name + '_t' + ext
		if not (saveInGDB or len(ext) > 0):
			grid_out_thick_file = os.path.join(os.path.dirname(out_name), os.path.basename(out_name)[0:11] + '_t')
		grid_out.save(grid_out_thick_file)
		grid_out_thick_file_name=os.path.basename(grid_out_thick_file)
		
		if arcpy.GetInstallInfo()['ProductName'] == 'Desktop':
			try:
				mxd = arcpy.mapping.MapDocument("CURRENT")
				dataFrame = arcpy.mapping.ListDataFrames(mxd, "*")[0]
				result = arcpy.MakeRasterLayer_management(grid_out_thick_file,grid_out_thick_file_name)
				layer = result.getOutput(0)
				arcpy.mapping.AddLayer(dataFrame, layer)
			except:
				str_message = 'Thickness raster has been saved but could not be added to the current project'
				arcpy.AddMessage(str_message)
		else:
			mxd = arcpy.mp.ArcGISProject("CURRENT")
			dataFrame = mxd.listMaps("*")[0]
			dataFrame.addDataFromPath(grid_out_thick_file)
		try:
			del grid_thick_after, grid_out, grid_out_thick_file, dataFrame
		except:
			pass
	
	# Saves the hillshade file (if selected) and tries to add it to the project
	if grid_hill_out == 'true':
		str_message = 'Saving hillshade file...'
		arcpy.AddMessage(str_message)
		arcpy.CheckOutExtension("3D")
		if arcpy.CheckOutExtension("3D")!="CheckedOut":
			arcpy.AddMessage("A 3D analyst license is needed for the hillshade")
		grid_out_hill_file = out_name + '_hshd' + ext
		if not (saveInGDB or len(ext) > 0):
			grid_out_hill_file = os.path.join(os.path.dirname(out_name), os.path.basename(out_name)[0:11] + '_h')
		arcpy.HillShade_3d(grid_dem_out, grid_out_hill_file,315, 45, 'NO_SHADOWS', 1)
		
		grid_out_hill_file_name=os.path.basename(grid_out_hill_file)

		if arcpy.GetInstallInfo()['ProductName'] == 'Desktop':
			try:
				mxd = arcpy.mapping.MapDocument("CURRENT")
				dataFrame = arcpy.mapping.ListDataFrames(mxd, "*")[0]
				result = arcpy.MakeRasterLayer_management(grid_out_hill_file,grid_out_hill_file_name)
				layer = result.getOutput(0)
				arcpy.mapping.AddLayer(dataFrame, layer)
			except:
				str_message = 'Hillshade raster has been saved but could not be added to the current project'
				arcpy.AddMessage(str_message)
		else:
			mxd = arcpy.mp.ArcGISProject("CURRENT")
			dataFrame = mxd.listMaps("*")[0]
			dataFrame.addDataFromPath(grid_out_hill_file)
		arcpy.CheckInExtension("3D")
		try:
			del grid_out_hill_file, dataFrame
		except:
			pass
	
	# Saves the key parameters in a table (creates the table if necessary)
	# The length of field names is limited to 10 characters in a dbf file
	if saveInGDB:
		summaryTableName = 'SLBL_results'
	else:
		summaryTableName = 'SLBL_results.dbf'
		
	summaryTable = os.path.join(os.path.dirname(grid_out_file),summaryTableName)
		
	if not arcpy.Exists(summaryTable):
		str_message = 'Creating summary table...'
		arcpy.AddMessage(str_message)
		arcpy.CreateTable_management(os.path.dirname(grid_out_file),summaryTableName)
		arcpy.AddField_management(summaryTable,"Name","TEXT")
		arcpy.AddField_management(summaryTable,"Volume_10e6m3","FLOAT")
		if saveInGDB:
			arcpy.AddField_management(summaryTable,"Reach_angle","FLOAT")
		else: #max length = 10
			arcpy.AddField_management(summaryTable,"Reach_angl","FLOAT")
		arcpy.AddField_management(summaryTable,"Max_thick","FLOAT")
		arcpy.AddField_management(summaryTable,"Avg_thick","FLOAT")
		arcpy.AddField_management(summaryTable,"Iterations","LONG")
		arcpy.AddField_management(summaryTable,"Cell_size","SHORT")
		arcpy.AddField_management(summaryTable,"Tolerance","FLOAT")
		arcpy.AddField_management(summaryTable,"Max_depth","FLOAT")
		arcpy.AddField_management(summaryTable,"Max_vol","FLOAT")
		arcpy.AddField_management(summaryTable,"Method","TEXT")
		if saveInGDB:
			arcpy.AddField_management(summaryTable,"Stop_criterion","FLOAT")
			arcpy.AddField_management(summaryTable,"Nb_neighbours","SHORT")
			arcpy.AddField_management(summaryTable,"Not_deepening","TEXT")
		else:
			arcpy.AddField_management(summaryTable,"Stop_crite","FLOAT")
			arcpy.AddField_management(summaryTable,"Nb_neighbo","SHORT")
			arcpy.AddField_management(summaryTable,"Not_deepen","TEXT")
		arcpy.AddField_management(summaryTable,"Inverse","TEXT")
		arcpy.AddField_management(summaryTable,"Date","DATE")
		
	str_message = 'Filling summary table...'
	arcpy.AddMessage(str_message)
	field_names = [f.name for f in arcpy.ListFields(summaryTable)]
	if not "Method" in field_names: #field added in a later version
		arcpy.AddField_management(summaryTable,"Method","TEXT")
	if not "Inverse" in field_names: #field added in a later version
		arcpy.AddField_management(summaryTable,"Inverse","TEXT")
	if not "Max_vol" in field_names: #field added in a later version
		arcpy.AddField_management(summaryTable,"Max_vol","FLOAT")
	if "Volume_Mm3" in field_names: #field name changed in a later version
		arcpy.AlterField_management(summaryTable, "Volume_Mm3", 'Volume_10e6m3')
	cur = arcpy.InsertCursor(summaryTable)
	row = cur.newRow()
	row.setValue('Name', os.path.basename(out_name))
	row.setValue('Volume_10e6m3', (volume/1000000))
	if saveInGDB:
		row.setValue('Reach_angle', scheidegger)
	else:
		row.setValue('Reach_angl', scheidegger)
	row.setValue('Max_thick', np.amax(grid_thickn))
	row.setValue('Avg_thick', np.mean(grid_thickn[grid_mask]))
	row.setValue('Iterations', iter)
	row.setValue('Cell_size', cellSize)
	row.setValue('Tolerance', tol)
	if np.isinf(maxt):
		row.setValue('Max_depth', -1)
	else:
		row.setValue('Max_depth', maxt)
	if np.isinf(maxv):
		row.setValue('Max_vol', -1)
	else:
		row.setValue('Max_vol', maxv)
	row.setValue('Method', criteria)
	if saveInGDB:
		row.setValue('Stop_criterion', stop)
		row.setValue('Nb_neighbours', int(nb_neigh))
		row.setValue('Not_deepening', not_deepen)
	else:
		row.setValue('Stop_crite', stop)
		row.setValue('Nb_neighbo', int(nb_neigh))
		row.setValue('Not_deepen', not_deepen)
	row.setValue('Inverse', inverse)
	row.setValue('Date', datetime.datetime.today())
	cur.insertRow(row)
	
	try:
		del grid_dem_out, grid_dem_after, cur, row
		arcpy.Delete_management(grid_mask_file)
	except:
		pass