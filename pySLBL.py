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
#import sys
from scipy.interpolate import griddata
import datetime
import inspect
import string, random


def SLBL(grid_dem,grid_mask,tol,maxt,maxv,z_min,planes=None):
	# initilizes thickness and difference grids
	s=grid_dem.shape
	grid_thickn = np.zeros(s)
	grid_diff = np.ones(s)

	# Creates a matrice to store the values of the neighbouring cells in the previous iteration
	mat_neigh = np.zeros(s)
	mat_neigh = np.expand_dims(mat_neigh,axis=2)
	if nb_neigh ==4:
		mat_neigh = np.tile(mat_neigh,(1,1,4))
	else: 
		mat_neigh = np.tile(mat_neigh,(1,1,8))

	# Creates a matrice where the proposed value and previous value are stored for comparison
	mat_comp = np.zeros(s)
	mat_comp = np.expand_dims(mat_comp,axis=2)
	mat_comp = np.tile(mat_comp,(1,1,2))

	# Initiate the slbl grid (altitude)
	grid_slbl = deepcopy(grid_dem)

	nb_iter = 0
	volume = 0.

	if np.isfinite(maxt):
		grid_maxt = grid_dem - maxt

	# The SLBL strarts here
	while np.amax(grid_diff)>stop and volume<maxv:
		nb_iter=nb_iter+1
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
		# limit to the lower altitude around the polygon
		if np.isfinite(z_min):
			mat_mean=np.maximum(mat_mean,z_min)
		# limit to the maximum thickness
		if np.isfinite(maxt):
			mat_mean=np.maximum(mat_mean,grid_maxt)
		if not planes is None:
			mat_mean=np.maximum(mat_mean,planes)
		
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
		
		if nb_iter%100==0:
			str_message = '{} iterations. Max diff is {:.3e} m, max thickness is {:.3f} m and volume is {:.6f} million m3'.format(nb_iter,np.amax(grid_diff),np.amax(grid_thickn),volume/1000000)
			arcpy.AddMessage(str_message)
	# The SLBL is finished
	
	return grid_slbl, grid_thickn, nb_iter

def define_extent(poly_extent,dem_extent,cellSize,src_mask = False):
	# Define the processing area from the polygon(s) extent and the DEM
	
	if src_mask != False:
		pt_ll = arcpy.PointGeometry(arcpy.Point(poly_extent.XMin,poly_extent.YMin),src_mask).projectAs(dem_extent.spatialReference).centroid
		pt_ur = arcpy.PointGeometry(arcpy.Point(poly_extent.XMax,poly_extent.YMax),src_mask).projectAs(dem_extent.spatialReference).centroid
		pt_ul = arcpy.PointGeometry(arcpy.Point(poly_extent.XMin,poly_extent.YMax),src_mask).projectAs(dem_extent.spatialReference).centroid
		pt_lr = arcpy.PointGeometry(arcpy.Point(poly_extent.XMax,poly_extent.YMin),src_mask).projectAs(dem_extent.spatialReference).centroid

		pt_min = arcpy.Point(min(pt_ll.X,pt_ul.X), min(pt_ll.Y,pt_lr.Y))
		pt_max = arcpy.Point(max(pt_lr.X,pt_ur.X), max(pt_ul.Y,pt_ur.Y))

		if verbose:
			str_message = 'Polygon src:{}'.format(src_mask.name)
			arcpy.AddMessage(str_message)
			str_message = 'Polygon: x:{}-{}, y:{}-{}, src:{}'.format(pt_min.X,pt_max.X,pt_min.Y,pt_max.Y,dem_extent.spatialReference.name)
			arcpy.AddMessage(str_message)
	else:
		pt_min = arcpy.Point(poly_extent.XMin, poly_extent.YMin)
		pt_max = arcpy.Point(poly_extent.XMax, poly_extent.YMax)

		if verbose:
			try:
				str_message = 'Polygon: x:{}-{}, y:{}-{}, src:{}'.format(poly_extent.XMin,poly_extent.XMax,poly_extent.YMin,poly_extent.YMax,poly_extent.spatialReference.name)
			except: #If the extent is created in the script, no SR can be defined
				str_message = 'Polygon: x:{}-{}, y:{}-{}, src:{}'.format(poly_extent.XMin,poly_extent.XMax,poly_extent.YMin,poly_extent.YMax,mask_desc.extent.spatialReference.name)
			arcpy.AddMessage(str_message)
	
	if verbose:
		str_message = 'DEM: x:{}-{}, y:{}-{}, src:{}'.format(dem_extent.XMin,dem_extent.XMax,dem_extent.YMin,dem_extent.YMax,dem_extent.spatialReference.name)
		arcpy.AddMessage(str_message)
	
	if (pt_min.X < dem_extent.XMin) or (pt_max.X > dem_extent.XMax) or (pt_min.Y < dem_extent.YMin) or (pt_max.Y > dem_extent.YMax):
		arcpy.AddError('The polygon layer should be completely contained by the DEM')
	else:
		# Looks for the minimum extent of the DEM that includes the polyon(s) and adds 2 pixels
		# if the coordinate of the polygon falls on a grid line, only one pixel remain around
		# the polygon on that side. One pixel is needed when the "not deepening" option
		# is activated
		xmin = deepcopy(dem_extent.XMin)
		while xmin < pt_min.X:
			xmin += cellSize
		xmin = xmin - (2*cellSize)
			
		xmax = deepcopy(dem_extent.XMax)
		while xmax > pt_max.X:
			xmax = xmax - cellSize
		xmax = xmax + (2*cellSize)
			
		ymin = deepcopy(dem_extent.YMin)
		while ymin < pt_min.Y:
			ymin += cellSize
		ymin = ymin - (2*cellSize)
			
		ymax = deepcopy(dem_extent.YMax)
		while ymax > pt_max.Y:
			ymax = ymax - cellSize
		ymax = ymax + (2*cellSize)
	
		if (xmin < dem_extent.XMin) or (xmax > dem_extent.XMax) or (ymin < dem_extent.YMin) or (ymax > dem_extent.YMax):
			arcpy.AddError('The DEM should be at least 2 pixels larger than the polygon(s)')
		
		if verbose:
			str_message = 'Processing extent: x:{}-{}, y:{}-{}, cellsize:{}, src:{}'.format(xmin,xmax,ymin,ymax,cellSize,dem_extent.spatialReference.name)
			arcpy.AddMessage(str_message)

	return arcpy.Extent(xmin, ymin, xmax, ymax)

def raster2numpy(ws,ext,mask_file,mask_desc,not_deepen,listValue):
	# This function creates the numpy grids for processing
	
	# Generate a random name for intermediate data (to avoid conflict if previous intermediate data weren't correctly deleted)
	lower_alphabet = string.ascii_lowercase
	random_string =''.join(random.choice(lower_alphabet) for i in range(5))
	grid_mask_file_temp = 'mask_temp_' + random_string
	if verbose:
		arcpy.AddMessage(u'grid_mask_file_temp initial = {}'.format(grid_mask_file_temp))
	grid_mask_file_temp = arcpy.ValidateTableName(grid_mask_file_temp,ws)
	grid_mask_file_temp = os.path.join(ws, grid_mask_file_temp + ext)
	if verbose:
		arcpy.AddMessage(u'grid_mask_file_temp validated = {}'.format(grid_mask_file_temp))
	
	# Generate a raster using the OID field as value
	oid_fieldname = mask_desc.OIDFieldName
	if verbose:
		arcpy.AddMessage('Processing extent: x:[{}:{}], y:[{}:{}]'.format(arcpy.env.extent.XMin,arcpy.env.extent.XMax,arcpy.env.extent.YMin,arcpy.env.extent.YMax))
		arcpy.AddMessage('Cellsize: {}'.format(arcpy.env.cellSize))
		arcpy.AddMessage('Coordinate system: {}'.format(arcpy.env.outputCoordinateSystem.name))
		arcpy.AddMessage(u'OID field name = {}'.format(oid_fieldname))
	arcpy.PolygonToRaster_conversion (mask_file, oid_fieldname, grid_mask_file_temp, "CELL_CENTER", "#", cellSize)

	if verbose:
		msg = '[%s]' % ', '.join(map(str, listValue))
		arcpy.AddMessage(msg)
		arcpy.AddMessage(u'path to Reclassify: {}'.format(os.path.abspath(inspect.getfile(arcpy.sa.Reclassify))))
		arcpy.AddMessage(u'Workspace: {}'.format(arcpy.env.workspace))
		arcpy.AddMessage(u'Scratch workspace: {}'.format(arcpy.env.scratchWorkspace))

	myRemapVal = arcpy.sa.RemapValue(listValue)
	grid_mask_file = arcpy.sa.Reclassify(grid_mask_file_temp, "Value", myRemapVal, "NODATA")
	if verbose:
		arcpy.AddMessage('grid_mask_file reclassified')
		desc = arcpy.Describe(grid_mask_file)
		layersource = os.path.join(desc.path, desc.name)
		arcpy.AddMessage(u'reclassified mask (intermediate) saved in: {}'.format(layersource))
	arcpy.Delete_management(grid_mask_file_temp)
	if verbose:
		arcpy.AddMessage('grid_mask_file_temp deleted')
	grid_mask = arcpy.RasterToNumPyArray (grid_mask_file,nodata_to_value=0)
	grid_mask = (grid_mask[:] > 0)
	if verbose:
		arcpy.AddMessage('grid_mask transformed to numpy array')
	
	lowerLeftArea= arcpy.Point(arcpy.env.extent.XMin,arcpy.env.extent.YMin)
	ncols = grid_mask_file.width
	nrows = grid_mask_file.height
	
	if verbose:
		arcpy.AddMessage('x_min = {}, y_min = {}, ncols = {}, nrows = {}'.format(arcpy.env.extent.XMin,arcpy.env.extent.YMin,ncols,nrows))
	
	# Convert the grid in numpy for computation
	grid_dem = arcpy.RasterToNumPyArray (grid_dem_file,lowerLeftArea,ncols,nrows,nodata_to_value=-9999)
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
		points = np.dstack((index[1],index[0]))
		values = grid_dem[mask].flatten()
		values = np.array(values)
		points = points.reshape(len(values),2)
		values = values.reshape(len(values),1)
		grid_y, grid_x = np.mgrid[0:s2[0]:1,0:s2[1]:1]
		grid_dem = griddata(points, values, (grid_x, grid_y), method='nearest')
		grid_dem = np.squeeze(grid_dem, axis=2)
		
	# Retrieves the minimum altitude around the landslide if this condition (not deepening) is set
	if not_deepen == 'true':
		grid_mask_border = arcpy.sa.Expand(grid_mask_file, 1, [1])
		grid_mask_border_np = arcpy.RasterToNumPyArray (grid_mask_border,lowerLeftArea,ncols,nrows,nodata_to_value=0)
		grid_mask_border_np = (grid_mask_border_np[:] > 0)
		z_min = np.nanmin(grid_dem[grid_mask_border_np])
		str_message = 'Minimum altitude in or around the landslide area {} m'.format(str(z_min))
		arcpy.AddMessage(str_message)
	else:
		z_min = np.nan

	arcpy.Delete_management(grid_mask_file)
	return grid_mask, grid_dem, z_min

def definetolerance(grid_mask, grid_dem):
	# This function defines the tolerance when the automatic mode is used
	
	# Create grids of x and y in a local reference system
	s2=grid_dem.shape
	grid_y, grid_x = np.mgrid[0:s2[0]*cellSize:1*cellSize,0:s2[1]*cellSize:1*cellSize]
	grid_y = np.flipud(grid_y)
	
	if verbose:
		str_message = 'Altitude of cell [0,0] = {}'.format(grid_dem[0,0])
		arcpy.AddMessage(str_message)
	
	# Finds the central point of the polygon
	x_mean = np.mean(grid_x[grid_mask])
	y_mean = np.mean(grid_y[grid_mask])
	
	if verbose:
		str_message = 'x_mean = {}, y_mean = {}'.format(x_mean,y_mean)
		arcpy.AddMessage(str_message)
		str_message = 'x: [{}-{}], y: [{}-{}]'.format(np.min(grid_x[grid_mask]),np.max(grid_x[grid_mask]),np.min(grid_y[grid_mask]),np.max(grid_y[grid_mask]))
		arcpy.AddMessage(str_message)
	
	# Vectorizes the grids for the points inside the polygon
	x = grid_x[grid_mask].flatten()
	y = grid_y[grid_mask].flatten()
	z = grid_dem[grid_mask].flatten()
	
	x = np.expand_dims(x,axis=1)
	y = np.expand_dims(y,axis=1)
	z = np.expand_dims(z,axis=1)
	
	# Fits a plan on the surface of the landslide
	# Inspired from: https://gist.github.com/RustingSword/e22a11e1d391f2ab1f2c
	XY1 = np.concatenate((x,y,np.ones(x.shape)),axis=1)
	if verbose:
		str_message = 'XY1({},{}) is {} with range {}-{}, z({},{}) is {} with range {}-{}'.format(XY1.shape[0],XY1.shape[1],type(XY1),np.min(XY1),np.max(XY1),z.shape[0],z.shape[1],type(z),np.min(z),np.max(z))
		arcpy.AddMessage(str_message)
	(a, b, c),resid,rank,s = np.linalg.lstsq(XY1, z, rcond=1e-10)
	
	# Calculates the normal vector of the plan
	normal = np.array([a[0], b[0], -1])
	nn = np.linalg.norm(normal)
	normal = normal / nn
	
	# Ensures that the vector points upwards
	if normal[2] < 0:
		normal = normal * -1
		
	if verbose:
		str_message = 'normal vector ({},{},{}), c = {}'.format(normal[0],normal[1],normal[2],c[0])
		arcpy.AddMessage(str_message)
		
		# https://se.mathworks.com/matlabcentral/answers/342134-how-can-i-find-strike-and-dip-of-a-plane-from-equation-of-a-plane
		n_e = np.array([0,0,1])
		dip = np.degrees(np.arccos(np.dot(normal,n_e)))
		dip_dir = np.degrees(np.arctan2(normal[0], normal[1]))
		
		str_message = 'dip dir = {}, dip {}'.format(dip_dir,dip)
		arcpy.AddMessage(str_message)

	if savefigs:
		# plot points and fitted surface
		from mpl_toolkits.mplot3d import Axes3D
		import matplotlib.pyplot as plt
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		z_plane = a[0]*x + b[0]*y + c[0]
		ax.scatter(x, y, z_plane, c='b', s=5)
		ax.scatter(x, y, z, c='r', s=5)
		plt.xlabel('X')
		plt.ylabel('Y')
		ax.set_zlabel('Z')
		ax.axis('equal')
		ax.axis('tight')
		save_name = validatename(figspath,'plane','.png')
		fig.savefig(os.path.join(figspath,save_name))
		plt.close()
	
	# finds the slope and y-intercept of the profile in the xy plan
	slope = normal[1]/normal[0]
	y0 = y_mean - (slope * x_mean)
		
	xmin = 0
	ymin = 0
	xmax = arcpy.env.extent.XMax - arcpy.env.extent.XMin
	ymax = arcpy.env.extent.YMax - arcpy.env.extent.YMin
	
	if verbose:
		str_message = 'xmax: {}, ymax: {}, y0: {}'.format(xmax,ymax,y0)
		arcpy.AddMessage(str_message)
	
	# finds the coordinates (xa,ya) and (xb,yb) formed by the intersection of
	# the profile with the study area
	y_xmin = xmin * slope + y0
	y_xmax = xmax * slope + y0
	x_ymin = (ymin - y0)/slope
	x_ymax = (ymax - y0)/slope
	
	if slope < 0:
		if x_ymax > xmin:
			xa = x_ymax
			ya = ymax
		else:
			xa = xmin
			ya = y_xmin
		if x_ymin < xmax:
			xb = x_ymin
			yb = ymin
		else:
			xb = xmax
			yb = y_xmax
	else:
		if x_ymin > xmin:
			xa = x_ymin
			ya = ymin
		else:
			xa = xmin
			ya = y_xmin
		if x_ymax < xmax:
			xb = x_ymax
			yb = ymax
		else:
			xb = xmax
			yb = y_xmax
	
	if verbose:
		str_message = 'profile: ({},{})-({},{})'.format(xa,ya,xb,yb)
		arcpy.AddMessage(str_message)
	
	#-- Interpolate the altitudes on the profile
	# Vectorizes the grids for all the points
	x = grid_x.flatten()
	y = grid_y.flatten()
	z = grid_dem.flatten()
	
	x = np.expand_dims(x,axis=1)
	y = np.expand_dims(y,axis=1)
	z = np.expand_dims(z,axis=1)
	
	# Generate points evenly spaced along the profile
	x_profile = np.linspace(xa,xb,num=100)
	y_profile = np.linspace(ya,yb,num=100)
	xy_profile = np.column_stack((x_profile,y_profile))
	
	l_profile = np.sqrt((x_profile-xa)**2 + (y_profile-ya)**2)
	
	XY = np.concatenate((x,y),axis=1)
	
	profil_int = griddata(XY, z, xy_profile, method='cubic')
	mask_int = griddata(XY, grid_mask.flatten(), xy_profile, method='nearest')
	
	if savefigs:
		# plot profile
		import matplotlib.pyplot as plt
		fig = plt.figure()
		ax = fig.gca()
		ax.plot(l_profile,profil_int,'b')
		ax.plot(l_profile[mask_int],profil_int[mask_int],'r')
		ax.axis('equal')
		ax.axis('tight')
		save_name = validatename(figspath,'profile','.png')
		fig.savefig(os.path.join(figspath,save_name))
		plt.close()

	# Calculate the tolerance using the altitude range of the profile (dZ) and its length (dL)
	dZ = np.max(profil_int[mask_int]) - np.min(profil_int[mask_int])
	dL = np.max(l_profile[mask_int]) - np.min(l_profile[mask_int])
	tol_max = 4 * (1-np.sqrt(2)) * dZ * (cellSize**2/dL**2)

	if verbose:
		arcpy.AddMessage('dZ={} ({}-{}), dL={}'.format(dZ,profil_int[0],profil_int[-1],dL))
		arcpy.AddMessage('tol max: {}'.format(tol_max))
	
	tol_inter = 2 * (1-np.sqrt(2)) * dZ * (cellSize**2/dL**2)
	tol_min = 0	
	
	return tol_inter, tol_min, tol_max

def fillsummarytable(summaryTable,out_basename,out_basename_validated,grid_thickn,grid_mask,nb_iter,tol):
	# Output the main results to the console
	str_message = 'SLBL computed in {} iterations'.format(nb_iter)
	arcpy.AddMessage(str_message)
	str_message = 'Maximum thickness: {:.3f} m'.format(np.amax(grid_thickn))
	arcpy.AddMessage(str_message)
	str_message = 'Average thickness: {:.3f} m'.format(np.mean(grid_thickn[grid_mask]))
	arcpy.AddMessage(str_message)
	volume = (np.sum(grid_thickn)*cellSize*cellSize)
	str_message = 'Total volume: {:.6f} million m3'.format(volume/1000000)
	arcpy.AddMessage(str_message)
	if volume < 251464.767769637:
		scheidegger = 30.9637565320735
	else:
		scheidegger = np.rad2deg(np.arctan(np.power(volume,-0.15666)*np.power(10,0.62419)))
	str_message = 'Angle of reach: {:.3f} degrees'.format(scheidegger)
	arcpy.AddMessage(str_message)
	
	# Saves the key parameters in the summary table (creates the table if necessary)
	# The length of field names is limited to 10 characters in a dbf file
	if not arcpy.Exists(summaryTable):
		str_message = 'Creating summary table...'
		arcpy.AddMessage(str_message)
		arcpy.CreateTable_management(ws,summaryTableName)
		arcpy.AddField_management(summaryTable,"Name","TEXT")
		if saveInGDB:
			arcpy.AddField_management(summaryTable,"Volume_10e6m3","FLOAT")
			arcpy.AddField_management(summaryTable,"Reach_angle","FLOAT")
		else: #max length = 10
			arcpy.AddField_management(summaryTable,"Volume","FLOAT")
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
	
	if arcpy.TestSchemaLock(summaryTable) == False:
		str_message = 'The table {} is locked and cannot be edited here. Make sure it is not being edited elsewhere'.format(summaryTable)
		arcpy.AddError(str_message)
		
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
	if str(out_basename_validated) == str(out_basename):
		row.setValue('Name', out_basename)
	else:
		row.setValue('Name', '{} ({})'.format(out_basename,out_basename_validated))
	if saveInGDB:
		row.setValue('Volume_10e6m3', (volume/1000000))
		row.setValue('Reach_angle', scheidegger)
	else:
		row.setValue('Volume', (volume/1000000))
		row.setValue('Reach_angl', scheidegger)
	row.setValue('Max_thick', np.amax(grid_thickn))
	row.setValue('Avg_thick', np.mean(grid_thickn[grid_mask]))
	row.setValue('Iterations', nb_iter)
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
	
def savegrids(grid_slbl,grid_thickn,ws,out_basename,ext):

	if extent_mode == 'Full extent of the DEM':
		lowerLeft= arcpy.Point(dem_extent.XMin,dem_extent.YMin)
		ncols = dem_desc.width
		nrows = dem_desc.height
		
		i_min = round((dem_extent.YMax - processing_extent.YMax)/cellSize)
		i_max = round((dem_extent.YMax - processing_extent.YMin)/cellSize)
		j_min = round((processing_extent.XMin - dem_extent.XMin)/cellSize)
		j_max = round((processing_extent.XMax - dem_extent.XMin)/cellSize)

		# Temporarly changes the processing extent to use the whole DEM extent
		arcpy.env.extent = dem_extent
		
		if verbose:
			str_message = 'i: {}-{}, j: {}-{}'.format(i_min,i_max,j_min,j_max)
			arcpy.AddMessage(str_message)
			str_message = 'lower left: {}, {}'.format(dem_extent.XMin,dem_extent.YMin)
			arcpy.AddMessage(str_message)
			str_message = 'ncols: {}, nrows: {}'.format(ncols,nrows)
			arcpy.AddMessage(str_message)
		
		# Saves the elevation file (if selected)
		if grid_mnt_out== 'true' or grid_hill_out == 'true':
			str_message = 'Saving elevation file...'
			arcpy.AddMessage(str_message)
			grid_dem_after = arcpy.RasterToNumPyArray (grid_dem_file,lowerLeft,ncols,nrows,nodata_to_value=-9999)
			grid_dem_after = grid_dem_after.astype(np.float32)
			grid_dem_after[i_min:i_max,j_min:j_max]=grid_slbl
			if verbose:
				str_message = 'grid_dem_after shape: {}'.format(grid_dem_after.shape)
				arcpy.AddMessage(str_message)
				str_message = 'lowerleft: {}'.format(lowerLeft)
				arcpy.AddMessage(str_message)
				str_message = 'grid slbl shape: {}'.format(grid_slbl.shape)
				arcpy.AddMessage(str_message)
		
			grid_dem_out = arcpy.NumPyArrayToRaster(grid_dem_after,lowerLeft,cellSize, value_to_nodata=-9999)
					
			if grid_mnt_out== 'true':
				grid_dem_out.save(os.path.join(ws,out_basename+ext))
				if verbose:
					str_message = 'grid_dem_out saved with extent {}'.format(grid_dem_out.extent)
					arcpy.AddMessage(str_message)
			
		# Saves the thickness file (if selected)
		if grid_diff_out== 'true':
			str_message = 'Saving thickness file...'
			arcpy.AddMessage(str_message)
			grid_thick_after=np.zeros(grid_dem_after.shape)
			grid_thick_after[i_min:i_max,j_min:j_max]=grid_thickn
			grid_out = arcpy.NumPyArrayToRaster(grid_thick_after,lowerLeft,cellSize)
			
			if not (saveInGDB or len(ext) > 0):
				grid_out_thick_file = os.path.join(ws, out_basename[0:11] + '_t')
			else:
				grid_out_thick_file = os.path.join(ws,out_basename + '_t' + ext)
			grid_out.save(grid_out_thick_file)
			grid_out_thick_file_name=os.path.basename(grid_out_thick_file)

	elif extent_mode == 'Clip around the polygon(s)':
		lowerLeft= arcpy.Point(arcpy.env.extent.XMin,arcpy.env.extent.YMin)
		
		# Converts the slbl grid to a raster and save it if selected
		if grid_mnt_out== 'true' or grid_hill_out == 'true':
			grid_dem_out = arcpy.NumPyArrayToRaster(grid_slbl,lowerLeft,cellSize, value_to_nodata=-9999)
			if grid_mnt_out== 'true':
				str_message = 'Saving elevation file...'
				arcpy.AddMessage(str_message)
				if verbose:
					str_message = u'name is {} with type {}, extension is {} with type {}'.format(out_basename,type(out_basename).__name__,ext,type(out_basename).__name__)
					arcpy.AddMessage(str_message)
				grid_dem_out.save(os.path.join(ws,out_basename + ext))
			
		# Saves the thickness file (if selected)
		if grid_diff_out== 'true':
			str_message = 'Saving thickness file...'
			arcpy.AddMessage(str_message)
			grid_out = arcpy.NumPyArrayToRaster(grid_thickn,lowerLeft,cellSize)
			if not (saveInGDB or len(ext) > 0):
				grid_out_thick_file_name = out_basename[0:11] + '_t'
			else:
				grid_out_thick_file_name = out_basename + '_t'
			grid_out_thick_file_name = validatename(ws,grid_out_thick_file_name,ext)
			grid_out_thick_file = os.path.join(ws, grid_out_thick_file_name + ext)
			grid_out.save(grid_out_thick_file)
			grid_out_thick_file_name=os.path.basename(grid_out_thick_file)
			
	# Add the Thickness file to the project
	if grid_diff_out== 'true':
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
		if not (saveInGDB or len(ext) > 0):
			grid_out_hill_file_name = out_basename[0:11] + '_h'
		else:
			grid_out_hill_file_name = out_basename + '_hshd'
		grid_out_hill_file_name = validatename(ws,grid_out_hill_file_name,ext)
		grid_out_hill_file = os.path.join(ws, grid_out_hill_file_name + ext)
		grid_out_hill_lyr = arcpy.sa.Hillshade(grid_dem_out,315, 45, 'NO_SHADOWS', 1)
		grid_out_hill_lyr.save(grid_out_hill_file)
		
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
		try:
			del grid_out_hill_file, dataFrame
		except:
			pass
		
	# Return to the processing extent in case the processing continues with other tolerances
	arcpy.env.extent = processing_extent
		
def validatename(ws,out_basename,ext):
	# Validation of the outputs name
	out_basename = str(out_basename)
	ext = str(ext)
	if ext != '.png':
		out_basename = arcpy.ValidateTableName(out_basename+ext,ws)
		if len(ext)>0:
			out_basename = out_basename[0:-len(ext)]
	if verbose:
		arcpy.AddMessage(u'basename = {}, with extension = {}'.format(out_basename,ext))
	if arcpy.Exists(os.path.join(ws,out_basename+ext)):
		i = 1
		out_basename_temp = out_basename + '_{}'.format(i)
		if verbose:
			arcpy.AddMessage(u'attempted basename = {}, with extension = {}'.format(out_basename_temp,ext))
		while arcpy.Exists(os.path.join(ws,out_basename_temp+ext)):
			i+=1
			out_basename_temp = out_basename + '_{}'.format(i)
			if verbose:
				arcpy.AddMessage(u'attempted basename = {}, with extension = {}'.format(out_basename_temp,ext))
		out_basename = out_basename_temp
	return out_basename

def limiting_planes(point_file,grid_dem,processing_extent,cellSize):
	
	try:
		arcpy.RecalculateFeatureClassExtent_management(point_file)
	except:
		pass
	point_desc = arcpy.Describe(point_file)
	try:
		# Retrieve the selected points (works if the input is a layer)
		point_set = point_desc.FIDSet
	except:
		#makes a layer first if the input is a file
		if arcpy.GetInstallInfo()['ProductName'] == 'Desktop':
			point_file = arcpy.mapping.Layer(point_file)
		else:
			#mask_file = arcpy.mp.Layer(mask_file)
			point_file = arcpy.MakeFeatureLayer_management(point_file)
		point_desc = arcpy.Describe(point_file)
		# Retrieve the selected points
		point_set = point_desc.FIDSet
	if len(point_set) > 0:
		point_set = point_set.split(';')
		point_set = [int(x) for x in point_set]

	if point_desc.shapeType != "Point":
		arcpy.AddError('The limiting planes layer must be a point layer')
	else:
		# Check if the points and DEM are in the same SRC. Reproject the points if necessary
		src_point = point_desc.extent.spatialReference.name
		src_dem = arcpy.env.outputCoordinateSystem.name
		if verbose:
			arcpy.AddMessage('SRC point: {}, SRC dem: {}, src are equal: {}'.format(src_point,src_dem,src_point == src_dem))
		if src_point != src_dem:
			# Generate a random name for intermediate data (to avoid conflict if previous intermediate data weren't correctly deleted)
			lower_alphabet = string.ascii_lowercase
			random_string =''.join(random.choice(lower_alphabet) for i in range(5))
			name_temp = 'point_temp_' + random_string
			if saveInGDB:
				point_file_temp = os.path.join(ws,name_temp)
			else:
				point_file_temp = os.path.join(ws,name_temp + '.shp')
			arcpy.management.Project(point_file, point_file_temp, dem_desc.spatialReference)
			arcpy.AddMessage('Point file reprojected in the spatial reference of the DEM')
			point_file = point_file_temp
			if arcpy.GetInstallInfo()['ProductName'] == 'Desktop':
				point_file = arcpy.mapping.Layer(point_file)
			else:
				point_file = arcpy.MakeFeatureLayer_management(point_file)
			point_desc = arcpy.Describe(point_file)
			
		oid_fieldname = point_desc.OIDFieldName
	
	#create empty plane matrix
	s2 = grid_dem.shape
	planes = np.ones(s2)
	planes = -np.inf * planes


	dip_fieldname = None
	dipdir_fieldname = None
	point_file_fields = [f.name for f in arcpy.ListFields(point_file)]
	dip_name_allowed = ['dip_angle','dip angle','angle','slope','dip']
	dipdir_name_allowed = ['dir','direction','azimut','azi','dipdir','dip_dir','dip_direction','dip direction']
	for name in point_file_fields:
		if name.lower() in dip_name_allowed:
			dip_fieldname = name
	for name in point_file_fields:
		if name.lower() in dipdir_name_allowed:
			dipdir_fieldname = name
	
	if dip_fieldname is None:
		arcpy.AddError('Dip field not found. Available fields {}. Accepted fieldnames: {}'.format(', '.join(point_file_fields),', '.join(dip_name_allowed)))
	if dipdir_fieldname is None:
		arcpy.AddError('Dip direction field not found. Available fields {}. Accepted fieldnames: {}'.format(', '.join(point_file_fields),', '.join(dipdir_name_allowed)))
	
	if verbose:
		arcpy.AddMessage('Dip field name: {}'.format(dip_fieldname))
		arcpy.AddMessage('Dip direction field name: {}'.format(dipdir_fieldname))
	
	fields = [oid_fieldname,'SHAPE@',dipdir_fieldname,dip_fieldname]
	
	rows = arcpy.da.SearchCursor(point_file,fields)
	if len(point_set) > 0:
		for row in rows:
			feat = row[0]
			if int(feat) in point_set:
				plane_temp = calculate_plane(row,s2,point_desc.hasZ,processing_extent,cellSize)
				planes[plane_temp>planes]=plane_temp[plane_temp>planes]
			else:
				pass
	else: #no selected points --> takes all
		for row in rows:
			plane_temp = calculate_plane(row,s2,point_desc.hasZ,processing_extent,cellSize)
			planes[plane_temp>planes]=plane_temp[plane_temp>planes]

	planes[planes==-np.inf]=np.nan
	planes = np.flipud(planes)
	
	arcpy.AddMessage('Planes constaint max={}, min={}, average={}'.format(np.max(planes),np.min(planes),np.mean(planes)))
	
	return planes
	
def calculate_plane(row,shape,hasZ,processing_extent,cellSize):
	
	for pnt in row[1]: #assuming no multipart...
		if hasZ:
			xy = np.array([[float(pnt.X), float(pnt.Y),float(pnt.Z)]])
		else:
			xy = np.array([[float(pnt.X), float(pnt.Y)]])
	x_pt = xy[0,0]-processing_extent.XMin
	y_pt = xy[0,1]-processing_extent.YMin
	if hasZ:
		z_pt = xy[0,2]
	else:
		grid_y, grid_x = np.mgrid[0:shape[0]*cellSize:cellSize,0:shape[1]*cellSize:cellSize]
		grid_x = grid_x.flatten()
		grid_y = np.flipud(grid_y)
		grid_y = grid_y.flatten()
		values = grid_dem.flatten()
		grid_xy = np.dstack((grid_x, grid_y))
		grid_xy = grid_xy.reshape(len(values),2)
		values = values.reshape(len(values),1)
		z_pt = griddata(grid_xy, values, (x_pt,y_pt), method='linear')
		if verbose:
			arcpy.AddMessage('x: {},y: {},z: {}'.format(x_pt,y_pt,z_pt))
	
	dip = np.radians(row[3])
	dipdir = np.radians(row[2])
	
	# normal vector (see: https://ocw.snu.ac.kr/sites/default/files/NOTE/4435.pdf)
	nx = np.sin(dip) * np.sin(dipdir)
	ny = np.sin(dip) * np.cos(dipdir)
	nz = np.cos(dip)
		
	a = -nx/nz
	b = -ny/nz
	c = (nx*x_pt+ny*y_pt+nz*z_pt)/nz
	
	#nx(x-x_pt) + ny(y-y_pt) + nz(z-z_pt) = 0
	
	grid_y, grid_x = np.mgrid[0:shape[0]*cellSize:cellSize,0:shape[1]*cellSize:cellSize]
	
	plane = a*grid_x + b*grid_y + c
	
	return plane

#-----------------------------------------------------------------------------
# Main script starts here

if __name__=="__main__":
	# Retrive the parameters from the GUI. There are two types of GUI:
	# Multiple areas or single area

	grid_dem_file = arcpy.GetParameterAsText(0)
	mask_file = arcpy.GetParameterAsText(1)
	tol_mode = arcpy.GetParameterAsText(2) # single value, auto, auto min,inter,max
	if tol_mode == 'Single value':
		tols = [float(arcpy.GetParameterAsText(3).replace(',','.'))]
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
	
	if os.path.isdir(arcpy.GetParameterAsText(10)): #gdb return false...
		# Multiple areas
		merge_poly = False
		ws = arcpy.GetParameterAsText(10)
		name_field = arcpy.GetParameterAsText(11)
		extent_mode = 'Clip around the polygon(s)'
		grid_mnt_out = arcpy.GetParameterAsText(12)
		grid_diff_out = arcpy.GetParameterAsText(13)
		grid_hill_out = arcpy.GetParameterAsText(14)
	else:
		# Single area
		merge_poly = True
		point_file = arcpy.GetParameterAsText(10)
		if len(point_file)==0:
			plane_constraint = False
		else:
			plane_constraint = True
		ws = os.path.dirname(arcpy.GetParameterAsText(11))
		out_basename = os.path.basename(arcpy.GetParameterAsText(11))
		extent_mode = arcpy.GetParameterAsText(12) # full extent of the DEM, Clip around the polygon(s)
		grid_mnt_out = 'true'
		grid_diff_out = arcpy.GetParameterAsText(13)
		grid_hill_out = arcpy.GetParameterAsText(14)
	
	#Debbuging tools
	verbose = False
	savefigs = False
	figspath = r'D:\Soft\pySLBL\Unit_tests'
	
	# Check the necessary extensions
	if arcpy.CheckOutExtension("Spatial")!="CheckedOut":
		arcpy.AddMessage("A Spatial analyst license is needed to run this tool")

	# Check if the results are saved in a database
	desc = arcpy.Describe(ws)
	if desc.workspaceType == 'LocalDatabase':
		saveInGDB = True
	else:
		saveInGDB = False
	
	# Define the summary table (created later if necessary)
	if saveInGDB:
		summaryTableName = 'SLBL_results'
	else:
		summaryTableName = 'SLBL_results.dbf'
		
	summaryTable = os.path.join(ws,summaryTableName)
	
	arcpy.env.workspace = ws

	grid_dem_file = arcpy.Raster(grid_dem_file)
	#grid_dem_lyr = os.path.basename(grid_dem_file)
	#arcpy.MakeRasterLayer_management (grid_dem_file, grid_dem_lyr, "", "", "1")

	# Convert the polygon features to a raster mask
	try:
		arcpy.RecalculateFeatureClassExtent_management(mask_file)
	except:
		pass
	mask_desc = arcpy.Describe(mask_file)
	try:
		# Retrieve the selected polygons (works if the input is a layer)
		Set = mask_desc.FIDSet
	except:
		#makes a layer first if the input is a file
		if arcpy.GetInstallInfo()['ProductName'] == 'Desktop':
			mask_file = arcpy.mapping.Layer(mask_file)
		else:
			#mask_file = arcpy.mp.Layer(mask_file)
			mask_file = arcpy.MakeFeatureLayer_management(mask_file)
		mask_desc = arcpy.Describe(mask_file)
		# Retrieve the selected polygons
		Set = mask_desc.FIDSet
	if len(Set) > 0:
		Set = Set.split(';')
		Set = [int(x) for x in Set]

	if mask_desc.shapeType != "Polygon":
		arcpy.AddError('The mask layer must be a polygon layer')
	else:
		# Retrive key parameters from the DEM that will be used for the raster mask
		# Desc is done on the referenced raster since the extent is otherwise 
		# returned in the SRC of the dataframe
		grid_dem_path = grid_dem_file.catalogPath
		dem_desc = arcpy.Describe(grid_dem_path)
		arcpy.env.outputCoordinateSystem = dem_desc.spatialReference
		cellSize = dem_desc.meanCellWidth
		# With mosaic datasets, the returned cellSize might be 0. --> Try to get
		# the cellsize from the first raster of the mosaic dataset
		if cellSize == 0:
			try:
				dem_desc_temp = arcpy.Describe(os.path.join(grid_dem_path,"raster.objectid=1"))
				cellSize = dem_desc_temp.meanCellWidth
			except:
				arcpy.AddError("Can't get the cell size of the DEM...")
		dem_extent = dem_desc.extent
		
		# Check if the mask and DEM are in the same SRC. Reproject the mask if necessary
		src_mask = mask_desc.extent.spatialReference.name
		src_dem = dem_desc.extent.spatialReference.name
		if verbose:
			arcpy.AddMessage('SRC mask: {}, SRC dem: {}, src are equal: {}'.format(src_mask,src_dem,src_mask == src_dem))
		if src_mask != src_dem:
			# Generate a random name for intermediate data (to avoid conflict if previous intermediate data weren't correctly deleted)
			lower_alphabet = string.ascii_lowercase
			random_string =''.join(random.choice(lower_alphabet) for i in range(5))
			name_temp = 'mask_temp_' + random_string
			if saveInGDB:
				mask_file_temp = os.path.join(ws,name_temp)
			else:
				mask_file_temp = os.path.join(ws,name_temp + '.shp')
			arcpy.management.Project(mask_file, mask_file_temp, dem_desc.spatialReference)
			arcpy.AddMessage('Mask file reprojected in the spatial reference of the DEM')
			mask_file = mask_file_temp
			if arcpy.GetInstallInfo()['ProductName'] == 'Desktop':
				mask_file = arcpy.mapping.Layer(mask_file)
			else:
				mask_file = arcpy.MakeFeatureLayer_management(mask_file)
			mask_desc = arcpy.Describe(mask_file)
			mask_reproject = True
		else:
			mask_reproject = False
		
		oid_fieldname = mask_desc.OIDFieldName

		if merge_poly == True:
			# Get the name and extension of the output
			if out_basename.find('.') != -1:
				ext='.' + out_basename.split('.')[-1]
				out_basename = '.'.join(os.path.basename(out_basename).split('.')[0:-1])
			else:
				ext=''
			fields = ['OID@','SHAPE@']
			rows = arcpy.da.SearchCursor(mask_file, fields)
			extents = []
			listValue = []
			if len(Set) > 0:
				for row in rows:
					feat = row[0]
					if feat in Set:
						extents.append(row[1].extent)
						listValue.append([feat,1])
					else:
						listValue.append([feat,0])
				listValue.append(['NoData',0]) #must be at the end of the list
				# Find the extent of the whole set of polygon. Inspired from:
				# https://community.esri.com/t5/python-questions/get-extent-object-for-features-returned-by-searchcursor/td-p/724402
				xmin = min([extent.XMin for extent in extents])
				xmax = max([extent.XMax for extent in extents])
				ymin = min([extent.YMin for extent in extents])
				ymax = max([extent.YMax for extent in extents])
				poly_extent = arcpy.Extent(xmin, ymin, xmax, ymax)
				processing_extent = define_extent(poly_extent,dem_extent,cellSize,src_mask = extents[0].spatialReference)
				arcpy.env.extent = processing_extent
				grid_mask, grid_dem, z_min = raster2numpy(ws,ext,mask_file,mask_desc,not_deepen,listValue)
				if plane_constraint:
					grid_planes = limiting_planes(point_file,grid_dem,processing_extent,cellSize)
				if tol_mode == 'Auto' or tol_mode == 'Auto min/inter/max':
					tols = definetolerance(grid_mask, grid_dem)
				if tol_mode == 'Auto':
					tols = [tols[0]]
				tol_nr = 0
				for tol in tols:
					if plane_constraint:
						grid_slbl, grid_thickn, nb_iter = SLBL(grid_dem,grid_mask,tol,maxt,maxv,z_min,planes=grid_planes)
					else:
						grid_slbl, grid_thickn, nb_iter = SLBL(grid_dem,grid_mask,tol,maxt,maxv,z_min)
					if tol_mode == 'Auto min/inter/max':
						if tol_nr == 0:
							out_basename_validated = validatename(ws,out_basename + '_inter',ext)
						if tol_nr == 1:
							out_basename_validated = validatename(ws,out_basename + '_min',ext)
						if tol_nr == 2:
							out_basename_validated = validatename(ws,out_basename + '_max',ext)
					else:
						out_basename_validated = validatename(ws,out_basename,ext)
					fillsummarytable(summaryTable,out_basename,out_basename_validated,grid_thickn,grid_mask,nb_iter,tol)
					if grid_mnt_out== 'true' or grid_hill_out == 'true' or grid_diff_out == 'true':
						savegrids(grid_slbl,grid_thickn,ws,out_basename_validated,ext)
					tol_nr += 1
			else:
				#no selected polygons --> takes all
				for row in rows:
					feat = row[0]
					listValue.append([feat,1])
				listValue.append(['NoData',0]) #must be at the end of the list
				poly_extent = mask_desc.extent
				processing_extent = define_extent(poly_extent,dem_extent,cellSize,src_mask = poly_extent.spatialReference)
				arcpy.env.extent = processing_extent
				grid_mask, grid_dem, z_min = raster2numpy(ws,ext,mask_file,mask_desc,not_deepen,listValue)
				if plane_constraint:
					grid_planes = limiting_planes(point_file,grid_dem,processing_extent,cellSize)
				if tol_mode == 'Auto' or tol_mode == 'Auto min/inter/max':
					tols = definetolerance(grid_mask, grid_dem)
				if tol_mode == 'Auto':
					tols = [tols[0]]
				tol_nr = 0
				for tol in tols:
					if plane_constraint:
						grid_slbl, grid_thickn, nb_iter = SLBL(grid_dem,grid_mask,tol,maxt,maxv,z_min,planes=grid_planes)
					else:
						grid_slbl, grid_thickn, nb_iter = SLBL(grid_dem,grid_mask,tol,maxt,maxv,z_min)
					if tol_mode == 'Auto min/inter/max':
						if tol_nr == 0:
							out_basename_validated = validatename(ws,out_basename + '_inter',ext)
						if tol_nr == 1:
							out_basename_validated = validatename(ws,out_basename + '_min',ext)
						if tol_nr == 2:
							out_basename_validated = validatename(ws,out_basename + '_max',ext)
					else:
						out_basename_validated = validatename(ws,out_basename,ext)
					fillsummarytable(summaryTable,out_basename,out_basename_validated,grid_thickn,grid_mask,nb_iter,tol)
					if grid_mnt_out== 'true' or grid_hill_out == 'true' or grid_diff_out == 'true':
						savegrids(grid_slbl,grid_thickn,ws,out_basename_validated,ext)
					tol_nr += 1
		else:
			if saveInGDB:
				ext=''
			else:
				ext='.tif'
			if len(Set) > 0:
				fields = ['OID@','SHAPE@',name_field]
				rows = arcpy.da.SearchCursor(mask_file, fields)
				for row in rows:
					feat = row[0]
					if feat in Set:
						str_message = 'Calculating SLBL for polygon {}'.format(row[2])
						arcpy.AddMessage(str_message)
						listValue = [[str(feat),1]]
						listValue.append(['NoData',0]) #must be at the end of the list
						poly_extent = row[1].extent
						processing_extent = define_extent(poly_extent,dem_extent,cellSize,src_mask = poly_extent.spatialReference)
						arcpy.env.extent = processing_extent
						clause = '"' + oid_fieldname + '" = ' + str(row[0])
						arcpy.AddMessage('Selecting feature ' +  str(row[0]))
						# Select the current polygon to avoid fail classification in the case of overlapping polygons
						arcpy.SelectLayerByAttribute_management(mask_file, "NEW_SELECTION", clause)
						grid_mask, grid_dem, z_min = raster2numpy(ws,ext,mask_file,mask_desc,not_deepen,listValue)
						if tol_mode == 'Auto' or tol_mode == 'Auto min/inter/max':
							tols = definetolerance(grid_mask, grid_dem)
						if tol_mode == 'Auto':
							tols = [tols[0]]
						tol_nr = 0
						for tol in tols:
							grid_slbl, grid_thickn, nb_iter = SLBL(grid_dem,grid_mask,tol,maxt,maxv,z_min)
							if tol_mode == 'Auto min/inter/max':
								if tol_nr == 0:
									out_basename_validated = validatename(ws,row[2] + '_inter',ext)
								if tol_nr == 1:
									out_basename_validated = validatename(ws,row[2] + '_min',ext)
								if tol_nr == 2:
									out_basename_validated = validatename(ws,row[2] + '_max',ext)
							else:
								out_basename_validated = validatename(ws,row[2],ext)
							str_message = u'Original name: {}, validated name: {}'.format(row[2],out_basename_validated)
							arcpy.AddMessage(str_message)
							fillsummarytable(summaryTable,row[2],out_basename_validated,grid_thickn,grid_mask,nb_iter,tol)
							if grid_mnt_out== 'true' or grid_hill_out == 'true' or grid_diff_out == 'true':
								savegrids(grid_slbl,grid_thickn,ws,out_basename_validated,ext)
							tol_nr += 1
				# Reselect the polygons that where selected before running the tool:
				clause = ''
				for feat in Set:
					if len(clause) > 0:
						clause = clause + ' OR "' + oid_fieldname + '" = ' + str(feat)
					else:
						clause = '"' + oid_fieldname + '" = ' + str(feat)
				arcpy.AddMessage(clause)
				arcpy.SelectLayerByAttribute_management(mask_file, "NEW_SELECTION", clause)
			else:
				#no selected polygons --> takes all
				fields = ['OID@','SHAPE@',name_field]
				rows = arcpy.da.SearchCursor(mask_file, fields)
				for row in rows:
					str_message = 'Calculating SLBL for polygon {}'.format(row[2])
					arcpy.AddMessage(str_message)
					feat = row[0]
					listValue = [[str(feat),1]]
					listValue.append(['NoData',0]) #must be at the end of the list
					poly_extent = row[1].extent
					processing_extent = define_extent(poly_extent,dem_extent,cellSize,src_mask = poly_extent.spatialReference)
					arcpy.env.extent = processing_extent
					clause = '"' + oid_fieldname + '" = ' + str(row[0])
					arcpy.AddMessage('Selecting feature ' +  str(row[0]))
					# Select the current polygon to avoid fail classification in the case of overlapping polygons
					arcpy.SelectLayerByAttribute_management(mask_file, "NEW_SELECTION", clause)
					grid_mask, grid_dem, z_min = raster2numpy(ws,ext,mask_file,mask_desc,not_deepen,listValue)
					if tol_mode == 'Auto' or tol_mode == 'Auto min/inter/max':
						tols = definetolerance(grid_mask, grid_dem)
					if tol_mode == 'Auto':
						tols = [tols[0]]
					tol_nr = 0
					for tol in tols:
						grid_slbl, grid_thickn, nb_iter = SLBL(grid_dem,grid_mask,tol,maxt,maxv,z_min)
						if tol_mode == 'Auto min/inter/max':
							if tol_nr == 0:
								out_basename_validated = validatename(ws,row[2] + '_inter',ext)
							if tol_nr == 1:
								out_basename_validated = validatename(ws,row[2] + '_min',ext)
							if tol_nr == 2:
								out_basename_validated = validatename(ws,row[2] + '_max',ext)
						else:
							out_basename_validated = validatename(ws,row[2],ext)
						str_message = u'Original name: {}, validated name: {}'.format(row[2],out_basename_validated)
						arcpy.AddMessage(str_message)
						fillsummarytable(summaryTable,row[2],out_basename_validated,grid_thickn,grid_mask,nb_iter,tol)
						if grid_mnt_out== 'true' or grid_hill_out == 'true' or grid_diff_out == 'true':
							savegrids(grid_slbl,grid_thickn,ws,out_basename_validated,ext)
						tol_nr += 1
				# clears the selection (returns to the original state)
				arcpy.SelectLayerByAttribute_management(mask_file, "CLEAR_SELECTION")

	if mask_reproject:
		arcpy.Delete_management(mask_file_temp)

	try:
		del grid_dem_out, grid_dem_after, cur, row
	except:
		pass