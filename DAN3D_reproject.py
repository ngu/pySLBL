# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# DAN3D_reproject.py
# Created on: 12.01.2021
# by: Pierrick.Nicolet@ngu.no
# ---------------------------------------------------------------------------
#
# This tool intends to shift back to their original position all the files
# saved by DAN3D. It scans the provided result folder and copies all the grids
# (“*.grd”), particles files (“parts*.txt”) as well as the summary files
# “finaloutput.txt” and “output.txt” to a specified output folder after
# shifting their coordinates back to their original position.
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


import os
import arcpy
import re

if __name__=="__main__":
	
	res_folder = arcpy.GetParameterAsText(0)
	header_file = arcpy.GetParameterAsText(1)
	output_folder = arcpy.GetParameterAsText(2)
	overwrite = arcpy.GetParameterAsText(3)
	
	with open(header_file,'r') as f:
		header_text = f.read()
	
	with open(header_file,'r') as f:
		header = f.readlines()
		x_shift = float(header[2].split(' ')[0])
		y_shift = float(header[3].split(' ')[0])
		width0 = float(header[1].split(' ')[0])
		height0 = float(header[1].split(' ')[1])
	
	os.walk(res_folder)
	
	file_list = os.listdir(res_folder)
	
	str_message = u'Reprojecting the files in folder {}...'.format(res_folder)
	arcpy.AddMessage(str_message)
		
	dim_warning = False

	for file in file_list:
		if file[-4:] == '.grd':
			source = os.path.join(res_folder, file)
			with open(source,"r") as f:
				if not dim_warning:
					header = f.readlines()
					width1 = float(re.split('\s|\t',header[1])[0])
					height1 = float(re.split('\s|\t',header[1])[1])
					if width1 != width0:
						str_message = 'The dimensions specified in the header file do not match the dimensions of the grids'
						arcpy.AddWarning(str_message)
						dim_warning = True
					elif height1 != height0:
						str_message = 'The dimensions specified in the header file do not match the dimensions of the grids'
						arcpy.AddWarning(str_message)
						dim_warning = True
			with open(source,"r") as f:
				lines = f.readlines()
				lines = lines[4:]
			target = os.path.join(output_folder, file)
			if overwrite != 'true':
				if not os.path.isfile(target):
					with open(target, "w") as f1:
						f1.writelines(header_text)
						f1.writelines(lines)
				else:
					str_message = u'Cannot overwrite file {} as overwrite is off'.format(target)
					arcpy.AddWarning(str_message)
			else:
				with open(target, "w") as f1:
					f1.writelines(header_text)
					f1.writelines(lines)
		if file[-4:] == '.txt':
			if file[:5] == 'parts':
				source = os.path.join(res_folder, file)
				with open(source,"r") as f:
					lines_corr = []
					lines = f.readlines()
					for line in lines:
						try:
							x0 = float(line.split('\t')[0])
							x1 = x0 + x_shift
							pos = line.find('\t')
							line_corr = str(x1) + line[pos:]
							y0 = float(line.split('\t')[1])
							y1 = y0 + y_shift
							pos = line_corr.find('\t')
							pos1 = line_corr.find('\t',pos+1)
							line_corr = line_corr[:pos+1] + str(y1) + line_corr[pos1:]
							lines_corr.append(line_corr)
						except:
							lines_corr.append(line)
				target = os.path.join(output_folder, file)
				if overwrite != 'true':
					if not os.path.isfile(target):
						with open(target, "w") as f1:
							f1.writelines(lines_corr)
					else:
						str_message = u'Cannot overwrite file {} as overwrite is off'.format(target)
						arcpy.AddWarning(str_message)
				else:
					with open(target, "w") as f1:
						f1.writelines(lines_corr)
			elif file[:11] == 'finaloutput':
				source = os.path.join(res_folder, file)
				with open(source,"r") as f:
					lines_corr = []
					lines = f.readlines()
					for line in lines:
						if line.split('\t')[0] == 'Final COM X-Position (m):':
							x0 = float(line.split('\t')[1])
							x1 = x0 + x_shift
							line_corr = line.replace(line.split('\t')[1],str(x1)+'\n')
							lines_corr.append(line_corr)
						elif line.split('\t')[0] == 'Final COM Y-Position (m):':
							y0 = float(line.split('\t')[1])
							y1 = y0 + y_shift
							line_corr = line.replace(line.split('\t')[1],str(y1)+'\n')
							lines_corr.append(line_corr)
						else:
							lines_corr.append(line)
				target = os.path.join(output_folder, file)
				if overwrite != 'true':
					if not os.path.isfile(target):
						with open(target, "w") as f1:
							f1.writelines(lines_corr)
					else:
						str_message = u'Cannot overwrite file {} as overwrite is off'.format(target)
						arcpy.AddWarning(str_message)
				else:
					with open(target, "w") as f1:
						f1.writelines(lines_corr)
			elif file[:6] == 'output':
				source = os.path.join(res_folder, file)
				with open(source,"r") as f:
					lines_corr = []
					lines = f.readlines()
					for line in lines:
						if line.split('\t')[0] == 'COM X-Position (m)=':
							x0 = float(line.split('\t')[1])
							x1 = x0 + x_shift
							line_corr = line.replace(line.split('\t')[1],str(x1)+'\n')
							lines_corr.append(line_corr)
						elif line.split('\t')[0] == 'COM Y-Position (m)=':
							y0 = float(line.split('\t')[1])
							y1 = y0 + y_shift
							line_corr = line.replace(line.split('\t')[1],str(y1)+'\n')
							lines_corr.append(line_corr)
						else:
							lines_corr.append(line)
				target = os.path.join(output_folder, file)
				if overwrite != 'true':
					if not os.path.isfile(target):
						with open(target, "w") as f1:
							f1.writelines(lines_corr)
					else:
						str_message = u'Cannot overwrite file {} as overwrite is off'.format(target)
						arcpy.AddWarning(str_message)
				else:
					with open(target, "w") as f1:
						f1.writelines(lines_corr)