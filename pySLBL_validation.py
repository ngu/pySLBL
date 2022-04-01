# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 14:14:33 2020

@author: Nicolet_Pierrick
"""

import arcpy
import os

class ToolValidator(object):
  """Class for validating a tool's parameter values and controlling
  the behavior of the tool's dialog."""

  def __init__(self):
    """Setup arcpy and the list of tool parameters."""
    self.params = arcpy.GetParameterInfo()

  def initializeParameters(self):
    """Refine the properties of a tool's parameters.  This method is
    called when the tool is opened."""
    lst=['4 neighbours, average','8 neighbours, average', '4 neighbours, min/max','8 neighbours, min/max']
    self.params[7].filter.list=lst
    self.params[7].value = '4 neighbours, average'

    lst=['Single value','Auto','Auto min/inter/max']
    self.params[2].filter.list=lst
    self.params[2].value = 'Single value'

    lst=['Full extent of the DEM','Clip around the polygon(s)']
    self.params[12].filter.list=lst
    self.params[12].value = 'Full extent of the DEM'

    return

  def updateParameters(self):
    """Modify the values and properties of parameters before internal
    validation is performed.  This method is called whenever a parameter
    has been changed."""


  def updateMessages(self):
    """Modify the messages created by internal validation for each tool
    parameter.  This method is called after internal validation."""
    if self.params[11].value:
      if arcpy.Exists(self.params[11].value.value):
        self.params[11].setErrorMessage("The feature class already exists")
      else:
        desc = arcpy.Describe(os.path.dirname(self.params[11].value.value))
        name = os.path.basename(self.params[11].value.value)
        if desc.workspaceType == 'LocalDatabase':
          self.params[11].clearMessage()
        else:
          if name.find('.') == -1:
            if len(name) > 13:
              self.params[11].setErrorMessage("The file name of an ESRI grid is limited to 13 characters. It is recommended to save in tiff (add '.tif' at the end of the name) or in a database, but you may also just change the name")
            else:
              self.params[11].setWarningMessage("It is recommended to save in tiff (add '.tif' at the end of the name) or in a database")
          else:
            self.params[11].clearMessage()

    
    if self.params[2].value == 'Single value':
      self.params[3].enabled = 1
      self.params[3].parameterType = 'Required'
    else:
      self.params[3].enabled = 0
      self.params[3].value = ''
      self.params[3].parameterType = 'Optional'
      self.params[3].clearMessage()

    if self.params[10].value:
      point_desc = arcpy.Describe(self.params[10].value)
      if point_desc.shapeType != "Point":
        self.params[10].setErrorMessage('The limiting planes layer must be a point layer')
      else:
        dip_fieldname = None
        dipdir_fieldname = None
        point_file_fields = [f.name for f in arcpy.ListFields(self.params[10].value)]
        dip_name_allowed = ['dip_angle','dip angle','angle','slope','dip']
        dip_name_allowed.reverse()
        dipdir_name_allowed = ['dir','direction','azimut','azi','dipdir','dip_dir','dip_direction','dip direction']
        dipdir_name_allowed.reverse()
        for name in point_file_fields:
          if name.lower() in dip_name_allowed:
            dip_fieldname = name
        for name in point_file_fields:
          if name.lower() in dipdir_name_allowed:
            dipdir_fieldname = name
        if dip_fieldname is None:
          self.params[10].setErrorMessage('Dip field not found. Accepted fieldnames: {}'.format(', '.join(dip_name_allowed)))
        elif dipdir_fieldname is None:
          self.params[10].setErrorMessage('Dip direction field not found. Accepted fieldnames: {}'.format(', '.join(dipdir_name_allowed)))
        else:
          self.params[10].clearMessage()
    else:
      self.params[10].clearMessage()

    if self.params[0].value and self.params[12].value == 'Full extent of the DEM':
      DEM = self.params[0].value
      desc = arcpy.Describe(DEM)
      xmin = desc.extent.XMin
      xmax = desc.extent.XMax
      ymin = desc.extent.YMin
      ymax = desc.extent.YMax
      res = desc.meanCellWidth
      ncol = (xmax - xmin)/res
      nrow = (ymax - ymin)/res
      if nrow*ncol > 10000000:
        self.params[12].setWarningMessage("The raster has {} rows and {} columns, which is quite large. It is recommanded to use the option 'Clip around the polygon(s)'".format(str(int(nrow)),str(int(ncol))))
      else:
        self.params[12].clearMessage()
    return