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
    self.params[6].filter.list=lst
    self.params[6].value = '4 neighbours, average'

    return

  def updateParameters(self):
    """Modify the values and properties of parameters before internal
    validation is performed.  This method is called whenever a parameter
    has been changed."""

  def updateMessages(self):
    """Modify the messages created by internal validation for each tool
    parameter.  This method is called after internal validation."""
    if self.params[2].value:
      desc = arcpy.Describe(os.path.dirname(self.params[2].value.value))
      name = os.path.basename(self.params[2].value.value)
      if desc.workspaceType == 'LocalDatabase':
        self.params[2].clearMessage()
      else:
        if name.find('.') == -1:
          if len(name) > 13:
            self.params[2].setErrorMessage("The file name of an ESRI grid is limited to 13 characters. Please change the name of the output, save in a database or in another format (for example tiff)")
          else:
            self.params[2].clearMessage()
        else:
          self.params[2].clearMessage()
    else:
      self.params[2].clearMessage()
    
    if self.params[0].value:
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
        self.params[0].setWarningMessage("The raster has {} rows and {} columns. As the mask will be converted to a raster using the same extent and resolution, it is recommended to use a smaller raster. The outputs will also use the same extent, while the SLBL will only be calculated on the masked area".format(str(int(nrow)),str(int(ncol))))
      else:
        self.params[0].clearMessage()
    return