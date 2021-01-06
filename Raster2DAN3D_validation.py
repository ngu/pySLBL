# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 14:14:33 2020

@author: Nicolet_Pierrick
"""

import arcpy

class ToolValidator(object):
  """Class for validating a tool's parameter values and controlling
  the behavior of the tool's dialog."""

  def __init__(self):
    """Setup arcpy and the list of tool parameters."""
    self.params = arcpy.GetParameterInfo()

  def initializeParameters(self):
    """Refine the properties of a tool's parameters.  This method is
    called when the tool is opened."""
    return

  def updateParameters(self):
    """Modify the values and properties of parameters before internal
    validation is performed.  This method is called whenever a parameter
    has been changed."""

  def updateMessages(self):
    """Modify the messages created by internal validation for each tool
    parameter.  This method is called after internal validation."""
    if self.params[1].value and self.params[0].value:
      DEM = self.params[0].value
      desc = arcpy.Describe(DEM)
      xmin_0 = desc.extent.XMin
      xmax_0 = desc.extent.XMax
      ymin_0 = desc.extent.YMin
      ymax_0 = desc.extent.YMax
      res_0 = desc.meanCellWidth

      DEPTH = self.params[1].value
      desc = arcpy.Describe(DEPTH)
      xmin_1 = desc.extent.XMin
      xmax_1 = desc.extent.XMax
      ymin_1 = desc.extent.YMin
      ymax_1 = desc.extent.YMax
      res_1 = desc.meanCellWidth
      
      if (xmin_0==xmin_1) and (xmax_0==xmax_1) and (ymin_0==ymin_1) and (ymax_0==ymax_1) and (res_0==res_1):
        self.params[1].clearMessage()
      else:
        self.params[1].setWarningMessage("The Source Depth raster should have the same extent and resolution as the DEM")
    else:
      self.params[1].clearMessage()

    if self.params[2].value and self.params[0].value:
      DEM = self.params[0].value
      desc = arcpy.Describe(DEM)
      xmin_0 = desc.extent.XMin
      xmax_0 = desc.extent.XMax
      ymin_0 = desc.extent.YMin
      ymax_0 = desc.extent.YMax
      res_0 = desc.meanCellWidth

      EROSION = self.params[2].value
      desc = arcpy.Describe(EROSION)
      xmin_1 = desc.extent.XMin
      xmax_1 = desc.extent.XMax
      ymin_1 = desc.extent.YMin
      ymax_1 = desc.extent.YMax
      res_1 = desc.meanCellWidth
      
      if (xmin_0==xmin_1) and (xmax_0==xmax_1) and (ymin_0==ymin_1) and (ymax_0==ymax_1) and (res_0==res_1):
        self.params[2].clearMessage()
      else:
        self.params[2].setWarningMessage("The Erosion raster should have the same extent and resolution as the DEM")
    else:
      self.params[2].clearMessage()

    return