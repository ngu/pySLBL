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
    self.params[11].filter.list=lst
    self.params[11].value = 'Full extent of the DEM'

    return

  def updateParameters(self):
    """Modify the values and properties of parameters before internal
    validation is performed.  This method is called whenever a parameter
    has been changed."""


  def updateMessages(self):
    """Modify the messages created by internal validation for each tool
    parameter.  This method is called after internal validation."""
    if self.params[10].value:
      if arcpy.Exists(self.params[10].value.value):
        self.params[10].setErrorMessage("The feature class already exists")
      else:
        desc = arcpy.Describe(os.path.dirname(self.params[10].value.value))
        name = os.path.basename(self.params[10].value.value)
        if desc.workspaceType == 'LocalDatabase':
          self.params[10].clearMessage()
        else:
          if name.find('.') == -1:
            if len(name) > 13:
              self.params[10].setErrorMessage("The file name of an ESRI grid is limited to 13 characters. It is recommended to save in tiff (add '.tif' at the end of the name) or in a database, but you may also just change the name")
            else:
              self.params[10].setWarningMessage("It is recommended to save in tiff (add '.tif' at the end of the name) or in a database")
          else:
            self.params[10].clearMessage()

    
    if self.params[2].value == 'Single value':
      self.params[3].enabled = 1
      self.params[3].parameterType = 'Required'
    else:
      self.params[3].enabled = 0
      self.params[3].value = ''
      self.params[3].parameterType = 'Optional'
      self.params[3].clearMessage()
    return