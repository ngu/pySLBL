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
    return