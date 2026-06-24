'''
This script executes the Fokker 100 example
'''

#IMPORTS
import numpy as np
from designTool.standard_airplane import standard_airplane
from designTool.geometry import geometry
from designTool.plots import plot_geometry
from designTool.aerodynamics import aerodynamics
import matplotlib.pyplot as plt

#=========================================

# SETUP

# Constants
ft2m = 0.3048
kt2ms = 0.514444
lb2N = 4.44822
gravity = 9.81

# Select airplane name from the standard_airplane function in designTool
airplane_name = 'Tomav'
#airplane_name = 'my_airplane'

#=========================================

# EXECUTION

# Load the airplane dictionary
airplane = standard_airplane(airplane_name)

# Execute the geometry module to compute all dimensions.
# This updates the airplane dictionary with new entries.
geometry(airplane)

plot_geometry(airplane, az1=0, az2=180, figname="Vista Frontal")

plot_geometry(airplane, az1=90, az2=0, figname="Vista Superior")

plot_geometry(airplane, az1=0, az2=90, figname="Vista Lateral")

plt.show()