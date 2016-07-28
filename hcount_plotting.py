#!/Users/josiehen/anaconda/bin/python
# ----------------------------------
from plotting_functions import *
from hsel_list import *
import numpy as np
import sys

# Subroutine
flush = sys.stdout.flush
def ffprint(string):
  print '%s' %(string)
  flush()

# Define variables
data = sys.argv[1]        # .dat file produced by count_hbond.py
system = sys.argv[2]      # System descriptor ie "zika_apo"
selection = sys.argv[3]   # Selection descriptor ie "LB3B4"

# Main routine, create a scatterplot/histogram of hydrogen bond count over time for a selection
bondlist = np.loadtxt(data)

nSteps = len(bondlist[:[)      # Give number of steps, will equal number of frames for a step of 1
ffprint('Number of steps: '%s' %(nSteps))

time = np.zeros(nSteps)
for i in range(nSteps):
  time[i] = i*2.0              # For 1000 frame step, convert frame time to ns
#print time[:]                 # Sanity check, should give time per frame in ns
#ffprint('bondlist[:]')
#print bondlist[:]             # Sanity check, should return data counted across row
scat_hist(time[:],bondlist[:],'k','Time (ns)','Hydrogen Bond Count', '%s.%s' %(system,selection), 'Hydrogen_Bonding')
