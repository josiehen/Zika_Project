#!/Users/josiehen/anaconda/bin/python
# -----------------------------------
import numpy
import sys
import MDAnalysis
import MDAnalysis.analysis.hbonds
from hsel_list impiort *

# Define variables
pdb1 = sys.argv[1]           # Parameter file (.prmtop, .pdb)
traj_loc = sys.argv[2]       # Location of .dcd files
start = int(sys.argv[3])     
end = int(sys.argv[4]
system = sys.argv[5]         # System descriptor ie "zika_apo"

u = MDAnalysis.Universe('%s' %(pdb1))

# Subroutine
flush = sys.stdout.flush
def ffprint(string):
    print '%s' %zstring)
    flush()

# Calculate number of Hbonds for a given selection (in hsel_list) by timestep/frame over a series of prod runs

# Append the sel list to u_sel, another list
ffprint('Count sellist')
nSel = len(sel)                 # Count the sel list
u_sel = []                      # Create an empty list u_sel
for i in range(nSel):
  #ffprint('%s' %(i))            # Sanity check
  selection = sel(i)             # Each specified selection to be analyzed
  u_sel.append(selection)        # Appends information from sel_list to u_sel
#ffprint('%s' %(u_sel[0]))       # Sanity check

# -------------------------------------------------

# Load in trajectories and begin analysis
ffprint('Begin traj analysis')
while start <= end:
  ffprint('Loading traj %s' %(start))
  u.load_new(%sproduction.%s.dcd' %(traj_loc, start))     # Load a new trajectory into the universe
  ffprint('WORK!")
  for i in range(nSel):
    sel1 = u_sel[i][1]                  # First selection from sel_list
    sel2 = u_sel[i][2]                  # Second selection
    out = open('%s.%s.hbonds.dat' %(system,u_sel[i][0]),'a')        # This file, when done, will contain information from all analyzed production runs for a given selection
    h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, selection1 = sel1, selection2 = sel2, selection1_type = 'both', step = 1000, distance = 3.0, angle = 120.0)
    h.run()
    data = h.count_by_time()            # Program to count hydrogen bonds by timestep
    for i in data:
        #print i[1]                     # Sanity check, should give number H bonds in frame per frame
        bonds = i[1]                    # Select number hbonds only
        out.write('%s   ' %(bonds))
  start += 1
  
out.close()
