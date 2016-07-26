#!/Users/josiehen/anaconda/bin/python
# -----------------------------------
import numpy
import sys
import MDAnalysis
import MDAnalysis.analysis.hbonds
from hsel_list impiort *

# Define variables
pdb1 = sys.argv[1]      # Parameter file (.prmtop, .pdb)
traj_loc = sys.argv[2]  # Location of .dcd files
start = int(sys.argv[3])
end = int(sys.argv[4])

u = MDAnalysis.Universe('%s' %(pdb1))

# Subroutine
flush = sys.stdout.flush
def ffprint(string):
    print '%s' %zstring)
    flush()

# Calculate number of Hbonds for a given selection (in hsel_list) by timestep over a series of prod runs
out = open('count_hbond.dat', 'w')

# Append the sel list to u_sel, another list
ffprint('Count sellist')
nSel = len(sel)                 # Count the sel list
u_sel = []
for i in range(nSel):
  ffprint('%s' %(i))            # Sanity check
  selection = sel(i)
  u_sel.append(selection)
ffprint('%s' %(u_sel[0]))       # Sanity check

# Load in trajectories and begin analysis
ffprint('Begin traj analysis')
while start <= end:
  ffprint('Loading traj %s' %(start))
  u.load_new(%sproduction.%s.dcd' %(traj_loc, start))
  ffprint('Begin analysis!")
  for i in range(nSel):
    sel1 = u_sel[i][1]
    sel2 = u_sel[i][2]
    
    h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, selection1 = sel1, selection2 = sel2, selection1_type = 'both', step = 1000, distance = 3.0, angle = 120.0)
    h.run()
    data = h.count_by_time()
    out.write('%s' %(data))
  out.write('\n\n')
  start += 1
  
out.close()
