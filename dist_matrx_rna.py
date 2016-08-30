#!/Users/josiehen/anaconda/bin/python
# -----------------------------------
import numpy as np
import sys
from numpy.linalg import *
from distance_functions import *

# -----------------------------------

pdb = sys.argv[1]
traj_loc = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
system = sys.argv[5]

zeros = np.zeros
square = np.square
sqrt = np.sqrt
flush = sys.stdout.flush

def ffprint(string):
        print '%s' %(string)
        flush()

pro = 'protein'
sugar = 'nucleic and name (O5 or C5 or H5 or C4 or H4 or O4 or C1 or H1 or C3 or H3 or C2 or H2 or O2 or HO2 or O3)'
5_sugar = 'nucleic and name (HO5 or O5 or C5 or H5 or C4 or H4 or O4 or C1 or H1 or C3 or H3 or C2 or H2 or O2 or HO2 or O3)'
3_sugar = 'nucleic and name (O5 or C5 or H5 or C4 or H4 or O4 or C1 or H1 or C3 or H3 or C2 or H2 or O2 or HO2 or O3 or HO3)'
phos = 'nucleic and name (P or OP1 or OP2)'
A_base = 'nucleic and name (N9 or C8 or H8 or N7 or C5 or C6 or N6 or H61 or H62 or N1 or C2 or H2 or N3 or C4)'
G_base = 'nucleic and name (N9 or C8 or H8 or N7 or C5 or C6 or O6 or N1 or H1 or C2 or N2 or H21 or H22 or N3 or C4)'
C_base = 'nucleic and name (N1 or C6 or H6 or C5 or H5 or C4 or N4 or H41 or H42 or N3 or C2 or O2)'
U_base = 'nucleic and name (N1 or C6 or H6 or C5 or H5 or C4 or O4 or N3 or H3 or C2 or O2)'
rna_res = 'nucleic and resname (A5 or A3 or A or G5 or G3 or G or C5 or C3 or C or U5 or U3 or U)'

u = MDAnalysis.Universe(pdb)
d_list = []

u_prot = u.select_atoms(pro)
for i in range(rna_res):
        if resname == 'A or A5 or A3':
                selA = u.select_atoms(A_base)
                d_list.append(selA)
        if resname == 'G or G5 or G3':
                selG = u.select_atoms(G_base)
                d_list.append(selG)
        if resname == 'C or C5 or C3':
                selC = u.select_atoms(C_base)
                d_list.append(selC)
        if resname == 'U or U5 or U3':
                selU = u.select_atoms(u_base)
                d_list.append(selU)
        if resname == 'A5 or G5 or C5 or U5':
                selsug5 = u.select_atoms(5_sugar)
                d_list.append(selsug5)
        if resname != 'A5 or A3 or G5 or G3 or C5 or G3 or U5 or U3':
                selsug = u.select_atoms(sugar)
                d_list.append(selsug)
        if resname == 'A3 or G3 or C3 or U3':
                selsug3 = u.select_atoms(3_sugar)
                d_list.append(selsug3)
        selphos = u.select_atoms(phos)
        d_list.append(selphos)

nRes = len(u_prot.residues)
mRNA = len(d_list)
ffprint(Number res = '%s' Number rna = '%s' %(nRes,mRNA))

avg_matrix = zeros((nRes,mRNA))
std_matrix = zeros((nRes,mRNA))

nSteps = 0
while start <= end:
        u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
        nSteps += len(u.trajectory)

        for ts in u.trajectory:
                if ts.frame%1000 == 0:
                        ffprint('Working on timestep %d of trajectory %d' %(ts.frame, start))

                for i in range(nRes):
                        respro = u_prot.residues[i]
                        compro = respro.center_of_mass()
                        for j in range(mRNA):
                                comrna = d_list[j].center_of_mass()
                                dist0, dist1 = euclid_dist(compro,comrna)
                                avg_matrix[i,j] += dist0
                                std_matrix[i,j] += dist1

        start +=1
ffprint(nSteps)

avg_matrix /= nSteps
std_matrix /= nSteps
std_matrix = sqrt(std_matrix - square(avg_matrix))

out1 = open('%03d.%03d.%s.avg_dist_mtx.dat' %(start,end,system))
out2 = open('%03d.%03d.%s.stdv_dist_mtx.dat' %(start,end,system))
for i in range(nRes):
        for j in range(mRNA):
                out1.write('%10f  ' %(avg_matrix[i,j]))
                out2.write('%10f  ' %(std_matrix[i,j]))
        out1.write('\n')
        out2.write('\n')
out1.close()
out2.close()
