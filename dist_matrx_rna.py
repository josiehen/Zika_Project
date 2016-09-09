#!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
# -----------------------------------
import MDAnalysis
import numpy as np
import sys
from numpy.linalg import *
from distance_functions import *

# -----------------------------------

pdb = sys.argv[1]               #Needs to be a .pdb with same number of atoms (used coordinates from frame 1, production 1)
traj_loc = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])

zeros = np.zeros
square = np.square
sqrt = np.sqrt
flush = sys.stdout.flush

def ffprint(string):
        print '%s' %(string)
        flush()

#Selections; shortened from saying 'or' inbetween each one
pro = 'protein'
sugar = "nucleic and name O5' H5'' C5' H5' C4' H4' O4' C1' H1' C3' H3' C2' H2' O2' HO2' O3' "
5_sugar = "nucleic and name HO5' H5'' O5' C5' H5' C4' H4' O4' C1' H1' C3' H3' C2' H2' O2' HO2' O3' "
3_sugar = "nucleic and name O5' C5' H5' H5'' C4' H4' O4' C1' H1' C3' H3' C2' H2' O2' HO2' O3' HO3' "
phos = 'nucleic and name P OP1 OP2'
A_base = 'nucleic and name N9 C8 H8 N7 C5 C6 N6 H61 H62 N1 C2 H2 N3 C4'
G_base = 'nucleic and name N9 C8 H8 N7 C5 C6 O6 N1 H1 C2 N2 H21 H22 N3 C4'
C_base = 'nucleic and name N1 C6 H6 C5 H5 C4 N4 H41 H42 N3 C2 O2'
U_base = 'nucleic and name N1 C6 H6 C5 H5 C4 O4 N3 H3 C2 O2'
rna = 'nucleic or resname A5 A3 A G5 G3 G C5 C3 C U5 U3 U'

u = MDAnalysis.Universe(pdb)

d_list = []
ressel = u.select_atoms(rna)
rna_res = ressel.n_residues      #Selects residues in nucleic - before was just a list
#Loop through RNA residues making selections. Changed a ton of 'if' statements to 'elif' or 'else'
for i in range(rna_res):
        resname = ressel.residues[i].resname
        if resname == 'A or A5 or A3':
                selA = u.select_atoms(A_base)
                d_list.append(selA)
        elif resname == 'G or G5 or G3':
                selG = u.select_atoms(G_base)
                d_list.append(selG)
        elif resname == 'C or C5 or C3':
                selC = u.select_atoms(C_base)
                d_list.append(selC)
        else:
                selU = u.select_atoms(U_base)
                d_list.append(selU)

        if resname == 'A5 or G5 or C5 or U5':
                selsug5 = u.select_atoms(5_sugar)
                d_list.append(selsug5)
        elif resname != 'A5 or A3 or G5 or G3 or C5 or G3 or U5 or U3':
                selsug = u.select_atoms(sugar)
                d_list.append(selsug)
        else: 
                selsug3 = u.select_atoms(3_sugar)
                d_list.append(selsug3)

        selphos = u.select_atoms(phos)
        d_list.append(selphos)

u_prot = u.select_atoms(pro)
nRes = len(u_prot.residues)
mRNA = len(d_list)            #d_list contains nucleic selections in order above

#Create empty numpy arrays filled with zeroes
avg_matrix = zeros((nRes,mRNA))
std_matrix = zeros((nRes,mRNA))

nSteps = 0
while start <= end:
        u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
        nSteps += len(u.trajectory)

        for ts in u.trajectory:
                if ts.frame%1000 == 0:
                        ffprint('Working on timestep %d of trajectory %d' %(ts.frame, start))
#For every timestep, calculate centers of mass for protein and nucleic residues
                for i in range(nRes):
                        respro = u_prot.residues[i]
                        compro = respro.center_of_mass()
                for i in range(mRNA):
                        comrna = d_list[i].center_of_mass()
#Here, make selections to calculate distances between COM and populate arrays
                for i in range(nRes):
                        for j in range(mRNA):
                                dist, dist2 = euclid_dist(compro,comrna)
                                avg_matrix[i,j] += dist
                                std_matrix[i,j] += dist2

        start +=1
ffprint(nSteps)

avg_matrix /= nSteps
std_matrix /= nSteps
std_matrix = sqrt(std_matrix - square(avg_matrix))

out1 = open('%03d.%03d.avg_dist_mtx.dat' %(start,end))
out2 = open('%03d.%03d.stdv_dist_mtx.dat' %(start,end))
for i in range(nRes):
        for j in range(mRNA):
                out1.write('%10f  ' %(avg_matrix[i,j]))
                out2.write('%10f  ' %(std_matrix[i,j]))
        out1.write('\n')
        out2.write('\n')
out1.close()
out2.close()
