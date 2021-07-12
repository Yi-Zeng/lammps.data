# -*- coding: utf-8 -*-
"""

This class is a collection of data file of LAMMPS 


"""


import io
import os
import shutil
import numpy as np
import math
import re
import copy
import codecs
import sklearn
#from scipy import stats
from sklearn.neighbors import NearestNeighbors
from itertools import islice
import matplotlib.pyplot as plt

from itertools import groupby
from operator import itemgetter

#title = ['atoms', 'bonds', 'angles', 'dihedrals']
title = ['atoms', 'bonds', 'angles', 'dihedrals', 'impropers']
coor  = ['x', 'y', 'z']


class Lammps(object):
    
    def writting_xyz_file(self, atoms_data, data_file_name):
        atom_matrix = []
        for i in range(int(len(atoms_data))):
            xyz = atoms_data[i,-3:]
            if int(atoms_data[i, 2]) == 1:
                atom_matrix.append('C %s %s %s \n'%(xyz[0], xyz[1], xyz[2]))
            else:
                atom_matrix.append('H %s %s %s \n'%(xyz[0], xyz[1], xyz[2]))              
                
        f = open(data_file_name, 'w')
        f.write('%s \n'%(int(len(atoms_data))))
        f.write('Atoms. Timestep: 0 \n')
        for i in range(int(len(atom_matrix))):
            f.write(atom_matrix[i])
        f.close()
    
    

    
        
    def writting_data_file(self, data_file_name, mol_data, extra_distance = 2):
        def write_matrix_to_file(Matrix, Matrix_type, f):      
            if Matrix_type == '0':       
                formatstr = '{:.0f} {:.0f} {:.0f} {:.0f}'         
                for i in range(3):
                    formatstr += ' {:f}'             
            else:
                formatstr = '{:.0f}' 
                for i in range(1,len(Matrix[0,:])):
                    formatstr += ' {:.0f}'
            formatstr += ' \n'          
            for i in range(len(Matrix)):
                f.write(formatstr.format(*Matrix[i]))
                

        coor_lim = np.zeros([2,3])

        coor_lim[0] = np.min(mol_data['atoms'][:,4:7], axis=0) - extra_distance
        coor_lim[1] = np.max(mol_data['atoms'][:,4:7], axis=0) + extra_distance
      
        f = open(data_file_name, 'w')
        f.write('LAMMPS Description\n')
        f.write('\n' + '\n')
        for i in range(len(title)):
            f.write('   ' + str(len(mol_data['%s'%(title[i])][:,0])) + '\t%s\n'%(title[i]))
        f.write('\n' + '\n')
        f.write('          2  atom types\n')
        f.write('          1  bond types\n')
        f.write('          1  angle types\n')
        f.write('          1  dihedral types\n')
        f.write('\n' + '\n')
        for i in range(len(coor)):
            f.write(' ' + str(coor_lim[0, i])+' '+str(coor_lim[1, i])+' %slo %shi\n'%(coor[i], coor[i]))
        f.write('\n')
        f.write('Masses\n')
        f.write('\n')
        f.write('  1 14.027\n')
        f.write('  2 15.035\n')
        f.write('\n')
        for i in range(len(title)):
            f.write('%s\n'%(title[i].title()))
            f.write('\n')
            write_matrix_to_file(mol_data['%s'%(title[i])], '%s'%i, f)
            f.write('\n')
        f.close()

    def writting_data_COMPASS(self, data_file_name, mol_data, extra_distance = 2):
        def write_matrix_to_file(Matrix, Matrix_type, f):      
            if Matrix_type == '0':       
                formatstr = '{:.0f} {:.0f} {:.0f} {:.0f}'         
                for i in range(3):
                    formatstr += ' {:f}'             
            else:
                formatstr = '{:.0f}' 
                for i in range(1,len(Matrix[0,:])):
                    formatstr += ' {:.0f}'
            formatstr += ' \n'          
            for i in range(len(Matrix)):
                f.write(formatstr.format(*Matrix[i]))
                

        coor_lim = np.zeros([2,3])

        coor_lim[0] = np.min(mol_data['atoms'][:,4:7], axis=0) - extra_distance
        coor_lim[1] = np.max(mol_data['atoms'][:,4:7], axis=0) + extra_distance
      
        f = open(data_file_name, 'w')
        f.write('LAMMPS Description\n')
        f.write('\n' + '\n')
        for i in range(len(title)):
            f.write('   ' + str(len(mol_data['%s'%(title[i])][:,0])) + '\t%s\n'%(title[i]))
        f.write('\n' + '\n')
        f.write('          2  atom types\n')
        f.write('          51  bond types\n')
        f.write('          97  angle types\n')
        f.write('          352  dihedral types\n')
        f.write('          10  improper types\n')
        f.write('\n' + '\n')
        for i in range(len(coor)):
            f.write(' ' + str(coor_lim[0, i])+' '+str(coor_lim[1, i])+' %slo %shi\n'%(coor[i], coor[i]))
        f.write('\n')
        f.write('Masses\n')
        f.write('\n')
        f.write('  1 12.01115\n')
        f.write('  2 1.00797\n')

        f.write('\n')
        for i in range(len(title)):
            f.write('%s\n'%(title[i].title()))
            f.write('\n')
            write_matrix_to_file(mol_data['%s'%(title[i])], '%s'%i, f)
            f.write('\n')
        f.close()
        
        
class datafile(Lammps):
    
    def __init__(self, file_name):
        self.name = file_name
    """    
    def readingDATA(self):
        with open(self.name) as f:
            temp = f.readlines()
        
        
        for i, val_i in enumerate(temp):
    """    
        
        
        
    def ideal_crystal(self, n_C, Nx, Ny, Nz, inverse_z = '0', swithch_xy = '0'):
        atoms_pos, n_mol = self.crs_seq1(n_C , Nx, Ny, Nz, inverse_z , swithch_xy )
        mol_data = self.mol_inf(atoms_pos, n_mol)
        return mol_data
    
    def output(self, mol_infor):
        #super().writting_data_file(self.name, mol_infor)
        super().writting_data_COMPASS(self.name, mol_infor)
        super().writting_xyz_file(mol_infor['atoms'], '%s.xyz'%self.name)
       
    
    def mol_inf(self, atoms_pos, n_mol):
        mol_inf = self.mol_data_array(n_mol)
        mol_inf['atoms'][:, 4:7] = atoms_pos
        return mol_inf
    
    def combine_data_file_with_pdb(self, n_C , Nx, Ny, Nz, pdb_file_name):
        """combining two data information from different alkanes to form a new data file   """
        atoms_pos_1, n_mol_1 = self.crs_seq1(n_C , Nx, Ny, Nz)
        atoms_pos_2 = pdbfile(pdb_file_name).data
        n_mol_2 =  pdbfile(pdb_file_name).n_mols       
        n_mol = np.append(n_mol_1, n_mol_2, axis = 0)
        atoms_pos = np.append(atoms_pos_1, atoms_pos_2, axis = 0)
        mol_info_mixture = self.mol_inf(atoms_pos, n_mol)

        return mol_info_mixture
    
    def combine_data_file(self, n_C , Nx, Ny, Nz, inverse_z, swithch_xy):
        """combining two data information from different alkanes to form a new data file   """
        atoms_pos, n_mol = self.crs_seq1(n_C[0] , Nx[0], Ny[0], Nz[0], inverse_z[0], swithch_xy[0])
        for i in range(1, len(n_C)):
            temp = self.crs_seq1(n_C[i] , Nx[i], Ny[i], Nz[i], inverse_z[i], swithch_xy[i])
            n_mol = np.append(n_mol, temp[1], axis = 0)
            atoms_pos = np.append(atoms_pos, temp[0], axis = 0)
        mol_info_mixture = self.mol_inf(atoms_pos, n_mol)
        return mol_info_mixture


    
    def straight_single_mol(self,n_C):
        dist = 2.54005432 # the distance between atom #1 and #3
        init_two_atoms = np.array([[  0.00,   0.00,   0.00],[  dist/2, 1.6, 0]])      
        mol_pos = np.zeros([n_C, 3])
        mol_pos[0:2] =init_two_atoms  
        for i in range(int(n_C//2 - 1)):
            mol_pos[2 * (i+1)] = init_two_atoms[0]
            mol_pos[2 * (i+1) + 1] = init_two_atoms[1]    
            mol_pos[2 * (i+1), 0] = mol_pos[2*i,0] + dist
            mol_pos[2 * (i+1) + 1, 0] = mol_pos[2*i+1,0] + dist
        return mol_pos
    
    def crs_seq1(self, n_C = 20, Nx = 2, Ny = 8, Nz = 20,  inverse_z = 0 , swithch_xy = 0 ):
        a = 3.9*1.122
        dy = a #distance between two molecules
        dz = math.sqrt(3)*a 
        n_mols_in_a_bulk = Ny*Nz
        n_atoms_in_a_bulk = n_mols_in_a_bulk * n_C        
        int_width = 4#*1.122
        n_mol = np.array([[n_C, n_mols_in_a_bulk*Nx]])      
        R_Vector = self.straight_single_mol(n_C)
        atoms_pos = np.zeros([n_mol[0,1]*n_C, 3])
        
        for i in range(0, n_mols_in_a_bulk):
            for j in range(0, n_C):
                atoms_pos[i * n_C + j, 0] = R_Vector[j, 0]

        mol_num = 0
        for iz in range(Nz):
            if iz % 2 == 0:
                for iy in range(Ny):
                    for j in range(n_C):
                        atoms_pos[mol_num*n_C+j, 1] = R_Vector[j, 1] + dy * iy
                        atoms_pos[mol_num*n_C+j, 2] = R_Vector[j, 2] +  0.5 * dz * iz
                    mol_num += 1
            else:
                for iy in range(Ny):
                    for j in range(n_C):
                        atoms_pos[mol_num*n_C+j, 1] = R_Vector[j, 1] + dy * (iy + 0.5) 
                        atoms_pos[mol_num*n_C+j, 2] = R_Vector[j, 2] +  0.5 * dz * iz
                    mol_num += 1
                    
                                        
        for n in range(1,Nx):
            for i in range(n_atoms_in_a_bulk):
                atoms_pos[n * n_atoms_in_a_bulk + i, 0] = atoms_pos[i, 0] + n * (R_Vector[-1,0] - R_Vector[0,0] + int_width)
                atoms_pos[n * n_atoms_in_a_bulk + i, 1] = atoms_pos[i, 1]   
                atoms_pos[n * n_atoms_in_a_bulk + i, 2] = atoms_pos[i, 2] 
        if  inverse_z == 1:  
            atoms_pos[:, 2] =  -atoms_pos[:, 2] - a

        
        if (swithch_xy == 1):
            atoms_pos[:,[0, 1]] =  atoms_pos[:,[1, 0]]

        return atoms_pos, n_mol
                    

        
        
                     
    def crs_seq2(self, n_C = 20, Nx = 2, Ny = 8, Nz = 20, a = 4, swithch_xz = '0'):
        """
        Ny = 6
        Nz = 4
        
        """
        dy = a #distance between two molecules
        dz = math.sqrt(3)*a 
        n_mols_in_a_bulk = Ny*Nz
        n_atoms_in_a_bulk = n_mols_in_a_bulk * n_C        
        int_width = 3.93*1.122
        n_mol = np.array([[n_C, n_mols_in_a_bulk*Nx* n_C]])      
        R_Vector = self.straight_single_mol(n_C)
        atom_pos = np.zeros([n_mol[0,1], 3])
        ####Atoms##
        for i in range(0, n_mols_in_a_bulk):
            for j in range(0, n_C):
                atom_pos[i * n_C + j, 0] = R_Vector[j, 0]
                                 
        Molecular_Laber = 0
        for n in range(Ny):
            for m in range(Nz):
                if (n == 0):
                    Np1 = 1
                else:
                    Np1 = 2
                if (m == 0):
                    Np2 = 1
                else:
                    Np2 = 2
                for p1 in range(Np1):
                    for p2 in range(Np2):
                        for s in range(2):
                            for j in range(n_C):
                                atom_pos[Molecular_Laber*n_C+j, 5] = R_Vector[j, 1] + ((-1)**p1) * n * dy + s * 0.5 * dy
                                atom_pos[Molecular_Laber*n_C+j, 6] = R_Vector[j, 2] + ((-1)**p2) * m * dz + s * 0.5 * dz
                            Molecular_Laber += 1
                            
        for n in range(1,n_bulks):
            for i in range(n_atoms_in_a_bulk):
                atom_pos[n * n_atoms_in_a_bulk + i, 4] = atom_pos[i, 4] + n * (R_Vector[-1,0] - R_Vector[0,0] + int_width)
                atom_pos[n * n_atoms_in_a_bulk + i, 5] = atom_pos[i, 5]   
                atom_pos[n * n_atoms_in_a_bulk + i, 6] = atom_pos[i, 6] 
    
        if (swithch_xz == '1'):
            atom_pos[:,[4,6]] =  atom_pos[:,[6,4]]
        else:
            pass
        return atom_pos
    
    def mol_data_fromMOLTEMPLATE(self):
        """
        This function reads the array of atoms, bonds, angles and dihedrals from MOLTEMPLATE datafile
        
        
        """      
        def readingMOLFILES(name):
            with open(name) as f:
                data = f.readlines()
            return np.loadtxt(data)
        files_list = ['Data Atoms', 'Data Bonds', 'Data Angles', 'Data Dihedrals', 'Data Impropers']
        #title_list = ['atoms', 'bonds', 'angles', 'dihedrals', 'impropers']
        mol_data = {}

        for i, val_i in enumerate(files_list):
            mol_data['%s'%(title[i])] = readingMOLFILES('%s'%(val_i))
        
        

        temp_arr = mol_data['atoms'][:,2]
        for i, val_i in enumerate(temp_arr):
            if val_i == float(6):
                mol_data['atoms'][i, 2] = 1
            elif val_i == float(13):
                mol_data['atoms'][i, 2] = 2       

        return mol_data
        
    def mol_data_array(self, n_mol):
        """
        This function creates the array of atoms, bonds, angles and dihedrals with the inforamtion of n_mol
        only for 
        
        
        """
        #n_mol = self.n_mols

        n_col = [7, 4, 5, 6] #the number of columns in 'atoms. bonds, angles and dihedrals'
        n_data = {}
        for i in range(len(title)):
            n_data['number of %s'%(title[i])] =  (np.zeros([len(n_mol)]).astype('int'),)
            for j in range(int(len(n_mol))):
                n_data['number of %s'%(title[i])][0][j] = n_mol[j, 1] * (n_mol[j, 0] - i) 
            n_data['number of %s'%(title[i])] += (n_col[i],) 
                
        mol_data ={}
        for i in range(len(title)):
            sum_var = sum(n_data['number of %s'%(title[i])][0])
            mol_data['%s'%(title[i])] = np.zeros([sum_var, n_data['number of %s'%(title[i])][1]])
            mol_data['%s'%(title[i])][:,0] = range(1, sum_var + 1)
        


        ############################ atoms part calcualtion
        mol_data['%s'%(title[0])][:, 2] = 1 # atoms type
        mol_data['%s'%(title[0])][:, 3] = 0
        
        cum_var = np.insert(np.cumsum(n_data['number of %s'%(title[0])][0]), 0, 0) # the cumsum value of atoms and 0 at the beginning
        cum_mol = np.insert(np.cumsum(n_mol[:, 1]), 0, 0)

            
        for i in range(len(n_mol)):
            for j in range(n_mol[i, 1]):
                 mol_data['%s'%(title[0])]\
                 [j*n_mol[i, 0]+cum_var[i]:n_mol[i, 0]*(j+1)+cum_var[i], 1] = j + 1 + cum_mol[i]
            
            mol_data['%s'%(title[0])][cum_var[i] : cum_var[i + 1] : n_mol[i, 0], 2] = 2
            mol_data['%s'%(title[0])][(cum_var[i] + n_mol[i, 0] - 1) : cum_var[i + 1]: n_mol[i, 0], 2] = 2
        #################################################################
        
        ############################ bonds, angles and dihedrals  calcualtion
        for k  in range(1, len(title)):
            mol_data['%s'%(title[k])][:, 1] = 1 # type of bonds, angles and dihedrals
            cum_vars = np.insert(np.cumsum(n_data['number of %s'%(title[k])][0]), 0, 0) # the cumsum values of bonds, angles and dishedrals and 0 at the beginning
            for i in range(len(n_mol)):
                for j in range(n_mol[i, 1]):
                    mol_data['%s'%(title[k])][j*(n_mol[i, 0]-k)+cum_vars[i]:(j+1)*(n_mol[i, 0]-k)+cum_vars[i], 2] \
                    = np.arange(n_mol[i, 0]*j + 1 + cum_var[i], n_mol[i, 0] * j + n_mol[i, 0] + cum_var[i] - k + 1)
            for i in range(3, len(mol_data['%s'%(title[k])][0])): # the details of bonds etc
                mol_data['%s'%(title[k])][:, i] = mol_data['%s'%(title[k])][:, i - 1] + 1 
        return mol_data
     
    
class xyzfile(Lammps):
     def __init__(self, file_name):
        self.name = file_name
        #self.type = data_type
        self.data = self.collect_pos_from_XYZ[0] 
        self.n_mols = self.collect_pos_from_XYZ[1] 
        self.atoms = self.collect_pos_from_XYZ[2]
        self.boun = self.xyz_boundary
        self.length = self.xyz_len
        self.timestep = self.collect_pos_from_XYZ[3]
        self.mols_char = self.mol_chara
        #self.pbc_data = self.pbc_modification()
        #self.select = self.picked_molecules_presentaion()
        
    @property
    def collect_pos_from_XYZ(self):    
        """
        Collect the atoms element and position from the xyz file
        return  all the atoms position in an array and the moleuclar type and their number
        """      
        with open(self.name) as f:
            lines = f.readlines()
        index_list = []
        for i in range(len(lines)):    
            if 'Atoms' in lines[i]:
                index_list.append(i)
        n_atoms = int(lines[0])
        n_fram = len(lines)//(n_atoms + 2)
        text = []
        atoms_pos = np.zeros([n_fram, n_atoms, 3])
        atoms_type = np.zeros([n_atoms]).astype('str')
        temp_i = 0
        if n_fram == 1:
            for j,s in enumerate(lines):
                temp_line = np.array(re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", s)).astype('float')
                if len(temp_line) == 3:
                    atoms_pos[0, temp_i] = temp_line
                    temp_i += 1
        for i in range(n_fram - 1):
            text.append(lines[index_list[i]])
            for j,s in enumerate(lines[(index_list[i] + 1) : (index_list[i + 1] - 1)]):
                temp_line = np.array(re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", s)).astype('float')
                if len(temp_line) == 4:
                    atoms_pos[-1, j] = temp_line[1:]
        for j,s in enumerate(lines[(index_list[-1] + 1) :]):
            text.append(lines[index_list[-1]])            
            temp_line = np.array(re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", s)).astype('float')
            atoms_type[j] = s[0]
            if len(temp_line) == 4:
                atoms_pos[-1, j] = temp_line[1:]
                atoms_type[j] = temp_line[0]
        atoms_num = np.array([]) 
        
        if atoms_type[0] == 'C':        
            for i in range(len(atoms_type)):            
                atoms_num = np.append(atoms_num, i) 
            len_shape = int(len(atoms_num)/2)
            label_arr = np.reshape(atoms_num, (len_shape,2)).astype(int)
            atoms_elem_num = np.unique((np.diff(np.reshape(atoms_num, (len_shape,2))) + 1),  True, False, True)
            label_list = (np.diff(np.reshape(atoms_num, (len_shape,2))) + 1).reshape([len_shape]).astype(int)
            seq = np.argsort(atoms_elem_num[1])
            n_mols = np.array([atoms_elem_num[0], atoms_elem_num[2], ])[:, seq].transpose().astype('int')
            atoms_array  = {}
            for i in range(len(n_mols)):
                atoms_array[i] = np.zeros([n_fram, len(np.where(label_list == n_mols[i,0])[0]), n_mols[i,0], 3])
                for j, val_j in enumerate(label_arr[np.where(label_list == n_mols[i, 0])[0]]):
                    atoms_array[i][:, j] = atoms_pos[:, val_j[0] : val_j[1] + 1]
        elif atoms_type[0] == '6.0':
            atoms_num_C = np.where(atoms_type == '6.0')[0]
            H_list = np.where(atoms_type == '13.0')[0]
            if len(H_list) > 500:
                H_index = np.array([])
                for i, val_i in enumerate(H_list[:200]):
                    if (H_list[i+1] == val_i + 1 and H_list[i+2] == val_i + 2) :
                        H_index = np.append(H_index, val_i) 
            atom_num_in_mol = int(H_index[1] + 3) #number of atoms in one molecule
            mol_num = len(lines)//atom_num_in_mol #number of molecules
            len_shape = len(atoms_num)//mol_num #number of carbon in one molecule
            label_arr = np.reshape(atoms_num, (mol_num, len_shape)).astype(int)
            atoms_elem_num = np.unique((2 + 1),  True, False, True)
            label_list = (np.diff(np.reshape(atoms_num, (len_shape,2))) + 1).reshape([len_shape]).astype(int)
            seq = np.argsort(atoms_elem_num[1])
            n_mols = np.array([atoms_elem_num[0], atoms_elem_num[2], ])[:, seq].transpose().astype('int')
            atoms_array  = {}
            for i in range(len(n_mols)):
                atoms_array[i] = np.zeros([n_fram, len(np.where(label_list == n_mols[i,0])[0]), n_mols[i,0], 3])
                for j, val_j in enumerate(label_arr[np.where(label_list == n_mols[i, 0])[0]]):
                    atoms_array[i][:, j] = atoms_pos[:, val_j[0] : val_j[1] + 1]
        return atoms_pos, n_mols, atoms_array, text
        

    def pbc_modification(self, data, timestep = 0):
        """
        Modify the coodinates of the this molecule with the periodic boudnary condiction
        data: the position information, (x,y,z)
        boun_len: the size of the simulation box of all the molecules.    
        """
        len_chara = np.abs(self.length[timestep] - 5)
        data_pbc = copy.deepcopy(data)
        for i in range(1, len(data)):
            line_arr = data_pbc[i] - data_pbc[i - 1]
            for j in range(3):
                if line_arr[j] > len_chara[j]:
                    data_pbc[i, j] -= self.length[timestep, j]
                elif line_arr[j] < - len_chara[j]:
                    data_pbc[i, j] += self.length[timestep, j]
        return data_pbc   
    
    @property
    def xyz_boundary(self):
        bound = np.zeros((2,) + self.data.min(1).shape)
        bound[0] = self.data.min(1)
        #self.data.min(1).reshape((1,) + self.data.min(1).shape).shape
        bound[1] = self.data.max(1)
        return bound
    @property
    def xyz_len(self):
        return self.data.max(1) - self.data.min(1)
    @property
    def distribution(self, num_bins = 20, order = 2): #order 0, 1, 2 mean direction x, y, z
        dist_data = ()
        for i in range(len(self.n_mols)):
            dist_data += (self.atoms[i].mean(-2),)
        dist_hist = np.zeros([len(self.atoms[0]), len(self.n_mols), num_bins, 2])
        
        for i in range(len(self.n_mols)):
            for j in range(len(dist_hist)):
               hist, bins = np.histogram(dist_data[i][j, :,order], num_bins)
               dist_hist[j, i, :, 0] = np.linspace(bins[0], bins[-1], len(hist))
               dist_hist[j, i, :, 1] = hist
        return dist_hist 
        


    def neighbor_mols(self, target_mol, block = 'off'):
            
        def distance_neigbor(i, target, mirro_arr):
            if i == 0:
                nei_ind = np.where(np.linalg.norm((target[0] - mirro_arr[:,:,-1,:]),axis = -1) < cut_d)
            elif i == 1:
                nei_ind = np.where(np.linalg.norm((target[-1] - mirro_arr[:,:,0,:]),axis = -1) < cut_d)
            nei_arr = np.zeros([len(nei_ind[0]), len(nei_ind)])
            for i in range(len(nei_ind)):
                nei_arr[:, i] = nei_ind[i]

            return nei_arr.astype('int')
            
            
            
        def pdb_atoms(data):   
            data_all = np.zeros((9,) + (data.shape))
            pbc_data = np.zeros_like(data)
            for i in range(len(data)):
                pbc_data[i] = self.pbc_modification(data[i])
            data_all[:] = pbc_data
            data_all[[1, 4, 6], :, :, 1] = pbc_data[:,:,1] - self.length[0, 1]
            data_all[[3, 5, 8], :, :, 1] = pbc_data[:,:,1] + self.length[0, 1]
            data_all[6:9, :, :, 2] = pbc_data[:,:,2] - self.length[0, 2]
            data_all[1:4, :, :, 2] = pbc_data[:,:,2] + self.length[0, 2]
            return data_all

        def atoms_pos_to_xyz_str(atom_arr):
            xyz_str = []

            for i in range(len(atom_arr)):
                if i == 0:
                    xyz_str += ['C %s %s %s \n'%(atom_arr[i, 0],atom_arr[i, 1],atom_arr[i, 2])]
                elif i == len(atom_arr) - 1 :
                    xyz_str += ['C %s %s %s \n'%(atom_arr[i, 0],atom_arr[i, 1],atom_arr[i, 2])]
                else:
                    xyz_str += ['H %s %s %s \n'%(atom_arr[i, 0],atom_arr[i, 1],atom_arr[i, 2])]

            return xyz_str
         
        cut_d = 10
        np_block = len(self.atoms[0][0])//154
        ind_arr = np.zeros([np_block, 2]).astype('int')
        for i in range(np_block):
            ind_arr[i, 0] = i * 154
            ind_arr[i, 1] = (i + 1) *154
        index = target_mol//154
        atoms_1 = self.atoms[0][0][ind_arr[index - 1, 0 ] : ind_arr[index - 1, 1 ] ]
        atoms_2 = self.atoms[0][0][ind_arr[index + 1, 0 ] : ind_arr[index + 1, 1 ] ]
        mirro_arr = np.zeros((2,9) + (atoms_1.shape))
        mirro_arr[0] = pdb_atoms(atoms_1)
        mirro_arr[1] = pdb_atoms(atoms_2)
        target_mol_pos = self.atoms[0][0][target_mol]
        mol_arr_index = ()
        for i in range(2):
            mol_arr_index += (distance_neigbor(i, target_mol_pos, mirro_arr[i]),)
        nei_mol_size = len(mol_arr_index[0]) + len(mol_arr_index[1])
        neighbor_mol_pos = np.zeros([nei_mol_size, 20, 3])
        numb = 0
        for i in range(2):
            for j, val_j in enumerate(mol_arr_index[i]):
                neighbor_mol_pos[numb] = mirro_arr[i, val_j[0], val_j[1]]
                numb += 1
        result = np.zeros((len(neighbor_mol_pos) + 1,) + (target_mol_pos.shape))
        result[:-1] =  neighbor_mol_pos
        result[-1] = target_mol_pos
        
        n_atoms = 20 * len(result)
        f1 = open('neighbor_mols.xyz','w')
        f1.write(str(int(n_atoms)) + '\n') #Write the number of atoms in first line
        f1.write('Atoms. Timestep: \n')
        for i in range(len(result)):
            for j in range(len(result[i])):
                f1.write(atoms_pos_to_xyz_str(result[i])[j])

        f1.close()
        return 0 
  



            
    def select(self, target_tuple, new_name, p_type = 'mol',fram = 0):   
        """Select certain molecules or atoms in xyz file to produce a new xyz file
        target_tuple: a tuple, each element has the label of molecules in each kind of n-alkanes.
        """
        
          
            
        def atoms_pos_to_xyz_str(atom_arr, target = 'on'):
            xyz_str = []
            for i in range(len(atom_arr)):
                if i == 0:
                    xyz_str += ['C %s %s %s \n'%(atom_arr[i, 0],atom_arr[i, 1],atom_arr[i, 2])]
                elif i == len(atom_arr) - 1 :
                    xyz_str += ['C %s %s %s \n'%(atom_arr[i, 0],atom_arr[i, 1],atom_arr[i, 2])]
                else:
                    xyz_str += ['H %s %s %s \n'%(atom_arr[i, 0],atom_arr[i, 1],atom_arr[i, 2])]
            return xyz_str
         
        if p_type == 'mol':
            n_atoms = 0
            for i in range(len(target_tuple)):
                n_atoms += self.n_mols[i, 0] * len(target_tuple[i])

            f1 = open('%s.xyz'%new_name,'w')
            f1.write(str(int(n_atoms)) + '\n') #Write the number of atoms in first line
            f1.write('Atoms. Timestep: \n')
            for i in range(len(target_tuple)):
                for j in target_tuple[i]:
                    for k in range(len(self.atoms[i][fram][j])):
                        f1.write(atoms_pos_to_xyz_str(self.atoms[i][fram][j])[k])
            f1.close()
            
        elif p_type == 'atom':
            target_list  = target_tuple[0]
            a = str(len(target_list))
            f = open('%s'%self.name)
            xyz = f.readlines()[2:]
            f.close()
            f1 = open('%s.xyz'%new_name,'w')
            f1.write(a + '\n') #Write the number of atoms in first line
            f1.write('Atoms. Timestep: \n')
            for k in range(len(target_list)):
                i = int(target_list[k]) - 1
                #for j in range(2 + (i-1) *n_atoms, 2+ i *n_atoms):
                f1.write(xyz[i])
                
            f1.close()
            
            
    def recolor_file(self):
        """change the color of the molecules
        """
        
        
        def atoms_pos_to_xyz_str_color(atom_arr, index):
            xyz_str = []
            color_list = ['C', 'O', 'H']

            xyz_str += ['%s %s %s %s \n'%(color_list[index], atom_arr[0],atom_arr[1],atom_arr[2])]
            return xyz_str
         
        n_atoms = 0
        for i in range(len(self.n_mols)):
            n_atoms += self.n_mols[i, 0] * self.n_mols[i, 1]

        f1 = open(self.name)
        lines = []
        for i in range(len(self.atoms[0])):
            lines.append(str(int(n_atoms)) + '\n') #Write the number of atoms in first line
            lines.append(self.timestep[i])
            for j in range(len(self.atoms)):
                for k in range(len(self.atoms[j][0])):
                    for l in range(len(self.atoms[j][0][k])):
                        lines.append(atoms_pos_to_xyz_str_color(self.atoms[j][i][k][l], j)[0])
        f1.close()
        f = open('%s_traj.xyz'%self.name, 'w')
        for l in lines:
            f.write('%s' % l)
        f.close()
        
    def readlist(self):
        with open(self.name) as f:
            lines = f.readlines()
        return lines
        
    def recolor_list(self):
        """change the color of the molecules
        """
        
        
        def atoms_pos_to_xyz_str_color(atom_arr, index):
            xyz_str = []
            color_list = ['C', 'O', 'H']

            xyz_str += ['%s %s %s %s \n'%(color_list[index], atom_arr[0],atom_arr[1],atom_arr[2])]
            return xyz_str
         
        n_atoms = 0
        for i in range(len(self.n_mols)):
            n_atoms += self.n_mols[i, 0] * self.n_mols[i, 1]

        f1 = open(self.name)
        lines = []
        for i in range(len(self.atoms[0])):
            lines.append(str(int(n_atoms)) + '\n') #Write the number of atoms in first line
            lines.append(self.timestep[i])
            for j in range(len(self.atoms)):
                for k in range(len(self.atoms[j][0])):
                    for l in range(len(self.atoms[j][0][k])):
                        lines.append(atoms_pos_to_xyz_str_color(self.atoms[j][i][k][l], j)[0])
        f1.close()
        return lines
    
    def inangle(v1, v2):
        """ Returns the angle in radians between vectors 'v1' and 'v2'::
    
                >>> angle_between((1, 0, 0), (0, 1, 0))
                1.5707963267948966
                >>> angle_between((1, 0, 0), (1, 0, 0))
                0.0
                >>> angle_between((1, 0, 0), (-1, 0, 0))
                3.141592653589793
        """
        v1_u = v1/ np.linalg.norm(v1)
        v2_u = v2/ np.linalg.norm(v2)
        return np.dot(v1_u, v2_u)           

    #@property    
    def mol_chara(self):
        '''Calculate the charateristics of all molecules in one XYZ file,distance between end atoms and the angle
        Cn: the array of number of carbon in molecules
        n_mol: the number of molecules for different n-alkanes
        '''
            
         
        def angle_between(v1, v2):
            inner = inangle(v1, v2)
            if len(v1) == 2:
                det = np.cross(v1, v2)
            else:
                #det = np.arcsin(np.linalg.norm(np.cross(v1, v2))/np.linalg.norm(v1)/np.linalg.norm(v2))/math.pi*180
                det = -1
            if det < 0:
                return inner
            else:
                return -inner

        def mol_chara_cal(data, timestep):
            """
            Calculte the order and length of a single  molecular chain.
            data is an array of  the position of all the atoms in this molcule
            
            
            return: an array with size [1,7] including x, y, z, distance between end atoms, alpha, beta, gamma
            """
            data_pbc = self.pbc_modification(data, timestep)
            #vec_arr = np.zeros([2, int((data.shape[0]-2)/2), 6])
            #vec_arr[0, :, :3] = data_pbc[2:][::2] - data_pbc[:-2][::2]
            #vec_arr[1, :, :3] = data_pbc[3:][::2] - data_pbc[:-3][::2]
            #head_end_arr = data[1:] - data[:-1]
            head_end_len = data_pbc[-2] - data_pbc[0]
            
            #ref_vec = [1, 0, 0]#the reference vector
            #for i in range(int(len(vec_arr))):
                #vec_arr[i, :, 3:] = mol_ang_projection(vec_arr[i,: , :3])
                #vec_arr[i, :, 3:] = mol_ang_2d_proj(vec_arr[i,: , :3])      
            x = np.mean(data_pbc[:,0])
            y = np.mean(data_pbc[:,1])
            z = np.mean(data_pbc[:,2])
            length = np.linalg.norm(head_end_len)/(len(data)-2)
            end_to_end_vec = head_end_len/np.linalg.norm(head_end_len)
            #end_to_end_vec = (data_pbc[-2] - data_pbc[0])/np.linalg.norm(data_pbc[-2] - data_pbc[0])
            #mol_angle = np.zeros([len(data), 5])
            #mol_angle[:,:3] = data
            #mol_angle[:, 3] = length
            #mol_angle[:, 4] = vec_arr[:, :, 3].mean()
            #mol_angle = np.append(np.array([x, y, z, length]), vec_arr[:, :, 3:].std(1).mean(0))
            mol_angle = np.append(np.array([x, y, z, length]), end_to_end_vec)
            return mol_angle
        
        
        def mol_ang_projection(data):
        
            ang_arr = np.zeros([len(data), 3])
            for i in range(len(data)):
                ang_arr[i, 0] = inangle([1, 0, 0], data[i])
                ang_arr[i, 1] = inangle([0, 1, 0], data[i])
                ang_arr[i, 2] = inangle([0, 0, 1], data[i])
            return ang_arr
        
        def mol_ang_2d_proj(data):
            ang_arr = np.zeros([len(data), 3])
            for i in range(len(data)):
                ang_arr[i, 0] = angle_between([1, 0], data[i][0:2])
                ang_arr[i, 1] = angle_between([1, 0], data[i][1:3])
                ang_arr[i, 2] = angle_between([1, 0], data[i][[2,0]])
            return ang_arr
        

        atoms_arr = self.atoms

        mols_chara = np.zeros([ len(self.atoms[0]), np.sum(self.n_mols[:, 1]), 7]) # 7 elements:x, y, z, length, alpha, beta, gamma
        cum_mol = np.insert(np.cumsum(self.n_mols[:, 1]), 0, 0)
        for i in range(len(self.n_mols)):       # molecular type                
            for j in range(len(self.atoms[0])):  #  time fram
                for k in range(self.n_mols[i, 1]):  # single molecule
                    mols_chara[j, cum_mol[i] + k ] = mol_chara_cal(atoms_arr[i][j][k], j)                 
        return mols_chara
    

    
    def sol_ratio(self, file = 0):
        """change the color of the molecules
        """
        
        
        def atoms_pos_to_xyz_str_color(atom_arr, index):
            xyz_str = []
            color_list = ['C', 'O', 'H']

            xyz_str += ['%s %s %s %s \n'%(color_list[index], atom_arr[0],atom_arr[1],atom_arr[2])]
            return xyz_str
        
        def nearest_neighbor():
            return 0



        chara_len = 1.15
        mol_sol_arr = self.mols_char
        n_mols = self.n_mols
        cum_mol = np.insert(np.cumsum(n_mols[:, 1]), 0, 0)
        sol_ratio = np.zeros([len(mol_sol_arr), len(n_mols), 4])
        for i in range(len(mol_sol_arr)):
            for j, val_j in enumerate(cum_mol[:-1]):
                #sol_ratio[i, j] = len(np.where(mols_chara[j][i][:, 4:].mean(-1) < 12)[0])/self.n_mols[j,1]
                #sol_ratio[i, j, 0] = len(np.where(mol_sol_arr[i, val_j:cum_mol[j + 1], 3] > chara_len)[0])/self.n_mols[j,1]
                sol_ratio[i, j, 1] = np.abs((mol_sol_arr[i, val_j:cum_mol[j + 1], 4:]*np.array([1,0,0])))[:, 0].mean()
                sol_ratio[i, j, 2] = np.abs((mol_sol_arr[i, val_j:cum_mol[j + 1], 4:]*np.array([0,1,0])))[:, 1].mean()
                sol_ratio[i, j, 3] = np.abs((mol_sol_arr[i, val_j:cum_mol[j + 1], 4:]*np.array([0,0,1])))[:, 2].mean()
        
        if file == 0:
            
            for i in range(len(mol_sol_arr)): #time fram
                straight_mols ={}
                tot_stra_mols = np.array([])
                for j, val_j in enumerate(cum_mol[:-1]):
                    straight_mols[j] = np.where(mol_sol_arr[i, val_j:cum_mol[j + 1], 3] > chara_len)[0] + val_j
                    tot_stra_mols = np.append(tot_stra_mols, straight_mols[j]).astype(int) 
                #neighbor_mols = {}
                for j in range(len(straight_mols)):
                    neighbor = np.array([])

                    for k in straight_mols[j]:
                        #cos_theta = np.sum(mol_sol_arr[i, k, 4:]* mol_sol_arr[i, tot_stra_mols, 4:], axis = -1)
                        dist = np.linalg.norm(mol_sol_arr[i, k, :3] - mol_sol_arr[i, tot_stra_mols, :3], axis = -1)
                        #neighbor_cos = np.where(abs(cos_theta) > 0.985)[0]
                        neighbor_dist = np.where(dist < 10)[0] 
                        #neighbor_k = np.intersect1d(neighbor_cos, neighbor_dist)
                        #neighbor_k = neighbor_dist
                        if len(neighbor_dist) > 3:
                            neighbor = np.append(neighbor, k)
                    neighbor = np.unique(neighbor)
                    #sol_mols = np.intersect1d(straight_mols, neighbor)
    
                    #sol_num = np.zeros([len(n_mols)])
                    """
                    for k in sol_mols:
                        for l in range(len(cum_mol) - 1):
                            if ((k < cum_mol[l + 1]) & (k >= cum_mol[l])):
                                sol_num[k] += 1
                    """
                    #neighbor_mols[j] = neighbor
                    sol_ratio[i, j, 0] = len(neighbor)/n_mols[j, 1]

                        
 
        elif file == 1:
            lines = []
            for i in range(len(mol_sol_arr)): #time fram
                straight_mols ={}
                tot_stra_mols = np.array([])
                for j, val_j in enumerate(cum_mol[:-1]):
                    straight_mols[j] = np.where(mol_sol_arr[i, val_j:cum_mol[j + 1], 3] > chara_len)[0] + val_j
                    tot_stra_mols = np.append(tot_stra_mols, straight_mols[j]).astype(int) 
                neighbor_mols = {}
                for j in range(len(straight_mols)):
                    neighbor = np.array([])

                    for k in straight_mols[j]:
                        #cos_theta = np.sum(mol_sol_arr[i, k, 4:]* mol_sol_arr[i, tot_stra_mols, 4:], axis = -1)
                        dist = np.linalg.norm(mol_sol_arr[i, k, :3] - mol_sol_arr[i, tot_stra_mols, :3], axis = -1)
                        #neighbor_cos = np.where(abs(cos_theta) > 0.985)[0]
                        neighbor_dist = np.where(dist < 10)[0] 
                        #neighbor_k = np.intersect1d(neighbor_cos, neighbor_dist)
                        #neighbor_k = neighbor_dist
                        if len(neighbor_dist) > 3:
                            neighbor = np.append(neighbor, k)
                    neighbor = np.unique(neighbor).astype(int)
                    neighbor_mols[j] = neighbor
                    sol_ratio[i, j, 0] = len(neighbor)/n_mols[j, 1]
                sol_num = np.zeros([len(neighbor_mols)])
                for j in range(len(sol_num)):
                    sol_num[j] = len(neighbor_mols[j])
                    

                if int(np.sum(sol_num * n_mols[:, 0])) != 0:
    
                    lines.append(str(int(np.sum(sol_num * n_mols[:, 0]))) + '\n') #Write the number of atoms in first line    
                    lines.append(self.timestep[i])
                    for j in range(len(neighbor_mols)):                        
                        for k in neighbor_mols[j]:
                            for l in range(n_mols[j, 0]):
                                lines.append(atoms_pos_to_xyz_str_color(self.atoms[j][i, k - n_mols[j, 1], l], j)[0])
    
            f = open('%s_sol_traj.xyz'%(self.name), 'w')
            for l in lines:
                f.write('%s' % (l))
            f.close()  
        return sol_ratio
        

    
    def to_dataFile(self, i):
        """Transfer atoms position to data file
        n_mol: number of molecules
        Cn: number of atoms in one molecule
        file_name: PDB file name
        This function is only applied for homogenous alkanes
        """
        data_filled = datafile('C20_C30').mol_inf(self.data[i], self.n_mols)
        super().writting_data_file('C20_C30.data', data_filled)
        

        
class pdbfile(Lammps):
    def __init__(self, file_name):
        self.name = file_name
        self.data = self.collect_pos_from_packmol_pdb()[0]
        self.n_mols = self.collect_pos_from_packmol_pdb()[1]
        
    def collect_pos_from_packmol_pdb(self):
        f = open(self.name)   
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        lines = f.readlines()
        f.close()
        n_mols = np.zeros([1,2])
        temp_array = np.zeros([len(lines[:-1]), len(re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", lines[0]))])
        for i,s in enumerate(lines[:-1]):
            temp_array[i] = np.array(re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", s)).astype('float') 
        atoms_pos = temp_array[:,-3:] 
        j = 0
        for i in range(1, len(temp_array)):
            if temp_array[i, 2] < temp_array[i -1, 2]:
                n_mols[j, 0] = int(temp_array[i-1, 1])
                n_mols[j, 1] = int(temp_array[i-1, 2])
                n_mols = np.append(n_mols, np.array([[0, 0]]), axis = 0)
                j += 1
        n_mols[j, 0] = int(temp_array[-1, 1])
        n_mols[j, 1] = int(temp_array[-1, 2])
        n_mols = n_mols.astype('int')             
        return atoms_pos, n_mols
    
    def collect_POS_CON_from_pdb_all(self):
        with open(self.name) as f:   
            lines = f.readlines()[3:]
        atoms = np.zeros([62, 6])
        atoms[:, 0] = np.arange(1,63)
        atoms[:, 1] = 1
        atoms[:20, 2] = 1
        atoms[20:, 2] = 2           
        n_atom = 0    
        for i in lines:
            a = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", i)
            if len(a) == 9:
                atoms[n_atom, 3:] = np.array(a).astype('float')[3:6]
                n_atom += 1
                
        
        return atoms
    
    
    def to_dataFile(self):
        data_filled = datafile('C20_C30').mol_inf(self.data, self.n_mols)
        super().writting_data_file('C20_C30.data', data_filled)
        
    def to_xyzFile(self):
        """
        Unfinished, need add the atoms based on differnet length of the molecules
        
        
        """
        atom_matrix = []
        atoms_data = self.data
        for i in range(int(len(atoms_data))):
            xyz = atoms_data[i,-3:]
            if atoms_data[i, 2] == 2:
                atom_matrix.append('C %s %s %s \n'%(xyz[0], xyz[1], xyz[2]))
            else:
                atom_matrix.append('H %s %s %s \n'%(xyz[0], xyz[1], xyz[2]))              
                
        f = open('%s.xyz'%self.name, 'w')
        f.write('%s \n'%(int(len(atoms_data))))
        f.write('Atoms. Timestep: 0 \n')
        for i in range(int(len(atom_matrix))):
            f.write(atom_matrix[i])
        f.close()
    
    
    
class logfile(Lammps):

    def __init__(self, file_name):
        self.name = file_name
        self.caption = self.read_from_log()[0]
        self.data = self.read_from_log()[1]
        self.ave_data = self.ave_log()
                                 
    
    
    def read_from_log(self):
        f = open(self.name)
        start_index  = f.read().index('Step ')
        f.seek(0,0)
        f.seek(start_index, 0)        
        lines = f.readlines()    
        f.close()
        for i in range(len(lines)):
            if ('Loop time of ') in lines[i]:
                end_index = i
        cap = re.findall('\w+', lines[0])
        data_log = np.loadtxt(lines[1:end_index])
        caption = {}
        for i in range(len(cap)):
            caption['%s'%(cap[i])] = i
        return caption, data_log
      
    def ave_log(self, step = 10000):
        """Get the mean value for 1 ns in LAMMPS cooling process.
        data: the raw log data from LAMMPS
        step: the number of steps to average the 1 ns time period
        
        """
        data = self.data
        ave_data = np.zeros([int(data.shape[0]/step),data.shape[1]])
        for i in range(len(ave_data)):
            ave_data[i] = np.mean(data[i*step: (i + 1)*step], axis=0)
        return ave_data

        
        
class fixave(Lammps):
    
    def __init__(self, file_name):
        self.name = file_name
        self.data = self.read_data()
    
    def read_data(self):
        with open(self.name) as f:
            lines = f.readlines()[2:]
        data = np.loadtxt(lines)
        return data
            
class space_ave(Lammps):
    """ This is a class for lammps data file got by command " fix ave/spatial"
    
    """
    

    def __init__(self, file_name):
        self.name = file_name
        self.data = self.read_data()
        self.ave_data = self.ave_data() # The data is averageing over time axis
        

    def read_data(self):
        with open(self.name) as f:
            lines = f.readlines()
        raw_data = lines[3:]
        n_chunk = int(re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", raw_data[0])[1])
        n_data = len(raw_data)//(n_chunk+1)
        del raw_data[0:len(raw_data):(n_chunk + 1)]
        temp_data = np.loadtxt(raw_data)[:, [1, 3]]
        temp = temp_data[:,1].reshape((n_chunk,n_data), order ='F')
        data = np.concatenate(((temp_data[:n_chunk,0])[:, np.newaxis], temp),axis=1)
        return data
    
    def ave_data(self):
        ave_data = np.zeros([len(self.data),2])
        ave_data[:,0] = self.data[:,0]
        ave_data[:,1] = np.mean(self.data[:,1:], axis = 1)
        return ave_data

                
    def NEMD_Correlation(self):
        temp_data = self.data
        para_arr = np.zeros([len(temp_data[0][1:]),5, 3])
        for i in range(1,len(temp_data[0])):
            Mean_Temp = temp_data[:,[0,i]]
            list1 = Mean_Temp[len(Mean_Temp)//20 : len(Mean_Temp)//20*8]
            list2 = Mean_Temp[len(Mean_Temp)//20*12 : len(Mean_Temp)//20*19]
            para_arr[i - 1,:,0] = np.array(stats.linregress(list1[:,0], list1[:,1]))
            para_arr[i - 1,:,1] = np.array(stats.linregress(list2[:,0], list2[:,1]))
            para_arr[i - 1,:,2] = np.abs(para_arr[i - 1,:,:2]).mean(1)
            x1 = np.linspace(list1[0,0],list1[-1,0],100)
            x2 = np.linspace(list2[0,0],list2[-1,0],100)
            y1 = para_arr[i - 1,0,0]*x1+para_arr[i - 1,1,0]
            y2 = para_arr[i - 1,0,1]*x2+para_arr[i - 1,1,1]
            #plt.figure()
            plt.plot(x1,y1,color = 'blue',linewidth = 2)
            plt.plot(x2,y2,color = 'blue',linewidth = 2)
            plt.scatter(Mean_Temp[:,0], Mean_Temp[:,1], s = 1)
            plt.xlabel('z ($\AA$)')
            plt.xlim()
            plt.ylabel('Temperature(K)')
            #plt.xticks([i*5 for i in range(0,2*n/5)])
            plt.grid()

        return para_arr
    
    
class TEMP_funs():
    
    def data(n_mol, data):
        result = np.zeros([n_mol, 3])
        for i, val_i in enumerate(data):
            result[i] =  np.array(re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", \
                  val_i)).astype('float')[4:7]
        return result
        
    
'''   
class phononfile(Lammps):

    def __init__(self, file_name):
        self.name = file_name
        self.data = self.read_from_log()[1]
        self.ave_data = self.ave_log()
                                 
    
    
    def read_from_log(self):
        f = open(self.name)
        lines = f.readlines() 
        f.close()
        index = np.zeros([1])
        for i, val_i in enumerate(lines):
            if val_i == '# qx\t qy \t qz \t\t Phi(q)\n':
                index = np.append(index,i)
        index += 1
        index = index.astype(int)
        for i in index[1:]:
            
        return caption, data_log
'''