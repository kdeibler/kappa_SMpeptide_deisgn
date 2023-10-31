############################################################################################
## some utilities useful for peptide generation, analysis, etc. Parisa Hosseinzadeh, 2018 ##
############################################################################################
from __future__ import print_function

import argparse
import glob

import math
import os
import sys
import csv
import itertools
from itertools import chain
from multiprocessing import Manager, Pool, freeze_support, cpu_count
import numpy as np
import rmsd
import re
import random

from pyrosetta import *
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.protocols.protein_interface_design.filters import HbondsToResidueFilter
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import NeighborhoodResidueSelector
from pyrosetta.rosetta.protocols.hbnet import UnsatSelector
from pyrosetta.rosetta.core.select.residue_selector import ResiduePDBInfoHasLabelSelector
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
from pyrosetta.rosetta.protocols.simple_ddg import DdgFilter
from rosetta import *
#from crosslinker_mover import PeptideCrosslinkerMover
from pyrosetta.bindings.utility import bind_method
from pyrosetta.rosetta.core.pack.rotamer_set import bb_independent_rotamers
from pyrosetta.rosetta.core.conformation import Residue

_DEBUG=False
_CHECK=False

#some basic necessary stuff
L_res=["GLY","ALA","CYS","MET","VAL","LEU","ILE","ASP","GLU","ASN","GLN","THR","SER","TYR","TRP","PHE","LYS","ARG","HIS","PRO"]
D_res=["GLY","DALA","DCYS","DMET","DVAL","DLEU","DILE","DASP","DGLU","DASN","DGLN","DTHR","DSER","DTYR","DTRP","DPHE","DLYS","DARG","DHIS","DPRO"]

#getting rotamers for a given residue
@bind_method(Residue)
def get_rotamers(self):
    if self.type().is_d_aa():
        rotamers = [pyrosetta.Vector1([-chi for chi in rotamer.chi()]) for rotamer in bb_independent_rotamers(self.type())]
    else:
        rotamers = [rotamer.chi() for rotamer in bb_independent_rotamers(self.type())]
    
    return rotamers

#random selection of a rotamer for finding rotamers
@bind_method(Residue)
def set_random_rotamer(self):
    rotamers = self.get_rotamers()
    one_rot=rotamers[random.randint(0,len(rotamers))]
    #self.set_chi(one_rot)
    if _DEBUG:
        print ("-----this is the random rotamer selected", one_rot,"out of total", len(rotamers))
    #setting the random rotamer
    for i in range(1,len(one_rot)+1):
        self.set_chi(i,one_rot[i])

#removing variant types (unfortunately necessary to add residues. Er)
def variant_remove(p):

    for ir in range(1,p.size()+1):
        if ( p.residue(ir).has_variant_type(core.chemical.UPPER_TERMINUS_VARIANT)):
            core.pose.remove_variant_type_from_pose_residue( p, core.chemical.UPPER_TERMINUS_VARIANT, ir)
        if ( p.residue(ir).has_variant_type(core.chemical.LOWER_TERMINUS_VARIANT)):
            core.pose.remove_variant_type_from_pose_residue( p, core.chemical.LOWER_TERMINUS_VARIANT, ir)
        if ( p.residue(ir).has_variant_type(core.chemical.CUTPOINT_LOWER)):
            core.pose.remove_variant_type_from_pose_residue( p, core.chemical.CUTPOINT_LOWER, ir)
        if ( p.residue(ir).has_variant_type(core.chemical.CUTPOINT_UPPER)):
            core.pose.remove_variant_type_from_pose_residue( p, core.chemical.CUTPOINT_UPPER, ir)

#n is number added, p is pose, c_ is boolean whether to add to c-term. n_ bool for n-term
#careful. Assumption here is that you're passing only one chain. a is residue
def extend(n, p, a, c_, n_):
    
    #first removing variatn types
    variant_remove(p)

    if (n_):
        for i in range(0, n):
            variant_remove(p)
            p.prepend_polymer_residue_before_seqpos(a, 1, True)
        if _DEBUG:
            p.dump_pdb('prepend_test.pdb')
        for res_no in range(1,n+1):
            p.set_omega(res_no,180.)
    if (c_):
        for i in range(0, n):
            variant_remove(p)
            p.append_residue_by_bond(a, True)
            if _DEBUG:
                p.dump_pdb('append_test.pdb')
        for res_no in range((p.size()-n)-1, p.size()+1):
            p.set_omega(res_no,180.)

    if (_DEBUG):
        p.dump_pdb('extended.pdb')


# bin sampler for added residues
# is pose and res is the residue that is being sampled and phi, psi are phi, psi values
def bin_sample(p,res, phi, psi):

    p.set_phi(res, phi)
    p.set_psi(res, psi)

    if (_DEBUG):
        p.dump_pdb('set_{}_{}.pdb'.format(phi,psi))

# cart-min relax
# is the chain number that has peptide
def relax(pose,c,cart):

    indeces=[]
    for resNo in range(1,pose.size()):
        if (pose.residue(resNo).chain() == c):
            indeces.append(resNo)
    my_score=pyrosetta.get_score_function()
    if cart:
        #set cart weigths
        my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded, 1.0)
        my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_angle, 1.0)
        my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_length, 1.0)
        my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_ring, 1.0)
        my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_torsion, 1.0)
        my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_proper, 1.0)
        my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_improper, 1.0)
    #set metal constraints
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.metalbinding_constraint, 1.0)
    frm = pyrosetta.rosetta.protocols.relax.FastRelax(my_score)
    frm.ramp_down_constraints(False)
    if cart:
        #set cart to be true
        frm.cartesian(True)
    else:
        frm.cartesian(False)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_jump(1,False)
    #setting bb false so that the phi_psi stays as I want
#mm.set_bb_true_range(indeces[0], indeces[-1])
    mm.set_chi_true_range(indeces[0], indeces[-1])
    frm.set_movemap(mm)
    frm.apply(pose)
    pyrosetta.rosetta.core.pose.remove_lower_terminus_type_from_pose_residue(pose, 1)
    if (_DEBUG):
        pose.dump_pdb('relax.pdb')


# cyclize mover
def cyclizer(p):

    pcm=protocols.cyclic_peptide.PeptideCyclizeMover()
    pcm.apply(p)

#get just the peptide, when peptide is chain c
def extractor (p_in,c):

    p=Pose()
    for resNo in range(1,p_in.size()):
        if (p_in.residue(resNo).chain() == c):
            p.append_residue_by_bond(p_in.residue(resNo),False)

    return p

# my hacky satisfier
# this only gives value of nearby atoms for carbonyls that are unsatisfied
def unsat_count(pose,chain):
    
    checker=0
    indeces=[]
    for resNo in range(1,pose.size()):
        if (pose.residue(resNo).chain() == chain):
            indeces.append(resNo)
    to_check = ','.join(map(str, indeces))
    index_sel=ResidueIndexSelector(to_check)
    unsats=UnsatSelector()
    unsats.set_consider_mainchain_only(False)
    unsats.set_scorefxn(pyrosetta.get_score_function())
    unsat_chain=AndResidueSelector(index_sel,unsats)
    subset=unsat_chain.apply(pose)
    all=[]
    for i in range(1,pose.size()):
        nearby=[]
        if (subset[i]):
            r1=pyrosetta.rosetta.core.conformation.Residue(pose.residue(i))
            for j in range(1,pose.size()):
                if (abs(j-i) > 1):
                    #r2=pyrosetta.rosetta.core.conformation.Residue(pose.residue(j))
                    for at in range (1,pose.residue(j).natoms()):
                        d_sq=r1.xyz('O').distance_squared(pose.xyz(core.id.AtomID(at,j)))
                        if (d_sq < 25.01):
                            nearby.append([i,j])
        if (len(nearby) > 30):
            all.append(i)
    
    return len(all)

# score
# p is pose and c is chain
def metrics(p,c):

    #you need to score so that Rosetta can give you the scores you want.
    scf=pyrosetta.get_score_function()
    scf(p)
    # getting score of the whole pose
    tot_energy_value=p.energies().total_energies()[core.scoring.score_type_from_name("total_score")]
    # getting score of the peptide chain
    pep_energy_value=0
    for resNo in range(1,p.size()+1):
        if (p.residue(resNo).chain() == c):
        # print (p.residue(resNo).seqpos())
            pep_energy_value+=p.energies().residue_total_energies(p.residue(resNo).seqpos())[core.scoring.score_type_from_name("total_score")]
    # getting shape complemetarity
    sc=ShapeComplementarityFilter()
    sc.jump_id(1)
    sc_value=sc.report_sm(p)
    #getting ddg_norepack
    ddg1=DdgFilter(0,pyrosetta.get_score_function(),1)
    ddg1.repack(0)
    ddg_value1=ddg1.report_sm(p)
    # getting ddg_repack
    ddg2=DdgFilter(0,pyrosetta.get_score_function(),1)
    ddg2.repack(0)
    ddg_value2=ddg2.report_sm(p)
    #buried unsat
    unsat_value=unsat_count(p,c)

    metrics=[("total energy",tot_energy_value),("peptide energy",pep_energy_value),("shape complementarity",sc_value),("ddg no_repack",ddg_value1), ("ddg repack",ddg_value2), ("number unsats",unsat_value)]

    return metrics

# csv writer
def add_to_score(metric,res_inf,res):
    
    everything=[]
    for info in res_inf:
        everything.append(info)
    for energy in metric:
        everything.append(energy[1])

    with open(r'scores_{}_cont.csv'.format(res), 'a') as f:
        writer = csv.writer(f, delimiter=' ')
        writer.writerow(everything)

# coordinate finder.this function gets x,y,z for an atom and returns an array corresponding to them
# ir is the residue number and ia is the atom name
def coord_find(p,ir,ia):

#name=re.compile('.HA ')
#if (len(re.findall(name, ia)) == 0):
    coord_xyz=p.xyz(core.id.AtomID(p.residue(ir).atom_index(ia),ir))
    coord_arr=[]
    x=coord_xyz[0]
    y=coord_xyz[1]
    z=coord_xyz[2]
    coord_arr.append(x)
    coord_arr.append(y)
    coord_arr.append(z)
    
    if (_DEBUG):
        print (ia, coord_arr)
    
    return coord_arr

# finds the center
def find_cent(A):
    
    sumA=[0,0,0]
    for i in range(len(A)):
        sumA[0]=sumA[0]+A[i][0]
        sumA[1]=sumA[1]+A[i][1]
        sumA[2]=sumA[2]+A[i][2]
    
    for i in range(3):
        sumA[i]=sumA[i]/len(A)
    
    if (_DEBUG):
        print ("========= find center function==========")
        print ("the elements are", '\n'
               ,A,'\n'
               ,"and the center is:",'\n'
               ,sumA)
    
    return sumA


# alignment (added because I am dumb with fold trees)
def transform(orig,rotd,resi, r_end):

    p_scaff=[]
    p_targ=[]
    for atom in range(1,orig.residue(resi).natoms()+1):
        #if (not (orig.residue(resi).atom_is_backbone(atom))):
        if (not orig.residue(resi).atom_is_hydrogen(atom)):
            p_scaff.append(coord_find(rotd,r_end,rotd.residue(r_end).atom_name(atom)))
            p_targ.append(coord_find(orig,resi,orig.residue(resi).atom_name(atom)))

    #step1: moving scaffold to the center
    T=find_cent(p_scaff)
    plusv=numeric.xyzVector_double_t(-1*T[0],-1*T[1],-1*T[2])
    #does not rotate
    noR=numeric.xyzMatrix_double_t.cols(1,0,0,0,1,0,0,0,1)
    rotd.apply_transform_Rx_plus_v(noR,plusv)
    if (_DEBUG):
        print ("============scaffold translation, step 1===============")
        print (T,noR)
        rotd.dump_pdb('translate_scaffold.pdb')

    #Step1': get the coordinates of target at the center
    T_targ=find_cent(p_targ)
    v_targ=numeric.xyzVector_double_t(-1*T_targ[0],-1*T_targ[1],-1*T_targ[2])
    orig.apply_transform_Rx_plus_v(noR,v_targ)
    if (_DEBUG):
        print ("============target translation, step 1===============")
        print (T_targ,noR)
        orig.dump_pdb('translate_target.pdb')

    #need to re-load the matrix now because the pose has changed
    p_scaff_new=[]
    p_targ_new=[]
    for atom in range(1,orig.residue(resi).natoms()+1):
        #if (not (orig.residue(resi).atom_is_backbone(atom))):
        if (not orig.residue(resi).atom_is_hydrogen(atom)):
            p_scaff_new.append(coord_find(rotd,r_end,rotd.residue(r_end).atom_name(atom)))
            p_targ_new.append(coord_find(orig,resi,orig.residue(resi).atom_name(atom)))


    #Step 2: get the rotation matrix
    #the magic of libraries
    #V=rmsd.kabsch(p_targ_new,p_scaff_new)
    semi_V=rmsd.kabsch(p_scaff_new,p_targ_new)
    V=np.linalg.inv(semi_V)
    if (_DEBUG):
        print ("============scaffold rotation, step 2===============")
        print ("the transformation matrix is",'\n', V)

    #Rotate the pose
    Rx=numeric.xyzMatrix_double_t.cols(V[0][0],V[1][0],V[2][0],V[0][1],V[1][1],V[2][1],V[0][2],V[1][2],V[2][2])
    noT=numeric.xyzVector_double_t(0,0,0)
    if (_DEBUG):
        print ("the old rmsd is:", rmsd.kabsch_rmsd(p_scaff, p_targ))

    #moving the pose
    rotd.apply_transform_Rx_plus_v(Rx,noT)
    if (_DEBUG):
        rotd.dump_pdb('rotate_scaffold.pdb')

    #Step3: translate the pose back to target (both the new and the original)
    scaff_trans=numeric.xyzVector_double_t(T_targ[0],T_targ[1],T_targ[2])
    rotd.apply_transform_Rx_plus_v(noR,scaff_trans)
    orig.apply_transform_Rx_plus_v(noR,scaff_trans)
    if (_DEBUG):
        rotd.dump_pdb('back_to_orig.pdb')
        print ("8888888", scaff_trans)

    #generating final set
    p_scaff_final=[]
    for atom in range(1,orig.residue(resi).natoms()+1):
        #if (not (orig.residue(resi).atom_is_backbone(atom))):
        if (not orig.residue(resi).atom_is_hydrogen(atom)):
            p_scaff_final.append(coord_find(rotd,r_end,rotd.residue(r_end).atom_name(atom)))
    if (_DEBUG):
        print ("the new rmsd is:", rmsd.kabsch_rmsd(p_scaff_final, p_targ))

# rotating a point around a line (given by l1 and l2) with a given degree
# line is given in the form of (a,b,c,u,v,w)
def rotater(point,line,ang):

    a=line[0]
    b=line[1]
    c=line[2]
    u=line[3]
    v=line[4]
    w=line[5]

    x=point[0]
    y=point[1]
    z=point[2]

    L=u*u+v*v+w*w

    new_x=((a*(v*v+w*w)-u*(b*v+c*w-u*x-v*y-w*z))*(1-math.cos(ang))+L*x*math.cos(ang)+(math.sqrt(L))*(-c*v+b*w-w*y+v*z)*math.sin(ang))/L
    
    new_y=((b*(u*u+w*w)-v*(a*u+c*w-u*x-v*y-w*z))*(1-math.cos(ang))+L*y*math.cos(ang)+(math.sqrt(L))*(c*u-a*w+w*x-u*z)*math.sin(ang))/L
    
    new_z=((c*(u*u+v*v)-w*(a*u+b*v-u*x-v*y-w*z))*(1-math.cos(ang))+L*z*math.cos(ang)+(math.sqrt(L))*(-b*u+a*v-v*x+u*y)*math.sin(ang))/L
    
    new_point=[new_x,new_y,new_z]

    return new_point

# find a line from two points
def line_finder(p1,p2):

    vec_x=p1[0]-p2[0]
    vec_y=p1[1]-p2[1]
    vec_z=p1[2]-p2[2]

    return [p1[0],p1[1],p2[2],vec_x,vec_y,vec_z]

# mutates, minimizes
def packer (p,resn,resi,c):

    indeces=[]
    for resNo in range(1,p.size()+1):
        if (p.residue(resNo).chain() == c):
            indeces.append(resNo)
    #step one is to add the residue you want
    mut=protocols.simple_moves.MutateResidue()
    mut.set_res_name(resn)
    mut.set_target(resi)
    mut.set_preserve_atom_coords(True)
    mut.apply(p)

    #step two is to minimize with a cart_min probably
    my_score=pyrosetta.get_score_function()
    #set cart weigths
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_angle, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_length, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_ring, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_torsion, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_proper, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_improper, 1.0)
    min=protocols.minimization_packing.MinMover()
    min.score_function(my_score)
    min.cartesian(True)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_jump(1,False)
    mm.set_bb_true_range(indeces[0], indeces[-1])
    mm.set_chi_true_range(indeces[0], indeces[-1])
    min.set_movemap(mm)
    min.apply(p)


