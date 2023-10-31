#############################################################
### This script allows for rough docking of pre-generated ###
### backbones to a residue of interest.                   ###
#############################################################

from pyrosetta import *
from rosetta import *
import sys, os, glob, math, random, subprocess
import numpy as np
import rosetta.protocols.rigid as rigid_moves
from rosetta.core.scoring.constraints import *
from rosetta.core.scoring.func import HarmonicFunc
from rosetta.protocols.cyclic_peptide import DeclareBond
import rmsd
import math

#global variables that need to be set as arguments upon submission
#arg1 => initial pdb complex to dock into
#arg2 => provide the location of pregenerated thioether backbones
comp=sys.argv[1]
thioether='thioether_6mers/{}.pdb'.format(sys.argv[2])
min_e=250

# initiating pyrosetta
init('-override_rsd_type_limit -ex1 -use_input_sc -ignore_unrecognized_res -mute all -extra_res_fa CYY.params CVV.params -detect_disulf false')

# defining score functions required
sfxn = create_score_function('ref2015')
sfxn_soft = create_score_function('ref2015_soft')
sfxn_mm_std = create_score_function('mm_std')

sfxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_dun_dev, 0.0)
sfxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_dun_rot, 0.0)
sfxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_dun_semi, 0.0)

# getting residue instances
create_residue = rosetta.core.conformation.ResidueFactory.create_residue
chm = rosetta.core.chemical.ChemicalManager.get_instance()
rts = chm.residue_type_set( 'fa_standard' )

# reading in the pdb file
pose_AB_complex = pose_from_pdb(comp)

#splitting by chain
AB_chains = pose_AB_complex.split_by_chain()

#defining peptide and protein. Assumes peptide is chain B
prot=Pose(AB_chains[1])
pep=Pose(AB_chains[2])

#reading in the thioether template
temp = pose_from_pdb(thioether)

#getting the chi values for the thioether residue
chis = []
for i in range(1, 5):
    # Adds chi values to last residue
    chis.append(temp.chi(i, temp.size()))
                
# makes the first residue CYY so that they can be overlaid later
pep.replace_residue(1, rosetta.core.conformation.ResidueFactory.create_residue( rts.name_map('CYY')), True)

#sets the chi for the last residue
for i in range(1, 5):
    pep.set_chi(i, 1, chis[i-1])

# transforming scaffold onto the peptide
transform(pep,temp, 1, temp.size())

# putting everything back together
final_complex=Pose()
new_pep=Pose()
variant_remove(prot)

#adding protein as chain A
for resNo in range(1,prot.size()+1):
    final_complex.append_residue_by_bond(prot.residue(resNo),False)
variant_remove(temp)

#adding scaffold as beginning of chain B
for resNo in range(1,temp.size()+1):
    if resNo == 1:
        new_pep.append_residue_by_jump(temp.residue(resNo),prot.size(),'','',True)
    else:
        new_pep.append_residue_by_bond(temp.residue(resNo),False)
variant_remove(pep)

#adding the peptide at the end as a tail
for resNo in range (2,pep.size()+1):
    new_pep.append_residue_by_bond(pep.residue(resNo),False)

# need to declare bond between the CYY and sth
new_pep.conformation().declare_chemical_bond(1, 'N', temp.size(), 'CE')

for resNo in range(1, new_pep.size()+1):
    if resNo == 1:
        final_complex.append_residue_by_jump(new_pep.residue(resNo),prot.size(),'','',True)
    else:
        final_complex.append_residue_by_bond(new_pep.residue(resNo),False)

#need to minimize
#setting up the movemap factory
movemap = MoveMap()
movemap.set_chi(True)
movemap.set_bb(True)
movemap.set_jump(True)
minmover_mm = pyrosetta.rosetta.protocols.minimization_packing.MinMover(movemap, sfxn_mm_std, 'linmin', 0.001, True)
minmover_beta = pyrosetta.rosetta.protocols.minimization_packing.MinMover(movemap, sfxn, 'linmin', 0.001, True)

#setting up the constraints
dist_cst = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(pyrosetta.rosetta.core.id.AtomID(final_complex.residue(prot.size()+temp.size()).atom_index('CE'), prot.size()+temp.size()), pyrosetta.rosetta.core.id.AtomID(final_complex[prot.size()+1].atom_index('N'), 1), pyrosetta.rosetta.core.scoring.func.HarmonicFunc(1.33, 0.01))
tors_cst = pyrosetta.rosetta.core.scoring.constraints.DihedralConstraint(pyrosetta.rosetta.core.id.AtomID(final_complex.residue(prot.size()+temp.size()).atom_index('CD'), prot.size()+temp.size()), pyrosetta.rosetta.core.id.AtomID(final_complex.residue(prot.size()+temp.size()).atom_index('CE'), prot.size()+temp.size()), pyrosetta.rosetta.core.id.AtomID(final_complex[prot.size()+1].atom_index('N'), 1), pyrosetta.rosetta.core.id.AtomID(final_complex[prot.size()+1].atom_index('CA'), 1), pyrosetta.rosetta.core.scoring.func.CircularHarmonicFunc(3.14, 0.005))
ang1_cst = pyrosetta.rosetta.core.scoring.constraints.AngleConstraint(pyrosetta.rosetta.core.id.AtomID(final_complex.residue(prot.size()+temp.size()).atom_index('CD'), prot.size()+temp.size()), pyrosetta.rosetta.core.id.AtomID(final_complex.residue(prot.size()+temp.size()).atom_index('CE'), prot.size()+temp.size()), pyrosetta.rosetta.core.id.AtomID(final_complex[prot.size()+1].atom_index('N'), 1), pyrosetta.rosetta.core.scoring.func.HarmonicFunc(2.01, 0.03))
ang2_cst = pyrosetta.rosetta.core.scoring.constraints.AngleConstraint(pyrosetta.rosetta.core.id.AtomID(final_complex.residue(prot.size()+temp.size()).atom_index('CE'), prot.size()+temp.size()), pyrosetta.rosetta.core.id.AtomID(final_complex[prot.size()+1].atom_index('N'), 1), pyrosetta.rosetta.core.id.AtomID(final_complex[prot.size()+1].atom_index('CA'), 1), pyrosetta.rosetta.core.scoring.func.HarmonicFunc(2.14, 0.03))
final_complex.add_constraint(dist_cst)
final_complex.add_constraint(tors_cst)
final_complex.add_constraint(ang1_cst)
final_complex.add_constraint(ang2_cst)

#applying minimization
for i in range(0, 10):
    #minmover_mm.apply(final_complex)
    minmover_beta.apply(final_complex)

# get a score and compared to original score cutoff at beginning of script
sc = sfxn(final_complex)
if sc < min_e:
    final_complex.dump_pdb('minimized_' + str(sc) + '.pdb')
