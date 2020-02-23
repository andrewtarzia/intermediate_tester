#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Perform ETKDG conformer analysis.

Author: Andrew Tarzia

Date Created: 23 Feb 2020

"""

from rdkit.Chem import AllChem as rdkit
import stk
from os import mkdir
from os.path import exists, join
from itertools import product
import json
import glob


def calculate_NN_distance(stk_mol):
    # Assumes that there are only two N atoms in the molecules
    for atom1, atom2 in product(stk_mol.atoms, repeat=2):
        if atom1.id == atom2.id:
            continue
        if atom1.atomic_number == 7 and atom2.atomic_number == 7:
            print(atom1, atom2)
            NN_dist = stk_mol.get_atom_distance(atom1.id, atom2.id)
            print(NN_dist)
            return NN_dist


def optimise_conformer(name, mol, solvent, output_file):
    print(f'doing xtb opt of {name}')
    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent

    print('doing XTB optimisation')
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'{name}_opt',
        gfn_version=2,
        num_cores=6,
        opt_level='extreme',
        charge=0,
        num_unpaired_electrons=0,
        max_runs=5,
        electronic_temperature=300,
        calculate_hessian=True,
        unlimited_memory=True,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )
    xtb_opt.optimize(mol=mol)
    mol.write(output_file)

    return mol


def build_conformers(mol, dir):
    # Build conformers with ETKDG.
    etkdg = rdkit.ETKDG()
    etkdg.randomSeed = 1000
    cids = rdkit.EmbedMultipleConfs(
        mol=mol,
        numConfs=200,
        params=etkdg
    )

    # Save each conformer into directory.
    for cid in cids:
        filename = join(dir, f'conformer_{cid}.mol')
        rdkit.MolToMolFile(mol, filename=filename, confId=cid)


def load_f_energy(file):
    with open(file, 'r') as f:
        for line in f.readlines():
            return line.rstrip()


def calculate_f_energy(name, mol, solvent, ey_file, fey_file):
    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent
    print(f'getting energy of {name}')
    xtb_energy = stk.XTBEnergy(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'{name}_ey',
        num_cores=6,
        charge=0,
        num_unpaired_electrons=0,
        electronic_temperature=300,
        unlimited_memory=True,
        calculate_free_energy=True,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )
    energy = xtb_energy.get_energy(mol)
    free_energy = xtb_energy.total_free_energies[mol]
    with open(ey_file, 'w') as f:
        f.write(str(energy))
    with open(fey_file, 'w') as f:
        f.write(str(free_energy))


def main():

    amine_files = [
        'ami1_opt_xtb.mol',
        'ami2_opt_xtb.mol',
        'ami3_opt_xtb.mol',
        'ami4_opt_xtb.mol',
    ]

    print(amine_files)
    results = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
    for ami in amine_files:
        res = {}
        amine = ami.replace('_opt_xtb.mol', '')
        ey_out_file = ami.replace('_opt_xtb.mol', 'conf_results.out')
        ami_dir = ami.replace('_opt_xtb.mol', '_confs')

        # Do calcs, only if the two energy files do not already exist.
        if exists(ey_out_file):
            with open(ey_out_file, 'r') as f:
                lst = json.load(f)
                results[amine] = lst
            print(results)
            input()
            continue

        # Read structure into rdkit.
        rdk_mol = rdkit.MolFromMolFile(ami)
        rdk_mol = rdkit.AddHs(rdk_mol)
        print(f'doing {amine}')
        # Write directory.
        if not exists(ami_dir):
            mkdir(ami_dir)

        # Get N conformers.
        print(f'building conformers of {amine}')
        build_conformers(mol=rdk_mol, dir=ami_dir)

        # For each conformer, iterate.
        print(f'getting energies of conformers of {amine}')
        for conf_file in glob.glob(join(ami_dir, 'conformer*mol')):
            if 'opt' in conf_file:
                continue
            print(conf_file)
            opt_conf_file = conf_file.replace('.mol', '_opt.mol')
            ey_file = conf_file.replace('.mol', '.ey')
            fey_file = conf_file.replace('.mol', '.fey')
            opt_ey_file = conf_file.replace('.mol', '_opt.ey')
            opt_fey_file = conf_file.replace('.mol', '_opt.fey')
            conf_id = conf_file.replace(ami_dir+'/', '')
            conf_id = conf_id.split('_')[1].replace('.mol', '')
            print(ey_file, conf_id)

            # Load in stk object.
            conf = stk.BuildingBlock.init_from_file(conf_file)
            print(conf)

            # Get xTB energy, optimise with xTB, and get energy again.
            # Read in opt and unopt energies.
            if not exists(ey_file) or not exists(fey_file):
                # Calculate energy.
                calculate_f_energy(
                    name=join(ami_dir, f'unoptey_{conf_id}'),
                    mol=conf,
                    solvent=('chcl3', 'verytight'),
                    ey_file=ey_file,
                    fey_file=fey_file
                )
            # Load energy.
            f_energy = load_f_energy(fey_file)

            if not exists(opt_ey_file) or not exists(opt_fey_file):
                if not exists(opt_conf_file):
                    # Optimise conformer.
                    optimise_conformer(
                        name=join(ami_dir, f'opt_{conf_id}'),
                        mol=conf,
                        solvent=('chcl3', 'verytight'),
                        output_file=opt_conf_file
                    )
                opt_conf = stk.BuildingBlock.init_from_file(
                    opt_conf_file
                )
                # Calculate energy.
                calculate_f_energy(
                    name=join(ami_dir, f'optey_{conf_id}'),
                    mol=opt_conf,
                    solvent=('chcl3', 'verytight'),
                    ey_file=opt_ey_file,
                    fey_file=opt_fey_file
                )
            # Load energy.
            opt_f_energy = load_f_energy(opt_fey_file)
            print(fey_file, opt_fey_file)
            print(f_energy)
            print(opt_f_energy)

            # Calculate N-N distance.
            conf = stk.BuildingBlock.init_from_file(
                conf_file,
                use_cache=False
            )
            NN_dist = calculate_NN_distance(conf)
            opt_conf = stk.BuildingBlock.init_from_file(
                opt_conf_file,
                use_cache=False
            )
            opt_NN_dist = calculate_NN_distance(opt_conf)
            res[conf_id] = {
                'f_energy': f_energy,
                'opt_f_energy': opt_f_energy,
                'NN_dis': NN_dist,
                'opt_NN_dis': opt_NN_dist,
            }
            print(res)

        # Save results to JSON file and results dict.
        results[amine] = res
        with open(ey_out_file, 'w') as f:
            json.dump(res, f)

        print('.................done.......................')