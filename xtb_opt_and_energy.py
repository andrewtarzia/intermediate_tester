#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Perform xTB optimisation and energy calculations on all structures.

Author: Andrew Tarzia

Date Created: 23 Feb 2020

"""

import stk
from os.path import exists
import glob


def optimize_structure(name, mol, solvent):
    print(f'doing xtb opt of {name}')
    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent

    print('doing XTB optimisation')
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'xtb_data/opt_{name}',
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

    return mol


def get_energy(name, mol, solvent):
    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent
    print(f'getting energy of {name}')
    xtb_energy = stk.XTBEnergy(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'xtb_data/ey_{name}',
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
    return energy, free_energy


def main():
    # Optimse and get energy of water.
    mol_file = 'water.mol'
    name = mol_file.replace('.mol', '')
    output_file = mol_file.replace('.mol', '_xtb.mol')
    ey_file = mol_file.replace('.mol', '_xtb.ey')
    fey_file = mol_file.replace('.mol', '_xtb.fey')

    if not exists(output_file):
        # Read molecule into stk.
        mol = stk.BuildingBlock.init_from_file(mol_file)
        # Optimise mol structure.
        mol = optimize_structure(
            name,
            mol,
            solvent=('chcl3', 'verytight')
        )
        # Output.
        mol.write(output_file)

    mol = stk.BuildingBlock.init_from_file(output_file)

    if not exists(ey_file) or not exists(fey_file):
        # Determine xtb energy.
        energy, free_energy = get_energy(
            name,
            mol,
            solvent=('chcl3', 'verytight')
        )
        # Output.
        with open(ey_file, 'w') as f:
            f.write(str(energy))
        with open(fey_file, 'w') as f:
            f.write(str(free_energy))

    # Do the rest.
    for mol_file in sorted(glob.glob('*a*_opt.mol')):
        name = mol_file.replace('.mol', '')
        output_file = mol_file.replace('.mol', '_xtb.mol')
        ey_file = mol_file.replace('.mol', '_xtb.ey')
        fey_file = mol_file.replace('.mol', '_xtb.fey')

        if not exists(output_file):
            # Read molecule into stk.
            mol = stk.BuildingBlock.init_from_file(mol_file)
            # Optimise mol structure.
            mol = optimize_structure(
                name,
                mol,
                solvent=('chcl3', 'verytight')
            )
            # Output.
            mol.write(output_file)

        mol = stk.BuildingBlock.init_from_file(output_file)

        if not exists(ey_file) or not exists(fey_file):
            # Determine xtb energy.
            energy, free_energy = get_energy(
                name,
                mol,
                solvent=('chcl3', 'verytight')
            )
            # Output.
            with open(ey_file, 'w') as f:
                f.write(str(energy))
            with open(fey_file, 'w') as f:
                f.write(str(free_energy))

        print('.................done.......................')


if __name__ == '__main__':
    main()
