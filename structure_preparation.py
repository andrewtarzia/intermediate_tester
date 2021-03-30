#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Perform optimisation and energy calculations on all structures.

Author: Andrew Tarzia

Date Created: 30 Mar 2021

"""

from os.path import exists
import pandas as pd

import stk
import stko

from reactions import landscape
from utilities import read_ey


def optimize_structure(name, mol, solvent):

    # OPLS FF
    print(f'doing opls opt of {name}')
    ff = stko.MacroModelForceField(
        macromodel_path='/home/atarzia/software/schrodinger_install',
        output_dir=f'FF_data/{name}_FFout',
        restricted=False,
    )
    mol = ff.optimize(mol)

    # OPLS MM conf search
    # MD process - run MD, collect N conformers, optimize each,
    # return lowest energy conformer
    print(f'doing MD opt of {name}')
    md = stk.MacroModelMD(
        macromodel_path='/home/atarzia/software/schrodinger_install',
        output_dir=f'MD_data/{name}_MDout',
        timeout=None,
        force_field=16,  # OPLS3e
        temperature=700,  # K
        conformers=1000,
        time_step=0.5,  # fs
        eq_time=10,  # ps
        simulation_time=5000,  # ps
    )
    mol = md.optimize(mol)

    # xTB opt.
    print(f'doing xtb opt of {name}')
    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent

    xtb_opt = stko.XTB(
        xtb_path='/home/atarzia/software/xtb-6.3.2/bin/xtb',
        output_dir=f'xtb_data/opt_{name}',
        gfn_version=2,
        num_cores=6,
        opt_level='extreme',
        max_runs=5,
        calculate_hessian=True,
        unlimited_memory=True,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )
    mol = xtb_opt.optimize(mol)

    return mol


def get_energy(name, mol, solvent):

    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent

    print(f'getting energy of {name}')
    xtb_energy = stko.XTBEnergy(
        xtb_path='/home/atarzia/software/xtb-6.3.2/bin/xtb',
        output_dir=f'xtb_data/ey_{name}',
        num_cores=6,
        unlimited_memory=True,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )
    energy = xtb_energy.get_energy(mol)
    return energy


def optimisation_sequence(mol_file):
    """
    Sequence for all molecules to undergo.

    """

    name = mol_file.replace('.mol', '')
    output_file = mol_file.replace('.mol', '_opt.mol')
    if not exists(output_file):
        # Read molecule into stk.
        mol = stk.BuildingBlock.init_from_file(mol_file)
        # Optimise mol structure.
        mol = optimize_structure(
            name=name,
            mol=mol,
            solvent=None,
        )
        # Output.
        mol.write(output_file)


def energy_calculation(mol_file, solvent):
    """
    Energy calculation sequence for all moelcules.

    """

    name = mol_file.replace('.mol', '')
    if solvent is None:
        ey_file = mol_file.replace('.mol', '_gas.ey')
    else:
        ey_file = mol_file.replace('.mol', f'_{solvent[0]}.ey')
    mol = stk.BuildingBlock.init_from_file(mol_file)

    if not exists(ey_file):
        # Determine xtb energy.
        energy = get_energy(
            name=name,
            mol=mol,
            solvent=solvent,
        )
        # Output.
        with open(ey_file, 'w') as f:
            f.write(str(energy))
    else:
        energy = read_ey(ey_file)

    return energy


def main():
    gas_energy_dict = {}
    chcl3_energy_dict = {}

    # Optimse and get energy of water and othe precursors.
    precursors = ['water', 'alde1', 'ami1', 'ami2', 'ami3', 'ami4']
    for p in precursors:
        print(f'doing {p}')
        optimisation_sequence(mol_file=f'{p}.mol')
        energy = energy_calculation(
            mol_file=f'{p}_opt.mol',
            solvent=None,
        )
        gas_energy_dict[p] = energy  # a.u.
        energy = energy_calculation(
            mol_file=f'{p}_opt.mol',
            solvent=('chcl3', 'verytight'),
        )
        chcl3_energy_dict[p] = energy  # a.u.
        print('..done..')

    print(gas_energy_dict, chcl3_energy_dict)

    # Do the rest.
    lscp = landscape()
    for intermediate in lscp:
        iname = intermediate['intname']
        print(f'doing {iname}')
        optimisation_sequence(mol_file=f'{iname}.mol')
        energy = energy_calculation(
            mol_file=f'{p}_opt.mol',
            solvent=None,
        )
        gas_energy_dict[p] = energy  # a.u.
        energy = energy_calculation(
            mol_file=f'{p}_opt.mol',
            solvent=('chcl3', 'verytight'),
        )
        chcl3_energy_dict[p] = energy  # a.u.
        print('..done..')

    df = pd.DataFrame({
        'names': [i for i in sorted(gas_energy_dict.keys())],
        'energy_gas': [
            gas_energy_dict[i] for i in sorted(gas_energy_dict.keys())
        ],
        'energy_chcl3': [
            chcl3_energy_dict[i]
            for i in sorted(gas_energy_dict.keys())
        ],
    })
    df.to_csv('extracted_xtb_energies.csv', index=False)


if __name__ == '__main__':
    main()
