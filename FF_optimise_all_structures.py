#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Perform FF optimisation and conformer search with OPLS.

Author: Andrew Tarzia

Date Created: 23 Feb 2020

"""

import stk
from os.path import exists
import glob


def settings():
    Settings = {
        'output_dir': None,
        'timeout': None,
        'force_field': 16,
        'temperature': 700,  # K
        'conformers': 200,  # change from 10000
        'time_step': 0.5,  # fs
        'eq_time': 10,  # ps
        # ps -- 50 ns changed from 100 ns
        'simulation_time': 5000,
        'maximum_iterations': 2500,
        'minimum_gradient': 0.05,
        'use_cache': False
    }
    return Settings


def optimize_structure(name, mol):
    Settings = settings()
    # restricted=False optimization with OPLS forcefield by default
    print(f'doing FF opt of {name}')
    ff = stk.MacroModelForceField(
        macromodel_path='/home/atarzia/software/schrodinger_install',
        output_dir=f'FF_data/{name}_FFout',
        restricted=False
    )
    # MD process - run MD, collect N conformers, optimize each,
    # return lowest energy conformer
    # print(f'doing MD opt of {name}')
    # md = stk.MacroModelMD(
    #     macromodel_path='/home/atarzia/software/schrodinger_install',
    #     output_dir=f'MD_data/{name}_MDout',
    #     timeout=Settings['timeout'],
    #     force_field=Settings['force_field'],
    #     temperature=Settings['temperature'],
    #     conformers=Settings['conformers'],
    #     time_step=Settings['time_step'],
    #     eq_time=Settings['eq_time'],
    #     simulation_time=Settings['simulation_time'],
    #     maximum_iterations=Settings['maximum_iterations'],
    #     minimum_gradient=Settings['minimum_gradient'],
    #     use_cache=Settings['use_cache']
    # )
    # seq = stk.Sequence(ff, md)
    ff.optimize(mol=mol)
    return mol


def main():
    for mol_file in sorted(glob.glob('*a*.mol')):
        name = mol_file.replace('.mol', '')
        if 'opt' in name or 'xtb' in name:
            continue
        output_file = mol_file.replace('.mol', '_opt.mol')
        if exists(output_file):
            continue

        # Read molecule into stk.
        mol = stk.BuildingBlock.init_from_file(mol_file)
        print(mol)

        # Optimise mol structure.
        mol = optimize_structure(name, mol)

        # Output.
        mol.write(output_file)


if __name__ == '__main__':
    main()
