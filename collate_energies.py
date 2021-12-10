#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Collate output of Gaussian single point energy calculations.

Author: Andrew Tarzia

Date Created: 25 Oct 2021

"""

import sys
import glob
import pandas as pd
import re


def xtb_coll_fn(name, directory, solvent, method):
    # This is kind of a double up because structure_prep already
    # collates these energies, just loading from that.
    data = pd.read_csv(f'{directory}/extracted_xtb_energies.csv')
    col_name = f'energy_{solvent}'
    row = data[data['names'] == name]
    energy = float(row[col_name].iloc[0])
    return energy


def g19_spe_coll_fn(name, directory, solvent, method):

    output_file = f'{directory}/{name}_opt_{solvent}_{method}.log'
    with open(output_file, 'r') as f:
        for line in f.readlines():
            if 'A.U. after' in line:
                l = line.rstrip()
                if method == 'pbe':
                    energy = l.split('SCF Done:  E(RPBE1PBE) =')
                elif method == 'mp2':
                    energy = l.split('SCF Done:  E(RHF) =')
                energy = energy[-1].split('A.U. after')[0]
                return float(energy)  # a.u.


def orca_spe_coll_fn(name, directory, solvent, method):

    output_file = f'{directory}/o_{name}_opt_{method}.out'
    nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
    with open(output_file, 'r') as f:
        checked = False
        for line in f.readlines():
            if method == 'B97-3c':
                # Make sure you get the final one after opt.
                if 'FINAL ENERGY EVALUATION AT THE STATIONARY' in line:
                    checked = True
            elif method == 'MP2':
                # Was single point.
                checked = True

            if 'FINAL SINGLE POINT ENERGY' in line and checked:
                string = nums.search(line.rstrip()).group(0)
    return float(string)  # a.u.


def main():
    if len(sys.argv) != 2:
        print('csv file to write energies to.')
        sys.exit()
    else:
        output_ey_file = sys.argv[1]

    structures = sorted(glob.glob('*_opt.mol'))
    print(f'there are {len(structures)} structures')
    structure_titles = [i.replace('_opt.mol', '') for i in structures]

    _energy_methods = {
        'xtbgas': {
            'collation_function': xtb_coll_fn,
            'directory': './',
            'solvent': 'gas',
            'method': 'xtb'
        },
        'xtbchcl3': {
            'collation_function': xtb_coll_fn,
            'directory': './',
            'solvent': 'chcl3',
            'method': 'xtb'
        },
        'ds_pbe_gas': {
            'collation_function': g19_spe_coll_fn,
            'directory': 'xtb_spe_DFT',
            'solvent': 'gas',
            'method': 'pbe'
        },
        'ds_pbe_dcm': {
            'collation_function': g19_spe_coll_fn,
            'directory': 'xtb_spe_DFT',
            'solvent': 'dcm',
            'method': 'pbe'
        },
        'ds_pbe_cfm': {
            'collation_function': g19_spe_coll_fn,
            'directory': 'xtb_spe_DFT',
            'solvent': 'cfm',
            'method': 'pbe'
        },
        'ds_mp2_gas': {
            'collation_function': g19_spe_coll_fn,
            'directory': 'xtb_spe_DFT',
            'solvent': 'gas',
            'method': 'mp2'
        },
        'ds_mp2_dcm': {
            'collation_function': g19_spe_coll_fn,
            'directory': 'xtb_spe_DFT',
            'solvent': 'dcm',
            'method': 'mp2'
        },
        'ds_mp2_cfm': {
            'collation_function': g19_spe_coll_fn,
            'directory': 'xtb_spe_DFT',
            'solvent': 'cfm',
            'method': 'mp2'
        },
        'ds_O_b97_gas': {
            'collation_function': orca_spe_coll_fn,
            'directory': 'orca_calculations',
            'solvent': 'gas',
            'method': 'B97-3c'
        },
        'ds_O_mp2_gas': {
            'collation_function': orca_spe_coll_fn,
            'directory': 'orca_calculations',
            'solvent': 'gas',
            'method': 'MP2'
        },
    }

    all_energies = {}
    for method in _energy_methods:
        all_energies[method] = {}
        for struct in structure_titles:
            fn = _energy_methods[method]['collation_function']
            energy = fn(
                name=struct,
                directory=_energy_methods[method]['directory'],
                solvent=_energy_methods[method]['solvent'],
                method=_energy_methods[method]['method'],
            )
            all_energies[method][struct] = energy

    with open(output_ey_file, 'w') as f:
        # Write top line.
        top_line = 'name'
        for method in all_energies:
            top_line += f',{method}'
        top_line += '\n'
        f.write(top_line)
        # Rows.
        for struct in structure_titles:
            line = struct
            for method in all_energies:
                energy = all_energies[method][struct]
                line += f',{energy}'
            line += '\n'
            f.write(line)


if __name__ == '__main__':
    main()
