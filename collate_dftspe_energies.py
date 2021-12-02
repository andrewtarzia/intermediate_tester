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


def get_energy(output_file, directory, method):

    with open(f'{directory}/{output_file}', 'r') as f:
        for line in f.readlines():
            if 'A.U. after' in line:
                l = line.rstrip()
                if method == 'pbe':
                    energy = l.split('SCF Done:  E(RPBE1PBE) =')
                elif method == 'mp2':
                    energy = l.split('SCF Done:  E(RHF) =')
                energy = energy[-1].split('A.U. after')[0]
                return float(energy)  # a.u.


def main():
    if len(sys.argv) != 2:
        print('csv file to write DFT energies to.')
        sys.exit()
    else:
        output_ey_file = sys.argv[1]

    structures = sorted(glob.glob('*_opt.mol'))
    directory = 'xtb_spe_DFT'

    _solvents = {
        'gas': None,
        'dcm': r'SCRF=(PCM,Solvent=Dichloromethane)',
        'cfm': r'SCRF=(PCM,Solvent=Chloroform)',
        # r'SCRF=(PCM,Solvent=DiMethylSulfoxide)'
    }

    _methods = {
        'pbe': ('PBE1PBE', 'Def2TZVP', 'EmpiricalDispersion=GD3'),
        'mp2': ('MP2', 'aug-cc-pVDZ', None),
    }

    all_energies = {}
    for meth in _methods:
        for solv in _solvents:
            for s in structures:
                name = s.replace('.mol', f'_{solv}_{meth}')
                output_file = (
                    s.replace('.mol', f'_{solv}_{meth}.log')
                )

                energy = get_energy(output_file, directory, meth)
                all_energies[name] = energy

    with open(output_ey_file, 'w') as f:
        f.write('name,energy\n')
        for e in all_energies:
            f.write(f'{e},{all_energies[e]}\n')


if __name__ == '__main__':
    main()
