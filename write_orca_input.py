#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Write input files for ORCA calculations.

Author: Andrew Tarzia

Date Created: 25 Oct 2021

"""

import os
import glob
import stk


def write_top_section(basename, method, solvent):

    if solvent is None:
        solvent_s = ''
    else:
        solvent_s = solvent

    if method == 'PBE0':
        string = (
            f'! DFT COPT RKS PBE0 def2-SVP D4 Grid6 NOFINALGRID '
            f'SlowConv TightSCF printbasis {solvent_s}'
            '\n\n'
            f'%base "{basename}"\n'
            '%maxcore 3000\n'
            '%scf\n   MaxIter 2000\nend\n'
        )
    elif method == 'B97-3c':
        string = (
            f'! DFT COPT B97-3c TightSCF printbasis {solvent_s} '
            'Grid6 NOFINALGRID SlowConv '
            '\n\n'
            f'%base "{basename}"\n'
            '%maxcore 3000\n'
            '%scf\n   MaxIter 2000\nend\n'
        )

    return string


def write_proc_section(np):

    return f'%pal\n   nprocs {np}\nend\n\n'


def write_molecule_section(prefix, method):

    charge = 0
    multi = 1

    return (
        f'* xyzfile {charge} {multi} o_{prefix}.xyz\n'
        f'# * xyzfile {charge} {multi} _o_{prefix}_{method}.xyz\n'
    )


def write_mp2_file(prefix, np, directory):

    basename = f'_o_{prefix}_MP2'
    infile = f'o_{prefix}_MP2.in'
    charge = 0
    multi = 1

    string = (
        '! RI-MP2 cc-pVTZ cc-pVTZ/C '
        'TIGHTSCF printbasis Grid6 NOFINALGRID SlowConv\n\n'
        f'%base "{basename}"\n'
        '%maxcore 3000\n'
        f'%pal\n   nprocs {np}\nend\n\n'
        f'* xyzfile {charge} {multi} _o_{prefix}_B97-3c.xyz\n'
    )

    with open(f'{directory}/{infile}', 'w') as f:
        f.write(string)


def write_input_file(
    prefix,
    np,
    directory,
    method,
    solvent,
):

    basename = f'_o_{prefix}_{method}'
    infile = f'o_{prefix}_{method}.in'

    string = write_top_section(
        basename=basename,
        method=method,
        solvent=solvent,
    )
    string += write_proc_section(np)
    string += write_molecule_section(prefix, method)

    with open(f'{directory}/{infile}', 'w') as f:
        f.write(string)


def main():

    num_proc = 32
    structures = glob.glob('*_opt.mol')
    directory = 'orca_calculations'
    if not os.path.exists(directory):
        os.mkdir(directory)

    _methods = {'b973c': 'B97-3c'}

    for s in structures:
        mol = stk.BuildingBlock.init_from_file(s)
        prefix = s.replace('.mol', '')
        mol.write(f'{directory}/o_{prefix}.xyz')

        write_input_file(
            prefix=prefix,
            np=num_proc,
            directory=directory,
            method=_methods['b973c'],
            solvent=None,
        )

        write_mp2_file(
            prefix=prefix,
            np=num_proc,
            directory=directory,
        )


if __name__ == '__main__':
    main()
