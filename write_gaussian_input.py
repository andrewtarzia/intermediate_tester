#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Write input files for Gaussian single point energy calculation.

Author: Andrew Tarzia

Date Created: 30 Mar 2021

"""

import glob
import stk


def write_top_section(name, np):

    new_name = name.split('.')[0]

    string = (
        f'%chk={new_name}.chk\n'
        f'%mem=60000mb\n'
        f'%nprocshared={np}\n\n'
    )

    return string


def write_defn_section(
    name,
    runtype,
    solvent=None,
    method='PBE1PBE',
    org_basis=None,
):

    if solvent is None:
        solvent_s = ''
    else:
        solvent_s = solvent

    new_name = name.split('.')[0]

    string = (
        f'#P {method}/{org_basis} '
        f'{runtype} '
        'SCF=(YQC,MaxCycle=900) '
        'int=(Grid=Superfinegrid) '
        'EmpiricalDispersion=GD3 '
        f'{solvent_s}\n\n'
        f'{new_name}\n\n'
    )

    return string


def write_molecule_section(struct, restart=False):

    charge = 0
    multiplicity = 1

    if restart:
        string = f'{charge} {multiplicity}\n\n'
    else:
        pos_mat = struct.get_position_matrix()
        string = f'{charge} {multiplicity}\n'
        for atom in struct.get_atoms():
            E = atom.__class__.__name__
            x, y, z = pos_mat[atom.get_id()]
            string += (
                f'{E} {round(x, 4)} {round(y, 4)} {round(z, 4)}\n'
            )
        string += '\n'

    return string


def run_file_top_line(name, np):

    string = (
        f'#PBS -N _{name}\n'
        f'#PBS -l select=1:ncpus={np}:mem=124gb\n'
        '#PBS -l walltime=72:00:00\n'
        'module load gaussian/g16-c01-avx\n\n'
        'cd $PBS_O_WORKDIR\n\n'
    )
    return string


def write_input_file(
    infile,
    struct,
    np,
    directory,
    runtype,
    method,
    org_basis,
    solvent,
):
    base_name = '_'+infile.replace('.in', '')

    string = write_top_section(
        name=base_name,
        np=np,
    )
    string += write_defn_section(
        name=base_name,
        runtype=runtype,
        method=method,
        solvent=solvent,
        org_basis=org_basis,
    )
    string += write_molecule_section(struct, restart=False)

    with open(f'{directory}/{infile}', 'w') as f:
        f.write(string)


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def main():

    num_proc = 32
    structures = glob.glob('*_opt.mol')
    directory = 'xtb_spe_DFT'

    _solvents = {
        'gas': None,
        'dcm': r'SCRF=(PCM,Solvent=Dichloromethane)',
        'cfm': r'SCRF=(PCM,Solvent=Chloroform)',
        # r'SCRF=(PCM,Solvent=DiMethylSulfoxide)'
    }

    run_lines = []
    for solv in _solvents:
        for s in structures:
            mol = stk.BuildingBlock.init_from_file(s)
            gau_inp = s.replace('.mol', f'_{solv}.gau')

            write_input_file(
                infile=gau_inp,
                struct=mol,
                np=num_proc,
                directory=directory,
                method='PBE1PBE',
                runtype='SP',
                org_basis='Def2TZVP',
                solvent=_solvents[solv],
            )
            run_line = (
                f"g16 < {gau_inp} > "
                f"{gau_inp.replace('.gau', '.log')}\n"
            )
            run_lines.append(run_line)

    # To add to a multiple run .sh file.
    for i, rls in enumerate(chunks(run_lines, 20)):
        with open(f'{directory}/spes_{i}.sh', 'w') as f:
            f.write(run_file_top_line(name=f's_{i}', np=num_proc))
            f.write('\n')
            for rl in rls:
                f.write(rl)


if __name__ == '__main__':
    main()
