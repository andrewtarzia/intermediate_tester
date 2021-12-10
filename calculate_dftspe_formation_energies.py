#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Calculate formation energies of intermediates from DFT energies.

Author: Andrew Tarzia

Date Created: 30 Mar 2021

"""

import sys
import pandas as pd

from reactions import landscape


def clean_extracted_dictionary(extracted_energy_file):
    extracted_ey = pd.read_csv(extracted_energy_file)
    extracted_ey['clean_name'] = [
        i.split('_opt')[0] for i in extracted_ey['name']
    ]
    extracted_ey['solvent'] = [
        i.split('_opt_')[1].split('_')[0]
        for i in extracted_ey['name']
    ]
    extracted_ey['method'] = [
        i.split('_opt_')[1].split('_')[1].replace('.gau', '')
        for i in extracted_ey['name']
    ]
    print(extracted_ey)
    return extracted_ey


def get_molecule_energy_from_frame(frame, molecule_name, key):
    t_frame = frame[frame['clean_name'] == molecule_name]
    if len(t_frame) != 1:
        raise ValueError(
            'Either no, or multiple instances of '
            f'{molecule_name} found in dataframe.'
        )

    return float(t_frame[key])


def main():
    if len(sys.argv) != 3:
        print(
            'csv file with DFT energies and prefix for output csv '
            'is required'
        )
        sys.exit()
    else:
        ey_file = sys.argv[1]
        output_prefix = sys.argv[2]

    _methods = ('pbe', 'mp2')
    _solvents = ('gas', 'dcm', 'cfm')
    _columns = (
        'intname', 'long-name', 'solvent', 'method', 'energy_au',
        'prodenergy_au', 'reactenergy_au', 'formenergy_au',
        'formenergypermine_au',
    )

    energy_frame = clean_extracted_dictionary(ey_file)

    lscp = landscape()

    for meth in _methods:
        m_frame = energy_frame[energy_frame['method'] == meth]
        for solv in _solvents:
            s_frame = m_frame[m_frame['solvent'] == solv]
            output_df = pd.DataFrame(columns=_columns)
            count = 1
            for intermediate in lscp:
                iname = intermediate['intname']
                lname = intermediate['long-name']

                # Get total energy.
                ey = get_molecule_energy_from_frame(
                    frame=s_frame,
                    molecule_name=iname,
                    key='energy',
                )
                # Sum product energies.
                prodey = sum(
                    [
                        get_molecule_energy_from_frame(
                            frame=s_frame,
                            # Modify name from xtb code.
                            molecule_name=i.replace('_opt_xtb', ''),
                            key='energy',
                        )
                        for i in intermediate['prod']
                    ]
                )
                # Sum reactant energies.
                reactey = sum(
                    [
                        get_molecule_energy_from_frame(
                            frame=s_frame,
                            # Modify name from xtb code.
                            molecule_name=i.replace('_opt_xtb', ''),
                            key='energy',
                        )
                        for i in intermediate['react']
                    ]
                )

                # Get reactant energies.
                formey = prodey - reactey
                formeyperimine = formey / intermediate['no.imine']
                if formey > 0:
                    print(intermediate, solv, formey, formeyperimine)

                # Add series to dataframe.
                output_df.loc[count] = [
                    iname, lname, solv, ey, prodey, reactey, formey,
                    formeyperimine,
                ]
                count += 1

            print(output_df)
            output_df.to_csv(
                (
                    f'{output_prefix}_{solv}_'
                    f'{meth}_formation_energies.csv'
                ),
                index=False,
            )


if __name__ == '__main__':
    main()
