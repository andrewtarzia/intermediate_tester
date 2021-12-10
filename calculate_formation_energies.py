#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Calculate formation energies and pathways.

Author: Andrew Tarzia

Date Created: 23 Feb 2020

"""

import sys
import pandas as pd

from reactions import landscape


def get_molecule_energy_from_frame(frame, molecule_name, method):
    by_name = frame[frame['name'] == molecule_name]
    try:
        by_method = float(by_name[method].iloc[0])
    except ValueError:
        raise ValueError(f'Issue with {molecule_name} and {method}')
    return by_method


def main():
    if len(sys.argv) != 4:
        print(
            'csv file with all DFT energies and output csv is required'
        )
    else:
        ey_file = sys.argv[1]
        formation_energy_file = sys.argv[2]
        per_imine_file = sys.argv[3]

    all_structure_energies = pd.read_csv(ey_file)
    print(all_structure_energies)

    _energy_methods = set(
        list(all_structure_energies.columns)
    )
    print(_energy_methods)
    # Get the landscape of intermediates.
    lscp = landscape()

    all_formation_energies = {}
    all_formation_energies_per_imine = {}
    for method in _energy_methods:
        all_formation_energies[method] = {}
        all_formation_energies_per_imine[method] = {}
        for intermediate in lscp:
            iname = intermediate['intname']

            # Sum product energies.
            prodey = [
                get_molecule_energy_from_frame(
                    frame=all_structure_energies,
                    molecule_name=i,
                    method=method,
                )
                for i in intermediate['prod']
            ]

            # Sum reactant energies.
            reactey = [
                get_molecule_energy_from_frame(
                    frame=all_structure_energies,
                    molecule_name=i,
                    method=method,
                )
                for i in intermediate['react']
            ]

            formey = sum(prodey) - sum(reactey)
            formeyperimine = formey / intermediate['no.imine']
            if formey > 0:
                prods = intermediate['prod']
                react = intermediate['react']
                print(
                    f'formation energy ({formey}) is positive for '
                    f'{iname} using {method}.\n'
                    f'products: {prods} with ey: {prodey}\n'
                    f'ractants: {react} with ey: {reactey}\n\n'
                )
            all_formation_energies[method][iname] = formey
            all_formation_energies_per_imine[method][iname] = (
                formeyperimine
            )

    with open(formation_energy_file, 'w') as f:
        # Write top line.
        top_line = 'name'
        for method in all_formation_energies:
            top_line += f',{method}'
        top_line += '\n'
        f.write(top_line)
        # Rows.
        for intermediate in lscp:
            iname = intermediate['intname']
            line = iname
            for method in all_formation_energies:
                energy = all_formation_energies[method][iname]
                line += f',{energy}'
            line += '\n'
            f.write(line)

    with open(per_imine_file, 'w') as f:
        # Write top line.
        top_line = 'name'
        for method in all_formation_energies:
            top_line += f',{method}'
        top_line += '\n'
        f.write(top_line)
        # Rows.
        for intermediate in lscp:
            iname = intermediate['intname']
            line = iname
            for method in all_formation_energies_per_imine:
                energy = (
                    all_formation_energies_per_imine[method][iname]
                )
                line += f',{energy}'
            line += '\n'
            f.write(line)


if __name__ == '__main__':
    main()
