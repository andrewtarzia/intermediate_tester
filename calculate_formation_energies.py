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


def main():
    if len(sys.argv) != 4:
        print(
            'csv file with DFT energies and output csv is required'
            ' - solv is also req: gas or chcl3.'
        )
    else:
        ey_file = sys.argv[1]
        output_name = sys.argv[2]
        solv = sys.argv[3]

    print(ey_file, output_name)
    extracted_ey = pd.read_csv(ey_file)
    extracted_ey['clean_name'] = [i for i in extracted_ey['names']]
    energy_dict = {
        i: j
        for i, j in zip(
            extracted_ey['clean_name'],
            extracted_ey[f'energy_{solv}'],
        )
    }
    print(energy_dict)
    lscp = landscape()
    df = pd.DataFrame({
        'intname': [i['intname'] for i in lscp],
        'long-name': [i['long-name'] for i in lscp],
    })

    energies = []
    prod_energies = []
    react_energies = []
    formation_energies = []
    formation_energies_per_imine = []
    for intermediate in lscp:
        iname = intermediate['intname']
        if intermediate['intname'] in energy_dict:
            # Get energies.
            ey = energy_dict[iname]
            prodey = sum(
                [
                    energy_dict[i.replace('_opt_xtb', '')]
                    for i in intermediate['prod']
                ]
            )
            reactey = sum(
                [
                    energy_dict[i.replace('_opt_xtb', '')]
                    for i in intermediate['react']
                ]
            )
            formey = prodey - reactey
            formeyperimine = formey / intermediate['no.imine']
            if formey > 0:
                input('he')

            energies.append(ey)
            prod_energies.append(prodey)
            react_energies.append(reactey)
            formation_energies.append(formey)
            formation_energies_per_imine.append(formeyperimine)
        else:
            raise ValueError()

    df['energy_au'] = energies
    df['prodenergy_au'] = prod_energies
    df['reactenergy_au'] = react_energies
    df['formenergy_au'] = formation_energies
    df['formenergypermine_au'] = formation_energies_per_imine
    df['formenergy_kjmol'] = [
        None if i is None else i*2625.5 for i in formation_energies
    ]
    df['formenergypermine_kjmol'] = [
        None if i is None else i*2625.5 for i in formation_energies_per_imine
    ]
    print(df)
    df.to_csv(output_name, index=False)


if __name__ == '__main__':
    main()
