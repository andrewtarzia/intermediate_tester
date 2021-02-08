#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Calculate and plot formation energies and pathways.

Author: Andrew Tarzia

Date Created: 23 Feb 2020

"""

import matplotlib.pyplot as plt
from reactions import reactions
from rdkit.Chem import AllChem as rdkit
from rdkit.Chem import Descriptors


def flat_line(ax, x, y, w=0, C='k', m='x', label=None):
    ax.plot([x - w, x, x + w], [y, y, y], c=C, lw=2)
    ax.scatter(
        x,
        y,
        marker=m,
        c=C,
        edgecolor='none',
        s=80,
        label=label
    )


def read_ey(file):

    with open(file, 'r') as f:
        lines = f.readlines()
        ey = float(lines[0].rstrip())

    return ey*2625.5


def plot_FE(
    X,
    Y,
    leg_info,
    names,
    same_sizers,
    ylbl,
    X_pos,
    title,
    ylim,
    amine=None
):

    fig, ax = plt.subplots(figsize=(8, 5))
    for ami in X:
        count1 = 0
        count2 = 0
        print(f'doing fe of {ami}')
        if amine is not None:
            if amine != ami:
                continue
        lab = leg_info[ami]['label']
        c = leg_info[ami]['c']
        for i in range(len(X[ami])):
            X_value = X[ami][i]
            Y_value = Y[ami][i] - min(Y[ami])
            print(title, names[ami][i], X_value, Y_value)
            if names[ami][i][-1] == '2':
                if count2 == 1:
                    flat_line(
                        ax,
                        x=X_value,
                        y=Y_value,
                        w=1,
                        C=c,
                        label=f'{lab}:imine'
                    )
                else:
                    flat_line(
                        ax,
                        x=X_value,
                        y=Y_value,
                        w=1,
                        C=c
                    )
                count2 += 1
            else:
                if count1 == 0:
                    flat_line(
                        ax,
                        x=X_value,
                        y=Y_value,
                        w=1,
                        C=c,
                        label=f'{lab}:aminal',
                        m='o'
                    )
                else:
                    flat_line(
                        ax,
                        x=X_value,
                        y=Y_value,
                        w=1,
                        C=c,
                        m='o'
                    )
                count1 += 1

    ax.legend(fontsize=16)
    ax.axhline(y=0, c='k', alpha=0.2, lw=2)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('cluster size', fontsize=16)
    ax.set_xlim(0, 55)
    ax.set_ylim(ylim)
    ax.set_ylabel(ylbl, fontsize=16)
    ax.set_xticks(list(X_pos.values()))
    ax.set_xticklabels(list(X_pos.keys()))
    fig.tight_layout()
    fig.savefig(title, dpi=720, bbox_inches='tight')


def get_mass(name):
    mol = rdkit.MolFromMolFile(name)
    # Add Hs.
    mol = rdkit.AddHs(mol)
    MW = Descriptors.ExactMolWt(mol)
    return MW


def main():
    distt = 5
    X_positions = {
        '[1+1]': distt*1,
        '[1+2]': distt*2,
        '[1+3]': distt*3,
        '[2+2]': distt*4,
        '[2+3]': distt*5,
        '[2+4]': distt*6,
        '[3+4]': distt*7,
        '[3+5]': distt*8,
        '[3+6]': distt*9,
        '[4+6]': distt*10
    }
    rxns = reactions()

    same_sizers = [
        'Al-2_Am1-4_1',
        'Al-2_Am1-4_2',
        'Al-2_Am2-4_1',
        'Al-2_Am2-4_2',
        'Al-2_Am3-4_1',
        'Al-2_Am3-4_2',
        'Al-2_Am4-4_1',
        'Al-2_Am4-4_2',
    ]
    leg_info = {
        'ami1': {
            'label': 'amine-1',
            'c': '#FF5733'
        },
        'ami2': {
            'label': 'amine-2',
            'c': '#FFC300'
        },
        'ami3': {
            'label': 'amine-3 (CC3)',
            'c': '#48C9B0'
        },
        'ami4': {
            'label': 'amine-4 (CC1)',
            'c': '#5499C7'
        }
    }

    X_values = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
    Y_values = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
    names = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
    Y_values_pI = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
    X_values_noaminal = {
        'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []
    }
    Y_values_noaminal = {
        'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []
    }
    Y_values_pI_noaminal = {
        'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []
    }
    names_noaminal = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}

    for rxn in rxns:
        # if rxn['ami'] != 'ami1':
        #     continue
        # Get Free energies (*.fey).
        prod_energies = [read_ey(f'{i}.ey') for i in rxn['prod']]
        react_energies = [read_ey(f'{i}.ey') for i in rxn['react']]
        # kJmol$^{-1}$
        FE = sum(prod_energies) - sum(react_energies)
        # kJmol$^{-1}$ per imine
        FE_pI = FE / rxn['no.imine']
        if FE > 0:
            input()
        X_values[rxn['ami']].append(X_positions[rxn['size']])
        Y_values[rxn['ami']].append(FE)
        Y_values_pI[rxn['ami']].append(FE_pI)
        names[rxn['ami']].append(rxn['long-name'])
        if rxn['long-name'].split('_')[-1] != '1':
            X_values_noaminal[rxn['ami']].append(
                X_positions[rxn['size']]
            )
            Y_values_noaminal[rxn['ami']].append(FE)
            Y_values_pI_noaminal[rxn['ami']].append(FE_pI)
            names_noaminal[rxn['ami']].append(rxn['long-name'])

        cage_file = f"{rxn['prod'][0]}.mol"
        mass = get_mass(cage_file)
        if rxn['long-name'].split('_')[-1] != '1':
            is_aminal = False
        else:
            is_aminal = True
        amine_energy = read_ey(f"{rxn['ami']}_opt_xtb.ey")
        aldehyde_energy = read_ey(f"{rxn['alde']}_opt_xtb.ey")
        water_energy = read_ey('water_opt_xtb.ey')
        cage_energy = read_ey(f"{rxn['prod'][0]}.ey")
        print(
            f"{rxn['long-name']},{rxn['size']},{mass},{is_aminal},"
            f"{cage_energy},"
            f"{rxn['ami']},{rxn['no_ami']},{amine_energy},"
            f"{rxn['alde']},{rxn['no_alde']},{aldehyde_energy},"
            f"{rxn['no.imine']},{rxn['no.imine']},{water_energy},"
            f"{sum(prod_energies)},{sum(react_energies)},"
            f"{FE},{FE_pI}"
        )

    plot_FE(
        X=X_values,
        Y=Y_values,
        leg_info=leg_info,
        names=names,
        same_sizers=same_sizers,
        ylbl='energy of formation [kJmol$^{-1}$]',
        # ylbl='free energy of formation [kJmol$^{-1}$]',
        X_pos=X_positions,
        title='fe_total.pdf',
        ylim=(-10, 700),
        amine=None
    )
    plot_FE(
        X=X_values,
        Y=Y_values_pI,
        leg_info=leg_info,
        names=names,
        same_sizers=same_sizers,
        ylbl='energy of formation\nper imine bond [kJmol$^{-1}$]',
        # ylbl='free energy of formation\nper imine bond [kJmol$^{-1}$]',
        X_pos=X_positions,
        title='fe_perimine.pdf',
        ylim=(-10, 110),
        amine=None
    )

    for ami in X_values:
        plot_FE(
            X=X_values,
            Y=Y_values,
            leg_info=leg_info,
            names=names,
            same_sizers=same_sizers,
            ylbl='energy of formation [kJmol$^{-1}$]',
            # ylbl='free energy of formation [kJmol$^{-1}$]',
            X_pos=X_positions,
            title=f'fe_total_{ami}.pdf',
            ylim=(-10, 700),
            amine=ami
        )
        plot_FE(
            X=X_values,
            Y=Y_values_pI,
            leg_info=leg_info,
            names=names,
            same_sizers=same_sizers,
            ylbl='energy of formation\nper imine bond [kJmol$^{-1}$]',
            # ylbl='free energy of formation\nper imine bond [kJmol$^{-1}$]',
            X_pos=X_positions,
            title=f'fe_perimine_{ami}.pdf',
            ylim=(-10, 110),
            amine=ami
        )
        plot_FE(
            X=X_values_noaminal,
            Y=Y_values_noaminal,
            leg_info=leg_info,
            names=names_noaminal,
            same_sizers=same_sizers,
            ylbl='energy of formation [kJmol$^{-1}$]',
            # ylbl='free energy of formation [kJmol$^{-1}$]',
            X_pos=X_positions,
            title=f'fe_total_{ami}_noaminal.pdf',
            ylim=(-10, 700),
            amine=ami
        )
        plot_FE(
            X=X_values_noaminal,
            Y=Y_values_pI_noaminal,
            leg_info=leg_info,
            names=names_noaminal,
            same_sizers=same_sizers,
            ylbl='energy of formation\nper imine bond [kJmol$^{-1}$]',
            # ylbl='free energy of formation\nper imine bond [kJmol$^{-1}$]',
            X_pos=X_positions,
            title=f'fe_perimine_{ami}_noaminal.pdf',
            ylim=(-10, 110),
            amine=ami
        )


if __name__ == '__main__':
    main()
