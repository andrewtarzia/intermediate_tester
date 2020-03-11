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


def flat_line(ax, x, y, w=0, C='k', m='x', label=None):
    ax.plot([x - w, x, x + w], [y, y, y], c=C, lw=2, label=label)
    ax.scatter(x, y, marker=m, c=C, edgecolor='none', s=80)


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
    amine=None
):

    fig, ax = plt.subplots(figsize=(8, 5))
    for ami in X:
        if amine is not None:
            if amine != ami:
                continue
        lab = leg_info[ami]['label']
        c = leg_info[ami]['c']
        for i in range(len(X[ami])):
            X_value = X[ami][i]
            Y_value = Y[ami][i] - min(Y[ami])
            print(X_value, Y_value)
            if i == 0:
                flat_line(
                    ax, x=X_value, y=Y_value,
                    w=1, C=c, label=lab
                )
            else:
                flat_line(ax, x=X_value, y=Y_value, w=1, C=c)
            if names[ami][i] in same_sizers:
                ax.text(
                    x=X_value-4,
                    y=Y_value-5,
                    s=names[ami][i]
                )

    ax.legend(fontsize=16)
    # ax.axhline(y=0, c='k', alpha=0.2, lw=2)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('cluster size', fontsize=16)
    ax.set_xlim(0, 30)
    # ax.set_ylim(-150, 150)
    ax.set_ylabel(ylbl, fontsize=16)
    ax.set_xticks(list(X_pos.values()))
    ax.set_xticklabels(list(X_pos.keys()))
    fig.tight_layout()
    fig.savefig(title, dpi=720, bbox_inches='tight')


def main():
    X_positions = {
        '[1+1]': 4*1,
        '[1+2]': 4*2,
        '[1+3]': 4*3,
        '[2+2]': 4*4,
        '[2+3]': 4*5,
        '[2+4]': 4*6,
        '[3+4]': 4*7,
        '[3+5]': 4*8,
        '[3+6]': 4*9,
        '[4+6]': 4*10
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

    for rxn in rxns:
        # if rxn['ami'] != 'ami1':
        #     continue
        # Get Free energies (*.fey).
        prod_energies = [read_ey(f'{i}.fey') for i in rxn['prod']]
        react_energies = [read_ey(f'{i}.fey') for i in rxn['react']]
        print(rxn)
        print(prod_energies)
        print(react_energies)
        # KJ/mol
        FE = sum(prod_energies) - sum(react_energies)
        print(FE)
        # KJ/mol per imine
        FE_pI = FE / rxn['no.imine']
        print(FE_pI)
        X_values[rxn['ami']].append(X_positions[rxn['size']])
        Y_values[rxn['ami']].append(FE)
        Y_values_pI[rxn['ami']].append(FE_pI)
        names[rxn['ami']].append(rxn['long-name'])

    plot_FE(
        X=X_values,
        Y=Y_values,
        leg_info=leg_info,
        names=names,
        same_sizers=same_sizers,
        ylbl='free energy of formation [kJ/mol]',
        X_pos=X_positions,
        title='fe_total.pdf',
        amine=None
    )
    plot_FE(
        X=X_values,
        Y=Y_values_pI,
        leg_info=leg_info,
        names=names,
        same_sizers=same_sizers,
        ylbl='free energy of formation\nper imine bond [kJ/mol]',
        X_pos=X_positions,
        title='fe_total.pdf',
        amine=None
    )

    for ami in X_values:
        plot_FE(
            X=X_values,
            Y=Y_values,
            leg_info=leg_info,
            names=names,
            same_sizers=same_sizers,
            ylbl='free energy of formation [kJ/mol]',
            X_pos=X_positions,
            title=f'fe_total_{ami}.pdf',
            amine=ami
        )
        plot_FE(
            X=X_values,
            Y=Y_values_pI,
            leg_info=leg_info,
            names=names,
            same_sizers=same_sizers,
            ylbl='free energy of formation\nper imine bond [kJ/mol]',
            X_pos=X_positions,
            title=f'fe_perimine_{ami}.pdf',
            amine=ami
        )


if __name__ == '__main__':
    main()
