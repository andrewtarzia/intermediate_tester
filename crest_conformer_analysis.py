#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Perform ETKDG conformer analysis.

Author: Andrew Tarzia

Date Created: 23 Feb 2020

"""

import matplotlib.pyplot as plt
import pandas as pd


def main():

    ami_def = {
        'Distance 17-2': 'ami1',
        'Distance 2-14': 'ami2',
        'Distance 2-12': 'ami3',
        'Distance 2-7': 'ami4'
    }

    data = pd.read_csv('NN_dist.csv')
    res = {
        'ami1': {'c': 'k', 'd': [], 'e': []},
        'ami2': {'c': 'b', 'd': [], 'e': []},
        'ami3': {'c': 'r', 'd': [], 'e': []},
        'ami4': {'c': 'green', 'd': [], 'e': []}
    }
    for i, row in data.iterrows():
        energy = row[row.first_valid_index()]
        dist = row[row.last_valid_index()]
        ami = ami_def[row.last_valid_index()]
        # print(row, energy, dist, ami)
        res[ami]['e'].append(float(energy))
        res[ami]['d'].append(float(dist))

    fig, ax = plt.subplots()
    for ami in res:
        ax.scatter(
            res[ami]['d'],
            [
                i-min(res[ami]['e'])
                for i in res[ami]['e']
            ],
            c=res[ami]['c'],
            label=ami,
            alpha=0.6,
            edgecolor='none',
            marker='o',
            s=80
        )

    ax.legend(fontsize=16)
    ax.axhline(y=0, c='k', alpha=0.2, lw=2)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel(r'N-N distance [$\mathrm{\AA}$]', fontsize=16)
    # ax.set_xlim(0, 30)
    ax.set_ylim(-5, 70)
    ax.set_ylabel(
        'free energy [kJ/mol]',
        fontsize=16
    )
    fig.tight_layout()
    fig.savefig(
        f'crest_conf_analysis.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


if __name__ == '__main__':
    main()
