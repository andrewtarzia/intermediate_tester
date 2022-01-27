#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module of utilities.

"""

from rdkit.Chem import AllChem as rdkit
from rdkit.Chem import Descriptors


def get_mass(name):
    mol = rdkit.MolFromMolFile(name)
    # Add Hs.
    mol = rdkit.AddHs(mol)
    MW = Descriptors.ExactMolWt(mol)
    return MW


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

    return ey


def interplot_xpos():
    distt = 5
    return {
        # '[1+1]': distt*1,
        '[1+2]': distt*1,
        '[1+3]': distt*2,
        # '[2+2]': distt*4,
        '[2+3]': distt*3,
        '[2+4]': distt*4,
        '[3+4]': distt*5,
        '[3+5]': distt*6,
        '[3+6]': distt*7,
        '[4+6]': distt*8,
    }


def interplot_ss():
    return [
        'Al-2_Am1-4_1',
        'Al-2_Am1-4_2',
        'Al-2_Am2-4_1',
        'Al-2_Am2-4_2',
        'Al-2_Am3-4_1',
        'Al-2_Am3-4_2',
        'Al-2_Am4-4_1',
        'Al-2_Am4-4_2',
    ]


def interplot_leg():
    return {
        'ami1': {
            'label': 'amine-1',
            'c': '#D81B60'
        },
        'ami2': {
            'label': 'amine-2',
            'c': '#1E88E5'
        },
        'ami3': {
            'label': 'amine-3 (CC3)',
            'c': '#FFC107'
        },
        'ami4': {
            'label': 'amine-4 (CC1)',
            'c': '#004D40'
        }
    }
