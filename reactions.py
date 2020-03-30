#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Define the reactions analysed.

Author: Andrew Tarzia

Date Created: 23 Feb 2020

"""


def reaction(lm, int_name, alde, ami, size, sts):

    dict = {
        'long-name': lm,
        'prod': [f'{int_name}_opt_xtb'] + ['water_xtb']*sts[2],
        'react': (
            [f'{alde}_opt_xtb']*sts[0] + [f'{ami}_opt_xtb']*sts[1]
        ),
        'no.imine': sts[2],
        'size': size,
        'alde': alde,
        'ami': ami
    }

    return dict


def reactions():
    """
    Define all reactions to get formation energy of.

    Products and reactants are given as the name of the .ey file.

    """

    reactions = []

    # Iterate over all amines and append reactions.
    for i in [1, 2, 3, 4]:

        reactions.append(reaction(
            f'Al-1_Am{i}-1_1',
            f'1_alde1_1_ami{i}_1',
            'alde1',
            f'ami{i}',
            '[1+1]',
            (1, 1, 1)
        ))
        reactions.append(reaction(
            f'Al-1_Am{i}-1_2',
            f'1_alde1_1_ami{i}_2',
            'alde1',
            f'ami{i}',
            '[1+1]',
            (1, 1, 1)
        ))
        reactions.append(reaction(
            f'Al-1_Am{i}-2_1',
            f'1_alde1_2_ami{i}_1',
            'alde1',
            f'ami{i}',
            '[1+2]',
            (1, 2, 2)
        ))
        reactions.append(reaction(
            f'Al-1_Am{i}-2_2',
            f'1_alde1_2_ami{i}_2',
            'alde1',
            f'ami{i}',
            '[1+2]',
            (1, 2, 2)
        ))
        reactions.append(reaction(
            f'Al-1_Am{i}-3_1',
            f'1_alde1_3_ami{i}_1',
            'alde1',
            f'ami{i}',
            '[1+3]',
            (1, 3, 3)
        ))
        reactions.append(reaction(
            f'Al-1_Am{i}-3_2',
            f'1_alde1_3_ami{i}_2',
            'alde1',
            f'ami{i}',
            '[1+3]',
            (1, 3, 3)
        ))
        reactions.append(reaction(
            f'Al-2_Am{i}-2_2',
            f'2_alde1_2_ami{i}_1',
            'alde1',
            f'ami{i}',
            '[2+2]',
            (2, 2, 4)
        ))
        reactions.append(reaction(
            f'Al-2_Am{i}-3_1',
            f'2_alde1_3_ami{i}_1',
            'alde1',
            f'ami{i}',
            '[2+3]',
            (2, 3, 5)
        ))
        reactions.append(reaction(
            f'Al-2_Am{i}-3_2',
            f'2_alde1_3_ami{i}_2',
            'alde1',
            f'ami{i}',
            '[2+3]',
            (2, 3, 5)
        ))
        reactions.append(reaction(
            f'Al-2_Am{i}-4_1',
            f'2_alde1_4_ami{i}_1',
            'alde1',
            f'ami{i}',
            '[2+4]',
            (2, 4, 6)
        ))
        reactions.append(reaction(
            f'Al-2_Am{i}-4_2',
            f'2_alde1_4_ami{i}_2',
            'alde1',
            f'ami{i}',
            '[2+4]',
            (2, 4, 6)
        ))
        reactions.append(reaction(
            f'Al-3_Am{i}-4_2',
            f'3_alde1_4_ami{i}_1',
            'alde1',
            f'ami{i}',
            '[3+4]',
            (3, 4, 8)
        ))
        reactions.append(reaction(
            f'Al-3_Am{i}-5_1',
            f'3_alde1_5_ami{i}_1',
            'alde1',
            f'ami{i}',
            '[3+5]',
            (3, 5, 9)
        ))
        reactions.append(reaction(
            f'Al-3_Am{i}-5_2',
            f'3_alde1_5_ami{i}_2',
            'alde1',
            f'ami{i}',
            '[3+5]',
            (3, 5, 9)
        ))
        reactions.append(reaction(
            f'Al-3_Am{i}-6_1',
            f'3_alde1_6_ami{i}_1',
            'alde1',
            f'ami{i}',
            '[3+6]',
            (3, 6, 9)
        ))
        reactions.append(reaction(
            f'Al-3_Am{i}-6_2',
            f'3_alde1_6_ami{i}_2',
            'alde1',
            f'ami{i}',
            '[3+6]',
            (3, 6, 9)
        ))
        reactions.append(reaction(
            f'Al-4_Am{i}-6_2',
            f'4_alde1_6_ami{i}_1',
            'alde1',
            f'ami{i}',
            '[4+6]',
            (4, 6, 12)
        ))

    return reactions
