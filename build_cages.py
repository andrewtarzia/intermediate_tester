#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Build the complete cage structures.

Author: Andrew Tarzia

Date Created: 23 Feb 2020

"""

import stk


def main():
    alde = stk.BuildingBlock.init_from_file(
        'alde1.mol', [stk.AldehydeFactory()]
    )
    ami1 = stk.BuildingBlock.init_from_file(
        'ami1.mol', [stk.PrimaryAminoFactory()]
    )
    cage1 = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix(
            building_blocks=(alde, ami1),
        )
    )
    cage1.write('4_alde1_6_ami1_1.mol')
    ami2 = stk.BuildingBlock.init_from_file(
        'ami2.mol', [stk.PrimaryAminoFactory()]
    )
    cage2 = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix(
            building_blocks=(alde, ami2),
        )
    )
    cage2.write('4_alde1_6_ami2_1.mol')
    ami3 = stk.BuildingBlock.init_from_file(
        'ami3.mol', [stk.PrimaryAminoFactory()]
    )
    cage3 = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix(
            building_blocks=(alde, ami3),
        )
    )
    cage3.write('4_alde1_6_ami3_1.mol')
    ami4 = stk.BuildingBlock.init_from_file(
        'ami4.mol', [stk.PrimaryAminoFactory()]
    )
    cage4 = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix(
            building_blocks=(alde, ami4),
        )
    )
    cage4.write('4_alde1_6_ami4_1.mol')


if __name__ == '__main__':
    main()
