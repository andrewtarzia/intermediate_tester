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
    alde = stk.BuildingBlock.init_from_file('alde1.mol', ['aldehyde'])
    ami1 = stk.BuildingBlock.init_from_file('ami1.mol', ['amine'])
    cage1 = stk.ConstructedMolecule(
        building_blocks=[alde, ami1],
        topology_graph=stk.cage.FourPlusSix()
    )
    cage1.write('4_alde1_6_ami1_1.mol')
    ami2 = stk.BuildingBlock.init_from_file('ami2.mol', ['amine'])
    cage2 = stk.ConstructedMolecule(
        building_blocks=[alde, ami2],
        topology_graph=stk.cage.FourPlusSix()
    )
    cage2.write('4_alde1_6_ami2_1.mol')
    ami3 = stk.BuildingBlock.init_from_file('ami3.mol', ['amine'])
    cage3 = stk.ConstructedMolecule(
        building_blocks=[alde, ami3],
        topology_graph=stk.cage.FourPlusSix()
    )
    cage3.write('4_alde1_6_ami3_1.mol')
    ami4 = stk.BuildingBlock.init_from_file('ami4.mol', ['amine'])
    cage4 = stk.ConstructedMolecule(
        building_blocks=[alde, ami4],
        topology_graph=stk.cage.FourPlusSix()
    )
    cage4.write('4_alde1_6_ami4_1.mol')


if __name__ == '__main__':
    main()
