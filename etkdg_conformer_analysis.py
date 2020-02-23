#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Perform ETKDG conformer analysis.

Author: Andrew Tarzia

Date Created: 23 Feb 2020

"""



def main():

    amine_files = [
        'ami1_opt_xtb.mol',
        'ami2_opt_xtb.mol',
        'ami3_opt_xtb.mol',
        'ami4_opt_xtb.mol',
    ]

    print(amine_files)
    results = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
    for ami in amine_files:
        amine = ami.replace('_opt_xtb.mol', '')
        ey_out_file = ami.replace('_opt_xtb.mol', 'conf_results.out')
        ami_dir = ami.replace('_opt_xtb.mol', '_confs')
