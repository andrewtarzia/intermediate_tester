#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Plot mass spectrum of three experiments.

Author: Andrew Tarzia

Date Created: 09 Feb 2022

"""

import matplotlib.pyplot as plt
import glob


def main():
    # To be run in mass_spec_image directory.
    files = glob.glob('*.txt')
    systems = {
        'DMHDA': '-',
        '+MPDA': '+CC13',
        '+monoamine': '+amine',
    }
    system_data = {}
    for sys in systems:
        sys_lbl = systems[sys]
        syss = []
        string = f'_{sys_lbl}_'
        for fi in files:
            if string in fi:
                syss.append(fi)
        system_data[sys] = syss

    y_pos = [2.2, 1.1, 0]
    _incs = {
        'Week1': 0,
        'Week2': 0.2,
        'Week3': 0.4,
        'Week4': 0.6,
        'Week6': 0.8,
    }
    _cs = {
        'Week1': '#648fff',
        'Week2': '#785ef0',
        'Week3': '#dc267f',
        'Week4': '#fe6100',
        'Week6': '#feb000',
    }

    fig, ax = plt.subplots(figsize=(8, 5))
    for i, sys in enumerate(systems):
        y_base = y_pos[i]
        print(y_base)
        sfiles = system_data[sys]
        print(sfiles)
        for j, fi in enumerate(sfiles):
            print(fi)
            week = fi.replace('.txt', '').replace('peaks', '')
            week = week.replace(f'_{systems[sys]}_', '')
            _inc = _incs[week]
            print(week, _inc)
            with open(fi, 'r') as f:
                lines = f.readlines()

            for li in lines:
                x = float(li.rstrip())
                print(li, x)
                xtoplot = [x, x]
                ytoplot = [y_base+_inc, y_base+_inc+0.16]
                ax.plot(
                    xtoplot,
                    ytoplot,
                    label=fi,
                    lw=2,
                    c=_cs[week],
                )

    peaks_picked = {
        # 305.2772,
        # 349.3114,
        '1+2': 415,
        # 460.4843,
        # 460.9885,
        # 476.5022,
        # 541.5114,
        '1+3': 540,
        # 573.5461,
        # 621.6327,
        # 649.6727,
        # 650.1745,
        '2+3': 666,
        # 919.8726,
        # 920.8757,
        # 921.8777,
        # 951.9074,
        '2+4': 793,
        '3+4': 919,
        '3+5': 1045,
        '3+6': 1189,
        '4+6': 1297,
    }
    for p in peaks_picked:
        x = peaks_picked[p]
        ax.axvline(x=x, lw=2, linestyle='--')
        ax.text(x=x, y=0.0, s=p)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('m/z', fontsize=16)
    ax.set_ylabel('count', fontsize=16)
    ax.set_xlim(300, 1400)
    ax.set_ylim(None, None)

    for y in y_pos:
        ax.axhline(y=y, c='k')
        for i in _incs:
            ax.axhline(y=y+_incs[i], c='k', alpha=0.2)

    fig.tight_layout()
    fig.savefig(
        'mass_specs.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()

if __name__ == '__main__':
    main()
