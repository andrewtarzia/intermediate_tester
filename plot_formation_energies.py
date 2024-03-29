#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Plot formation energies and pathways.

Author: Andrew Tarzia

Date Created: 30 Mar 2021

"""

import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

from reactions import landscape
from utilities import (
    flat_line, interplot_leg, interplot_ss, interplot_xpos
)


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
        print(f'doing plot_fe of {ami}')
        if amine is not None:
            if amine != ami:
                continue
        lab = leg_info[ami]['label']
        c = leg_info[ami]['c']
        for i in range(len(X[ami])):
            X_value = X[ami][i]
            Y_value = Y[ami][i] - min(Y[ami])
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
    ax.set_xlim(0, 45)
    ax.set_ylim(ylim)
    ax.set_ylabel(ylbl, fontsize=16)
    ax.set_xticks(list(X_pos.values()))
    ax.set_xticklabels(list(X_pos.keys()))
    fig.tight_layout()
    fig.savefig(title, dpi=720, bbox_inches='tight')
    plt.close()


def manu_fig2a(
    X,
    Y,
    leg_info,
    names,
    ylbl,
    X_pos,
    title,
    ylim,
    amine=None
):

    fig, ax = plt.subplots(figsize=(8, 5))
    for ami in X:
        count2 = 0
        print(f'doing manu_fig2a of {ami}')
        if amine is not None:
            if amine != ami:
                continue
        lab = leg_info[ami]['label']
        c = leg_info[ami]['c']
        for i in range(len(X[ami])):
            X_value = X[ami][i]
            Y_value = Y[ami][i] - min(Y[ami])
            if names[ami][i][-1] == '2':
                if count2 == 1:
                    flat_line(
                        ax,
                        x=X_value,
                        y=Y_value,
                        w=1,
                        C=c,
                        m='o',
                        label=f'{lab}'
                    )
                else:
                    flat_line(
                        ax,
                        x=X_value,
                        y=Y_value,
                        w=1,
                        m='o',
                        C=c
                    )
                count2 += 1

    ax.legend(fontsize=16)
    ax.axhline(y=0, c='k', alpha=0.2, lw=2)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('cluster size', fontsize=16)
    ax.set_xlim(0, 45)
    ax.set_ylim(ylim)
    ax.set_ylabel(ylbl, fontsize=16)
    ax.set_xticks(list(X_pos.values()))
    ax.set_xticklabels(list(X_pos.keys()))
    fig.tight_layout()
    fig.savefig(title, dpi=720, bbox_inches='tight')
    plt.close()


def manu_si_table(
    X,
    Y,
    names,
):

    print('======== for table ========')
    for ami in X:
        print(f'doing fe of {ami}')
        for i in range(len(X[ami])):
            X_value = X[ami][i]
            Y_value = Y[ami][i] - min(Y[ami])
            print(names[ami][i], X_value, Y_value)
    print('===========================')

def manu_fig2b_plot_FE_withaminal(
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

    fig, ax = plt.subplots(figsize=(8, 4))
    for ami in X:
        count1 = 0
        count2 = 0
        print(f'doing manu_fig2b_waminal of {ami}')
        if amine is not None:
            if amine != ami:
                continue
        lab = leg_info[ami]['label']
        c = leg_info[ami]['c']
        for i in range(len(X[ami])):
            X_value = X[ami][i]
            Y_value = Y[ami][i] - min(Y[ami])
            if names[ami][i][-1] == '2':
                if count2 == 1:
                    flat_line(
                        ax,
                        x=X_value,
                        y=Y_value,
                        w=1,
                        C=c,
                        m='o',
                        label=f'{lab}:imine'
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
                    )
                else:
                    flat_line(
                        ax,
                        x=X_value,
                        y=Y_value,
                        w=1,
                        C=c,

                    )
                count1 += 1

    ax.legend(fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('cluster size', fontsize=16)
    ax.set_xlim(0, 45)
    ax.set_ylim(ylim)
    ax.set_ylabel(ylbl, fontsize=16)
    ax.set_xticks(list(X_pos.values()))
    ax.set_xticklabels(list(X_pos.keys()))
    fig.tight_layout()
    fig.savefig(title, dpi=720, bbox_inches='tight')
    plt.close()


def main():
    if len(sys.argv) != 3:
        print('csv file with formation energies and output prefix')
        sys.exit()
    else:
        fe_csv = sys.argv[1]
        fe_pi_csv = sys.argv[2]

    figures_directory = 'figures'
    if not os.path.exists(figures_directory):
        os.mkdir(figures_directory)

    fe_data = pd.read_csv(fe_csv)
    fe_pi_data = pd.read_csv(fe_pi_csv)

    _energy_methods = set(
        i for i in list(fe_data.columns)
        if i != 'name'
    )
    # Get plot properties.
    xpos = interplot_xpos()
    leg_info = interplot_leg()
    same_sizers = interplot_ss()

    lscp = landscape()

    for method in _energy_methods:
        output_prefix = method
        names = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
        X_v = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
        Y_v = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
        Y_v_pI = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
        names_noam = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
        X_v_noam = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
        Y_v_noam = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
        Y_v_pI_noam = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}

        for intermediate in lscp:
            iname = intermediate['intname']
            spec_fe_df = fe_data[fe_data['name'] == iname]
            FE = float(spec_fe_df[method].iloc[0]) * 2625.5
            spec_fe_pi_df = fe_pi_data[fe_pi_data['name'] == iname]
            FE_pI = float(spec_fe_pi_df[method].iloc[0]) * 2625.5
            X_v[intermediate['ami']].append(xpos[intermediate['size']])
            Y_v[intermediate['ami']].append(FE)
            Y_v_pI[intermediate['ami']].append(FE_pI)
            names[intermediate['ami']].append(
                intermediate['long-name']
            )
            if intermediate['long-name'].split('_')[-1] != '1':
                X_v_noam[intermediate['ami']].append(
                    xpos[intermediate['size']]
                )
                Y_v_noam[intermediate['ami']].append(FE)
                Y_v_pI_noam[intermediate['ami']].append(FE_pI)
                names_noam[intermediate['ami']].append(
                    intermediate['long-name']
                )

        plot_FE(
            X=X_v,
            Y=Y_v,
            leg_info=leg_info,
            names=names,
            same_sizers=same_sizers,
            ylbl='formation energy [kJ mol$^{-1}$]',
            # ylbl='free energy of formation [kJ mol$^{-1}$]',
            X_pos=xpos,
            title=f'{figures_directory}/{output_prefix}_fe_total.pdf',
            ylim=(-10, 700),
            amine=None
        )
        plot_FE(
            X=X_v,
            Y=Y_v_pI,
            leg_info=leg_info,
            names=names,
            same_sizers=same_sizers,
            ylbl='formation energy\nper imine [kJ mol$^{-1}$]',
            X_pos=xpos,
            title=f'{figures_directory}/{output_prefix}_fe_perimine.pdf',
            ylim=(-10, 110),
            amine=None
        )
        manu_fig2a(
            X=X_v,
            Y=Y_v,
            leg_info=leg_info,
            names=names,
            ylbl='formation energy [kJ mol$^{-1}$]',
            X_pos=xpos,
            title=(
                f'{figures_directory}/fig2a_{output_prefix}_fe_total.pdf'
            ),
            ylim=(0, None),
            amine=None
        )
        if output_prefix == 'xtbgas':
            manu_si_table(
                X=X_v,
                Y=Y_v,
                names=names,
            )

        for ami in X_v:
            plot_FE(
                X=X_v,
                Y=Y_v,
                leg_info=leg_info,
                names=names,
                same_sizers=same_sizers,
                ylbl='formation energy [kJ mol$^{-1}$]',
                # ylbl='free energy of formation [kJ mol$^{-1}$]',
                X_pos=xpos,
                title=(
                    f'{figures_directory}/{output_prefix}_fe_total_{ami}'
                    '.pdf'
                ),
                ylim=(-10, 700),
                amine=ami,
            )
            plot_FE(
                X=X_v,
                Y=Y_v_pI,
                leg_info=leg_info,
                names=names,
                same_sizers=same_sizers,
                ylbl='formation energy\nper imine [kJ mol$^{-1}$]',
                # ylbl=(
                #   'free energy of formation\nper imine bond '
                #   '[kJ mol$^{-1}$]'
                # ),
                X_pos=xpos,
                title=(
                    f'{figures_directory}/{output_prefix}_fe_perimine_'
                    f'{ami}.pdf'
                ),
                ylim=(-10, 110),
                amine=ami
            )
            if output_prefix == 'ds_O_b97_gas':
                ylimmax = 180
            elif output_prefix == 'xtbgas':
                ylimmax = 600
            else:
                ylimmax = None
            manu_fig2b_plot_FE_withaminal(
                X=X_v,
                Y=Y_v,
                leg_info=leg_info,
                names=names,
                same_sizers=same_sizers,
                ylbl='formation energy [kJ mol$^{-1}$]',
                X_pos=xpos,
                title=(
                    f'{figures_directory}/manu_fig2b{output_prefix}'
                    f'_fe_total_{ami}.pdf'
                ),
                ylim=(-10, ylimmax),
                amine=ami
            )
            manu_fig2b_plot_FE_withaminal(
                X=X_v,
                Y=Y_v_pI,
                leg_info=leg_info,
                names=names,
                same_sizers=same_sizers,
                ylbl='formation energy\nper imine [kJ mol$^{-1}$]',
                # ylbl=(
                #   'free energy of formation\nper imine bond '
                #   '[kJ mol$^{-1}$]'
                # ),
                X_pos=xpos,
                title=(
                    f'{figures_directory}/manu_{output_prefix}_fe_'
                    f'perimine_{ami}.pdf'
                ),
                ylim=(0, 45),
                amine=ami
            )
            plot_FE(
                X=X_v_noam,
                Y=Y_v_noam,
                leg_info=leg_info,
                names=names_noam,
                same_sizers=same_sizers,
                ylbl='formation energy [kJ mol$^{-1}$]',
                # ylbl='free energy of formation [kJ mol$^{-1}$]',
                X_pos=xpos,
                title=(
                    f'{figures_directory}/{output_prefix}_fe_total_{ami}_'
                    'noaminal.pdf'
                ),
                ylim=(-10, 700),
                amine=ami
            )
            plot_FE(
                X=X_v_noam,
                Y=Y_v_pI_noam,
                leg_info=leg_info,
                names=names_noam,
                same_sizers=same_sizers,
                ylbl='formation energy\nper imine [kJ mol$^{-1}$]',
                # ylbl='free energy of formation\nper imine bond [kJ mol$^{-1}$]',
                X_pos=xpos,
                title=(
                    f'{figures_directory}/{output_prefix}_fe_perimine_'
                    f'{ami}_noaminal.pdf'
                ),
                ylim=(-10, 110),
                amine=ami
            )


if __name__ == '__main__':
    main()
