#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Perform ETKDG conformer analysis.

Author: Andrew Tarzia

Date Created: 23 Feb 2020

"""

from rdkit.Chem import AllChem as rdkit
import stk
from os import mkdir
from os.path import exists, join
from itertools import product
import json
import numpy as np
import glob
import matplotlib.pyplot as plt


def calculate_NN_distance(stk_mol):
    # Assumes that there are only two N atoms in the molecules
    for atom1, atom2 in product(stk_mol.atoms, repeat=2):
        if atom1.id == atom2.id:
            continue
        if atom1.atomic_number == 7 and atom2.atomic_number == 7:
            print(atom1, atom2)
            NN_dist = stk_mol.get_atom_distance(atom1.id, atom2.id)
            print(NN_dist)
            return NN_dist


def optimise_conformer(name, mol, solvent, output_file):
    print(f'doing xtb opt of {name}')
    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent

    print('doing XTB optimisation')
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'{name}_opt',
        gfn_version=2,
        num_cores=6,
        opt_level='extreme',
        charge=0,
        num_unpaired_electrons=0,
        max_runs=5,
        electronic_temperature=300,
        calculate_hessian=True,
        unlimited_memory=True,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )
    xtb_opt.optimize(mol=mol)
    mol.write(output_file)

    return mol


def build_conformers(mol, dir):
    # Build conformers with ETKDG.
    etkdg = rdkit.ETKDG()
    etkdg.randomSeed = 1000
    cids = rdkit.EmbedMultipleConfs(
        mol=mol,
        numConfs=500,
        params=etkdg
    )

    # Save each conformer into directory.
    for cid in cids:
        filename = join(dir, f'conformer_{cid}.mol')
        rdkit.MolToMolFile(mol, filename=filename, confId=cid)


def load_f_energy(file):
    with open(file, 'r') as f:
        for line in f.readlines():
            return line.rstrip()


def calculate_f_energy(name, mol, solvent, ey_file, fey_file):
    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent
    print(f'getting energy of {name}')
    xtb_energy = stk.XTBEnergy(
        xtb_path='/home/atarzia/software/xtb-190806/bin/xtb',
        output_dir=f'{name}_ey',
        num_cores=6,
        charge=0,
        num_unpaired_electrons=0,
        electronic_temperature=300,
        unlimited_memory=True,
        calculate_free_energy=True,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )
    energy = xtb_energy.get_energy(mol)
    free_energy = xtb_energy.total_free_energies[mol]
    with open(ey_file, 'w') as f:
        f.write(str(energy))
    with open(fey_file, 'w') as f:
        f.write(str(free_energy))


def main():

    amines = {
        'ami1': 'CC(C)[C@H](N)[C@@H](N)C(C)C',
        'ami2': 'CC(C)(N)CN',
        'ami3': 'N[C@@H]1CCCC[C@H]1N',
        'ami4': 'NCCN'
    }

    print(amines)
    results = {'ami1': [], 'ami2': [], 'ami3': [], 'ami4': []}
    for amine in amines:
        res = {}
        ey_out_file = f'{amine}_conf_results.out'
        ami_dir = f'{amine}_confs'

        # Do calcs, only if the two energy files do not already exist.
        if exists(ey_out_file):
            with open(ey_out_file, 'r') as f:
                lst = json.load(f)
                results[amine] = lst
            continue

        # Read structure into rdkit.
        rdk_mol = rdkit.MolFromSmiles(amines[amine])
        rdk_mol = rdkit.AddHs(rdk_mol)
        print(f'doing {amine}')
        # Write directory.
        if not exists(ami_dir):
            mkdir(ami_dir)

        # Get N conformers.
        print(f'building conformers of {amine}')
        build_conformers(mol=rdk_mol, dir=ami_dir)

        # For each conformer, iterate.
        print(f'getting energies of conformers of {amine}')
        for conf_file in glob.glob(join(ami_dir, 'conformer*mol')):
            if 'opt' in conf_file:
                continue
            print(conf_file)
            opt_conf_file = conf_file.replace('.mol', '_opt.mol')
            ey_file = conf_file.replace('.mol', '.ey')
            fey_file = conf_file.replace('.mol', '.fey')
            opt_ey_file = conf_file.replace('.mol', '_opt.ey')
            opt_fey_file = conf_file.replace('.mol', '_opt.fey')
            conf_id = conf_file.replace(ami_dir+'/', '')
            conf_id = conf_id.split('_')[1].replace('.mol', '')
            print(ey_file, conf_id)

            # Load in stk object.
            conf = stk.BuildingBlock.init_from_file(conf_file)
            print(conf)

            # Get xTB energy, optimise with xTB, and get energy again.
            # Read in opt and unopt energies.
            if not exists(ey_file) or not exists(fey_file):
                # Calculate energy.
                calculate_f_energy(
                    name=join(ami_dir, f'unoptey_{conf_id}'),
                    mol=conf,
                    solvent=('chcl3', 'verytight'),
                    ey_file=ey_file,
                    fey_file=fey_file
                )
            # Load energy.
            f_energy = load_f_energy(fey_file)

            if not exists(opt_ey_file) or not exists(opt_fey_file):
                if not exists(opt_conf_file):
                    # Optimise conformer.
                    optimise_conformer(
                        name=join(ami_dir, f'opt_{conf_id}'),
                        mol=conf,
                        solvent=('chcl3', 'verytight'),
                        output_file=opt_conf_file
                    )
                opt_conf = stk.BuildingBlock.init_from_file(
                    opt_conf_file
                )
                # Calculate energy.
                calculate_f_energy(
                    name=join(ami_dir, f'optey_{conf_id}'),
                    mol=opt_conf,
                    solvent=('chcl3', 'verytight'),
                    ey_file=opt_ey_file,
                    fey_file=opt_fey_file
                )
            # Load energy.
            opt_energy = load_f_energy(opt_ey_file)
            opt_f_energy = load_f_energy(opt_fey_file)
            print(fey_file, opt_fey_file)
            print(f_energy)
            print(opt_f_energy)

            # Calculate N-N distance.
            conf = stk.BuildingBlock.init_from_file(
                conf_file,
                use_cache=False
            )
            NN_dist = calculate_NN_distance(conf)
            opt_conf = stk.BuildingBlock.init_from_file(
                opt_conf_file,
                use_cache=False
            )
            opt_NN_dist = calculate_NN_distance(opt_conf)
            res[conf_id] = {
                'f_energy': f_energy,
                'opt_f_energy': opt_f_energy,
                'opt_energy': opt_energy,
                'NN_dis': NN_dist,
                'opt_NN_dis': opt_NN_dist,
            }

        # Save results to JSON file and results dict.
        results[amine] = res
        with open(ey_out_file, 'w') as f:
            json.dump(res, f)

        print('.................done.......................')

    # Plot.
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

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    figh1, axsh = plt.subplots(4, 1, figsize=(8, 10), sharex=True)
    # Remove horizontal space between axes
    figh1.subplots_adjust(hspace=0)
    figh2, axsh2 = plt.subplots(4, 1, figsize=(8, 10), sharex=True)
    # Remove horizontal space between axes
    figh2.subplots_adjust(hspace=0)

    fig1 = plt.figure(figsize=(8, 8))
    ax_s1 = plt.axes(rect_scatter)
    ax_s1.tick_params(
        axis='both', which='major', labelsize=16,
        direction='in', top=True, right=True)
    ax_h1x = plt.axes(rect_histx)
    ax_h1x.tick_params(
        axis='both', which='major', labelsize=16,
        direction='in', labelbottom=False
    )
    ax_h1y = plt.axes(rect_histy)
    ax_h1y.tick_params(
        axis='both', which='major', labelsize=16,
        direction='in', labelleft=False
    )
    xlim = (2.2, 4.2)
    ylim = (-3, 22)
    binwidth_x = 0.1
    binwidth_y = 1
    for i, ami in enumerate(results):
        fig = plt.figure(figsize=(8, 8))
        ax_s = plt.axes(rect_scatter)
        ax_s.tick_params(
            axis='both', which='major', labelsize=16,
            direction='in', top=True, right=True
        )
        ax_hx = plt.axes(rect_histx)
        ax_hx.tick_params(
            axis='both', which='major', labelsize=16,
            direction='in', labelbottom=False
        )
        ax_hy = plt.axes(rect_histy)
        ax_hy.tick_params(
            axis='both', which='major', labelsize=16,
            direction='in', labelleft=False
        )
        lab = leg_info[ami]['label']
        c = leg_info[ami]['c']
        X1 = [float(results[ami][i]['NN_dis']) for i in results[ami]]
        X2 = [
            float(results[ami][i]['opt_NN_dis']) for i in results[ami]
        ]
        Y1 = [float(results[ami][i]['f_energy']) for i in results[ami]]
        Y1 = [2625.5*(i-min(Y1)) for i in Y1]
        Y2 = [
            float(results[ami][i]['opt_f_energy'])
            for i in results[ami]
        ]
        Y2 = [2625.5*(i-min(Y2)) for i in Y2]

        bins_x = np.arange(xlim[0], xlim[1] + binwidth_x, binwidth_x)
        bins_y = np.arange(ylim[0], ylim[1] + binwidth_y, binwidth_y)
        ax_hx.hist(
            X2,
            bins=bins_x,
            alpha=0.6,
            edgecolor='k',
            facecolor=c
        )
        ax_hy.hist(
            Y2,
            bins=bins_y,
            orientation='horizontal',
            alpha=0.6,
            edgecolor='k',
            facecolor=c
        )

        ax_hx.set_xlim(xlim)
        ax_hy.set_ylim(ylim)

        ax_h1x.hist(
            X2,
            bins=bins_x,
            alpha=0.6,
            edgecolor='k',
            facecolor=c
        )
        ax_h1y.hist(
            Y2,
            bins=bins_y,
            orientation='horizontal',
            alpha=0.6,
            edgecolor='k',
            facecolor=c
        )

        ax_h1x.set_xlim(xlim)
        ax_h1y.set_ylim(ylim)

        # Unoptimised.
        # ax.scatter(
        #     X1,
        #     Y1,
        #     c=c,
        #     alpha=1,
        #     edgecolor='none',
        #     marker='o',
        #     s=80,
        #     label=lab
        # )
        # Optimised.
        ax_s.scatter(
            X2,
            Y2,
            c=c,
            alpha=0.6,
            edgecolor='none',
            marker='o',
            s=80,
            label=lab
        )
        ax_s1.scatter(
            X2,
            Y2,
            c=c,
            alpha=0.6,
            edgecolor='none',
            marker='o',
            s=80,
            label=lab
        )

        ax_s.legend(fontsize=16)
        ax_s.axhline(y=0, c='k', alpha=0.2, lw=2)
        ax_s.tick_params()
        ax_s.set_xlabel(r'N-N distance [$\mathrm{\AA}$]', fontsize=16)
        ax_s.set_xlim(xlim)
        ax_s.set_ylim(ylim)
        ax_s.set_ylabel(
            'free energy [kJ/mol]',
            fontsize=16
        )
        fig.tight_layout()
        fig.savefig(
            f'{ami}_etkdg_conf_analysis.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()

        # Only want N-N distances for conformers within 20 kJ/mol.
        final_xs = []
        for x, y in zip(X2, Y2):
            if y < 10:
                final_xs.append(x)

        print(ami, len(final_xs), len(X2))

        axsh[i].hist(
            final_xs,
            bins=bins_x,
            alpha=0.6,
            edgecolor='k',
            facecolor=c,
            density=True,
            label=lab
        )
        axsh[i].tick_params(axis='both', which='major', labelsize=16)
        axsh[i].set_xlim(xlim)
        axsh[i].set_ylabel('frequency', fontsize=16)

        axsh2[i].hist(
            Y2,
            bins=bins_y,
            alpha=0.6,
            edgecolor='k',
            facecolor=c,
            density=True,
            label=lab
        )
        axsh2[i].tick_params(axis='both', which='major', labelsize=16)
        axsh2[i].set_ylabel('frequency', fontsize=16)

    ax_s1.legend(fontsize=16)
    ax_s1.axhline(y=0, c='k', alpha=0.2, lw=2)
    ax_s1.tick_params(axis='both', which='major', labelsize=16)
    ax_s1.set_xlabel(r'N-N distance [$\mathrm{\AA}$]', fontsize=16)
    ax_s1.set_xlim(xlim)
    ax_s1.set_ylim(ylim)
    ax_s1.set_ylabel(
        'free energy [kJ/mol]',
        fontsize=16
    )
    fig1.tight_layout()
    fig1.savefig(
        f'etkdg_conf_analysis.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()

    figh1.legend(fontsize=16)
    axsh[3].set_xlabel(r'N-N distance [$\mathrm{\AA}$]', fontsize=16)
    figh1.tight_layout()
    figh1.savefig(
        f'etkdg_conf_DD_analysis_hist.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()

    figh2.legend(fontsize=16)
    axsh2[3].set_xlabel('free energy [kJ/mol]', fontsize=16)
    figh2.tight_layout()
    figh2.savefig(
        f'etkdg_conf_E_analysis_hist.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


if __name__ == '__main__':
    main()
