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
import json
import numpy as np
import glob
import matplotlib.pyplot as plt


def calculate_NN_distance(stk_mol):
    # Assumes that there are only two N atoms in the molecules
    N_atom_ids = [i.id for i in stk_mol.atoms if i.atomic_number == 7]
    NN_dist = stk_mol.get_atom_distance(N_atom_ids[0], N_atom_ids[1])
    print(NN_dist)
    return NN_dist


def get_dihedral(pt1, pt2, pt3, pt4):
    """
    Calculate the dihedral between four points.

    Uses Praxeolitic formula --> 1 sqrt, 1 cross product

    Output in range (-pi to pi).

    From: https://stackoverflow.com/questions/20305272/
    dihedral-torsion-angle-from-four-points-in-cartesian-
    coordinates-in-python
    (new_dihedral(p))

    """
    p0 = np.asarray(pt1)
    p1 = np.asarray(pt2)
    p2 = np.asarray(pt3)
    p3 = np.asarray(pt4)

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def calculate_NCCN_dihedral(stk_mol):
    # Assumes that there are only two N atoms in the molecules
    N_atom_ids = [i.id for i in stk_mol.atoms if i.atomic_number == 7]
    print(N_atom_ids)
    # Get all bonds from N atoms.
    all_N_bonds = []
    for bond in stk_mol.bonds:
        if bond.atom1.id in N_atom_ids or bond.atom2.id in N_atom_ids:
            all_N_bonds.append(bond)

    # Find N bonds to C atoms, that are also connected.
    connected_C_atom_ids = [None, None]
    for b in all_N_bonds:
        non_N = b.atom1.id if b.atom2.id in N_atom_ids else b.atom2.id
        N_id = b.atom2.id if b.atom2.id in N_atom_ids else b.atom1.id
        for b2 in all_N_bonds:
            if connected_C_atom_ids[0] is not None:
                break
            if b == b2:
                continue
            non_N2 = (
                b2.atom1.id
                if b2.atom2.id in N_atom_ids else b2.atom2.id
            )
            N2_id = (
                b2.atom2.id
                if b2.atom2.id in N_atom_ids else b2.atom1.id
            )
            # Find a bond between non_N and non_N2 - thats the dihedral
            # Else, move on.
            for bond in stk_mol.bonds:
                if connected_C_atom_ids[0] is not None:
                    break
                cond1 = bond.atom1.id in [non_N, non_N2]
                cond2 = bond.atom2.id in [non_N, non_N2]
                if cond1 and cond2:
                    print('found!', bond, b, b2)
                    # preserve order
                    if N_id == N_atom_ids[0]:
                        connected_C_atom_ids[0] = non_N
                        connected_C_atom_ids[1] = non_N2
                    elif N2_id == N_atom_ids[0]:
                        connected_C_atom_ids[0] = non_N2
                        connected_C_atom_ids[1] = non_N

    print(N_atom_ids, connected_C_atom_ids)
    pts = [
        i for i in stk_mol.get_atom_positions(
            atom_ids=N_atom_ids+connected_C_atom_ids
        )
    ]
    NCCN_dihed = get_dihedral(
        pt1=pts[0],
        pt2=pts[2],
        pt3=pts[3],
        pt4=pts[1],
    )

    print(NCCN_dihed)
    return NCCN_dihed


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
        'ami4': 'NCCN',
        'ami1b': 'CC(C)[C@H](N)[C@@H](/N=C/c1ccccc1)C(C)C',
        'ami2b': 'CC(C)(N)C/N=C/c1ccccc1',
        'ami3b': 'N[C@@H]1CCCC[C@H]1/N=C/c1ccccc1',
        'ami4b': 'NCC/N=C\\c1ccccc1'
    }

    print(amines)
    results = {
        'ami1': [], 'ami2': [], 'ami3': [], 'ami4': [],
        'ami1b': [], 'ami2b': [], 'ami3b': [], 'ami4b': []
    }
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
            NCCN_dihed = calculate_NCCN_dihedral(conf)
            opt_conf = stk.BuildingBlock.init_from_file(
                opt_conf_file,
                use_cache=False
            )
            opt_NN_dist = calculate_NN_distance(opt_conf)
            opt_NCCN_dihed = calculate_NCCN_dihedral(opt_conf)
            res[conf_id] = {
                'f_energy': f_energy,
                'opt_f_energy': opt_f_energy,
                'opt_energy': opt_energy,
                'NN_dis': NN_dist,
                'opt_NN_dis': opt_NN_dist,
                'NCCN_dihed': NCCN_dihed,
                'opt_NCCN_dihed': opt_NCCN_dihed,
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
        },
        'ami1b': {
            'label': 'amine-1 + benzaldehyde',
            'c': '#FFEAE5'
        },
        'ami2b': {
            'label': 'amine-2 + benzaldehyde',
            'c': '#FFF8DF'
        },
        'ami3b': {
            'label': 'amine-3 + benzaldehyde',
            'c': '#E8F8F5'
        },
        'ami4b': {
            'label': 'amine-4 + benzaldehyde',
            'c': '#EAF2F8'
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
    figh1d, axshd = plt.subplots(4, 1, figsize=(8, 10), sharex=True)
    # Remove horizontal space between axes
    figh1d.subplots_adjust(hspace=0)
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
    xlimd = (-5, 185)
    ylim = (-3, 22)
    binwidth_x = 0.1
    binwidth_xd = 5
    binwidth_y = 1
    for i, ami in enumerate(results):
        if i > 3:
            continue
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
        X2 = [
            float(results[ami][i]['opt_NN_dis']) for i in results[ami]
        ]
        NC2 = [
            abs(float(results[ami][i]['opt_NCCN_dihed']))
            for i in results[ami]
        ]

        Y1 = [float(results[ami][i]['f_energy']) for i in results[ami]]
        Y1 = [2625.5*(i-min(Y1)) for i in Y1]
        Y2 = [
            float(results[ami][i]['opt_energy'])
            for i in results[ami]
        ]
        Y2 = [2625.5*(i-min(Y2)) for i in Y2]
        ami_b = f'{ami}b'
        cb = leg_info[ami_b]['c']
        labb = leg_info[ami_b]['label']
        Xb2 = [
            float(results[ami_b][i]['opt_NN_dis'])
            for i in results[ami_b]
        ]
        NCb2 = [
            abs(float(results[ami_b][i]['opt_NCCN_dihed']))
            for i in results[ami_b]
        ]
        Yb1 = [
            float(results[ami_b][i]['f_energy'])
            for i in results[ami_b]
        ]
        Yb1 = [2625.5*(i-min(Yb1)) for i in Yb1]
        Yb2 = [
            float(results[ami_b][i]['opt_energy'])
            for i in results[ami_b]
        ]
        Yb2 = [2625.5*(i-min(Yb2)) for i in Yb2]

        bins_x = np.arange(xlim[0], xlim[1] + binwidth_x, binwidth_x)
        bins_xd = np.arange(
            xlimd[0],
            xlimd[1] + binwidth_xd,
            binwidth_xd
        )
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
            'energy [kJmol$^{-1}$]',
            fontsize=16
        )
        fig.tight_layout()
        fig.savefig(
            f'{ami}_etkdg_conf_analysis.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()

        # Only want N-N distances for conformers within 10 kJmol$^{-1}$.
        final_xs = []
        final_dihed_xs = []
        for x, y, xd in zip(X2, Y2, NC2):
            if y < 10:
                final_xs.append(x)
                final_dihed_xs.append(xd)
        final_xbs = []
        final_dihed_xbs = []
        for x, y, xd in zip(Xb2, Yb2, NCb2):
            if y < 10:
                final_xbs.append(x)
                final_dihed_xbs.append(xd)

        print(ami, len(final_xs), len(X2))
        print(ami_b, len(final_xbs), len(Xb2))

        axsh[i].hist(
            final_xs,
            bins=bins_x,
            alpha=0.6,
            edgecolor='k',
            linewidth=1.2,
            facecolor=c,
            density=True,
            label=lab
        )
        axsh[i].hist(
            final_xbs,
            bins=bins_x,
            alpha=0.7,
            edgecolor='k',
            linewidth=1.2,
            facecolor=cb,
            density=True,
            label=labb
        )
        axsh[i].tick_params(axis='both', which='major', labelsize=16)
        axsh[i].set_xlim(xlim)
        axsh[i].set_ylabel('frequency', fontsize=16)

        axshd[i].hist(
            final_dihed_xs,
            bins=bins_xd,
            alpha=0.6,
            edgecolor='k',
            linewidth=1.2,
            facecolor=c,
            density=True,
            label=lab
        )
        axshd[i].hist(
            final_dihed_xbs,
            bins=bins_xd,
            alpha=0.7,
            edgecolor='k',
            linewidth=1.2,
            facecolor=cb,
            density=True,
            label=labb
        )
        axshd[i].tick_params(axis='both', which='major', labelsize=16)
        axshd[i].set_xlim(xlimd)
        axshd[i].set_ylabel('frequency', fontsize=16)

        axsh2[i].hist(
            Y2,
            bins=bins_y,
            alpha=0.6,
            edgecolor='k',
            linewidth=1.2,
            facecolor=c,
            density=True,
            label=lab
        )
        axsh2[i].hist(
            Yb2,
            bins=bins_y,
            alpha=0.7,
            edgecolor='k',
            linewidth=1.2,
            facecolor=cb,
            density=True,
            label=labb
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
        'energy [kJmol$^{-1}$]',
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
    axsh2[3].set_xlabel('energy [kJmol$^{-1}$]', fontsize=16)
    figh2.tight_layout()
    figh2.savefig(
        f'etkdg_conf_E_analysis_hist.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()

    figh1d.legend(fontsize=16, loc=1)
    axshd[3].set_xlabel(r'NCCN dihedral [$^{\circ}$]', fontsize=16)
    figh1d.tight_layout()
    figh1d.savefig(
        f'etkdg_conf_DI_analysis_hist.pdf',
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


if __name__ == '__main__':
    main()
