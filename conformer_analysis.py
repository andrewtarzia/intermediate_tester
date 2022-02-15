#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Perform conformer analysis of small molecules.

Author: Andrew Tarzia

Date Created: 07 Dec 2020

"""

from rdkit.Chem import AllChem as rdkit
import sys
import stk
# from scipy.stats import gaussian_kde
import stko
from os import mkdir
from os.path import exists, join
import json
import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean


def result_key_defn():

    return {
        'f_energy': {
            'label': 'relative energy [kJmol$^{-1}$]',
            'lim': (0, 100),
            'width': 5,
        },
        'NN_dis': {
            'label': r'N-N distance [$\mathrm{\AA}$]',
            'lim': (2.5, 4),
            'width': 0.05,
        },
        'NCCN_dihed': {
            'label': r'NCCN dihedral [$^{\circ}$]',
            'lim': (0, 200),
            'width': 5,
        },
    }


def leg_info():
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
        },
        'ami1b': {
            'label': 'amine-1 + benzaldehyde',
            'c': '#33DBFF'
        },
        'ami2b': {
            'label': 'amine-2 + benzaldehyde',
            'c': '#003CFF'
        },
        'ami3b': {
            'label': 'amine-3 + benzaldehyde',
            'c': '#C94861'
        },
        'ami4b': {
            'label': 'amine-4 + benzaldehyde',
            'c': '#C78254'
        }
    }


def calculate_NN_distance(stk_mol):
    position_matrix = stk_mol.get_position_matrix()
    # Assumes that there are only two N atoms in the molecules
    N_atom_ids = [
        i.get_id() for i in stk_mol.get_atoms()
        if i.get_atomic_number() == 7
    ]
    distance = euclidean(
        u=position_matrix[N_atom_ids[0]],
        v=position_matrix[N_atom_ids[1]],
    )
    print(distance)
    return float(distance)


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
    smarts = '[#7]-[#6]-[#6]-[#7]'
    # Find torsions.
    rdkit_mol = stk_mol.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)
    # Calculate torsional angle for all imines.
    torsion_list = []
    for atom_ids in rdkit_mol.GetSubstructMatches(
        rdkit.MolFromSmarts(smarts),
    ):
        torsion = get_dihedral(
            pt1=tuple(
                stk_mol.get_atomic_positions(atom_ids[0])
            )[0],
            pt2=tuple(
                stk_mol.get_atomic_positions(atom_ids[1])
            )[0],
            pt3=tuple(
                stk_mol.get_atomic_positions(atom_ids[2])
            )[0],
            pt4=tuple(
                stk_mol.get_atomic_positions(atom_ids[3])
            )[0]
        )
        torsion_list.append(abs(torsion))
    if len(torsion_list) != 1:
        raise ValueError(
            f'something wrong with the torsions in {stk_mol}'
        )

    return torsion_list[0]


def build_conformers(mol, dir, solvent=None):
    nconfs = 500
    # Build conformers with ETKDG.
    etkdg = rdkit.ETKDGv3()
    etkdg.randomSeed = 1000
    etkdg.numThreads = 4
    etkdg.useRandomCoords = True
    # etkdg.pruneRmsThresh = 0.05
    print('building conformers')
    cids = rdkit.EmbedMultipleConfs(
        mol=mol,
        numConfs=nconfs,
        params=etkdg
    )
    print(
        f'extracted {len(cids)} conformers out of {nconfs}'
    )

    # Save each conformer into directory.
    stk_conformers = {}
    for cid in cids:
        filename = join(dir, f'conformer_{cid}.mol')
        rdkit.MolToMolFile(mol, filename=filename, confId=cid)

        print(f'doing xtb opt of {cid}')
        if solvent is None:
            solvent_str = None
            solvent_grid = 'normal'
        else:
            solvent_str, solvent_grid = solvent

        conf = stk.BuildingBlock.init_from_file(filename)
        xtb_opt = stko.XTB(
            xtb_path='/home/atarzia/software/xtb-6.3.2/bin/xtb',
            output_dir=f'{dir}/{cid}_opt',
            gfn_version=2,
            num_cores=6,
            opt_level='normal',
            charge=0,
            num_unpaired_electrons=0,
            max_runs=5,
            electronic_temperature=300,
            calculate_hessian=True,
            unlimited_memory=True,
            solvent=solvent_str,
            solvent_grid=solvent_grid
        )
        conf = xtb_opt.optimize(mol=conf)
        conf.write(join(dir, f'conformer_{cid}.mol'))
        stk_conformers[cid] = conf
    return stk_conformers


def load_opls_conformers(
    amine,
    temp_mol,
    structure_dir,
    outputdir,
    solvent=None
):
    conformer_xyzs = sorted(glob.glob(
        f'{structure_dir}/{amine}*c.xyz'
    ))
    print(
        f'extracting {len(conformer_xyzs)} conformers'
    )

    # Save each conformer into directory.
    stk_conformers = {}
    for confxyz in conformer_xyzs:
        cid = int(
            confxyz.split('_')[1].replace('c.xyz', '')
        )
        filename = join(outputdir, f'conformer_{cid}.mol')

        print(f'doing xtb opt of {cid}')
        if solvent is None:
            solvent_str = None
            solvent_grid = 'normal'
        else:
            solvent_str, solvent_grid = solvent

        conf = temp_mol.with_structure_from_file(confxyz)

        xtb_opt = stko.XTB(
            xtb_path='/home/atarzia/software/xtb-6.3.2/bin/xtb',
            output_dir=f'{outputdir}/{cid}_opt',
            gfn_version=2,
            num_cores=6,
            opt_level='normal',
            charge=0,
            num_unpaired_electrons=0,
            max_runs=5,
            electronic_temperature=300,
            calculate_hessian=True,
            unlimited_memory=True,
            solvent=solvent_str,
            solvent_grid=solvent_grid
        )
        conf = xtb_opt.optimize(mol=conf)
        conf.write(filename)
        stk_conformers[cid] = conf
    return stk_conformers


def load_f_energy(file):
    with open(file, 'r') as f:
        for line in f.readlines():
            return float(line.rstrip())


def calculate_f_energy(name, mol, solvent, ey_file, fey_file):

    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent

    print(f'getting energy of {name}')
    xtb_energy = stko.XTBEnergy(
        xtb_path='/home/atarzia/software/xtb-6.3.2/bin/xtb',
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


def run_etkdg_workflow(amines):

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
                results[amine] = json.load(f)
            save_lowest_energy_structures(
                results=results[amine],
                name=amine,
                dir=ami_dir,
            )
            continue

        # Read structure into rdkit.
        rdk_mol = rdkit.MolFromSmiles(amines[amine])
        rdk_mol = rdkit.AddHs(rdk_mol)
        print(f'doing {amine}')
        # Write directory.
        if not exists(ami_dir):
            mkdir(ami_dir)

        # Get N conformers.
        print(f'building and optimising conformers of {amine}')
        stk_conformers = build_conformers(mol=rdk_mol, dir=ami_dir)

        # For each conformer, iterate.
        print(f'getting properties of conformers of {amine}')
        for cid in stk_conformers:
            conf = stk_conformers[cid]
            ey_file = join(ami_dir, f'conformer_{cid}.ey')
            fey_file = join(ami_dir, f'conformer_{cid}.fey')
            print(cid, ey_file, fey_file)

            # Get xTB energy, optimise with xTB, and get energy again.
            # Read in opt and unopt energies.
            if not exists(ey_file) or not exists(fey_file):
                # Calculate energy.
                calculate_f_energy(
                    name=join(ami_dir, f'{cid}'),
                    mol=conf,
                    solvent=None,
                    ey_file=ey_file,
                    fey_file=fey_file
                )
            # Load energy.
            f_energy = load_f_energy(fey_file)
            print(fey_file, f_energy)

            # Calculate N-N distance and NCCN dihedral.
            NN_dist = calculate_NN_distance(conf)
            if NN_dist > 3:
                print(cid, NN_dist, amine)
            NCCN_dihed = calculate_NCCN_dihedral(conf)
            res[cid] = {
                'f_energy': f_energy,
                'NN_dis': NN_dist,
                'NCCN_dihed': NCCN_dihed,
            }

        # Save results to JSON file and results dict.
        results[amine] = res
        with open(ey_out_file, 'w') as f:
            json.dump(res, f)
        print('.................done.......................')

    return results


def run_opls3e_workflow(amines, structure_dir):

    results = {
        'ami1': [], 'ami2': [], 'ami3': [], 'ami4': [],
        'ami1b': [], 'ami2b': [], 'ami3b': [], 'ami4b': []
    }
    for amine in amines:
        if 'b' == amine[-1]:
            continue
        res = {}
        ey_out_file = f'{amine}_opls3e_conf_results.out'
        ami_dir = f'{amine}_opls3e_confs'

        # Do calcs, only if the two energy files do not already exist.
        if exists(ey_out_file):
            with open(ey_out_file, 'r') as f:
                results[amine] = json.load(f)
            continue

        # Read a temparary structure into stk.
        temp_mol = stk.BuildingBlock.init_from_file(
            f'{structure_dir}/{amine}.mol'
        )
        print(f'doing {amine}')
        # Write directory.
        if not exists(ami_dir):
            mkdir(ami_dir)

        # Get N conformers.
        print(f'building and optimising conformers of {amine}')
        stk_conformers = load_opls_conformers(
            amine=amine,
            temp_mol=temp_mol,
            structure_dir=structure_dir,
            outputdir=ami_dir,
        )

        # For each conformer, iterate.
        print(f'getting properties of conformers of {amine}')
        for cid in stk_conformers:
            conf = stk_conformers[cid]
            ey_file = join(ami_dir, f'conformer_{cid}.ey')
            fey_file = join(ami_dir, f'conformer_{cid}.fey')
            print(cid, ey_file, fey_file)

            # Get xTB energy, optimise with xTB, and get energy again.
            # Read in opt and unopt energies.
            if not exists(ey_file) or not exists(fey_file):
                # Calculate energy.
                calculate_f_energy(
                    name=join(ami_dir, f'{cid}'),
                    mol=conf,
                    solvent=None,
                    ey_file=ey_file,
                    fey_file=fey_file
                )
            # Load energy.
            f_energy = load_f_energy(fey_file)
            print(fey_file, f_energy)

            # Calculate N-N distance and NCCN dihedral.
            NN_dist = calculate_NN_distance(conf)
            NCCN_dihed = calculate_NCCN_dihedral(conf)
            res[cid] = {
                'f_energy': f_energy,
                'NN_dis': NN_dist,
                'NCCN_dihed': NCCN_dihed,
            }

        # Save results to JSON file and results dict.
        results[amine] = res
        with open(ey_out_file, 'w') as f:
            json.dump(res, f)
        print('.................done.......................')

    return results


def density_of_all_amines(results, filename, ignore_b):

    # For each property.
    result_keys = result_key_defn()
    for rkey in result_keys:
        rkeyinfo = result_keys[rkey]
        if rkey == 'NN_dis':
            # bottoms = [12, 8, 4, 0]
            density = True
            ylab = 'frequency'
            # ylim = 24
        elif rkey == 'NCCN_dihed':
            # bottoms = [0.6, 0.5, 0.15, 0]
            density = True
            ylab = 'frequency'
            # ylim = 0.7
        else:
            # bottoms = [0, 0, 0, 0]
            density = True
            ylab = 'frequency'
            # ylim = None

        fig, ax = plt.subplots(8, 1, sharex=True, figsize=(8, 10))
        # Remove horizontal space between axes
        fig.subplots_adjust(hspace=0)

        # For each amine.
        amine_names = [
            'ami1', 'ami1b', 'ami2', 'ami2b',
            'ami3', 'ami3b', 'ami4', 'ami4b',
        ]
        for i, ami in enumerate(amine_names):
            if ignore_b:
                if ami[-1] == 'b':
                    continue
            result_dict = results[ami]

            # Get data.
            data = {
                cid: (
                    result_dict[cid]['f_energy'],
                    result_dict[cid][rkey]
                )
                for cid in result_dict
            }

            eydata = {cid: data[cid][0] for cid in data}
            # Set to relative energies.
            eydata = {
                cid: eydata[cid]-min(eydata.values())
                for cid in eydata
            }
            # Set to kJ/mol.
            eydata = {cid: eydata[cid]*2625.5 for cid in eydata}

            if rkey == 'f_energy':
                xdata = [eydata[cid] for cid in eydata]
            else:
                xdata = [
                    data[cid][1] for cid in data
                ]
            print(ami, rkey)
            print(min(xdata), max(xdata))

            # Plot data.
            # density = gaussian_kde(
            #     xdata,
            #     bw_method=0.3,
            # )
            xwidth = rkeyinfo['width']
            # xs = np.linspace(
            #     rkeyinfo['lim'][0] - xwidth,
            #     rkeyinfo['lim'][1] + xwidth,
            #     1000,
            # )
            # ax.plot(
            #     xs, density(xs),
            #     linewidth=2.,
            #     color=leg_info()[ami]['c'],
            #     alpha=1.0,
            #     label=leg_info()[ami]['label'],
            # )

            xbins = np.arange(
                rkeyinfo['lim'][0] - xwidth,
                rkeyinfo['lim'][1] + xwidth,
                xwidth
            )
            ax[i].hist(
                x=xdata,
                bins=xbins,
                density=density,
                # bottom=bottoms[i],
                histtype='stepfilled',
                # histtype='step',  fill=False,
                stacked=True,
                linewidth=1.,
                facecolor=leg_info()[ami]['c'],
                alpha=1.0,
                # color=leg_info()[ami]['c'],
                # color='white',
                edgecolor='k',
                label=leg_info()[ami]['label'],
            )
            ax[i].set_yticks([])
            # ax[i].set_ylim(-1, 1)

            ax[i].set_xlim(rkeyinfo['lim'])
            ax[i].tick_params(left=False, bottom=False)
            # ax.set_ylim(0, ylim)
        ax[7].tick_params(labelsize=16, bottom=True)
        ax[7].set_xlabel(rkeyinfo['label'], fontsize=16)
        ax[7].set_ylabel(ylab, fontsize=16)
        fig.legend(fontsize=16, ncol=2)
        fig.tight_layout()
        fig.savefig(
            f'{filename}_{rkey}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


def manuscript_density_of_all_amines(results, filename):

    # For each property.
    result_keys = result_key_defn()
    for rkey in result_keys:
        rkeyinfo = result_keys[rkey]
        if rkey == 'NN_dis':
            # bottoms = [12, 8, 4, 0]
            density = True
            ylab = 'frequency'
            # ylim = 24
        elif rkey == 'NCCN_dihed':
            # bottoms = [0.6, 0.5, 0.15, 0]
            density = True
            ylab = 'frequency'
            # ylim = 0.7
        else:
            continue

        fig, ax = plt.subplots(4, 1, sharex=True, figsize=(8, 5))
        # Remove horizontal space between axes
        fig.subplots_adjust(hspace=0)

        # For each amine.
        for i, ami in enumerate(results):
            if ami[-1] == 'b':
                continue
            result_dict = results[ami]

            # Get data.
            data = {
                cid: (
                    result_dict[cid]['f_energy'],
                    result_dict[cid][rkey]
                )
                for cid in result_dict
            }

            eydata = {cid: data[cid][0] for cid in data}
            # Set to relative energies.
            eydata = {
                cid: eydata[cid]-min(eydata.values())
                for cid in eydata
            }
            # Set to kJ/mol.
            eydata = {cid: eydata[cid]*2625.5 for cid in eydata}

            if rkey == 'f_energy':
                xdata = [eydata[cid] for cid in eydata]
            else:
                # Only want lowest 10 kJ/mol.
                xdata = [
                    data[cid][1] for cid in data
                    if eydata[cid] < 10
                ]

            print(min(xdata), max(xdata))

            # Plot data.
            # density = gaussian_kde(
            #     xdata,
            #     bw_method=0.3,
            # )
            xwidth = rkeyinfo['width']
            # xs = np.linspace(
            #     rkeyinfo['lim'][0] - xwidth,
            #     rkeyinfo['lim'][1] + xwidth,
            #     1000,
            # )
            # ax.plot(
            #     xs, density(xs),
            #     linewidth=2.,
            #     color=leg_info()[ami]['c'],
            #     alpha=1.0,
            #     label=leg_info()[ami]['label'],
            # )

            xbins = np.arange(
                rkeyinfo['lim'][0] - xwidth,
                rkeyinfo['lim'][1] + xwidth,
                xwidth
            )
            ax[i].hist(
                x=xdata,
                bins=xbins,
                density=density,
                # bottom=bottoms[i],
                histtype='stepfilled',
                # histtype='step',  fill=False,
                stacked=True,
                linewidth=1.,
                facecolor=leg_info()[ami]['c'],
                alpha=1.0,
                # color=leg_info()[ami]['c'],
                # color='white',
                edgecolor='k',
                label=leg_info()[ami]['label'],
            )
            ax[i].set_yticks([])
            # ax[i].set_ylim(-1, 1)

            ax[i].set_xlim(rkeyinfo['lim'])
            ax[i].tick_params(left=False, bottom=False)
            # ax.set_ylim(0, ylim)
        ax[3].tick_params(labelsize=16, bottom=True)
        ax[3].set_xlabel(rkeyinfo['label'], fontsize=16)
        ax[3].set_ylabel(ylab, fontsize=16)
        fig.legend(fontsize=16, ncol=2)
        fig.tight_layout()
        fig.savefig(
            f'{filename}_{rkey}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


def compare_etkdg_opls(etkdg, opls):

    # For each property.
    result_keys = result_key_defn()
    for rkey in result_keys:
        rkeyinfo = result_keys[rkey]
        if rkey == 'NN_dis':
            bottoms = [12, 0]
            density = True
            ylab = 'frequency'
            ylim = 24
        elif rkey == 'NCCN_dihed':
            bottoms = [0.1, 0]
            density = True
            ylab = 'frequency'
            ylim = 0.2
        else:
            continue

        fig, ax = plt.subplots(figsize=(8, 5))
        # For amine 1 only.
        etkdg_dict = etkdg['ami1']
        opls_dict = opls['ami1']

        # Get data.
        etkdg_data = {
            cid: (
                etkdg_dict[cid]['f_energy'],
                etkdg_dict[cid][rkey]
            )
            for cid in etkdg_dict
        }
        opls_data = {
            cid: (
                opls_dict[cid]['f_energy'],
                opls_dict[cid][rkey]
            )
            for cid in opls_dict
        }

        # Show all.
        etkdg_xdata = [etkdg_data[cid][1] for cid in etkdg_data]
        opls_xdata = [opls_data[cid][1] for cid in opls_data]

        # Plot data.
        xwidth = rkeyinfo['width']
        xbins = np.arange(
            rkeyinfo['lim'][0] - xwidth,
            rkeyinfo['lim'][1] + xwidth,
            xwidth
        )
        hist, bin_edges = np.histogram(
            a=etkdg_xdata,
            bins=xbins,
            density=density,
        )
        ax.bar(
            bin_edges[:-1],
            hist,
            bottom=bottoms[0],
            align='edge',
            alpha=1.0,
            width=xwidth,
            color='#FF5733',
            edgecolor='k',
            label='ETKDGv3',
        )

        hist, bin_edges = np.histogram(
            a=opls_xdata,
            bins=xbins,
            density=density,
        )
        ax.bar(
            bin_edges[:-1],
            hist,
            bottom=bottoms[1],
            align='edge',
            alpha=1.0,
            width=xwidth,
            color='#33DBFF',
            edgecolor='k',
            label='Macromodel/OPLS3e',
        )

        # ax.axvline(x=min_etkdg_s, c='#FF5733', lw=2, linestyle='--')
        # ax.axvline(x=min_opls_s, c='#33DBFF', lw=2, linestyle='--')

        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlim(rkeyinfo['lim'])
        ax.set_ylim(0, ylim)
        ax.set_xlabel(rkeyinfo['label'], fontsize=16)
        ax.set_ylabel(ylab, fontsize=16)
        if bottoms[0] != 0:
            ax.set_yticks([])
        ax.tick_params(left=False)
        ax.legend(fontsize=16, ncol=2)
        fig.tight_layout()
        fig.savefig(
            f'etkdgvsopls_{rkey}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


def scatter_of_all_conformers(
    result_dict,
    b_result_dict,
    filename,
    c,
    cb,
    plot_benzalde=True,
):

    result_keys = result_key_defn()

    for rkey in result_keys:
        if rkey == 'f_energy':
            continue

        rkeyinfo = result_keys[rkey]
        fig, ax = plt.subplots(figsize=(8, 5))

        if plot_benzalde:
            data = {
                cid: (
                    b_result_dict[cid]['f_energy'],
                    b_result_dict[cid][rkey]
                )
                for cid in b_result_dict
            }
            ydata = [data[cid][0] for cid in data]
            # Set to relative energies.
            ydata = [i-min(ydata) for i in ydata]
            # Set to kJ/mol.
            ydata = [i*2625.5 for i in ydata]
            xdata = [data[cid][1] for cid in data]
            ax.scatter(
                xdata,
                ydata,
                c=cb,
                alpha=1.0,
                edgecolor='white',
                marker='o',
                s=60,
                label='+ benzaldehyde',
            )

        data = {
            cid: (result_dict[cid]['f_energy'], result_dict[cid][rkey])
            for cid in result_dict
        }
        ydata = [data[cid][0] for cid in data]
        # Set to relative energies.
        ydata = [i-min(ydata) for i in ydata]
        # Set to kJ/mol.
        ydata = [i*2625.5 for i in ydata]
        xdata = [data[cid][1] for cid in data]
        ax.scatter(
            xdata,
            ydata,
            c=c,
            alpha=1.0,
            edgecolor='white',
            marker='o',
            s=60,
            label='free',
        )

        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlim(rkeyinfo['lim'])
        ax.set_ylim(result_keys['f_energy']['lim'])
        ax.set_xlabel(rkeyinfo['label'], fontsize=16)
        ax.set_ylabel(
            result_keys['f_energy']['label'],
            fontsize=16
        )
        ax.legend(fontsize=16)
        fig.tight_layout()
        fig.savefig(
            f'{filename}_{rkey}.pdf',
            dpi=720,
            bbox_inches='tight'
        )
        plt.close()


def save_lowest_energy_structures(results, name, dir):

    # Lowest 10 conformers by energy.
    low_e_conformers = [
        x
        for _, x in sorted(
            zip(
                [results[i]['f_energy'] for i in results],
                [i for i in results]
            ),
            key=lambda pair: pair[0]
        )
    ][:10]
    struct_files = [
        f'{dir}/conformer_{i}.mol' for i in low_e_conformers
    ]
    for i, s in enumerate(struct_files):
        mol = stk.BuildingBlock.init_from_file(s)
        mol.write(f'{name}_lowe_{i}.mol')


def main():
    if len(sys.argv) != 2:
        print(
            'structure_dir for OPLS3e structures must be given as'
            ' command line argument'
        )
        sys.exit()
    else:
        structure_dir = sys.argv[1]

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

    etkdg_results = run_etkdg_workflow(amines)
    for ami in amines:
        print(f'{ami} has {len(etkdg_results[ami])} etkdg conformers')

    if exists(structure_dir):
        opls3e_results = run_opls3e_workflow(amines, structure_dir)
    else:
        opls3e_results = None

    # Plot.
    for ami in amines:
        if ami[-1] == 'b':
            continue
        etkdg = etkdg_results[ami]
        etkdgb = etkdg_results[ami+'b']
        scatter_of_all_conformers(
            etkdg, etkdgb,
            filename=f'{ami}_etkdgall',
            c=leg_info()[ami]['c'],
            cb=leg_info()[ami+'b']['c'],
        )
        if opls3e_results is not None:
            opls3e = opls3e_results[ami]
            opls3eb = opls3e_results[ami+'b']
            scatter_of_all_conformers(
                opls3e, opls3eb,
                filename=f'{ami}_opls3eall',
                c=leg_info()[ami]['c'],
                cb=leg_info()[ami+'b']['c'],
                plot_benzalde=False,
            )

    density_of_all_amines(
        etkdg_results,
        filename='topetkdg_',
        ignore_b=False,
    )
    manuscript_density_of_all_amines(
        etkdg_results,
        filename='manu_etkdg_',
    )
    if opls3e_results is not None:
        density_of_all_amines(
            opls3e_results,
            filename='topopls3e',
            ignore_b=True,
        )
        compare_etkdg_opls(etkdg_results, opls3e_results)


if __name__ == '__main__':
    main()
