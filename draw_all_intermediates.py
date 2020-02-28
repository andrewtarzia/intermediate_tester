#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Draw all intermediates in rdkit.

Author: Andrew Tarzia

Date Created: 28 Feb 2020

"""

import stk
from rdkit.Chem import AllChem as rdkit
import glob
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import DrawingOptions


def draw_mol_to_svg(mol, filename):
    """
    Draw a single molecule to an SVG file with transparent BG.

    """
    # change BG to transperent
    # (https://sourceforge.net/p/rdkit/mailman/message/31637105/)
    o = DrawingOptions()
    o.bgColor = None
    # Use copy of molecule to avoid changing instance of mol.
    new_mol = rdkit.MolFromMolBlock(rdkit.MolToMolBlock(mol))
    rdkit.Compute2DCoords(new_mol)
    Draw.MolToFile(
        new_mol, filename, fitImage=True, imageType='svg', options=o
    )


def main():
    for mol_file in sorted(glob.glob('*a*.mol')):
        name = mol_file.replace('.mol', '')
        if 'opt' in name or 'xtb' in name:
            continue

        # Read molecule into stk.
        mol = stk.BuildingBlock.init_from_file(mol_file)
        print(mol)
        smi = rdkit.MolToSmiles(mol.to_rdkit_mol())
        print(smi)
        rdk_mol = rdkit.MolFromSmiles(smi)

        draw_mol_to_svg(rdk_mol, f'{name}.svg')


if __name__ == '__main__':
    main()
