# intermediate_tester
Build and optimise cage intermediates

# Order to run intermediate analysis:
* Before running anything, run build_cages.py and manually build precursors!
* Run:
    * draw_all_intermediates.py
        * produces 'name'.svg image of structure from rdkit

    * structure_preparation.py
        * Runs opls3e conformer search and optimisation -> '{name}_ff.mol', '{name}_md.mol'
        * Then runs xTB optimisation -> '{name}_opt.{mol, ey}' [multiple ey files for multiple solvents]

    * write_gaussian_input.py
        * writes into xtb_spe_DFT/
        * writes input files for single point energy calculations

    * RUN single-point energy calculations through Gaussian

    * write_orca_input.py
        * writes into orca_calculations/
        * writes input files for B97-3c optimisations and single point energy calculations

    * RUN single-point energy and optimisation calculations through ORCA

    * An important note:
        * The higher-level methods for performing single-point energy evaluations and geometry optimisations showed method dependance
        * There were stability and convergence issues with the aug basis sets for MP2 calculations
        * The PBE1PBE single-points on xtb optimised structures and B97-3c energies from B97-3c optimised structures showed different relative reaction energies
        * Therefore, only the GFN2-xTB gas-phase energies were used in the manuscript

    * Collate all energies from all calculation sets using collate_energies.py
        * One command-line argument: output_ey_file: .csv file to save all energies too
        * edit this script to add new methods

    * calculate_formation_energies.py
        * produces two csv files with formation energies of all species for all methods.
        * arguments:
            * ey_file : all_total_energies.csv (from above)
            * formation_energy_file : file to output all formation energies too
            * per_imine_file : file to output all formation energies per imine bonds formed too

    * plot_formation_energies.py fe_csv_file fe_pi_csv_file
        * fe_csv_file: file containing formation energies to plot
        * fe_pi_csv_file: file containing formation energies per imine to plot
        * produces plots in figures/ directory:
            * fe_perimine.pdf
            * fe_total.pdf
            * fe_perimine_ami{i}.pdf

# Conformer analysis:
* conformer_analysis.py will run conformer generation with the ETKDGv3 algorithm in `rdkit` plot the distribution of geometrical properties
    * Amines are defined in the script for the molecules to analyse by their SMILES string
    * All conformers optimised with GFN2-xTB (version 6.3.2)
    * All energies and free energies calculated with GFN2-xTB (version 6.3.2)
    * If conformer generation with macromodel/maestro is performed manually, these can be linked in for analysis (a CLI option for the script to this directory allows this)
        * All conformers optimised with GFN2-xTB (version 6.3.2)
    * The 10 lowest energy conformers will be saved as {amine}\_lowe\_{i}.mol in the working directory
