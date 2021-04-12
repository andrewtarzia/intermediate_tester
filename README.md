# intermediate_tester
Build and optimise cage intermediates

# Order:
* Before running anything, run build_cages.py and manually build precursors!
* Run:
    * draw_all_intermediates.py
        * produces 'name'.svg image of structure from rdkit
    * structure_preparation.py
        * Runs opls3e conformer search and optimisation -> '{name}_ff.mol', '{name}_md.mol'
        * Then runs xTB optimisation -> '{name}_opt.{mol, ey}' [multiple ey files for multiple solvents]
    * calculate_formation_energies.py
        * produces XX.csv with formation energies of all species
    * plot_formation_energies.py
        * produces plots:
            * fe_perimine.pdf
            * fe_total.pdf
            * fe_perimine_ami{i}.pdf

# Conformer analysis:
* conformer_analysis.py will run conformer generation with the ETKDG algorithm in `rdkit` plot the distribution of geometrical properties
    * All energies and free energies calculated performed with GFN2-xTB (version 6.3.2)
    * If conformer generation with macromodel/maestro is performed manually, these can be linked in for analysis
    * The 10 lowest energy conformers will be saved as {amine}_lowe_{i}.mol in the working directory
