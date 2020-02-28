# intermediate_tester
Build and optimise cage intermediates

# Order:
* Before running anything, run build_cages.py and manually build precursors!
* Run:
    * draw_all_intermediates.py
        * produces 'name'.svg image of structure from rdkit
    * FF_optimise_all_structures.py
        * produces 'name'_opt.mol
    * xtb_opt_and_energy.py
        * produces 'name'_opt_xtb.{mol, ey, fey}
    * calculate_formation_energies.py
        * produces plots:
            * fe_perimine.pdf
            * fe_total.pdf
            * fe_perimine_ami{i}.pdf
