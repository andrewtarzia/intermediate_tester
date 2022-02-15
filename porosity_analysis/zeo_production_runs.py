import os
import matplotlib.pyplot as plt


def convert_zsa_to_xyz(file):
    """
    Convert .zsa coordinates into XYZ file for visualisation.

    """
    with open(file, 'r') as f:
        data = f.readlines()

    for i, j in enumerate(data):
        if 'color red' in j:
            red_mention = i

    greens = data[1:red_mention]
    reds = data[red_mention+1:]

    n_atoms = len(greens) + len(reds)
    xyz_file = file.replace('.zsa', '_z.xyz')

    with open(xyz_file, 'w') as f:
        f.write(f'{n_atoms}\nWritten by Andrew Tarzia!\n')
        for g in greens:
            id = 'H'
            D = g.rstrip().replace('{', '').replace('}', '')
            x, y, z = [
                i for i in D.replace('point', '').split(' ') if i
            ]
            f.write(f'{id} {x} {y} {z}\n')
        for g in reds:
            id = 'P'
            D = g.rstrip().replace('{', '').replace('}', '')
            x, y, z = [
                i for i in D.replace('point', '').split(' ') if i
            ]
            f.write(f'{id} {x} {y} {z}\n')


def get_values(lines, file):
    values = {}
    line_of_int = lines[0]
    nl = line_of_int.replace(file, '').replace('@  ', '').rstrip()
    nl = [i for i in nl.split(' ') if i]
    values['volume'] = float(nl[1])
    values['density'] = float(nl[3])
    values['asa_m2g-1'] = float(nl[9])
    values['nasa_m2g-1'] = float(nl[-1])

    return values


def main():

    crystals = [
        'CC21-alpha_100K_publ.cif',
    ]
    zeo_path = '/home/atarzia/software/zeo++-0.3/network'
    probes = [1.0, 1.2, 1.3, 1.4, 1.55, 1.6, 1.7, 1.82]
    samplins = [50000]

    results = {}
    for cryst in crystals:
        results[cryst] = {}
        for probe in probes:
            results[cryst][probe] = {}
            for sampl in samplins:
                output = cryst.replace('.cif', f'_{probe}_{sampl}.out')
                cmd = (
                    f'{zeo_path} -ha -sa {probe} {probe} {sampl} '
                    f'{output} {cryst}'
                )
                if not os.path.exists(output):
                    os.system(cmd)

                with open(output, 'r') as f:
                    lines = f.readlines()
                values = get_values(lines, output)
                results[cryst][probe][sampl] = values

                # Do zsa.
                output = cryst.replace(
                    '.cif', f'_{probe}_{sampl}_zsa.out'
                )
                cmd = (
                    f'{zeo_path} -ha -zsa {probe} {probe} {sampl} '
                    f'{output} {cryst}'
                )
                if not os.path.exists(output):
                    os.system(cmd)
                    convert_zsa_to_xyz(output)

    print(results)
    for cryst in results:
        for probe in results[cryst]:
            for sampl in results[cryst][probe]:
                da = results[cryst][probe][sampl]
                asa = da['asa_m2g-1']
                nasa = da['nasa_m2g-1']
                print(
                    f'{cryst}, {probe}, {sampl}: '
                    f'ASA={asa}, NASA={nasa}'
                )

    fig, ax = plt.subplots(figsize=(8, 5))

    markers = ['o', 'X', 'P', 'D', '<', '>']
    for i, res in enumerate(results):
        sres = results[res]
        for probe in sres:
            da = results[res][probe]
            print(da)
            nasas = [da[i]['nasa_m2g-1'] for i in da]
            asas = [da[i]['asa_m2g-1'] for i in da]
            print(nasas, asas)
            if nasas[-1] != 0.0:
                ax.scatter(
                    [i for i in da],
                    # [((sa-nasas[-1])/nasas[-1])*100 for sa in nasas],
                    [sa for sa in nasas],
                    c='k',
                    marker=markers[i],
                    label=f'{res}-{probe}-NASA',
                    edgecolor='k',
                    s=50,
                )
            else:
                ax.scatter(
                    [i for i in da],
                    # [((sa-nasas[-1])) for sa in nasas],
                    [sa for sa in nasas],
                    c='k',
                    marker=markers[i],
                    label=f'{res}-{probe}-NASA',
                    edgecolor='k',
                    s=50,
                )
            if asas[-1] != 0.0:
                ax.scatter(
                    [i for i in da],
                    # [((asa-asas[-1])/asas[-1])*100 for asa in asas],
                    [sa for sa in asas],
                    c='r',
                    marker=markers[i],
                    label=f'{res}-{probe}-ASA',
                    edgecolor='k',
                    s=50,
                )
            else:
                ax.scatter(
                    [i for i in da],
                    # [((sa-asas[-1])) for sa in asas],
                    [sa for sa in asas],
                    c='r',
                    marker=markers[i],
                    label=f'{res}-{probe}-ASA',
                    edgecolor='k',
                    s=50,
                )

    ax.legend(fontsize=16)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('sampling', fontsize=16)
    ax.set_ylabel('SA [m2/g]', fontsize=16)

    fig.tight_layout()
    fig.savefig('SA_vs_probe.pdf', dpi=720, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    main()
