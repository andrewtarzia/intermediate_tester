import os


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


def get_res_values(lines, file):
    values = {}
    line_of_int = lines[0]
    nl = line_of_int.replace(file, '').replace('@  ', '').rstrip()
    nl = [i for i in nl.split(' ') if i]
    print(nl)
    values['di'] = float(nl[0])
    values['df'] = float(nl[1])
    values['dif'] = float(nl[2])

    return values


def main():

    crystals = [
        'CC21-alpha_100K_publ.cif',
    ]
    zeo_path = '/home/atarzia/software/zeo++-0.3/network'
    probes = [1.0, 1.55, 1.82]
    samplins = [1000, 5000, 10000, 20000, 30000, 50000, 100000, 200000]

    results = {}
    for cryst in crystals:
        resoutput = cryst.replace('.cif', '.res')
        cmd = (
            f'{zeo_path} -ha -res '
            f'{resoutput} {cryst}'
        )
        if not os.path.exists(resoutput):
            os.system(cmd)

        with open(resoutput, 'r') as f:
            lines = f.readlines()
        values = get_res_values(lines, resoutput)
        print(f'{cryst}: {values}')

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


if __name__ == '__main__':
    main()
