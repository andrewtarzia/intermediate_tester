#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Perform simple analysis and peak finding of mass spec data.

Author: Andrew Tarzia

Date Created: 17 Mar 2021

"""

import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks


def get_peaks(xdata, ydata, height):
    peaks, _ = find_peaks(ydata, height=height)
    x = xdata[peaks]
    return x


def get_clean_data(files, suffix):
    data_dict = {}
    total_counts = []
    for file in files:
        data = pd.read_csv(
            file,
            names=['#Point', 'X(Thompsons)', 'Y(Counts)'],
            header=None,
            index_col=['#Point']
        ).iloc[2:]
        data['x'] = data['X(Thompsons)'].astype(float)
        data['y'] = data['Y(Counts)'].astype(float)

        total_counts.append(sum(list(data['y'])))
        data_dict[file.replace(suffix, '')] = data

    # max_counts = max(total_counts)
    # count_ratios = [i/max_counts for i in total_counts]

    for i, d in enumerate(data_dict):
        # data_dict[d]['yrelative'] = ((
        # data_dict[d]['y']/sum(list(data_dict[d]['y']))
        # )*count_ratios[i] )
        # Use internal relative because that way we can just look at
        # at least 10 percent hit and mark peaks. Faulty experiment
        # is faulty..
        data_dict[d]['yrelative'] = (
            (data_dict[d]['y']/sum(list(data_dict[d]['y'])))
        )

    return data_dict


def plot_on_top(data_dict, out_name):
    fig, ax = plt.subplots(figsize=(8, 5))
    for wk in data_dict:
        ax.plot(
            data_dict[wk]['x'], data_dict[wk]['yrelative'], label=wk
        )
    ax.legend(fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('m/z', fontsize=16)
    ax.set_ylabel('relative ion counts', fontsize=16)
    ax.set_xlim(None, 1500)
    ax.set_title(out_name, fontsize=16)
    fig.tight_layout()
    fig.savefig(f'ontop_{out_name}.pdf', dpi=720, bbox_inches='tight')


def plot_peaks(data_dict, peaks, out_name):
    fig, ax = plt.subplots(figsize=(8, 5))
    for idx, wk in enumerate(data_dict):
        ax.plot(
            data_dict[wk]['x'],
            data_dict[wk]['yrelative']+(idx*0.2),
            label=wk,
        )
        ax.scatter(
            peaks[wk], [(idx*0.2)+0.1 for i in peaks[wk]], marker='*'
        )
    ax.legend(fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('m/z', fontsize=16)
    ax.set_ylabel('relative ion counts', fontsize=16)
    ax.set_xlim(250, 1500)
    ax.set_title(out_name, fontsize=16)
    fig.tight_layout()
    fig.savefig(f'peaks_{out_name}.pdf', dpi=720, bbox_inches='tight')


def save_peaks(peaks, out_name):

    for wk in peaks:
        string = '\n'.join([str(i) for i in peaks[wk]])
        with open(f'peaks_{out_name}_{wk}.txt', 'w') as f:
            f.write(string)


def main():
    dirs = {
        'iPr': {
            'name': '-',
            'files': [
                'Week1.CSV', 'Week2_2.CSV', 'Week3_2.CSV',
                'Week4.CSV', 'Week6.CSV'
            ],
        },
        'iPr & CC13': {
            'name': '+CC13',
            'files': [
                'Week1.CSV', 'Week2_2.CSV', 'Week3.CSV',
                'Week4.CSV', 'Week6.CSV'
            ],
        },
        'iPr & amine': {
            'name': '+amine',
            'files': [
                'Week1.CSV', 'Week2_2.CSV', 'Week3.CSV',
                'Week4.CSV', 'Week6.CSV'
            ],
        },
    }

    height = 0.01
    for dir in dirs:
        print(dir)
        full_files = [dir+'/'+i for i in dirs[dir]['files']]
        wks = get_clean_data(
            files=full_files,
            suffix='.CSV',
        )
        wks_updated = {
            i.replace(dir+'/', '').split('_')[0]: wks[i] for i in wks
        }
        plot_on_top(wks_updated, dirs[dir]['name'])

        # Get peaks and save.
        peaks = {}
        for idx, wk in enumerate(wks_updated):
            x = get_peaks(
                xdata=wks_updated[wk]['x'],
                ydata=wks_updated[wk]['yrelative'],
                height=height,
            )
            print(len(x), 'peaks found in', wk)
            peaks[wk] = list(x)

        plot_peaks(wks_updated, peaks, dirs[dir]['name'])
        save_peaks(peaks, dirs[dir]['name'])


if __name__ == '__main__':
    main()
