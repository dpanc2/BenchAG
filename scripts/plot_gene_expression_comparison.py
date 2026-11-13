#!/usr/bin/env python3
"""Compare the three gene-expression metric distributions from saved evaluations."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evaluation', action='append', required=True, help='Label=/path/to/evaluation; repeat for each species/model')
    parser.add_argument('--out', required=True, help='Output PNG')
    parser.add_argument('--title', default='Gene-level expression: independent TPM benchmark')
    args = parser.parse_args()
    metrics = ['cross_gene_raw', 'cross_gene_normalized', 'cross_track_normalized']
    colors = ['#e76f51', '#2a9d8f', '#457b9d', '#9b5de5']
    fig, ax = plt.subplots(figsize=(8, 5.2))
    handles = []
    width = .7 / len(args.evaluation)
    for k, item in enumerate(args.evaluation):
        label, directory = item.split('=', 1)
        scores = pd.read_csv(Path(directory) / 'correlations.csv')
        summary = pd.read_csv(Path(directory) / 'summary.csv')
        color = colors[k % len(colors)]
        legend_label = f'{label} ({summary.n_genes.iloc[0]:,} genes, {summary.n_tracks.iloc[0]} tracks)'
        handles.append(Line2D([0], [0], marker='o', color=color, linewidth=0, label=legend_label))
        for j, metric in enumerate(metrics, 1):
            values = scores.loc[scores.metric == metric, 'pearson_r'].dropna().to_numpy()
            if not len(values):
                raise ValueError(f'{label}: no correlations for {metric}')
            pos = j + (k - (len(args.evaluation) - 1) / 2) * width
            if len(values) > 1 and np.ptp(values) > 1e-12:
                violin = ax.violinplot(values, positions=[pos], widths=width*.9, showextrema=False)
                violin['bodies'][0].set_facecolor(color)
                violin['bodies'][0].set_alpha(.65)
            ax.scatter(pos, values.mean(), color=color, s=30, zorder=3, edgecolor='white', linewidth=.5)
            ax.annotate(f'{values.mean():.3f}', (pos, values.mean()), xytext=(-5 if k == 0 else 5, 14 if k == 0 else 3),
                        textcoords='offset points', ha='right' if k == 0 else 'left', fontsize=10)
    ax.set_xticks([1, 2, 3], ['Raw\nCross-gene', 'Normalized\nCross-gene', 'Normalized\nCross-track'])
    ax.set_ylabel('Pearson r per track or gene')
    ax.set_ylim(-1.02, 1.06)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_title(args.title, fontsize=12, pad=12)
    ax.axhline(0, color='#dddddd', linewidth=.8, zorder=0)
    fig.legend(handles=handles, loc='lower center', frameon=False, ncol=min(2, len(handles)), fontsize=9)
    fig.tight_layout(rect=(0, .07, 1, 1))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
