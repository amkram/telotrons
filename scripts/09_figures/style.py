"""Shared paths, colours, fonts and helpers of the figure scripts."""

from pathlib import Path
import csv
import json
import os

import matplotlib

# render without a display (must precede the pyplot import)
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

PACKAGE = Path(os.environ['TELO_PACKAGE'])
for directory in ['figures', 'qa']:
    (PACKAGE / directory).mkdir(parents=True, exist_ok=True)

# sizes in inches: one millimetre, and the 180 mm figure width
MM = 1 / 25.4
WIDTH = 180 * MM

# text, rule, F/R repeat and TSD colours
INK = '#222222'
MUTED = '#666666'
RULE = '#aaaaaa'
F_COLOR = '#EF6B35'
R_COLOR = '#2A78D6'
TSD_COLOR = '#CFEDE3'

# lineages in display order, with their labels and colours
CLADES = [
    'haptophyte',
    'cryptophyte',
    'dinoflagellate',
    'chrysophyte',
    'apicomplexa',
    'nematode',
    'stramenopile',
]
NAMES = [
    'Haptophytes',
    'Cryptophyte',
    'Dinoflagellates',
    'Chrysophyte',
    'Eimeria',
    'Mesorhabditis',
    'Other stramenopile',
]
COLORS = ['#6F4C9B', '#0072B2', '#D55E00', '#E69F00', '#009E73', '#CC79A7', '#888888']

# the six classified architectures and their colours
ARCH = ['R', 'F', 'R-F', 'F-R', 'R-L-F', 'F-L-R']
ARCH_COLORS = {
    'R': R_COLOR,
    'F': F_COLOR,
    'R-F': '#9C728E',
    'F-R': '#5C9E8C',
    'R-L-F': '#8F79B6',
    'F-L-R': '#B2CBAA',
}

# taxon names italicised at the start of a label, longest first so a binomial wins over its genus
TAXA = sorted(
    [
        'Pristionchus',
        'Aureococcus',
        'Fragilariopsis',
        'Mesorhabditis spiculigera',
        'Pristionchus pacificus',
        'Baffinella frigidus',
        'Diacronema lutheri',
        'Emiliania huxleyi',
        'Phaeodactylum tricornutum',
        'Durusdinium trenchii',
        'Toxoplasma gondii',
        'Haptolina',
        'Symbiodinium natans',
        'Cyclospora cayetanensis',
        'Eimeria acervulina',
        'Eimeria tenella',
        'Phaeocystis',
        'Emiliania',
        'Symbiodinium',
        'Eimeria',
        'Mesorhabditis',
    ],
    key=len,
    reverse=True,
)

plt.rcParams.update(
    {
        'font.family': 'DejaVu Sans',
        'font.size': 7.5,
        'axes.labelsize': 7.5,
        'axes.titlesize': 7.5,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7.5,
        'legend.fontsize': 7,
        'axes.linewidth': 0.6,
        'lines.linewidth': 1.1,
        # keep text editable in PDF and SVG output
        'pdf.fonttype': 42,
        'svg.fonttype': 'none',
        'text.color': INK,
        'axes.labelcolor': INK,
        'xtick.color': INK,
        'ytick.color': INK,
        'savefig.facecolor': 'white',
    }
)


def read(name, folder='tables'):
    """Rows of a package table."""
    with (PACKAGE / folder / name).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def read_inputs(name):
    """Rows of a curated input under inputs/."""
    return read(name, 'inputs')


def save(fig, stem):
    for extension in ['pdf', 'svg', 'png']:
        fig.savefig(PACKAGE / 'figures' / f'{stem}.{extension}', dpi=450)
    plt.close(fig)


def write_qa(name, data):
    """A layout or content record under qa/."""
    (PACKAGE / 'qa' / name).write_text(json.dumps(data, indent=2) + '\n')


def axes_style(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color(RULE)
    ax.tick_params(width=0.6, length=2.5, pad=2)


def panel(fig, label, x, y):
    fig.text(x, y, label, fontsize=10, fontweight='bold', va='top')


def display_name(label):
    """Italicise the taxon name that starts the label, in mathtext."""
    label = label.replace('(n=', '(n = ')
    for taxon in TAXA:
        if label.startswith(taxon):
            rest = label[len(taxon) :].replace('MAG_', 'MAG ')
            return r'$\it{' + taxon.replace(' ', r'\ ') + r'}$' + rest
    return label.replace('MAG_', 'MAG ')
