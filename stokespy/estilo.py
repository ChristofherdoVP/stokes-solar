"""Estilo das figuras do artigo: tons legíveis em impressão e vírgula decimal."""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

CORES = {"I": "#222222", "Q": "#1f5fa8", "U": "#2a9d55", "V": "#c0392b",
         "cinza": "#8a8a8a", "destaque": "#d98a00"}


def _virgula(x, pos=None):
    s = f"{x:g}"
    return s.replace("-", "−").replace(".", ",")


VIRGULA = FuncFormatter(_virgula)


COLUNA = 3.35   # largura de uma coluna da RBEF, em polegadas (8,5 cm)
PAGINA = 7.0    # largura de duas colunas


def setup(base=10):
    """base=8 para figuras do artigo (largura de coluna); 10 para o suplementar."""
    mpl.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": base, "axes.labelsize": base, "axes.titlesize": base,
        "legend.fontsize": base - 1, "xtick.labelsize": base - 1, "ytick.labelsize": base - 1,
        "lines.linewidth": 1.2 if base < 10 else 1.4, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "mathtext.fontset": "dejavuserif",
        "font.family": "DejaVu Serif", "axes.unicode_minus": True,
    })


def virgula(ax, eixos="xy"):
    """Aplica vírgula decimal aos eixos indicados."""
    if "x" in eixos:
        ax.xaxis.set_major_formatter(VIRGULA)
    if "y" in eixos:
        ax.yaxis.set_major_formatter(VIRGULA)


def salvar(fig, nome, pasta="../figuras"):
    import os
    os.makedirs(pasta, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{pasta}/{nome}.{ext}")
    plt.show()
