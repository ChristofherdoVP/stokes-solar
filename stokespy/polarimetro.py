"""Modulação, demodulação e eficiência polarimétrica (Seção 4 do artigo)."""

import numpy as np

from .polarizacao import mueller_polarizer, mueller_retarder


def modulation_matrix(angles, Delta, analyzer=0.0):
    """Matriz O (n x 4) de um retardador girante seguido de polarizador fixo.

    angles   : orientações do eixo rápido do retardador (rad)
    Delta    : retardo (rad); pi/2 para lâmina de quarto de onda
    analyzer : ângulo do polarizador analisador (rad)
    Cada linha é a primeira linha de M_ana @ M_mod,i, normalizada para que o
    termo em I valha 1 (um analisador ideal transmitiria metade da luz).
    """
    Mana = mueller_polarizer(analyzer)
    O = np.array([(Mana @ mueller_retarder(t, Delta))[0] for t in np.atleast_1d(angles)])
    return O / O[:, :1]


def demodulation_matrix(O):
    """Matriz D = (O^T O)^-1 O^T (pseudoinversa)."""
    return np.linalg.pinv(O)


def efficiency(O):
    """Eficiências (eps_I, eps_Q, eps_U, eps_V) de del Toro Iniesta & Collados (2000)."""
    D = demodulation_matrix(O)
    n = O.shape[0]
    return 1.0 / np.sqrt(n * np.sum(D**2, axis=1))


def simulate_measurements(S, O, photons=1e6, rng=None):
    """Intensidades medidas com ruído de fótons (Poisson).

    S       : vetor de Stokes normalizado (I = 1)
    photons : número médio de fótons por exposição para I = 1
    Devolve as intensidades em unidades de I.
    """
    rng = np.random.default_rng(rng)
    expected = photons * (O @ S)
    return rng.poisson(expected) / photons
