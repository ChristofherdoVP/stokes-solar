"""stokespy: funções didáticas de espectropolarimetria solar.

Acompanha o artigo "Espectropolarimetria solar: dos parâmetros de Stokes ao
campo magnético fotosférico, com Python".

Convenções (Seção 2 do artigo):
    * vetor de Stokes S = (I, Q, U, V)
    * E_x = a_x cos(wt - phi_x), E_y = a_y cos(wt - phi_y), delta = phi_y - phi_x
    * V > 0 para 0 < delta < pi (rotação anti-horária vista pelo observador, IAU)
    * comprimentos de onda em angstrom, campo magnético em gauss
"""

from .polarizacao import (
    stokes_from_field, stokes_from_ellipse, degree_of_polarization,
    mueller_polarizer, mueller_retarder, rotation, rotate_element,
)
from .polarimetro import modulation_matrix, demodulation_matrix, efficiency, simulate_measurements
from .zeeman import (
    C_ZEEMAN, lande_g, g_effective, zeeman_splitting, doppler_width,
    voigt_faraday, propagation_matrix, me_synthesis, LINHAS,
)
from .campo_fraco import wfa_blos, wfa_btrans

__all__ = [
    "stokes_from_field", "stokes_from_ellipse", "degree_of_polarization",
    "mueller_polarizer", "mueller_retarder", "rotation", "rotate_element",
    "modulation_matrix", "demodulation_matrix", "efficiency", "simulate_measurements",
    "C_ZEEMAN", "lande_g", "g_effective", "zeeman_splitting", "doppler_width",
    "voigt_faraday", "propagation_matrix", "me_synthesis", "LINHAS",
    "wfa_blos", "wfa_btrans",
]
