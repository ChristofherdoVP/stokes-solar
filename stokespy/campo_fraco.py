"""Aproximação de campo fraco (Seção 7 do artigo)."""

import numpy as np

from .zeeman import C_ZEEMAN


def wfa_blos(lam, I, V, lam0, g, mask=None):
    """B_LOS (G) por mínimos quadrados: V = -C B_LOS dI/dlambda.

    Aceita perfis 1D (n_lambda) ou cubos (..., n_lambda).
    mask : seleção booleana de pontos espectrais (por exemplo, sem linhas telúricas)
    """
    dI = np.gradient(I, lam, axis=-1)
    if mask is not None:
        dI, V = dI[..., mask], V[..., mask]
    C = C_ZEEMAN * g * lam0**2
    return -np.sum(V * dI, axis=-1) / (C * np.sum(dI**2, axis=-1))


def wfa_btrans(lam, I, Q, U, lam0, G):
    """B_perp (G) e azimute chi (rad) pela aproximação de segunda ordem.

    Q = -(1/4) G (C' lam0^2 B_perp)^2 cos(2chi) d2I/dlambda2, com C' = 4,67e-13.
    G : fator de Landé transversal (g_ef^2 para um tripleto normal).
    """
    d2I = np.gradient(np.gradient(I, lam, axis=-1), lam, axis=-1)
    k = -0.25 * G * (C_ZEEMAN * lam0**2)**2
    den = k * np.sum(d2I**2, axis=-1)
    Bq2 = np.sum(Q * d2I, axis=-1) / den   # B_perp^2 cos 2chi
    Bu2 = np.sum(U * d2I, axis=-1) / den   # B_perp^2 sin 2chi
    return np.hypot(Bq2, Bu2) ** 0.5, 0.5 * np.arctan2(Bu2, Bq2)
