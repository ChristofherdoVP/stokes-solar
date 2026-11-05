"""Efeito Zeeman e síntese de Milne-Eddington (Seções 5 e 6 do artigo).

Unidades: comprimento de onda em angstrom, campo em gauss, velocidade em km/s.
"""

import numpy as np
from scipy.special import wofz

C_ZEEMAN = 4.6686e-13  # A^-1 G^-1: Delta lambda_B = C g lambda0^2 B
C_LUZ = 2.99792458e5   # km/s
K_B = 1.380649e-23     # J/K
UMA = 1.66053907e-27   # kg

# Linhas usadas no artigo: lambda0 (A), g_ef, massa atômica (u)
LINHAS = {
    "Fe I 630,15 nm": dict(lam0=6301.50, g=1.67, massa=55.845),
    "Fe I 630,25 nm": dict(lam0=6302.49, g=2.50, massa=55.845),
    "Fe I 617,3 nm":  dict(lam0=6173.33, g=2.50, massa=55.845),
    "Fe I 1564,8 nm": dict(lam0=15648.5, g=3.00, massa=55.845),
}


def lande_g(S, L, J):
    """Fator de Landé de um nível em acoplamento LS."""
    if J == 0:
        return 0.0
    return 1.5 + (S * (S + 1) - L * (L + 1)) / (2 * J * (J + 1))


def g_effective(gu, gl, Ju, Jl):
    """Fator de Landé efetivo de uma transição."""
    return 0.5 * (gu + gl) + 0.25 * (gu - gl) * (Ju * (Ju + 1) - Jl * (Jl + 1))


def zeeman_splitting(lam0, g, B):
    """Separação Zeeman Delta lambda_B (A)."""
    return C_ZEEMAN * g * lam0**2 * B


def doppler_width(lam0, T, xi=0.0, massa=55.845):
    """Largura Doppler (A) para temperatura T (K) e microturbulência xi (km/s)."""
    v_th = np.sqrt(2 * K_B * T / (massa * UMA)) / 1e3
    return lam0 / C_LUZ * np.sqrt(v_th**2 + xi**2)


def voigt_faraday(a, v):
    """Perfis de Voigt H(a, v) e Faraday-Voigt F(a, v), com w(z) = H + 2iF."""
    w = wofz(v + 1j * a)
    return w.real, 0.5 * w.imag


def propagation_matrix(lam, lam0, g, B, gamma, chi, vlos, dlD, eta0, a):
    """Matriz de propagação K (n_lambda x 4 x 4), normalizada pelo contínuo.

    Tripleto equivalente com fator de Landé efetivo g. Ângulos em rad.
    Na convenção adotada, V tem o sinal de -B cos(gamma) dI/dlambda
    (Seção 7 do artigo).
    """
    v = (lam - lam0 - lam0 * vlos / C_LUZ) / dlD
    vB = zeeman_splitting(lam0, g, B) / dlD
    phi_p, psi_p = voigt_faraday(a, v)
    phi_b, psi_b = voigt_faraday(a, v + vB)   # componente sigma deslocada para o azul
    phi_r, psi_r = voigt_faraday(a, v - vB)   # componente sigma deslocada para o vermelho

    sg2, cg = np.sin(gamma)**2, np.cos(gamma)
    c2x, s2x = np.cos(2 * chi), np.sin(2 * chi)

    def coef(p, b, r):
        I = 0.5 * (p * sg2 + 0.5 * (b + r) * (1 + cg**2))
        Ql = 0.5 * (p - 0.5 * (b + r)) * sg2
        return I, Ql * c2x, Ql * s2x, 0.5 * (r - b) * cg

    eI, eQ, eU, eV = coef(phi_p, phi_b, phi_r)
    _, rQ, rU, rV = coef(psi_p, psi_b, psi_r)

    n = lam.size
    K = np.zeros((n, 4, 4))
    K[:, 0, 0] = K[:, 1, 1] = K[:, 2, 2] = K[:, 3, 3] = 1 + eta0 * eI
    K[:, 0, 1] = K[:, 1, 0] = eta0 * eQ
    K[:, 0, 2] = K[:, 2, 0] = eta0 * eU
    K[:, 0, 3] = K[:, 3, 0] = eta0 * eV
    K[:, 1, 2], K[:, 2, 1] = eta0 * rV, -eta0 * rV
    K[:, 1, 3], K[:, 3, 1] = -eta0 * rU, eta0 * rU
    K[:, 2, 3], K[:, 3, 2] = eta0 * rQ, -eta0 * rQ
    return K


def me_synthesis(lam, lam0, g, B, gamma, chi, vlos, dlD, eta0, a, S0, S1, mu=1.0):
    """Perfis de Stokes emergentes de uma atmosfera de Milne-Eddington.

    I(0) = (S0 + mu S1 K^-1) e0  (forma matricial de Unno-Rachkovsky).
    Devolve um array (4, n_lambda) com I, Q, U, V.
    """
    K = propagation_matrix(lam, lam0, g, B, gamma, chi, vlos, dlD, eta0, a)
    e0 = np.array([1.0, 0.0, 0.0, 0.0])
    sol = np.linalg.solve(K, np.broadcast_to(e0, (lam.size, 4))[..., None])[..., 0]
    out = mu * S1 * sol
    out[:, 0] += S0
    return out.T
