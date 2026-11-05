"""Parâmetros de Stokes e matrizes de Mueller (Seções 2 e 3 do artigo)."""

import numpy as np


def stokes_from_field(ax, ay, delta):
    """Vetor de Stokes de luz totalmente polarizada.

    ax, ay : amplitudes das componentes do campo elétrico
    delta  : diferença de fase phi_y - phi_x (rad)
    """
    return np.array([ax**2 + ay**2,
                     ax**2 - ay**2,
                     2 * ax * ay * np.cos(delta),
                     2 * ax * ay * np.sin(delta)])


def stokes_from_ellipse(psi, chi, p=1.0, I=1.0):
    """Vetor de Stokes a partir do azimute psi e da elipticidade chi (rad)."""
    return np.array([I,
                     p * I * np.cos(2 * chi) * np.cos(2 * psi),
                     p * I * np.cos(2 * chi) * np.sin(2 * psi),
                     p * I * np.sin(2 * chi)])


def degree_of_polarization(S):
    """Graus de polarização total, linear e circular."""
    I, Q, U, V = S
    p = np.sqrt(Q**2 + U**2 + V**2) / I
    return p, np.sqrt(Q**2 + U**2) / I, np.abs(V) / I


def mueller_polarizer(theta):
    """Polarizador linear ideal com eixo de transmissão no ângulo theta (rad)."""
    c, s = np.cos(2 * theta), np.sin(2 * theta)
    return 0.5 * np.array([[1, c, s, 0],
                           [c, c * c, c * s, 0],
                           [s, c * s, s * s, 0],
                           [0, 0, 0, 0]])


def mueller_retarder(theta, Delta):
    """Retardador linear com eixo rápido em theta e retardo Delta (rad).

    Forma consistente com delta = phi_y - phi_x (Seção 3.2 do artigo).
    """
    c, s = np.cos(2 * theta), np.sin(2 * theta)
    cd, sd = np.cos(Delta), np.sin(Delta)
    return np.array([[1, 0, 0, 0],
                     [0, c * c + s * s * cd, c * s * (1 - cd), s * sd],
                     [0, c * s * (1 - cd), s * s + c * c * cd, -c * sd],
                     [0, -s * sd, c * sd, cd]])


def rotation(alpha):
    """Rotação do sistema de referência por alpha (rad)."""
    c, s = np.cos(2 * alpha), np.sin(2 * alpha)
    return np.array([[1, 0, 0, 0],
                     [0, c, s, 0],
                     [0, -s, c, 0],
                     [0, 0, 0, 1]])


def rotate_element(M0, theta):
    """Matriz de um elemento girado por theta: R(-theta) M0 R(theta)."""
    return rotation(-theta) @ M0 @ rotation(theta)
