import numpy as np, stokespy as sp
lam0, g = 6302.49, 2.5
lam = lam0 + np.linspace(-0.6, 0.6, 1201)
dlD = sp.doppler_width(lam0, 6000, 1.0)
par = dict(vlos=0.0, dlD=dlD, eta0=10.0, a=0.1, S0=0.3, S1=0.7)

def test_luz_nao_polarizada():
    S = 0.5*np.array([1,0,0,0])
    assert sp.degree_of_polarization(np.array([1,0,0,0.]))[0] == 0

def test_p_menor_que_1():
    rng = np.random.default_rng(0)
    for _ in range(100):
        S = sp.stokes_from_ellipse(*rng.uniform(-1,1,2), p=rng.uniform())
        M = sp.mueller_retarder(*rng.uniform(0,3,2)) @ sp.mueller_polarizer(rng.uniform(0,3))
        S2 = M @ S
        if S2[0] > 1e-12: assert sp.degree_of_polarization(S2)[0] <= 1 + 1e-9

def test_rotacao_consistente():
    for t in np.linspace(0, np.pi, 7):
        assert np.allclose(sp.rotate_element(sp.mueller_retarder(0, 1.1), t), sp.mueller_retarder(t, 1.1))
        assert np.allclose(sp.rotate_element(sp.mueller_polarizer(0), t), sp.mueller_polarizer(t))

def test_exemplo_lamina():
    S = sp.mueller_retarder(np.pi/4, np.pi/2) @ sp.mueller_polarizer(0) @ np.array([1,0,0,0.])
    assert np.allclose(S, [0.5, 0, 0, -0.5])

def test_otimo_132():
    O = sp.modulation_matrix(np.deg2rad([15.12, -15.12, 51.69, -51.69]), np.deg2rad(131.8))
    assert np.allclose(sp.efficiency(O)[1:], 1/np.sqrt(3), atol=1e-3)

def test_B_zero():
    S = sp.me_synthesis(lam, lam0, g, 0.0, 0.5, 0.3, **par)
    assert np.abs(S[1:]).max() < 1e-12

def test_campo_fraco_recupera():
    for gam in [0.0, 0.6, 2.5]:
        S = sp.me_synthesis(lam, lam0, g, 20.0, gam, 0.0, **par)
        b = sp.wfa_blos(lam, S[0], S[3], lam0, g)
        assert abs(b - 20*np.cos(gam)) < 0.01*20

def test_lamina_continua_simetria_QU():
    S = sp.me_synthesis(lam, lam0, g, 1000.0, np.pi/2, 0.0, **par)
    assert np.abs(S[2]).max() < 1e-12 and np.abs(S[3]).max() < 1e-12
    assert np.allclose(S[1], S[1][::-1], atol=1e-6)
