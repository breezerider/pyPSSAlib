import numpy as np
import pypssalib as m


def test_version():
    assert m.__version__ == "0.1.0-dev0"


def test_method():
    pssa = m.pSSAlib(m.SSA.SPDM)
    assert pssa.method == m.SSA.SPDM

    pssa.method = m.SSA.PDM
    assert pssa.method == m.SSA.PDM


def test_run_tcs_dm(monkeypatch):
    monkeypatch.setenv("GSL_RNG_SEED", "01062024")
    expected = np.loadtxt("tests/tcs_dm.dat", dtype=int)

    I_0 = 50
    HK = 100
    RR = 100

    pssa = m.pSSAlib(m.SSA.DM)
    actual = pssa.sample_testcase_population(m.Testcase.TCS, [0.2, 6.5, 0.6, 0.2, 4.0, 4.0, I_0, HK, 0, RR, 0], 100, 10)
    actual = actual.squeeze()
    assert (expected == actual).all()


def test_run_sbd_pdm(monkeypatch):
    monkeypatch.setenv("GSL_RNG_SEED", "01062024")
    expected = np.loadtxt("tests/sbd_pdm.dat", dtype=int)

    pssa = m.pSSAlib(m.SSA.PDM)
    actual = pssa.sample_testcase_trajectory(
        m.Testcase.SingleBirthDeath, [1, 1, 10, 0], 2100, time_start=2000, time_step=1
    )
    actual = actual.squeeze()
    assert (expected == actual).all()


def test_run_ca_spdm(monkeypatch):
    monkeypatch.setenv("GSL_RNG_SEED", "01062024")
    expected = np.loadtxt("tests/ca_spdm.dat", dtype=int)

    pssa = m.pSSAlib(m.SSA.SPDM)
    actual = pssa.sample_testcase_trajectory(
        m.Testcase.ColloidalAggregation, [2, 2.1, 0.1, 1.0, 0.01, 0.1, 15, 0, 0], 3100, time_start=3000, time_step=1
    )
    actual = actual.squeeze()
    assert (expected == actual).all()


def test_run_clc_pssacr(monkeypatch):
    monkeypatch.setenv("GSL_RNG_SEED", "01062024")
    samples = 10
    expected = []
    for i in range(samples):
        expected.append(np.loadtxt(f"tests/clc_pssacr_{i}.dat", dtype=int))
    expected = np.stack(expected).squeeze()

    system_size = 50
    reaction_rates = np.repeat(1, system_size)
    omega = 1
    initial_pops = np.repeat(1, system_size)

    pssa = m.pSSAlib(m.SSA.PSSACR)
    actual = pssa.sample_testcase_trajectory(
        m.Testcase.CyclicLinearChain,
        [system_size, *reaction_rates, omega, *initial_pops],
        1000,
        time_start=990,
        time_step=1,
        samples=samples,
    )

    actual = actual.squeeze()
    assert (expected == actual).all()


def test_homoreaction_pdf(monkeypatch):
    expected = np.loadtxt("tests/homoreaction_pdf.dat", dtype=float)

    omega = 1.0
    k1 = 0.016
    k2 = 10.0

    N = np.arange(0, 100, 1)
    analytical_pdf = m.homoreaction_pdf(N, k1, k2, omega)

    actual = np.stack((N, analytical_pdf), axis=1, dtype=float)
    assert np.isclose(actual, expected).all()
