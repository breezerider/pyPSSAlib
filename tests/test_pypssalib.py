import numpy as np
import pypssalib as m
import pytest


def test_version():
    assert m.__version__ == "0.1.0-dev4"


def test_method():
    pssa = m.pSSAlib(m.SSA.SPDM)
    assert pssa.method == m.SSA.SPDM

    pssa.method = m.SSA.PDM
    assert pssa.method == m.SSA.PDM


def test_run_invalid_testcase():
    pssa = m.pSSAlib(m.SSA.SPDM)
    with pytest.raises(RuntimeError) as excinfo:
        pssa.sample_testcase_population(m.Testcase.none, [1, 1, 10, 0], 100.0, 1)
    assert "Failed to intialize test case model" in str(excinfo.value)
    pssa.sample_testcase_population(m.Testcase.SingleBirthDeath, [1, 1, 10, 0], 100.0, 1)


def test_run_invalid_params():
    pssa = m.pSSAlib(m.SSA.SPDM)
    with pytest.raises(RuntimeError) as excinfo:
        pssa.sample_testcase_population(m.Testcase.SingleBirthDeath, None, 100.0, 1)
    assert "Parameters must be a float vector" in str(excinfo.value)
    pssa.sample_testcase_population(m.Testcase.SingleBirthDeath, [1, 1, 10, 0], 100.0, 1)


def test_run_invalid_time():
    pssa = m.pSSAlib(m.SSA.SPDM)
    # time_start
    with pytest.raises(RuntimeError) as excinfo:
        pssa.sample_testcase_trajectory(m.Testcase.SingleBirthDeath, [1, 1, 10, 0], 100.0, time_start=100.0)
    assert "Time interval must include at least one time point" in str(excinfo.value)
    pssa.sample_testcase_trajectory(m.Testcase.SingleBirthDeath, [1, 1, 10, 0], 100.0, time_start=99.0, time_step=1)
    # time_step
    with pytest.raises(RuntimeError) as excinfo:
        pssa.sample_testcase_trajectory(
            m.Testcase.SingleBirthDeath, [1, 1, 10, 0], 100.0, time_start=99.0, time_step=0.0
        )
    assert "Time interval must include at least one time point" in str(excinfo.value)
    pssa.sample_testcase_trajectory(m.Testcase.SingleBirthDeath, [1, 1, 10, 0], 100.0, time_start=99.0, time_step=1)


def test_run_invalid_samples():
    pssa = m.pSSAlib(m.SSA.SPDM)
    with pytest.raises(RuntimeError) as excinfo:
        pssa.sample_testcase_population(m.Testcase.SingleBirthDeath, [1, 1, 10, 0], 100.0, 0)
    assert "Number of samples must be a positive integer" in str(excinfo.value)
    pssa.sample_testcase_population(m.Testcase.SingleBirthDeath, [1, 1, 10, 0], 100.0, 1)


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


def test_homoreaction_pdf():
    expected = np.loadtxt("tests/homoreaction_pdf.dat", dtype=float)

    omega = 1.0
    k1 = 0.016
    k2 = 10.0

    N = np.arange(0, 100, 1)
    analytical_pdf = m.homoreaction_pdf(N, k1, k2, omega)

    actual = np.stack((N, analytical_pdf), axis=1, dtype=float)
    assert np.isclose(actual, expected).all()


class CATest:
    k1_gen = 1.0
    k11_asc = 2.1
    k11_dis = 0.1
    k1_deg = 0.01
    k2_deg = 0.1
    omega = 15
    S1_0 = 0
    S2_0 = 0

    # parameter vector
    @classmethod
    def k(cls, S1_0=0, S2_0=0):
        return [2, cls.k1_gen, cls.k11_asc, cls.k11_dis, cls.k1_deg, cls.k2_deg, cls.omega, S1_0, S2_0]

    # analytical ODEs
    @classmethod
    def dsdt(cls, s):
        return np.array(
            [
                cls.k1_gen - 2.0 * s[0] * s[0] * cls.k11_asc + 2.0 * s[1] * cls.k11_dis - cls.k1_deg * s[0],
                s[0] * s[0] * cls.k11_asc - s[1] * cls.k11_dis - cls.k2_deg * s[1],
            ]
        )

    # analytical Jacobian
    @classmethod
    def jacobian(cls, s):
        return np.array(
            [
                [-4.0 * s[0] * cls.k11_asc - cls.k1_deg, 2.0 * cls.k11_dis],
                [2.0 * s[0] * cls.k11_asc, -cls.k11_dis - cls.k2_deg],
            ]
        )

    # analytical Jacobian
    @classmethod
    def lyapunovQ(cls, s):
        # stoichiometry matrix
        S = np.array(
            [
                [1, 0],
                [-2, 1],
                [2, -1],
                [-1, 0],
                [0, -1],
            ]
        )
        F = np.array([cls.k1_gen, cls.k11_asc * s[0] ** 2, cls.k11_dis * s[1], cls.k1_deg * s[0], cls.k2_deg * s[1]])

        return cls.omega * S.T.dot(np.diag(F)).dot(S)


def test_odes():
    # create simulator instance
    pssa = m.pSSAlib()

    nx, ny = 5, 5
    sx, sy = np.linspace(10, 100, nx, dtype=int), np.linspace(10, 100, ny, dtype=int)
    SX, SY = np.meshgrid(sx, sy, indexing='ij')
    for i in range(nx):
        for j in range(ny):
            # number of molecules
            S = np.array([SX[i, j], SY[i, j]])
            # species concentration
            s = S / CATest.omega

            expected = CATest.dsdt(s)

            actual = pssa.odes(m.Testcase.ca, CATest.k(), S)

            assert np.isclose(actual, expected).all()


def test_jacobian():
    # create simulator instance
    pssa = m.pSSAlib()

    nx, ny = 5, 5
    sx, sy = np.linspace(10, 100, nx, dtype=int), np.linspace(10, 100, ny, dtype=int)
    SX, SY = np.meshgrid(sx, sy, indexing='ij')
    for i in range(nx):
        for j in range(ny):
            # number of molecules
            S = np.array([SX[i, j], SY[i, j]])
            # species concentration
            s = S / CATest.omega

            expected = CATest.jacobian(s)

            actual = pssa.jacobian(m.Testcase.ca, CATest.k(), S)

            assert np.isclose(actual, expected).all()


def test_lyapunovQ():
    # create simulator instance
    pssa = m.pSSAlib()

    nx, ny = 5, 5
    sx, sy = np.linspace(10, 100, nx, dtype=int), np.linspace(10, 100, ny, dtype=int)
    SX, SY = np.meshgrid(sx, sy, indexing='ij')
    for i in range(nx):
        for j in range(ny):
            # number of molecules
            S = np.array([SX[i, j], SY[i, j]])
            # species concentration
            s = S / CATest.omega

            expected = CATest.lyapunovQ(s)

            actual = pssa.lyapunovQ(m.Testcase.ca, CATest.k(), S)

            print(expected)
            print(actual)

            assert np.isclose(actual, expected).all()


def test_print_ca():
    # create simulator instance
    pssa = m.pSSAlib()
    expected = """'[ColloidalAggregation]' model:

Compartment volume: 15

Initial species population:
[S0] = 1\t[S1] = 2\t

Reaction network:
[R0]: 0 * [ ] -1-> 1 * [S0]
[R1]: 2 * [S0] -2.1-> 1 * [S1]
[R2]: 1 * [S1] -0.1-> 2 * [S0]
[R3]: 1 * [S0] -0.01-> 0 * [ ]
[R4]: 1 * [S1] -0.1-> 0 * [ ]
"""

    actual = pssa.print_test_case(m.Testcase.ca, CATest.k(S1_0=1, S2_0=2))

    assert expected == actual


def test_print_ed():
    # create simulator instance
    pssa = m.pSSAlib()
    expected = """'[EnzymaticDegration]' model:

Compartment volume: 15

Initial species population:
[mRNA] = 1\t[Complex] = 2\t[Enzyme] = 3\t[Protein] = 4\t

Reaction network:
[mRNA_Production]: 0 * [ ] -2.25-> 1 * [mRNA]
[Protein_Translation]: 1 * [mRNA] -4-> 1 * [mRNA] + 1 * [Protein]
[Complex_Formation]: 1 * [Enzyme] + 1 * [Protein] -100-> 1 * [Complex]
[Complex_Detachment]: 1 * [Complex] -1-> 1 * [Enzyme] + 1 * [Protein]
[Complex_Degraded]: 1 * [Complex] -1-> 1 * [Enzyme]
[mRNA_Degradation]: 1 * [mRNA] -10-> 0 * [ ]
"""

    actual = pssa.print_test_case(m.Testcase.ed, [2.25, 4, 100, 1, 1, 10, 15, 1, 2, 3, 4])

    assert expected == actual
