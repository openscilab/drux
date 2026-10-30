"""Tests for the Korsmeyer-Peppas model implementation in drux package."""

from pytest import raises
from numpy import isclose
from re import escape
from drux import KorsmeyerPeppasModel

TEST_CASE_NAME = "Korsmeyer-Peppas model tests"
M0, k, n = 0.01, 0.1, 0.6
SIM_DURATION, SIM_TIME_STEP = 1000, 10
RELATIVE_TOLERANCE = 1e-2


def test_korsmeyer_peppas_parameters():
    model = KorsmeyerPeppasModel(M0=M0, k=k, n=n)
    assert model._parameters.M0 == M0
    assert model._parameters.k == k
    assert model._parameters.n == n


def test_korsmeyer_peppas_default_initial_amount():
    model = KorsmeyerPeppasModel(k=k, n=n)
    assert model._parameters.M0 == 0


def test_invalid_parameters():
    with raises(ValueError, match=escape("Initial amount of drug in the solution (M0) must be non-negative.")):
        KorsmeyerPeppasModel(M0=-M0, k=k, n=n).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Release rate constant (k) must be non-negative.")):
        KorsmeyerPeppasModel(M0=M0, k=-k, n=n).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Release exponent (n) must be positive.")):
        KorsmeyerPeppasModel(M0=M0, k=k, n=-n).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Release exponent (n) must be positive.")):
        KorsmeyerPeppasModel(M0=M0, k=k, n=0).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)


def test_repr():
    model = KorsmeyerPeppasModel(M0=M0, k=k, n=n)
    repr_str = repr(model)
    assert repr_str == f"drux.KorsmeyerPeppasModel(k={k}, n={n}, M0={M0})"


def test_korsmeyer_peppas_simulation():  # Reference: https://www.ptfarm.pl/pub/File/Acta_Poloniae/2010/3/217.pdf
    model = KorsmeyerPeppasModel(M0=M0, k=k, n=n)
    profile = model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    actual_release = [M0 + (k * t**n) for t in range(0, SIM_DURATION + SIM_TIME_STEP, SIM_TIME_STEP)]
    assert all(isclose(p, a, rtol=RELATIVE_TOLERANCE) for p, a in zip(profile, actual_release))


def test_korsmeyer_peppas_simulation_errors():
    model = KorsmeyerPeppasModel(M0=M0, k=k, n=n)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=-100, time_step=10)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=100, time_step=-10)

    with raises(ValueError, match="Time step cannot be greater than duration"):
        model.simulate(duration=10, time_step=20)


def test_korsmeyer_peppas_plot1():
    model = KorsmeyerPeppasModel(M0=M0, k=k, n=n)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    fig, ax = model.plot()
    assert fig is not None
    assert ax is not None
    assert ax.get_title() == model._plot_parameters["title"]
    assert ax.get_xlabel() == model._plot_parameters["xlabel"]
    assert ax.get_ylabel() == model._plot_parameters["ylabel"]
    assert [text.get_text() for text in ax.get_legend().get_texts()] == [model._plot_parameters["label"]]


def test_korsmeyer_peppas_plot2():
    model = KorsmeyerPeppasModel(M0=M0, k=k, n=n)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    fig, ax = model.plot(title="test-title", xlabel="test-xlabel", ylabel="test-ylabel", label="test-label")
    assert fig is not None
    assert ax is not None
    assert ax.get_title() == "test-title"
    assert ax.get_xlabel() == "test-xlabel"
    assert ax.get_ylabel() == "test-ylabel"
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["test-label"]


def test_korsmeyer_peppas_plot_error():
    model = KorsmeyerPeppasModel(M0=M0, k=k, n=n)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.plot()

    model._time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model._release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.plot()


def test_korsmeyer_peppas_release_rate():  # Reference: https://www.wolframalpha.com/input?i=get+the+derivative+of+M0+%2B+k+*+t%5En+with+respect+to+t
    model = KorsmeyerPeppasModel(M0=M0, k=k, n=n)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    rate = model.get_release_rate()
    actual_rate = [k * n * t**(n - 1)
                   for t in range(SIM_TIME_STEP, SIM_DURATION + SIM_TIME_STEP, SIM_TIME_STEP)]  # not defined at t=0
    # skip the first points to avoid the numerical derivative inaccuracies near t=0
    assert all(isclose(r, a, rtol=RELATIVE_TOLERANCE) for r, a in zip(rate[10:], actual_rate[9:]))


def test_korsmeyer_peppas_release_rate_error():
    model = KorsmeyerPeppasModel(M0=M0, k=k, n=n)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.get_release_rate()

    model._time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model._release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.get_release_rate()


def test_korsmeyer_peppas_time_for_release():  # Reference: https://www.wolframalpha.com/input?i=solve+for+t+in+0.01+%2B+0.1*t%5E0.6+%3D+0.5*%280.01+%2B+0.1*1000%5E0.6%29
    model = KorsmeyerPeppasModel(M0=M0, k=k, n=n)
    model.simulate(duration=SIM_DURATION, time_step=1)
    tx = model.time_for_release(0.5 * model._release_profile[-1])
    assert isclose(tx, 314.1, rtol=1e-2)


def test_korsmeyer_peppas_time_for_release_error():
    model = KorsmeyerPeppasModel(M0=M0, k=k, n=n)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.time_for_release(0.5)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match="Target release must be non-negative."):
        model.time_for_release(-0.1)

    with raises(ValueError, match="Target release exceeds maximum release of the simulated duration."):
        model.time_for_release(model._release_profile[-1] + 0.1)
