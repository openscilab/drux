"""Tests for the Higuchi model implementation in drux package."""

from pytest import raises
from numpy import isclose
from math import sqrt
from re import escape
from drux import HiguchiModel

TEST_CASE_NAME = "Higuchi model tests"
D, C0, CS = 1e-6, 1, 0.5
SIM_DURATION, SIM_TIME_STEP = 1000, 10
RELATIVE_TOLERANCE = 1e-2


def test_higuchi_parameters():
    model = HiguchiModel(D=D, c0=C0, cs=CS)
    assert model._parameters.D == D
    assert model._parameters.c0 == C0
    assert model._parameters.cs == CS


def test_invalid_parameters():
    with raises(ValueError, match=escape("Diffusivity (D) must be positive.")):
        HiguchiModel(D=-D, c0=C0, cs=CS).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Initial drug concentration (c0) must be positive.")):
        HiguchiModel(D=D, c0=-C0, cs=CS).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Solubility (cs) must be positive.")):
        HiguchiModel(D=D, c0=C0, cs=-CS).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Solubility (cs) must be lower or equal to initial concentration (c0).")):
        HiguchiModel(D=D, c0=0.5, cs=1).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)


def test_higuchi_simulation():  # Reference: https://www.sciencedirect.com/science/article/abs/pii/S0022354915333037
    model = HiguchiModel(D=D, c0=C0, cs=CS)
    profile = model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    actual_release = [sqrt(D * t * (2 * C0 - CS) * CS) for t in range(0, 1001, 10)]
    assert all(isclose(p, a, rtol=RELATIVE_TOLERANCE) for p, a in zip(profile, actual_release))


def test_higuchi_simulation_errors():
    model = HiguchiModel(D=D, c0=C0, cs=CS)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=-100, time_step=10)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=100, time_step=-10)

    with raises(ValueError, match="Time step cannot be greater than duration"):
        model.simulate(duration=10, time_step=20)


def test_higuchi_plot1():
    model = HiguchiModel(D=D, c0=C0, cs=CS)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    fig, ax = model.plot()
    assert fig is not None
    assert ax is not None
    assert ax.get_title() == model._plot_parameters["title"]
    assert ax.get_xlabel() == model._plot_parameters["xlabel"]
    assert ax.get_ylabel() == model._plot_parameters["ylabel"]
    assert [text.get_text() for text in ax.get_legend().get_texts()] == [model._plot_parameters["label"]]


def test_higuchi_plot2():
    model = HiguchiModel(D=D, c0=C0, cs=CS)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    fig, ax = model.plot(title="test-title", xlabel="test-xlabel", ylabel="test-ylabel", label="test-label")
    assert fig is not None
    assert ax is not None
    assert ax.get_title() == "test-title"
    assert ax.get_xlabel() == "test-xlabel"
    assert ax.get_ylabel() == "test-ylabel"
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["test-label"]


def test_higuchi_plot_error():
    model = HiguchiModel(D=D, c0=C0, cs=CS)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.plot()

    model._time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model._release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.plot()


def test_higuchi_release_rate():  # Reference: https://www.wolframalpha.com/input?i=get+the+derivative+of+sqrt%28D*C_s*%282*C_0-C_s%29*t%29+with+respect+to+t
    model = HiguchiModel(D=D, c0=C0, cs=CS)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    rate = model.get_release_rate()
    actual_rate = [sqrt(D * t * (2 * C0 - CS) * CS) / (2 * t)
                   for t in range(1, 1001, 10)]  # not defined at t=0, set to 0
    # skip first point to avoid near zero division issues
    assert all(isclose(r, a, rtol=1e-2) for r, a in zip(rate[10:], actual_rate[10:]))


def test_higuchi_release_rate_error():
    model = HiguchiModel(D=D, c0=C0, cs=CS)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.get_release_rate()

    model._time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model._release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.get_release_rate()


def test_repr():
    model = HiguchiModel(D=D, c0=C0, cs=CS)
    repr_str = repr(model)
    assert repr_str == f"drux.HiguchiModel(D={D}, c0={C0}, cs={CS})"


def test_higuchi_time_for_release():  # Reference: https://www.wolframalpha.com/input?i=solve+for+t+in+sqrt%2810%5E%28-6%29*0.5*%282*1.5-0.5%29*t%29+%3D+0.5*sqrt%2810%5E%28-6%29*0.5*%282*1.5-0.5%29*1000%29
    model = HiguchiModel(D=D, c0=C0, cs=CS)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    tx = model.time_for_release(0.5 * model._release_profile[-1])
    assert isclose(tx, 250.0, rtol=1e-2)


def test_higuchi_time_for_release_error():
    model = HiguchiModel(D=D, c0=C0, cs=CS)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.time_for_release(0.5)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match="Target release must be non-negative."):
        model.time_for_release(-0.1)

    with raises(ValueError, match="Target release exceeds maximum release of the simulated duration."):
        model.time_for_release(min(model._release_profile[-1] + 0.1, 1))
