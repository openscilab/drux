"""Tests for the Higuchi model implementation in drux package."""

from pytest import raises
from unittest import mock
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
    assert model.params.D == D
    assert model.params.c0 == C0
    assert model.params.cs == CS


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


@mock.patch("matplotlib.pyplot.subplots")
def test_higuchi_plot(mock_subplots: mock.MagicMock):
    mock_subplots.return_value = (mock.MagicMock(), mock.MagicMock())
    model = HiguchiModel(D=D, c0=C0, cs=CS)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    fig, ax = model.plot()
    assert fig is not None
    assert ax is not None
    mock_subplots.assert_called_once()


def test_higuchi_plot_error():
    model = HiguchiModel(D=D, c0=C0, cs=CS)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.plot()

    model.time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model.release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.plot()

    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    with mock.patch.dict('sys.modules', {'matplotlib': None}):
        with raises(ImportError, match="Matplotlib is required for plotting but not installed."):
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

    model.time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model.release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.get_release_rate()


def test_higuchi_time_for_release():  # Reference: https://www.wolframalpha.com/input?i=solve+for+t+in+sqrt%2810%5E%28-6%29*0.5*%282*1.5-0.5%29*t%29+%3D+0.5*sqrt%2810%5E%28-6%29*0.5*%282*1.5-0.5%29*1000%29
    model = HiguchiModel(D=D, c0=C0, cs=CS)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    tx = model.time_for_release(0.5 * model.release_profile[-1])
    assert isclose(tx, 250.0, rtol=1e-2)


def test_higuchi_time_for_release_error():
    model = HiguchiModel(D=D, c0=C0, cs=CS)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.time_for_release(0.5)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match="Target release must be non-negative."):
        model.time_for_release(-0.1)

    with raises(ValueError, match="Target release exceeds maximum release of the simulated duration."):
        model.time_for_release(min(model.release_profile[-1] + 0.1, 1))
