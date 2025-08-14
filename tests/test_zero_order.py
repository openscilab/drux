"""Tests for the Zero-order model implementation in drux package."""

from pytest import raises
from unittest import mock
from numpy import isclose
from math import sqrt
from re import escape
from drux import ZeroOrderModel

TEST_CASE_NAME = "Zero-order model tests"
M0, k0 = 0.01, 0.1
SIM_DURATION, SIM_TIME_STEP = 1000, 10
RELATIVE_TOLERANCE = 1e-2


def test_zer_order_parameters():
    model = ZeroOrderModel(M0=M0, k0=k0)
    assert model.params.M0 == M0
    assert model.params.k0 == k0

def test_invalid_parameters():
    with raises(ValueError, match=escape("Initial amount of drug in the solution (M0) must be non-negative.")):
        ZeroOrderModel(M0=-M0, k0=k0).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Release rate (k0) must be non-negative.")):
        ZeroOrderModel(M0=M0, k0=-k0).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)


def test_zero_order_simulation():  # Reference: https://europepmc.org/article/pmc/3425064
    model = ZeroOrderModel(M0=M0, k0=k0)
    profile = model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    actual_release = [M0 + (k0 * t) for t in range(0, 1001, 10)]
    assert all(isclose(p, a, rtol=RELATIVE_TOLERANCE) for p, a in zip(profile, actual_release))


def test_zero_order_simulation_errors():
    model = ZeroOrderModel(M0=M0, k0=k0)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=-100, time_step=10)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=100, time_step=-10)

    with raises(ValueError, match="Time step cannot be greater than duration"):
        model.simulate(duration=10, time_step=20)


@mock.patch("matplotlib.pyplot.subplots")
def test_zero_order_plot(mock_subplots: mock.MagicMock):
    mock_subplots.return_value = (mock.MagicMock(), mock.MagicMock())
    model = ZeroOrderModel(M0=M0, k0=k0)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    fig, ax = model.plot()
    assert fig is not None
    assert ax is not None
    mock_subplots.assert_called_once()


def test_zero_order_plot_error():
    model = ZeroOrderModel(M0=M0, k0=k0)

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


def test_zero_order_release_rate():  # Reference: https://www.wolframalpha.com/input?i=get+the+derivative+of+M0+%2B+%28k0+*+t%29+with+respect+to+t
    model = ZeroOrderModel(M0=M0, k0=k0)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    rate = model.get_release_rate()
    actual_rate = [k0 for _ in range(0, 1001, 10)]
    assert all(isclose(r, a, rtol=1e-2) for r, a in zip(rate, actual_rate))


def test_zero_order_release_rate_error():
    model = ZeroOrderModel(M0=M0, k0=k0)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.get_release_rate()

    model.time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model.release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.get_release_rate()


def test_zero_order_time_for_release():  # Reference: https://www.wolframalpha.com/input?i=solve+for+t+in+0.01+%2B+%280.1*t%29+%3D+0.5*+%280.01+%2B+%280.1*1000%29%29
    model = ZeroOrderModel(M0=M0, k0=k0)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    tx = model.time_for_release(0.5 * model.release_profile[-1])
    assert isclose(tx, 499.95, rtol=1e-2)


def test_zero_order_time_for_release_error():
    model = ZeroOrderModel(M0=M0, k0=k0)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.time_for_release(0.5)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match="Target release must be non-negative."):
        model.time_for_release(-0.1)

    with raises(ValueError, match="Target release exceeds maximum release of the simulated duration."):
        model.time_for_release(model.release_profile[-1] + 0.1)
