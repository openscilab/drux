"""Tests for the Weibull model implementation in drux package."""

from pytest import raises
from numpy import isclose
from re import escape
from math import exp
from drux import WeibullModel

TEST_CASE_NAME = "Weibull model tests"
M, a, b = 1, 0.095, 0.7
SIM_DURATION, SIM_TIME_STEP = 100, 1
RELATIVE_TOLERANCE = 1e-1


def test_weibull_parameters():
    model = WeibullModel(M=M, a=a, b=b)
    assert model._parameters.M == M
    assert model._parameters.a == a
    assert model._parameters.b == b


def test_invalid_parameters():
    with raises(ValueError, match=escape("Entire releasable amount of drug (M) must be non-negative.")):
        WeibullModel(M=-M, a=a, b=b).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Scale parameter (a) must be positive.")):
        WeibullModel(M=M, a=-a, b=b).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Shape parameter (b) must be positive.")):
        WeibullModel(M=M, a=a, b=-b).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)


def test_repr():
    model = WeibullModel(M=M, a=a, b=b)
    repr_str = repr(model)
    assert repr_str == f"drux.WeibullModel(M={M}, a={a}, b={b})"


def test_weibull_simulation():  # Reference: https://www.mdpi.com/2073-4360/13/17/2897
    model = WeibullModel(M=M, a=a, b=b)
    profile = model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    actual_release = [M * (1 - exp(-a * t ** b)) for t in range(0, SIM_DURATION+SIM_TIME_STEP, SIM_TIME_STEP)]
    assert all(isclose(p, r, rtol=RELATIVE_TOLERANCE) for p, r in zip(profile, actual_release))


def test_weibull_simulation_errors():
    model = WeibullModel(M=M, a=a, b=b)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=-100, time_step=10)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=100, time_step=-10)

    with raises(ValueError, match="Time step cannot be greater than duration"):
        model.simulate(duration=10, time_step=20)


def test_weibull_plot1():
    model = WeibullModel(M=M, a=a, b=b)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    fig, ax = model.plot()
    assert fig is not None
    assert ax is not None
    assert ax.get_title() == model._plot_parameters["title"]
    assert ax.get_xlabel() == model._plot_parameters["xlabel"]
    assert ax.get_ylabel() == model._plot_parameters["ylabel"]
    assert [text.get_text() for text in ax.get_legend().get_texts()] == [model._plot_parameters["label"]]


def test_weibull_plot2():
    model = WeibullModel(M=M, a=a, b=b)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    fig, ax = model.plot(title="test-title", xlabel="test-xlabel", ylabel="test-ylabel", label="test-label")
    assert fig is not None
    assert ax is not None
    assert ax.get_title() == "test-title"
    assert ax.get_xlabel() == "test-xlabel"
    assert ax.get_ylabel() == "test-ylabel"
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["test-label"]


def test_weibull_plot_error():
    model = WeibullModel(M=M, a=a, b=b)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.plot()

    model._time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model._release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.plot()


def test_weibull_release_rate():  # Reference: https://www.wolframalpha.com/input?i=get+the+derivative+of+%28M+*+%281+-+exp%28-+%28%28t+%2F+r%29+**+beta%29%29%29%29+with+respect+to+t
    small_timestep = 0.001  # smaller time step for better accuracy in numerical derivative
    small_duration = 10
    model = WeibullModel(M=M, a=a, b=b)
    model.simulate(duration=small_duration, time_step=small_timestep)
    rate = model.get_release_rate().tolist()
    import numpy as np
    actual_rate = [
        (M*b*a) * t**(b - 1) * exp(-a * t**b)
        for t in np.arange(small_timestep, small_duration + small_timestep, small_timestep)
    ]  # avoid t=0 to prevent division by zero
    assert all(isclose(r, ar, rtol=1e-1) for r, ar in zip(rate[2:], actual_rate[1:]))  # skip first two points due to numerical derivative inaccuracies


def test_weibull_release_rate_error():
    model = WeibullModel(M=M, a=a, b=b)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.get_release_rate()

    model._time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model._release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.get_release_rate()


def test_weibull_time_for_release():  # Reference: https://www.wolframalpha.com/input?i=solve+for+t+in+%281+-+exp%28-0.095+*+t+**+0.7%29%29+%3D+0.8*0.908
    model = WeibullModel(M=M, a=a, b=b)
    model.simulate(duration=SIM_DURATION, time_step=1)
    tx = model.time_for_release(0.8 * model._release_profile[-1])
    assert isclose(tx, 42, rtol=1e-2)


def test_weibull_time_for_release_error():
    model = WeibullModel(M=M, a=a, b=b)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.time_for_release(0.5)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match="Target release must be non-negative."):
        model.time_for_release(-0.1)

    with raises(ValueError, match="Target release exceeds maximum release of the simulated duration."):
        model.time_for_release(model._release_profile[-1] + 0.1)
