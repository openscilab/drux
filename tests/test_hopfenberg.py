"""Tests for the Hopfenberg model implementation in drux package."""

from pytest import raises
from numpy import isclose
from re import escape
from drux import HopfenbergModel

TEST_CASE_NAME = "Hopfenberg model tests"
M, k0, c0, a0, n = 1, 0.00067, 0.0374, 3.51, 2
SIM_DURATION, SIM_TIME_STEP = 100, 1
RELATIVE_TOLERANCE = 1e-1


def test_hopfenberg_parameters():
    model = HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=n)
    assert model._parameters.M == M
    assert model._parameters.k0 == k0
    assert model._parameters.c0 == c0
    assert model._parameters.a0 == a0
    assert model._parameters.n == n


def test_invalid_parameters():
    with raises(ValueError, match=escape("Entire releasable amount of drug (M) must be non-negative.")):
        HopfenbergModel(M=-M, k0=k0, c0=c0, a0=a0, n=n).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Erosion rate constant (k0) must be non-negative")):
        HopfenbergModel(M=M, k0=-k0, c0=c0, a0=a0, n=n).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Initial drug concentration (c0) must be positive.")):
        HopfenbergModel(M=M, k0=k0, c0=-c0, a0=a0, n=n).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Initial radius or half-thickness (a0) must be positive.")):
        HopfenbergModel(M=M, k0=k0, c0=c0, a0=-a0, n=n).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match=escape("Geometry factor (n) must be 1 (slab), 2 (cylinder), or 3 (sphere).")):
        HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=4).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)


def test_repr():
    model = HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=n)
    repr_str = repr(model)
    assert repr_str == f"drux.HopfenbergModel(M={M}, k0={k0}, c0={c0}, a0={a0}, n={n})"


def test_hopfenberg_simulation():  # Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC3500559/
    model = HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=n)
    profile = model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    actual_release = [M * (1 - (1 - (k0 * t) / (c0 * a0))**n) for t in range(0, SIM_DURATION+SIM_TIME_STEP, SIM_TIME_STEP)]
    assert all(isclose(p, r, rtol=RELATIVE_TOLERANCE) for p, r in zip(profile, actual_release))


def test_hopfenberg_simulation_errors():
    model = HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=n)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=-100, time_step=10)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=100, time_step=-10)

    with raises(ValueError, match="Time step cannot be greater than duration"):
        model.simulate(duration=10, time_step=20)


def test_hopfenberg_plot1():
    model = HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=n)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    fig, ax = model.plot()
    assert fig is not None
    assert ax is not None
    assert ax.get_title() == model._plot_parameters["title"]
    assert ax.get_xlabel() == model._plot_parameters["xlabel"]
    assert ax.get_ylabel() == model._plot_parameters["ylabel"]
    assert [text.get_text() for text in ax.get_legend().get_texts()] == [model._plot_parameters["label"]]


def test_hopfenberg_plot2():
    model = HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=n)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    fig, ax = model.plot(title="test-title", xlabel="test-xlabel", ylabel="test-ylabel", label="test-label")
    assert fig is not None
    assert ax is not None
    assert ax.get_title() == "test-title"
    assert ax.get_xlabel() == "test-xlabel"
    assert ax.get_ylabel() == "test-ylabel"
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["test-label"]


def test_hopfenberg_plot_error():
    model = HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=n)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.plot()

    model._time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model._release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.plot()


def test_hopfenberg_release_rate():
    small_timestep = 0.001  # smaller time step for better accuracy in numerical derivative
    small_duration = 10
    model = HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=n)
    model.simulate(duration=small_duration, time_step=small_timestep)
    rate = model.get_release_rate().tolist()
    import numpy as np
    # Ref: https://www.wolframalpha.com/input?i=get+the+derivative+of+%28M+*+%281+-+%281+-+%28k0+*+t%29+%2F+%28c0+*+a0%29%29**n%29%29+with+respect+to+t
    actual_rate = [
        k0*M*n*((1-((k0*t)/(a0*c0)))**(n-1))/(a0*c0)
        for t in np.arange(small_timestep, small_duration + small_timestep, small_timestep)
    ]  # avoid t=0 to prevent division by zero
    assert all(isclose(r, ar, rtol=1e-1) for r, ar in zip(rate[2:], actual_rate[1:]))  # skip first two points due to numerical derivative inaccuracies


def test_hopfenberg_release_rate_error():
    model = HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=n)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.get_release_rate()

    model._time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model._release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.get_release_rate()


def test_hopfenberg_time_for_release():
    model = HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=n)
    model.simulate(duration=SIM_DURATION, time_step=1)
    tx = model.time_for_release(0.8 * model._release_profile[-1])
    # Ref: https://www.wolframalpha.com/input?i=solve+for+t+in+%281+-+%281+-+%280.00067*+t%29+%2F+%280.0374+*+3.51%29%29**2%29+%3D+0.8*0.7602750594732341
    assert isclose(tx, 74, rtol=1e-2)


def test_hopfenberg_time_for_release_error():
    model = HopfenbergModel(M=M, k0=k0, c0=c0, a0=a0, n=n)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.time_for_release(0.5)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match="Target release must be non-negative."):
        model.time_for_release(-0.1)

    with raises(ValueError, match="Target release exceeds maximum release of the simulated duration."):
        model.time_for_release(model._release_profile[-1] + 0.1)
