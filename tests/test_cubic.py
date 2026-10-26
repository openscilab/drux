"""Tests for the cubic model implementation in drux package."""

from pytest import raises
from numpy import isclose
from re import escape
from drux import CubicModel

TEST_CASE_NAME = "Cubic model tests"
M0, k1, k2, k3 = 0.01, 0.1, 0.00005, 0.00000005
SIM_DURATION, SIM_TIME_STEP = 1000, 10
RELATIVE_TOLERANCE = 1e-2


def test_cubic_parameters():
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=k3)
    assert model._parameters.M0 == M0
    assert model._parameters.k1 == k1
    assert model._parameters.k2 == k2
    assert model._parameters.k3 == k3


def test_cubic_default_initial_amount():
    model = CubicModel(k1=k1, k2=k2, k3=k3)
    assert model._parameters.M0 == 0


def test_invalid_parameters():
    with raises(ValueError, match=escape("Initial amount of drug in the solution (M0) must be non-negative.")):
        CubicModel(M0=-M0, k1=k1, k2=k2, k3=k3).simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)


def test_cubic_negative_coefficients():  # k1, k2 and k3 are not restricted, e.g. k3 < 0 gives a concave end of profile
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=-k3)
    profile = model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    actual_release = [M0 + (k1 * t) + (k2 * t**2) - (k3 * t**3)
                      for t in range(0, SIM_DURATION + SIM_TIME_STEP, SIM_TIME_STEP)]
    assert all(isclose(p, a, rtol=RELATIVE_TOLERANCE) for p, a in zip(profile, actual_release))

    model = CubicModel(M0=M0, k1=k1, k2=-k2, k3=k3)
    profile = model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    actual_release = [M0 + (k1 * t) - (k2 * t**2) + (k3 * t**3)
                      for t in range(0, SIM_DURATION + SIM_TIME_STEP, SIM_TIME_STEP)]
    assert all(isclose(p, a, rtol=RELATIVE_TOLERANCE) for p, a in zip(profile, actual_release))


def test_repr():
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=k3)
    repr_str = repr(model)
    assert repr_str == f"drux.CubicModel(k1={k1}, k2={k2}, k3={k3}, M0={M0})"


def test_cubic_simulation():  # Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC8912140/
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=k3)
    profile = model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    actual_release = [M0 + (k1 * t) + (k2 * t**2) + (k3 * t**3)
                      for t in range(0, SIM_DURATION + SIM_TIME_STEP, SIM_TIME_STEP)]
    assert all(isclose(p, a, rtol=RELATIVE_TOLERANCE) for p, a in zip(profile, actual_release))


def test_cubic_simulation_errors():
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=k3)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=-100, time_step=10)

    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=100, time_step=-10)

    with raises(ValueError, match="Time step cannot be greater than duration"):
        model.simulate(duration=10, time_step=20)


def test_cubic_plot1():
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=k3)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    fig, ax = model.plot()
    assert fig is not None
    assert ax is not None
    assert ax.get_title() == model._plot_parameters["title"]
    assert ax.get_xlabel() == model._plot_parameters["xlabel"]
    assert ax.get_ylabel() == model._plot_parameters["ylabel"]
    assert [text.get_text() for text in ax.get_legend().get_texts()] == [model._plot_parameters["label"]]


def test_cubic_plot2():
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=k3)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    fig, ax = model.plot(title="test-title", xlabel="test-xlabel", ylabel="test-ylabel", label="test-label")
    assert fig is not None
    assert ax is not None
    assert ax.get_title() == "test-title"
    assert ax.get_xlabel() == "test-xlabel"
    assert ax.get_ylabel() == "test-ylabel"
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["test-label"]


def test_cubic_plot_error():
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=k3)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.plot()

    model._time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model._release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.plot()


def test_cubic_release_rate():  # Reference: https://www.wolframalpha.com/input?i=get+the+derivative+of+M0+%2B+k1+*+t+%2B+k2+*+t%5E2+%2B+k3+*+t%5E3+with+respect+to+t
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=k3)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)
    rate = model.get_release_rate()
    actual_rate = [k1 + (2 * k2 * t) + (3 * k3 * t**2) for t in range(0, SIM_DURATION + SIM_TIME_STEP, SIM_TIME_STEP)]
    assert all(isclose(r, a, rtol=RELATIVE_TOLERANCE) for r, a in zip(rate, actual_rate))


def test_cubic_release_rate_error():
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=k3)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.get_release_rate()

    model._time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    model._release_profile = [0.0]
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.get_release_rate()


def test_cubic_time_for_release():  # Reference: https://www.wolframalpha.com/input?i=solve+for+t+in+0.01+%2B+0.1*t+%2B+0.00005*t%5E2+%2B+0.00000005*t%5E3+%3D+0.5*%280.01+%2B+0.1*1000+%2B+0.00005*1000%5E2+%2B+0.00000005*1000%5E3%29
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=k3)
    model.simulate(duration=SIM_DURATION, time_step=1)
    tx = model.time_for_release(0.5 * model._release_profile[-1])
    assert isclose(tx, 650.6, rtol=1e-2)


def test_cubic_time_for_release_error():
    model = CubicModel(M0=M0, k1=k1, k2=k2, k3=k3)

    with raises(ValueError, match=escape("No simulation data available. Run simulate() first.")):
        model.time_for_release(0.5)
    model.simulate(duration=SIM_DURATION, time_step=SIM_TIME_STEP)

    with raises(ValueError, match="Target release must be non-negative."):
        model.time_for_release(-0.1)

    with raises(ValueError, match="Target release exceeds maximum release of the simulated duration."):
        model.time_for_release(model._release_profile[-1] + 0.1)
