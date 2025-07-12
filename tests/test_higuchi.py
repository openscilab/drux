"""Tests for the Higuchi model implementation in drux package."""

from pytest import raises
from unittest import mock
from numpy import isclose
from math import sqrt
from drux import HiguchiModel

TEST_CASE_NAME = "Higuchi model tests"


def test_higuchi_parameters():
    model = HiguchiModel(D=1e-6, c0=1.5, cs=0.5)
    assert model.params.D == 1e-6
    assert model.params.c0 == 1.5
    assert model.params.cs == 0.5


def test_invalid_parameters():
    with raises(ValueError, match="Diffusivity \\(D\\) must be positive."):
        HiguchiModel(D=-1e-6, c0=1.5, cs=0.5).simulate(duration=1000, time_step=10)
    
    with raises(ValueError, match="Initial drug concentration \\(c0\\) must be positive."):
        HiguchiModel(D=1e-6, c0=-1.5, cs=0.5).simulate(duration=1000, time_step=10)
    
    with raises(ValueError, match="Solubility \\(cs\\) must be positive."):
        HiguchiModel(D=1e-6, c0=1.5, cs=-0.5).simulate(duration=1000, time_step=10)
    
    with raises(ValueError, match="Solubility \\(cs\\) must be lower or equal to initial concentration \\(c0\\)."):
        HiguchiModel(D=1e-6, c0=0.5, cs=1.5).simulate(duration=1000, time_step=10)


def test_higuchi_simulation(): # Reference: https://www.sciencedirect.com/science/article/abs/pii/S0022354915333037
    D, A, Cs = 1e-6, 1.5, 0.5
    model = HiguchiModel(D=D, c0=A, cs=Cs)
    profile = model.simulate(duration=1000, time_step=10)
    actual_release = [sqrt(D * t * (2 * A - Cs) * Cs) for t in range(0, 1001, 10)]
    assert all(isclose(p, a, rtol=1e-3) for p, a in zip(profile, actual_release))


def test_higuchi_simulation_errors():
    model = HiguchiModel(D=1e-6, c0=1.5, cs=0.5)
    
    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=-100, time_step=10)
    
    with raises(ValueError, match="Duration and time step must be positive values"):
        model.simulate(duration=100, time_step=-10)
    
    with raises(ValueError, match="Time step cannot be greater than duration"):
        model.simulate(duration=10, time_step=20)


@mock.patch("matplotlib.pyplot.subplots")
def test_higuchi_plot(mock_subplots: mock.MagicMock):
    mock_subplots.return_value = (mock.MagicMock(), mock.MagicMock())
    model = HiguchiModel(D=1e-6, c0=1.5, cs=0.5)
    model.simulate(duration=1000, time_step=10)
    
    fig, ax = model.plot()
    assert fig is not None
    assert ax is not None
    mock_subplots.assert_called_once()


def test_higuchi_plot_error():
    model = HiguchiModel(D=1e-6, c0=1.5, cs=0.5)
    
    with raises(ValueError, match="No simulation data available. Run simulate\\(\\) first."):
        model.plot()
    
    model.time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    model.release_profile = [0.0]  # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.plot()
    
    model.simulate(duration=1000, time_step=10)
    with mock.patch.dict('sys.modules', {'matplotlib': None}):
        with raises(ImportError, match="Matplotlib is required for plotting but not installed."):
            model.plot()


def test_higuchi_release_rate():
    model = HiguchiModel(D=1e-6, c0=1.5, cs=0.5)
    model.simulate(duration=1000, time_step=10)
    rate = model.get_release_rate()
    assert True # TODO: need a reference for rate based on a sequence of release values


def test_higuchi_release_rate_error():
    model = HiguchiModel(D=1e-6, c0=1.5, cs=0.5)
    
    with raises(ValueError, match="No simulation data available. Run simulate\\(\\) first."):
        model.get_release_rate()

    model.time_points = [0]  # manually set time points to simulate error (TODO: it will be caught with prior errors)
    model.release_profile = [0.0]  # manually set a too short profile to simulate error (TODO: it will be caught with prior errors)
    with raises(ValueError, match="Release profile is too short to calculate release rate."):
        model.get_release_rate()


def test_higuchi_time_for_release(): # Reference: https://www.wolframalpha.com/input?i=solve+for+t+in+sqrt%2810%5E%28-6%29*0.5*%282*1.5-0.5%29*t%29+%3D+0.5*sqrt%2810%5E%28-6%29*0.5*%282*1.5-0.5%29*1000%29
    model = HiguchiModel(D=1e-6, c0=1.5, cs=0.5)
    model.simulate(duration=1000, time_step=10)
    tx = model.time_for_release(0.5 * model.release_profile[-1])
    assert tx == 250.0


def test_higuchi_time_for_release_error():
    model = HiguchiModel(D=1e-6, c0=1.5, cs=0.5)
    
    with raises(ValueError, match="No simulation data available. Run simulate\\(\\) first."):
        model.time_for_release(0.5)
    model.simulate(duration=1000, time_step=10)
    
    with raises(ValueError, match="Target release must be between 0 and 1."):
        model.time_for_release(-0.1)
    
    with raises(ValueError, match="Target release must be between 0 and 1."):
        model.time_for_release(2.0)
    
    with raises(ValueError, match="Target release exceeds maximum release of the simulated duration."):
        model.time_for_release(min(model.release_profile[-1] + 0.1, 1))
