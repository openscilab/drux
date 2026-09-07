"""Tests for the CurveFit class implementation in drux package."""

from pytest import raises
from numpy import isclose, array
from math import exp, sqrt
from re import escape
from drux import CurveFit, ZeroOrderModel, FirstOrderModel, HiguchiModel, WeibullModel, HopfenbergModel

TEST_CASE_NAME = "CurveFit tests"
LONG_DURATION, LONG_TIME_STEP = 1000, 10
SHORT_DURATION, SHORT_TIME_STEP = 100, 1
LONG_TIME = list(range(0, LONG_DURATION + LONG_TIME_STEP, LONG_TIME_STEP))
SHORT_TIME = list(range(0, SHORT_DURATION + SHORT_TIME_STEP, SHORT_TIME_STEP))
RELATIVE_TOLERANCE = 1e-2
ABSOLUTE_TOLERANCE = 1e-6

# Parameters and release profiles below are the ones used in the corresponding model tests
ZERO_ORDER_M0, ZERO_ORDER_K0 = 0.01, 0.1
FIRST_ORDER_M0, FIRST_ORDER_K = 0.1, 0.003
HIGUCHI_D, HIGUCHI_C0, HIGUCHI_CS = 1e-6, 1, 0.5
WEIBULL_M, WEIBULL_A, WEIBULL_B = 1, 0.095, 0.7
HOPFENBERG_M, HOPFENBERG_K0, HOPFENBERG_C0, HOPFENBERG_A0, HOPFENBERG_N = 1, 0.00067, 0.0374, 3.51, 2

ZERO_ORDER_RELEASE = [ZERO_ORDER_M0 + (ZERO_ORDER_K0 * t) for t in LONG_TIME]
FIRST_ORDER_RELEASE = [FIRST_ORDER_M0 * (1 - exp(-FIRST_ORDER_K * t)) for t in LONG_TIME]
HIGUCHI_RELEASE = [sqrt(HIGUCHI_D * (2 * HIGUCHI_C0 - HIGUCHI_CS) * HIGUCHI_CS * t) for t in LONG_TIME]
WEIBULL_RELEASE = [WEIBULL_M * (1 - exp(-WEIBULL_A * t ** WEIBULL_B)) for t in SHORT_TIME]
HOPFENBERG_RELEASE = [
    HOPFENBERG_M * (1 - (1 - (HOPFENBERG_K0 * t) / (HOPFENBERG_C0 * HOPFENBERG_A0))**HOPFENBERG_N) for t in SHORT_TIME]

# Zero-mean alternating perturbation, kept deterministic so that the fit is reproducible
NOISE = [0.5 * (-1)**i for i in range(len(LONG_TIME))]
NOISY_ZERO_ORDER_RELEASE = [r + noise for r, noise in zip(ZERO_ORDER_RELEASE, NOISE)]
# A value that is exactly representable, so that the total sum of squares of the profile is exactly zero
CONSTANT_RELEASE_VALUE = 0.5
CONSTANT_RELEASE = [CONSTANT_RELEASE_VALUE for _ in LONG_TIME]


def test_curve_fit_parameters():
    fitter = CurveFit(
        model_name="hopfenberg",
        time=SHORT_TIME,
        release_profile=HOPFENBERG_RELEASE,
        known_parameters={"n": HOPFENBERG_N})
    assert fitter._model_name == "hopfenberg"
    assert fitter._model_class == HopfenbergModel
    assert fitter._time.tolist() == SHORT_TIME
    assert fitter._release_profile.tolist() == HOPFENBERG_RELEASE
    assert fitter._known_parameters == {"n": HOPFENBERG_N}
    assert fitter._parameter_names == ["M", "k0", "c0", "a0", "n"]
    assert fitter._free_parameters == ["M", "k0", "c0", "a0"]
    assert fitter._fit_result is None


def test_default_known_parameters():
    fitter = CurveFit(model_name="zero_order", time=LONG_TIME, release_profile=ZERO_ORDER_RELEASE)
    assert fitter._known_parameters == {}
    assert fitter._parameter_names == ["M0", "k0"]
    assert fitter._free_parameters == ["M0", "k0"]


def test_invalid_parameters():
    with raises(ValueError, match=escape("Unknown model 'not_a_model'. Available models: ")):
        CurveFit(model_name="not_a_model", time=LONG_TIME, release_profile=ZERO_ORDER_RELEASE)

    with raises(ValueError, match=escape("Time and release profile must have the same length.")):
        CurveFit(model_name="zero_order", time=LONG_TIME, release_profile=ZERO_ORDER_RELEASE[:-1])

    with raises(ValueError, match=escape("At least 2 data points are required for fitting.")):
        CurveFit(model_name="zero_order", time=[0], release_profile=[ZERO_ORDER_M0])

    with raises(ValueError, match=escape(
            "Unknown parameter(s) ['a', 'b'] for model 'zero_order'. Valid parameters: ['M0', 'k0'].")):
        CurveFit(
            model_name="zero_order",
            time=LONG_TIME,
            release_profile=ZERO_ORDER_RELEASE,
            known_parameters={"a": WEIBULL_A, "b": WEIBULL_B})

    with raises(ValueError, match=escape("All parameters are known; there is nothing left to fit.")):
        CurveFit(
            model_name="zero_order",
            time=LONG_TIME,
            release_profile=ZERO_ORDER_RELEASE,
            known_parameters={"M0": ZERO_ORDER_M0, "k0": ZERO_ORDER_K0})


def test_repr():
    fitter = CurveFit(model_name="zero_order", time=LONG_TIME, release_profile=ZERO_ORDER_RELEASE)
    assert repr(fitter) == "drux.CurveFit(zero_order, not fitted)"

    fitter.fit()
    assert repr(fitter) == "drux.CurveFit(zero_order: M0=0.0100, k0=0.1000, R^2=1.0000)"


def test_zero_order_fit():  # Reference: https://europepmc.org/article/pmc/3425064
    fitter = CurveFit(model_name="zero_order", time=LONG_TIME, release_profile=ZERO_ORDER_RELEASE)
    result = fitter.fit()
    assert result.model_name == "zero_order"
    assert isclose(result.parameters["M0"], ZERO_ORDER_M0, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.parameters["k0"], ZERO_ORDER_K0, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.r_squared, 1.0, atol=ABSOLUTE_TOLERANCE)
    assert isinstance(result.model, ZeroOrderModel)
    assert isclose(result.model._parameters.M0, ZERO_ORDER_M0, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.model._parameters.k0, ZERO_ORDER_K0, rtol=RELATIVE_TOLERANCE)


def test_zero_order_fit_with_known_parameters():
    fitter = CurveFit(
        model_name="zero_order",
        time=LONG_TIME,
        release_profile=ZERO_ORDER_RELEASE,
        known_parameters={"M0": ZERO_ORDER_M0})
    result = fitter.fit()
    assert result.parameters["M0"] == ZERO_ORDER_M0
    assert isclose(result.parameters["k0"], ZERO_ORDER_K0, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.r_squared, 1.0, atol=ABSOLUTE_TOLERANCE)


def test_first_order_fit():  # Reference: https://europepmc.org/article/pmc/3425064
    fitter = CurveFit(model_name="first_order", time=LONG_TIME, release_profile=FIRST_ORDER_RELEASE)
    result = fitter.fit()
    assert result.model_name == "first_order"
    assert isclose(result.parameters["M0"], FIRST_ORDER_M0, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.parameters["k"], FIRST_ORDER_K, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.r_squared, 1.0, atol=ABSOLUTE_TOLERANCE)
    assert isinstance(result.model, FirstOrderModel)


def test_first_order_fit_with_known_parameters():
    fitter = CurveFit(
        model_name="first_order",
        time=LONG_TIME,
        release_profile=FIRST_ORDER_RELEASE,
        known_parameters={"M0": FIRST_ORDER_M0})
    result = fitter.fit()
    assert result.parameters["M0"] == FIRST_ORDER_M0
    assert isclose(result.parameters["k"], FIRST_ORDER_K, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.r_squared, 1.0, atol=ABSOLUTE_TOLERANCE)


def test_higuchi_fit():  # c0 and cs are known, since only the D * (2*c0 - cs) * cs product is identifiable
    fitter = CurveFit(
        model_name="higuchi",
        time=LONG_TIME,
        release_profile=HIGUCHI_RELEASE,
        known_parameters={"c0": HIGUCHI_C0, "cs": HIGUCHI_CS})
    result = fitter.fit()
    assert result.model_name == "higuchi"
    assert result.parameters["c0"] == HIGUCHI_C0
    assert result.parameters["cs"] == HIGUCHI_CS
    assert isclose(result.parameters["D"], HIGUCHI_D, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.r_squared, 1.0, atol=ABSOLUTE_TOLERANCE)
    assert isinstance(result.model, HiguchiModel)


def test_weibull_fit():  # Reference: https://www.mdpi.com/2073-4360/13/17/2897
    fitter = CurveFit(model_name="weibull", time=SHORT_TIME, release_profile=WEIBULL_RELEASE)
    result = fitter.fit()
    assert result.model_name == "weibull"
    assert isclose(result.parameters["M"], WEIBULL_M, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.parameters["a"], WEIBULL_A, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.parameters["b"], WEIBULL_B, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.r_squared, 1.0, atol=ABSOLUTE_TOLERANCE)
    assert isinstance(result.model, WeibullModel)


def test_weibull_fit_with_known_parameters():
    fitter = CurveFit(
        model_name="weibull",
        time=SHORT_TIME,
        release_profile=WEIBULL_RELEASE,
        known_parameters={"M": WEIBULL_M})
    result = fitter.fit()
    assert result.parameters["M"] == WEIBULL_M
    assert isclose(result.parameters["a"], WEIBULL_A, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.parameters["b"], WEIBULL_B, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.r_squared, 1.0, atol=ABSOLUTE_TOLERANCE)


def test_fit_interleaved_known_parameters():  # a known parameter declared between two free ones
    fitter = CurveFit(
        model_name="weibull",
        time=SHORT_TIME,
        release_profile=WEIBULL_RELEASE,
        known_parameters={"a": WEIBULL_A})
    result = fitter.fit()
    assert fitter._free_parameters == ["M", "b"]
    assert list(result.parameters) == ["M", "a", "b"]  # the model declaration order is kept
    assert result.parameters["a"] == WEIBULL_A
    assert isclose(result.parameters["M"], WEIBULL_M, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.parameters["b"], WEIBULL_B, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.r_squared, 1.0, atol=ABSOLUTE_TOLERANCE)


# n is a geometry factor, so it is known in every test; c0 and a0 are known, since only their product is identifiable
def test_hopfenberg_fit():  # Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC3500559/
    fitter = CurveFit(
        model_name="hopfenberg",
        time=SHORT_TIME,
        release_profile=HOPFENBERG_RELEASE,
        known_parameters={"c0": HOPFENBERG_C0, "a0": HOPFENBERG_A0, "n": HOPFENBERG_N})
    result = fitter.fit()
    assert result.model_name == "hopfenberg"
    assert result.parameters["c0"] == HOPFENBERG_C0
    assert result.parameters["a0"] == HOPFENBERG_A0
    assert result.parameters["n"] == HOPFENBERG_N
    assert isclose(result.parameters["M"], HOPFENBERG_M, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.parameters["k0"], HOPFENBERG_K0, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.r_squared, 1.0, atol=ABSOLUTE_TOLERANCE)
    assert isinstance(result.model, HopfenbergModel)


def test_hopfenberg_fit_with_known_parameters():
    fitter = CurveFit(
        model_name="hopfenberg",
        time=SHORT_TIME,
        release_profile=HOPFENBERG_RELEASE,
        known_parameters={"M": HOPFENBERG_M, "c0": HOPFENBERG_C0, "a0": HOPFENBERG_A0, "n": HOPFENBERG_N})
    result = fitter.fit()
    assert result.parameters["M"] == HOPFENBERG_M
    assert isclose(result.parameters["k0"], HOPFENBERG_K0, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.r_squared, 1.0, atol=ABSOLUTE_TOLERANCE)


def test_fit_initial_guess():
    fitter = CurveFit(
        model_name="higuchi",
        time=LONG_TIME,
        release_profile=HIGUCHI_RELEASE,
        known_parameters={"c0": HIGUCHI_C0, "cs": HIGUCHI_CS})
    result = fitter.fit(initial_guess=[HIGUCHI_D])
    assert isclose(result.parameters["D"], HIGUCHI_D, rtol=RELATIVE_TOLERANCE)
    assert isclose(result.r_squared, 1.0, atol=ABSOLUTE_TOLERANCE)


def test_fit_bounds():
    upper_k0 = ZERO_ORDER_K0 / 2  # the release rate is capped below its actual value
    fitter = CurveFit(model_name="zero_order", time=LONG_TIME, release_profile=ZERO_ORDER_RELEASE)
    result = fitter.fit(initial_guess=[ZERO_ORDER_M0, upper_k0], bounds=([1e-10, 1e-10], [1.0, upper_k0]))
    assert isclose(result.parameters["k0"], upper_k0, rtol=RELATIVE_TOLERANCE)
    assert result.parameters["M0"] <= 1.0
    assert result.r_squared < 1.0


def test_fit_noisy_release_profile():
    fitter = CurveFit(model_name="zero_order", time=LONG_TIME, release_profile=NOISY_ZERO_ORDER_RELEASE)
    result = fitter.fit()
    assert isclose(result.parameters["k0"], ZERO_ORDER_K0, rtol=RELATIVE_TOLERANCE)
    assert 0.99 < result.r_squared < 1.0


def test_fit_constant_release_profile():
    fitter = CurveFit(model_name="zero_order", time=LONG_TIME, release_profile=CONSTANT_RELEASE)
    result = fitter.fit()
    assert isclose(result.parameters["M0"], CONSTANT_RELEASE_VALUE, rtol=RELATIVE_TOLERANCE)
    assert result.r_squared == 0.0  # the total sum of squares is zero, so R^2 is undefined


def test_fit_sequence_types():
    list_result = CurveFit(model_name="zero_order", time=LONG_TIME, release_profile=ZERO_ORDER_RELEASE).fit()
    tuple_result = CurveFit(
        model_name="zero_order",
        time=tuple(LONG_TIME),
        release_profile=tuple(ZERO_ORDER_RELEASE)).fit()
    array_result = CurveFit(
        model_name="zero_order",
        time=array(LONG_TIME),
        release_profile=array(ZERO_ORDER_RELEASE)).fit()
    assert list_result.parameters == tuple_result.parameters == array_result.parameters


def test_fitted_model_simulation():
    fitter = CurveFit(model_name="weibull", time=SHORT_TIME, release_profile=WEIBULL_RELEASE)
    result = fitter.fit()
    profile = result.model.simulate(duration=SHORT_DURATION, time_step=SHORT_TIME_STEP)
    assert all(isclose(p, r, rtol=RELATIVE_TOLERANCE) for p, r in zip(profile, WEIBULL_RELEASE))


def test_get_result():
    fitter = CurveFit(model_name="zero_order", time=LONG_TIME, release_profile=ZERO_ORDER_RELEASE)
    result = fitter.fit()
    assert fitter.get_result() is result

    upper_k0 = ZERO_ORDER_K0 / 2
    bounded_result = fitter.fit(
        initial_guess=[ZERO_ORDER_M0, upper_k0],
        bounds=([1e-10, 1e-10], [1.0, upper_k0]))
    assert fitter.get_result() is bounded_result  # the latest fit result is kept


def test_get_result_error():
    fitter = CurveFit(model_name="zero_order", time=LONG_TIME, release_profile=ZERO_ORDER_RELEASE)

    with raises(ValueError, match=escape("No fit result available. Run fit() first.")):
        fitter.get_result()
