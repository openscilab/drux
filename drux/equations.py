# -*- coding: utf-8 -*-
"""Centralized equation and validation registry for drug release models."""

from math import exp, sqrt

# ---------------------------------------------------------------------------
# Internal registries – accessed via helper functions, not directly.
# ---------------------------------------------------------------------------
_EQUATIONS = {}
_VALIDATIONS = {}


# ---------------------------------------------------------------------------
# Decorators
# ---------------------------------------------------------------------------

def equation(name):
    """Register a simulation equation for the given model name.

    The decorated function must have the signature ``(p, t) -> float``
    where *p* is a parameters dataclass instance.
    """
    def decorator(fn):
        _EQUATIONS[name] = fn
        return fn
    return decorator


def validation(name):
    """
    Register validation rules for the given model name.

    The decorated function must return a dict of validation rules in the form::

        {
            "param_name": {"check": callable, "error": str},
            ...
        }

    Cross-parameter validations should include ``"cross_param": True``.
    """
    def decorator(fn):
        _VALIDATIONS[name] = fn()
        return fn
    return decorator


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def get_equation(name):
    """Return the simulation equation for *name*, or raise KeyError."""
    return _EQUATIONS[name]


def get_validation(name):
    """Return the validation rules for *name*, or raise KeyError."""
    return _VALIDATIONS[name]


def list_equations():
    """Return a list of registered equation names."""
    return list(_EQUATIONS.keys())


# ---------------------------------------------------------------------------
# Simulation equations  (p, t) -> float
# ---------------------------------------------------------------------------

@equation("zero_order")
def zero_order(p, t):
    """Zero-order kinetics: M(t) = M0 + k0 * t."""
    return p.M0 + p.k0 * t


@equation("first_order")
def first_order(p, t):
    """First-order kinetics: M(t) = M0 * (1 - exp(-k * t))."""
    return p.M0 * (1 - exp(-p.k * t))


@equation("higuchi")
def higuchi(p, t):
    """Higuchi model: M(t) = sqrt(D * (2*c0 - cs) * cs * t)."""
    return sqrt(p.D * (2 * p.c0 - p.cs) * p.cs * t)


@equation("weibull")
def weibull(p, t):
    """Weibull model: M(t) = M * (1 - exp(-a * t^b))."""
    return p.M * (1 - exp(-p.a * t ** p.b))


@equation("hopfenberg")
def hopfenberg(p, t):
    """Hopfenberg model: M(t) = M * (1 - (1 - k0*t/(c0*a0))^n)."""
    return p.M * (1 - (1 - (p.k0 * t) / (p.c0 * p.a0)) ** p.n)


# ---------------------------------------------------------------------------
# Validation rules – co-located with equations
# ---------------------------------------------------------------------------

@validation("zero_order")
def zero_order_validation():
    """Rules for validating the zero-order model."""
    from .messages import ERROR_ZERO_ORDER_RELEASE_RATE, ERROR_ZERO_ORDER_INITIAL_AMOUNT
    return {
        "k0": {"check": lambda v: v >= 0, "error": ERROR_ZERO_ORDER_RELEASE_RATE},
        "M0": {"check": lambda v: v >= 0, "error": ERROR_ZERO_ORDER_INITIAL_AMOUNT},
    }


@validation("first_order")
def first_order_validation():
    """Rules for validating the first-order model."""
    from .messages import ERROR_FIRST_ORDER_RELEASE_RATE, ERROR_FIRST_ORDER_INITIAL_AMOUNT
    return {
        "k": {"check": lambda v: v >= 0, "error": ERROR_FIRST_ORDER_RELEASE_RATE},
        "M0": {"check": lambda v: v >= 0, "error": ERROR_FIRST_ORDER_INITIAL_AMOUNT},
    }


@validation("higuchi")
def higuchi_validation():
    """Rules for validating the Higuchi model."""
    from .messages import (
        ERROR_INVALID_DIFFUSION,
        ERROR_INVALID_CONCENTRATION,
        ERROR_INVALID_SOLUBILITY,
        ERROR_SOLUBILITY_HIGHER_THAN_CONCENTRATION,
    )
    return {
        "D": {"check": lambda v: v > 0, "error": ERROR_INVALID_DIFFUSION},
        "c0": {"check": lambda v: v > 0, "error": ERROR_INVALID_CONCENTRATION},
        "cs": {"check": lambda v: v > 0, "error": ERROR_INVALID_SOLUBILITY},
        "cs_vs_c0": {
            "check": lambda p: p.cs <= p.c0,
            "error": ERROR_SOLUBILITY_HIGHER_THAN_CONCENTRATION,
            "cross_param": True,
        },
    }


@validation("weibull")
def weibull_validation():
    """Rules for validating the Weibull model."""
    from .messages import (
        ERROR_RELEASABLE_AMOUNT,
        ERROR_WEIBULL_SCALE_PARAMETER,
        ERROR_WEIBULL_SHAPE_PARAMETER,
    )
    return {
        "M": {"check": lambda v: v >= 0, "error": ERROR_RELEASABLE_AMOUNT},
        "a": {"check": lambda v: v > 0, "error": ERROR_WEIBULL_SCALE_PARAMETER},
        "b": {"check": lambda v: v > 0, "error": ERROR_WEIBULL_SHAPE_PARAMETER},
    }


@validation("hopfenberg")
def hopfenberg_validation():
    """Rules for validating the Hopfenberg model."""
    from .messages import (
        ERROR_RELEASABLE_AMOUNT,
        ERROR_INVALID_EROSION_CONSTANT,
        ERROR_INVALID_CONCENTRATION,
        ERROR_INVALID_INITIAL_RADIUS,
        ERROR_INVALID_GEOMETRY_FACTOR,
    )
    return {
        "M": {"check": lambda v: v >= 0, "error": ERROR_RELEASABLE_AMOUNT},
        "k0": {"check": lambda v: v >= 0, "error": ERROR_INVALID_EROSION_CONSTANT},
        "c0": {"check": lambda v: v > 0, "error": ERROR_INVALID_CONCENTRATION},
        "a0": {"check": lambda v: v > 0, "error": ERROR_INVALID_INITIAL_RADIUS},
        "n": {
            "check": lambda v: v in (1, 2, 3),
            "error": ERROR_INVALID_GEOMETRY_FACTOR,
        },
    }
