# -*- coding: utf-8 -*-
"""Centralized model registry for drug release models."""

from math import exp, sqrt

from .messages import (
    ERROR_ZERO_ORDER_RELEASE_RATE,
    ERROR_ZERO_ORDER_INITIAL_AMOUNT,
    ERROR_FIRST_ORDER_RELEASE_RATE,
    ERROR_FIRST_ORDER_INITIAL_AMOUNT,
    ERROR_INVALID_DIFFUSION,
    ERROR_INVALID_CONCENTRATION,
    ERROR_INVALID_SOLUBILITY,
    ERROR_SOLUBILITY_HIGHER_THAN_CONCENTRATION,
    ERROR_RELEASABLE_AMOUNT,
    ERROR_WEIBULL_SCALE_PARAMETER,
    ERROR_WEIBULL_SHAPE_PARAMETER,
    ERROR_INVALID_EROSION_CONSTANT,
    ERROR_INVALID_INITIAL_RADIUS,
    ERROR_INVALID_GEOMETRY_FACTOR,
)


# ---------------------------------------------------------------------------
# Internal registries
# ---------------------------------------------------------------------------
_EQUATIONS = {}
_VALIDATIONS = {}
_PARAMS = {}


# ---------------------------------------------------------------------------
# Registration helpers
# ---------------------------------------------------------------------------

def equation(name):
    """Register a simulation equation for the given model name.

    The decorated function must have the signature ``(p, t) -> float``
    where *p* is a parameters namespace instance.
    """
    def decorator(fn):
        _EQUATIONS[name] = fn
        return fn
    return decorator


def register_validation(name, rules):
    """Register validation rules for the given model name.

    :param name: Model name key.
    :param rules: Dict of validation rules, e.g.::

        {
            "param_name": {"check": callable, "error": str},
            ...
        }

    Cross-parameter validations should include ``"cross_param": True``.
    """
    _VALIDATIONS[name] = rules


def register_params(name, params):
    """Register parameter metadata for the given model name.

    :param name: Model name key.
    :param params: OrderedDict-style dict of parameter specs, e.g.::

        {
            "k0": {"type": float, "description": "...", "unit": "...", "default": None},
            ...
        }
    """
    _PARAMS[name] = params


# ---------------------------------------------------------------------------
# Public accessors
# ---------------------------------------------------------------------------

def get_model_config(name):
    """Return the full model configuration dict for *name*.

    :returns: ``{"params": ..., "equation": ..., "validation": ...}``
    :raises KeyError: if *name* is not registered.
    """
    return {
        "params": _PARAMS[name],
        "equation": _EQUATIONS[name],
        "validation": _VALIDATIONS[name],
    }


def get_equation(name):
    """Return the simulation equation for *name*."""
    return _EQUATIONS[name]


def get_validation(name):
    """Return the validation rules for *name*."""
    return _VALIDATIONS[name]


def get_params(name):
    """Return the parameter metadata for *name*."""
    return _PARAMS[name]


def list_models():
    """Return a list of registered model names."""
    return list(_EQUATIONS.keys())


# ---------------------------------------------------------------------------
# Zero-order model
# ---------------------------------------------------------------------------

register_params("zero_order", {
    "k0": {
        "type": float,
        "description": "Zero-order release rate constant",
        "unit": "mg/s",
        "default": None,
    },
    "M0": {
        "type": float,
        "description": "Initial amount of drug in the solution",
        "unit": "mg",
        "default": 0,
    },
})


@equation("zero_order")
def zero_order(p, t):
    """Zero-order kinetics: M(t) = M0 + k0 * t."""
    return p.M0 + p.k0 * t


register_validation("zero_order", {
    "k0": {"check": lambda v: v >= 0, "error": ERROR_ZERO_ORDER_RELEASE_RATE},
    "M0": {"check": lambda v: v >= 0, "error": ERROR_ZERO_ORDER_INITIAL_AMOUNT},
})


# ---------------------------------------------------------------------------
# First-order model
# ---------------------------------------------------------------------------

register_params("first_order", {
    "k": {
        "type": float,
        "description": "First-order release rate constant",
        "unit": "1/s",
        "default": None,
    },
    "M0": {
        "type": float,
        "description": "Entire releasable amount of drug",
        "unit": "mg",
        "default": 1,
    },
})


@equation("first_order")
def first_order(p, t):
    """First-order kinetics: M(t) = M0 * (1 - exp(-k * t))."""
    return p.M0 * (1 - exp(-p.k * t))


register_validation("first_order", {
    "k": {"check": lambda v: v >= 0, "error": ERROR_FIRST_ORDER_RELEASE_RATE},
    "M0": {"check": lambda v: v >= 0, "error": ERROR_FIRST_ORDER_INITIAL_AMOUNT},
})


# ---------------------------------------------------------------------------
# Higuchi model
# ---------------------------------------------------------------------------

register_params("higuchi", {
    "D": {
        "type": float,
        "description": "Drug diffusivity in the polymer carrier",
        "unit": "cm^2/s",
        "default": None,
    },
    "c0": {
        "type": float,
        "description": "Initial drug concentration",
        "unit": "mg/cm^3",
        "default": None,
    },
    "cs": {
        "type": float,
        "description": "Drug solubility in the polymer",
        "unit": "mg/cm^3",
        "default": None,
    },
})


@equation("higuchi")
def higuchi(p, t):
    """Higuchi model: M(t) = sqrt(D * (2*c0 - cs) * cs * t)."""
    return sqrt(p.D * (2 * p.c0 - p.cs) * p.cs * t)


register_validation("higuchi", {
    "D": {"check": lambda v: v > 0, "error": ERROR_INVALID_DIFFUSION},
    "c0": {"check": lambda v: v > 0, "error": ERROR_INVALID_CONCENTRATION},
    "cs": {"check": lambda v: v > 0, "error": ERROR_INVALID_SOLUBILITY},
    "cs_vs_c0": {
        "check": lambda p: p.cs <= p.c0,
        "error": ERROR_SOLUBILITY_HIGHER_THAN_CONCENTRATION,
        "cross_param": True,
    },
})


# ---------------------------------------------------------------------------
# Weibull model
# ---------------------------------------------------------------------------

register_params("weibull", {
    "M": {
        "type": float,
        "description": "Entire releasable amount of drug",
        "unit": "mg",
        "default": 1,
    },
    "a": {
        "type": float,
        "description": "Scale factor",
        "unit": "dimensionless",
        "default": None,
    },
    "b": {
        "type": float,
        "description": "Shape factor",
        "unit": "dimensionless",
        "default": None,
    },
})


@equation("weibull")
def weibull(p, t):
    """Weibull model: M(t) = M * (1 - exp(-a * t^b))."""
    return p.M * (1 - exp(-p.a * t ** p.b))


register_validation("weibull", {
    "M": {"check": lambda v: v >= 0, "error": ERROR_RELEASABLE_AMOUNT},
    "a": {"check": lambda v: v > 0, "error": ERROR_WEIBULL_SCALE_PARAMETER},
    "b": {"check": lambda v: v > 0, "error": ERROR_WEIBULL_SHAPE_PARAMETER},
})


# ---------------------------------------------------------------------------
# Hopfenberg model
# ---------------------------------------------------------------------------

register_params("hopfenberg", {
    "M": {
        "type": float,
        "description": "Entire releasable amount of drug",
        "unit": "mg",
        "default": 1,
    },
    "k0": {
        "type": float,
        "description": "Erosion rate constant",
        "unit": "mg/(mm^2·s)",
        "default": None,
    },
    "c0": {
        "type": float,
        "description": "Initial drug concentration in the matrix",
        "unit": "mg/mm^3",
        "default": None,
    },
    "a0": {
        "type": float,
        "description": "Initial radius or half-thickness of the device",
        "unit": "mm",
        "default": None,
    },
    "n": {
        "type": "int",
        "description": "Geometry factor (1=slab, 2=cylinder, 3=sphere)",
        "unit": "dimensionless",
        "default": None,
    },
})


@equation("hopfenberg")
def hopfenberg(p, t):
    """Hopfenberg model: M(t) = M * (1 - (1 - k0*t/(c0*a0))^n)."""
    return p.M * (1 - (1 - (p.k0 * t) / (p.c0 * p.a0)) ** p.n)


register_validation("hopfenberg", {
    "M": {"check": lambda v: v >= 0, "error": ERROR_RELEASABLE_AMOUNT},
    "k0": {"check": lambda v: v >= 0, "error": ERROR_INVALID_EROSION_CONSTANT},
    "c0": {"check": lambda v: v > 0, "error": ERROR_INVALID_CONCENTRATION},
    "a0": {"check": lambda v: v > 0, "error": ERROR_INVALID_INITIAL_RADIUS},
    "n": {
        "check": lambda v: v in (1, 2, 3),
        "error": ERROR_INVALID_GEOMETRY_FACTOR,
    },
})
