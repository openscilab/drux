# -*- coding: utf-8 -*-
"""Drux parameters and constants."""

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

DRUX_VERSION = "0.3"

MODELS_REGISTRY = {
    "zero_order": {
        "params": {
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
        },
        "equation": lambda p, t: p.M0 + p.k0 * t,
        "validation": {
            "k0": {"check": lambda v: v >= 0, "error": ERROR_ZERO_ORDER_RELEASE_RATE},
            "M0": {"check": lambda v: v >= 0, "error": ERROR_ZERO_ORDER_INITIAL_AMOUNT},
        },
    },
    "first_order": {
        "params": {
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
        },
        "equation": lambda p, t: p.M0 * (1 - exp(-p.k * t)),
        "validation": {
            "k": {"check": lambda v: v >= 0, "error": ERROR_FIRST_ORDER_RELEASE_RATE},
            "M0": {
                "check": lambda v: v >= 0,
                "error": ERROR_FIRST_ORDER_INITIAL_AMOUNT,
            },
        },
    },
    "higuchi": {
        "params": {
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
        },
        "equation": lambda p, t: sqrt(p.D * (2 * p.c0 - p.cs) * p.cs * t),
        "validation": {
            "D": {"check": lambda v: v > 0, "error": ERROR_INVALID_DIFFUSION},
            "c0": {"check": lambda v: v > 0, "error": ERROR_INVALID_CONCENTRATION},
            "cs": {"check": lambda v: v > 0, "error": ERROR_INVALID_SOLUBILITY},
            "cs_vs_c0": {
                "check": lambda p: p.cs <= p.c0,
                "error": ERROR_SOLUBILITY_HIGHER_THAN_CONCENTRATION,
                "cross_param": True,
            },
        },
    },
    "weibull": {
        "params": {
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
            "M": {
                "type": float,
                "description": "Entire releasable amount of drug",
                "unit": "mg",
                "default": 1,
            },
        },
        "equation": lambda p, t: p.M * (1 - exp(-p.a * t**p.b)),
        "validation": {
            "M": {"check": lambda v: v >= 0, "error": ERROR_RELEASABLE_AMOUNT},
            "a": {"check": lambda v: v > 0, "error": ERROR_WEIBULL_SCALE_PARAMETER},
            "b": {"check": lambda v: v > 0, "error": ERROR_WEIBULL_SHAPE_PARAMETER},
        },
    },
    "hopfenberg": {
        "params": {
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
            "M": {
                "type": float,
                "description": "Entire releasable amount of drug",
                "unit": "mg",
                "default": 1,
            },
        },
        "equation": lambda p, t: p.M * (1 - (1 - (p.k0 * t) / (p.c0 * p.a0)) ** p.n),
        "validation": {
            "M": {"check": lambda v: v >= 0, "error": ERROR_RELEASABLE_AMOUNT},
            "k0": {"check": lambda v: v >= 0, "error": ERROR_INVALID_EROSION_CONSTANT},
            "c0": {"check": lambda v: v > 0, "error": ERROR_INVALID_CONCENTRATION},
            "a0": {"check": lambda v: v > 0, "error": ERROR_INVALID_INITIAL_RADIUS},
            "n": {
                "check": lambda v: v in (1, 2, 3),
                "error": ERROR_INVALID_GEOMETRY_FACTOR,
            },
        },
    },
}
