# -*- coding: utf-8 -*-
"""Drux parameters and constants."""

from .equations import get_equation, get_validation

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
        "equation": get_equation("zero_order"),
        "validation": get_validation("zero_order"),
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
        "equation": get_equation("first_order"),
        "validation": get_validation("first_order"),
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
        "equation": get_equation("higuchi"),
        "validation": get_validation("higuchi"),
    },
    "weibull": {
        "params": {
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
        },
        "equation": get_equation("weibull"),
        "validation": get_validation("weibull"),
    },
    "hopfenberg": {
        "params": {
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
        },
        "equation": get_equation("hopfenberg"),
        "validation": get_validation("hopfenberg"),
    },
}
