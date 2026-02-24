# -*- coding: utf-8 -*-
"""Generic utilities for Drux."""

from dataclasses import field, make_dataclass
from .params import MODELS_REGISTRY


def create_parameters_dataclass(model_name: str):
    """
    Create a dataclass from MODELS_REGISTRY configuration.

    :param model_name: Name of the model in MODELS_REGISTRY
    :return: Dataclass type for model parameters
    """
    config = MODELS_REGISTRY[model_name]
    fields = []

    for param_name, param_info in config["params"].items():
        param_type = param_info["type"]
        default_value = param_info["default"]

        if default_value is not None:
            fields.append((param_name, param_type, field(default=default_value)))
        else:
            fields.append((param_name, param_type))

    class_name = f"{model_name.title().replace('_', '')}Parameters"
    return make_dataclass(class_name, fields)
