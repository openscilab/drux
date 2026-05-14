# -*- coding: utf-8 -*-
"""Drux Higuchi model implementation."""

from .base_model import create_model_class

HiguchiModel = create_model_class(
    model_name="higuchi",
    class_name="HiguchiModel",
    label="Higuchi Model",
    docstring="Simulator for the Higuchi drug release model using analytical expressions based on concentration conditions.",
)
