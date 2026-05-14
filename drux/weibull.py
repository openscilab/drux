# -*- coding: utf-8 -*-
"""Drux Weibull model implementation."""

from .base_model import create_model_class

WeibullModel = create_model_class(
    model_name="weibull",
    class_name="WeibullModel",
    label="Weibull Model",
    docstring="Simulator for the Weibull drug release model using analytical expressions based on concentration conditions.",
)
