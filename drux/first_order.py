# -*- coding: utf-8 -*-
"""Drux first-order model implementation."""

from .base_model import create_model_class

FirstOrderModel = create_model_class(
    model_name="first_order",
    class_name="FirstOrderModel",
    label="First-Order Model",
    docstring="Simulator for the first-order drug release model.",
)
