# -*- coding: utf-8 -*-
"""Drux zero-order model implementation."""

from .base_model import create_model_class

ZeroOrderModel = create_model_class(
    model_name="zero_order",
    class_name="ZeroOrderModel",
    label="Zero-Order Model",
    docstring="Simulator for the zero-order drug release model.",
)
