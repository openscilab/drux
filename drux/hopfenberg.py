# -*- coding: utf-8 -*-
"""Drux Hopfenberg model implementation."""

from .base_model import create_model_class

HopfenbergModel = create_model_class(
    model_name="hopfenberg",
    class_name="HopfenbergModel",
    label="Hopfenberg Model",
    docstring="Simulator for the Hopfenberg drug release model for surface-eroding polymers.",
)
