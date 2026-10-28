# -*- coding: utf-8 -*-
"""Drux modules."""

from .params import DRUX_VERSION
from .higuchi import HiguchiModel, HiguchiParameters
from .zero_order import ZeroOrderModel, ZeroOrderParameters
from .first_order import FirstOrderModel, FirstOrderParameters
from .weibull import WeibullModel, WeibullParameters
from .hopfenberg import HopfenbergModel, HopfenbergParameters
from .quadratic import QuadraticModel, QuadraticParameters
from .cubic import CubicModel, CubicParameters
from .curve_fit import CurveFit, FitResult

__version__ = DRUX_VERSION
