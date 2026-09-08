<div align="center">
<img src="https://github.com/openscilab/drux/raw/main/otherfiles/logo.png" width="350">
    <h1>Drux: Drug Release Analysis Framework</h1>
    <br/>
    <a href="https://badge.fury.io/py/drux"><img src="https://badge.fury.io/py/drux.svg" alt="PyPI version"></a>
    <a href="https://www.python.org/"><img src="https://img.shields.io/badge/built%20with-Python3-green.svg" alt="built with Python3"></a>
    <a href="https://codecov.io/gh/openscilab/drux"><img src="https://codecov.io/gh/openscilab/drux/branch/dev/graph/badge.svg?token=5O41J3XX2L"></a>
    <a href="https://github.com/openscilab/drux"><img alt="GitHub repo size" src="https://img.shields.io/github/repo-size/openscilab/drux"></a>
    <a href="https://discord.gg/8Rf6bGBtse"><img src="https://img.shields.io/discord/1064533716615049236.svg" alt="Discord Channel"></a>
</div>

----------


## Overview
<p align="justify">
Drux is a Python-based framework for simulating drug release profiles using mathematical models. It offers a reproducible and extensible platform to model, analyze, and visualize time-dependent drug release behavior, making it ideal for pharmaceutical research and development. By combining simplicity with scientific rigor, Drux provides a robust foundation for quantitative analysis of drug delivery kinetics.
</p>
<table>
    <tr>
        <td align="center">PyPI Counter</td>
        <td align="center">
            <a href="https://pepy.tech/projects/drux">
                <img src="https://static.pepy.tech/badge/drux">
            </a>
        </td>
    </tr>
    <tr>
        <td align="center">Github Stars</td>
        <td align="center">
            <a href="https://github.com/openscilab/drux">
                <img src="https://img.shields.io/github/stars/openscilab/drux.svg?style=social&label=Stars">
            </a>
        </td>
    </tr>
</table>
<table>
    <tr> 
        <td align="center">Branch</td>
        <td align="center">main</td>
        <td align="center">dev</td>
    </tr>
    <tr>
        <td align="center">CI</td>
        <td align="center">
            <img src="https://github.com/openscilab/drux/actions/workflows/test.yml/badge.svg?branch=main">
        </td>
        <td align="center">
            <img src="https://github.com/openscilab/drux/actions/workflows/test.yml/badge.svg?branch=dev">
            </td>
    </tr>
</table>
<table>
    <tr> 
        <td align="center">Code Quality</td>
        <td align="center"><a href="https://www.codefactor.io/repository/github/openscilab/drux"><img src="https://www.codefactor.io/repository/github/openscilab/drux/badge" alt="CodeFactor"></a></td>
        <td align="center"><a href="https://app.codacy.com/gh/openscilab/drux/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade"><img src="https://app.codacy.com/project/badge/Grade/06ed95529d284c81a846205baa1f4c6a"></a></td>
    </tr>
</table>


## Installation

### PyPI
- Check [Python Packaging User Guide](https://packaging.python.org/installing/)
- Run `pip install drux==0.4`
### Source code
- Download [Version 0.4](https://github.com/openscilab/drux/archive/v0.4.zip) or [Latest Source](https://github.com/openscilab/drux/archive/dev.zip)
- Run `pip install .`

## Supported Models
### Zero-Order
The Zero-Order model describes a constant rate of drug release over time. According to this model, the cumulative amount of drug released at time $t$ is given by:

$$
M_t = M_0 + k_0 t
$$

where:
- $M_t (mg)$ is the cumulative absolute amount of drug released at time $t$.
- $M_0 (mg)$ is the initial amount of drug in the system. $M_0$ defaults to zero in this model.
- $k_0 (\frac{mg}{s})$ is the zero-order release rate constant.

#### Applications
1. Tablets with extended release
2. Transdermal Patches
3. Implantable Device
4. Intraocular Implants
5. Infusion Systems

### First-Order
The first-order drug release model describes a process where the rate of drug release is proportional to the remaining amount of drug in the system. According to this model, the cumulative amount of drug released at time $t$ is given by:

$$
M_t = M_0 (1 - e^{-kt})
$$

where:
- $M_t (mg)$ is the cumulative absolute amount of drug released at time $t$.
- $M_0 (mg)$ is entire releasable amount of drug (the asymptotic maximum).
- $k (\frac{1}{s})$ is the first-order release rate constant.

#### Applications
1. Immediate-release tablets and capsules
2. Liquid drug formulations (oral solutions, intravenous injections)
3. Controlled-release matrix systems
4. Elastomeric infusion pumps

### Higuchi
The Higuchi model describes the release of a drug from a matrix system, where the drug diffuses through a porous medium.
The Higuchi equation addressed important aspects of drug transport and release from planar
devices. According to this model, the cumulative amount of drug released at time $t$ is given by:

$$
M_t = \sqrt{D(2c_0 - c_s)c_st}
$$

where:
- $M_t (\frac{mg}{cm^2})$ is the cumulative absolute amount of drug released at time $t$
- $D ({\frac{cm^2}{s}})$ is the drug diffusivity in the polymer carrier
- $c_0 (\frac{mg}{cm^3})$ is the initial drug concentration (total concentration of drug in the matrix)
- $c_s (\frac{mg}{cm^3})$ is the solubility of the drug in the polymer (carrier)

⚠️ The Higuchi model assumes that $c_0 \ge c_s$

#### Applications
1. Matrix Tablets
2. Hydrophilic polymer matrices
3. Controlled - Release Microspheres
4. Semisolid Systems
5. Implantable Drug delivery systems

### Weibull
Weibull model is an empirical model used to describe drug release kinetics from various pharmaceutical dosage forms. It is characterized by a flexible empirical equation that captures various release kinetics. 
According to this model, the cumulative amount of drug released at time $t$ is given by:

$$
M_t = M \left(1 - e^{-at^b}\right)
$$

where:
- $M_t (mg)$ is the cumulative absolute amount of drug released at time $t$
- $M (mg)$ is the total amount of drug released at infinite time
- $a$ is the scale parameter (related to the release rate)
- $b$ is the shape parameter (indicates the release mechanism)

#### Applications
1. Controlled-release modeling
2. Dissolution profiling
3. Comparative studies
4. Vivo predictions

### Hopfenberg
The Hopfenberg model describes drug release from surface-eroding polymers with a constant surface area. 
It is particularly useful for modeling drug release from biodegradable polymers where the drug is uniformly distributed throughout the matrix.
According to this model, the cumulative amount of drug released at time $t$ is given by:

$$
M_t = M_{\infty} \left(1 - \left(1 - \frac{k_0t}{c_0 a_0}\right)^n\right)
$$

where:
- $M_t (mg)$ is the cumulative absolute amount of drug released at time $t
- $M_{\infty} (mg)$ is the total amount of drug released at infinite time
- $k_0 (\frac{mg}{mm^2 s})$ is the surface erosion rate constant
- $c_0 (\frac{mg}{mm^3})$ is the initial drug concentration in the polymer matrix
- $a_0 (mm)$ is the initial radius (for spheres) or half-thickness (for slabs) of the device
- $n$ is the geometry-dependent exponent (1 for slabs, 2 for cylinders, 3 for spheres)

#### Applications
1. Biodegradable Implants
2. Surface-eroding drug delivery systems
3. Transdermal Patches
4. Injectable depots

## Curve Fitting
The models above calculate a release profile from known parameters. The `CurveFit` class does the opposite operation. It calculates the parameters of a model from measured experimental data. `CurveFit` uses the non-linear least squares method. This method decreases the difference between the measured release and the calculated release.

`CurveFit` gives the quality of the fit as the coefficient of determination:

$$
R^2 = 1 - \frac{\sum_i \left(M_i - \hat{M}_i\right)^2}{\sum_i \left(M_i - \bar{M}\right)^2}
$$

where:
- $M_i (mg)$ is the measured drug release at time $t_i$
- $\hat{M}_i (mg)$ is the drug release that the model calculates at time $t_i$
- $\bar{M} (mg)$ is the mean value of the measured drug release

An $R^2$ value near 1 shows a good fit. A low value shows that the model does not describe the data correctly.

You can fit all the models above. This table gives the name and the parameters of each model:

| Model        | `model_name`   | Parameters                 |
| ------------ | -------------- | -------------------------- |
| Zero-Order   | `zero_order`   | `M0`, `k0`                 |
| First-Order  | `first_order`  | `M0`, `k`                  |
| Higuchi      | `higuchi`      | `D`, `c0`, `cs`            |
| Weibull      | `weibull`      | `M`, `a`, `b`              |
| Hopfenberg   | `hopfenberg`   | `M`, `k0`, `c0`, `a0`, `n` |

The `fit` method starts from a value of 1 for each unknown parameter, and keeps all the parameters positive. If the result is not good, give your own values in the `initial_guess` and `bounds` arguments. These arguments contain only the unknown parameters, in the order of the table. The `get_result` method gives the result of the last fit again.

The `fit` method returns a `FitResult` object, where:
- `model_name` is the name of the fitted model
- `parameters` is a dictionary of the known and the calculated parameter values
- `r_squared` is the fit's coefficient of determination
- `model` is a model object with the calculated parameters. You can simulate and plot this model.

### Known Parameters
If you know the value of a parameter, give it in the `known_parameters` argument. Use the parameter names in the table above. `CurveFit` keeps these values constant and calculates only the other parameters.

⚠️ All the parameters of a model cannot always be calculated from a release profile. Some models use a group of parameters only as one product or one ratio. Many different combinations of these parameters give the same value for this group. Thus, the fit has no unique solution, and the calculated values can be different from the true physical values.  To prevent this problem, give at least **two** parameters of each group in the `known_parameters` argument. Two of the models have such a group:
- The Higuchi model uses $D$, $c_0$ and $c_s$ only in the product $D(2c_0 - c_s)c_s$. Give two of these three parameters, for example `known_parameters={"c0": 1, "cs": 0.5}`.
- The Hopfenberg model uses $k_0$, $c_0$ and $a_0$ only in the ratio $\frac{k_0}{c_0 a_0}$. Give two of these three parameters, and also the geometry factor $n$, for example `known_parameters={"c0": 0.0374, "a0": 3.51, "n": 2}`. Also, $n$ is always assumed to be known because the user should know the shape of the device.

### Applications
1. Analysis of dissolution test data
2. Estimation of the release rate constants
3. Selection of the correct model for an experiment
4. Comparison of different formulations
5. Quality control of production batches

## Usage
### Zero-Order Model
```python
from drux import ZeroOrderModel
model = ZeroOrderModel(k0=0.1, M0=0)
model.simulate(duration=1000, time_step=10)
model.plot(show=True)
```
<img src="https://github.com/openscilab/drux/raw/main/otherfiles/zero_order_plot.png" alt="Zero-order Plot">

### First-Order Model
```python
from drux import FirstOrderModel
model = FirstOrderModel(k=0.003, M0=0.1)
model.simulate(duration=1000, time_step=10)
model.plot(show=True)
```
<img src="https://github.com/openscilab/drux/raw/main/otherfiles/first_order_plot.png" alt="First-order Plot">

### Higuchi Model
```python
from drux import HiguchiModel
model = HiguchiModel(D=1e-6, c0=1, cs=0.5)
model.simulate(duration=1000, time_step=10)
model.plot(show=True)
```
<img src="https://github.com/openscilab/drux/raw/main/otherfiles/higuchi_plot.png" alt="Higuchi Plot">

### Weibull Model

```python
from drux import WeibullModel
model = WeibullModel(M=1, a=0.095, b=0.7)
model.simulate(duration=100, time_step=1)
model.plot(show=True)
```
<img src="https://github.com/openscilab/drux/raw/main/otherfiles/weibull_plot.png" alt="Weibull Plot">

### Hopfenberg Model

```python
from drux import HopfenbergModel
model = HopfenbergModel(M=1, k0=0.00067, c0=0.0374, a0=3.51, n=2)
model.simulate(duration=100, time_step=1)
model.plot(show=True)
```
<img src="https://github.com/openscilab/drux/raw/main/otherfiles/hopfenberg_plot.png" alt="Hopfenberg Plot">

### Curve Fitting

```python
from drux import CurveFit
time = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
release_profile = [0, 0.35, 0.48, 0.57, 0.63, 0.69, 0.73, 0.77, 0.80, 0.83, 0.85]
fitter = CurveFit(model_name="weibull", time=time, release_profile=release_profile)
result = fitter.fit()
print(result.parameters)
print(result.r_squared)
result.model.simulate(duration=100, time_step=1)
result.model.plot(show=True)
```
<img src="https://github.com/openscilab/drux/raw/main/otherfiles/curve_fit_plot.png" alt="Curve Fit Plot">


## Issues & bug reports

Just fill an issue and describe it. We'll check it ASAP! or send an email to [drux@openscilab.com](mailto:drux@openscilab.com "drux@openscilab.com"). 

- Please complete the issue template

You can also join our discord server

<a href="https://discord.gg/8Rf6bGBtse">
  <img src="https://img.shields.io/discord/1064533716615049236.svg?style=for-the-badge" alt="Discord Channel">
</a>


## References
<blockquote>1- T. Higuchi, "Rate of release of medicaments from ointment bases containing drugs in suspension," <i>Journal of Pharmaceutical Sciences</i>, vol. 50, no. 10, pp. 874–875, 1961.</blockquote>
<blockquote>2- D. R. Paul, "Elaborations on the Higuchi model for drug delivery," <i>International Journal of Pharmaceutics</i>, vol. 418, no. 1, pp. 13–17, 2011.</blockquote>
<blockquote>3- R. T. Medarametla, K. V. Gopaiah, J. N. Suresh Kumar, G. Anand Babu, M. Shaggir, G. Raghavendra, D. Naveen Reddy, and B. Venkamma, "Drug Release Kinetics and Mathematical Models," <i>International Journal of Science and Research Methodology</i>, vol. 27, no. 9, pp. 12–19, Sep. 2024.</blockquote>
<blockquote>4- R. Vaju and K. V. Murthy, "Development and validation of new discriminative dissolution method for carvedilol tablets," <i>Indian Journal of Pharmaceutical Sciences</i>, vol. 73, no. 5, pp. 527–536, Sep. 2011.</blockquote>
<blockquote>5- S. Dash, "Kinetic modeling on drug release from controlled drug delivery systems," <i>Acta Poloniae Pharmaceutica</i>, 2010.</blockquote>
<blockquote>6- K. H. Ramteke, P. A. Dighe, A. R. Kharat, S. V. Patil, <i>Mathematical models of drug dissolution: A review</i>, <i>Sch. Acad. J. Pharm.</i>, vol. 3, no. 5, pp. 388-396, 2014.</blockquote>
<blockquote>7- C. Corsaro, G. Neri, A. M. Mezzasalma, and E. Fazio, "Weibull modeling of controlled drug release from Ag-PMA nanosystems," <i>Polymers</i>, vol. 13, no. 17, p. 2897, 2021.</blockquote>
<blockquote>8- H. B. Hopfenberg, "Controlled release from erodible slabs, cylinders, and spheres," in <i>Controlled Release Polymeric Formulations</i>, ACS Symposium Series, vol. 33, pp. 26–32, 1976.</blockquote>
<blockquote>9- H. V. Chavda, M. S. Patel, and C. N. Patel, "Preparation and in vitro evaluation of guar gum based triple-layer matrix tablet of diclofenac sodium," <i>Research in Pharmaceutical Sciences</i>, vol. 7, no. 1, pp. 57–64, Jan. 2012.</blockquote>


## Show your support
### Star this repo

Give a ⭐️ if this project helped you!

### Donate to our project
If you do like our project and we hope that you do, can you please support us? Our project is not and is never going to be working for profit. We need the money just so we can continue doing what we do ;-) .			

<a href="https://openscilab.com/#donation" target="_blank"><img src="https://github.com/openscilab/drux/raw/main/otherfiles/donation.png" height="90px" width="270px" alt="Drux Donation"></a>
