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
- Run `pip install drux==0.3`
### Source code
- Download [Version 0.3](https://github.com/openscilab/drux/archive/v0.3.zip) or [Latest Source](https://github.com/openscilab/drux/archive/dev.zip)
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

## Show your support
### Star this repo

Give a ⭐️ if this project helped you!

### Donate to our project
If you do like our project and we hope that you do, can you please support us? Our project is not and is never going to be working for profit. We need the money just so we can continue doing what we do ;-) .			

<a href="https://openscilab.com/#donation" target="_blank"><img src="https://github.com/openscilab/drux/raw/main/otherfiles/donation.png" height="90px" width="270px" alt="Drux Donation"></a>
