# -*- coding: utf-8 -*-
"""Setup module."""
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


def get_requires() -> list:
    """Read requirements.txt."""
    requirements = open("requirements.txt", "r").read()
    return list(filter(lambda x: x != "", requirements.split()))


def read_description() -> str:
    """Read README.md and CHANGELOG.md."""
    try:
        with open("README.md") as r:
            description = "\n"
            description += r.read()
        with open("CHANGELOG.md") as c:
            description += "\n"
            description += c.read()
        return description
    except Exception:
        return '''Drux is a Python-based framework for simulating drug release profiles using mathematical models.
                It offers a reproducible and extensible platform to model, analyze, and visualize time-dependent drug release behavior, making it ideal for pharmaceutical research and development.
                By combining simplicity with scientific rigor, Drux provides a robust foundation for quantitative analysis of drug delivery kinetics.
               '''


setup(
    name='drux',
    packages=[
        'drux', ],
    version='0.3',
    description='Drux: Drug Release Analysis Framework',
    long_description=read_description(),
    long_description_content_type='text/markdown',
    author='Drux Development Team',
    author_email='drux@openscilab.com',
    url='https://github.com/openscilab/drux',
    download_url='https://github.com/openscilab/drux/tarball/v0.3',
    keywords="drug-release drug-delivery mathematical-modeling simulation kinetics",
    project_urls={
            'Source': 'https://github.com/openscilab/drux',
    },
    install_requires=get_requires(),
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Manufacturing',
        'Intended Audience :: Science/Research',
        'Topic :: Education',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    license='MIT',
)
