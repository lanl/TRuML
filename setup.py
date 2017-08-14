"""setup.py: setuptools control"""


import re
from setuptools import setup

version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('truml/truml.py').read(),
    re.M
    ).group(1)

with open('README.rst', 'rb') as f:
    long_desc = f.read().decode('utf-8')

setup(name='TRuML',
      packages=['truml'],
      entry_points={'console_scripts': ['truml = truml.truml:main']},
      version=version,
      description='Translation between BNGL models and Kappa models',
      long_description=long_desc,
      author='Ryan Suderman',
      author_email='ryants@lanl.gov',
      install_requires=['pyparsing', 'deepdiff', 'networkx', 'nose'])