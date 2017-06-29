from setuptools import setup

setup(name='rbconvert',
      version='0.1',
      description='Conversion of BNGL models to Kappa models',
      author='Ryan Suderman',
      package=['rbconvert'],
      install_requires=['pyparsing', 'deepdiff', 'networkx', 'nose'])