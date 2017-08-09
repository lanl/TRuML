from setuptools import setup

setup(name='truml',
      version='0.1',
      description='Interconversion between BNGL models and Kappa models',
      author='Ryan Suderman',
      package=['truml'],
      install_requires=['pyparsing', 'deepdiff', 'networkx', 'nose'])