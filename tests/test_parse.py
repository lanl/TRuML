import os
import sys
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))

from rbconvert.parse_bngl import *

from nose.tools import raises

def test_valid_bonds():
	num_bond = Bond(3)
	wild_bond = Bond(-1,w=True)
	any_bond = Bond(-3,a=True)

	assert num_bond.wild == False
	assert num_bond.any == False
	assert wild_bond.write_as_bngl() == r'!+'
	assert any_bond.write_as_bngl() == r'!?'
	assert wild_bond.write_as_kappa() == r'!_'
	assert any_bond.write_as_kappa() == r'?'

@raises(AssertionError)
def test_invalid_bond_0():
	Bond(-3)

@raises(ValueError)
def test_invalid_bond_1():
	Bond('hello')

def test_sites():
	s0 = Site('bound',s='state',b=Bond(0))
	s1 = Site('unbound',b=Bond(-1,a=True))

	assert s0.write_as_bngl() == r'bound~state!0'
	assert s0.write_as_kappa() == r'bound~state!0'
	assert s1.write_as_bngl() == r'unbound!?'
	assert s1.write_as_kappa() == r'unbound?'

def test_molec_def():
	m0 = MoleculeDef('Molec',{'site0':['a','b'],'site1':[]})

	assert m0.write_as_bngl() == r'Molec(site0~a~b,site1)'
	assert m0.write_as_kappa() == r'%agent: Molec(site0~a~b,site1)'