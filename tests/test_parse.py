import os
import sys
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))

from rbconvert.parse_bngl import *

from nose.tools import raises

def test_valid_bonds():
	num_bond = Bond('123')
	wild_bond = Bond(-1,w=True)
	any_bond = Bond(-3,a=True)

	assert num_bond.wild == False
	assert num_bond.any == False
	assert wild_bond.write_as_bngl() == r'!+'
	assert any_bond.write_as_bngl() == r'!?'
	assert wild_bond.write_as_kappa() == r'!_'
	assert any_bond.write_as_kappa() == r'?'

@raises(AssertionError)
def test_invalid_bonds():
	Bond(-3)