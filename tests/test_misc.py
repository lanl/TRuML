import os
import sys
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))

from rbconvert.parse_bngl import *

from nose.tools import raises

class TestMisc:

	@classmethod
	def setup_class(cls):
		cls.bond0 = Bond(1)
		cls.bond1 = Bond(-1,w=True)
		cls.bond2 = Bond(-2,w=True)

	@classmethod
	def teardown_class(cls):
		pass

	def test_bond_equality(self):
		assert self.bond0 == Bond(1)
		assert self.bond0 != Bond(2)
		assert self.bond0 != self.bond1
		assert self.bond1 == self.bond2