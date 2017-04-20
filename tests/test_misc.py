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

		cls.pattern = 'A(x!1).A(x!1,y!2).B(y!2)'
		cls.pattern2 = 'BSA(DNP!+,DNP!+,DNP!1,DNP).IgE(Fab!1,Fab!2).BSA(DNP!2,DNP!3).IgE(Fab!3,Fab)'

	@classmethod
	def teardown_class(cls):
		pass

	def test_bond_equality(self):
		assert self.bond0 == Bond(1)
		assert self.bond0 != Bond(2)
		assert self.bond0 != self.bond1
		assert self.bond1 == self.bond2

	def test_graph_builder(self):
		graph = BNGLReader.parse_cpattern(self.pattern)._build_graph()
		assert len(graph.nodes()) == 7
		assert len(graph.edges()) == 8 # bond edges count twice
		assert graph.node[1]['state'] == ''
		assert graph.node[4]['bond'] == 'b'

		graph2 = BNGLReader.parse_cpattern(self.pattern2)._build_graph()
		assert len(graph2.nodes()) == 14
		assert len(graph2.edges()) == 16
		assert graph2.node[1]['bond'] == 'w'