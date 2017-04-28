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
		cls.pattern3 = 'A(x!1,y!2).A(x!2,y!3).A(x!3,y!4).A(x!4,y!1)'

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
		assert len(graph.nodes()) == 3
		assert len(graph.edges()) == 2
		assert graph.node[1]['name'] == 'A:xb_yb'
		assert graph[0][1]['name'] == 'x-x'

		graph2 = BNGLReader.parse_cpattern(self.pattern2)._build_graph()
		assert len(graph2.nodes()) == 4
		assert len(graph2.edges()) == 3

	def test_automorphism_counter(self):
		assert BNGLReader.parse_cpattern(self.pattern).automorphisms() == 1
		assert BNGLReader.parse_cpattern(self.pattern3).automorphisms() == 8

	@raises(NotConvertedException)
	def test_not_converted(self):
		BNGLReader.parse_cpattern(self.pattern2).automorphisms()