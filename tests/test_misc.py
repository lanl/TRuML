from nose.tools import raises
from .context import objects
from .context import readers
from .context import rbexceptions


class TestMisc:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        cls.bond0 = objects.Bond(1)
        cls.bond1 = objects.Bond(-1, w=True)
        cls.bond2 = objects.Bond(-2, w=True)

        cls.m0 = objects.Molecule('A', [])
        cls.m1 = objects.Molecule('B', [])

        cls.pattern = 'A(x!1).A(x!1,y!2).B(y!2)'
        cls.pattern2 = 'BSA(DNP!+,DNP!+,DNP!1,DNP).IgE(Fab!1,Fab!2).BSA(DNP!2,DNP!3).IgE(Fab!3,Fab)'
        cls.pattern3 = 'A(x!1,y!2).A(x!2,y!3).A(x!3,y!4).A(x!4,y!1)'

    @classmethod
    def teardown_class(cls):
        pass

    def test_bond_equality(self):
        assert self.bond0 == objects.Bond(1)
        assert self.bond0 != objects.Bond(2)
        assert self.bond0 != self.bond1
        assert self.bond1 == self.bond2

    def test_molecule_ordering(self):
        assert sorted([self.m1, self.m0]) == [self.m0, self.m1]

    def test_graph_builder(self):
        graph = readers.BNGLReader.parse_cpattern(self.pattern)._build_graph()
        assert len(graph.nodes()) == 3
        assert len(graph.edges()) == 2
        assert graph.node[1]['name'] == 'A:xb_yb'
        assert graph[0][1]['name'] == 'x-x'

        graph2 = readers.BNGLReader.parse_cpattern(self.pattern2)._build_graph()
        assert len(graph2.nodes()) == 4
        assert len(graph2.edges()) == 3

    def test_automorphism_counter(self):
        assert readers.BNGLReader.parse_cpattern(self.pattern).automorphisms() == 1
        assert readers.BNGLReader.parse_cpattern(self.pattern3).automorphisms() == 8

    @raises(rbexceptions.NotConvertedException)
    def test_not_converted(self):
        readers.BNGLReader.parse_cpattern(self.pattern2).automorphisms()
