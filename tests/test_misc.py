from nose.tools import raises
from .context import objects
from .context import readers
from .context import rbexceptions
from .context import utils


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
        cls.m2 = objects.Molecule('A', [objects.Site("s", 0, b=objects.Bond(1))])
        cls.m3 = objects.Molecule('A', [objects.Site("s", 0, b=objects.Bond(1)), objects.Site("t", 1, b=objects.Bond(2))])
        cls.m4 = objects.Molecule('B', [objects.Site("s", 0, b=objects.Bond(2))])

        cls.pattern = 'A(x!1).A(x!1,y!2).B(y!2)'
        cls.pattern2 = 'BSA(DNP!+,DNP!+,DNP!1,DNP).IgE(Fab!1,Fab!2).BSA(DNP!2,DNP!3).IgE(Fab!3,Fab)'
        cls.pattern3 = 'A(x!1,y!2).A(x!2,y!3).A(x!3,y!4).A(x!4,y!1)'

        cls.pattern4 = 'A(x!1),A(x!1,y!2),B(x!2),C(y)'
        cls.pattern4 = 'A(x!1),A(x!1,y!2),C(y!2)'

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

    def test_bound_to(self):
        assert not self.m3.bound_to(self.bond0)
        assert self.m3.bound_to(self.m2)
        assert not self.m3.bound_to(self.m0)

    def test_adj_list(self):
        al0 = utils.build_adj_list([self.m2, self.m3, self.m4])
        assert al0 == [[1], [0, 2], [1]]
        al1 = utils.build_adj_list([self.m2, self.m3, self.m0, self.m4])
        assert al1 == [[1], [0, 3], [], [1]]

    def test_get_comps(self):
        cmps0 = utils.get_connected_components([self.m2, self.m3, self.m4])
        assert cmps0 == [[self.m2, self.m3, self.m4]]
        cmps1 = utils.get_connected_components([self.m2, self.m3, self.m0, self.m4])
        assert cmps1 == [[self.m2, self.m3, self.m4], [self.m0]]
