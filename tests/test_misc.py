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

        cls.sd0 = objects.SiteDef('DNP')
        cls.sd1 = objects.SiteDef('Fab')
        cls.md0 = objects.MoleculeDef('A', [objects.SiteDef('x'), objects.SiteDef('y'), objects.SiteDef('s'),
                                            objects.SiteDef('t')], {'x': 'x', 'y': 'y', 's': 's', 't': 't'})
        cls.md1 = objects.MoleculeDef('B', [objects.SiteDef('y'), objects.SiteDef('s')], {'y': 'y', 's': 's'})
        cls.md2 = objects.MoleculeDef('BSA', [cls.sd0, cls.sd0, cls.sd0, cls.sd0],
                                      {'DNP0': 'DNP', 'DNP1': 'DNP', 'DNP2': 'DNP', 'DNP3': 'DNP'})
        cls.md3 = objects.MoleculeDef('IgE', [cls.sd1, cls.sd1], {'Fab0': 'Fab', 'Fab2': 'Fab'})
        cls.md4 = objects.MoleculeDef('C', [], {})
        cls.md5 = objects.MoleculeDef('D', [], {})

        cls.s0 = objects.Site("s", 0, b=objects.Bond(1))
        cls.s1 = objects.Site("t", 1, b=objects.Bond(2))
        cls.s2 = objects.Site("s", 0, b=objects.Bond(2))
        cls.s3 = objects.Site("t", 0, b=objects.Bond(2))

        cls.m0 = objects.Molecule('A', [], cls.md0)
        cls.m1 = objects.Molecule('B', [], cls.md1)
        cls.m2 = objects.Molecule('A', [cls.s0], cls.md0)
        cls.m3 = objects.Molecule('A', [cls.s0, cls.s1], cls.md0)
        cls.m4 = objects.Molecule('B', [cls.s2], cls.md1)
        cls.m5 = objects.Molecule('A', [objects.Site('x', 0, b=objects.Bond(1))], cls.md0)
        cls.m6 = objects.Molecule('B', [objects.Site('y', 0, b=objects.Bond(1))], cls.md1)
        cls.m7 = objects.Molecule('A', [objects.Site('x', 0)], cls.md0)
        cls.m8 = objects.Molecule('B', [objects.Site('y', 0)], cls.md1)
        cls.m9 = objects.Molecule('C', [], cls.md4)
        cls.m10 = objects.Molecule('D', [], cls.md5)

        cls.pattern = 'A(x!1).A(x!1,y!2).B(y!2)'
        cls.pattern2 = 'BSA(DNP!+,DNP!+,DNP!1,DNP).IgE(Fab!1,Fab!2).BSA(DNP!2,DNP!3).IgE(Fab!3,Fab)'
        cls.pattern3 = 'A(x!1,y!2).A(x!2,y!3).A(x!3,y!4).A(x!4,y!1)'
        cls.pattern4 = 'A(x!4).A(x!4,y!1).B(y!1)'
        cls.pattern5 = 'IgE(Fab!4,Fab!3).BSA(DNP!4,DNP!+,DNP,DNP!+).BSA(DNP!3,DNP!2).IgE(Fab,Fab!2)'

        cls.rate0 = objects.Rate(3)
        cls.rate1 = objects.Rate('rate')
        cls.rate2 = objects.Rate(objects.Expression(['rate', '/', '2']))

        cls.sd4 = objects.SiteDef('b', [])
        cls.md4 = objects.MoleculeDef('B', [cls.sd4, cls.sd4], {'b0': 'b', 'b1': 'b'})
        cls.m8 = objects.Molecule('B', [objects.Site('b', 0)], cls.md4)
        cls.m10 = objects.Molecule('B', [objects.Site('b', 0, b=objects.Bond(1))], cls.md4)
        cls.p6 = objects.CPattern([cls.m8])
        cls.p8 = objects.CPattern([cls.m10, cls.m10])
        cls.rule6 = objects.Rule(objects.CPatternList([cls.p6, cls.p6]), objects.CPatternList([cls.p8]), cls.rate0)

    @classmethod
    def teardown_class(cls):
        pass

    def test_bond_equality(self):
        assert self.bond0 == objects.Bond(1)
        assert self.bond0 != objects.Bond(2)
        assert self.bond0 != self.bond1
        assert self.bond1 == self.bond2

    def test_site_equality(self):
        assert self.s0 != self.s2
        assert self.s1 == self.s3

    def test_pattern_equality(self):
        pp = readers.BNGLReader.parse_cpattern(self.pattern, [self.md0, self.md1])
        pp4 = readers.BNGLReader.parse_cpattern(self.pattern4, [self.md0, self.md1])
        assert pp.is_isomorphic(pp4)
        pp2 = readers.BNGLReader.parse_cpattern(self.pattern2, [self.md2, self.md3])
        pp5 = readers.BNGLReader.parse_cpattern(self.pattern5, [self.md2, self.md3])
        assert pp2.is_isomorphic(pp5)

    def test_molecule_ordering(self):
        assert sorted([self.m1, self.m0]) == [self.m0, self.m1]

    def test_graph_builder(self):
        graph = readers.BNGLReader.parse_cpattern(self.pattern, [self.md0, self.md1])._build_graph()
        assert len(graph.nodes()) == 3
        assert len(graph.edges()) == 2
        assert graph.node[1]['name'] == 'A:xb_yb'
        assert graph[0][1]['name'] == 'x-x'

        graph2 = readers.BNGLReader.parse_cpattern(self.pattern2, [self.md2, self.md3])._build_graph()
        assert len(graph2.nodes()) == 4
        assert len(graph2.edges()) == 3

    def test_automorphism_counter(self):
        assert readers.BNGLReader.parse_cpattern(self.pattern, [self.md0, self.md1]).automorphisms() == 1
        assert readers.BNGLReader.parse_cpattern(self.pattern3, [self.md0, self.md1]).automorphisms() == 8
        assert objects.CPattern([objects.PlaceHolderMolecule()]).automorphisms() == 1

    @raises(rbexceptions.NotConvertedException)
    def test_not_converted(self):
        readers.BNGLReader.parse_cpattern(self.pattern2, [self.md2, self.md3, self.md0]).automorphisms()

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

    def test_contains_variable(self):
        assert not self.rate0.contains_variable('rate')
        assert self.rate1.contains_variable('rate')
        assert self.rate2.contains_variable('rate')

    def test_factors(self):
        print self.rule6._unique_reactant_indices()
        assert self.rule6._unique_reactant_indices() == {0: 2}
        assert self.rule6.rate_factor(True) == 0.5
        assert self.rule6.rate_factor(False) == 2.0