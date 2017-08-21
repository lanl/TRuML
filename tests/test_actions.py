from .context import objects
from .context import readers

class TestAction:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        cls.wild_bond = objects.Bond(-1, w=True)

        cls.s2 = objects.Site('site2', 1, b=cls.wild_bond)
        cls.s3 = objects.Site('site0', 0, s='state')

        cls.sd0 = objects.SiteDef('site0', ['state', 'state2'])
        cls.sd1 = objects.SiteDef('site1', [])
        cls.sd2 = objects.SiteDef('site2', [])
        cls.sd3 = objects.SiteDef('a', ['0', '1', '2'])
        cls.sd4 = objects.SiteDef('b', [])
        cls.sd5 = objects.SiteDef('c', [])

        cls.md2 = objects.MoleculeDef('Test2', [cls.sd3, cls.sd3, cls.sd3, cls.sd3, cls.sd4, cls.sd5, cls.sd5],
                                      {'a0': 'a', 'a1': 'a', 'a2': 'a', 'a3': 'a', 'b': 'b', 'c0': 'c', 'c1': 'c'},
                                      hss=True)
        cls.md3 = objects.MoleculeDef('A', [cls.sd3, cls.sd3, cls.sd3], {'a0': 'a', 'a1': 'a', 'a2': 'a'})
        cls.md4 = objects.MoleculeDef('B', [cls.sd4, cls.sd4], {'b0': 'b', 'b1': 'b'})
        cls.md5 = objects.MoleculeDef('C', [], {})

        cls.m7 = objects.Molecule('A', [objects.Site('a', 0, b=objects.Bond(-1, w=True)), objects.Site('a', 1)], cls.md3)
        cls.m8 = objects.Molecule('B', [objects.Site('b', 0)], cls.md4)
        cls.m9 = objects.Molecule('A', [objects.Site('a', 0, b=objects.Bond(-1, w=True)),
                                        objects.Site('a', 1, b=objects.Bond(1))], cls.md3)
        cls.m10 = objects.Molecule('B', [objects.Site('b', 0, b=objects.Bond(1))], cls.md4)

        cls.p2 = objects.CPattern([objects.Molecule('A', [], cls.md3)])
        cls.p3 = objects.CPattern([objects.Molecule('B', [], cls.md4)])
        cls.p4 = objects.CPattern([objects.Molecule('C', [], cls.md5)])
        cls.p5 = objects.CPattern([cls.m7])
        cls.p6 = objects.CPattern([cls.m8])
        cls.p7 = objects.CPattern([cls.m9, cls.m10])

        cls.rate0 = objects.Rate(3)

        cls.rule0 = objects.Rule([cls.p2], [cls.p3], cls.rate0)
        cls.rule1 = objects.Rule([cls.p5, cls.p6], [cls.p7], cls.rate0)
        cls.rule2 = objects.Rule([cls.p4, cls.p2, cls.p3], [cls.p3, cls.p4], cls.rate0)

    @classmethod
    def teardown_class(cls):
        pass

    def test_build_mol_map(self):
        mmap0 = readers.BNGLReader()._build_mol_map(self.rule0.lhs, self.rule0.rhs)
        assert mmap0[0] is None
        assert len(mmap0.keys()) == 1

        mmap1 = readers.BNGLReader()._build_mol_map(self.rule1.lhs, self.rule1.rhs)
        assert mmap1[0] == 0
        assert mmap1[1] == 1
        assert len(mmap1.keys()) == 2

        mmap2 = readers.BNGLReader()._build_mol_map(self.rule2.lhs, self.rule2.rhs)
        assert mmap2[0] == 1
        assert mmap2[1] is None
        assert mmap2[2] == 0

    def test_action_parse(self):
        pass

    def test_action_apply(self):
        pass