from .context import objects
from .context import readers


class TestAction:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        cls.wild_bond = objects.Bond(-1, w=True)

        cls.s3 = objects.Site('site0', 0, s='state')

        cls.sd0 = objects.SiteDef('site0', ['state', 'state2'])
        cls.sd3 = objects.SiteDef('a', ['0', '1', '2'])
        cls.sd4 = objects.SiteDef('b', [])
        cls.sd5 = objects.SiteDef('c', [])

        cls.md2 = objects.MoleculeDef('Test2', [cls.sd3, cls.sd3, cls.sd3, cls.sd3, cls.sd4, cls.sd5, cls.sd5],
                                      {'a0': 'a', 'a1': 'a', 'a2': 'a', 'a3': 'a', 'b': 'b', 'c0': 'c', 'c1': 'c'},
                                      hss=True)
        cls.md3 = objects.MoleculeDef('A', [cls.sd3, cls.sd3, cls.sd3], {'a0': 'a', 'a1': 'a', 'a2': 'a'})
        cls.md4 = objects.MoleculeDef('B', [cls.sd4, cls.sd4], {'b0': 'b', 'b1': 'b'})
        cls.md5 = objects.MoleculeDef('C', [cls.sd0], {'site0': 'site0'})

        cls.m0 = objects.Molecule('C', [objects.Site('site0', 0, s='state')], cls.md5)

        cls.m7 = objects.Molecule('A', [objects.Site('a', 0, b=objects.Bond(-1, w=True)),
                                        objects.Site('a', 1)], cls.md3)
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
        cls.rule3 = objects.Rule([cls.p3, cls.p5, cls.p6], [cls.p7, cls.p4], cls.rate0)
        cls.rule4 = objects.Rule([cls.p7], [cls.p6], cls.rate0)

    @classmethod
    def teardown_class(cls):
        pass

    def test_interface_map(self):
        assert self.p2[0].interface_diff_map(self.p3[0]) == dict()
        idm0 = self.m7.interface_diff_map(self.m9)
        assert len(idm0) == 1
        assert 1 in idm0.keys()
        assert idm0[1] == (None, ('a', None, objects.Bond(1), True))

    def test_build_mol_map(self):
        r0_lhs_mols = [x for cp in self.rule0.lhs for x in cp.molecule_list]
        r0_rhs_mols = [x for cp in self.rule0.rhs for x in cp.molecule_list]
        mmap0 = readers.BNGLReader()._build_mol_map(r0_lhs_mols, r0_rhs_mols)
        assert mmap0[0] is None
        assert len(mmap0.keys()) == 1

        r1_lhs_mols = [x for cp in self.rule1.lhs for x in cp.molecule_list]
        r1_rhs_mols = [x for cp in self.rule1.rhs for x in cp.molecule_list]
        mmap1 = readers.BNGLReader()._build_mol_map(r1_lhs_mols, r1_rhs_mols)
        assert mmap1[0] == 0
        assert mmap1[1] == 1
        assert len(mmap1.keys()) == 2

        r2_lhs_mols = [x for cp in self.rule2.lhs for x in cp.molecule_list]
        r2_rhs_mols = [x for cp in self.rule2.rhs for x in cp.molecule_list]
        mmap2 = readers.BNGLReader()._build_mol_map(r2_lhs_mols, r2_rhs_mols)
        assert mmap2[0] == 1
        assert mmap2[1] is None
        assert mmap2[2] == 0

    def test_action_parse(self):
        a0 = readers.BNGLReader()._build_actions(self.rule0.lhs, self.rule0.rhs)
        assert len(a0) == 2
        assert isinstance(a0[0], objects.Degradation)
        assert isinstance(a0[1], objects.Synthesis)
        assert a0[0].mol_index == 0
        assert isinstance(a0[1].cpattern, objects.CPattern)
        assert len(a0[1].cpattern.molecule_list) == 1

        a1 = readers.BNGLReader()._build_actions(self.rule1.lhs, self.rule1.rhs)
        assert len(a1) == 2
        assert isinstance(a1[0], objects.Binding)
        assert isinstance(a1[1], objects.Binding)

        a2 = readers.BNGLReader()._build_actions(self.rule1.rhs, self.rule1.lhs)
        assert len(a2) == 2
        assert isinstance(a2[0], objects.Unbinding)
        assert isinstance(a2[1], objects.Unbinding)

        a3 = readers.BNGLReader()._build_actions(self.rule3.lhs, self.rule3.rhs)
        print a3
        assert len(a3) == 4

    def test_synth_action_apply(self):
        synth = objects.Synthesis(self.p3)
        rhs = synth.apply(self.rule0.lhs)
        assert len(rhs) == 2

    def test_statechange_action_apply(self):
        sc = objects.StateChange(0, 0, ('site0', 'state', 'state2'))
        rhs = sc.apply([objects.CPattern([self.m0])])
        assert len(rhs) == 1
        assert rhs[0][0].sites[0].state == 'state2'

    def test_degradation_action_apply(self):
        d = objects.Degradation(0)
        rhs = d.apply(self.rule4.lhs)
        assert len(rhs) == 1
        assert rhs[0][0].name == 'B'
