from .context import objects
from .context import utils


class TestAction:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        cls.wild_bond = objects.Bond(-1, w=True)

        cls.s3 = objects.Site('site0', 0, s='state')
        cls.s4 = objects.Site('site', 0, s='state')

        cls.sd0 = objects.SiteDef('site0', ['state', 'state2'])
        cls.sd1 = objects.SiteDef('site', ['state', 'state2'])
        cls.sd3 = objects.SiteDef('a', ['0', '1', '2'])
        cls.sd4 = objects.SiteDef('b', [])
        cls.sd5 = objects.SiteDef('c', [])

        cls.md2 = objects.MoleculeDef('Test2', [cls.sd3, cls.sd3, cls.sd3, cls.sd3, cls.sd4, cls.sd5, cls.sd5],
                                      {'a0': 'a', 'a1': 'a', 'a2': 'a', 'a3': 'a', 'b': 'b', 'c0': 'c', 'c1': 'c'},
                                      hss=True)
        cls.md3 = objects.MoleculeDef('A', [cls.sd3, cls.sd3, cls.sd3], {'a0': 'a', 'a1': 'a', 'a2': 'a'})
        cls.md4 = objects.MoleculeDef('B', [cls.sd4, cls.sd4], {'b0': 'b', 'b1': 'b'})
        cls.md5 = objects.MoleculeDef('C', [cls.sd0], {'site0': 'site0'})
        cls.md6 = objects.MoleculeDef('D', [cls.sd1, cls.sd1], {'site0': 'site', 'site1': 'site'})

        cls.m0 = objects.Molecule('C', [objects.Site('site0', 0, s='state')], cls.md5)
        cls.m1 = objects.Molecule('D', [cls.s4], cls.md6)
        cls.m2 = objects.Molecule('C', [objects.Site('site0', 0, s='state2', b=objects.Bond(1))], cls.md5)
        cls.m3 = objects.Molecule('D', [objects.Site('site', 0, s='state', b=objects.Bond(1))], cls.md6)
        cls.p0 = objects.CPattern([cls.m0])
        cls.p1 = objects.CPattern([cls.m1])
        cls.p8 = objects.CPattern([cls.m2, cls.m3])

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
        cls.rule2 = objects.Rule([cls.p4, cls.p2, cls.p3], [cls.p4, cls.p3], cls.rate0)
        cls.rule3 = objects.Rule(objects.CPatternList([cls.p3, cls.p5, cls.p6]), objects.CPatternList([cls.p7, cls.p4]), cls.rate0)
        cls.rule4 = objects.Rule(objects.CPatternList([cls.p7]), objects.CPatternList([cls.p6]), cls.rate0, delmol=True)
        cls.rule5 = objects.Rule(objects.CPatternList([cls.p0, cls.p1]), objects.CPatternList([cls.p8]), cls.rate0)

    @classmethod
    def teardown_class(cls):
        pass

    def test_interface_map(self):
        assert self.p2[0].interface_diff_map(self.p3[0]) == dict()
        idm0 = self.m7.interface_diff_map(self.m9)
        assert len(idm0) == 1
        s0 = objects.Site('a', 1)
        assert s0 in idm0.keys()
        assert idm0[s0] == (-1, objects.Bond(1))

    def test_build_mol_map(self):
        r0_lhs_mols = utils.flatten_pattern([self.p2])
        r0_rhs_mols = utils.flatten_pattern([self.p3])
        mmap0 = objects.Rule._build_mol_map(r0_lhs_mols, r0_rhs_mols)
        assert mmap0[0] is None
        assert len(mmap0.keys()) == 1

        r1_lhs_mols = utils.flatten_pattern([self.p5, self.p6])
        r1_rhs_mols = utils.flatten_pattern([self.p7])
        mmap1 = objects.Rule._build_mol_map(r1_lhs_mols, r1_rhs_mols)
        assert mmap1[0] == 0
        assert mmap1[1] == 1
        assert len(mmap1.keys()) == 2

        r2_lhs_mols = utils.flatten_pattern([self.p4, self.p2, self.p3])
        r2_rhs_mols = utils.flatten_pattern([self.p4, self.p3])
        mmap2 = objects.Rule._build_mol_map(r2_lhs_mols, r2_rhs_mols)
        assert mmap2[0] == 0
        assert mmap2[1] is None
        assert mmap2[2] == 1

    def test_action_parse(self):
        a0 = self.rule0._build_actions()
        assert len(a0) == 2
        assert isinstance(a0[1], objects.Degradation)
        assert isinstance(a0[0], objects.Synthesis)
        assert a0[1].mol_index == 0
        assert isinstance(a0[0].cpattern, objects.CPattern)
        assert len(a0[0].cpattern.molecule_list) == 1
        a1 = self.rule1._build_actions()
        assert len(a1) == 2
        assert isinstance(a1[0], objects.BondChange)
        assert isinstance(a1[1], objects.BondChange)
        a2 = self.rule1._build_actions(True)
        assert len(a2) == 2
        assert isinstance(a2[0], objects.BondChange)
        assert isinstance(a2[1], objects.BondChange)

        a3 = self.rule3._build_actions()
        assert len(a3) == 4

    def test_synth_action_apply(self):
        synth = objects.Synthesis(self.p3)
        rhs = synth.apply(self.rule0.lhs)
        assert len(rhs) == 1
        assert len(rhs[0]) == 3

    def test_statechange_action_apply(self):
        lhs = [objects.CPattern([self.m1])]
        conv_lhs = []
        for l in lhs:
            conv_lhs.append(l.convert())
        sc = objects.StateChange(0, self.s4, 'state2', self.m1.mdef)
        rhss = [sc.apply(c) for c in conv_lhs]
        assert len(rhss) == len(conv_lhs)
        for rhs in rhss:
            assert len(rhs) == 1
            assert rhs[0][0][0].sites[0].state == 'state2'

    def test_degradation_action_apply(self):
        d = objects.Degradation(0)
        rhs = d.apply(self.rule4.lhs)
        assert len(rhs) == 1
        assert rhs[0][0][0].name == 'B'

    def test_binding_action_apply(self):
        b0 = objects.BondChange(0, objects.Site('a', 1), objects.Bond(1), self.md3)
        b1 = objects.BondChange(1, objects.Site('b', 0), objects.Bond(1), self.md4)
        rhs0 = b0.apply(b1.apply(self.rule1.lhs)[0])
        assert len(rhs0) == 1
        assert rhs0[0][0][0].sites[1].bond == objects.Bond(1)
        assert rhs0[0][0][1].sites[0].bond == objects.Bond(1)
        rhs1 = b1.apply(b0.apply(self.rule1.lhs)[0])
        assert len(rhs1) == 1
        assert rhs1[0][0][0].sites[1].bond == objects.Bond(1)
        assert rhs1[0][0][1].sites[0].bond == objects.Bond(1)

    def test_unbinding_action_apply(self):
        u0 = objects.BondChange(0, objects.Site('a', 1, b=objects.Bond(1)), None, self.md3)
        u1 = objects.BondChange(1, objects.Site('b', 0, b=objects.Bond(1)), None, self.md4)
        rhs0 = u0.apply(u1.apply(self.rule1.rhs)[0])
        assert len(rhs0) == 1
        assert len(rhs0[0]) == 2
        assert rhs0[0][0][0].sites[1].bond is None
        assert rhs0[0][1][0].sites[0].bond is None
        rhs1 = u1.apply(u0.apply(self.rule1.rhs)[0])
        assert len(rhs1) == 1
        assert len(rhs1[0]) == 2
        assert rhs1[0][0][0].sites[1].bond is None
        assert rhs1[0][1][0].sites[0].bond is None
        m = objects.MultiAction([u0, u1])
        rhs2 = m.apply(self.rule1.rhs)
        assert len(rhs2) == 1
        assert len(rhs2[0]) == 2
        assert rhs2[0][0][0].sites[1].bond is None
        assert rhs2[0][1][0].sites[0].bond is None

    def test_simultaneous_bond_state_apply(self):
        bs0 = objects.BondAndStateChange(0, objects.Site('site0', 0, s='state'), 'state2', objects.Bond(1), self.md5)
        bs1 = objects.BondChange(1, self.s4, objects.Bond(1), self.md6)
        acts = objects.MultiAction([bs0, bs1])
        rhs = acts.apply(self.rule5.lhs)
        assert len(rhs) == 1
        assert len(rhs[0]) == 1
        assert rhs[0][0][0].sites[0].state == 'state2'
        assert rhs[0][0][0].sites[0].bond == objects.Bond(1)

    def test_deg_unbinding_parse_apply(self):
        acts = self.rule4._build_actions()
        assert len(acts) == 2
        rhs = acts.apply(self.rule4.lhs)
        assert len(rhs) == 1
        assert len(rhs[0]) == 1
        assert rhs[0][0][0].sites[0].bond is None
