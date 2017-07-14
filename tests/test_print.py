from rbconvert.objects import *
from nose.tools import raises

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


class TestPrint:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        cls.num_bond = Bond(3)
        cls.wild_bond = Bond(-1, w=True)
        cls.any_bond = Bond(-3, a=True)

        cls.s0 = Site('site0', 0, s='state', b=cls.num_bond)
        cls.s1 = Site('site1', 0, b=cls.any_bond)
        cls.s2 = Site('site2', 1, b=cls.wild_bond)
        cls.s3 = Site('site0', 0, s='state')

        cls.sd0 = SiteDef('site0', ['state', 'state2'])
        cls.sd1 = SiteDef('site1', [])
        cls.sd2 = SiteDef('site2', [])
        cls.sd3 = SiteDef('a', ['0', '1', '2'])
        cls.sd4 = SiteDef('b', [])
        cls.sd5 = SiteDef('c', [])

        cls.md0 = MoleculeDef('Test0', [cls.sd0, cls.sd1, cls.sd2],
                              {'site0': 'site0', 'site1': 'site1', 'site2': 'site2'})
        cls.md1 = MoleculeDef('Test1', [cls.sd1], {'site1': 'site1'})
        cls.md2 = MoleculeDef('Test2', [cls.sd3, cls.sd3, cls.sd3, cls.sd3, cls.sd4, cls.sd5, cls.sd5],
                              {'a0': 'a', 'a1': 'a', 'a2': 'a', 'a3': 'a', 'b': 'b', 'c0': 'c', 'c1': 'c'}, hss=True)
        cls.md3 = MoleculeDef('A', [cls.sd3, cls.sd3, cls.sd3], {'a0': 'a', 'a1': 'a', 'a2': 'a'})
        cls.md4 = MoleculeDef('B', [cls.sd4, cls.sd4], {'b0': 'b', 'b1': 'b'})

        cls.m0 = Molecule('Test0', [cls.s0, cls.s2])
        cls.m1 = Molecule('Test1', [cls.s1])
        cls.m2 = Molecule('Test2', [Site('a', 0)])
        cls.m3 = Molecule('Test2', [Site('a', 0), Site('a', 1)])
        cls.m4 = Molecule('Test2', [Site('a', 0), Site('a', 1, b=Bond(1)), Site('c', 2)])
        cls.m5 = Molecule('Test2', [Site('a', 0), Site('a', 1, b=Bond(1)), Site('a', 2), Site('b', 3)])
        cls.m6 = Molecule('Test2', [Site('a', 0), Site('a', 1, b=Bond(1)), Site('a', 2, Bond(2)), Site('b', 3)])
        cls.m7 = Molecule('A', [Site('a', 0, b=Bond(-1, w=True)), Site('a', 1)])
        cls.m8 = Molecule('B', [Site('b', 0)])
        cls.m9 = Molecule('A', [Site('a', 0, b=Bond(-1, w=True)), Site('a', 1, b=Bond(1))])
        cls.m10 = Molecule('B', [Site('b', 0, b=Bond(1))])

        cls.p0 = CPattern([cls.m0, Molecule('Test0', [cls.s0])])
        cls.p1 = CPattern([Molecule('Test0', [cls.s3, cls.s2]), Molecule('Test0', [cls.s3])])
        cls.p2 = CPattern([Molecule('A', [])])
        cls.p3 = CPattern([Molecule('B', [])])
        cls.p4 = CPattern([cls.m2, cls.m0])
        cls.p5 = CPattern([cls.m7])
        cls.p6 = CPattern([cls.m8])
        cls.p7 = CPattern([cls.m9, cls.m10])

        cls.i0 = InitialCondition(cls.p0, 10)
        # implement functionality to print initial condition as kappa/bngl expression
        cls.i1 = InitialCondition(cls.p0, Expression(['x', '+', '10']), False)

        cls.par0 = Parameter('rate', 1e6)

        expr0 = Expression(['ln', '(', '10', ')', '+', 'x', '-', '356'])
        cls.func0 = Function('rate_expr', expr0)

        cls.rate0 = Rate(3)
        cls.rate1 = Rate(expr0)
        cls.rate2 = Rate('rate')
        cls.rate3 = Rate(5, intra=True)

        cls.rule0 = Rule([cls.p2], [cls.p3], cls.rate0)
        cls.rule1 = Rule([cls.p2], [cls.p3], cls.rate1, False, cls.rate2)
        cls.rule2 = Rule([cls.p2], [cls.p3], cls.rate0, True, cls.rate2)
        cls.rule3 = Rule([cls.p5, cls.p6], [cls.p7], cls.rate0)

        cls.obs0 = Observable("Obs0", [cls.p3], 'm')
        cls.obs1 = Observable("Obs1", [cls.p2, cls.p3], 's')
        cls.obs2 = Observable("Obs2", [CPattern([cls.m2])], 'm')

    @classmethod
    def teardown_class(cls):
        pass

    def test_valid_bonds(self):
        assert self.num_bond.wild == False
        assert self.num_bond.any == False
        assert self.wild_bond.write_as_bngl() == r'!+'
        assert self.any_bond.write_as_bngl() == r'!?'
        assert self.wild_bond.write_as_kappa() == r'!_'
        assert self.any_bond.write_as_kappa() == r'?'

    @raises(AssertionError)
    def test_invalid_bond_0(self):
        Bond(-3)

    @raises(ValueError)
    def test_invalid_bond_1(self):
        Bond('hello')

    def test_sites(self):
        assert self.s0.write_as_bngl() == r'site0~state!3'
        assert self.s0.write_as_kappa() == r'site0~state!3'
        assert self.s1.write_as_bngl() == r'site1!?'
        assert self.s1.write_as_kappa() == r'site1?'

    def test_molec_def(self):
        assert self.md0.write_as_bngl() == r'Test0(site0~state~state2,site1,site2)'
        assert self.md0.write_as_kappa() == r'%agent: Test0(site0~state~state2,site1,site2)'
        assert self.md2.write_as_bngl() == r'Test2(a~0~1~2,a~0~1~2,a~0~1~2,a~0~1~2,b,c,c)'
        assert self.md2.write_as_kappa() == r'%agent: Test2(a0~0~1~2,a1~0~1~2,a2~0~1~2,a3~0~1~2,b,c0,c1)'

    def test_molecules(self):
        assert self.m0.write_as_bngl() == r'Test0(site0~state!3,site2!+)'
        assert self.m0.write_as_kappa() == r'Test0(site0~state!3,site2!_)'
        assert self.m1.write_as_bngl() == r'Test1(site1!?)'
        assert self.m1.write_as_kappa() == r'Test1(site1?)'

    def test_site_renaming(self):
        km2 = self.m2.convert(self.md2)
        km3 = self.m3.convert(self.md2)
        km4 = self.m4.convert(self.md2)
        km5 = self.m5.convert(self.md2)
        km6 = self.m6.convert(self.md2)
        assert len(km2) == 4    # 4 choose 1
        assert len(km3) == 48   # (4 choose 2) * 2^2 * 2^1
        assert len(km4) == 192  # km3 * 2 * 2^1
        assert len(km5) == 24   # (4 choose 2) * (2 choose 1) * 2
        assert len(km6) == 144  # (4 choose 1) * (3 choose 1) * (2 choose 1) * 3 * 2

    def test_patterns(self):
        assert self.p0.write_as_bngl() == r'Test0(site0~state!3,site2!+).Test0(site0~state!3)'
        assert self.p0.write_as_kappa() == r'Test0(site0~state!3,site2!_),Test0(site0~state!3)'
        assert self.p4.write_as_bngl() == r'Test2(a).Test0(site0~state!3,site2!+)'
        kp4 = self.p4.convert([self.md0, self.md2])
        assert len(kp4) == 4
        print sorted([k.write_as_kappa() for k in kp4])
        assert sorted([k.write_as_kappa() for k in kp4]) == [r'Test2(a0),Test0(site0~state!3,site2!_)',
                                                             r'Test2(a1),Test0(site0~state!3,site2!_)',
                                                             r'Test2(a2),Test0(site0~state!3,site2!_)',
                                                             r'Test2(a3),Test0(site0~state!3,site2!_)']

    def test_init_conditions(self):
        assert self.i0.write_as_bngl() == 'Test0(site0~state!3,site2!+).Test0(site0~state!3) 10'
        assert self.i0.write_as_kappa() == r'%init: 10 Test0(site0~state!3,site2!_),Test0(site0~state!3)'
        assert self.i1.write_as_bngl() == 'Test0(site0~state!3,site2!+).Test0(site0~state!3) x+10'
        assert self.i1.write_as_kappa() == r"%init: 'x'+10 Test0(site0~state!3,site2!_),Test0(site0~state!3)"

    # TODO write more to check function map functionality
    def test_pars_and_funcs(self):
        assert self.par0.write_as_bngl() == r'rate 1000000.0'
        assert self.par0.write_as_kappa() == r"%var: 'rate' 1000000.0"
        assert self.func0.write_as_bngl() == r'rate_expr()=ln(10)+x-356'
        assert self.func0.write_as_kappa() == r"%obs: 'rate_expr' [log](10)+'x'-356"
        assert self.func0.write_as_kappa(as_obs=False) == r"%var: 'rate_expr' [log](10)+'x'-356"

    def test_rate(self):
        assert self.rate0.write_as_bngl() == r'3'
        assert self.rate0.write_as_kappa() == r'3'
        assert self.rate1.write_as_bngl() == r'ln(10)+x-356'
        assert self.rate1.write_as_kappa() == r"[log](10)+'x'-356"
        assert self.rate2.write_as_bngl() == r'rate'
        assert self.rate2.write_as_kappa() == r"'rate'"
        assert self.rate3.write_as_bngl() == r'5'
        assert self.rate3.write_as_kappa() == r'0 {5}'

    def test_rule(self):
        assert self.rule0.write_as_bngl() == r'A() -> B() 3'
        assert self.rule0.write_as_kappa() == r'A() -> B() @ 3'
        assert self.rule1.write_as_bngl() == r'A() -> B() ln(10)+x-356'
        assert self.rule1.write_as_kappa() == r"A() -> B() @ [log](10)+'x'-356"
        assert self.rule2.write_as_bngl() == r'A() <-> B() 3,rate'
        assert self.rule2.write_as_kappa() == r"A() <-> B() @ 3,'rate'"

    # TODO CHECK this
    def test_rule_expansion(self):
        x = self.rule3.convert([self.md3, self.md4], [self.md3, self.md4])
        print x
        print len(x)
        assert len(self.rule3.convert([self.md3, self.md4], [self.md3, self.md4])) == 18

    def test_obs(self):
        assert self.obs0.write_as_bngl() == r'Molecules Obs0 B()'
        assert self.obs0.write_as_kappa([self.md4]) == r"%obs: 'Obs0' |B()|"
        assert self.obs1.write_as_bngl() == r'Species Obs1 A() B()'
        assert self.obs1.write_as_kappa([self.md4, self.md3]) == r"%obs: 'Obs1' |A()|+|B()|"
        assert self.obs2.write_as_bngl() == r'Molecules Obs2 Test2(a)'
        assert self.obs2.write_as_kappa([self.md2]) == r"%obs: 'Obs2' |Test2(a0)|+|Test2(a1)|+|Test2(a2)|+|Test2(a3)|"

    @raises(Exception)
    def test_invalid_obs(self):
        Observable("InvalidType", [self.p2], 'f')
