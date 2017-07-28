from .context import readers
from .context import objects


class TestParseKappa:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        cls.mdef0 = "%agent: M(x,y,z~0~1)"

        cls.mol0 = "M(x,y,z~0)"
        cls.mol1 = "M(x,y,z~0!_)"
        cls.mol2 = "M+_2-985798f(x?,y,z~0)"

        cls.cp0 = "A(),B(x!1),C(y!1,z~s?)"

        cls.init0 = "%init: 10 A(x)"
        cls.init1 = "%init: 10 + 'x' B(),C()"

        cls.expr0 = "10 + 'x'"
        cls.expr1 = "[log] 100 / [max] 10 100 - [int] 7.342"

        cls.rule0 = 'A(x),B(x) -> A(x!1),B(x!1) @ 1'
        cls.rule1 = "'rule' %s" % cls.rule0
        cls.rule2 = "A(a~0),B(y) <-> A(a~1),B(y) @ %s {1}, 0.1 {10}" % cls.expr0
        cls.rule3 = "A(x),B(x) <-> A(x!1),B(x!1) @ %s {0}, 0.01" % cls.expr1
        cls.rule4 = "A(x),B(x) <-> A(x!1),B(x!1) @ 1, 0.1"

    def test_rule_parse(self):
        rule0s = readers.KappaReader.parse_rule(self.rule0)
        assert len(rule0s) == 1
        assert not rule0s[0].rev
        rule1s = readers.KappaReader.parse_rule(self.rule1)
        assert rule1s[0].lhs == rule0s[0].lhs
        assert rule1s[0].label == 'rule'
        rule2s = readers.KappaReader.parse_rule(self.rule2)
        assert len(rule2s) == 4
        for r in rule2s:
            assert not r.rev
        rule3s = readers.KappaReader.parse_rule(self.rule3)
        assert len(rule3s) == 3
        rule4s = readers.KappaReader.parse_rule(self.rule4)
        assert len(rule4s) == 1
        assert rule4s[0].rev

    def test_cpattern_parse(self):
        pcp0 = readers.KappaReader.parse_cpattern(self.cp0)
        assert len(pcp0.molecule_list) == 3
        assert pcp0.molecule_list[2].sites[0].bond.num == 1
        assert not pcp0.molecule_list[0].sites
        assert pcp0.molecule_list[2].sites[1].state == 's'

    def test_init_parse(self):
        assert readers.KappaReader.parse_init(self.init0).write_as_kappa() == "%init: 10 A(x)"
        assert readers.KappaReader.parse_init(self.init1).write_as_kappa() == "%init: 10+'x' B(),C()"

    def test_eq_parse(self):
        assert readers.KappaReader.parse_alg_expr(self.expr0).asList() == ['10', '+', "'x'"]
        assert readers.KappaReader.parse_alg_expr(self.expr1).asList() == \
               ['[log]', '100', '/', '[max]', '10', '100', '-', '[int]', '7.342']

    def test_mdef_parse(self):
        assert readers.KappaReader.parse_mtype(self.mdef0).write_as_kappa() == "%agent: M(x,y,z~0~1)"

    def test_mol_parse(self):
        assert readers.KappaReader.parse_molecule(self.mol0).write_as_kappa() == "M(x,y,z~0)"
        pmol1 = readers.KappaReader.parse_molecule(self.mol1)
        assert pmol1.write_as_kappa() == "M(x,y,z~0!_)"
        assert pmol1.sites[2].bond.wild
        pmol2 = readers.KappaReader.parse_molecule(self.mol2)
        assert pmol2.write_as_kappa() == "M+_2-985798f(x?,y,z~0)"
        assert pmol2.sites[0].bond.any
        assert pmol2.name == "M+_2-985798f"

    def test_vars_parse(self):
        kr = readers.KappaReader()
        kr.lines = ["%var: 'a' 3", "%var: 'b' 3 + 'a'", "%var: 'c' |C(x!_,y~state?)|", "%var: 'd' |A()| + 'b'"]
        model = kr.parse()
        assert len(model.functions) == 1
        assert model.parameters[0].name == 'a'
        assert model.parameters[1].name == 'b'
        assert isinstance(model.parameters[1].value, objects.Expression)
        cmd = objects.MoleculeDef("C", [objects.SiteDef('x'), objects.SiteDef('y', ['state', 'state2'])], {'x': 'x', 'y': 'y'})
        assert model.observables[0].write_as_kappa([cmd]) == "%obs: 'c' |C(x!_,y~state?)|"
        assert len(model.observables) == 2
        print model.observables[0].name
        print model.observables[1].name
        assert model.observables[1].name == "anon_obs0"


class TestParseBNGL:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        cls.mdef0 = "Molecule(site0,site1~0~P~PP)"
        cls.mdef1 = "Mol(a,b~0~1,c~a~b,b~0~1,c~a~b)"
        cls.mol0 = "Mol(sa,sb!+,sc!3,sd~0!?)"
        cls.mol1 = "Mol(a,b~0,b~1)"
        cls.init0 = cls.mol0 + ' 100'
        cls.init1 = cls.mol0 + '\t(x+3)/k'
        cls.obs0 = "Molecules Mol0 " + cls.mol0
        cls.obs1 = "Species Mol1 " + cls.mol0
        cls.param0 = "kcat 1"
        cls.param1 = "kp=km/kd/(NA*V)"
        cls.rule0 = "A(a) + B(b)<->A(a!1).B(b!1) kp,km"
        # intermolecular rate
        cls.rule1 = "A(a~r)+B(b,c!1).C(c!1) -> A(a~r!2).B(b!2,c!1).C(c!1) kp / log10(10)"
        # intramolecular rule
        cls.rule2 = "A(a~s).B(b,c!1).C(c!1) -> A(a~s!2).B(b!2,c!1).C(c!1) kp/log10(10)"
        cls.rule3 = "K(s!1).S(k!1,active~0!?) -> K(s!1).S(k!1,active~P!?) kcat + 1"

    @classmethod
    def teardown_class(cls):
        pass

    def test_mdef_parse(self):
        assert readers.BNGLReader.parse_mtype(self.mdef0).write_as_bngl() == self.mdef0
        md1 = readers.BNGLReader.parse_mtype(self.mdef1)
        md1.site_name_map['b0'] = 'b'
        md1.site_name_map['b1'] = 'b'
        print md1.write_as_bngl()
        assert md1.write_as_bngl() == "Mol(a,b~0~1,c~a~b,b~0~1,c~a~b)"

    def test_mol_parse(self):
        assert readers.BNGLReader.parse_molecule(self.mol0).write_as_bngl() == self.mol0
        mol1 = readers.BNGLReader.parse_molecule(self.mol1)
        mol1.write_as_bngl() == "Mol(a,b~0,b~1)"

    def test_init_parse(self):
        assert readers.BNGLReader.parse_init(self.init0).write_as_bngl() == self.mol0 + ' 100.0'
        assert readers.BNGLReader.parse_init(self.init1).write_as_bngl() == self.mol0 + ' (x+3)/k'

    def test_obs_parse(self):
        assert readers.BNGLReader.parse_obs(self.obs0).write_as_bngl() == self.obs0
        assert readers.BNGLReader.parse_obs(self.obs1).write_as_bngl() == self.obs1

    def test_params_parse(self):
        assert readers.BNGLReader.parse_param(self.param0).write_as_bngl() == self.param0
        assert readers.BNGLReader.parse_param(self.param1).write_as_bngl() == "kp km/kd/(NA*V)"

    def test_rule_parse(self):
        prule0 = readers.BNGLReader.parse_rule(self.rule0)
        assert prule0.rev is True
        assert prule0.write_as_bngl() == "A(a)+B(b) <-> A(a!1).B(b!1) kp,km"
        prule1 = readers.BNGLReader.parse_rule(self.rule1)
        assert prule1.write_as_bngl() == "A(a~r)+B(b,c!1).C(c!1) -> A(a~r!2).B(b!2,c!1).C(c!1) kp/log10(10)"
        # assert prule1.write_as_kappa() == "A(a~r),B(b,c!1),C(c!1) -> A(a~r!2),B(b!2,c!1),C(c!1) @ 'kp'/([log](10)/[log](10))"
        prule2 = readers.BNGLReader.parse_rule(self.rule2)
        assert prule2.rate.intra_binding is True
        assert prule2.write_as_bngl() == self.rule2
        # assert prule2.write_as_kappa() == "A(a~s),B(b,c!1),C(c!1) -> A(a~s!2),B(b!2,c!1),C(c!1) @ 0 {'kp'/([log](10)/[log](10))}"
        prule3 = readers.BNGLReader.parse_rule(self.rule3)
        assert prule3.rate.intra_binding is False
        assert prule3.write_as_bngl() == "K(s!1).S(k!1,active~0!?) -> K(s!1).S(k!1,active~P!?) kcat+1"
        # assert prule3.write_as_kappa() == "K(s!1),S(k!1,active~0?) -> K(s!1),S(k!1,active~P?) @ 'kcat'+1"
