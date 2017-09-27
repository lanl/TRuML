from nose.tools import raises
from .context import readers
from .context import objects
from .context import rbexceptions


class TestParseKappa:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        cls.mdef0 = "%agent: M(x,y,z~0~1)"
        cls.mdef1 = "%agent: M+_2-985798f(x, y, z~0~1)"

        cls.mdef2 = "%agent: A(x,a~0~1)"
        cls.mdef3 = "%agent: B(x,y,site~state~state2)"
        cls.mdef4 = "%agent: C(y,z~t~s)"

        cls.mds = [cls.mdef2, cls.mdef3, cls.mdef4]

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
        cls.rule3 = "'label with space' A(x),B(x) <-> A(x!1),B(x!1) @ %s {0}, 0.01" % cls.expr1
        cls.rule4 = "A(x),B(x) <-> A(x!1),B(x!1) @ 1, 0.1"
        cls.rule5 = " <-> A(x) @ 'rate', 'rate'"
        cls.rule6 = "B(site~state!_) -> @ [log] 3"

        cls.obs0 = "%obs: 'ste5 dimerized' |Ste5(ste5!1),Ste5(ste5!1)|"

    def test_rule_parse(self):
        mds = [readers.KappaReader.parse_mtype(x) for x in self.mds]
        rule0s = readers.KappaReader.parse_rule(self.rule0, mds)
        assert len(rule0s) == 1
        assert not rule0s[0].rev
        rule1s = readers.KappaReader.parse_rule(self.rule1, mds)
        assert len(rule1s[0].lhs) == len(rule0s[0].lhs)
        for i in range(len(rule1s[0].lhs)):
            assert rule1s[0].lhs[i].is_isomorphic(rule0s[0].lhs[i])
        assert rule1s[0].label == 'rule'
        rule2s = readers.KappaReader.parse_rule(self.rule2, mds)
        assert len(rule2s) == 4
        for r in rule2s:
            assert not r.rev
        rule3s = readers.KappaReader.parse_rule(self.rule3, mds)
        assert len(rule3s) == 3
        assert rule3s[0].label == 'label with space'
        rule4s = readers.KappaReader.parse_rule(self.rule4, mds)
        assert len(rule4s) == 1
        assert rule4s[0].rev
        rule5s = readers.KappaReader.parse_rule(self.rule5, mds)
        assert rule5s[0].lhs == []
        assert rule5s[0].rev
        rule6s = readers.KappaReader.parse_rule(self.rule6, mds)
        assert rule6s[0].rhs == []
        assert rule6s[0].delmol

    def test_cpattern_parse(self):
        pmdef2 = readers.KappaReader.parse_mtype(self.mdef2)
        pmdef3 = readers.KappaReader.parse_mtype(self.mdef3)
        pmdef4 = readers.KappaReader.parse_mtype(self.mdef4)
        pcp0 = readers.KappaReader.parse_cpatterns(self.cp0, [pmdef2, pmdef3, pmdef4])
        assert len(pcp0) == 2
        assert len(pcp0[1].molecule_list) == 2
        assert len(pcp0[0].molecule_list) == 1
        assert pcp0[1].molecule_list[1].sites[0].bond.num == 1
        assert not pcp0[0].molecule_list[0].sites
        assert pcp0[1].molecule_list[1].sites[1].state == 's'

    def test_init_parse(self):
        pmdef2 = readers.KappaReader.parse_mtype(self.mdef2)
        assert readers.KappaReader.parse_init(self.init0, [pmdef2])[0].write_as_kappa() == "%init: 10 A(x)"
        bdef = objects.MoleculeDef('B', [], {})
        cdef = objects.MoleculeDef('C', [], {})
        assert readers.KappaReader.parse_init(self.init1, [bdef, cdef])[0].write_as_kappa() == "%init: 10+'x' B()"
        assert readers.KappaReader.parse_init(self.init1, [bdef, cdef])[1].write_as_kappa() == "%init: 10+'x' C()"

    def test_eq_parse(self):
        assert readers.KappaReader.parse_alg_expr(self.expr0).asList() == ['10', '+', 'x']
        assert readers.KappaReader.parse_alg_expr(self.expr1).asList() == \
               ['[log]', '100', '/', '[max]', '10', '100', '-', '[int]', '7.342']

    def test_mdef_parse(self):
        assert readers.KappaReader.parse_mtype(self.mdef0).write_as_kappa() == "%agent: M(x,y,z~0~1)"

    def test_mol_parse(self):
        pmdef0 = readers.KappaReader.parse_mtype(self.mdef0)
        assert readers.KappaReader.parse_molecule(self.mol0, [pmdef0]).write_as_kappa() == "M(x,y,z~0)"
        pmol1 = readers.KappaReader.parse_molecule(self.mol1, [pmdef0])
        assert pmol1.write_as_kappa() == "M(x,y,z~0!_)"
        assert pmol1.sites[2].bond.wild
        pmdef1 = readers.KappaReader.parse_mtype(self.mdef1)
        pmol2 = readers.KappaReader.parse_molecule(self.mol2, [pmdef1])
        assert pmol2.write_as_kappa() == "M+_2-985798f(x?,y,z~0)"
        assert pmol2.sites[0].bond.any
        assert pmol2.name == "M+_2-985798f"

    def test_vars_parse(self):
        kr = readers.KappaReader()
        kr.lines = ["%agent: C(x, y~state~state2)", "%agent: A()", "%agent: Ste5(ste5, ste4)", "%agent: Ste4(ste5)",
                    "%var: 'a' 3", "%var: 'b' 3 + 'a'", "%var: 'c' |C(x!_,y~state?)|", "%var: 'd' |A()| + 'b'",
                    "%obs: 'membrane Ste5' |Ste5(ste4!1),Ste4(ste5!1)|", "%var: 'combo' 'membrane Ste5' / 'a'", self.obs0]
        model = kr.parse()
        assert model.molecules[0].name == 'C'
        print model.molecules[0].site_name_map
        assert model.molecules[0].site_name_map == {'x': 'x', 'y': 'y'}
        assert len(model.functions) == 2
        assert model.functions[1].name == 'combo'
        assert model.parameters[0].name == 'a'
        assert model.parameters[1].name == 'b'
        assert set(model.parameters[1].value.atom_list) == {'3', '+', 'a'}
        assert isinstance(model.parameters[1].value, objects.Expression)
        assert model.observables[0].write_as_kappa() == "%obs: 'c' |C(x!_,y~state?)|"
        assert len(model.observables) == 4
        assert model.observables[1].name == "anon_obs0"
        assert model.observables[-2].name == "membrane Ste5"
        assert len(model.observables[-1].cpatterns) == 1


class TestParseBNGL:
    def __init__(self):
        pass

    @classmethod
    def setup_class(cls):
        cls.mdef0 = "Molecule(site0,site1~0~P~PP)"
        cls.mdef1 = "Mol(a,b~0~1,c~a~b,b~0~1,c~a~b)"
        cls.mdef2 = "Mol(sa,sb,sc,sd~0~1)"
        cls.mdef3 = "Mol(a,b~0~1,sb~0~1)"
        cls.mdef4 = "A(a~r~s)"
        cls.mdef5 = "B(b,c)"
        cls.mdef6 = "C(c)"
        cls.mdef7 = "K(s)"
        cls.mdef8 = "S(k, active~0~P)"
        cls.mds = [cls.mdef4, cls.mdef5, cls.mdef6, cls.mdef7, cls.mdef8]
        cls.mol0 = "Mol(sa,sb!+,sc!3,sd~0!?)"
        cls.mol1 = "Mol(a,b~0,b~1)"
        cls.mol2 = "Mol(a,b~?!+)"
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
        cls.rule4 = "A() <-> 0 rate, rate"
        cls.rule5 = "0 -> B(x) 4"
        cls.rule6 = "bdeg: B(x!+) -> 0 kdeg DeleteMolecules"

    @classmethod
    def teardown_class(cls):
        pass

    def test_mdef_parse(self):
        assert readers.BNGLReader.parse_mtype(self.mdef0).write_as_bngl() == self.mdef0
        md1 = readers.BNGLReader.parse_mtype(self.mdef1)
        md1.site_name_map['b0'] = 'b'
        md1.site_name_map['b1'] = 'b'
        assert md1.write_as_bngl() == "Mol(a,b~0~1,c~a~b,b~0~1,c~a~b)"

    def test_mol_parse(self):
        pmdef2 = readers.BNGLReader.parse_mtype(self.mdef2)
        assert readers.BNGLReader.parse_molecule(self.mol0, [pmdef2]).write_as_bngl() == self.mol0
        pmdef3 = readers.BNGLReader.parse_mtype(self.mdef3)
        mol1 = readers.BNGLReader.parse_molecule(self.mol1, [pmdef3])
        mol1.write_as_bngl() == "Mol(a,b~0,b~1)"
        mol2 = readers.BNGLReader.parse_molecule(self.mol2, [pmdef3])
        mol2.write_as_bngl() == "Mol(a,b!+)"

    def test_init_parse(self):
        pmd2 = readers.BNGLReader.parse_mtype(self.mdef2)
        assert readers.BNGLReader.parse_init(self.init0, [pmd2]).write_as_bngl() == self.mol0 + ' 100.0'
        assert readers.BNGLReader.parse_init(self.init1, [pmd2]).write_as_bngl(
            {'x': 'x', 'k': 'k'}) == self.mol0 + ' (x+3)/k'

    def test_obs_parse(self):
        pmdef2 = readers.BNGLReader.parse_mtype(self.mdef2)
        pmdef3 = readers.BNGLReader.parse_mtype(self.mdef3)
        assert readers.BNGLReader.parse_obs(self.obs0, [pmdef2]).write_as_bngl({"Mol0": "Mol0"}) == self.obs0
        assert readers.BNGLReader.parse_obs(self.obs1, [pmdef3]).write_as_bngl({"Mol1": "Mol1"}) == self.obs1

    def test_params_parse(self):
        assert readers.BNGLReader.parse_param(self.param0).write_as_bngl({"kcat": "kcat"}) == self.param0
        namespace = {"kp": "kp", "km": "km", "kd": "kd", "NA": "NA", "V": "V"}
        assert readers.BNGLReader.parse_param(self.param1).write_as_bngl(namespace) == "kp km/kd/(NA*V)"

    def test_rule_parse(self):
        mds = [readers.BNGLReader.parse_mtype(x) for x in self.mds]
        prule0 = readers.BNGLReader.parse_rule(self.rule0, mds)
        assert prule0.rev is True
        assert prule0.write_as_bngl({"kp": "kp", "km": "km"}) == "A(a)+B(b) <-> A(a!1).B(b!1) kp,km"
        prule1 = readers.BNGLReader.parse_rule(self.rule1, mds)
        assert prule1.write_as_bngl({"kp": "kp"}) == "A(a~r)+B(b,c!1).C(c!1) -> A(a~r!2).B(b!2,c!1).C(c!1) kp/log10(10)"
        assert prule1.write_as_kappa() == "A(a~r),B(b,c!1),C(c!1) -> A(a~r!2),B(b!2,c!1),C(c!1) @ 'kp'/([log](10)/[log](10))"
        prule2 = readers.BNGLReader.parse_rule(self.rule2, mds)
        assert prule2.rate.intra_binding is True
        assert prule2.write_as_bngl({"kp": "kp"}) == self.rule2
        assert prule2.write_as_kappa() == "A(a~s),B(b,c!1),C(c!1) -> A(a~s!2),B(b!2,c!1),C(c!1) @ 0 {'kp'/([log](10)/[log](10))}"
        prule3 = readers.BNGLReader.parse_rule(self.rule3, mds)
        assert prule3.rate.intra_binding is False
        assert prule3.write_as_bngl({"kcat": "kcat"}) == "K(s!1).S(k!1,active~0!?) -> K(s!1).S(k!1,active~P!?) kcat+1"
        assert prule3.write_as_kappa() == "K(s!1),S(k!1,active~0?) -> K(s!1),S(k!1,active~P?) @ 'kcat'+1"
        prule4 = readers.BNGLReader.parse_rule(self.rule4, mds)
        assert len(prule4.rhs) == 0
        assert prule4.rate.write_as_bngl({"rate": "rate2"}) == 'rate2'
        assert prule4.delmol
        prule5 = readers.BNGLReader.parse_rule(self.rule5, mds)
        assert len(prule5.lhs) == 0
        assert len(prule5.rhs) == 1
        prule6 = readers.BNGLReader.parse_rule(self.rule6, mds)
        assert prule6.label == 'bdeg'
        assert prule6.rate.rate == 'kdeg'
        assert prule6.delmol

    @raises(rbexceptions.NotCompatibleException)
    def test_invalid_rule_rate(self):
        lines = ['begin molecule types', 'A(x)', 'end molecule types',
                 'begin observables', 'Molecules A A().A()', 'end observables',
                 'begin functions', 'f() A / 2', 'end functions',
                 'begin reaction rules', 'A(x)+A(x) -> A(x!1).A(x!1) f()', 'end reaction rules']

        br = readers.BNGLReader()
        br.lines = lines
        m = br.parse()
        m.write_as_kappa('test.ka', True)
