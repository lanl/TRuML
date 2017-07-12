from rbconvert.objects import *
from rbconvert.readers import *

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


class TestParseKappa:
    def __init__(self):
        pass


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
        assert BNGLReader.parse_mtype(self.mdef0).write_as_bngl() == self.mdef0
        md1 = BNGLReader.parse_mtype(self.mdef1)
        md1.site_name_map['b0'] = 'b'
        md1.site_name_map['b1'] = 'b'
        print md1.write_as_bngl()
        assert md1.write_as_bngl() == "Mol(a,b~0~1,c~a~b,b~0~1,c~a~b)"

    def test_mol_parse(self):
        assert BNGLReader.parse_molecule(self.mol0).write_as_bngl() == self.mol0
        mol1 = BNGLReader.parse_molecule(self.mol1)
        mol1.write_as_bngl() == "Mol(a,b~0,b~1)"

    def test_init_parse(self):
        assert BNGLReader.parse_init(self.init0).write_as_bngl() == self.mol0 + ' 100.0'
        assert BNGLReader.parse_init(self.init1).write_as_bngl() == self.mol0 + ' (x+3)/k'

    def test_obs_parse(self):
        assert BNGLReader.parse_obs(self.obs0).write_as_bngl() == self.obs0
        assert BNGLReader.parse_obs(self.obs1).write_as_bngl() == self.obs1

    def test_params_parse(self):
        assert BNGLReader.parse_param(self.param0).write_as_bngl() == self.param0
        assert BNGLReader.parse_param(self.param1).write_as_bngl() == "kp km/kd/(NA*V)"

    def test_rule_parse(self):
        prule0 = BNGLReader.parse_rule(self.rule0)
        assert prule0.rev is True
        assert prule0.write_as_bngl() == "A(a)+B(b) <-> A(a!1).B(b!1) kp,km"
        prule1 = BNGLReader.parse_rule(self.rule1)
        assert prule1.write_as_bngl() == "A(a~r)+B(b,c!1).C(c!1) -> A(a~r!2).B(b!2,c!1).C(c!1) kp/log10(10)"
        # assert prule1.write_as_kappa() == "A(a~r),B(b,c!1),C(c!1) -> A(a~r!2),B(b!2,c!1),C(c!1) @ 'kp'/([log](10)/[log](10))"
        prule2 = BNGLReader.parse_rule(self.rule2)
        assert prule2.rate.intra_binding is True
        assert prule2.write_as_bngl() == self.rule2
        # assert prule2.write_as_kappa() == "A(a~s),B(b,c!1),C(c!1) -> A(a~s!2),B(b!2,c!1),C(c!1) @ 0 {'kp'/([log](10)/[log](10))}"
        prule3 = BNGLReader.parse_rule(self.rule3)
        assert prule3.rate.intra_binding is False
        assert prule3.write_as_bngl() == "K(s!1).S(k!1,active~0!?) -> K(s!1).S(k!1,active~P!?) kcat+1"
        # assert prule3.write_as_kappa() == "K(s!1),S(k!1,active~0?) -> K(s!1),S(k!1,active~P?) @ 'kcat'+1"
