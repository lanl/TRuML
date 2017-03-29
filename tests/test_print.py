import os
import sys
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))

from rbconvert.parse_bngl import *

from nose.tools import raises

class TestPrint:

	@classmethod
	def setup_class(cls):
		cls.num_bond = Bond(3)
		cls.wild_bond = Bond(-1,w=True)
		cls.any_bond = Bond(-3,a=True)

		cls.s0 = Site('site0',s='state',b=cls.num_bond)
		cls.s1 = Site('site1',b=cls.any_bond)
		cls.s2 = Site('site2',b=cls.wild_bond)
		cls.s3 = Site('site0',s='state')
		
		cls.md0 = MoleculeDef('Test0',{'site0':['state','state2'],'site1':[]},{'site0':'site0','site1':'site1'})
		cls.md1 = MoleculeDef('Test1',{'site1':[]},{'site1':'site1'})
		cls.md2 = MoleculeDef('Test2',{'a':[],'a':[],'a':[],'a':[],'b':[],'c':[],'c':[]},{'a0':'a','a1':'a','a2':'a','a3':'a','b':'b','c0':'c','c1':'c'},hss=True)

		cls.m0 = Molecule('Test0',[cls.s0,cls.s2])
		cls.m1 = Molecule('Test1',[cls.s1])
		cls.m2 = Molecule('Test2',[Site('a')])
		cls.m3 = Molecule('Test2',[Site('a'),Site('a')])
		cls.m4 = Molecule('Test2',[Site('a'),Site('a',b=Bond(1)),Site('c')])
		cls.m5 = Molecule('Test2',[Site('a'),Site('a',b=Bond(1)),Site('a'),Site('b')])
		cls.m6 = Molecule('Test2',[Site('a'),Site('a',b=Bond(1)),Site('a',Bond(2)),Site('b')])

		cls.p0 = CPattern([cls.m0,Molecule('Test0',[cls.s0])])
		cls.p1 = CPattern([Molecule('Test0',[cls.s3,cls.s2]),Molecule('Test0',[cls.s3])])
		cls.p2 = CPattern([Molecule('A',[])])
		cls.p3 = CPattern([Molecule('B',[])])

		cls.i0 = InitialCondition(cls.p0,10)
		# implement functionality to print initial condition as kappa/bngl expression
		# cls.i1 = InitialCondition(cls.p0,"x+10")

		cls.par0 = Parameter('rate',1e6)

		expr0 = Expression(['ln', '(', '10', ')', '+', 'x', '-', '356'])
		cls.func0 = Function('rate_expr',expr0)

		cls.rate0 = Rate(3)
		cls.rate1 = Rate(expr0)
		cls.rate2 = Rate('rate')
		cls.rate3 = Rate(5,intra=True)

		cls.rule0 = Rule([cls.p2],[cls.p3],cls.rate0)
		cls.rule1 = Rule([cls.p2],[cls.p3],cls.rate1,False,cls.rate2)
		cls.rule2 = Rule([cls.p2],[cls.p3],cls.rate0,True,cls.rate2)

		cls.obs0 = Observable("Obs0",[cls.p3],'m')
		cls.obs1 = Observable("Obs1",[cls.p2,cls.p3],'s')

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
		assert self.md0.write_as_bngl() == r'Test0(site0~state~state2,site1)'
		assert self.md0.write_as_kappa() == r'%agent: Test0(site0~state~state2,site1)'

	def test_molecules(self):
		assert self.m0.write_as_bngl() == r'Test0(site0~state!3,site2!+)'
		assert self.m0.write_as_kappa(self.md0)[0] == r'Test0(site0~state!3,site2!_)'
		assert self.m1.write_as_bngl() == r'Test1(site1!?)'
		assert self.m1.write_as_kappa(self.md1)[0] == r'Test1(site1?)'

	def test_site_renaming(self):
		km2 = self.m2.write_as_kappa(self.md2)
		km3 = self.m3.write_as_kappa(self.md2)
		km4 = self.m4.write_as_kappa(self.md2)
		km5 = self.m5.write_as_kappa(self.md2)
		km6 = self.m6.write_as_kappa(self.md2)
		assert len(km2) == 4 # 4 choose 1
		assert len(km3) == 6 # 4 choose 2
		assert len(km4) == 24 # (4 choose 2) * (2 choose 1) * (2 choose 1)
		assert len(km5) == 12 # (4 choose 2) * (2 choose 1)
		assert len(km6) == 24 # (4 choose 1) * (3 choose 1) * (2 choose 1)

	def test_patterns(self):
		assert self.p0.write_as_bngl() == r'Test0(site0~state!3,site2!+).Test0(site0~state!3)'
		# assert self.p0.write_as_kappa([self.md0]) == r'Test0(site0~state!3,site2!_),Test0(site0~state!3)'

	def test_init_conditions(self):
		assert self.i0.write_as_bngl() == 'Test0(site0~state!3,site2!+).Test0(site0~state!3) 10'
		# assert self.i0.write_as_kappa([self.md0]) == r'%init: 10 Test0(site0~state!3,site2!_),Test0(site0~state!3)'
		# assert self.i1.write_as_bngl() == r'Test0(site0~state!3,site2!+).Test0(site0~state!3)\tx+10'
		# assert self.i1.write_as_kappa() == r"%init: 'x' + 10 Test0(site0~state!3,site2!+).Test0(site0~state!3)"

	#TODO write more to check function map functionality
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
		# assert self.rule0.write_as_kappa() == r'A() -> B() @ 3'
		assert self.rule1.write_as_bngl() == r'A() -> B() ln(10)+x-356'
		# assert self.rule1.write_as_kappa() == r"A() -> B() @ [log](10)+'x'-356"
		assert self.rule2.write_as_bngl() == r'A() <-> B() 3,rate'
		# assert self.rule2.write_as_kappa() == r"A() <-> B() @ 3,'rate'"

	def test_obs(self):
		assert self.obs0.write_as_bngl() == r'Molecules Obs0 B()'
		# assert self.obs0.write_as_kappa() == r"%obs: 'Obs0' |B()|"
		assert self.obs1.write_as_bngl() == r'Species Obs1 A() B()'
		# assert self.obs1.write_as_kappa() == r"%obs: 'Obs1' |A()|+|B()|"

	@raises(Exception)
	def test_invalid_obs(self):
		Observable("InvalidType",[cls.p2],'f')

