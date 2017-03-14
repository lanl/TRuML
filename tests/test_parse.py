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

		cls.md0 = MoleculeDef('Molec',{'site0':['a','b'],'site1':[]})

		cls.m0 = Molecule('Test0',[cls.s0,cls.s2])
		cls.m1 = Molecule('Test1',[cls.s1])

		cls.s3 = Site('site0',s='state')
		cls.p0 = Pattern([cls.m0,Molecule('Test0',[cls.s0])])
		cls.p1 = Pattern([Molecule('Test0',[cls.s3,cls.s2]),Molecule('Test0',[cls.s3])])

		cls.i0 = InitialCondition(cls.p0,10)
		# implement functionality to print initial condition as kappa/bngl expression
		# cls.i1 = InitialCondition(cls.p0,"x+10")

		cls.par0 = Parameter('rate',1e6)

		cls.expr0 = Expression('rate_expr',['ln', '(', '10', ')', '+', 'x', '-', '356'])

		cls.rule0 = Rule([Pattern([Molecule('A',[])])],[Pattern([Molecule('B',[])])],Rate(10))
		cls.rule1 = Rule([Pattern([Molecule('A',[])])],[Pattern([Molecule('B',[])])],Rate(cls.expr0),False,'x')
		cls.rule2 = Rule([Pattern([Molecule('A',[])])],[Pattern([Molecule('B',[])])],Rate(10),True,Rate(cls.par0))

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
		assert self.md0.write_as_bngl() == r'Molec(site0~a~b,site1)'
		assert self.md0.write_as_kappa() == r'%agent: Molec(site0~a~b,site1)'

	def test_molecules(self):
		assert self.m0.write_as_bngl() == r'Test0(site0~state!3,site2!+)'
		assert self.m0.write_as_kappa() == r'Test0(site0~state!3,site2!_)'
		assert self.m1.write_as_bngl() == r'Test1(site1!?)'
		assert self.m1.write_as_kappa() == r'Test1(site1?)'

	def test_patterns(self):
		assert self.p0.write_as_bngl() == r'Test0(site0~state!3,site2!+).Test0(site0~state!3)'
		assert self.p0.write_as_kappa() == r'Test0(site0~state!3,site2!_),Test0(site0~state!3)'

	def test_init_conditions(self):
		assert self.i0.write_as_bngl() == 'Test0(site0~state!3,site2!+).Test0(site0~state!3)\t10'
		assert self.i0.write_as_kappa() == r'%init: 10 Test0(site0~state!3,site2!_),Test0(site0~state!3)'
		# assert self.i1.write_as_bngl() == r'Test0(site0~state!3,site2!+).Test0(site0~state!3)\tx+10'
		# assert self.i1.write_as_kappa() == r"%init: 'x' + 10 Test0(site0~state!3,site2!+).Test0(site0~state!3)"

	#TODO write more to check function map functionality
	def test_pars_and_exprs(self):
		assert self.par0.write_as_bngl() == r'rate 1000000.0'
		assert self.par0.write_as_kappa() == r"%var: 'rate' 1000000.0"
		assert self.expr0.write_as_bngl() == r'rate_expr=ln(10)+x-356'
		assert self.expr0.write_as_kappa() == r"%obs: 'rate_expr' [log](10)+'x'-356"
		assert self.expr0.write_as_kappa(as_obs=False) == r"%var: 'rate_expr' [log](10)+'x'-356"

	def test_rule(self):
		assert self.rule0.write_as_bngl() == r'A() -> B() 10'
		assert self.rule0.write_as_kappa() == r'A() -> B() @ 10'
		# assert self.rule1.write_as_bngl() == r'A() -> B() 10'
		# assert self.rule1.write_as_kappa() == r'A() -> B() @ 10'
		assert self.rule2.write_as_bngl() == r'A() <-> B() 10,rate'
		print self.rule2.write_as_kappa(),'\n',r"A() <-> B() 10,'rate'"
		assert self.rule2.write_as_kappa() == r"A() <-> B() @ 10,'rate'"



