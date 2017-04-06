import re
from pyparsing import Literal,CaselessLiteral,Word,Combine,Group,Optional,\
    ZeroOrMore,Forward,nums,alphas
from deepdiff import DeepDiff
from itertools import product, combinations
from copy import deepcopy

class NotAMoleculeException(Exception):
	def __init__(self,s):
		print "%s is not a molecule"%s

class NotCompatibleException(Exception):
	def __init__(self,s):
		print s

# full definition
class MoleculeDef:
	def __init__(self,n,st,snm,hss=False):
		self.name = n
		self.sites = st
		# map of Kappa name to BNGL name
		self.site_name_map = snm
		self.has_site_symmetry = hss
		self._invert_site_name_map()

	# map of BNGL name to list of Kappa names (includes 1-1 mappings)
	def _invert_site_name_map(self):
		self.inv_site_name_map = {}
		for k in self.site_name_map.keys():
			v = self.site_name_map[k] 
			if v not in self.inv_site_name_map.keys():
				self.inv_site_name_map[v] = [k]
			else:
				self.inv_site_name_map[v].append(k)

	def _all_site_states(self,is_bngl=True):
		ss = []
		if not is_bngl:
			k_track = {s:list(reversed(sorted(self.inv_site_name_map[s]))) for s in self.inv_site_name_map.keys()}
		# ordered easy testing
		for k,v in self.sites:
			site_name = None
			if is_bngl:
				site_name = k 
			else:
				site_name = k_track[k].pop()
			if not v:
				ss.append(site_name)
			else:
				ss.append('%s~%s'%(site_name,'~'.join(v)))
		return ','.join(ss)

	def write_as_bngl(self):
		return "%s(%s)"%(self.name,self._all_site_states())

	def write_as_kappa(self):
		return "%%agent: %s(%s)"%(self.name,self._all_site_states(False))

	def __repr__(self):
		sites_string = ','.join(['%s->%s'%(k,v) for k,v in self.sites])
		return "MoleculeDef(name: %s, sites: %s)"%(self.name,sites_string)

class Molecule:
	def __init__(self,name,sites):
		self.name = name
		self.sites = sites # list of Sites
		self.site_dict = {s.name:s for s in self.sites}
	
	def write_as_bngl(self):
		return self.name + '(' + ','.join([s.write_as_bngl() for s in self.sites]) + ')'

	# returns list of molecule strings
	# check to see how many of possible symm sites are present in pattern
	def write_as_kappa(self,mdef):

		def kappa_string(ss):
			return self.name + '(' + ','.join(sorted(ss)) + ')' # sorted for consistent testing

		un_site_names = set([s.name for s in self.sites])
		un_configs_per_site = {s:{} for s in un_site_names}
		for s in self.sites:
			if s not in un_configs_per_site[s.name]:
				un_configs_per_site[s.name][s] = 0
			un_configs_per_site[s.name][s] += 1

		def rename_sites(names,site):
			return tuple([Site(name,s=site.state,b=site.bond) for name in names])

		k_configs = {}
		for sn in un_configs_per_site.keys():
			k_configs[sn] = []
			k_sn_names = set(mdef.inv_site_name_map[sn])
			cur_combs = []
			for s,n in un_configs_per_site[sn].iteritems():
				if len(cur_combs) == 0:
					cur_combs = [rename_sites(names,s) for names in combinations(k_sn_names,n)]
				else:
					tmp_combs = []
					for cc in cur_combs:
						rem_names = k_sn_names - set(map(lambda l: l.name, cc))
						new_combs = [rename_sites(names,s) for names in combinations(rem_names,n)]
						for nc in new_combs:
							tmp_combs.append(cc+nc)
					cur_combs = tmp_combs
			k_configs[sn] = cur_combs

		k_prod = list(product(*k_configs.values()))
		return [kappa_string([e.write_as_kappa() for t in tt for e in t]) for tt in k_prod]



	def __repr__(self):
		return 'Molecule(name: %s, sites: %s)'%(self.name,', '.join([str(x) for x in self.sites]))

class Site:
	def __init__(self,n,s=None,b=None):
		self.name = n
		self.state = s
		self.bond = b

	def _add_state_bond(self,s,kappa=False):
		if self.state is not None:
			s += '~%s'%self.state
		if self.bond is not None and kappa:
			s += self.bond.write_as_kappa()
		elif self.bond is not None:
			s += self.bond.write_as_bngl()
		return s

	def write_as_bngl(self):
		return self._add_state_bond(self.name)

	# for testing only (check for site symmetries occurs in Molecule class)
	def write_as_kappa(self):
		return self._add_state_bond(self.name,True)

	def __eq__(self,other):
		if isinstance(other,self.__class__):
			return (self.name == other.name and self.state == other.state and self.bond == other.bond)
		return False

	def __ne__(self,other):
		return not self == other

	def __hash__(self):
		return hash((self.name,self.state,self.bond))

	def __repr__(self):
		return 'Site(name: %s, state: %s, bond: %s)'%(self.name,self.state,self.bond)

# TODO make sure to include ! operator
class Bond:
	def __init__(self,n,w=False,a=False):
		self.wild = w
		self.num = int(n) # negative numbers indicate absence of specific bond, positive numbers will override w and a
		self.any = a
		if self.num < 0:
			self.num = -1
			assert(self.wild or self.any)
		assert(not (self.wild and self.any))

	def write_as_bngl(self):
		s = ''
		if self.num >= 0:
			s = '!%s'%self.num
		if self.wild:
			s = '!+'
		elif self.any:
			s = '!?'
		return s

	def write_as_kappa(self):
		s = ''
		if self.wild:
			s = '!_'
		elif self.any:
			s = '?'
		else:
			s = '!%s'%self.num
		return s

	def __eq__(self,other):
		if isinstance(other,self.__class__):
			num_check = False
			if (self.num < 0 and other.num < 0) or (self.num == other.num):
				num_check = True
			return (self.wild == other.wild and self.any == other.any and num_check)
		return False

	def __ne__(self,other):
		return not self == other

	def __hash__(self):
		return hash((self.num,self.wild,self.any))

	def __repr__(self):
		b_string = str(self.num)
		if self.wild:
			b_string = 'wild'
		elif self.any:
			b_string = 'any'
		return b_string

# CPattern is defined as a pattern (list of molecules) that are attached
class CPattern:
	def __init__(self,ml):
		self.molecule_list = ml

	def write_as_bngl(self):
		return '.'.join([m.write_as_bngl() for m in self.molecule_list])

	def write_as_kappa(self,mdefs): # mdefs is all the moleculedefs corresponding to molecules in the molecule list
		k_str_mol_list = []
		for m in self.molecule_list:
			for md in mdefs:
				if m.name == md.name:
					k_str_mol_list.append(m.write_as_kappa(md))
		k_patterns = product(*k_str_mol_list)
		# remove doubles, preserve molecule order
		seen = set()
		un_k_patterns = []
		for pat in k_patterns:
			s_pat = tuple(sorted(pat))
			if s_pat not in seen:
				seen.add(s_pat)
				un_k_patterns.append(','.join(pat))
		return un_k_patterns

	def __repr__(self):
		return '\n'.join([str(x) for x in self.molecule_list])

class InitialCondition:
	def __init__(self,s,a,ain=True):
		self.species = s
		self.amount = a # can be number or expression
		self.amount_is_number = ain

	def write_as_bngl(self):
		amount = self.amount if self.amount_is_number else self.amount.write_as_bngl()
		return '%s %s'%(self.species.write_as_bngl(),amount)

	# if there are symmetries, the initial condition amount is divided evenly among species
	def write_as_kappa(self,mdefs):
		symm_strings = self.species.write_as_kappa(mdefs)
		num_symm = len(symm_strings)
		if num_symm == 1:
			amount = self.amount if self.amount_is_number else self.amount.write_as_kappa()
		else:
			amount = self.amount/float(num_symm) if self.amount_is_number else self.amount.write_as_kappa()
		res = []
		for s in symm_strings:
			res.append('%%init: %s %s'%(amount,s))
		return res

	def __repr__(self):
		return "Init(species: %s, quantity: %s)"%(self.species,self.amount)

class Parameter:
	def __init__(self,n,v):
		self.name = n # string
		self.value = v # number

	def write_as_bngl(self):
		return '%s %s'%(self.name,self.value)

	def write_as_kappa(self):
		return '%%var: \'%s\' %s'%(self.name, self.value)

	def __repr__(self):
		return "Parameter(name: %s, value: %s)"%(self.name,self.value)

# special parsing required for 'if', 'log' functions
# can implement conversion of certain types of values (e.g. log10(x) to log(x)/log(10))
class Expression:
	def __init__(self,atom_list):
		self.atom_list = atom_list # list from parse_math_expr listing (in order) operators, values, variables

	def write_as_bngl(self):
		return ''.join(self.atom_list)

	def write_as_kappa(self):
		expr = ''

		i = 0
		while (i < len(self.atom_list)):
			a = self.atom_list[i]
			if a in bngl_to_kappa_func_map.keys():
				trig_func_match = re.compile('sinh|cosh|tanh|asinh|acosh|atanh')
				if re.match('log',a) or re.match(trig_func_match,a):
					expr += bngl_to_kappa_func_map[a](self.atom_list[i+2])
					i += 4
				else:
					expr += bngl_to_kappa_func_map[a]
			elif re.match('[A-Za-z]',a):
				expr += '\'%s\''%a
			else:
				expr += a
			i += 1

		return expr

	def __repr__(self):
		return "Expression(expr: %s)"%self.write_as_bngl()

class Function:
	def __init__(self,name,expr):
		self.name = name
		self.expr = expr

	def write_as_bngl(self):
		return '%s()=%s'%(self.name,self.expr.write_as_bngl())

	def write_as_kappa(self,as_obs=True):
		dec = 'obs' if as_obs else 'var'
		return "%%%s: '%s' %s"%(dec,self.name,self.expr.write_as_kappa())

	def __repr__(self):
		return "Function(name: %s, expr: %s"%(self.name, self.expr)

class Rate:
	def __init__(self,r,intra=False): # can be string (including a single number) or Expression
		self.rate = r
		self.intra_binding = intra

	def write_as_bngl(self):
		try:
			return self.rate.write_as_bngl()
		except AttributeError:
			return str(self.rate)

	def write_as_kappa(self):
		try:
			rate_string = self.rate.write_as_kappa()
		except AttributeError:
			rate_string =  str(self.rate) if _is_number(self.rate) else "'%s'"%self.rate
		return rate_string if not self.intra_binding else '0 {%s}'%rate_string

	def __repr__(self):
		return "Rate: %s"%self.rate

#TODO implement check for rate as raw number before writing
class Rule:
	# lhs, rhs are lists of CPatterns, rate/rev_rate are Rates, rev is bool (true for reversible rules), 
	# amb_mol is boolean denoting a rule that (in Kappa) has ambiguous molecularity
	def __init__(self,lhs,rhs,rate,rev=False,rev_rate=None):
		self.lhs = lhs
		self.rhs = rhs
		self.rate = rate
		self.rev = rev
		self.arrow = '->' if not rev else '<->'
		self.rev_rate = None if not rev else rev_rate # rev overrides rev_rate

	def write_as_bngl(self):
		lhs_string = '+'.join([p.write_as_bngl() for p in self.lhs])
		rhs_string = '+'.join([p.write_as_bngl() for p in self.rhs])
		if self.rev:
			rate_string = self.rate.write_as_bngl() + ',' + self.rev_rate.write_as_bngl()
		else:
			rate_string = self.rate.write_as_bngl()
		return '%s %s %s %s'%(lhs_string,self.arrow,rhs_string,rate_string)

	def write_as_kappa(self,mdefs):

		# possible actions  (root object is list (rule.lhs and rule.rhs))
		#  - iterable_item_added/removed (binding, unbinding)
		#  - type_changes (binding, unbinding)
		#  - value_changes (state change)

		def kappa_rule_string(lhss,ars,rhss,rs):
			return '%s %s %s @ %s'%(lhss,ars,rhss,rs)

		lhs_strings = product(*[p.write_as_kappa(mdefs) for p in self.lhs])


		# lhs_string = ','.join([p.write_as_kappa() for p in self.lhs])
		# rhs_string = ','.join([p.write_as_kappa() for p in self.rhs])
		# if self.rev:
		# 	rate_string = self.rate.write_as_kappa() + ',' + self.rev_rate.write_as_kappa()
		# else:
		# 	rate_string = self.rate.write_as_kappa()
		# return '%s %s %s @ %s'%(lhs_string,self.arrow,rhs_string,rate_string)

	def __repr__(self):
		if not self.rev:
			return "Rule(lhs: %s, rhs: %s, rate: %s)"%(self.lhs,self.rhs,self.rate)
		else:
			return "Rule(lhs: %s, rhs: %s, rate: %s, rev_rate: %s)"%(self.lhs,self.rhs,self.rate,self.rev_rate)

class Observable:
	def __init__(self,n,ps,t='m'):
		self.name = n
		if re.match('[sS]$',t):
			self.type = 'Species'
		elif re.match('[mM]$',t):
			self.type = 'Molecules'
		else:
			raise Exception("not a valid observable type: %s"%t)
		self.cpatterns = ps # a list of CPatterns

	def write_as_bngl(self):
		return "%s %s %s"%(self.type,self.name,' '.join([p.write_as_bngl() for p in self.cpatterns]))

	def write_as_kappa(self,mdefs):
		if self.type == 'Species':
			print "Kappa does not have a Species-like observable; printing '%s' as Molecules-like observable"%self.name

		obs_strs = []
		for p in self.cpatterns:
			# sorted for determinism (testing)
			kos = '+'.join(sorted(['|%s|'%x for x in p.write_as_kappa(mdefs)]))
			obs_strs.append(kos)

		obs = '+'.join(obs_strs)
		return '%%obs: \'%s\' %s'%(self.name,obs)

	def __repr__(self):
		return "Obs(name: %s, pattern: %s)"%(self.name, ' '.join(self.cpatterns))

# TODO check types for add_* functions
class Model:
	def __init__(self):
		self.molecules = []
		self.initial_cond = []
		self.observables = []
		self.functions = []
		self.rules = []
		self.parameters = []

	def write_as_bngl(self,file_name):
		s = 'begin model\n\nbegin parameters\n\n'
		for p in self.parameters:
			s += '\t%s\n'%p.write_as_bngl()
		s += '\nend parameters\n\n'
		s += 'begin molecule types\n\n'
		for m in self.molecules:
			s += '\t%s\n'%m.write_as_bngl()
		s += '\nend molecule types\n\n'
		s += 'begin initial conditions\n\n'
		for i in self.initial_cond:
			s += '\t%s\n'%i.write_as_bngl()
		s += '\nend initial conditions\n\n'
		s += 'begin observables\n\n'
		for o in self.observables:
			s += '\t%s\n'%o.write_as_bngl()
		s += '\nend observables\n\n'
		s += 'begin functions\n\n'
		for f in self.functions:
			s += '\t%s\n'%f.write_as_bngl()
		s += '\nend functions\n\n'
		s += 'begin reaction rules'
		for r in self.rules:
			s += '\t%s\n'%r.write_as_bngl()
		s += 'end reaction rules\n\n'
		s += 'end model\n'

		f = open('%s'%file_name)
		f.write(s)
		f.close()

	# check for rules with molecular ambiguity
	def write_as_kappa(self,func_as_obs=False):
		s = ''
		for m in self.molecules:
			s += '%s\n'%m.write_as_kappa()
		s += '\n'
		for p in self.parameters:
			s += '%s\n'%p.write_as_kappa()
		s += '\n'
		for o in self.observables:
			s += '%s\n'%o.write_as_kappa()
		s += '\n'
		for f in self.functions:
			s += '%s\n'%f.write_as_kappa(func_as_obs) # defaults to printing all "functions" as observables
		s += '\n'
		for i in self.initial_cond:
			s += '%s\n'%i.write_as_kappa()
		s += '\n'
		for r in self.rules:
			s += '%s\n'%r.write_as_kappa()
		s += '\n'

	def add_molecule(self,mol):
		self.molecules.append(mol)

	def add_init(self,init):
		self.initial_cond.append(init)

	def add_obs(self,obs):
		self.observables.append(obs)

	def add_func(self,func):
		self.functions.append(func)

	def add_rule(self,rule):
		self.rules.append(rule)

	def add_parameter(self,param):
		self.parameters.append(param)

class Reader:
	def __init__(self,file_name):
		f = open(file_name)
		d = f.readlines()
		f.close()
		self.lines = d
		self.file_name = file_name

# ignores perturbation and action commands
class KappaReader(Reader):
	pass

# ignores action commands
class BNGLReader(Reader):
	def __init__(self,file_name):
		super().__init__(file_name)
		self.is_def_block = False
		self.is_init_block = False
		self.is_param_block = False
		self.is_rule_block = False
		self.is_obs_block = False
		self.is_func_block = False

	# assumes that the file has the molecule types block before the rules block
	def parse(self):
		# to accommodate line-continuation characters
		cur_line = ''
		model = Model(self.file_name)
		for i,l in enumerate(self.lines):
			if re.match('begin parameters',l):
				self.is_param_block = True
				continue
			elif re.match('end parameters'):
				self.is_param_block = False
				continue
			elif re.match('begin molecule types',l):
				self.is_def_block = True
				continue
			elif re.match('end molecule types',l):
				self.is_def_block = False
				continue
			elif re.match('begin seed species',l):
				self.is_init_block = True
				continue
			elif re.match('end seed species',l):
				self.is_init_block = False
				continue
			elif re.match('begin observables',l):
				self.is_obs_block = True
				continue
			elif re.match('end observables',l):
				self.is_obs_block = False
				continue
			elif re.match('begin functions',l):
				self.is_func_block = True
				continue
			elif re.match('end functions',l):
				self.is_func_block = False
				continue
			elif re.match('begin reaction rules',l):
				self.is_rule_block = True
				continue
			elif re.match('end reaction rules',l):
				self.is_rule_block = False
				continue

			# determines presence of line continuation, file cannot have backslashes in other contexts
			if re.search("\\\\\s*$",l):
				# saves current line, stripping trailing and leading whitespace, continues to subsequent line
				cur_line += re.sub('\\\\','',l.strip())
				continue
			else:
				cur_line += l.strip()
				if self.is_param_block:
					model.add_parameter(self.parse_param(cur_line))
				elif self.is_def_block:
					model.add_molecule(self.parse_mtype(cur_line))
				elif self.is_init_block:
					model.add_init(self.parse_init(cur_line))
				elif self.is_obs_block:
					model.add_obs(self.parse_obs(cur_line))
				elif self.is_func_block:
					model.add_func(self.parse_func(cur_line))
				elif self.is_rule_block:
					model.add_rule(self.parse_rule(cur_line))
				else:
					continue
				cur_line = ''

		return model

	@staticmethod
	def parse_bond(b):
		if re.match('\+',b):
			return Bond(-1,w=True)
		elif re.match('\?',b):
			return Bond(-1,a=True)
		elif b.isdigit():
			return Bond(b)
		else:
			raise ValueError("Illegal bond: %s"%b)

	@staticmethod
	def parse_mtype(line):
		psplit = re.split('\(',line.strip())
		name = psplit[0]
		
		site_name_map = {} # maps distinct site names for tracking conversion to kappa if needed from equivalent site names in BNGL

		sites = re.split(',',psplit[1].strip(')'))
		site_defs = []
		site_name_counter = {}
		has_site_symmetry = False
		for s in sites:
			site_split = re.split('~',s)
			site_name = site_split[0]
			site_defs.append((site_name,[] if len(site_split) == 1 else site_split[1:]))
			if site_name in site_name_counter.keys():
				site_name_counter[site_name] += 1
				if not has_site_symmetry:
					has_site_symmetry = True
			else:
				site_name_counter[site_name] = 1

		for sn in site_name_counter.keys():
			if site_name_counter[sn] == 1:
				site_name_counter.pop(sn)
				site_name_map[sn] = sn

		for sn in site_name_counter.keys():
			while site_name_counter[sn] > 0:
				site_name_map[sn + str(site_name_counter[sn]-1)] = sn
				site_name_counter[sn] -= 1

		return MoleculeDef(name, site_defs, site_name_map, has_site_symmetry)

	@classmethod
	def parse_molecule(cls,line):
		sline = line.strip()
		msplit = re.split('\(',sline)
		mname = msplit[0]
		if not re.match('[A-Za-z]\w*\(.*\)\s*$',sline):
			raise NotAMoleculeException(sline)
		sites = re.split(',',msplit[1].strip(')'))
		if not sites[0]:
			return Molecule(mname,[])
		site_list = []
		for s in sites:
			if '~' in s:
				tsplit = re.split('~',s)
				name = tsplit[0]
				if '!' in s:
					bsplit = re.split('!',tsplit[1])
					bond = cls.parse_bond(bsplit[1])
					site_list.append(Site(name,s=bsplit[0],b=bond))
				else:
					site_list.append(Site(name,s=tsplit[1]))
			else:
				if '!' in s:
					bsplit = re.split('!',s)
					name = bsplit[0]
					bond = cls.parse_bond(bsplit[1])
					site_list.append(Site(name,b=bond))
				else:
					site_list.append(Site(s))
		return Molecule(mname,site_list)

	# TODO implement parsing for expression (need to identify variables for conversion to kappa syntax)
	@classmethod
	def parse_init(cls,line):
		isplit = re.split('\s+',line.strip())
		spec = cls.parse_cpattern(isplit[0])
		amount = ' '.join(isplit[1:])
		amount_is_number = _is_number(amount)
		p_amount = float(amount) if amount_is_number else Expression(cls.parse_math_expr(amount))
		return InitialCondition(spec,p_amount,amount_is_number)

	@classmethod
	def parse_cpattern(cls,line):
		ssplit = re.split('(?<=\))\.',line.strip())
		m_list = []
		for s in ssplit:
			m_list.append(cls.parse_molecule(s))
		return CPattern(m_list)

	@classmethod
	def parse_obs(cls,line):
		osplit = re.split('\s+',line.strip())
		otype = osplit[0][0]
		oname = osplit[1]
		oCPattern = [cls.parse_cpattern(p) for p in osplit[2:]]
		return Observable(oname,oCPattern,otype)

	@staticmethod
	def parse_param(line):
		sline = line.strip()
		s_char = ''
		for x in sline:
			if re.match('\s',x) or re.match('=',x):
				s_char = x
				break
		psplit = re.split(s_char,sline)
		pname = psplit[0]
		pexpr = s_char.join(psplit[1:])
		return Parameter(pname,pexpr)

	# assumes that pattern mapping is left to right and that there is 
	# only 1 component on either side of the rule (doesn't make sense to 
	# have components that aren't operated on).  The change will be from
	# a Site with bond = None to a Site with a Bond object containing a 
	# link to another Molecule in the same component
	@staticmethod
	def _has_intramolecular_binding(lhs_cp,rhs_cp):
		d = DeepDiff(lhs_cp,rhs_cp)
		try:
			changed = d.get('type_changes').keys()
		except AttributeError:
			return False
		num_changed_bonds = 0
		for c in changed:
			if re.search('sites\[.\]\.bond$',c):
				num_changed_bonds += 1
		return num_changed_bonds == 2

	# TODO parse rule label
	@classmethod
	def parse_rule(cls,line):
		sline = line.strip()
		rhs = ''
		lhs = ''
		is_reversible = True if re.search('<->',sline) else False
		parts = re.split('->',sline)
		lhs_cpatterns = [cls.parse_cpattern(x) for x in re.split('(?<!!)\+',parts[0].rstrip('<'))]
		rem = [x.strip() for x in re.split('(?<!!)\+',parts[1].strip())]
		# if the rule is an unbinding rule or has a '+' in its rate expression
		if len(rem) > 1:
			one_past_final_mol_index = 0
			for i,t in enumerate(rem):
				try:
					cls.parse_cpattern(t)
				except NotAMoleculeException:
					one_past_final_mol_index = i
					break
			last_split = re.split('\s+',rem[one_past_final_mol_index])
 			mol,first_rate_part = last_split[0], ' '.join(last_split[1:])
			rhs_cpatterns = [cls.parse_cpattern(x) for x in (rem[:one_past_final_mol_index] + [mol])]
			rate_string = first_rate_part + '+' + '+'.join(rem[one_past_final_mol_index+1:])
			if is_reversible:
				rate0,rate1 = re.split(',',rate_string)
				return Rule(lhs_cpatterns,rhs_cpatterns,cls.parse_rate(rate0),is_reversible,cls.parse_rate(rate1))
			else:
				return Rule(lhs_cpatterns,rhs_cpatterns,cls.parse_rate(rate_string),is_reversible)
		else:
			rem_parts = re.split('(?<!!)\s+',parts[1].strip())
			rhs_cpatterns = [cls.parse_cpattern(rem_parts[0])]
			is_intra_l_to_r = False
			if len(lhs_cpatterns) == 1 and len(rhs_cpatterns) == 1:
				is_intra_l_to_r = cls._has_intramolecular_binding(lhs_cpatterns[0],rhs_cpatterns[0])
			rate_string = ' '.join(rem_parts[1:])
			if is_reversible:
				is_intra_r_to_l = cls._has_intramolecular_binding(rhs_cpatterns[0],lhs_cpatterns[0])
				rate0_string,rate1_string = re.split(',',rate_string)
				rate0 = cls.parse_rate(rate0_string,is_intra_l_to_r)
				rate1 = cls.parse_rate(rate1_string,is_intra_r_to_l)
				return Rule(lhs_cpatterns,rhs_cpatterns,rate0,is_reversible,rate1)
			else:
				rate0 = cls.parse_rate(rate_string,is_intra_l_to_r)
				return Rule(lhs_cpatterns,rhs_cpatterns,rate0,is_reversible)

	@classmethod
	def parse_rate(cls,rs,is_intra=False):
		rss = rs.strip()
		expr = cls.parse_math_expr(rss)
		if len(expr) > 1:
			return Rate(Expression(expr),is_intra)
		else:
			return Rate(rs,is_intra)

	# needs to identify other user-defined functions + stuff in parse_math_expr
	@classmethod
	def parse_func(cls,line):
		sline = line.strip()
		s_char = ''
		for x in sline:
			if re.match('\s',x) or re.match('=',s):
				s_char = x
				break
		name,func = re.split(s_char,sline)
		if re.search('\(.\)',name): # a variable in between the parentheses means the function is local (not Kappa compatible)
			raise NotCompatibleException("Kappa functions cannot accommodate local functions:\n\t%s\n"%sline)
		p_func = cls.parse_math_expr(func)
		return Expression(name,p_func.asList())

	# needs to be able to identify built in functions, numbers, variables, (previously defined functions?)
	# functions are an alphanumeric string starting with a letter; they are preceded by an operator or parenthesis and encompass something in parentheses
	# parameters are also alphanumeric strings starting with a letter; they are preceded by operators or parentheses and succeeded by operators
	@staticmethod
	def parse_math_expr(line):

		point = Literal( "." )
		e = CaselessLiteral( "E" )
		fnumber = Combine( Word( "+-"+nums, nums ) + 
		                   Optional( point + Optional( Word( nums ) ) ) +
		                   Optional( e + Word( "+-"+nums, nums ) ) )
		ident = Word(alphas, alphas+nums+"_$")

		plus  = Literal( "+" )
		minus = Literal( "-" )
		mult  = Literal( "*" )
		div   = Literal( "/" )
		lpar  = Literal( "(" )
		rpar  = Literal( ")" )
		addop  = plus | minus
		multop = mult | div
		expop = Literal( "^" )
		pi    = CaselessLiteral( "PI" )

		expr = Forward()
		atom = (Optional("-") + ( pi | e | fnumber | ident + lpar + expr + rpar | ident) | ( lpar + expr + rpar )) 

		# by defining exponentiation as "atom [ ^ factor ]..." instead of "atom [ ^ atom ]...", we get right-to-left exponents, instead of left-to-righ
		# that is, 2^3^2 = 2^(3^2), not (2^3)^2.
		factor = Forward()
		factor << atom + ZeroOrMore( ( expop + factor ))

		term = factor + ZeroOrMore( ( multop + factor ))
		expr << term + ZeroOrMore( ( addop + term ))
		pattern = expr

		return pattern.parseString(line.strip())

def _is_number(n):
	try:
		float(n)
	except (ValueError, AttributeError):
		return False
	else:
		return True


# KaSim has: inf, [mod]
# BNGL has: hyperbolic funcs, inverse trig, abs
# Likely compatible: max, min, rint, if, sum, avg, boolean expressions

bngl_to_kappa_func_map = {
	
	'_pi': '[pi]',
	'_e': '[exp](1)',
	'sin': '[sin]',
	'cos': '[cos]',
	'tan': '[tan]',
	'sqrt': '[sqrt]',
	'time': '[T]',
	'log2': lambda s: '([log](%s)/[log](2))'%s,
	'log10': lambda s: '([log](%s)/[log](10))'%s,
	'ln': '[log]',
	'exp': '[exp]',
	'sinh': lambda s: '([exp](%s) - [exp](-%s))/2.0'%(s,s),
	'cosh': lambda s: '([exp](%s) + [exp](-%s))/2.0'%(s,s),
	'tanh': lambda s: '([exp](%s) - [exp](-%s))/([exp](%s) + [exp](-%s))'%(s,s,s,s),
	'asinh': lambda s: '[log](%s + [sqrt](%s^2 + 1))'%(s,s),
	'acosh': lambda s: '[log](%s + [sqrt](%s^2 - 1))'%(s,s),
	'atanh': lambda s: '0.5*[log]((1+%s)/(1-%s))'%(s,s)
}

# TO CHECK:
# calculate arctrig or trigh
# rint is round
bngl_other_builtin_funcs = set(['asin','acos','atan','rint','abs','if','min','max','sum','avg'])
bngl_binary_operators = set(['==','!=','>=','<='])

# not sure how boolean expressions work in Kappa
# overlap_binary_operators = set(['*','/','^','+','-','&&','||','<','>'])

kappa_other_builtin_funcs = set(['[Tsim]','[Tmax]','[E]','[E-]','[Emax]','[pp]','inf'])
# [int] is floor
kappa_other_operators = set(['[int]','[mod]','[max]','[min]','[not]','[true]','[false]'])
