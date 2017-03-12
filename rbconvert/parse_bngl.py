import re
from pyparsing import Literal,CaselessLiteral,Word,Combine,Group,Optional,\
    ZeroOrMore,Forward,nums,alphas

class NotAMoleculeException(Exception):
	def __init__(self,s):
		print "%s is not a molecule"%s

class NotCompatibleException(Exception):
	def __init__(self,s):
		print s

class MoleculeDef:
	def __init__(self,n,sd):
		self.name = n
		self.sites = sd

	def _all_site_states(self):
		ss = []
		# ordered for deterministic behavior
		for k in sorted(self.sites.keys()):
			v = self.sites[k]
			if not v:
				ss.append(k)
			else:
				ss.append('%s~%s'%(k,'~'.join(v)))
		return ','.join(ss)

	def write_as_bngl(self):
		return "%s(%s)"%(self.name,self._all_site_states())

	def write_as_kappa(self):
		return "%%agent: %s(%s)"%(self.name,self._all_site_states())

class Molecule:
	def __init__(self,name,sites):
		self.name = name
		self.sites = sites # list of Sites
	
	def write_as_bngl(self):
		return self.name + '(' + ','.join([s.write_as_bngl() for s in self.sites]) + ')'

	def write_as_kappa(self):
		return self.name + '(' + ','.join([s.write_as_kappa() for s in self.sites]) + ')'		

	def __repr__(self):
		return 'Molecule(name: %s, sites: %s)'%(self.name,', '.join([str(x) for x in self.sites]))

class Site:
	def __init__(self,n,s=None,b=None):
		self.name = n
		self.state = s
		self.bond = b

	def write_as_bngl(self):
		s = self.name
		if self.state is not None:
			s += '~%s'%self.state
		if self.bond is not None:
			s += self.bond.write_as_bngl()
		return s

	def write_as_kappa(self):
		s = self.name
		if self.state is not None:
			s += '~%s'%self.state
		if self.bond is not None:
			s += self.bond.write_as_kappa()
		return s

	def __repr__(self):
		return 'Site(name: %s, state: %s, bond: %s)'%(self.name,self.state,self.bond)

# TODO make sure to include ! operator
class Bond:
	def __init__(self,n,w=False,a=False):
		self.wild = w
		self.num = int(n) # negative numbers indicate absence of specific bond, positive numbers will override w and a
		self.any = a
		if self.num < 0:
			assert(self.wild or self.any)

	def write_as_bngl(self):
		s = ''
		if self.wild:
			s = '!+'
		elif self.any:
			s = '!?'
		else:
			s = '!%s'%self.num
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

	def __repr__(self):
		b_string = str(self.num)
		if self.wild:
			b_string = 'wild'
		elif self.any:
			b_string = 'any'
		return b_string

class Pattern:
	def __init__(self,ml):
		self.molecule_list = ml

	def write_as_bngl(self):
		return '.'.join([m.write_as_bngl() for m in self.molecule_list])

	def write_as_kappa(self):
		return ','.join([m.write_as_kappa() for m in self.molecule_list])

	def __repr__(self):
		return '\n'.join([str(x) for x in self.molecule_list])

class InitialCondition:
	def __init__(self,s,a,ais=False):
		self.species = s
		self.amount = a # can be number or parameter string
		self.amount_is_parameter = ais

	def write_as_bngl(self):
		return '%s\t%s'%(self.species.write_as_bngl(),self.amount)

	def write_as_kappa(self):
		amount = self.amount if not self.amount_is_parameter else "'%s'"%self.amount
		return '%%init: %s %s'%(amount,self.species.write_as_kappa())

class Parameter:
	def __init__(self,n,v):
		self.name = n # string
		self.value = v # number

	def write_as_bngl(self):
		return '%s %s'%(self.name,self.value)

	def write_as_kappa(self):
		return '%%var: \'%s\' %s'%(self.name, self.value)

# special parsing required for 'if', 'log' functions
# can implement conversion of certain types of values (e.g. log10(x) to log(x)/log(10))
class Expression:
	def __init_(self,n,atom_list):
		self.name = n # assigned label (includes parentheses and stuff inside)
		self.atom_list = atom_list # list from _parse_math_expr listing (in order) operators, values, variables

	def write_as_bngl(self):
		return '%s'%''.join(atom_list)

	def write_as_kappa(self,as_obs=True):
		s = '%%%s: \'%s\' '%(self.name, 'var' if not as_obs else 'obs')
		expr = ''

		i = 0
		while (i < len(self.atom_list)):
			a = self.atom_list[i]
			if re.match('\w+',a) and a not in bngl_builtin_funcs:
				expr += '\'%s\''%a
			elif a in bngl_to_kappa_func_map.keys():
				trig_func_match = re.compile('sinh|cosh|tanh|asinh|acosh|atanh')
				if re.match('log',a) or re.match(trig_func_match,a):
					expr += bngl_to_kappa_func_map[a](self.atom_list[i+2])
					i += 4
				else:
					expr += bngl_to_kappa_func_map[a]
			else:
				expr += a
			i += 1

		return s + expr + '\n'

#TODO implement check for rate as raw number before writing
class Rule:
	# lhs, rhs are lists of Patterns, rate is parameter, number or expression, rev is bool (true for reversible rules)
	def __init__(self,lhs,rhs,rate,rev):
		self.lhs = lhs
		self.rhs = rhs
		self.rate = rate
		self.rev = rev
		self.arrow = '->' if not rev else '<->'

	def write_as_bngl(self):
		lhs_string = '+'.join([p.write_as_bngl() for p in self.lhs])
		rhs_string = '+'.join([p.write_as_bngl() for p in self.rhs])
		return '%s %s %s %s'%(lhs_string,self.arrow,rhs_string,rate.write_as_bngl())

	def write_as_kappa(self):
		lhs_string = ','.join([p.write_as_kappa() for p in self.lhs])
		rhs_string = ','.join([p.write_as_kappa() for p in self.rhs])
		return '%s %s %s @ %s'%(lhs_string,self.arrow,rhs_string,rate.write_as_kappa())


class Observable:
	def __init__(self,n,ps,t):
		self.name = n
		self.type = t
		assert(re.match('[sS]$',self.type) or re.match('[mM]$',self.type))
		self.patterns = ps # a list of patterns

	def write_as_bngl():
		return "%s %s %s"%(self.type,self.name,' '.join(self.patterns))

	def write_as_kappa():
		print "Kappa does not distinguish between BNGL \'Molecules\' and \'Species\' observables\nKappa's observables are analogous to the Molecules observable type"
		obs = ['+'.join(['|%s|'%p.write_as_kappa() for p in self.patterns])]
		return '%%obs: \'%s\' %s'%(self.name,obs)

# class Function:
# 	def __init__(self,n,f):
# 		self.name = n
# 		self.function = f

# 	def write_as_bngl():
# 		pass

# 	def write_as_kappa():
# 		pass

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
			s += '%s\n'%f.write_as_kappa(func_as_obs)
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
					model.add_parameter(_parse_param(cur_line))
				elif self.is_def_block:
					model.add_molecule(_parse_mtype(cur_line))
				elif self.is_init_block:
					model.add_init(_parse_init(cur_line))
				elif self.is_obs_block:
					model.add_obs(_parse_obs(cur_line))
				elif self.is_func_block:
					model.add_func(_parse_func(cur_line))
				elif self.is_rule_block:
					model.add_rule(_parse_rule(cur_line))
				else:
					continue
				cur_line = ''

		return model

	def _is_number(n):
		try:
			float(n)
		except ValueError:
			return False
		else:
			return True

	def _parse_bond(b):
		if re.match('\+',b):
			return Bond(-1,w=True)
		elif re.match('\?',b):
			return Bond(-1,a=True)
		elif b.isdigit():
			return Bond(b)
		else:
			raise ValueError("Illegal bond: %s"%b)

	def _parse_mtype(line):
		psplit = re.split('\(',line.strip())
		name = psplit[0]
		site_dict = {} # site name: list of possible states
		sites = re.split(',',psplit[1].strip(')'))
		for s in sites:
			site_split = re.split('~',s)
			if len(site_split) == 1:
				site_dict[site_split[0]] = []
			else:
				site_dict[site_split[0]] = site_split[1:]
		return MoleculeDef(name, site_dict)

	def _parse_molecule(line):
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
					bond = _parse_bond(bsplit[1])
					site_list.append(Site(name,s=bsplit[0],b=bond))
				else:
					site_list.append(Site(name,s=tsplit[1]))
			else:
				if '!' in s:
					bsplit = re.split('!',s)
					name = bsplit[0]
					bond = _parse_bond(bsplit[1])
					site_list.append(Site(name,b=bond))
				else:
					site_list.append(Site(s))
		return Molecule(mname,site_list)

	def _parse_init(line):
		isplit = split('\s+',line.strip())
		spec = _parse_species(isplit[0])
		amount = isplit[1]
		return InitialCondition(spec,amount,not _is_number(amount))

	def _parse_pattern(line):
		ssplit = re.split('(?<=\))\.',line.strip())
		m_list = []
		for s in ssplit:
			m_list.append(_parse_molecule(s))
		return Pattern(m_list)

	def _parse_obs(line):
		osplit = split('\s+',line.strip())
		otype = osplit[0][0]
		oname = osplit[1]
		opattern = [_parse_pattern(p) for p in osplit[2:]]
		return Observable(oname,opattern,otype)

	def _parse_param(line):
		sline = line.strip()
		s_char = ''
		for x in sline:
			if re.match('\s',x) or re.match('=',x):
				s_char = x
				break
		psplit = split(s_char,sline)
		pname = psplit[0]
		pexpr = s_char.join(psplit[1:])
		return Parameter(pname,pexpr)

	# TODO parse rule label
	def _parse_rule(line):
		sline = line.strip()
		rhs = ''
		lhs = ''
		is_reversible = True if re.search('<->',sline) else False
		parts = re.split('->',sline)
		lhs_patterns = [_parse_pattern(x) for x in re.split('(?<!!)\+',parts[0].rstrip('<'))]
		rem = re.split('(?<!!)\+',parts[1].strip())
		if len(rem) > 1:
			one_past_final_mol_index = 0
			for i,t in enumerate(rem):
				try:
					_parse_pattern(t)
				except NotAMoleculeException:
					one_past_final_mol_index = i
			mol,first_rate_part = re.split('\s+',rem[one_past_final_mol_index])
			rhs_patterns = [_parse_pattern(x) for x in (rem[:one_past_final_mol_index] + [mol])]
			rate = first_rate_part + '+'.join(rem[one_past_final_mol_index+1:])
		else:
			rem_parts = re.split('(?<!!)\s+',parts[1])
			rhs_patterns = list(_parse_pattern(rem_parts[0]))
			rate = ' '.join(rem_parts[1:])
		return lhs_patterns,rhs_patterns,rate
		return Rule(lhs_patterns,rhs,rate,is_reversible)

	# needs to identify other user-defined functions + stuff in parse_math_expr
	def _parse_func(line):
		sline = line.strip()
		s_char = ''
		for x in sline:
			if re.match('\s',x) or re.match('=',s):
				s_char = x
				break
		name,func = re.split(s_char,sline)
		if re.search('\(.\)',name): # a variable in between the parentheses means the function is local (not Kappa compatible)
			raise NotCompatibleException("Kappa functions cannot accommodate local functions:\n\t%s\n"%sline)
		p_func = _parse_math_expr(func)
		return Expression(name,p_func.asList())

	# needs to be able to identify built in functions, numbers, variables, (previously defined functions?)
	# functions are an alphanumeric string starting with a letter; they are preceded by an operator or parenthesis and encompass something in parentheses
	# parameters are also alphanumeric strings starting with a letter; they are preceded by operators or parentheses and succeeded by operators
	def _parse_math_expr(line):

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
	'log2': lambda s: '[log](%s)/[log](2)'%s,
	'log10': lambda s: '[log](%s)/[log](10)'%s,
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

