"""truml.readers: module containing classes for reading BNGL and Kappa model files"""

from deepdiff import DeepDiff
from objects import *

import itertools as it
import logging
import pyparsing as pp
import re
import utils


class Reader(object):
    """Inherited class for reading rule-based modeling files."""

    def __init__(self, file_name):
        """
        Reader initialization function

        Parameters
        ----------
        file_name : str
            Rule-based model file
        """
        self.file_name = file_name
        if self.file_name is not None:
            try:
                f = open(file_name)
                self.lines = f.readlines()
                f.close()
                logging.info("Read file %s" % self.file_name)
            except IOError:
                logging.error("Cannot find model file %s" % file_name)
                raise rbexceptions.NoModelsException("Cannot find model file %s" % file_name)
        else:
            self.lines = []

    @staticmethod
    def ignore_line(l):
        return l == '' or re.match('\s*\n', l)


# ignores perturbation and action commands
class KappaReader(Reader):
    """Reader for Kappa model files"""

    def __init__(self, file_name=None):
        """
        Kappa initialization function

        Parameters
        ----------
        file_name : str
        """
        super(KappaReader, self).__init__(file_name)
        self.lines = list(it.ifilterfalse(self.ignore_line, [re.sub('//.*$', '', l).strip() for l in self.lines]))
        # var_dict keeps track of read variables and observables and what type they are
        # Variables can be constant values (c), patterns (p), or dynamic expressions (d)
        self.var_dict = {}
        self.num_anon_pats = 0

    def get_agents(self):
        """Get molecule/agent definitions"""
        m = Model()
        for l in self.lines:
            logging.debug("Parsing: %s" % l.strip())
            if re.match('%agent', l):
                m.add_molecule_def(self.parse_mtype(l))
        return m

    def parse(self):
        # First get agent definitions
        model = self.get_agents()

        if len(model.molecules) == 0:
            logging.error("No molecules found")
            raise rbexceptions.NoMoleculesException("Did not find any molecules.  "
                                                    "Check file and TRuML command-line arguments")

        for i, l in enumerate(self.lines):
            if re.match('%agent', l):
                continue

            logging.debug("Parsing: %s" % l.strip())

            if re.match('%init', l):
                inits = self.parse_init(l, model.molecules)
                for init in inits:
                    model.add_init(init)
            elif re.match('%var', l) or re.match('%obs', l):
                match = re.match("%(obs|var):\s*('.*?')\s*(.*)$", l)
                name = match.group(2).strip("'")

                bname = name
                if re.search('\W', name):
                    bname = re.sub('\W', '_', name)
                    logging.warning(
                        "Exact conversion of observable '%s' to BNGL is not possible.  Renamed to '%s'" % (
                            name, bname))

                if bname in model.convert_namespace.values():
                    rebname = bname + '_'
                    logging.warning(
                        "Name '%s' already exists due to inexact conversion.  Renamed to '%s'" % (bname, rebname))
                    bname = rebname

                model.convert_namespace[name] = bname
                expr_list = KappaReader.parse_alg_expr(match.group(3).strip())

                if self.var_contains_pattern(expr_list):
                    if len(expr_list) == 1:
                        model.add_obs(Observable(name, self.parse_cpatterns(expr_list[0].strip('|'), model.molecules)))
                        self.var_dict[name] = 'p'
                    else:
                        pat_dict, subst_expr_list = self.get_var_patterns(expr_list)
                        for p in pat_dict.keys():
                            model.add_obs(Observable(p, self.parse_cpatterns(pat_dict[p].strip('|'), model.molecules)))
                            self.var_dict[p] = 'p'
                        model.add_func(Function(name, Expression(subst_expr_list)))
                        self.var_dict[name] = 'd'
                elif self.var_is_dynamic_no_pat(expr_list):
                    model.add_func(Function(name, Expression(expr_list)))
                    self.var_dict[name] = 'd'
                else:
                    model.add_parameter(Parameter(name, Expression(expr_list)))
                    self.var_dict[name] = 'c'
            elif re.search('@', l):
                rules = self.parse_rule(l, model.molecules)
                for rule in rules:
                    model.add_rule(rule)

        logging.info("Parsed Kappa model file %s" % self.file_name)
        return model

    def var_is_dynamic_no_pat(self, expr_list):
        for atom in expr_list:
            if re.match('\[T', atom):
                return True
            for k, v in self.var_dict.iteritems():
                if re.match("%s" % k, atom) and (v == 'd' or v == 'p'):
                    return True
        return False

    @staticmethod
    def var_contains_pattern(expr_list):
        for atom in expr_list:
            if re.match('\|', atom):
                return True
        return False

    def get_var_patterns(self, expr_list):
        pat_dict = {}
        new_expr_list = []
        for atom in expr_list:
            if re.match('\|', atom):
                anon = "anon_obs%s" % self.num_anon_pats
                self.num_anon_pats += 1
                pat_dict[anon] = atom
                new_expr_list.append(anon)
            else:
                new_expr_list.append(atom)
        return pat_dict, new_expr_list

    @staticmethod
    def parse_init(line, mdefs):
        sline = re.split('\s+', line)
        amount = ' '.join(sline[1:-1])
        patterns = KappaReader.parse_cpatterns(sline[-1], mdefs)
        amount_is_number = True if is_number(amount) else False
        if not amount_is_number:
            amount = Expression(KappaReader.parse_alg_expr(amount))
        return [InitialCondition(pattern, amount, amount_is_number) for pattern in patterns]

    @staticmethod
    def _get_sitedef(s):
        return SiteDef(s.name, s.states[:])

    @staticmethod
    def _get_molecdef(s):
        site_defs = s.sites[:]
        site_name_map = {x.name: x.name for x in site_defs}
        return MoleculeDef(s.name, site_defs, site_name_map)

    @staticmethod
    def parse_mtype(line):
        mol_string = re.sub('\s*%agent:\s*', '', line)

        name = pp.Word(pp.alphas, bodyChars=pp.alphanums+"_+-") | pp.Word("_", bodyChars=pp.alphanums+"_+-")

        lpar = pp.Suppress(pp.Literal("("))
        rpar = pp.Suppress(pp.Literal(")"))
        lbrace = pp.Suppress(pp.Literal("{"))
        rbrace = pp.Suppress(pp.Literal("}"))
        comma = pp.Suppress(pp.Literal(","))

        site = name.setResultsName('name') + \
               pp.Optional(lbrace + pp.Group(name + pp.ZeroOrMore(pp.Optional(comma) + name)).setResultsName('states') +
                           rbrace)
        site.setParseAction(KappaReader._get_sitedef)
        agent = name.setResultsName('name') + lpar + \
                pp.Group(pp.Optional(site + pp.ZeroOrMore(pp.Optional(comma) + site))).setResultsName('sites') + rpar
        agent.setParseAction(KappaReader._get_molecdef)

        return agent.parseString(mol_string)[0]

    @staticmethod
    def parse_molecule(mstr, mdefs):
        smstr = mstr.strip()
        msplit = re.split('\(', smstr)
        mname = msplit[0]
        mtype = None
        for mdef in mdefs:
            if mdef.name == mname:
                mtype = mdef
        if mtype is None:
            raise rbexceptions.UnknownMoleculeTypeException(mname)

        if not re.match('[A-Za-z][-+\w]*\(.*\)\s*$', smstr):
            raise rbexceptions.NotAMoleculeException(smstr)
        sites = re.split(',', msplit[1].strip(')'))
        if not sites[0]:
            return Molecule(mname, [], mtype)
        site_list = []
        for i in range(len(sites)):
            s = sites[i]
            if '~' in s:
                tsplit = re.split('~', s)
                name = tsplit[0]
                if '!' in s:
                    bsplit = re.split('!', tsplit[1])
                    bond = Bond(-1, w=True) if re.match('_', bsplit[1]) else Bond(int(bsplit[1]))
                    site_list.append(Site(name, i, s=bsplit[0], b=bond))
                elif re.search('\?$', s):
                    bond = Bond(-1, a=True)
                    site_list.append(Site(name, i, s=tsplit[1].strip('?'), b=bond))
                else:
                    site_list.append(Site(name, i, s=tsplit[1]))
            else:
                if '!' in s:
                    bsplit = re.split('!', s)
                    name = bsplit[0]
                    bond = Bond(-1, w=True) if re.match('_', bsplit[1]) else Bond(int(bsplit[1]))
                    site_list.append(Site(name, i, b=bond))
                elif re.search('\?$', s):
                    bond = Bond(-1, a=True)
                    site_list.append(Site(s.strip('?'), i, b=bond))
                else:
                    site_list.append(Site(s, i))
        return Molecule(mname, site_list, mtype)

    @staticmethod
    def parse_cpatterns(s, mdefs):
        mol_list = []
        in_par = 0
        cur_mol = ''
        for c in s:
            if re.match('\(', c):
                in_par += 1
            elif re.match('\)', c):
                in_par -= 1
            if re.match(',', c) and in_par == 0:
                mol_list.append(KappaReader.parse_molecule(cur_mol, mdefs))
                cur_mol = ''
                continue
            cur_mol += c
        mol_list.append(KappaReader.parse_molecule(cur_mol, mdefs))
        conn_cmps = utils.get_connected_components(mol_list)
        return CPatternList([CPattern(c) for c in conn_cmps])

    @staticmethod
    def parse_rule(line, mdefs):
        reversible = False
        rule_str = line

        label = None
        label_match = re.match("\s*'(.*?)'", line)
        if label_match:
            label = label_match.group(1).strip("'")
            rule_str = line[label_match.end():].strip()

        rule_cps = re.split('->', rule_str)

        # Check to see if the rule is reversible
        if re.search('<$', rule_cps[0]):
            reversible = True
            rule_cps[0] = rule_cps[0].strip('<')

        if rule_cps[0].strip() == '':
            lhs_patts = []
        else:
            lhs_patts = KappaReader.parse_cpatterns(rule_cps[0], mdefs)

        rhs_cps = re.split('@', rule_cps[1])
        if rhs_cps[0].strip() == '':
            rhs_patts = []
        else:
            rhs_patts = KappaReader.parse_cpatterns(rhs_cps[0].strip(), mdefs)

        n_lhs_mols = sum([p.num_molecules() for p in lhs_patts])
        n_rhs_mols = sum([p.num_molecules() for p in rhs_patts])
        delmol = n_lhs_mols > n_rhs_mols
        if delmol:
            logging.debug("Rule '%s' is a degradation rule" % line)

        if reversible:
            rate_cps = re.split(',', rhs_cps[1])
            if re.search('{', rate_cps[0]) and re.search('{', rate_cps[1]):
                frate_cps = re.split('{', rate_cps[0].strip())
                rrate_cps = re.split('{', rate_cps[1].strip())
                inter_frate = Rate(Expression(KappaReader.parse_alg_expr(frate_cps[0].strip()).asList()))
                intra_frate = Rate(Expression(KappaReader.parse_alg_expr(frate_cps[1].strip()).asList()), True)
                inter_rrate = Rate(Expression(KappaReader.parse_alg_expr(rrate_cps[0].strip()).asList()))
                intra_rrate = Rate(Expression(KappaReader.parse_alg_expr(rrate_cps[1].strip('}').strip()).asList()),
                                   True)
                inter_frule = Rule(lhs_patts, rhs_patts, inter_frate, label=label, delmol=delmol)
                intra_frule = Rule(lhs_patts, rhs_patts, intra_frate, label=label, delmol=delmol)
                inter_rrule = Rule(rhs_patts, lhs_patts, inter_rrate, label=label, delmol=delmol)
                intra_rrule = Rule(rhs_patts, lhs_patts, intra_rrate, label=label, delmol=delmol)
                return [inter_frule, intra_frule, inter_rrule, intra_rrule]
            elif re.search('{', rate_cps[0]):
                frate_cps = re.split('{', rate_cps[0].strip())
                inter_frate = Rate(Expression(KappaReader.parse_alg_expr(frate_cps[0].strip()).asList()))
                intra_frate = Rate(Expression(KappaReader.parse_alg_expr(frate_cps[1].strip()).asList()), True)
                rrate = Rate(Expression(KappaReader.parse_alg_expr(rate_cps[1].strip()).asList()))
                inter_frule = Rule(lhs_patts, rhs_patts, inter_frate, label=label, delmol=delmol)
                intra_frule = Rule(lhs_patts, rhs_patts, intra_frate, label=label, delmol=delmol)
                rrule = Rule(rhs_patts, lhs_patts, rrate, label=label, delmol=delmol)
                return [inter_frule, intra_frule, rrule]
            elif re.search('{', rate_cps[1]):
                rrate_cps = re.split('{', rate_cps[1].strip())
                inter_rrate = Rate(Expression(KappaReader.parse_alg_expr(rrate_cps[0].strip()).asList()))
                intra_rrate = Rate(Expression(KappaReader.parse_alg_expr(rrate_cps[1].strip('}').strip()).asList()),
                                   True)
                frate = Rate(Expression(KappaReader.parse_alg_expr(rate_cps[0].strip()).asList()))
                inter_rrule = Rule(rhs_patts, lhs_patts, inter_rrate, label=label, delmol=delmol)
                intra_rrule = Rule(rhs_patts, lhs_patts, intra_rrate, label=label, delmol=delmol)
                frule = Rule(lhs_patts, rhs_patts, frate, label=label, delmol=delmol)
                return [inter_rrule, intra_rrule, frule]
            else:
                frate = Rate(Expression(KappaReader.parse_alg_expr(rate_cps[0].strip()).asList()))
                rrate = Rate(Expression(KappaReader.parse_alg_expr(rate_cps[1].strip()).asList()))
                return [Rule(lhs_patts, rhs_patts, frate, reversible, rrate, label=label, delmol=delmol)]
        else:
            rate = Rate(Expression(KappaReader.parse_alg_expr(rhs_cps[1].strip())))
            return [Rule(lhs_patts, rhs_patts, rate, label=label, delmol=delmol)]

    @staticmethod
    def parse_alg_expr(estr):
        point = pp.Literal(".")
        e = pp.CaselessLiteral("E")
        fnumber = pp.Combine(pp.Word("+-" + pp.nums, pp.nums) +
                             pp.Optional(point + pp.Optional(pp.Word(pp.nums))) +
                             pp.Optional(e + pp.Word("+-" + pp.nums, pp.nums)))

        # infix operators
        plus = pp.Literal("+")
        minus = pp.Literal("-")
        mult = pp.Literal("*")
        div = pp.Literal("/")
        mod = pp.Literal("[mod]")
        lpar = pp.Literal("(")
        rpar = pp.Literal(")")
        expop = pp.Literal("^")

        addop = plus | minus
        multop = mult | div | mod

        # constants
        inf = pp.Literal("inf")
        pi = pp.Literal("[pi]")
        events = pp.Literal("[E]")
        null_events = pp.Literal("[E-]")
        event_limit = pp.Literal("[Emax]")
        time = pp.Literal("[T]")
        cpu_time = pp.Literal("[Tsim]")
        time_limit = pp.Literal("[Tmax]")
        plot_points = pp.Literal("[pp]")

        constant = inf | pi | events | null_events | event_limit | time | cpu_time | time_limit | plot_points

        # variables
        variable = pp.QuotedString("'")

        bondvalue = pp.Word(pp.nums) | pp.Literal("_")
        kbond = pp.Literal("?") | pp.Combine(pp.Literal("!") + bondvalue)
        kstate = pp.Word("~", pp.alphanums)
        site = pp.Combine(pp.Word(pp.alphas, pp.alphanums + "_") + pp.Optional(kstate) + pp.Optional(kbond))
        siteList = pp.delimitedList(site, delim=',', combine=True)
        mol = pp.Combine(pp.Word(pp.alphas, pp.alphanums + "_") + lpar + (pp.Empty() ^ siteList) + rpar)
        molList = pp.delimitedList(mol, delim=',', combine=True)
        # patterns
        pattern = pp.Combine(
            pp.Literal("|") + molList + pp.Literal("|"))

        # unary functions (one arg)
        logfunc = pp.Literal("[log]")
        expfunc = pp.Literal("[exp]")
        sinfunc = pp.Literal("[sin]")
        cosfunc = pp.Literal("[cos]")
        tanfunc = pp.Literal("[tan]")
        sqrtfunc = pp.Literal("[sqrt]")
        floorfunc = pp.Literal("[int]")

        unary_one_funcs = logfunc | expfunc | sinfunc | cosfunc | tanfunc | sqrtfunc | floorfunc

        # unary functions (two args)
        maxfunc = pp.Literal("[max]")
        minfunc = pp.Literal("[min]")

        unary_two_funcs = maxfunc | minfunc

        expr = pp.Forward()
        atom = (pp.Optional("-") + (
            constant | variable | fnumber | lpar + expr + rpar | unary_one_funcs + expr | unary_two_funcs + expr + expr | pattern))

        factor = pp.Forward()
        factor << atom + pp.ZeroOrMore((expop + factor))

        term = factor + pp.ZeroOrMore((multop + factor))
        expr << term + pp.ZeroOrMore((addop + term))
        fullExpr = expr

        return fullExpr.parseString(estr.strip())


# ignores action commands
class BNGLReader(Reader):
    """Reader for BNGL model files"""

    def __init__(self, file_name=None):
        """
        BNGLReader initialization function

        Parameters
        ----------
        file_name : str
        """
        super(BNGLReader, self).__init__(file_name)
        self.lines = self.condense_line_continuation(
            it.ifilterfalse(self.ignore_line, [re.sub('#.*$', '', l).strip() for l in self.lines]))

        self.is_def_block = False
        self.is_init_block = False
        self.is_param_block = False
        self.is_rule_block = False
        self.is_obs_block = False
        self.is_func_block = False

    @staticmethod
    def condense_line_continuation(lines):
        condensed_lines = []
        cur_line = ''
        for l in lines:
            if re.search("\\\\\s*$", l):
                # Saves current line, stripping trailing and leading whitespace, continues to subsequent line
                cur_line += re.sub('\\\\', '', l.strip())
                continue
            else:
                cur_line += l.strip()
                condensed_lines.append(cur_line)
                cur_line = ''
        return condensed_lines

    def get_molecule_types(self):
        m = Model()
        for l in self.lines:
            if re.match('begin molecule types', l):
                logging.debug("Entering molecule types block")
                self.is_def_block = True
                continue
            elif re.match('end molecule types', l):
                logging.debug("Leaving molecule types block")
                self.is_def_block = False
                break
            if self.is_def_block:
                m.add_molecule_def(self.parse_mtype(l))
        return m

    def parse(self):
        """
        Function to parse BNGL model files

        This function assumes that the file has the molecule types block before the rules block
        """
        model = self.get_molecule_types()

        if len(model.molecules) == 0:
            logging.error("No molecules found")
            raise rbexceptions.NoMoleculesException("Did not find any molecules.  "
                                                    "Check file and TRuML command-line arguments")

        for i, l in enumerate(self.lines):

            if re.match('\s*\n', l):
                continue

            logging.debug("Parsing %s" % l.strip())

            if re.match('begin\s+parameters', l):
                logging.debug("Entering parameter block")
                self.is_param_block = True
                continue
            elif re.match('end\s+parameters', l):
                logging.debug("Leaving parameter block")
                self.is_param_block = False
                continue
            elif re.match('begin\s+(seed\s+)?species', l):
                logging.debug("Entering initial conditions block")
                self.is_init_block = True
                continue
            elif re.match('end\s+(seed\s+)?species', l):
                logging.debug("Leaving initial conditions block")
                self.is_init_block = False
                continue
            elif re.match('begin\s+observables', l):
                logging.debug("Entering observables block")
                self.is_obs_block = True
                continue
            elif re.match('end\s+observables', l):
                logging.debug("Leaving observables block")
                self.is_obs_block = False
                continue
            elif re.match('begin\s+functions', l):
                logging.debug("Entering functions block")
                self.is_func_block = True
                continue
            elif re.match('end\s+functions', l):
                logging.debug("Leaving functions block")
                self.is_func_block = False
                continue
            elif re.match('begin\s+reaction\s+rules', l):
                logging.debug("Entering rules block")
                self.is_rule_block = True
                continue
            elif re.match('end\s+reaction\s+rules', l):
                logging.debug("Leaving rules block")
                self.is_rule_block = False
                continue

            if self.is_param_block:
                model.add_parameter(self.parse_param(l))
            elif self.is_init_block:
                model.add_init(self.parse_init(l, model.molecules))
            elif self.is_obs_block:
                obs = self.parse_obs(l, model.molecules)
                if obs is not None:
                    model.add_obs(self.parse_obs(l, model.molecules))
            elif self.is_func_block:
                model.add_func(self.parse_func(l))
            elif self.is_rule_block:
                model.add_rule(self.parse_rule(l, model.molecules))
            else:
                continue

        logging.info("Parsed BNGL model file: %s" % self.file_name)
        return model

    @staticmethod
    def parse_bond(b):
        """
        Function that parses bonds

        Parameters
        ----------
        b : str
            BNGL string that represents a bond

        Returns
        -------
        Bond
            Converts BNGL string to Bond instance. Raises ValueError if the string
            is malformed.
        """
        if re.match('\+', b):
            return Bond(-1, w=True)
        elif re.match('\?', b):
            return Bond(-1, a=True)
        elif b.isdigit():
            return Bond(b)
        else:
            raise ValueError("Illegal bond: %s" % b)

    @staticmethod
    def parse_mtype(line):
        """
        Function that parses molecule type definitions

        Parameters
        ----------
        line : str
            Line from BNGL file that represents a molecule type definition

        Returns
        -------
        MoleculeDef
            Builds MoleculeDef
        """
        psplit = re.split('\(', line.strip())
        name = psplit[0]

        site_name_map = {}  # tracks conversion to kappa by mapping BNGL site names to Kappa site namess

        sites = re.split(',', psplit[1].strip(')'))
        site_defs = []
        site_name_counter = {}
        has_site_symmetry = False
        for s in sites:
            site_split = re.split('~', s)
            site_name = site_split[0]
            site_defs.append(SiteDef(site_name, [] if len(site_split) == 1 else site_split[1:]))
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
                site_name_map[sn + str(site_name_counter[sn] - 1)] = sn
                site_name_counter[sn] -= 1

        return MoleculeDef(name, site_defs, site_name_map, has_site_symmetry)

    @classmethod
    def parse_molecule(cls, mstr, mdefs):
        """
        Function that parses molecules.

        Parameters
        ----------
        mstr : str
            String in BNGL file that represents a single molecule
        mdefs : list of MoleculeDef instances

        Returns
        -------
        Molecule
            Builds a Molecule or raises a NotAMoleculeException
        """
        smstr = mstr.strip()
        msplit = re.split('\(', smstr)
        mname = msplit[0]
        mtype = None
        for mdef in mdefs:
            if mdef.name == mname:
                mtype = mdef
        if mtype is None:
            raise rbexceptions.UnknownMoleculeTypeException(mname)

        if not re.match('[A-Za-z]\w*\(.*\)\s*$', smstr):
            raise rbexceptions.NotAMoleculeException(smstr)
        sites = re.split(',', msplit[1].strip(')'))
        if not sites[0]:
            return Molecule(mname, [], mtype)
        site_list = []
        for i in range(len(sites)):
            s = sites[i]
            if '~' in s:
                tsplit = re.split('~', s)
                name = tsplit[0]
                if '!' in s:
                    bsplit = re.split('!', tsplit[1])
                    bond = cls.parse_bond(bsplit[1])
                    state = None if re.match('\?', bsplit[0]) else bsplit[0]
                    site_list.append(Site(name, i, s=state, b=bond))
                else:
                    state = None if re.match('\?', tsplit[1]) else tsplit[1]
                    site_list.append(Site(name, i, s=state))
            else:
                if '!' in s:
                    bsplit = re.split('!', s)
                    name = bsplit[0]
                    bond = cls.parse_bond(bsplit[1])
                    site_list.append(Site(name, i, b=bond))
                else:
                    site_list.append(Site(s, i))
        return Molecule(mname, site_list, mtype)

    @classmethod
    def parse_init(cls, line, mdefs):
        """
        Function that parses initial conditions

        Parameters
        ----------
        line : str
            Line in BNGL file that represents an initial condition

        Returns
        -------
        InitialCondition
        """
        isplit = re.split('\s+', line.strip())
        spec = cls.parse_cpattern(isplit[0], mdefs)
        amount = ' '.join(isplit[1:])
        amount_is_number = is_number(amount)
        p_amount = float(amount) if amount_is_number else Expression(cls.parse_math_expr(amount))
        return InitialCondition(spec, p_amount, amount_is_number)

    @classmethod
    def parse_cpattern(cls, pstr, mdefs):
        """
        Function that parses patterns connected by the '.' operator

        Parameters
        ----------
        pstr : str
            String in BNGL file that represents a pattern
        mdefs : list of MoleculeDef instance

        Returns
        -------
        CPattern
        """
        spstr = re.split('(?<=\))\.', pstr.strip())
        m_list = []
        for s in spstr:
            m_list.append(cls.parse_molecule(s, mdefs))
        return CPattern(m_list)

    @classmethod
    def parse_obs(cls, line, mdefs):
        """
        Function that parses observables

        Parameters
        ----------
        line : str
            Line in BNGL file that represents an observable

        Returns
        -------
        Observable
        """
        osplit = re.split('\s+', line.strip())
        otype = osplit[0][0]
        oname = osplit[1]
        if re.search('==', line) and re.match('[sS]', otype):
            logging.warning("Observable '%s' cannot be converted to Kappa due to '==' syntax" % line.strip())
        else:
            oCPattern = CPatternList([cls.parse_cpattern(p, mdefs) for p in osplit[2:]])
            return Observable(oname, oCPattern, otype)

    @staticmethod
    def parse_param(line):
        """
        Function that parses parameters

        Parameters
        ----------
        line : str
            Line in BNGL file that represents a parameter

        Returns
        -------
        Parameter
        """
        sline = line.strip()
        if re.search('=', sline):
            s_char = '='
        else:
            s_char = re.compile('\s+')
        psplit = re.split(s_char, sline)
        pname = psplit[0].strip()
        pexpr = ' '.join(psplit[1:]).strip()
        if re.search('[*/+-]', pexpr) or re.match('[A-Za-z]', pexpr):
            pval = Expression(BNGLReader.parse_math_expr(pexpr))
        else:
            pval = pexpr
        return Parameter(pname, pval)

    # assumes that pattern mapping is left to right and that there is
    # only 1 component on either side of the rule (doesn't make sense to
    # have components that aren't operated on).  The change will be from
    # a Site with bond = None to a Site with a Bond object containing a
    # link to another Molecule in the same component
    @staticmethod
    def _has_intramolecular_binding(lhs_cp, rhs_cp):
        """
        Function that determines whether or not there is intramolecular binding

            This assumes that pattern mapping is left to right and that there is
            only 1 component on either side of the rule, since it doesn't make sense to
            have components that aren't operated on.  The change will be from
            a Site with bond = None to a Site with a Bond object containing a
            link to another Molecule in the same component

        Parameters
        ----------
        lhs_cp : CPattern
            Rule's left-hand side (the reactants)
        rhs_cp : CPattern
            Rule's right-hand side (the products)

        Returns
        -------
        bool
            True if the number of bonds changed is two (intramolecular bond formation),
            False if not
        """
        d = DeepDiff(lhs_cp, rhs_cp)
        try:
            changed = d.get('type_changes').keys()
        except AttributeError:
            return False
        num_changed_bonds = 0
        for c in changed:
            if re.search('sites\[.\]\.bond$', c):
                num_changed_bonds += 1
        return num_changed_bonds == 2

    @classmethod
    def parse_rule(cls, line, mdefs):
        """
        Function that parses rules

        Parameters
        ----------
        line : str
            Line in BNGL file that represents a reaction rule

        Returns
        -------
        Rule
        """
        sline = line.strip()

        sp = re.split(':', sline)
        label = None
        if len(sp) > 1:
            label = sp[0].strip()
            sline = sp[1].strip()

        is_reversible = True if re.search('<->', sline) else False
        parts = re.split('->', sline)

        lhs = re.split('(?<!!)\+', parts[0].rstrip('<'))
        if len(lhs) == 1 and lhs[0].strip() == '0':
            lhs_cpatterns = CPatternList([])
        else:
            lhs_cpatterns = CPatternList([cls.parse_cpattern(x, mdefs) for x in lhs])
        rem = [x.strip() for x in re.split('(?<!!)\+', parts[1].strip())]

        def del_mol_warning(s):
            if not re.search('DeleteMolecules', s):
                logging.warning(
                    "Degradation rule '%s' will remove full complexes.  This cannot be exactly translated into Kappa" % sline)
                logging.warning(
                    "Writing this rule in Kappa will only remove the matched pattern and could result in side effects")

        if len(rem) > 1:
            one_past_final_mol_index = 0
            for i, t in enumerate(rem):
                try:
                    cls.parse_cpattern(t, mdefs)
                except rbexceptions.NotAMoleculeException:
                    one_past_final_mol_index = i
                    break
            last_split = re.split('\s+', rem[one_past_final_mol_index])
            mol, first_rate_part = last_split[0], ' '.join(last_split[1:])

            if re.match('0', rem[0]):
                rhs_cpatterns = CPatternList([])
                del_mol_warning(rem[-1])
                rem[-1] = re.sub('\s*DeleteMolecules', '', rem[-1])
                delmol = True
            else:
                rhs_cpatterns = CPatternList([cls.parse_cpattern(x, mdefs) for x in (rem[:one_past_final_mol_index] + [mol])])
                n_lhs_mols = sum([p.num_molecules() for p in lhs_cpatterns])
                n_rhs_mols = sum([p.num_molecules() for p in rhs_cpatterns])
                delmol = n_lhs_mols > n_rhs_mols
                if delmol:
                    del_mol_warning(rem[-1])
                rem[-1] = re.sub('\s*DeleteMolecules', '', rem[-1])

            if len(rem[one_past_final_mol_index + 1:]) == 0:
                rate_string = first_rate_part
            else:
                rate_string = first_rate_part + '+' + '+'.join(rem[one_past_final_mol_index + 1:])

            if is_reversible:
                rate0, rate1 = re.split(',\s*', rate_string)
                return Rule(lhs_cpatterns, rhs_cpatterns, cls.parse_rate(rate0), is_reversible, cls.parse_rate(rate1),
                            label=label, delmol=delmol)
            else:
                return Rule(lhs_cpatterns, rhs_cpatterns, cls.parse_rate(rate_string), is_reversible, label=label,
                            delmol=delmol)
        else:
            rem_parts = re.split('(?<!!)\s+', parts[1].strip())

            if re.match('0', rem_parts[0]):
                rhs_cpatterns = CPatternList([])
                del_mol_warning(rem_parts[-1])
                rem_parts[-1] = re.sub('\s*DeleteMolecules', '', rem_parts[-1])
                delmol = True
            else:
                rhs_cpatterns = CPatternList([cls.parse_cpattern(rem_parts[0], mdefs)])
                n_lhs_mols = sum([p.num_molecules() for p in lhs_cpatterns])
                n_rhs_mols = sum([p.num_molecules() for p in rhs_cpatterns])
                delmol = n_lhs_mols > n_rhs_mols
                if delmol:
                    del_mol_warning(rem_parts[-1])
                rem_parts[-1] = re.sub('\s*DeleteMolecules', '', rem_parts[-1])

            rem_parts = [x for x in rem_parts if x != '']

            is_intra_l_to_r = False
            if len(lhs_cpatterns) == 1 and len(rhs_cpatterns) == 1:
                is_intra_l_to_r = cls._has_intramolecular_binding(lhs_cpatterns[0], rhs_cpatterns[0])
            rate_string = ' '.join(rem_parts[1:])
            if is_reversible:
                is_intra_r_to_l = False
                if len(lhs_cpatterns) == 1 and len(rhs_cpatterns) == 1:
                    is_intra_r_to_l = cls._has_intramolecular_binding(rhs_cpatterns[0], lhs_cpatterns[0])
                rate_split = re.split(',\s*', rate_string)
                rate0 = cls.parse_rate(rate_split[0], is_intra_l_to_r)
                rate1 = cls.parse_rate(rate_split[1], is_intra_r_to_l)
                return Rule(lhs_cpatterns, rhs_cpatterns, rate0, is_reversible, rate1, label=label, delmol=delmol)
            else:
                rate0 = cls.parse_rate(rate_string, is_intra_l_to_r)
                return Rule(lhs_cpatterns, rhs_cpatterns, rate0, is_reversible, label=label, delmol=delmol)

    @classmethod
    def parse_rate(cls, rs, is_intra=False):
        """
        Function for parsing rates

        Parameters
        ----------
        rs : str
            String in BNGL file corresponding to rate
        is_intra : bool
            True if the rate string is an intramolecular association rate, False otherwise

        Returns
        -------
        Rate
        """
        rss = rs.strip()
        expr = cls.parse_math_expr(rss)
        if len(expr) > 1:
            return Rate(Expression(expr), is_intra)
        else:
            return Rate(rs, is_intra)

    # needs to identify other user-defined functions + stuff in parse_math_expr
    @classmethod
    def parse_func(cls, line):
        """
        Function to parse BNGL functions

        Parameters
        ----------
        line : str
            Line in BNGL file that represents a function (i.e. a dynamic quantity)

        Returns
        -------
        Expression
        """
        sline = line.strip()
        if re.search('=', sline):
            s_char = '='
        else:
            s_char = ' '
        ssplit = re.split(s_char, sline)
        name = ssplit[0].strip()
        func = s_char.join(ssplit[1:]).strip()
        if re.search('\(.+\)', name):  # a variable in between the parentheses means the function is local
            raise rbexceptions.NotCompatibleException(
                "Kappa functions cannot accommodate BNGL local functions:\n\t%s\n" % sline)
        p_func = cls.parse_math_expr(func)
        return Function(name, Expression(p_func.asList()))

    # needs to be able to identify built in functions, numbers, variables, (previously defined functions?)
    # functions are an alphanumeric string starting with a letter; they are preceded by an operator or parenthesis and encompass something in parentheses
    # parameters are also alphanumeric strings starting with a letter; they are preceded by operators or parentheses and succeeded by operators
    @staticmethod
    def parse_math_expr(estr):
        """
        Function to parse algebraic expressions

        Parameters
        ----------
        estr : str
            String in BNGL file corresponding to an algebraic expression

        Returns
        -------
        list
            List of algebraic tokens, including functions, variables, numbers, and operators
        """

        point = pp.Literal(".")
        e = pp.CaselessLiteral("E")
        fnumber = pp.Combine(pp.Word("+-" + pp.nums, pp.nums) +
                             pp.Optional(point + pp.Optional(pp.Word(pp.nums))) +
                             pp.Optional(e + pp.Word("+-" + pp.nums, pp.nums)))
        ident = pp.Word(pp.alphas, pp.alphas + pp.nums + "_$")

        plus = pp.Literal("+")
        minus = pp.Literal("-")
        mult = pp.Literal("*")
        div = pp.Literal("/")
        lpar = pp.Literal("(")
        rpar = pp.Literal(")")
        addop = plus | minus
        multop = mult | div
        expop = pp.Literal("^")
        lowerPi = pp.Literal("pi")
        upperPi = pp.Literal("PI")
        pi = lowerPi | upperPi

        func = pp.Combine(ident + lpar + rpar)

        expr = pp.Forward()
        atom = (pp.Optional("-") + (pi ^ e ^ fnumber ^ ident + lpar + expr + rpar ^ func ^ ident) ^ (lpar + expr + rpar))

        # by defining exponentiation as "atom [ ^ factor ]..." instead of "atom [ ^ atom ]...", we get right-to-left exponents, instead of left-to-righ
        # that is, 2^3^2 = 2^(3^2), not (2^3)^2.
        factor = pp.Forward()
        factor << atom + pp.ZeroOrMore((expop + factor))

        term = factor + pp.ZeroOrMore((multop + factor))
        expr << term + pp.ZeroOrMore((addop + term))
        pattern = expr

        return pattern.parseString(estr.strip())
