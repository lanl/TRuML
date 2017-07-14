import rbexceptions
from deepdiff import DeepDiff
from objects import *
import pyparsing as pp


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
            f = open(file_name)
            d = f.readlines()
            f.close()
            self.lines = d
        else:
            self.lines = []


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
        # var_dict keeps track of read variables and what type they are
        # Variables can be constant values (c), patterns (p), or dynamic expressions (d)
        self.var_dict = {}
        self.num_anon_pats = 0

    def parse(self):
        cur_line = ''
        model = Model()
        for i, l in enumerate(self.lines):
            if re.search("\\\\\s*$", l):
                # Saves current line, stripping trailing and leading whitespace, continues to subsequent line
                cur_line += re.sub('\\\\', '', l.strip())
                continue
            else:
                cur_line += l.strip()
                if re.match('%init', cur_line):
                    model.add_init(self.parse_init(cur_line))
                elif re.match('%agent', cur_line):
                    model.add_molecule(self.parse_mtype(cur_line))
                elif re.match('%var', cur_line):
                    scur_line = re.split('\s+', cur_line)
                    name = scur_line[1].strip("'")
                    expr_list = KappaReader.parse_alg_expr(' '.join(scur_line[2:])).asList()
                    if self.var_contains_pattern(expr_list):
                        if len(expr_list) == 1:
                            model.add_obs(Observable(name, [self.parse_cpattern(expr_list[0].strip('|'))]))
                        else:
                            pat_dict, subst_expr_list = self.get_var_patterns(expr_list)
                            for p in pat_dict.keys():
                                model.add_obs(Observable(p, [self.parse_cpattern(pat_dict[p].strip('|'))]))
                            model.add_func(Function(name, Expression(subst_expr_list)))
                    elif self.var_is_dynamic_no_pat(expr_list):
                        model.add_func(Function(name, Expression(expr_list)))
                    else:
                        model.add_parameter(Parameter(name, Expression(expr_list)))

                elif re.match('%obs', cur_line):
                    model.add_obs(self.parse_obs(cur_line))
                elif re.search('@', cur_line):
                    model.add_rule(self.parse_rule(cur_line))
                else:
                    continue
                cur_line = ''

        return model

    def var_is_dynamic_no_pat(self, expr_list):
        for atom in expr_list:
            if re.match('\[T', atom) or self.var_contains_pattern(expr_list):
                return True
            for k, v in self.var_dict.iteritems():
                if re.match("'%s'" % k, atom) and (v == 'd' or v == 'p'):
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
    def parse_init(line):
        sline = re.split('\s+', line)
        amount = ' '.join(sline[1:-1])
        pattern = KappaReader.parse_cpattern(sline[-1])
        amount_is_number = True if is_number(amount) else False
        if not amount_is_number:
            amount = Expression(KappaReader.parse_alg_expr(amount))
        return InitialCondition(pattern, amount, amount_is_number)

    @staticmethod
    def parse_mtype(line):
        sline = re.split('\s+', line.strip())

        psplit = re.split('\(', sline[1])
        name = psplit[0]

        site_name_map = {}  # tracks conversion to kappa by mapping BNGL site names to Kappa site namess

        sites = re.split(',', psplit[1].strip(')'))
        site_defs = []
        for s in sites:
            site_split = re.split('~', s)
            site_name = site_split[0]
            site_defs.append(SiteDef(site_name, [] if len(site_split) == 1 else site_split[1:]))

        return MoleculeDef(name, site_defs, site_name_map, False)

    @staticmethod
    def parse_molecule(mstr):
        smstr = mstr.strip()
        msplit = re.split('\(', smstr)
        mname = msplit[0]
        if not re.match('[A-Za-z][-+\w]*\(.*\)\s*$', smstr):
            raise rbexceptions.NotAMoleculeException(smstr)
        sites = re.split(',', msplit[1].strip(')'))
        if not sites[0]:
            return Molecule(mname, [])
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
        return Molecule(mname, site_list)

    @staticmethod
    def parse_cpattern(s):
        mol_list = []
        in_par = 0
        cur_mol = ''
        for c in s:
            if re.match('\(', c):
                in_par += 1
            elif re.match('\)', c):
                in_par -= 1
            if re.match(',', c) and in_par == 0:
                mol_list.append(KappaReader.parse_molecule(cur_mol))
                cur_mol = ''
                continue
            cur_mol += c
        mol_list.append(KappaReader.parse_molecule(cur_mol))
        return CPattern(mol_list)

    @staticmethod
    def parse_rule(line):
        pass

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
        variable = pp.Combine(pp.Literal("'") + pp.Word(pp.alphanums + "_") + pp.Literal("'"))

        # patterns
        pattern = pp.Combine(pp.Literal("|") + pp.Word(pp.alphas, pp.alphanums + "_") + lpar + (pp.Empty() ^ pp.CharsNotIn(")(")) + rpar + pp.Literal("|"))

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
        super(Reader, self).__init__(file_name)
        self.is_def_block = False
        self.is_init_block = False
        self.is_param_block = False
        self.is_rule_block = False
        self.is_obs_block = False
        self.is_func_block = False

    # TODO implement as simple grammar
    def parse(self):
        """
        Function to parse BNGL model files

        This function assumes that the file has the molecule types block before the rules block
        """
        cur_line = ''  # used for line continuation
        model = Model()
        for i, l in enumerate(self.lines):
            if re.match('begin parameters', l):
                self.is_param_block = True
                continue
            elif re.match('end parameters'):
                self.is_param_block = False
                continue
            elif re.match('begin molecule types', l):
                self.is_def_block = True
                continue
            elif re.match('end molecule types', l):
                self.is_def_block = False
                continue
            elif re.match('begin seed species', l):
                self.is_init_block = True
                continue
            elif re.match('end seed species', l):
                self.is_init_block = False
                continue
            elif re.match('begin observables', l):
                self.is_obs_block = True
                continue
            elif re.match('end observables', l):
                self.is_obs_block = False
                continue
            elif re.match('begin functions', l):
                self.is_func_block = True
                continue
            elif re.match('end functions', l):
                self.is_func_block = False
                continue
            elif re.match('begin reaction rules', l):
                self.is_rule_block = True
                continue
            elif re.match('end reaction rules', l):
                self.is_rule_block = False
                continue

            # determines presence of line continuation, file cannot have backslashes in other contexts
            if re.search("\\\\\s*$", l):
                # saves current line, stripping trailing and leading whitespace, continues to subsequent line
                cur_line += re.sub('\\\\', '', l.strip())
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
    def parse_molecule(cls, mstr):
        """
        Function that parses molecules.

        Parameters
        ----------
        mstr : str
            String in BNGL file that represents a single molecule

        Returns
        -------
        Molecule
            Builds a Molecule or raises a NotAMoleculeException
        """
        smstr = mstr.strip()
        msplit = re.split('\(', smstr)
        mname = msplit[0]
        if not re.match('[A-Za-z]\w*\(.*\)\s*$', smstr):
            raise rbexceptions.NotAMoleculeException(smstr)
        sites = re.split(',', msplit[1].strip(')'))
        if not sites[0]:
            return Molecule(mname, [])
        site_list = []
        for i in range(len(sites)):
            s = sites[i]
            if '~' in s:
                tsplit = re.split('~', s)
                name = tsplit[0]
                if '!' in s:
                    bsplit = re.split('!', tsplit[1])
                    bond = cls.parse_bond(bsplit[1])
                    site_list.append(Site(name, i, s=bsplit[0], b=bond))
                else:
                    site_list.append(Site(name, i, s=tsplit[1]))
            else:
                if '!' in s:
                    bsplit = re.split('!', s)
                    name = bsplit[0]
                    bond = cls.parse_bond(bsplit[1])
                    site_list.append(Site(name, i, b=bond))
                else:
                    site_list.append(Site(s, i))
        return Molecule(mname, site_list)

    # TODO implement parsing for expression (need to identify variables for conversion to kappa syntax)
    @classmethod
    def parse_init(cls, line):
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
        spec = cls.parse_cpattern(isplit[0])
        amount = ' '.join(isplit[1:])
        amount_is_number = is_number(amount)
        p_amount = float(amount) if amount_is_number else Expression(cls.parse_math_expr(amount))
        return InitialCondition(spec, p_amount, amount_is_number)

    @classmethod
    def parse_cpattern(cls, pstr):
        """
        Function that parses patterns connected by the '.' operator

        Parameters
        ----------
        pstr : str
            String in BNGL file that represents a pattern

        Returns
        -------
        CPattern
        """
        spstr = re.split('(?<=\))\.', pstr.strip())
        m_list = []
        for s in spstr:
            m_list.append(cls.parse_molecule(s))
        return CPattern(m_list)

    @classmethod
    def parse_obs(cls, line):
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
        oCPattern = [cls.parse_cpattern(p) for p in osplit[2:]]
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
        s_char = ''
        for x in sline:
            if re.match('\s', x) or re.match('=', x):
                s_char = x
                break
        psplit = re.split(s_char, sline)
        pname = psplit[0]
        pexpr = s_char.join(psplit[1:])
        return Parameter(pname, pexpr)

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

    # TODO parse rule label, change so that lhs and rule 'action' is returned
    @classmethod
    def parse_rule(cls, line):
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
        rhs = ''
        lhs = ''
        is_reversible = True if re.search('<->', sline) else False
        parts = re.split('->', sline)
        lhs_cpatterns = [cls.parse_cpattern(x) for x in re.split('(?<!!)\+', parts[0].rstrip('<'))]
        rem = [x.strip() for x in re.split('(?<!!)\+', parts[1].strip())]
        # if the rule is an unbinding rule or has a '+' in its rate expression
        if len(rem) > 1:
            one_past_final_mol_index = 0
            for i, t in enumerate(rem):
                try:
                    cls.parse_cpattern(t)
                except rbexceptions.NotAMoleculeException:
                    one_past_final_mol_index = i
                    break
            last_split = re.split('\s+', rem[one_past_final_mol_index])
            mol, first_rate_part = last_split[0], ' '.join(last_split[1:])
            rhs_cpatterns = [cls.parse_cpattern(x) for x in (rem[:one_past_final_mol_index] + [mol])]
            rate_string = first_rate_part + '+' + '+'.join(rem[one_past_final_mol_index + 1:])
            if is_reversible:
                rate0, rate1 = re.split(',', rate_string)
                return Rule(lhs_cpatterns, rhs_cpatterns, cls.parse_rate(rate0), is_reversible, cls.parse_rate(rate1))
            else:
                return Rule(lhs_cpatterns, rhs_cpatterns, cls.parse_rate(rate_string), is_reversible)
        else:
            rem_parts = re.split('(?<!!)\s+', parts[1].strip())
            rhs_cpatterns = [cls.parse_cpattern(rem_parts[0])]
            is_intra_l_to_r = False
            if len(lhs_cpatterns) == 1 and len(rhs_cpatterns) == 1:
                is_intra_l_to_r = cls._has_intramolecular_binding(lhs_cpatterns[0], rhs_cpatterns[0])
            rate_string = ' '.join(rem_parts[1:])
            if is_reversible:
                is_intra_r_to_l = cls._has_intramolecular_binding(rhs_cpatterns[0], lhs_cpatterns[0])
                rate0_string, rate1_string = re.split(',', rate_string)
                rate0 = cls.parse_rate(rate0_string, is_intra_l_to_r)
                rate1 = cls.parse_rate(rate1_string, is_intra_r_to_l)
                return Rule(lhs_cpatterns, rhs_cpatterns, rate0, is_reversible, rate1)
            else:
                rate0 = cls.parse_rate(rate_string, is_intra_l_to_r)
                return Rule(lhs_cpatterns, rhs_cpatterns, rate0, is_reversible)

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
        s_char = ''
        for x in sline:
            if re.match('\s', x) or re.match('=', x):
                s_char = x
                break
        name, func = re.split(s_char, sline)
        if re.search('\(.\)',
                     name):  # a variable in between the parentheses means the function is local (not Kappa compatible)
            raise rbexceptions.NotCompatibleException("Kappa functions cannot accommodate local functions:\n\t%s\n" % sline)
        p_func = cls.parse_math_expr(func)
        return Expression(name, p_func.asList())

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
        pi = pp.CaselessLiteral("PI")

        expr = pp.Forward()
        atom = (pp.Optional("-") + (pi | e | fnumber | ident + lpar + expr + rpar | ident) | (lpar + expr + rpar))

        # by defining exponentiation as "atom [ ^ factor ]..." instead of "atom [ ^ atom ]...", we get right-to-left exponents, instead of left-to-righ
        # that is, 2^3^2 = 2^(3^2), not (2^3)^2.
        factor = pp.Forward()
        factor << atom + pp.ZeroOrMore((expop + factor))

        term = factor + pp.ZeroOrMore((multop + factor))
        expr << term + pp.ZeroOrMore((addop + term))
        pattern = expr

        return pattern.parseString(estr.strip())
