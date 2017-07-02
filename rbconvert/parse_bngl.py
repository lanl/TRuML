import re
from pyparsing import Literal, CaselessLiteral, Word, Combine, Optional, \
    ZeroOrMore, Forward, nums, alphas
from deepdiff import DeepDiff
import itertools as it
import networkx as nx
import networkx.algorithms.isomorphism as iso


class NotAMoleculeException(Exception):
    """Raised when a string is expected to, but does not, conform to molecule syntax"""

    def __init__(self, s):
        print "%s is not a molecule" % s


class NotCompatibleException(Exception):
    """Raised when a BNGL string or model cannot be converted to Kappa"""

    def __init__(self, s):
        print s


class NotConvertedException(Exception):
    """Raised when a string (if required) has not been converted to Kappa compatible syntax"""

    def __init__(self):
        print "Must convert object due to identically named sites"

class SiteDef:
    """A site definition composed of a name and a finite set of states"""

    def __init__(self, n, ss):
        self.name = n
        self.state_list = ss

    def _write(self):
        if self.state_list:
            return "%s~%s"%(self.name, '~'.join(self.state_list))
        else:
            return self.name

    def write_as_bngl(self):
        return self._write()

    def write_as_kappa(self):
        return self._write()

    def __repr__(self):
        return "SiteDef(%s: %s)"%(self.name, ','.join(self.state_list))


class MoleculeDef:
    """A BNGL molecule type or Kappa agent signature"""

    def __init__(self, n, sds, snm, hss=False):
        """
        MoleculeDef initialization function that calculates inverse site name mapping.

        Parameters
        ----------
        n : str
            Alphanumeric string that identifies the molecule type
        st : list
            List of SiteDef instances
        snm : dict
            Dictionary that maps Kappa site names to BNGL site names
        hss : bool
            True if BNGL definition has identically named sites

        """
        self.name = n
        self.sites = sds
        self.site_name_map = snm
        self.has_site_symmetry = hss
        self._invert_site_name_map()

    def _invert_site_name_map(self):
        """Builds a dictionary of BNGL site names to a list of Kappa site names"""
        self.inv_site_name_map = {}
        for k in self.site_name_map.keys():
            v = self.site_name_map[k]
            if v not in self.inv_site_name_map.keys():
                self.inv_site_name_map[v] = [k]
            else:
                self.inv_site_name_map[v].append(k)

    def convert(self):
        """Converts MoleculeDef to use Kappa-compatible site names if necessary"""
        ss = []
        snm = {}
        k_track = {s: list(reversed(sorted(self.inv_site_name_map[s]))) for s in self.inv_site_name_map.keys()}
        for s in self.sites:
            name = k_track[s.name].pop()
            ss.append(SiteDef(name, s.state_list))
            snm[name] = s.name
        return MoleculeDef(self.name, ss, snm)

    def _write(self, is_bngl=True):
        """
        Writes MoleculeDef as string

        Parameters
        ----------
        is_bngl : bool
            If True, writes as BNGL string, else Kappa string

        Returns
        -------
        str
            MoleculeDef as BNGL molecule type or Kappa agent signature
        """
        if is_bngl or not self.has_site_symmetry:
            ss = [s.write_as_bngl() for s in self.sites]
        else:
            md = self.convert()
            ss = [s.write_as_kappa() for s in md.sites]
        return "%s(%s)" % (self.name, ','.join(ss))

    def write_as_bngl(self):
        """Writes MoleculeDef as BNGL string"""
        return self._write()

    def write_as_kappa(self):
        """Writes MoleculeDef as Kappa string"""
        return "%%agent: %s"%self._write(is_bngl=False)

    def __repr__(self):
        sites_string = ','.join(['%s->%s' % (k, v) for k, v in self.sites])
        return "MoleculeDef(name: %s, sites: %s)" % (self.name, sites_string)


class Molecule:
    """
    An individual molecule/agent inside a pattern.
        Note that it may not contain all sites listed in the associated
        MoleculeDef class
    """

    def __init__(self, name, sites):
        """
        Molecule initialization function. Sites are sorted by a predefined index

        Parameters
        ----------
        name : str
            The name of the molecule type
        sites : list
            A list of Sites that appear in the pattern

        """
        self.name = name
        self.sites = sorted(sites, key=lambda s: s.index)  # list of Sites

    def _node_name(self):
        """
        Provides a unique label for the Molecule based on its sites' states and bonds

        Returns
        -------
        str
            Unique string identifier based on bonding and sites' states

        """
        sstrs = []
        for s in self.sites:
            sstr = s._site_plus_state()
            if s.bond is not None:
                if s.bond.wild:
                    sstr += '?'
                elif s.bond.any:
                    sstr += '!'
                else:
                    sstr += 'b'
            else:
                sstr += 'u'
            sstrs.append(sstr)

        return self.name + ':' + '_'.join(sorted(sstrs))

    def has_identical_sites(self):
        """
        Checks to see if the Molecule has identically named sites

        Returns
        -------
        bool
            True if the molecule has sites with the same name, False otherwise
        """
        sn_set = set([s.name for s in self.sites])
        return len(sn_set) != len(self.sites)

    def _check_overlap(self, ms):
        pass

    def _enumerate_site(self, site_name, index, mdef, state=False):
        """
        Provides a list of sites in bound and unbound bond states, optionally
        enumerating all potential states as well

        Parameters
        ----------
        site_name : str
            Kappa-compatible name of site
        index : int
            Index of site in molecule
        mdef : MoleculeDef
        state : bool
            If True, enumerate site states as well as bond state, otherwise
            only bond states.

        Returns
        -------
        list
            List of Sites
        """
        ss = []
        if state:
            for sn, ss in mdef.sites:
                ss.append(Site(site_name))
        else:
            ss.append(Site(site_name, index, b=Bond(-1, w=True)))
            ss.append(Site(site_name, index, b=None))
            return ss

    # TODO
    # consider cases with potential for double counting events.  For example, if expanding a BNGL
    # molecule with identical sites, direct expansion (i.e. enumerating the renamed sites) may not
    # be appropriate.  There may be overlap in the rules, so we must enumerate binding or general
    # site state possibilities in this case
    def convert(self, mdef):
        """
        Converts a molecule that may have multiple identically named sites
        into a list of molecules compatible with Kappa semantics that
        enumerate all molecule configurations compatible with the original
        BNGL molecule configuration

        Parameters
        ----------
        mdef : MoleculeDef

        Returns
        -------
        list
            Returns a list of Molecules compatible with Kappa semantics,
            where each Molecule is a unique configuration to accommodate
            any identically named sites.
        """
        un_site_names = set([s.name for s in self.sites])
        un_configs_per_site = {s: {} for s in un_site_names}
        for i in range(len(self.sites)):
            s = self.sites[i]
            if s not in un_configs_per_site[s.name]:
                un_configs_per_site[s.name][s] = 0
            un_configs_per_site[s.name][s] += 1

        def rename_sites(names, site):
            return tuple([Site(name, site.index, s=site.state, b=site.bond) for name in names])

        # Check for the possibility of overlapping patterns
        possible_overlap = {k: False for k in un_configs_per_site.keys()}
        for k in un_configs_per_site.keys():
            num_identical_sites = len(mdef.inv_site_name_map[k])
            if num_identical_sites > 1 and k in un_configs_per_site.keys():
                num_present_sites = sum(un_configs_per_site[k].values())
                if num_identical_sites > num_present_sites > 1:
                    possible_overlap[k] = True
                    break

        k_configs = {}
        for sn in un_configs_per_site.keys():
            k_configs[sn] = []
            k_sn_names = set(mdef.inv_site_name_map[sn])
            cur_combs = []
            for s, n in un_configs_per_site[sn].iteritems():
                if len(cur_combs) == 0:
                    cur_combs = [rename_sites(names, s) for names in it.combinations(k_sn_names, n)]
                else:
                    tmp_combs = []
                    for cc in cur_combs:
                        rem_names = k_sn_names - set(map(lambda l: l.name, cc))
                        new_combs = [rename_sites(names, s) for names in it.combinations(rem_names, n)]
                        for nc in new_combs:
                            tmp_combs.append(cc + nc)
                    cur_combs = tmp_combs
            k_configs[sn] = cur_combs
            if possible_overlap[sn]:

                # generate all possible configurations of site sn
                # for remaining (unseen) sites, enumerate configurations

                for idx in range(len(cur_combs),len(mdef.inv_site_name_map[k])):
                    new_combs = gen_configs(idx)


                bound_configs = [Site(sn, i, state, 'w')]
                unbound_configs = []

        k_prod = list(it.product(*k_configs.values()))

        return [Molecule(self.name, [e for t in tt for e in t]) for tt in k_prod]

    def _write(self, bngl=True):
        """
        Writes the Molecule in Kappa or BNGL syntax

        Parameters
        ----------
        bngl : bool
            Writes the string as BNGL compatible if True, otherwise writes as Kappa

        Returns
        -------
        str
            String representation of Molecule
        """
        ss = []
        for s in self.sites:
            if bngl:
                ss.append(s.write_as_bngl())
            else:
                ss.append(s.write_as_kappa())
        return '%s(%s)' % (self.name, ','.join(ss))

    def write_as_bngl(self):
        """Writes Molecule as BNGL string"""
        return self._write()

    # returns list of molecule strings
    # check to see how many of possible symm sites are present in pattern
    def write_as_kappa(self):
        """Writes Molecule as Kappa string and checks for conversion if necessary"""
        if len(set([s.name for s in self.sites])) < len(self.sites):
            raise NotConvertedException
        return self._write(False)

    def __repr__(self):
        return 'Molecule(name: %s, sites: %s)' % (self.name, ', '.join([str(x) for x in self.sites]))


class Site:
    """
    A subunit of a Molecule that can engage in binding and exist
    various states
    """

    def __init__(self, n, i, s=None, b=None):
        """
        Site initialization function

        Parameters
        ----------
        n : str
            Site name
        i : int
            The index of this site in the Molecule
        s : str
            The site's state (could be None)
        b : Bond
            The site's binding state
        """
        self.name = n
        self.index = i
        self.state = s
        self.bond = b

    def _site_plus_state(self):
        """Builds Kappa/BNGL-compatible string composed of site's name state"""
        state = '' if self.state is None else '~%s' % self.state
        return self.name + state

    def _write(self, kappa=False):
        """
        Builds site string

        Parameters
        ----------
        kappa : bool
            If True, then the Site is written in Kappa syntax, BNGL if False

        Returns
        -------
        str
            Site written as BNGL or Kappa string
        """
        s = self._site_plus_state()
        if self.bond is not None and kappa:
            s += self.bond.write_as_kappa()
        elif self.bond is not None:
            s += self.bond.write_as_bngl()
        return s

    def write_as_bngl(self):
        """Write Site as BNGL string"""
        return self._write()

    def write_as_kappa(self):
        """Write Site as Kappa string"""
        return self._write(True)

    def __eq__(self, other):
        """
        Check for equality with another site

        Parameters
        ----------
        other : Site
            Another Site for equality checking

        Returns
        -------
        bool
            True if Sites are equivalent in terms of name, state, and bond
        """
        if isinstance(other, self.__class__):
            return self.name == other.name and self.state == other.state and self.bond == other.bond
        return False

    def __ne__(self, other):
        """Check for inequality with another Site"""
        return not self == other

    def __hash__(self):
        return hash((self.name, self.state, self.bond))

    def __repr__(self):
        return 'Site(name: %s, state: %s, bond: %s)' % (self.name, self.state, self.bond)


class Bond:
    """Represents a bond of any type in a Site"""

    def __init__(self, n, w=False, a=False):
        """
        Bond initialization function

        Parameters
        ----------
        n : int
            Bond number.  Negative numbers indicate the absence of a specific bond.
            Positive numbers override the other arguments.
        w : bool
            True denotes a 'wild' bond.
        a : bool
            True denotes 'any' bond, including the absence of a bond
        """
        self.wild = w
        self.num = int(n)  # negative numbers indicate absence of specific bond, positive numbers will override w and a
        self.any = a
        if self.num < 0:
            self.num = -1
            assert (self.wild or self.any)
        assert (not (self.wild and self.any))

    def write_as_bngl(self):
        """Write bond as BNGL string"""
        s = ''
        if self.num >= 0:
            s = '!%s' % self.num
        if self.wild:
            s = '!+'
        elif self.any:
            s = '!?'
        return s

    def write_as_kappa(self):
        """Write bond as Kappa string"""
        s = ''
        if self.wild:
            s = '!_'
        elif self.any:
            s = '?'
        else:
            s = '!%s' % self.num
        return s

    def __eq__(self, other):
        """
        Check for equality with another Bond

        Parameters
        ----------
        other : Bond

        Returns
        -------
        bool
            True if both bond numbers are positive or if both bonds are
            wild or any.
        """
        if isinstance(other, self.__class__):
            num_check = False
            if (self.num < 0 and other.num < 0) or (self.num == other.num):
                num_check = True
            return (self.wild == other.wild and self.any == other.any and num_check)
        return False

    def __ne__(self, other):
        """Check for inequality with another Bond"""
        return not self == other

    def __hash__(self):
        return hash((self.num, self.wild, self.any))

    def __repr__(self):
        b_string = str(self.num)
        if self.wild:
            b_string = 'wild'
        elif self.any:
            b_string = 'any'
        return b_string


# CPattern is defined as a pattern (list of molecules) that are attached
class CPattern:
    """
    Defines a pattern as a list of Molecules

        In BNGL the pattern is defined as molecules joined by the '.' operator.  This
        means that a pattern matches either a single molecule or a complex of multiple
        molecules.  This is in contrast to Kappa, where all patterns are joined by the
        ',' operator and may match multiple molecules or complexes.
    """

    def __init__(self, ml):
        """
        CPattern initialization function

        Parameters
        ----------
        ml : list
            List of Molecules that are a part of the pattern
        """
        self.molecule_list = ml

    def _build_graph(self):
        """
        Builds a graph representation of the CPattern

        Returns
        -------
        Graph
        """
        g = nx.Graph()

        # Put molecules and sites into graph as nodes where nodes are named
        # with the molecule name and a sorted list of sites (and states if
        # present) without bonds.  Edges are named with sorted site names
        for i, mol in enumerate(self.molecule_list):
            nstr = mol._node_name()
            g.add_node(i, name=nstr)

        # enter bonds into graph
        for i in range(len(self.molecule_list) - 1):
            for j in range(i + 1, len(self.molecule_list)):
                for s in self.molecule_list[i].sites:
                    if s.bond is None:
                        break
                    for s2 in self.molecule_list[j].sites:
                        if s2.bond is None:
                            break
                        if s.bond == s2.bond and s.bond.num >= 0:
                            bond_name = '-'.join(sorted([s.name, s2.name]))
                            g.add_edge(i, j, name=bond_name)

        return g

    # returns a list of graphs with all permutations of identical nodes
    @staticmethod
    def _permute(g):
        """
        Finds all the permutations of a graph's identical nodes

        Parameters
        ----------
        g : Graph
            An input graph to permute
        Returns
        -------
        list
            List of Graph instances with all permutations of identically named nodes
        """
        # dict of all (node name, node list) pairs
        node_name_dict = {}
        for node in g.nodes():
            name = g.node[node]['name']
            if name not in node_name_dict.keys():
                node_name_dict[name] = [node]
            else:
                node_name_dict[name].append(node)

        # Remove nodes with unique names
        for n in node_name_dict.keys():
            if len(node_name_dict[n]) == 1:
                node_name_dict.pop(n)
            else:
                node_name_dict[n] = sorted(node_name_dict[n])

        # Find all permutations of all node types (possible node mappings)
        per_node_tuples = []
        for n in node_name_dict.keys():
            op = node_name_dict[n]  # Original permutation of node type n
            perms = list(it.permutations(op))  # All permutations of node type n

            # Convert each permutation into a list of tuples.  The list is a
            # mapping from the original ordering to a permutation.  Tuples is a
            # list of these lists of tuples for each permutation of node type n
            tuples = [list(it.izip(op, p)) for p in perms]
            per_node_tuples.append(tuples)

        # First take the Cartesian product over the list of lists of
        # permutations for each node type producing all combinations of
        # permutations over all node types.  Each result is a tuple of node
        # type permutations (a tuple of lists of tuples), so flatten the
        # highest tuple to a single list of tuples.  Convert it to a
        # dictionary for the actual node relabeling
        relabeling_dicts = [dict(it.chain(*x)) for x in it.product(*per_node_tuples)]

        return [nx.relabel_nodes(g, rd) for rd in relabeling_dicts]

    # since conversion to Kappa removes identical site names, automorphisms only exist
    # on the level of the molecule
    def automorphisms(self):
        """
        Determines the number of automorphisms in a CPattern

        Notes
        -----
        Only performed after conversion to Kappa compatible site names.  Thus
        the count will only include automorphisms on the level of molecules

        Returns
        -------
        int
            The number of automorphisms in this CPattern
        """
        # Check to make sure conversion to Kappa compatible site names has occurred
        for m in self.molecule_list:
            if m.has_identical_sites():
                raise NotConvertedException
        # If all molecules are unique, exit with count 1.  Otherwise calculate the
        # number of automorphisms
        if len([m._node_name() for m in self.molecule_list]) == len(set([m.name for m in self.molecule_list])):
            return 1
        else:
            g = self._build_graph()
            am = 0
            # Matching functions that take the node name string (node type) into account
            nm = iso.categorical_node_match('name', '')
            em = iso.categorical_edge_match('name', '')
            for gp in self._permute(g):
                is_iso = iso.is_isomorphic(g, gp, edge_match=em, node_match=nm)
                equal_edges = set(g.edges()) == set(gp.edges())
                # Automorphisms require that a permutation of a graph's nodes
                # is both isomorphic to the original graph and that the
                # permutation's edges are preserved from the original graph
                if is_iso and equal_edges:
                    am += 1
            return am

    # returns a list of objects with renamed sites
    def convert(self, mdefs):
        """
        Converts a CPattern to a list of Kappa compatible CPatterns

        Parameters
        ----------
        mdefs : list
            List of MoleculeDefs corresponding to Molecules in the CPattern

        Returns
        -------
        list
            List of unique Kappa compatible CPatterns that correspond to the
            original CPattern
        """
        k_str_mol_list = []
        for m in self.molecule_list:
            for md in mdefs:
                if m.name == md.name:
                    k_str_mol_list.append(m.convert(md))
        k_patterns = it.product(*k_str_mol_list)
        # Remove doubles and preserve molecule order
        seen = set()
        un_k_patterns = []
        for pat in k_patterns:
            s_pat = tuple(sorted(pat))
            if s_pat not in seen:
                seen.add(s_pat)
                un_k_patterns.append(CPattern(pat))
        return un_k_patterns

    def _write(self, bngl=True):
        """
        Writes the CPattern in either BNGL or Kappa syntax

        Parameters
        ----------
        bngl : bool
            If True, write as BNGL string, otherwise Kappa

        Returns
        -------
        str
        """
        cps = []
        for m in self.molecule_list:
            if bngl:
                cps.append(m.write_as_bngl())
            else:
                cps.append(m.write_as_kappa())
        joiner = '.' if bngl else ','
        return joiner.join(cps)

    def write_as_bngl(self):
        """Write the CPattern as a BNGL string"""
        return self._write()

    def write_as_kappa(self):
        """Write the CPattern as a Kappa string"""
        return self._write(False)

    def __repr__(self):
        return '\n'.join([str(x) for x in self.molecule_list])


class InitialCondition:
    """Initial conditions for seeding simulation"""

    def __init__(self, s, a, ain=True):
        """
        Initial condition initialization function

        Parmameters
        -----------
        s : CPattern or Molecule
        a : Expression or float
        ain : bool
            True if a is a number, False if not
        """
        self.species = s
        self.amount = a
        self.amount_is_number = ain

    def convert(self, mdefs):
        """
        Converts species to Kappa compatible species

        Parameters
        ----------
        mdefs : MoleculeDef

        Returns
        -------
        list
            List of InitialCondition instances
        """
        ss = self.species.convert(mdefs)
        lss = len(ss)
        if lss == 1:
            amount = self.amount
        else:
            amount = self.amount / float(lss) if self.amount_is_number else Expression(
                ['('] + self.amount.expr + [')', '/', lss])
        return [InitialCondition(s, amount, self.amount_is_number) for s in ss]

    def write_as_bngl(self):
        """Write as BNGL string"""
        amount = self.amount if self.amount_is_number else self.amount.write_as_bngl()
        return '%s %s' % (self.species.write_as_bngl(), amount)

    # if there are symmetries, the initial condition amount is divided evenly among species
    def write_as_kappa(self):
        """
        Write as Kappa string

        If there are symmetries resulting from identically named sites then
        the initial condition quantity is divided evenly among the species
        """
        amount = self.amount if self.amount_is_number else self.amount.write_as_kappa()
        return '%%init: %s %s' % (amount, self.species.write_as_kappa())

    def __repr__(self):
        return "Init(species: %s, quantity: %s)" % (self.species, self.amount)


class Parameter:
    """Defines a constant value"""

    def __init__(self, n, v):
        """
        Parameter initialization function

        Parameters
        ----------
        n : str
            Parameter name
        v : float or int
            Parameter value
        """
        self.name = n
        self.value = v

    def write_as_bngl(self):
        """Writes Parameter as BNGL string"""
        return '%s %s' % (self.name, self.value)

    def write_as_kappa(self):
        """Writes Parameter as Kappa string"""
        return '%%var: \'%s\' %s' % (self.name, self.value)

    def __repr__(self):
        return "Parameter(name: %s, value: %s)" % (self.name, self.value)


# special parsing required for 'if', 'log' functions
# can implement conversion of certain types of values (e.g. log10(x) to log(x)/log(10))
class Expression:
    """Defines an algebraic expression"""

    def __init__(self, atom_list):
        """
        Expression initialization function

        Parameters
        ----------
        atom_list : list
            List of tokens from parse_math_expr function.  Ordered as operators,
            values, variables
        """
        self.atom_list = atom_list  # list from parse_math_expr listing (in order) operators, values, variables

    def write_as_bngl(self):
        """Writes Expression as BNGL string"""
        return ''.join(self.atom_list)

    def write_as_kappa(self):
        """Writes Expression as Kappa string"""
        expr = ''

        i = 0
        while (i < len(self.atom_list)):
            a = self.atom_list[i]
            if a in bngl_to_kappa_func_map.keys():
                trig_func_match = re.compile('sinh|cosh|tanh|asinh|acosh|atanh')
                if re.match('log', a) or re.match(trig_func_match, a):
                    expr += bngl_to_kappa_func_map[a](self.atom_list[i + 2])
                    i += 4
                else:
                    expr += bngl_to_kappa_func_map[a]
            elif re.match('[A-Za-z]', a):
                expr += '\'%s\'' % a
            else:
                expr += a
            i += 1

        return expr

    def __repr__(self):
        return "Expression(expr: %s)" % self.write_as_bngl()


class Function:
    """Defines dynamic quantities (i.e. those that depend on Observables)"""

    def __init__(self, name, expr):
        """
        Function initialization function

        Parameters
        ----------
        name : str
            Function name
        expr : Expression
            Function expression
        """
        self.name = name
        self.expr = expr

    def write_as_bngl(self):
        """Writes function as BNGL string"""
        return '%s()=%s' % (self.name, self.expr.write_as_bngl())

    def write_as_kappa(self, as_obs=True):
        """
        Writes function as Kappa string

        Parameters
        ----------
        as_obs : bool
            If True, writes functions as observables, otherwise as variables
        """
        dec = 'obs' if as_obs else 'var'
        return "%%%s: '%s' %s" % (dec, self.name, self.expr.write_as_kappa())

    def __repr__(self):
        return "Function(name: %s, expr: %s" % (self.name, self.expr)


class Rate:
    """Defines a Rule's Rate"""

    def __init__(self, r, intra=False):  # can be string (including a single number) or Expression
        """
        Rate initialization function

        Parameters
        ----------
        r : str or float or Expression
            Number, string or algebraic expression defining a rate law
        intra : bool
            True if the rate governs an intramolecular reaction, False otherwise
        """
        self.rate = r
        self.intra_binding = intra

    def write_as_bngl(self):
        """Write Rate as BNGL string"""
        try:
            return self.rate.write_as_bngl()
        except AttributeError:
            return str(self.rate)

    def write_as_kappa(self):
        """Write Rate as Kappa string"""
        try:
            rate_string = self.rate.write_as_kappa()
        except AttributeError:
            rate_string = str(self.rate) if _is_number(self.rate) else "'%s'" % self.rate
        return rate_string if not self.intra_binding else '0 {%s}' % rate_string

    def __repr__(self):
        return "Rate: %s" % self.rate


# TODO implement check for rate as raw number before writing
# TODO add check for automorphisms in CPatterns
class Rule:
    """Defines a rule"""

    # lhs, rhs are lists of CPatterns, rate/rev_rate are Rates, rev is bool (true for reversible rules),
    # amb_mol is boolean denoting a rule that (in Kappa) has ambiguous molecularity
    def __init__(self, lhs, rhs, rate, rev=False, rev_rate=None):
        """
        Rule initialization function

        Parameters
        ----------
        lhs : list
            List of CPatterns
        rhs : list
            List of CPatterns
        rate : Rate
            Rate for the lhs -> rhs reaction
        rev : bool
            True if the rule is reversible, False otherwise (default)
        rev_rate : Rate
            Rate for the reverse rhs -> lhs reaction if present
        """
        self.lhs = lhs
        self.rhs = rhs
        self.rate = rate
        self.rev = rev
        self.arrow = '->' if not rev else '<->'
        self.rev_rate = None if not rev else rev_rate  # rev overrides rev_rate

    def write_as_bngl(self):
        """Writes the rule as a BNGL string"""

        lhs_string = '+'.join([p.write_as_bngl() for p in self.lhs])
        rhs_string = '+'.join([p.write_as_bngl() for p in self.rhs])
        if self.rev:
            rate_string = self.rate.write_as_bngl() + ',' + self.rev_rate.write_as_bngl()
        else:
            rate_string = self.rate.write_as_bngl()
        return '%s %s %s %s' % (lhs_string, self.arrow, rhs_string, rate_string)

    def write_as_kappa(self, mdefs):
        """Writes the rule as a Kappa string"""

        # possible actions  (root object is list (rule.lhs and rule.rhs))
        #  - iterable_item_added/removed (binding, unbinding)
        #  - type_changes (binding, unbinding)
        #  - value_changes (state change)

        def kappa_rule_string(lhss, ars, rhss, rs):
            return '%s %s %s @ %s' % (lhss, ars, rhss, rs)

        lhs_strings = it.product(*[p.write_as_kappa(mdefs) for p in self.lhs])

        pass

    # lhs_string = ','.join([p.write_as_kappa() for p in self.lhs])
    # rhs_string = ','.join([p.write_as_kappa() for p in self.rhs])
    # if self.rev:
    # 	rate_string = self.rate.write_as_kappa() + ',' + self.rev_rate.write_as_kappa()
    # else:
    # 	rate_string = self.rate.write_as_kappa()
    # return '%s %s %s @ %s'%(lhs_string,self.arrow,rhs_string,rate_string)

    def __repr__(self):
        if not self.rev:
            return "Rule(lhs: %s, rhs: %s, rate: %s)" % (self.lhs, self.rhs, self.rate)
        else:
            return "Rule(lhs: %s, rhs: %s, rate: %s, rev_rate: %s)" % (self.lhs, self.rhs, self.rate, self.rev_rate)


class Observable:
    """Defines observables"""

    def __init__(self, n, ps, t='m'):
        """
        Observable initialization functions

        Parameters
        ----------
        n : str
            Observable name
        ps : list
            List of CPatterns
        t : str
            Must be 'm' or 's'. Only 'm' is compatible with Kappa
        """
        self.name = n
        if re.match('[sS]$', t):
            self.type = 'Species'
        elif re.match('[mM]$', t):
            self.type = 'Molecules'
        else:
            raise Exception("not a valid observable type: %s" % t)
        self.cpatterns = ps  # a list of CPatterns

    def write_as_bngl(self):
        """Writes Observable as BNGL string"""
        return "%s %s %s" % (self.type, self.name, ' '.join([p.write_as_bngl() for p in self.cpatterns]))

    def write_as_kappa(self, mdefs):
        """Writes Observable as Kappa string"""
        if self.type == 'Species':
            print "Kappa does not have a Species-like observable; printing '%s' as Molecules-like observable" % self.name

        obs_strs = []
        for p in self.cpatterns:
            # sorted for determinism (testing)
            kos = '+'.join(sorted(['|%s|' % x.write_as_kappa() for x in p.convert(mdefs)]))
            obs_strs.append(kos)

        obs = '+'.join(obs_strs)
        return '%%obs: \'%s\' %s' % (self.name, obs)

    def __repr__(self):
        return "Obs(name: %s, pattern: %s)" % (self.name, ' '.join(self.cpatterns))


class Model:
    """Object for rule-based model"""

    def __init__(self):
        """Model initialization function"""
        self.molecules = []
        self.initial_cond = []
        self.observables = []
        self.functions = []
        self.rules = []
        self.parameters = []

    def write_as_bngl(self, file_name):
        """Writes Model as BNGL file"""
        s = 'begin model\n\nbegin parameters\n\n'
        for p in self.parameters:
            s += '\t%s\n' % p.write_as_bngl()
        s += '\nend parameters\n\n'
        s += 'begin molecule types\n\n'
        for m in self.molecules:
            s += '\t%s\n' % m.write_as_bngl()
        s += '\nend molecule types\n\n'
        s += 'begin initial conditions\n\n'
        for i in self.initial_cond:
            s += '\t%s\n' % i.write_as_bngl()
        s += '\nend initial conditions\n\n'
        s += 'begin observables\n\n'
        for o in self.observables:
            s += '\t%s\n' % o.write_as_bngl()
        s += '\nend observables\n\n'
        s += 'begin functions\n\n'
        for f in self.functions:
            s += '\t%s\n' % f.write_as_bngl()
        s += '\nend functions\n\n'
        s += 'begin reaction rules'
        for r in self.rules:
            s += '\t%s\n' % r.write_as_bngl()
        s += 'end reaction rules\n\n'
        s += 'end model\n'

        f = open('%s' % file_name)
        f.write(s)
        f.close()

    # check for rules with molecular ambiguity
    def write_as_kappa(self, func_as_obs=False):
        """
        Writes Model as Kappa file

        Parameters
        ----------
        func_as_obs : bool
            If True, writes functions as observables, otherwise as variables
        """
        s = ''
        for m in self.molecules:
            s += '%s\n' % m.write_as_kappa()
        s += '\n'
        for p in self.parameters:
            s += '%s\n' % p.write_as_kappa()
        s += '\n'
        for o in self.observables:
            s += '%s\n' % o.write_as_kappa()
        s += '\n'
        for f in self.functions:
            s += '%s\n' % f.write_as_kappa(func_as_obs)  # defaults to printing all "functions" as observables
        s += '\n'
        for i in self.initial_cond:
            s += '%s\n' % i.write_as_kappa()
        s += '\n'
        for r in self.rules:
            s += '%s\n' % r.write_as_kappa()
        s += '\n'

    def add_molecule(self, mol):
        """
        Adds a Molecule to Model

        Parameters
        ----------
        mol : Molecule
        """
        self.molecules.append(mol)

    def add_init(self, init):
        """
        Adds an InitialCondition to Model

        Parameters
        ----------
        init : InitialCondition
        """
        self.initial_cond.append(init)

    def add_obs(self, obs):
        """
        Adds an Observable to Model

        Parameters
        ----------
        obs : Observable
        """
        self.observables.append(obs)

    def add_func(self, func):
        """
        Adds a Function to Model

        Parameters
        ----------
        func : Function
        """
        self.functions.append(func)

    def add_rule(self, rule):
        """
        Adds a Rule to Model

        Parameters
        ----------
        rule : Rule
        """
        self.rules.append(rule)

    def add_parameter(self, param):
        """
        Adds a Parameter to Model

        Parameters
        ----------
        param : Parameter
        """
        self.parameters.append(param)


class Reader:
    """Inherited class for reading rule-based modeling files."""

    def __init__(self, file_name):
        """
        Reader initialization function

        Parameters
        ----------
        file_name : str
            Rule-based model file
        """
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
    """Reader for BNGL model files"""

    def __init__(self, file_name):
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
        model = Model(self.file_name)
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
            raise NotAMoleculeException(smstr)
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
        amount_is_number = _is_number(amount)
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
                except NotAMoleculeException:
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
            raise NotCompatibleException("Kappa functions cannot accommodate local functions:\n\t%s\n" % sline)
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

        point = Literal(".")
        e = CaselessLiteral("E")
        fnumber = Combine(Word("+-" + nums, nums) +
                          Optional(point + Optional(Word(nums))) +
                          Optional(e + Word("+-" + nums, nums)))
        ident = Word(alphas, alphas + nums + "_$")

        plus = Literal("+")
        minus = Literal("-")
        mult = Literal("*")
        div = Literal("/")
        lpar = Literal("(")
        rpar = Literal(")")
        addop = plus | minus
        multop = mult | div
        expop = Literal("^")
        pi = CaselessLiteral("PI")

        expr = Forward()
        atom = (Optional("-") + (pi | e | fnumber | ident + lpar + expr + rpar | ident) | (lpar + expr + rpar))

        # by defining exponentiation as "atom [ ^ factor ]..." instead of "atom [ ^ atom ]...", we get right-to-left exponents, instead of left-to-righ
        # that is, 2^3^2 = 2^(3^2), not (2^3)^2.
        factor = Forward()
        factor << atom + ZeroOrMore((expop + factor))

        term = factor + ZeroOrMore((multop + factor))
        expr << term + ZeroOrMore((addop + term))
        pattern = expr

        return pattern.parseString(estr.strip())


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
    'log2': lambda s: '([log](%s)/[log](2))' % s,
    'log10': lambda s: '([log](%s)/[log](10))' % s,
    'ln': '[log]',
    'exp': '[exp]',
    'sinh': lambda s: '([exp](%s) - [exp](-%s))/2.0' % (s, s),
    'cosh': lambda s: '([exp](%s) + [exp](-%s))/2.0' % (s, s),
    'tanh': lambda s: '([exp](%s) - [exp](-%s))/([exp](%s) + [exp](-%s))' % (s, s, s, s),
    'asinh': lambda s: '[log](%s + [sqrt](%s^2 + 1))' % (s, s),
    'acosh': lambda s: '[log](%s + [sqrt](%s^2 - 1))' % (s, s),
    'atanh': lambda s: '0.5*[log]((1+%s)/(1-%s))' % (s, s)
}

# TO CHECK:
# calculate arctrig or trigh
# rint is round
bngl_other_builtin_funcs = set(['asin', 'acos', 'atan', 'rint', 'abs', 'if', 'min', 'max', 'sum', 'avg'])
bngl_binary_operators = set(['==', '!=', '>=', '<='])

# not sure how boolean expressions work in Kappa
# overlap_binary_operators = set(['*','/','^','+','-','&&','||','<','>'])

kappa_other_builtin_funcs = set(['[Tsim]', '[Tmax]', '[E]', '[E-]', '[Emax]', '[pp]', 'inf'])
# [int] is floor
kappa_other_operators = set(['[int]', '[mod]', '[max]', '[min]', '[not]', '[true]', '[false]'])
