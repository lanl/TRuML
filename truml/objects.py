import re
import itertools as it
import logging
import networkx as nx
import networkx.algorithms.isomorphism as iso
import rbexceptions


class SiteDef:
    """A site definition composed of a name and a finite set of states"""

    def __init__(self, n, ss=[]):
        self.name = n
        self.state_list = ss

    def _write(self):
        if self.state_list:
            return "%s~%s" % (self.name, '~'.join(self.state_list))
        else:
            return self.name

    def write_as_bngl(self):
        return self._write()

    def write_as_kappa(self):
        return self._write()

    def __repr__(self):
        if not self.state_list:
            return "SiteDef(%s)" % self.name
        else:
            return "SiteDef(%s: %s)" % (self.name, ','.join(self.state_list))


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
        return "%%agent: %s" % self._write(is_bngl=False)

    def __repr__(self):
        return "MoleculeDef(name: %s, sites: %s)" % (self.name, self.sites)


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

    def _enumerate_site(self, site_name, index, mdef, need_state=False):
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
        need_state : bool
            If True, enumerate site states as well as bond state, otherwise
            only bond states.

        Returns
        -------
        list
            List of Sites
        """
        ss = []
        if need_state:
            for s in mdef.sites:
                if site_name == s.name:
                    for state in s.state_list:
                        ss.append(Site(site_name, index, state, b=Bond(-1, w=True)))
                        ss.append(Site(site_name, index, state, b=None))
                    break
        else:
            ss.append(Site(site_name, index, b=Bond(-1, w=True)))
            ss.append(Site(site_name, index, b=None))
        return ss

    @staticmethod
    def _site_state_present(ss):
        """
        Checks to see if a Site in a list of Site instances has a defined state

        Parameters
        ----------
        ss : list
            List of Site instances

        Returns
        -------
        bool
            True if at least one site has a defined state, False otherwise

        """
        for s in ss:
            if s.state:
                return True
        return False

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

        def rename_site(name, site):
            return Site(name, site.index, s=site.state, b=site.bond)

        def rename_sites(names, site):
            return tuple([rename_site(name, site) for name in names])

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

            if possible_overlap[sn]:
                need_state = self._site_state_present(un_configs_per_site[sn])
                indices = range(sum(un_configs_per_site[sn].values()), len(mdef.inv_site_name_map[k]))

                for idx in indices:
                    possible_sites = self._enumerate_site(sn, idx, mdef, need_state)
                    tmp_combs = []
                    for cc in cur_combs:
                        rem_names = k_sn_names - set(map(lambda l: l.name, cc))
                        new_combs = [rename_site(x, y) for x, y in it.product(rem_names, possible_sites)]
                        for nc in new_combs:
                            tmp_combs.append(cc + (nc,))
                    cur_combs = tmp_combs

            k_configs[sn] = cur_combs

        k_prod = it.product(*k_configs.values())

        return sorted([Molecule(self.name, [e for t in tt for e in sorted(t)]) for tt in k_prod])

    def bound_to(self, other):
        if not isinstance(other, self.__class__):
            return False
        else:
            for s in self.sites:
                for t in other.sites:
                    if s.bond is None or t.bond is None:
                        continue
                    if s.bond.num == t.bond.num and s.bond.num > -1:
                        return True
            return False

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
            raise rbexceptions.NotConvertedException
        return self._write(False)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.name == other.name and frozenset(self.sites) == frozenset(other.sites)
        else:
            return False

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        return self.write_as_bngl() < other.write_as_bngl()

    def __hash__(self):
        return hash((self.name, frozenset(self.sites)))

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

    def __lt__(self, other):
        return self.write_as_bngl() < other.write_as_bngl()

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
            True if both bond numbers are equal or if both bonds are
            wild or any.
        """
        if isinstance(other, self.__class__):
            num_check = False
            if (self.num < 0 and other.num < 0) or (self.num == other.num):
                num_check = True
            return self.wild == other.wild and self.any == other.any and num_check
        return False

    def __ne__(self, other):
        """Check for inequality with another Bond"""
        return not self == other

    def __hash__(self):
        return hash((-1 if self.num < 0 else self.num, self.wild, self.any))

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

    def num_molecules(self):
        """Determines the number of molecules in the pattern"""
        return len(self.molecule_list)

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
                raise rbexceptions.NotConvertedException
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
            List of MoleculeDef instances corresponding to Molecule instances in the CPattern

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
        return [CPattern(pat) for pat in k_patterns]
        # # Remove doubles and preserve molecule order
        # seen = set()
        # un_k_patterns = []
        # for pat in k_patterns:
        #     if pat not in seen:
        #         seen.add(pat)
        #         un_k_patterns.append(CPattern(pat))
        # return un_k_patterns

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

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return frozenset(self.molecule_list) == frozenset(other.molecule_list)
        else:
            return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(frozenset(self.molecule_list))

    def __repr__(self):
        return "CPattern(" + '--'.join([str(x) for x in self.molecule_list]) + ")"


# TODO prevent dynamic quantities from being used as the initial amount
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

    def write_as_bngl(self, namespace):
        """Write as BNGL string"""
        amount = self.amount if self.amount_is_number else self.amount.write_as_bngl(namespace)
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
        Parameter initialization function.  The value can be an integer or a constant
        algebraic expression (cannot contain dynamic quantities such as observables)

        Parameters
        ----------
        n : str
            Parameter name
        v : float or int or Expression
            Parameter value
        """
        self.name = n
        self.value = v

    def write_as_bngl(self, namespace):
        """Writes Parameter as BNGL string"""
        bname = namespace[self.name]
        val = self.value.write_as_bngl(namespace) if isinstance(self.value, Expression) else self.value
        return '%s %s' % (bname, val)

    def write_as_kappa(self):
        """Writes Parameter as Kappa string"""
        val = self.value.write_as_kappa() if isinstance(self.value, Expression) else self.value
        return '%%var: \'%s\' %s' % (self.name, val)

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
            List of tokens from parse_alg_expr functions.
        """
        self.atom_list = atom_list

    def write_as_bngl(self, namespace):
        """Writes Expression as BNGL string"""
        conv_atom_list = []
        for atom in self.atom_list:
            if atom in namespace.keys():
                conv_atom_list.append(namespace[atom])
            else:
                conv_atom_list.append(atom)
        return ''.join(conv_atom_list)

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
        return "Expression(expr: %s)" % ''.join(self.atom_list)


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

    def write_as_bngl(self, namespace):
        """Writes function as BNGL string"""
        bname = namespace[self.name]
        return '%s()=%s' % (bname, self.expr.write_as_bngl(namespace))

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

    def write_as_bngl(self, namespace):
        """Write Rate as BNGL string"""
        try:
            return self.rate.write_as_bngl(namespace)
        except AttributeError:
            if isinstance(self.rate, str):
                return namespace[self.rate]
            else:
                return str(self.rate)

    def write_as_kappa(self):
        """Write Rate as Kappa string"""
        try:
            rate_string = self.rate.write_as_kappa()
        except AttributeError:
            rate_string = str(self.rate) if is_number(self.rate) else "'%s'" % self.rate
        return rate_string if not self.intra_binding else '0 {%s}' % rate_string

    def __repr__(self):
        return "Rate: %s" % self.rate


class Rule:
    """Defines a rule"""

    # lhs, rhs are lists of CPatterns, rate/rev_rate are Rates, rev is bool (true for reversible rules),
    # amb_mol is boolean denoting a rule that (in Kappa) has ambiguous molecularity
    def __init__(self, lhs, rhs, rate, rev=False, rev_rate=None, label=None, delmol=False):
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
        label : str
            Label identifying the rule
        delmol : bool
            If the rule governs degradation, this allows fragments to persist if unseen
            bonds can be broken by the rule's application
        """
        self.lhs = lhs
        self.rhs = rhs
        self.rate = rate
        self.rev = rev
        self.arrow = '->' if not rev else '<->'
        self.rev_rate = None if not rev else rev_rate  # rev overrides rev_rate
        self.label = label
        self.delmol = delmol

    def convert(self, lhs_mdefs, rhs_mdefs):
        """
        Converts a Rule to a Kappa-compatible naming scheme

        Parameters
        ----------
        lhs_mdefs : list
            List of MoleculeDef instances corresponding to CPatterns on the
            Rule's left-hand side
        rhs_mdefs : list
            List of MoleculeDef instances corresponding to CPatterns on the
            Rule's right-hand side
        Returns
        -------
        list
            List of Rule instances
        """
        k_lhs, k_rhs = [], []
        for cp in self.lhs:
            k_lhs.append(cp.convert(lhs_mdefs))
        for cp in self.rhs:
            k_rhs.append(cp.convert(rhs_mdefs))
        all_lhs = it.product(*k_lhs)
        all_rhs = it.product(*k_rhs)

        z = zip(all_lhs, all_rhs)  # order in lhs and rhs conversions are preserved
        rs = []
        for l, r in z:
            rs.append(Rule(l, r, self.rate, self.rev, self.rev_rate))
        return list(set(rs))

    def write_as_bngl(self, namespace, dot=False):
        """Writes the rule as a BNGL string"""
        if not self.lhs:
            lhs_string = '0'
        elif dot:
            lhs_string = '.'.join([p.write_as_bngl() for p in self.lhs])
        else:
            lhs_string = '+'.join([p.write_as_bngl() for p in self.lhs])
        if not self.rhs:
            rhs_string = '0'
        elif dot:
            rhs_string = '.'.join([p.write_as_bngl() for p in self.rhs])
        else:
            rhs_string = '+'.join([p.write_as_bngl() for p in self.rhs])
        if self.rev:
            rate_string = self.rate.write_as_bngl(namespace) + ',' + self.rev_rate.write_as_bngl(namespace)
        else:
            rate_string = self.rate.write_as_bngl(namespace)
        delmol_string = '' if not self.delmol else " DeleteMolecules"
        return '%s %s %s %s%s' % (lhs_string, self.arrow, rhs_string, rate_string, delmol_string)

    def write_as_kappa(self):
        """Writes the rule as a Kappa string"""

        # possible actions  (root object is list (rule.lhs and rule.rhs))
        #  - iterable_item_added/removed (binding, unbinding)
        #  - type_changes (binding, unbinding)
        #  - value_changes (state change)
        lhs_string, rhs_string = '', ''
        if self.lhs:
            lhs_string = ','.join([p.write_as_kappa() for p in self.lhs])
        if self.rhs:
            rhs_string = ','.join([p.write_as_kappa() for p in self.rhs])
        if self.rev:
            rate_string = self.rate.write_as_kappa() + ',' + self.rev_rate.write_as_kappa()
        else:
            rate_string = self.rate.write_as_kappa()
        return '%s %s %s @ %s' % (lhs_string, self.arrow, rhs_string, rate_string)

    def __eq__(self, other):
        """Tests for Rule equality with another Rule.  Only considers LHS and RHS patterns
        and reversibility."""
        if isinstance(other, self.__class__):
            lhs_check = frozenset(self.lhs) == frozenset(other.lhs)
            rhs_check = frozenset(self.rhs) == frozenset(other.rhs)
            rev_check = self.rev == other.rev
            return lhs_check and rhs_check and rev_check
        else:
            return False

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self.lhs, self.rhs, self.rev, self.label))

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

    def write_as_bngl(self, namespace):
        """Writes Observable as BNGL string"""
        bname = namespace[self.name]
        return "%s %s %s" % (self.type, bname, ' '.join([p.write_as_bngl() for p in self.cpatterns]))

    def write_as_kappa(self, mdefs):
        """Writes Observable as Kappa string"""
        if self.type == 'Species':
            print "Kappa does not have a Species-like observable; printing '%s' as Molecules-like observable" % self.name

        obs_strs = []
        for p in self.cpatterns:
            # sorted for determinism (testing)
            kos = '+'.join(sorted(['|%s|' % x.write_as_kappa() for x in set(p.convert(mdefs))]))
            obs_strs.append(kos)

        obs = '+'.join(obs_strs)
        return '%%obs: \'%s\' %s' % (self.name, obs)

    def __repr__(self):
        return "Obs(name: %s, pattern: %s)" % (self.name, ' '.join([str(cp) for cp in self.cpatterns]))


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

        # Observable, Function, and Parameter instances may need to be renamed when converting from Kappa to BNGL
        self.convert_namespace = dict()  # Tracks names

    def write_as_bngl(self, file_name, dnp):
        """Writes Model as BNGL file"""
        s = 'begin model\n\nbegin parameters\n\n'
        for p in self.parameters:
            s += '\t%s\n' % p.write_as_bngl(self.convert_namespace)
        s += '\nend parameters\n\n'
        s += 'begin molecule types\n\n'
        for m in self.molecules:
            s += '\t%s\n' % m.write_as_bngl()
        s += '\nend molecule types\n\n'
        s += 'begin seed species\n\n'
        for i in self.initial_cond:
            s += '\t%s\n' % i.write_as_bngl(self.convert_namespace)
        s += '\nend seed species\n\n'
        s += 'begin observables\n\n'
        for o in self.observables:
            s += '\t%s\n' % o.write_as_bngl(self.convert_namespace)
        s += '\nend observables\n\n'
        s += 'begin functions\n\n'
        for f in self.functions:
            s += '\t%s\n' % f.write_as_bngl(self.convert_namespace)
        s += '\nend functions\n\n'
        s += 'begin reaction rules\n\n'
        for r in self.rules:
            s += '\t%s\n' % r.write_as_bngl(self.convert_namespace)
            if dnp:
                s += '\t%s\n' % r.write_as_bngl(dnp)
        s += '\nend reaction rules\n\n'
        s += 'end model\n'

        wf = open(file_name, 'w')
        wf.write(s)
        wf.close()

        logging.info("Wrote model to file (BNGL): %s" % file_name)

    # check for rules with molecular ambiguity
    def write_as_kappa(self, file_name, func_as_obs):
        """
        Writes Model as Kappa file

        Parameters
        ----------
        file_name : str
            Name of file to write
        func_as_obs : bool
            If True, writes functions as observables, otherwise as variables
        """
        s = ''
        for m in self.molecules:
            s += '%s\n' % m.write_as_kappa()
        if len(self.molecules) > 0:
            s += '\n'
        for p in self.parameters:
            s += '%s\n' % p.write_as_kappa()
        if len(self.parameters) > 0:
            s += '\n'
        for o in self.observables:
            s += '%s\n' % o.write_as_kappa()
        if len(self.observables) > 0:
            s += '\n'
        for f in self.functions:
            s += '%s\n' % f.write_as_kappa(func_as_obs)  # defaults to printing all "functions" as observables
        if len(self.functions) > 0:
            s += '\n'
        for i in self.initial_cond:
            s += '%s\n' % i.write_as_kappa()
        if len(self.initial_cond) > 0:
            s += '\n'
        for r in self.rules:
            s += '%s\n' % r.write_as_kappa()
        if len(self.rules) > 0:
            s += '\n'

        wf = open(file_name, 'w')
        wf.write(s)
        wf.close()

        logging.info("Wrote model to file (Kappa): %s" % file_name)

    def add_molecule_def(self, mol):
        """
        Adds a MoleculeDef to Model

        Parameters
        ----------
        mol : MoleculeDef
        """
        logging.debug("Added a molecule definition to the model: %s" % mol)
        self.molecules.append(mol)

    def add_init(self, init):
        """
        Adds an InitialCondition to Model

        Parameters
        ----------
        init : InitialCondition
        """
        logging.debug("Added an initial condition to the model: %s" % init)
        self.initial_cond.append(init)

    def add_obs(self, obs):
        """
        Adds an Observable to Model

        Parameters
        ----------
        obs : Observable
        """
        logging.debug("Added an observable to the model: %s" % obs)
        self.observables.append(obs)

    def add_func(self, func):
        """
        Adds a Function to Model

        Parameters
        ----------
        func : Function
        """
        logging.debug("Added a function to the model: %s" % func)
        self.functions.append(func)

    def add_rule(self, rule):
        """
        Adds a Rule to Model

        Parameters
        ----------
        rule : Rule
        """
        logging.debug("Added a rule to the model: %s" % rule)
        self.rules.append(rule)

    def add_parameter(self, param):
        """
        Adds a Parameter to Model

        Parameters
        ----------
        param : Parameter
        """
        logging.debug("Added a parameter to the model: %s" % param)
        self.parameters.append(param)


def is_number(n):
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
