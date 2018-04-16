"""Classes representing the semantics of rule-based modeling languages"""


import re
import itertools as it
import logging
import networkx as nx
import networkx.algorithms.isomorphism as iso
import rbexceptions
import utils

from copy import deepcopy
from math import factorial


class SiteDef:
    """A site definition composed of a name and a finite set of states"""

    def __init__(self, n, ss=[]):
        self.name = n
        self.state_list = ss

    def write_as_bngl(self):
        if self.state_list:
            return "%s~%s" % (self.name, '~'.join(self.state_list))
        else:
            return self.name

    def write_as_kappa(self):
        if self.state_list:
            return "%s{%s}" % (self.name, ','.join(self.state_list))
        else:
            return self.name

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
        sds : list
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
        if is_bngl:
            ss = [s.write_as_bngl() for s in self.sites]
        elif not is_bngl and not self.has_site_symmetry:
            ss = [s.write_as_kappa() for s in self.sites]
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


class MoleculeTemplate(object):
    def __init__(self):
        self.name = ''
        self.sites = []
        self.mdef = None

    @staticmethod
    def is_placeholder():
        return NotImplementedError("Function must be implemented in subclass")

    def write_as_bngl(self):
        return NotImplementedError("Function must be implemented in subclass")

    def write_as_kappa(self):
        return NotImplementedError("Function must be implemented in subclass")

    def has_same_interface(self, other):
        return NotImplementedError("Function must be implemented in subclass")

    def bound_to(self, other):
        return NotImplementedError("Function must be implemented in subclass")

    def _node_name(self):
        return NotImplementedError("Function must be implemented in subclass")

    def convert(self):
        return NotImplementedError("Function must be implemented in subclass")

    def has_identical_sites(self):
        return NotImplementedError("Function must be implemented in subclass")


class PlaceHolderMolecule(MoleculeTemplate):
    def __init__(self):
        super(PlaceHolderMolecule, self).__init__()

    @staticmethod
    def write_as_kappa():
        return '.'

    @staticmethod
    def write_as_bngl():
        return None

    def __eq__(self, other):
        return isinstance(other, PlaceHolderMolecule)

    @staticmethod
    def is_placeholder():
        return True

    def has_same_interface(self, other):
        return True

    def bound_to(self, other):
        return False

    def _node_name(self):
        return '.'

    def convert(self):
        return [self]

    def has_identical_sites(self):
        return False

    def __repr__(self):
        return "PlaceHolderMolecule"


class Molecule(MoleculeTemplate):
    """
    An individual molecule/agent inside a pattern.
        Note that it may not contain all sites listed in the associated
        MoleculeDef class
    """

    def __init__(self, name, sites, md):
        """
        Molecule initialization function. Sites are sorted by a predefined index

        Parameters
        ----------
        name : str
            The name of the molecule type
        sites : list
            A list of Sites that appear in the pattern
        md : MoleculeDef
            The corresponding MoleculeDef for this Molecule

        """
        super(Molecule, self).__init__()
        self.name = name
        self.sites = sorted(sites, key=lambda s: s.index)  # list of Sites
        self.mdef = md

    @staticmethod
    def is_placeholder():
        return False

    def _node_name(self):
        """
        Provides a unique label for the Molecule based on its sites' states and bonds.

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

    def _enumerate_site(self, site_name, index, need_state=False):
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
            for s in self.mdef.sites:
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

    def convert(self):
        """
        Converts a molecule that may have multiple identically named sites
        into a list of molecules compatible with Kappa semantics that
        enumerate all molecule configurations compatible with the original
        BNGL molecule configuration

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
                un_configs_per_site[s.name][s] = []
            un_configs_per_site[s.name][s].append(s.index)

        def rename_site(name, site, index):
            return Site(name, index, s=site.state, b=site.bond)

        def rename_sites(names, site, idcs):
            return tuple([rename_site(name, site, index) for name, index in zip(names, idcs)])

        # Check for the possibility of overlapping patterns
        possible_overlap = {k: False for k in un_configs_per_site.keys()}
        for k in un_configs_per_site.keys():
            num_identical_sites = len(self.mdef.inv_site_name_map[k])
            if num_identical_sites > 1 and k in un_configs_per_site.keys():
                num_present_sites = sum([len(idcs) for idcs in un_configs_per_site[k].values()])
                if num_identical_sites > num_present_sites >= 1:
                    possible_overlap[k] = True
                    break

        k_configs = {}
        for sn in un_configs_per_site.keys():
            k_configs[sn] = []
            k_sn_names = set(self.mdef.inv_site_name_map[sn])
            cur_combs = []

            for s, idcs in un_configs_per_site[sn].iteritems():
                if len(cur_combs) == 0:
                    cur_combs = [rename_sites(names, s, idcs) for names in it.combinations(k_sn_names, len(idcs))]
                else:
                    tmp_combs = []
                    for cc in cur_combs:
                        rem_names = k_sn_names - set(map(lambda l: l.name, cc))
                        new_combs = [rename_sites(names, s, idcs) for names in it.combinations(rem_names, len(idcs))]
                        for nc in new_combs:
                            tmp_combs.append(cc + nc)
                    cur_combs = tmp_combs

            if possible_overlap[sn]:
                need_state = self._site_state_present(un_configs_per_site[sn])
                num_rem_sites = len(self.mdef.inv_site_name_map[k]) - sum([len(idcs) for idcs in un_configs_per_site[sn].values()])
                indices = range(len(self.sites), len(self.sites) + num_rem_sites)

                for idx in indices:
                    possible_sites = self._enumerate_site(sn, idx, need_state)
                    tmp_combs = []
                    for cc in cur_combs:
                        rem_names = k_sn_names - set(map(lambda l: l.name, cc))
                        new_combs = [rename_site(x, y, idx) for x, y in it.product(rem_names, possible_sites)]
                        for nc in new_combs:
                            tmp_combs.append(cc + (nc,))
                    cur_combs = tmp_combs

            k_configs[sn] = cur_combs

        k_prod = it.product(*k_configs.values())

        return sorted([Molecule(self.name, [e for t in tt for e in sorted(t)], self.mdef) for tt in k_prod])

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

    def has_same_interface(self, other):
        if isinstance(other, PlaceHolderMolecule):
            return True
        elif isinstance(other, self.__class__):
            return self.name == other.name and \
                   sorted([s.name for s in self.sites]) == sorted([s.name for s in other.sites])
        else:
            return False

    @staticmethod
    def _diff_quant(d):
        if d == (-1, -1):
            return 0
        elif d[0] == -1:
            return 1
        elif d[1] == -1:
            return 2
        else:
            return 3

    def interface_diff_map(self, other):
        imap = dict()
        if len(self.sites) != len(other.sites):
            return None

        used_other_idcs = []
        for s in self.sites:
            imap[s] = None
            possible = []
            # Map to most similar site
            for j, t in enumerate(other.sites):
                diff = s.diff(t)
                # Check to see if there is any difference between sites s and t
                if t.name == s.name and j not in used_other_idcs and diff == (-1, -1):
                    used_other_idcs.append(j)
                    imap[s] = (-1, -1)
                    possible = []
                    break
                elif t.name == s.name and j not in used_other_idcs:
                    possible.append((j, self._diff_quant(diff)))
            if imap[s] is None and len(possible) == 0:
                return None  # there is no match for site s in other.sites
            elif imap[s] is None:
                sd = sorted(possible, key=lambda l: l[1])
                idx = sd[0][0]
                imap[s] = s.diff(other.sites[idx])
                used_other_idcs.append(idx)
        if len(used_other_idcs) < len(other.sites):
            return None  # there are unmatched sites in other
        else:
            return {k: v for k, v in imap.iteritems() if v != (-1, -1)}

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
            raise rbexceptions.NotConvertedException(self.write_as_bngl())
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

    @staticmethod
    def is_placeholder():
        return False


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
        return self.name if self.state is None else '%s~%s' % (self.name, self.state)

    def write_as_bngl(self):
        """
        Builds site string in BNGL syntax

        Returns
        -------
        str
            Site written as BNGL string
        """
        state = '' if self.state is None else '~?' if self.state == 'WILD' else '~%s' % self.state
        bond = '' if self.bond is None else self.bond.write_as_bngl()
        return self.name + state + bond

    def write_as_kappa(self):
        """
        Builds site string in Kappa syntax

        Returns
        -------
        str
            Site written as Kappa string
        """
        state = '' if self.state is None else '{#}' if self.state == 'WILD' else '{%s}' % self.state
        bond = '[.]' if self.bond is None else self.bond.write_as_kappa()
        return self.name + state + bond

    def diff(self, other):
        """
        Provides a 2-tuple describing the difference in state and bond.  Not symmetric

        Parameters
        ----------
        other : Site
            A product site in a corresponding Molecule

        Returns
        -------
        tuple
            Contains 2 elements containing information about site state and site bond, respectively.
            The first element is the other site state, and the second element is the other site Bond.
        """
        diff_tuple = [-1, -1]
        if self.state != other.state:
            diff_tuple[0] = other.state
        if self.bond != other.bond:
            diff_tuple[1] = other.bond
        return tuple(diff_tuple)

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
        if self.num >= 0:
            return '!%s' % self.num
        if self.wild:
            return '!+'
        elif self.any:
            return '!?'
        return ''

    def write_as_kappa(self):
        """Write bond as Kappa string"""
        if self.wild:
            return '[_]'
        elif self.any:
            return '[#]'
        else:
            return '[%s]' % self.num
        return '[.]'

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
        self.placeholder = (len(self.molecule_list) == 1) and isinstance(self.molecule_list[0], PlaceHolderMolecule)

    def __getitem__(self, item):
        return self.molecule_list[item]

    def __len__(self):
        return len(self.molecule_list)

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
                for s in self[i].sites:
                    if s.bond is None:
                        continue
                    for s2 in self[j].sites:
                        if s2.bond is None:
                            continue
                        if s.bond.num == s2.bond.num and s.bond.num >= 0:
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
        Determines the number of automorphisms in a :class:`~truml.objects.CPattern`

        Only performed after conversion to Kappa compatible site names.  Thus
        the count will only include automorphisms on the level of molecules

        Returns
        -------
        int
            The number of automorphisms in this CPattern
        """
        if self.placeholder:
            return 1

        # Check to make sure conversion to Kappa compatible site names has occurred
        for m in self.molecule_list:
            if m.has_identical_sites():
                raise rbexceptions.NotConvertedException(m.write_as_bngl())
        # If all molecules are unique, exit with count 1.  Otherwise calculate the
        # number of automorphisms
        if len([m._node_name() for m in self.molecule_list if not m.is_placeholder()]) == \
                len(set([m.name for m in self.molecule_list if not m.is_placeholder()])):
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
    def convert(self):
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
            k_str_mol_list.append(m.convert())
        k_patterns = list(it.product(*k_str_mol_list))
        for pat in k_patterns:
            if len(utils.get_connected_components(list(pat))) > 1 and self.num_molecules() > 1:
                logging.warning("CPattern components are not explicitly connected in '%s'" % self.write_as_bngl())
                logging.warning(
                    "The above pattern may be attempting to detect polymers which is not possible in Kappa")
                logging.warning(
                    "Please specify explicit binding if possible"
                )
                raise rbexceptions.NotCompatibleException(
                    "Pattern '%s' cannot be converted to Kappa-compatible syntax" % self.write_as_bngl())

        return [CPattern(pat) for pat in k_patterns]

    def is_isomorphic(self, other):
        if isinstance(other, self.__class__):
            if self.num_molecules() != other.num_molecules():
                return False
            else:
                sg = self._build_graph()
                og = other._build_graph()
                nm = iso.categorical_node_match('name', '')
                em = iso.categorical_edge_match('name', '')
                return iso.is_isomorphic(sg, og, node_match=nm, edge_match=em)
        return False

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
                if m.is_placeholder():
                    continue
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
        return "CPattern(" + '--'.join([str(x) for x in self.molecule_list]) + ")"


class CPatternList:
    """List of CPatterns"""

    def __init__(self, cps):
        self.cpatterns = cps

    def __getitem__(self, item):
        return self.cpatterns[item]

    def __len__(self):
        return len(self.cpatterns)

    def append(self, cp):
        self.cpatterns.append(cp)

    def convert(self):
        c_cps = []
        for cp in self.cpatterns:
            c_cps.append(cp.convert())
        return list(it.imap(lambda p: CPatternList(list(p)), it.product(*c_cps)))

    def write_as_bngl(self, dot):
        all_placeholder = True
        for cp in self.cpatterns:
            all_placeholder = all_placeholder and cp.placeholder

        if all_placeholder:
            return '0'

        valid_cpatterns = [cp for cp in self.cpatterns if not cp.placeholder]
        joiner = '.' if dot else '+'
        return joiner.join([cp.write_as_bngl() for cp in valid_cpatterns])

    def write_as_kappa(self):
        return ','.join([cp.write_as_kappa() for cp in self.cpatterns])

    def __str__(self):
        return "[%s]" % ','.join([str(cp) for cp in self.cpatterns])

    def __repr__(self):
        return str(self)


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

    def convert(self):
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
        ss = self.species.convert()
        lss = len(ss)
        if lss == 1:
            amount = self.amount
        else:
            amount = self.amount / float(lss) if self.amount_is_number else Expression(
                ['('] + self.amount.expr + [')', '/', lss])
        return [InitialCondition(s, amount, self.amount_is_number) for s in ss]

    def write_as_bngl(self, namespace=dict()):
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

    def write_as_bngl(self, namespace=dict()):
        """Writes Parameter as BNGL string"""
        bname = namespace[self.name] if self.name in namespace.keys() else self.name
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

    def write_as_bngl(self, namespace=dict()):
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
        expr = []

        i = 0
        while (i < len(self.atom_list)):
            a = self.atom_list[i]
            if a in bngl_to_kappa_func_map.keys():
                trig_func_match = re.compile('sinh|cosh|tanh|asinh|acosh|atanh')
                if re.match('log', a) or re.match(trig_func_match, a):
                    expr.append(bngl_to_kappa_func_map[a](self.atom_list[i + 2]))
                    i += 4
                else:
                    expr.append(bngl_to_kappa_func_map[a])
            elif re.match('[A-Za-z]', a):
                expr.append('\'%s\'' % a)
            elif re.match('[+-]', a):
                expr.append(' %s ' % a)
            else:
                expr.append(a)
            i += 1

        return ''.join(expr)

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

    def write_as_bngl(self, namespace=dict()):
        """Writes function as BNGL string"""
        bname = namespace[self.name] if self.name in namespace.keys() else self.name
        return '%s=%s' % (bname, self.expr.write_as_bngl(namespace))

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


# TODO check to make sure rate doesn't have parentheses
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

    def contains_variable(self, var):
        """
        Checks to see if the Rate contains a particular named variable.

        Parameters
        ----------
        var : str
            String corresponding to the name of an Observable, Function, or Parameter instance
        Returns
        -------
        bool
            True if Rate involves the named variable, False otherwise
        """
        if isinstance(self.rate, float):
            return False
        elif isinstance(self.rate, str):
            return self.rate == var
        elif isinstance(self.rate, Expression):
            return var in self.rate.atom_list

    def write_as_bngl(self, namespace=dict()):
        """Write Rate as BNGL string"""
        try:
            return self.rate.write_as_bngl(namespace)
        except AttributeError:
            if isinstance(self.rate, str) and self.rate in namespace.keys():
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
        lhs : CPatternList
            The reactants
        rhs : CPatternList
            The products
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

        self.lhs_mols = utils.flatten_pattern(self.lhs)
        self.rhs_mols = utils.flatten_pattern(self.rhs)

        # This is set to None once an appropriate mapping is built
        self.mol_map = self._build_mol_map(self.lhs_mols, self.rhs_mols)

        if self._placeholder_check():
            self.insert_placeholders()

    def _placeholder_check(self):
        if len(self.lhs_mols) != len(self.rhs_mols):
            return True
        elif None in self.mol_map.values():
            return True
        elif set(range(len(self.rhs_mols))) != set(self.mol_map.values()):
            return True
        self.mol_map = None
        return False

    def insert_placeholders(self):
        """Inserts PlaceHolderMolecule instances into rules originally written in BNGL for conversion to Kappa"""

        lmol2cp = utils.flatten_pattern_todict(self.lhs)
        rmol2cp = utils.flatten_pattern_todict(self.rhs)

        lhs_list, rhs_list = [], []
        for li in sorted(self.mol_map.keys()):
            lhs_list.append(self.lhs_mols[li])
            if self.mol_map[li] is not None:
                rhs_list.append(self.rhs_mols[self.mol_map[li]])
            else:
                rhs_list.append(PlaceHolderMolecule())
                # Update Molecule to CPattern mapping to accommodate new Molecule (and CPattern)
                if len(rmol2cp.keys()) < li + 1:
                    rmol2cp[li] = li
                else:
                    for ri in reversed(sorted(rmol2cp.keys())):
                        if ri >= li:
                            rmol2cp[ri + 1] = rmol2cp[ri] + 1
                        if li == ri:
                            rmol2cp[ri] = ri
                        if ri < li:
                            break

        mapped_rhs_idcs = set(it.ifilterfalse(lambda l: l is None, self.mol_map.values()))
        unmapped_rhs_idcs = set(range(len(self.rhs_mols))) - mapped_rhs_idcs

        maxi = len(lmol2cp)
        for ri in unmapped_rhs_idcs:
            lhs_list.append(PlaceHolderMolecule())

            # Update Molecule to CPattern mapping to accommodate new Molecule (and CPattern)
            lmol2cp[maxi] = maxi
            maxi += 1

            rhs_list.append(self.rhs_mols[ri])

        lhs_cps, rhs_cps = {v: [] for v in set(lmol2cp.values())}, {v: [] for v in set(rmol2cp.values())}
        for i in lmol2cp.keys():
            lhs_cps[lmol2cp[i]].append(lhs_list[i])
        for i in rmol2cp.keys():
            rhs_cps[rmol2cp[i]].append(rhs_list[i])

        self.lhs = CPatternList([CPattern(lhs_cps[k]) for k in sorted(lhs_cps.keys())])
        self.rhs = CPatternList([CPattern(rhs_cps[k]) for k in sorted(rhs_cps.keys())])
        self.lhs_mols = utils.flatten_pattern(self.lhs)
        self.rhs_mols = utils.flatten_pattern(self.rhs)
        # Sanity check
        assert len(self.lhs_mols) == len(self.rhs_mols)
        self.mol_map = None

    def convert(self):
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
        logging.debug("Attempting to convert rule: %s" % self.write_as_bngl(from_bngl=True))
        rs = []
        converted_lhss = self.lhs.convert()

        for conv_lhs in converted_lhss:
            converted_rhss = self._build_actions().apply(conv_lhs)
            for conv_rhs in converted_rhss:
                rs.append(Rule(conv_lhs, conv_rhs, self.rate, rev=self.rev, rev_rate=self.rev_rate, label=self.label,
                               delmol=self.delmol))

        un_rules = [rs[0]]
        for rule in rs[1:]:
            is_isomorphic = False
            for un_rule in un_rules:
                if rule.is_isomorphic(un_rule):
                    is_isomorphic = True
                    break
            if not is_isomorphic:
                un_rules.append(rule)

        logging.debug("Converted into %s Kappa rule(s):\n\t%s" % (len(un_rules), '\n\t'.join([x.write_as_kappa() for x in un_rules])))

        return un_rules

    def is_isomorphic(self, other):
        if isinstance(other, self.__class__):
            if not (len(self.lhs) == len(other.lhs) and len(self.rhs) == len(other.rhs)):
                return False

            for i in range(len(self.lhs)):
                if not self.lhs[i].is_isomorphic(other.lhs[i]):
                    return False

            for i in range(len(self.rhs)):
                if not self.rhs[i].is_isomorphic(other.rhs[i]):
                    return False

            return True
        else:
            return False

    def _build_actions(self, rev=False):
        """Builds a list of Action instances corresponding to the differences between the reactant
        list of Molecule instances and the product list of Molecule instances"""
        action_list = []
        if rev:
            lhs_mols = self.rhs_mols
            rhs_mols = self.lhs_mols
        else:
            lhs_mols = self.lhs_mols
            rhs_mols = self.rhs_mols

        to_synth = []
        for i in range(len(lhs_mols)):
            if isinstance(rhs_mols[i], PlaceHolderMolecule):
                action_list.append(Degradation(i))
                continue
            elif isinstance(lhs_mols[i], PlaceHolderMolecule):
                to_synth.append(i)
                continue

            smap = lhs_mols[i].interface_diff_map(rhs_mols[i])
            mdef = lhs_mols[i].mdef
            for k in smap.keys():
                diff = smap[k]
                if diff[0] != -1 and diff[1] != -1:
                    action_list.append(BondAndStateChange(i, k, diff[0], diff[1], mdef))
                else:
                    if diff[0] != -1:
                        action_list.append(StateChange(i, k, diff[0], mdef))
                    if diff[1] != -1:
                        action_list.append(BondChange(i, k, diff[1], mdef))

        if len(to_synth) > 0:
            conn_cmps = utils.get_connected_components([rhs_mols[i] for i in to_synth])
            for c in conn_cmps:
                action_list.append(Synthesis(CPattern(c)))

        return MultiAction(action_list)

    @staticmethod
    def _build_mol_map(lhs, rhs):
        """Builds a map between Molecule instances on the LHS and RHS of a rule"""
        mmap = dict()
        used_rhs_mol_idcs = []
        for i, lm in enumerate(lhs):
            mmap[i] = None
            for j, rm in enumerate(rhs):
                if lm.has_same_interface(rm) and j not in used_rhs_mol_idcs:
                    mmap[i] = j
                    used_rhs_mol_idcs.append(j)
                    break
        return mmap

    def _unique_reactant_indices(self):
        actual_cpattern_idcs = [i for i in range(len(self.lhs.cpatterns)) if not self.lhs.cpatterns[i].placeholder]
        if len(actual_cpattern_idcs) > 0:
            unique_reactant_first_idcs = {0: 1}
        else:
            return None
        for l1 in actual_cpattern_idcs[1:]:
            unique = True
            for u in unique_reactant_first_idcs:
                if self.lhs.cpatterns[l1].is_isomorphic(self.lhs.cpatterns[u]):
                    unique = False
                    unique_reactant_first_idcs[u] += 1
                    break
            if unique:
                unique_reactant_first_idcs[l1] = 1
        return unique_reactant_first_idcs

    def rate_factor(self, b2k):
        auto_factor = 1.0
        for lhs_patt in self.lhs.cpatterns:
            auto_factor *= lhs_patt.automorphisms()

        reactant_counts = self._unique_reactant_indices()
        if reactant_counts:
            multiple_reactant_factor = 1
            for k, v in reactant_counts.iteritems():
                multiple_reactant_factor *= factorial(v)

            factor = auto_factor * multiple_reactant_factor

            return 1.0 / factor if b2k else factor
        return 1.0

    def write_as_bngl(self, namespace=dict(), dot=False, from_bngl=False):
        """Writes the rule as a BNGL string"""
        lhs_string = self.lhs.write_as_bngl(dot)
        rhs_string = self.rhs.write_as_bngl(dot)
        rfactor_string = ''
        if not from_bngl:
            rfactor = self.rate_factor(b2k=False)
            rfactor_string = '' if rfactor == 1.0 else ' * %s' % rfactor

        if self.rev:
            rate_string = self.rate.write_as_bngl(namespace) + rfactor_string + ',' + \
                          self.rev_rate.write_as_bngl(namespace) + rfactor_string
        else:
            rate_string = self.rate.write_as_bngl(namespace) + rfactor_string
        delmol_string = '' if not self.delmol else " DeleteMolecules"
        return '%s %s %s %s%s' % (lhs_string, self.arrow, rhs_string, rate_string, delmol_string)

    def write_as_kappa(self, from_kappa=False):
        """Writes the rule as a Kappa string"""

        # possible actions  (root object is list (rule.lhs and rule.rhs))
        #  - iterable_item_added/removed (binding, unbinding)
        #  - type_changes (binding, unbinding)
        #  - value_changes (state change)

        lhs_string = self.lhs.write_as_kappa()
        rhs_string = self.rhs.write_as_kappa()
        rfactor_string = ''
        if not from_kappa:
            rfactor = self.rate_factor(b2k=True)
            rfactor_string = '' if rfactor == 1.0 else ' * %s' % rfactor

        if self.rev:
            rate_string = self.rate.write_as_kappa() + rfactor_string + ',' + self.rev_rate.write_as_kappa() + rfactor_string
        else:
            rate_string = self.rate.write_as_kappa() + rfactor_string
        return '%s %s %s @ %s' % (lhs_string, self.arrow, rhs_string, rate_string)

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
        ps : CPatternList
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

    @staticmethod
    def _calc_factor(cp):
        # assumes conversion to Kappa-compatible syntax
        f = 1
        for mol in cp:
            un_site_names = set([mol.mdef.site_name_map[s.name] for s in mol.sites])
            un_configs_per_site = {s: {} for s in un_site_names}
            for site in mol.sites:
                bngl_site_name = mol.mdef.site_name_map[site.name]
                bngl_site = Site(bngl_site_name, site.index, s=site.state, b=site.bond)
                if bngl_site not in un_configs_per_site[bngl_site_name]:
                    un_configs_per_site[bngl_site_name][bngl_site] = 1
                else:
                    un_configs_per_site[bngl_site_name][bngl_site] += 1
            m_int_symm = 1
            for d in un_configs_per_site.values():
                for n in d.values():
                    m_int_symm *= factorial(n)
            f *= m_int_symm
        return f

    def _ftos(self, cp):
        f = self._calc_factor(cp)
        return '' if f == 1 else '%s*' % str(f)

    def write_as_bngl(self, namespace=dict()):
        """Writes Observable as BNGL string"""
        bname = namespace[self.name] if self.name in namespace.keys() else self.name
        return "%s %s %s" % (self.type, bname, ' '.join([p.write_as_bngl() for p in self.cpatterns]))

    def write_as_kappa(self):
        """Writes Observable as Kappa string"""
        if self.type == 'Species':
            logging.warning(
                "Kappa does not have a Species-like observable; printing '%s' as Molecules-like observable" % self.name)

        obs_strs = []
        for p in self.cpatterns:
            kps = p.convert()
            un_kps = [kps[0]]
            for kp in kps[1:]:
                is_isomorphic = False
                for un_kp in un_kps:
                    if un_kp.is_isomorphic(kp):
                        is_isomorphic = True
                        break
                if not is_isomorphic:
                    un_kps.append(kp)

            # sorted for determinism (testing)
            kos = '+'.join(sorted(['%s|%s|' % (self._ftos(x), x.write_as_kappa()) for x in un_kps]))
            obs_strs.append(kos)

        obs = '+'.join(obs_strs)
        return '%%obs: \'%s\' %s' % (self.name, obs)

    def __repr__(self):
        return "Obs(name: %s, pattern: %s)" % (self.name, ' '.join([str(cp) for cp in self.cpatterns]))


class Model:
    """Object for rule-based model"""

    def __init__(self, bngl):
        """Model initialization function"""
        self.molecules = []
        self.initial_cond = []
        self.observables = []
        self.functions = []
        self.rules = []
        self.parameters = []
        self.bngl = bngl  # True if reading BNGL, False if reading Kappa

        # Observable, Function, and Parameter instances may need to be renamed when converting from Kappa to BNGL
        self.convert_namespace = dict()  # Tracks names

        # Observable or Function names that cannot be translated to Kappa.
        self.invalid_names = set()

    def write_as_bngl(self, file_name, dnp):
        """Writes Model as BNGL file"""
        logging.debug("Writing model to BNGL file: '%s'" % file_name)

        s = 'begin model\n\nbegin parameters\n\n'
        for p in self.parameters:
            logging.debug("Writing parameter %s to file" % p)
            s += '\t%s\n' % p.write_as_bngl(self.convert_namespace)
        s += '\nend parameters\n\n'
        s += 'begin molecule types\n\n'
        for m in self.molecules:
            logging.debug("Writing molecule type %s to file" % m)
            s += '\t%s\n' % m.write_as_bngl()
        s += '\nend molecule types\n\n'
        s += 'begin seed species\n\n'
        for i in self.initial_cond:
            logging.debug("Writing initial condition %s to file" % i)
            s += '\t%s\n' % i.write_as_bngl(self.convert_namespace)
        s += '\nend seed species\n\n'
        s += 'begin observables\n\n'
        for o in self.observables:
            logging.debug("Writing observable %s to file" % o)
            s += '\t%s\n' % o.write_as_bngl(self.convert_namespace)
        s += '\nend observables\n\n'
        s += 'begin functions\n\n'
        for f in self.functions:
            logging.debug("Writing function %s to file" % f)
            s += '\t%s\n' % f.write_as_bngl(self.convert_namespace)
        s += '\nend functions\n\n'
        s += 'begin reaction rules\n\n'
        for r in self.rules:
            logging.debug("Writing rule %s to file" % r)
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
        logging.debug("Writing model to Kappa file: '%s'" % file_name)

        s = ''
        for m in self.molecules:
            logging.debug("Writing molecule type %s to file" % m)
            km = m.convert()
            s += '%s\n' % km.write_as_kappa()
        if len(self.molecules) > 0:
            s += '\n'

        for p in self.parameters:
            logging.debug("Writing parameter %s to file" % p)
            s += '%s\n' % p.write_as_kappa()
        if len(self.parameters) > 0:
            s += '\n'

        for o in self.observables:
            try:
                logging.debug("Writing observable %s to file" % o)
                s += '%s\n' % o.write_as_kappa()
            except rbexceptions.NotCompatibleException:
                self.invalid_names.add(o.name)
                logging.warning("Incompatible observable '%s' not converted to Kappa" % o.write_as_bngl())
        if len(self.observables) > 0:
            s += '\n'

        for f in self.functions:
            valid_func = True
            for atom in f.expr.atom_list:
                if atom in self.invalid_names:
                    self.invalid_names.add(f.name)
                    valid_func = False
                    logging.warning("Incompatible function '%s' not converted to Kappa" % f.write_as_bngl())
            if valid_func:
                logging.debug("Writing function %s to file" % f)
                s += '%s\n' % f.write_as_kappa(func_as_obs)  # defaults to printing all "functions" as observables
        if len(self.functions) > 0:
            s += '\n'

        for i in self.initial_cond:
            kis = i.convert()
            logging.debug("Writing initial condition %s to file" % i)
            s += '%s\n' % '\n'.join([ki.write_as_kappa() for ki in kis])
        if len(self.initial_cond) > 0:
            s += '\n'

        for r in self.rules:
            try:
                for inv in self.invalid_names:
                    if r.rate.contains_variable(inv):
                        logging.critical("Rule's rate contains an incompatible variable")
                        raise rbexceptions.NotCompatibleException(
                            "Rate '%s' contains an incompatible variable" % r.rate.write_as_bngl())
                    elif r.rev_rate:
                        if r.rev_rate.contains_variable(inv):
                            logging.critical("Rule's reverse rate contains an incompatible variable")
                            raise rbexceptions.NotCompatibleException(
                                "Reverse rate '%s' contains an incompatible variable" % r.rev_rate.write_a_bngl())
                krs = r.convert()
                logging.debug("Writing rule %s to file" % r)
                s += '%s\n' % '\n'.join([kr.write_as_kappa() for kr in krs])
            except rbexceptions.NotCompatibleException as nce:
                logging.critical("Cannot convert rule '%s' to Kappa" % r.write_as_bngl())
                raise nce
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


# when reading from BNGL files with identical site names, this schema will have Action objects
# that contain site names corresponding to the BNGL format.  Thus, the molecule definition
# needs to be included to detect and rewrite patterns that have converted to Kappa-compatible
# site names
class Action(object):
    """
    Abstract class that defines an action that when applied to a CPatternList,
    results in a list of CPatternList instances
    """
    def __init__(self):
        pass

    def apply(self, cps):
        NotImplementedError("apply is not implemented")


class BondChange(Action):
    """Action subclass that defines an bond action on a site"""
    def __init__(self, i, s, nb, md):
        super(BondChange, self).__init__()
        self.mol_index = i
        self.site = s
        self.new_bond = nb
        self.molecule_def = md

    def apply(self, cps):
        cps_copy = deepcopy(cps)
        mols = utils.flatten_pattern(cps_copy)
        applications = []
        for s in mols[self.mol_index].sites:
            try:
                bname = self.molecule_def.site_name_map[s.name]
            except KeyError:
                bname = s.name
            tmp_site = Site(bname, s.index, s=s.state, b=s.bond)
            if tmp_site == self.site:
                mols_copy = deepcopy(mols)
                mols_copy[self.mol_index].sites[s.index].bond = self.new_bond
                new_cps = CPatternList([CPattern(x) for x in utils.get_connected_components(mols_copy)])
                applications.append(new_cps)
        return applications

    def __str__(self):
        return "BondChange(%s, %s, %s)" % (self.mol_index, self.site, self.new_bond)

    def __repr__(self):
        return str(self)


class StateChange(Action):
    """Action subclass that defines a change in a Site instance's state"""
    def __init__(self, i, s, ns, md):
        super(StateChange, self).__init__()
        self.mol_index = i
        self.site = s
        self.new_state = ns
        self.molecule_def = md

    def apply(self, cps):
        cps_copy = deepcopy(cps)
        mols = utils.flatten_pattern(cps_copy)
        applications = []
        for s in mols[self.mol_index].sites:
            try:
                bname = self.molecule_def.site_name_map[s.name]
            except KeyError:
                bname = s.name
            tmp_site = Site(bname, s.index, s=s.state, b=s.bond)
            if tmp_site == self.site:
                mols_copy = deepcopy(mols)
                mols_copy[self.mol_index].sites[s.index].state = self.new_state
                new_cps = CPatternList([CPattern(x) for x in utils.get_connected_components(mols_copy)])
                applications.append(new_cps)

        return applications

    def __str__(self):
        return "StateChange(%s, %s, %s)" % (self.mol_index, self.site, self.new_state)

    def __repr__(self):
        return str(self)


class Degradation(Action):
    """Action subclass that defines the removal of a Molecule instance"""
    def __init__(self, i):
        super(Degradation, self).__init__()
        self.mol_index = i

    @staticmethod
    def _filter_explicit_bonds(bond):
        if bond is None:
            return False
        elif bond.num < 0:
            return False
        else:
            return True

    def apply(self, cps):
        cps_copy = deepcopy(cps)
        mols = utils.flatten_pattern(cps_copy)
        mols.pop(self.mol_index)
        return [CPatternList([CPattern(x) for x in utils.get_connected_components(mols)])]

    def __str__(self):
        return "Degradation(%s)" % self.mol_index

    def __repr__(self):
        return str(self)


class Synthesis(Action):
    """Action subclass that defines the addition of a CPattern instance"""
    def __init__(self, cp):
        super(Synthesis, self).__init__()
        self.cpattern = cp

    def apply(self, cps):
        cps_copy = deepcopy(cps)
        cps_copy.append(self.cpattern)
        return [CPatternList(cps_copy)]

    def __str__(self):
        return "Synthesis(%s)" % self.cpattern

    def __repr__(self):
        return str(self)


class BondAndStateChange(Action):
    def __init__(self, i, s, ns, nb, md):
        super(BondAndStateChange, self).__init__()
        self.mol_index = i
        self.site = s
        self.new_state = ns
        self.new_bond = nb
        self.molecule_def = md

    def apply(self, cps):
        cps_copy = deepcopy(cps)
        mols = utils.flatten_pattern(cps_copy)
        applications = []
        for s in mols[self.mol_index].sites:
            try:
                bname = self.molecule_def.site_name_map[s.name]
            except KeyError:
                bname = s.name
            tmp_site = Site(bname, s.index, s=s.state, b=s.bond)
            if tmp_site == self.site:
                mols_copy = deepcopy(mols)
                mols_copy[self.mol_index].sites[s.index].state = self.new_state
                mols_copy[self.mol_index].sites[s.index].bond = self.new_bond
                new_cps = CPatternList([CPattern(x) for x in utils.get_connected_components(mols_copy)])
                applications.append(new_cps)

        return applications

    def __str__(self):
        return "BondAndStateChange(%s, %s, %s, %s)" % (self.mol_index, self.site, self.new_state, self.new_bond)

    def __repr__(self):
        return str(self)


class MultiAction(Action):
    """Class that contains a list of Action instances to be applied to a CPattern or list of CPattern instances"""
    def __init__(self, ls):
        super(MultiAction, self).__init__()
        self.action_list = self._order_actions(ls)

    @staticmethod
    def _order_actions(ls):
        ordered_action_list = []
        for action in ls:
            if isinstance(action, Degradation):
                ordered_action_list.append(action)
            else:
                ordered_action_list.insert(0, action)
        return ordered_action_list

    def apply(self, cps):
        cpss = [deepcopy(cps)]  # list of CPatternList instances
        for action in self.action_list:
            cpsss = []
            for cpsi in cpss:
                cpsss.append(action.apply(cpsi))  # 3-D list
            cpss = [j for i in cpsss for j in i]
        return cpss

    def __len__(self):
        return len(self.action_list)

    def __getitem__(self, item):
        if isinstance(item, (int, slice)):
            return self.action_list[item]
        else:
            raise TypeError

    def __str__(self):
        return "MultiAction\n\t%s" % '\n\t'.join([str(action) for action in self.action_list])

    def __repr__(self):
        return str(self)


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
