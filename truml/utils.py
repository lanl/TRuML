"""truml.utils: miscellaneous support functions"""


def build_adj_list(mols):
    """
    Builds an adjacency list from a list of molecules

    Parameters
    ----------
    mols : list of Molecule instances

    Returns
    -------
    Adjacency list
    """
    al = []
    for i, m0 in enumerate(mols):
        bound = []
        for j, m1 in enumerate(mols):
            if i != j and m0.bound_to(m1):
                bound.append(j)
        al.append(bound)
    return al


def dfs(al, m, v):
    """
    Depth first search (recursive)

    Parameters
    ----------
    al : Adjacency list
    m : Node list
    v : Current node index

    Returns
    -------

    """
    m[v] = 1
    for x in al[v]:
        if m[x] == 0:
            dfs(al, m, x)
    return None


def get_connected_components(mols):
    """
    Determines explicitly connected components in a list of molecules

    Parameters
    ----------
    mols : list of Molecule instances

    Returns
    -------
    list of list of Molecule instances
    """
    al = build_adj_list(mols)
    visited = [0] * len(al)
    cmps = []
    for x in range(len(al)):
        if visited[x] == 1:
            continue
        prev_visited = list(visited)
        dfs(al, visited, x)
        cur_comp = []
        for i in range(len(visited)):
            if visited[i] == 1 and prev_visited[i] == 0:
                cur_comp.append(mols[i])
        cmps.append(cur_comp)
    return cmps


def flatten_pattern(cps):
    return [m for cp in cps for m in cp.molecule_list]
