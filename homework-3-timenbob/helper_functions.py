"""Helper functions for HW3"""
import numpy as np
from copy import deepcopy
from matplotlib.axes import Axes
import matplotlib.pyplot as plt

class Node:
    def __init__(
        self,
        name: str,
        left: "Node",
        left_distance: float,
        right: "Node",
        right_distance: float,
        confidence: float = None,
    ):
        """A node in a binary tree produced by neighbor joining algorithm.

        Parameters
        ----------
        name: str
            Name of the node.
        left: Node
            Left child.
        left_distance: float
            The distance to the left child.
        right: Node
            Right child.
        right_distance: float
            The distance to the right child.
        confidence: float
            The confidence level of the split determined by the bootstrap method.
            Only used if you implement Bonus Problem 1.

        Notes
        -----
        The current public API needs to remain as it is, i.e., don't change the
        names of the properties in the template, as the tests expect this kind
        of structure. However, feel free to add any methods/properties/attributes
        that you might need in your tree construction.

        """
        self.name = name
        self.left = left
        self.left_distance = left_distance
        self.right = right
        self.right_distance = right_distance
        self.confidence = confidence


def neighbor_joining(distances: np.ndarray, labels: list) -> Node:
    """The Neighbor-Joining algorithm.

    For the same results as in the later test dendrograms;
    add new nodes to the end of the list/matrix and
    in case of ties, use np.argmin to choose the joining pair.

    Parameters
    ----------
    distances: np.ndarray
        A 2d square, symmetric distance matrix containing distances between
        data points. The diagonal entries should always be zero; d(x, x) = 0.
    labels: list
        A list of labels corresponding to entries in the distances matrix.
        Use them to set names of nodes.

    Returns
    -------
    Node
        A root node of the neighbor joining tree.

    """
    D = distances.copy()
    nodes = {name: Node(name, None, 0, None, 0) for name in labels}

    while len(D) > 2:
        
        n = len(labels)
        
        R = np.sum(D, axis=1)

        Q = np.full_like(D, np.inf)#isto tabelo kot D samo da povsod inf
        for i in range(n):#naredimo Q
            for j in range(i + 1, n):
                Q[i, j] = (n - 2) * D[i, j] - R[i] - R[j]

        i, j = np.unravel_index(np.argmin(Q), Q.shape)#..unravel_index najde index
        if i > j:
            i, j = j, i

        #Distance from the pair members to the new node
        Li = 0.5 * D[i, j] + (R[i] - R[j]) / (2 * (n - 2))
        Lj = D[i, j] - Li

        #novo vozlisce
        new_name = f"({labels[i]},{labels[j]})"
        new_node = Node(
            name=new_name,
            left=nodes[labels[i]],
            left_distance=float(Li),
            right=nodes[labels[j]],
            right_distance=float(Lj),
        )

        # Distance of the other taxa from the new node
        new_row = []
        for k in range(n):
            if k not in (i, j):
                d_new = 0.5 * (D[i, k] + D[j, k] - D[i, j])
                new_row.append(d_new)

        #fukne stran i j  stolpec in vrstico
        D = np.delete(D, (i, j), axis=0)
        D = np.delete(D, (i, j), axis=1)
        #urinemo notr nove
        D = np.vstack((D, new_row))
        new_row.append(0)
        D = np.column_stack((D, new_row))

        new_label = new_name
        keep = [k for k in range(len(labels)) if k not in (i, j)]
        labels = [labels[k] for k in keep] + [new_label]
        
        nodes[new_label] = new_node

    i, j = 0, 1
    root_name = f"({labels[i]},{labels[j]})"
    root_distance = D[i, j] / 2
    root = Node(
        name=root_name,
        left=nodes[labels[i]],
        left_distance=float(root_distance),
        right=nodes[labels[j]],
        right_distance=float(root_distance),
    )

    return root

def plot_nj_tree(tree: Node, ax: Axes = None, **kwargs) -> None:
    """
    Plot a Neighbor-Joining phylogenetic tree horizontally.

    Parameters
    ----------
    tree : Node
        Root node of the tree.
    ax : Axes
        Matplotlib Axes object for plotting.
    kwargs :
        colors : dict (optional)
            Mapping from group name → color.
        taxonomies_groups : dict (optional)
            Mapping from group name → list of leaf names.

    Returns
    -------
    Axes
        The matplotlib Axes with the plotted tree.
    """
    colors = kwargs.get("colors", None)
    taxa_grupe = kwargs.get("taxa_grupe", None)
    scale = kwargs.get("scale", 1.0)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 15))

    # da lahko v visin ok zgleda
    y_positions = {}
    current_y = 0

    def assign_y(node):
        nonlocal current_y
        if node.left is None and node.right is None:#ko nimamo vec nic sam pusti isti y
            y_positions[node.name] = current_y
            current_y += 1
            return y_positions[node.name]
        left_y = assign_y(node.left)
        right_y = assign_y(node.right)
        y_positions[node.name] = (left_y + right_y) / 2 #sticisca damo umes
        return y_positions[node.name]

    assign_y(tree)

    def draw_node(node, x):
        name = node.name
        color = "black"


        if taxa_grupe and colors:
            for group, members in taxa_grupe.items():
                name1 = name.split("|")[0]
                if name1 in members:
                    color = colors[group]
                    #print('bl')
                    break

        # leaf node
        if node.left is None and node.right is None:
            ax.text(x + 0.2, y_positions[name], name,
                    va="center", fontsize=10, color=color)
            return

        # left branch
        left_x = x + node.left_distance * scale
        left_y = y_positions[node.left.name]
        ax.plot([x, left_x], [left_y, left_y], color="black")
        draw_node(node.left, left_x)

        # right branch
        right_x = x + node.right_distance * scale
        right_y = y_positions[node.right.name]
        ax.plot([x, right_x], [right_y, right_y], color="black")
        draw_node(node.right, right_x)

        # vertical connector
        ax.plot([x, x], [left_y, right_y], color="black")
        
        if hasattr(node, "confidence") and node.confidence is not None and (node.left or node.right):
            ax.text(x + 0.05, y_positions[name], f"{node.confidence:.2f}",
                fontsize=8, color="gray", va="bottom")
            
    draw_node(tree, 0)

 
    ax.set_ylim(-0.5, current_y - 0.5)
    ax.set_xlabel("Genetic distance")
    ax.set_ylabel("")
    ax.set_title("Neighbor-Joining Phylogenetic Tree")


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False, labelleft=False)


    x_min, x_max = ax.get_xlim()
    bar_len = (x_max - x_min) * 0.05
    ax.plot([0, bar_len], [-0.8, -0.8], color="black", lw=2)
    ax.text(bar_len / 2, -1.1, f"{bar_len:.1f}", ha="center", va="top", fontsize=9)

    if taxa_grupe and colors:
        handles = [plt.Line2D([0], [0], color=c, marker='o', lw=0, label=g)
                   for g, c in colors.items()]
        ax.legend(handles=handles, loc='upper left', fontsize=9, frameon=False)

    return ax


def _find_a_parent_to_node(tree: Node, node: Node) -> tuple:
    """Utility function for reroot_tree"""
    stack = [tree]

    while len(stack) > 0:

        current_node = stack.pop()
        if node.name == current_node.left.name:
            return current_node, "left"
        elif node.name == current_node.right.name:
            return current_node, "right"

        stack += [
            n for n in [current_node.left, current_node.right] if n.left is not None
        ]

    return None


def _remove_child_from_parent(parent_node: Node, child_location: str) -> None:
    """Utility function for reroot_tree"""
    setattr(parent_node, child_location, None)
    setattr(parent_node, f"{child_location}_distance", 0.0)


def reroot_tree(original_tree: Node, outgroup_node: Node) -> Node:

    #potreboval sem veliko pomoci is strani llm da sem to naredil

    """A function to create a new root and invert a tree accordingly.

    This function reroots tree with nodes in original format. If you
    added any other relational parameters to your nodes, these parameters
    will not be inverted! You can modify this implementation or create
    additional functions to fix them.

    Parameters
    ----------
    original_tree: Node
        A root node of the original tree.
    outgroup_node: Node
        A Node to set as an outgroup (already included in a tree).
        Find it by it's name and then use it as parameter.

    Returns
    -------
    Node
        Inverted tree with a new root node.
    """
    tree = deepcopy(original_tree)

    stars, lokacije = _find_a_parent_to_node(tree, outgroup_node)
    if lokacije == "left":
        razdalja = stars.left_distance
    else: 
        razdalja = stars.right_distance
        
    _remove_child_from_parent(stars, lokacije)

    nov_koren = Node("nov_koren", stars, razdalja / 2, outgroup_node, razdalja / 2)
    trenutno_vozlisce = stars

    while tree != trenutno_vozlisce:
        stars, lokacije = _find_a_parent_to_node(tree, trenutno_vozlisce)

        if lokacije == "left":
            razdalja = stars.left_distance
        else: 
            razdalja = stars.right_distance
            
        _remove_child_from_parent(stars, lokacije)

        prazna_stran = "left" if trenutno_vozlisce.left is None else "right"
        
        if prazna_stran == "left":
            trenutno_vozlisce.left_distance = razdalja
            trenutno_vozlisce.left = stars
        else:
            trenutno_vozlisce.right_distance = razdalja  
            trenutno_vozlisce.right = stars

        if tree.name == stars.name:
            break
        trenutno_vozlisce = stars

    if lokacije == "left":
        lokacije2 = "right"
        razdalja_drugega = stars.right_distance
        drugi_otrok = stars.right
    else:
        lokacije2 = "left" 
        razdalja_drugega = stars.left_distance
        drugi_otrok = stars.left
    
    if prazna_stran == "left":
        trenutno_vozlisce.left_distance = razdalja_drugega + razdalja
        trenutno_vozlisce.left = drugi_otrok
    else:
        trenutno_vozlisce.right_distance = razdalja_drugega + razdalja
        trenutno_vozlisce.right = drugi_otrok

    return nov_koren

def sort_children_by_leaves(tree: Node) -> None:
    """Sort the children of a tree by their corresponding number of leaves.

    The tree can be changed inplace.

    Paramteres
    ----------
    tree: Node
        The root node of the tree.

    """
    def pomoc(tree):
        if tree.left is None and tree.right is None:
            return 1
        else:
            return (pomoc(tree.left) if tree.left else 0) + (pomoc(tree.right) if tree.right else 0) 

    if tree.left and tree.right:
        #poglej kok ma listov
        l=pomoc(tree.left)
        r=pomoc(tree.right)
        if r < l:#levo un ka ma mn potomcev
            l=tree.left
            lp=tree.left_distance
            r=tree.right
            rp=tree.right_distance
            
            tree.left = r
            tree.right = l
            tree.left_distance = rp
            tree.right_distance = lp
        sort_children_by_leaves(tree.left)
        sort_children_by_leaves(tree.right)


def plot_nj_tree_radial(tree: Node, ax: Axes = None, **kwargs) -> None:
    """A function for plotting neighbor joining phylogeny dendrogram
    with a radial layout.

    Parameters
    ----------
    tree: Node
        The root of the phylogenetic tree produced by `neighbor_joining(...)`.
    ax: Axes
        A matplotlib Axes object which should be used for plotting.
    kwargs
        Feel free to replace/use these with any additional arguments you need.

    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>>
    >>> tree = neighbor_joining(distances)
    >>> fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    >>> plot_nj_tree_radial(tree=tree, ax=ax)
    >>> fig.savefig("example_radial.png")

    """

    #V veliki meri sem se pomagal s llm!
    
    colors = kwargs.get("colors", None)
    taxa_grupe = kwargs.get("taxa_grupe", None)
    scale = kwargs.get("scale", 1.0)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    
 
    linewidth = kwargs.get('lw', kwargs.get('linewidth', 1))
    barva_vej = kwargs.get('color', 'black')
    
 
    listi = {}
    koordinate = {}
    koti_klina = {}
    velikosti_klinov = {}
    
    def listi_fun(node):
        if node in listi:
            return listi[node]
            
        if node.left is None and node.right is None:
            listi[node] = 1
        else:
            levo_st = listi_fun(node.left) if node.left else 0
            desno_st = listi_fun(node.right) if node.right else 0
            listi[node] = levo_st + desno_st
        
        return listi[node]
    
    def dobi_dolzino_veje(node, stars):
        if stars.left == node:
            return stars.left_distance * scale
        elif stars.right == node:
            return stars.right_distance * scale
        return 0.0
    
    def dodeli_koordinate(node, stars=None, zacetni_kot=0.0, velikost_klina=2*np.pi):
        
        if stars is None:
            koordinate[node] = (0.0, 0.0)
            koti_klina[node] = 0.0
            velikosti_klinov[node] = 2 * np.pi
        else:
            dolzina_veje = dobi_dolzino_veje(node, stars)
            kot = zacetni_kot + velikost_klina / 2
            stars_x, stars_y = koordinate[stars]
            
            x = stars_x + dolzina_veje * np.cos(kot)
            y = stars_y + dolzina_veje * np.sin(kot)
            koordinate[node] = (x, y)
            

            ax.plot([stars_x, x], [stars_y, y], color=barva_vej, linewidth=linewidth)
            
            if node.left is None and node.right is None:
                barva_tocke = barva_vej
                if hasattr(node, 'name') and node.name and taxa_grupe and colors:
                    for skupina, clani in taxa_grupe.items():
                        osnovno_ime = node.name.split("|")[0]
                        if osnovno_ime in clani:
                            barva_tocke = colors[skupina]
                            break
                
                ax.plot(x, y, 'o', markersize=3, color=barva_tocke)
                
                if hasattr(node, 'name') and node.name:
                    kot_oznake = kot
                    x_oznake = x + 0.15 * np.cos(kot_oznake)
                    y_oznake = y + 0.15 * np.sin(kot_oznake)
                    
                    if np.pi/2 <= kot_oznake <= 3*np.pi/2:
                        ha = 'right'
                    else:
                        ha = 'left'
                    
                    barva_oznake = "black"
                    if taxa_grupe and colors:
                        for skupina, clani in taxa_grupe.items():
                            osnovno_ime = node.name.split("|")[0]
                            if osnovno_ime in clani:
                                barva_oznake = colors[skupina]
                                break
                    
                    ax.text(x_oznake, y_oznake, node.name, 
                           ha=ha, va='center', fontsize=7, color=barva_oznake)

        otroci = []
        if node.left is not None:
            otroci.append(node.left)
        if node.right is not None:
            otroci.append(node.right)
        
        if otroci:
            skupno_otroci_listi = sum(listi_fun(otrok) for otrok in otroci)
            trenutni_kot = zacetni_kot
            
            for otrok in otroci:
                velikost_klina_otroka = (listi_fun(otrok) / skupno_otroci_listi) * velikost_klina
                dodeli_koordinate(otrok, node, trenutni_kot, velikost_klina_otroka)
                trenutni_kot += velikost_klina_otroka
    

    skupno_listi = listi_fun(tree)
    

    dodeli_koordinate(tree)
    

    ax.set_aspect('equal')
    ax.axis('off')

    if koordinate:
        vse_x = [x for x, y in koordinate.values()]
        vse_y = [y for x, y in koordinate.values()]
        max_obseg = max(max(abs(x) for x in vse_x), max(abs(y) for y in vse_y)) * 1.3
        ax.set_xlim(-max_obseg, max_obseg)
        ax.set_ylim(-max_obseg, max_obseg)

    return ax



























    

    























