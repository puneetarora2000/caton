import numpy as np
import matplotlib.delaunay as sd
import networkx as nx

# make sure orders don't get messed up!!

class graph_info(object):
    def __init__(self,channel_graph,blowup):
        """channel_graph MUST be in terms of channels 0,...,N"""
        self.ChannelGraph = channel_graph
        
        N = channel_graph.number_of_nodes()
        self.IsComplete = (channel_graph.number_of_edges() == N*(N-1)/2)
        
        ChannelGraph2 = self.ChannelGraph.copy()
        for Node in ChannelGraph2.nodes_iter():
            ChannelGraph2.add_edge(Node,Node)

        self.Graph = blowup_graph(ChannelGraph2,blowup)
        for Node in self.Graph.nodes_iter():
            self.Graph.remove_edge(Node,Node)	

        self.Cliques,self.Seps = get_cliques_and_seps(self.Graph)
        self.Rank = max([len(clique) for clique in self.Cliques])+1
        self.NParams = self.Graph.number_of_edges()/2 + self.Graph.number_of_nodes()
        
    def concentration(self,S):
        if self.IsComplete:
            return np.linalg.inv(S)
        else:
            return sparse_conc(S,self.Cliques,self.Seps,self.Graph)
    
def restrict_to_graph(Arr,G):
    """
    G doesn't have an edge i->j (for i != j), then Out[i,j] = 0

    Parameters
    -------
    G : networkx Graph
    Arr : ndarray (2D square)

    Returns
    ------
    Out : ndarray (2D square)
    """    
    Restrictor  = np.identity(Arr.shape[0],dtype=np.float32)
    Restrictor  += np.array(nx.adj_matrix(G),dtype=np.float32)
    return Arr*Restrictor

def triangle_graph(Locs):
    """
    Returns a graph giving the Delaunay triangulation of a set of two-dimensional points    

    Parameters
    -------
    Locs : ndarray, or anything that gets cast into ndarray upon np.array(Locs)
        Locs.shape = (n,2), where n is the number of points


    Returns
    ------
    out : networkx Graph
    """
    Locs = np.array(Locs)
    if len(Locs)==1:
        return nx.complete_graph(1)
    else:
        Triangulation = sd.Triangulation(Locs[:,0],Locs[:,1])
        return nx.Graph(Triangulation.node_graph())

def path_graph(Locs):
    """
    Returns a line graph for a set of one-dimensional locations.    

    Parameters
    -------
    Locs : ndarray, or anything that gets cast into ndarray upon np.array(Locs)

    Returns
    ------
    out : networkx Graph
    """
    Locs = np.array(Locs).flatten()
    if len(Locs)==1:
        return nx.complete_graph(1)    
    else:
        Graph = nx.Graph()
        SortInds = np.argsort(Locs)
        for src,targ in zip(SortInds[1:],SortInds[:-1]):
            Graph.add_edge(src,targ)
        return Graph
                
    
def get_cliques_and_seps(Graph):
    """
    Gets cliques and separators for the input graph

    Parameters
    -------
    Graph : networkx Graph

    Returns
    ------
    Cliques : list of lists of node labels from Graph
    Seps : list of lists of node labels from Graph
    """

    ##### Construct Clique Graph and order the Cliques to satisfy running intersection property
    CliqueGraph = nx.clique.make_max_clique_graph(Graph)
    CliqueDic = dict([(node,clique) for (node,clique) in zip(CliqueGraph,nx.clique.find_cliques(Graph))])
    ### I used a dict because graph returned by make_max_clique_graph has nodes 1,2,..., which correspond to
    ### indices 0,1,... from list returned by find_cliques
    for src,targ in CliqueGraph.edges_iter():
        CliqueGraph[src][targ]['weight'] = -len(set(CliqueDic[src]).intersection(set(CliqueDic[targ])))
    # for src,trg in CliqueGraph.edges_iter():
        # assert cg.get_edge_data(src,trg) != 0
    ### made sure nodes of CliqueGraph correspond (in order) to Cliques
    if len(CliqueDic) == 1: return [CliqueDic[1]],[] # if graph has one node, mst has no edges so we get no graph
    
    MST = nx.Graph(nx.mst(CliqueGraph)) # Minimal weight spanning tree for CliqueGraph
    CliqueTree = nx.dfs_tree(MST) # MST turned into directed graph for traversal
    OrderedCliques = nx.search.dfs_preorder(CliqueTree) # Ordering of CliqueTree
   
    Seps = []     
    CliqueList = []
    
    CliqueList.append(CliqueDic[OrderedCliques[0]])
    for Node in OrderedCliques[1:]:
        ### First node is root of tree and has no predecessors. Other nodes of one predecessor each
        ### since this is a directed tree graph
        Seps.append(list(set(CliqueDic[Node]).intersection(set(CliqueDic[CliqueTree.predecessors(Node)[0]]))))        
        CliqueList.append(CliqueDic[Node])

    return CliqueList, Seps

def blowup_graph(G,scale):
    mat = nx.to_numpy_matrix(G)
    bigmat =  np.repeat(np.repeat(mat,scale,axis=0),scale,axis=1)
    return nx.from_numpy_matrix(bigmat)

def sparse_conc(S,Cliques,Seps,G):
    """
    Finds the concentration matrix for sample covariance matrix S, with the sparsity pattern given by G. 

    Parameters
    -------
    S : positive-definite D-by-D matrix.
    Cliques : list of lists of integers between 0 and D
    Seps : list of lists of integers between 0 and D
    G : networkx graph with nodes 0...D-1

    Returns
    ------
    Conc : positive-definite D-by-D matrix that has zeros
    wherever the original graph was missing an edge.
    """    
    Conc = np.zeros_like(S)
    for Clique in Cliques:
        Conc[np.ix_(Clique,Clique)] += np.linalg.inv(S[np.ix_(Clique,Clique)])
    for Sep in Seps:
        Conc[np.ix_(Sep,Sep)] -= np.linalg.inv(S[np.ix_(Sep,Sep)])
    return restrict_to_graph(Conc,G)