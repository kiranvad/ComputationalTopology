'''
This is a helper function to compute representative cycles using dionysus.
This function has the following functions implemented:

    FindRepresentativeCycles   : Which determines representative cycles of any particular bar
    ComputeConnectedComponents : Computes connected components for each bar above a specified threshold
    from_hfp_to_connectedcomps : Computes connected components of a given bar
    
(c) Copyright Kiran Vaddi 01-2020

'''
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import networkx as nx
import time

__all__ = ["from_hfp_to_connectedcomps", "FindRepresentativeCycles", "ComputeConnectedComponents"]

class homcycles(object):
    
    def __init__(self, hfp, filtration, num_vertices, thresh=0.0):
        """
        A class for various operations using homology cycles:
        
        Parameters:
        -----------
                hfp : Homology peristenece class from dionysus2
                filtration : filtration from dionysus2
                num_vertices : number of 0-dimensional simplices in the data (typically number of data points)
                thresh : death threshold for computation
        """
        self.thresh = thresh
        self.hfp = hfp
        self.filtration = filtration
        self.num_vertices = num_vertices

    def FindRepresentativeCycles(self):
        '''
        This function can be used to collect representative cycles of the bars in a persistence
        diagram.
        This function takes a homology persistence class from dionysus as a input along with
        an optional threshold value in terms of the lifetime of the bar.
        The output is a numpy array of representative cycles. Each row in the array is one cycle per bar,
        with last three columns of each row reprsenting: dimension, birth and death of the cycle.
        Rest of the columns in the array are vertices in the cycle.

        All the dimensional cycles in the homology persistence are computed by default.
        '''
        
        repcycles = []
        ccverts = []
        for i,c in enumerate(self.hfp):
            if len(c)>0:
                birth = self.filtration[i].data 
                dim = self.filtration[i].dimension()-1
                deaths =[]
                vertices = []
                for x in c:
                    simplex = self.filtration[x.index]
                    deaths.append(self.filtration[x.index].data)
                    for s in simplex:
                        vertices.append(s)
                death = np.max(np.unique(np.asarray(deaths)))
                vertices = np.unique(np.asarray(vertices))
                vertices = np.append(vertices,[dim, birth, death])
                if birth-death>self.thresh:
                    repcycles.append(vertices)
        return repcycles
    
    def from_hfp_to_connectedcomps(self,index):
        """
        Computes connected components upto a death of a 0-cycle.
        Parameters:
        -----------
                hfp : Homology peristenece class from dionysus2
                filtration : filtration from dionysus2
                num_vertices : number of 0-dimensional simplices in the data (typically number of data points)
                index : filtration index of death of 0-cycle

        Returns:
        --------
        List of connected components as a set
        """
        total_iters = np.arange(self.num_vertices,index+1)
        adj = np.eye(self.num_vertices)
        for i,doi in enumerate(total_iters):
            current_vertices = list(self.filtration[doi])
            if doi==index:
                repvertex = current_vertices[1]
            adj[current_vertices[0],current_vertices[1]] = 1;
        adj = np.maximum( adj, adj.T)
        graph = csr_matrix(adj)
        n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)
        G = nx.from_scipy_sparse_matrix(graph)
        concomps = nx.node_connected_component(G, repvertex)
        # print('Vertex %d is connected to'%(repvertex),concomps)
        return concomps

    def ComputeConnectedComponents(self):
        """
        Computes connected components of bars above a threshold.
        Parameters:
        -----------
                hfp          : Homology peristenece class from dionysus2
                filtration   : filtration from dionysus2
                num_vertices : number of 0-dimensional simplices in the data (typically number of data points)
                thresh       : Threshold of death time of bars (default, 0.0)

        Returns:
        --------
        List of connected components as a dictornary, one for each bar. Dictornary has following keys:
                'cc'    : connected conponents as a list
                'birth' : birth of the bar
                'death' : death of the bar
                'index' : Index into the filtration where the bar is ended (birth of its representative 1-cycles)
        """
        all_concomps = []
        for i,c in enumerate(self.hfp):
            if len(c)>0:
                dim = self.filtration[i].dimension()-1
                if dim==0:
                    death = self.filtration[i].data 
                    births =[]
                    for x in c:
                        simplex = self.filtration[x.index]
                        births.append(self.filtration[x.index].data)
                    birth = np.max(np.unique(np.asarray(births)))
                    if death-birth>self.thresh:
                        print('Computing Bar %d which was killed at %0.3f'%(i,death))
                        start = time.time()
                        concomps={}
                        concomps['cc'] = self.from_hfp_to_connectedcomps(i)
                        concomps['birth'] = birth
                        concomps['death'] = death
                        concomps['index'] = i
                        all_concomps.append(concomps.copy())
                        end = time.time()
                        print('Elapsed Time is : ',end-start)
        return all_concomps