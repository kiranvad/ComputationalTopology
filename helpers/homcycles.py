'''
This is a helper function to compute representative cycles using dionysus.
This function has the following functions implemented:

    findConnectedComponents : Which determines vertices in any particular bar of zero-dimension
    FindRepresentativeCycles: Which determines representative cycles of any particular bar
    
(c) Copyright Kiran Vaddi 01-2020

'''
import dionysus as dion
import numpy as np

def findConnectedComponents(start_index):
    '''
    A helper function to compute connected component vertices.
    This function takes index of any filtration where a 0-cycle death is detected and
    returns all the vetices that are in the bar corresponding to the 0-cycle that was dead.
    '''
    nodes = list(filtration[start_index])
    pairs = []
    vertices = []
    while not any(p == hfp.unpaired for p in pairs):
        nodes_temp = []
        vertices.append(nodes)
        for n in nodes:
            # print(n, hfp.pair(n))
            pairs.append(hfp.pair(n))
            if hfp.pair(n) != hfp.unpaired:
                nodes_temp.append(list(filtration[hfp.pair(n)]))
        nodes_temp = np.unique(np.asarray(nodes_temp))
        nodes = np.setdiff1d(nodes_temp, nodes).tolist()
    vertices = np.unique(np.hstack(vertices))
    # print(vertices)
    return vertices

def FindRepresentativeCycles(HomPer, thresh=0.0):
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
    for i,c in enumerate(hfp):
        if len(c)>0:
            birth = filtration[i].data 
            dim = filtration[i].dimension()-1
            deaths =[]
            vertices = []
            for x in c:
                simplex = filtration[x.index]
                deaths.append(filtration[x.index].data)
                for s in simplex:
                    vertices.append(s)
            death = np.max(np.unique(np.asarray(deaths)))
            vertices = np.unique(np.asarray(vertices))
            vertices = np.append(vertices,[dim, birth, death])
            if birth-death>thresh:
                repcycles.append(vertices)
                if dim==0:
                    ccverts.append(findConnectedComponents(i))                
    return repcycles, ccverts
    