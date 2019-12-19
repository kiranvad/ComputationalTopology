"""
A set of helper functions for computational topology repository
"""

import numpy as np
import matplotlib.pyplot as plt




def plot_a_bar(p, q, c='b', linestyle='-'):
    plt.plot([p[0], q[0]], [p[1], q[1]], c=c, linestyle=linestyle, linewidth=1)

def plot_barcode(diagrams, dimensions=all):
    """
    # An example usage:
    from ripser import ripser
    from sklearn import datasets
    data = datasets.make_circles(n_samples=100)[0] + 5 * datasets.make_circles(n_samples=100)[0]
    dgms = ripser(data)['dgms']
    plot_barcode(dgms,dimensions=0)    

    # This is adapted from : https://github.com/tmelorc/persim/blob/master/persim/visuals.py
    """
    if not isinstance(diagrams, list):
    # Must have diagrams as a list for processing downstream
        diagrams = [diagrams]
    if isinstance(dimensions,int):
        dimensions = dimensions+1
    else: 
        dimensions = len(diagrams)
        
    for dim in range(dimensions):
        number_of_bars = len(diagrams[dim])
        
        number_of_bars_fin = 0
        number_of_bars_inf = 0 

        if number_of_bars > 0:
            ax  = plt.gca()
            for i in range(number_of_bars):
                birth = [diagrams[dim][i, 0], i]
                death = [diagrams[dim][i, 1], i]
                maximum_death = np.nanmax(diagrams[dim][diagrams[dim] < 1E308].flatten())
                if np.isinf(death[0]):
                    number_of_bars_inf += 1
                    plot_a_bar(birth, [1.05 * maximum_death, i])
                    plt.scatter([1.05 * maximum_death], [i], c='b', s=10, marker='>')
                else:
                    number_of_bars_fin += 1
                    plot_a_bar(birth, death)
            title = "%d-dimensional bars: %d finite, %d infinite" % (dim, number_of_bars_fin, number_of_bars_inf)
            ax.set_title(title, fontsize=10)
            plt.yticks([])
            plt.show()
    
    
    
    
    
    
    