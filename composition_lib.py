import copy  
import numpy as np

# update composition 
def update_composition(comp0, ele1, ele2):
    """
    This is how I can generate the matrix with all different values
    """

    # deepcopy so that you don't have a single composition in the entire matrix
    # this allows comp to have the same values as comp0 but will not update comp0
    comp = copy.deepcopy(comp0)
    
    comp.update({ele1[0] : round(ele1[1],5)})
    comp.update({ele2[0] : round(ele2[1],5)})

    return comp
    
# generate composition 
def compositions_matrix(composition0, element1, element2):
    """
    This is how we take the varied elements and use them to generate a matrix with varied compositions
    """
    # take dicts of elements and generate linspaces
    x_vals = np.linspace(element1["start"], element1["end"], element1["length"])
    y_vals = np.linspace(element2["start"], element2["end"], element2["length"])
    
    # makes our matrix w/ a list comprehension (nested)
    composition_matr = [[update_composition(composition0, (element1["name"],x), (element2["name"], y)) for y in y_vals] for x in x_vals]
    return composition_matr