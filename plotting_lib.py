import matplotlib.pyplot as plt
import os
import numpy as np 
from composition_lib import compositions_matrix
from TC_lib import matrix_TC_caller, single_TC_caller

# for plotting/saving
def subtitle(composition0):
    dependent_element = "Fe"
    plotName = dependent_element
    for key, value in sorted(composition0.items(), key=lambda item: item[1], reverse=True):
        plotName = plotName+"-%s%s" % (round(value*100,3), key)
    return plotName

# file name generator
def file_name(composition0, element1={}, element2={}):
    if len(element1) == 0:
        return "point_calc"
    title = "Vary " + element1["name"] + "-" + element2["name"] + "__(" + str(element1["start"]) + "-" + str(element1["end"]) + "-" + str(element1["length"]) + ")(" + str(element2["start"]) + "-" + str(element2["end"]) + "-" + str(element2["length"]) + ")"
    return title

# Generate Data and Save
def gen_and_save(composition0, tests, element1={}, element2={}, temps={}, overwrite=False):
    folder = subtitle(composition0) 
    if not os.path.exists(folder):
        os.mkdir(folder)
    folder += "/"
    filename = folder + file_name(composition0, element1=element1, element2=element2)

    if "fr" or "csc" or "hcs" or "meta_del_ferrite" or "laves" in tests:
        tests.append("printability")

    if "gamma_prime" or "apbe" in tests:
        tests.append("phase_frac_and_apbe")

    if "dg_diff" or "strength" in tests:
        tests.append("strength_and_df") 

    # remove duplicates
    tests = list(dict.fromkeys(tests))

    """
    MOVE OVERWRITE FUNCTIONALITY HERE!!!!!!

    """


    if len(element1) == 0:
        data = single_TC_caller(tests, composition0, temps)
    else:
        comp_matr = compositions_matrix(composition0, element1, element2)
        data = matrix_TC_caller(tests, comp_matr, temps)

    if "printability" in tests:
        fr, csc, bcc, laves = data["fr"], data["csc"], data["BCC_frac"], data["laves_frac"]
        save_mat_fr = np.asarray(fr)
        save_mat_csc = np.asarray(csc)
        save_mat_bcc = np.asarray(bcc)
        save_mat_laves = np.asarray(laves)
        np.savez_compressed(filename + "_printability", fr=save_mat_fr, csc=save_mat_csc, bcc=save_mat_bcc, laves=save_mat_laves)
            
    if "stable_del_ferrite" in tests:
        del_ferrite = data["del_ferrite"]
        save_mat_del_ferrite = np.asarray(del_ferrite)
        np.savez_compressed(filename + "_del_ferrite", del_ferrite=save_mat_del_ferrite) 

    if "asp" in tests:
        asp = data["asp"]
        save_mat_asp = np.asarray(asp)
        np.savez_compressed(filename + "_asp", asp=save_mat_asp)

    if "phase_frac_and_apbe" in tests:
        gamma_prime, apbe = data["gamma_prime_mole_fraction"], data["apbe"]
        save_mat_gamma = np.asarray(gamma_prime)
        save_mat_apbe = np.asarray(apbe)
        np.savez_compressed(filename + "_phase_frac_and_apbe", gamma_prime=save_mat_gamma, apbe=save_mat_apbe)

    if "strength_and_df" in tests:
        dg_diff, strength = data["dg_diff"], data["strength"]
        save_mat_dg_diff = np.asarray(dg_diff)
        save_mat_strength = np.asarray(strength)
        np.savez_compressed(filename+"_strength_and_df", dg_diff=save_mat_dg_diff, strength=save_mat_strength)

# Load to use 
def load_and_use(composition0, tests, element1={}, element2={}, temps={}):
    folder = subtitle(composition0) 
    if not os.path.exists(folder):
        os.mkdir(folder)
    folder += "/"
    filename = folder + file_name(composition0, element1=element1, element2=element2)

    data = {} 
    if "printability" or "hcs" or "fr" or "csc" or "meta_del_ferrite" or "laves" in tests:
        res = np.load(filename + "_printability.npz")
        # fr, csc, bcc, laves = res["fr"], res["csc"], res["bcc"], res["laves"]
        data.update( dict(res) )# {"fr": fr, "csc": csc, "bcc":bcc, "laves":laves})
    if "stable_del_ferrite" in tests:
        res = np.load(filename + "_del_ferrite.npz")
        # del_ferrite = res["del_ferrite"]
        data.update( dict(res) )# {"del_ferrite":del_ferrite} )
    if "asp" in tests:
        res = np.load(filename + "_asp.npz")
        # asp = res["asp"]
        data.update( dict(res) )# {"asp" : asp} )
    if "phase_frac_and_apbe" or "gamma_prime" or "apbe" in tests:
        res = np.load(filename + "_phase_frac_and_apbe.npz")
        # gamma_prime, apbe = res["gamma_prime"], res["apbe"]
        data.update( dict(res) )# {"gamma_prime" : gamma_prime, "apbe" : apbe} )
    if "strength_and_df" or "dg_diff" or "strength" in tests:
        res = np.load(filename + "_strength_and_df.npz")
        # dg_diff, strength = res["dg_diff"], res["strength"]
        data.update( dict(res) )# {"dg_diff": dg_diff, "strength" : strength} )
    
    return data

# Plotter
def plotter(composition0, tests, element1, element2, temps, manual=False):
    subtitle_ = subtitle(composition0) 
    folder = subtitle + "/Plots" 
    if not os.path.exists(folder):
        os.mkdir(folder)
    folder += "/"
    x_vals = np.linspace(element1["start"]*100, element1["end"]*100, element1["length"])
    y_vals = np.linspace(element2["start"]*100, element2["end"]*100, element2["length"])
    
    X, Y = np.meshgrid(x_vals, y_vals) 
    
    data = load_and_use(composition0, tests, element1, element2, temps) 
    
    title = tests[0] 
    for i in tests[1:]: 
        title += " + " + i
    title += " Contours"
    
    figsize_ = 10
    fig, ax = plt.subplots(figsize=(figsize_, figsize_)) 
    
    # hcs
    def hcs(fr, csc):
        return fr * csc / 100 

    # update the tests list
    if "printability" in tests:
        tests.remove("printability")
        tests.append("fr")
        tests.append("csc")
        tests.append("hcs")
        tests.append("meta_del_ferrite")
        tests.append("laves")
    if "phase_frac_and_apbe" in tests:
        tests.remove("phase_frac_and_apbe")
        tests.append("gamma_prime")
        tests.append("apbe") 
    if "strength_and_df" in tests:
        tests.remove("strength_and_df")
        tests.append("dg_diff")
        tests.append("strength")        
    
    # remove duplicates
    tests = list(dict.fromkeys(tests))

    color_dict = { # need to have better colors yikes
        'fr' : 'green',
        'csc' : '#cc5500',
        'hcs' : 'k',
        'meta_del_ferrite' : 'blue',
        'laves' : 'red',
        'stable_del_ferrite' : 'blue',
        'asp' : 'red',
        'gamma_prime' : 'green',
        'apbe' : '#cc5500',
        'dg_diff' : 'red',
        'strength' : 'blue'
    }

    fmt_dict = { # check significant digits, these are hard-coded!
        'fr' : "fr: %1.0f",
        'csc' : 'csc: %1.1f',
        'hcs' : 'hcs: %1.2f',
        'meta_del_ferrite' : 'del-ferrite (AP): %1.2f',
        'laves' : 'laves: %1.3f',
        'stable_del_ferrite' : 'del-ferrite (PS): %1.3f',
        'asp' : 'asp: %1.0f',
        'gamma_prime' : 'gamma\': %1.0f',
        'apbe' : 'APBE: %1.0f',
        'dg_diff' : 'dg_diff: %1.0f',
        'strength' : 'strength: %1.0f'       
    }

    data_dict = { # this is a sad necessity, a weakness in my naming convention :(
        'fr' : "fr",
        'csc' : 'csc',
        'hcs' : 'hcs,
        'meta_del_ferrite' : 'bcc',
        'laves' : 'laves',
        'stable_del_ferrite' : 'del-ferrite',
        'asp' : 'asp',
        'gamma_prime' : 'gamma_prime,
        'apbe' : 'apbe',
        'dg_diff' : 'dg_diff',
        'strength' : 'strength' 
    }

    def contoured(test):
        contour_font = 20
        contour = ax.contour(X, Y, data[data_dict[test]], colors=color_dict[test], linestyles='dashed')
        fmt = fmt_dict[test]
        ax.clabel(contour, inline=True, fontsize=contour_font, fmt=fmt, manual=manual)

    for test in tests:
        contoured(test)

    title_font = 22
    subtitle_font = 18   
    ax_font = 22
    tick_font = 15
    fig.suptitle(title, y=0.945, fontsize=title_font)
    ax.set_title(subtitle_, fontsize=subtitle_font)
    ax.set_xlabel(element1["name"] + ' (wt %)', fontsize=ax_font)
    ax.set_ylabel(element2["name"] + ' (wt %)', fontsize=ax_font)
    ax.tick_params(axis='x',labelsize=tick_font)
    ax.tick_params(axis='y',labelsize=tick_font)
    ax.grid(axis='both',alpha=0.5)
    plt.savefig(folder+title+file_name(composition0, element1, element2)+".png")
    fig.show() 

# Runner 
def run(composition0, tests, element1={}, element2={}, temps={}, manual=False, overwrite=False):

    gen_and_save(composition0, tests, element1=element1, element2=element2, temps=temps, overwrite=overwrite) 

    # if single point, print out a nice table and save the nice table
    if len(element1) == 0:
        data = load_and_use(compositions_matrix, tests, temps=temps)

        """
        ADD PRETTY TABLE STUFF HERE
        """
        return None
    
    # else plot the stuff
    plotter(composition0, tests, element1, element2, temps, manual=manual)

 