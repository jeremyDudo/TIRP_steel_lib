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
def gen_and_save(composition0, tests, element1={}, element2={}, temps={}):
    folder = subtitle(composition0) 
    if not os.path.exists(folder):
        os.mkdir(folder)
    folder += "/"
    filename = folder + file_name(composition0, element1=element1, element2=element2)

    if len(element1) == 0:
        data = single_TC_caller(tests, composition0, temps)
    else:
        comp_matr = compositions_matrix(composition0, element1, element2)
        data = matrix_TC_caller(tests, comp_matr, temps)

    if "printability" in tests:
        fr, csc, bcc, laves = data["FR"], data["CSC"], data["BCC_frac"], data["laves_frac"]
        save_mat_fr = np.asarray(fr)
        save_mat_csc = np.asarray(csc)
        save_mat_bcc = np.asarray(bcc)
        save_mat_laves = np.asarray(laves)
        np.savez_compressed(filename + "_printability", fr=save_mat_fr, csc=save_mat_csc, bcc=save_mat_bcc, laves=save_mat_laves)
            
    if "stable_del_ferrite" in tests:
        del_ferrite = data["del_ferrite"]
        save_mat_del_ferrite = np.asarray(del_ferrite)
        np.savez_compressed(filename + "_del_ferrite", del_ferrite=save_mat_del_ferrite) 

    if "ASP" in tests:
        asp = data["asp"]
        save_mat_asp = np.asarray(asp)
        np.savez_compressed(filename + "_asp", asp=save_mat_asp)

    if "phase_frac_and_APBE" in tests:
        gammaPrime, apbe = data["gammaPrime_mole_fraction"], data["APBEamount"]
        save_mat_gamma = np.asarray(gammaPrime)
        save_mat_apbe = np.asarray(apbe)
        np.savez_compressed(filename + "_phase_frac_and_APBE", gammaPrime=save_mat_gamma, apbe=save_mat_apbe)

    if "strength_and_DF" in tests:
        dG_diff, strength = data["dG_diff"], data["strength"]
        save_mat_dG_diff = np.asarray(dG_diff)
        save_mat_strength = np.asarray(strength)
        np.savez_compressed(filename+"_strength_and_DF", dG_diff=save_mat_dG_diff, strength=save_mat_strength)

# Load to use 
def load_and_use(composition0, tests, element1={}, element2={}, temps={}):
    folder = subtitle(composition0) 
    if not os.path.exists(folder):
        os.mkdir(folder)
    folder += "/"
    filename = folder + file_name(composition0, element1=element1, element2=element2)

    data = {} 
    if "printability" in tests:
        res = np.load(filename + "_printability.npz")
        fr, csc, bcc, laves = res["fr"], res["csc"], res["bcc"], res["laves"]
        data.update( {"fr": fr, "csc": csc, "bcc":bcc, "laves":laves})
    if "stable_del_ferrite" in tests:
        res = np.load(filename + "_del_ferrite.npz")
        del_ferrite = res["del_ferrite"]
        data.update( {"del_ferrite":del_ferrite} )
    if "ASP" in tests:
        res = np.load(filename + "_asp.npz")
        asp = res["asp"]
        data.update( {"asp" : asp} )
    if "phase_frac_and_APBE" in tests:
        res = np.load(filename + "_phase_frac_and_APBE.npz")
        gammaPrime, apbe = res["gammaPrime"], res["apbe"]
        data.update( {"gammaPrime" : gammaPrime, "apbe" : apbe} )
    if "strength_and_DF" in tests:
        res = np.load(filename + "_strength_and_DF.npz")
        dG_diff, strength = res["dG_diff"], res["strength"]
        data.update( {"dG_diff": dG_diff, "strength" : strength} )
    
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
    
    contour_font = 20
    if "FR" in tests: 
        fr_contours = ax.contour(X, Y, fr, colors='green', linestyles='dashdot')
        fr_fmt = '%1.0f'
        ax.clabel(fr_contours, inline=True, fontsize=contour_font, fmt=fr_fmt, manual=manual)
    
    if "CSC" in tests:
        csc_contours = ax.contour(X, Y, csc, colors='#cc5500', linestyles='dashed')
        csc_fmt = '%1.1f'
        ax.clabel(csc_contours, inline=True, fontsize=contour_font, fmt=csc_fmt, manual=manual)
    
    if "HCS" in tests:
        hcs_contours = ax.contour(X, Y, HCS(fr, csc), colors='k', linestyles='dashdot')
        hcs_fmt = '%1.2f'
        ax.clabel(hcs_contours, inline=True, fontsize=contour_font, fmt=hcs_fmt, manual=manual)        
    
    if "delta-Ferrite(at print)" in tests:
        bcc_contours = ax.contour(X, Y, bcc, colors='blue', linestyles='dashed')#, alpha=.85)
        bcc_fmt = '%1.2f'# '(AP)del-F: %1.1f'
        ax.clabel(bcc_contours, inline=True, fontsize=contour_font, fmt=bcc_fmt, manual=manual)        
    
    if "Laves" in tests:
        lav_contours = ax.contour(X, Y, laves, colors='red', linestyles='dashdot')
        lav_fmt = '%1.3f'# 'Fe2Ti_Laves: %1.1f'
        ax.clabel(lav_contours, inline=True, fontsize=contour_font, fmt=lav_fmt, manual=manual)        
    
    if "Stable del-Ferrite" in tests:
        meta_contours = ax.contour(X, Y, del_ferrite, colors='blue', linestyles='dashed', alpha=.45)
        meta_fmt = '%1.3f'
        ax.clabel(meta_contours, inline=True, fontsize=contour_font, fmt=meta_fmt, manual=manual)
            
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
    plt.savefig(folder+title+file_name(composition0, element1, element2, "printability")+".png")
    fig.show() 

# Runner 
def run():
    pass 

# Printer
def printer():
    pass 