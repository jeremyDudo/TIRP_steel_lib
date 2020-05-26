import matplotlib.pyplot as plt
import os
import numpy as np 
from composition_lib import compositions_matrix
from TC_lib import TC_caller
from prettytable import PrettyTable
import sys

# for plotting/saving
def subtitle(composition0):
    """
    Currently Hard-Coded for Fe systems!
    Makes the subtitle for the plots, which is also used to name the folder within "Results"

    Ex: Fe-16.0Ni-13.0Cr-3.0Ti-1.3Mo-0.3V-0.2Al-0.03C
    """
    dependent_element = "Fe"
    plotName = dependent_element
    for key, value in sorted(composition0.items(), key=lambda item: item[1], reverse=True):
        plotName = plotName+"-%s%s" % (round(value*100,3), key)
    return plotName

# file name generator
def file_name(composition0, element1={}, element2={}):
    """
    For single points: names npz files as point_calc_test_name.npz

    For matrix variations: names npz files as Vary Element1-Element2__(element1Start-element1Finish-element1Length)(element2Start-element2Finish-element2Length)_test_name.npz

    These are for readability at a glance (so user knows if the data has been generated for the test)
    """
    if len(element1) == 0:
        return "point_calc"
    title = "Vary " + element1["name"] + "-" + element2["name"] + "__(" + str(element1["start"]) + "-" + str(element1["end"]) + "-" + str(element1["length"]) + ")(" + str(element2["start"]) + "-" + str(element2["end"]) + "-" + str(element2["length"]) + ")"
    return title

# Generate Data and Save
def gen_and_save(composition0, testss, element1={}, element2={}, temps={}, overwrite=False, disp=False):
    """
    This is the big caller for this script
    Takes the variables defined in run.py and determines if single or matrix calcs

    This calls our TC_lib funcs, and then writes the data to .npz compressed files for later use
    """

    # Save our data in a results folder for cleanness
    folder = "Results"
    if not os.path.exists(folder):
        os.mkdir(folder)
    # Differentiate tests based on compositions
    subtitle_ = subtitle(composition0)
    folder += "/" + subtitle_
    if not os.path.exists(folder):
        os.mkdir(folder)
    folder += "/"

    # This is how we name our individual tests within the composition folder
    filename = folder + file_name(composition0, element1=element1, element2=element2)
    tests = [test for test in testss]

    # all of these tests are generated in the printability call in the TC_calc conditional matching "printability"
    if "fr" in tests or "csc" in tests or "hcs" in tests or "meta_del_ferrite" in tests or "laves" in tests:
        tests.append("printability")

    # both of theses tests are generated at once in "phase_frac_and_apbe"
    if "gamma_prime" in tests or "apbe" in tests:
        tests.append("phase_frac_and_apbe")

    # again, the TC_lib matches "strength_and_df", but the user wants "dg_diff" or "strength" so this is a compromise
    if "dg_diff" in tests or "strength" in tests:
        tests.append("strength_and_df") 

    # remove duplicates
    tests = list(dict.fromkeys(tests))

    # for each test, give filename a different extension for differentiability
    filename_dict = {
        "printability" : filename + "_printability",
        "stable_del_ferrite" : filename + "_del_ferrite",
        "asp" : filename + "_asp",
        "phase_frac_and_apbe" : filename + "_phase_frac_and_apbe",
        "strength_and_df" : filename+"_strength_and_df"
    }

    # sorry for Ugly comprehension, but need to fliter list down to only those being called in the TC_caller & filename_dict for the loop to work
    # we want only the tests that TC_lib cares about so the loop is faster
    tests = [k for k in tests if 'printability' in k or 'stable_del_ferrite' in k or 'asp' in k or 'phase_frac_and_apbe' in k or 'strength_and_df' in k]

    # removing elements from the list you loop over kept crashing so had to do it in two steps
    remover = []
    for test in tests:
        # check if file has been saved and if we don't want to overwrite
        if os.path.exists(filename_dict[test]+".npz") and not overwrite: 
            # if we have data we don't need to overwrite, don't run the test again
            remover.append(test)
    # step 2 to remove elements that we don't need to (re)generate
    tests = [test for test in tests if test not in remover]

    # check if a singlepoint by seeing if we are testing ranges of elements
    if len(element1) == 0:
        data = TC_caller(tests, composition0, temps, disp=disp)
    else:
        comp_matr = compositions_matrix(composition0, element1, element2)
        data = TC_caller(tests, comp_matr, temps)

    # pull all data out of generated data and put into the save file with nice dict keys!
    if "printability" in tests:
        fr, csc, bcc, laves = data["fr"], data["csc"], data["BCC_frac"], data["laves_frac"]
        save_mat_fr = np.asarray(fr)
        save_mat_csc = np.asarray(csc)
        save_mat_bcc = np.asarray(bcc)
        save_mat_laves = np.asarray(laves)
        np.savez_compressed(filename + "_printability", fr=save_mat_fr, csc=save_mat_csc, meta_del_ferrite=save_mat_bcc, laves=save_mat_laves)
            
    if "stable_del_ferrite" in tests:
        del_ferrite = data["del_ferrite"]
        save_mat_del_ferrite = np.asarray(del_ferrite)
        np.savez_compressed(filename + "_del_ferrite", stable_del_ferrite=save_mat_del_ferrite) 

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
    """
    Effectively the other half of "gen_and_save"
    Takes same arguments minus "overwrite" b/c you aren't going to overwrite anything while loading

    This just locates the files saved by "gen_and_save", loads it, and updates it all into one dictionary with the same keys
    """
    # find file that was saved
    folder = "Results"
    subtitle_ = subtitle(composition0)
    folder += "/" + subtitle_ 
    folder += "/"
    filename = folder + file_name(composition0, element1=element1, element2=element2)

    # commented out lines show what handles & info is being pushed into the data dict
    # Please do not remove, I think it is helpful to see what is being moved where
    data = {} 

    if "printability" in tests or "hcs" in tests or "fr" in tests or "csc" in tests or "meta_del_ferrite" in tests or "laves" in tests:
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
    if "phase_frac_and_apbe" in tests or "gamma_prime" in tests or "apbe" in tests:
        res = np.load(filename + "_phase_frac_and_apbe.npz")
        # gamma_prime, apbe = res["gamma_prime"], res["apbe"]
        data.update( dict(res) )# {"gamma_prime" : gamma_prime, "apbe" : apbe} )
    if "strength_and_df" in tests or "dg_diff" in tests or "strength" in tests:
        res = np.load(filename + "_strength_and_df.npz")
        # dg_diff, strength = res["dg_diff"], res["strength"]
        data.update( dict(res) )# {"dg_diff": dg_diff, "strength" : strength} )
    
    return data

# Plotter
def plotter(composition0, tests, element1, element2, temps, manual=False):
    """
    Calls our data generator/loader
    Unpacks all of our data
    Plots our data
    """

    # We want the plots to go in the same folder as the rest of the data, but in a subfolder "plots"
    # The conditionals are absolutely redundant at this point, might be deleted for speed?
    folder = "Results"
    if not os.path.exists(folder):
        os.mkdir(folder)
    subtitle_ = subtitle(composition0)
    folder += "/" + subtitle_
    if not os.path.exists(folder):
        os.mkdir(folder)
    folder += "/Plots" 
    if not os.path.exists(folder):
        os.mkdir(folder)
    folder += "/"
    
    # Convert our element dictionaries into np.linspaces 
    x_vals = np.linspace(element1["start"]*100, element1["end"]*100, element1["length"])
    y_vals = np.linspace(element2["start"]*100, element2["end"]*100, element2["length"])
    
    # for contour plotting
    X, Y = np.meshgrid(x_vals, y_vals) 
    
    # load the data
    data = load_and_use(composition0, tests, element1, element2, temps) 
    
    # title is to help show what tests are run
    title = tests[0] 
    for i in tests[1:]: 
        title += " + " + i
    title += " contours"

    # open a plot
    figsize_ = 10
    fig, ax = plt.subplots(figsize=(figsize_, figsize_)) 
    
    # hcs
    def hcs(fr, csc):
        return fr * csc / 100 

    # update the tests list
    # the big names are the names of the functions that calculate more than a single thing --
    # we don't want redundencies in the plotting, and we want the tests to only have the names 
    # of plots that we can add
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

    # for contours, so we can have a single function
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

    # for contours, so we can have a single function
    # the lables are a maybe, Greg doesn't like them
    fmt_dict = { # check significant digits, these are hard-coded!
        'fr' : "fr: %1.0f",
        'csc' : 'csc: %1.1f',
        'hcs' : 'hcs: %1.2f',
        'meta_del_ferrite' : 'del_ferrite (AP): %1.2f',
        'laves' : 'laves: %1.3f',
        'stable_del_ferrite' : 'del_ferrite (PS): %1.3f',
        'asp' : 'asp: %1.0f',
        'gamma_prime' : 'gamma\': %1.0f',
        'apbe' : 'APBE: %1.0f',
        'dg_diff' : 'dg_diff: %1.0f',
        'strength' : 'strength: %1.0f'       
    }

    # function-ize our plotting so that it is Easy
    def contoured(test):

        # legible font
        contour_font = 20

        # hcs is not a listed calc, but is a function of two listed calcs, check for it (this is definitely slow and should likely be fixed in the TC_lib)
        if test == "hcs":
            Z = hcs(data['fr'], data['csc'])
        else:
            Z = data[test]
        
        # actually generate and plot the contour
        contour = ax.contour(X, Y, Z, colors=color_dict[test], linestyles='dashed')
        fmt = fmt_dict[test]
        ax.clabel(contour, inline=True, fontsize=contour_font, fmt=fmt, manual=manual)

    # plot all of the contours
    for test in tests:
        contoured(test)

    # Format the rest of the plot

    # legible font sizes for a ppt
    title_font = 22
    subtitle_font = 18   
    ax_font = 22
    tick_font = 15
    # supertitle with the tests 
    fig.suptitle(title, y=0.945, fontsize=title_font)
    # subtitle with the composition
    ax.set_title(subtitle_, fontsize=subtitle_font)
    # name lables with our varied elements
    ax.set_xlabel(element1["name"] + ' (wt %)', fontsize=ax_font)
    ax.set_ylabel(element2["name"] + ' (wt %)', fontsize=ax_font)
    # set size of ticks
    ax.tick_params(axis='x',labelsize=tick_font)
    ax.tick_params(axis='y',labelsize=tick_font)
    # pale grid for readability
    ax.grid(axis='both',alpha=0.5)
    # save
    plt.savefig(folder+file_name(composition0, element1, element2)+title+".png")
    fig.show() 

# Runner 
def run(composition0, tests, element1={}, element2={}, temps={}, manual=False, overwrite=False, disp=False):
    """
    To be called in run.py
    Takes parameters defined in run.py
    Runs the tests you want!

    disp is for single point scheil calculations if you want to see a Scheil Curve
        it is disabled for matrix viewing as you will have way too many sheil curves
    """
    # make sure data is generated & saved
    gen_and_save(composition0, tests, element1=element1, element2=element2, temps=temps, overwrite=overwrite, disp=disp) 

    # if single point, print out a nice table and save the nice table
    if len(element1) == 0:
        # load the data
        data = load_and_use(composition0, tests, temps=temps)
        
        # again hcs is a function of two named tests, so ugly conditionals
        def hcs(fr, csc):
            return fr*csc/100
        if "hcs" in tests:
            data.update( {"hcs": hcs(data['fr'], data['csc'])} )

        # make a pretty table
        table = PrettyTable() 

        # headers of the columns are the names of each test
        field_names = tests

        # format the column headers nicely
        field_names.sort() 
        data_dict = { 
            'fr' : "FR",
            'csc' : 'CSC',
            'hcs' : 'HCS',
            'meta_del_ferrite' : 'Del-Ferrite (At Print)',
            'laves' : 'Laves',
            'stable_del_ferrite' : 'Del-Ferrite (Post Solutionizing)',
            'asp' : 'ASP',
            'gamma_prime' : 'Gamma\'',
            'apbe' : 'APBE',
            'dg_diff' : 'dG_diff',
            'strength' : 'Strength' 
        }
        field_names = [data_dict[name] for name in field_names]
        table.field_names = field_names

        # row with the numerical data
        row = [data[key] for key in field_names] 
        table.add_row(row)

        # just calling the test our single_point file
        title_ = subtitle(composition0)
        print(table.get_string(title=title_))

        # save the table to a text file
        # HEY: THIS FILE IS ALWAYS OVERWRITTEN 
        table_txt = table.get_string()
        folder = "Results"
        if not os.path.exists(folder):
            os.mkdir(folder)
        folder += "/" + title_
        if not os.path.exists(folder):
            os.mkdir(folder)
        with open(folder+ '/single_point.txt','w') as file:
            file.write(table_txt)
        return None
    
    # else plot the stuff
    plotter(composition0, tests, element1, element2, temps, manual=manual)
    return None
 