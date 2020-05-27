from plotting_lib import run 

dependent_element = "Fe"

# Fe is dependent element
# The values are wt. fraction of each element
composition = {
    "C": 0.0003,
    "Cr": .13,
    "Ni": .16,
    "Mo": .013,
    "Ti": 0.030,
    "V": 0.003,
    "Al": 0.002,
    #"B": 0.01/100
    "O": 0.015/100
} 

# define the two elements that define the range of values that are being varied
# wt fraction values
el1 = {
    "name" : "Cr",
    "start" : 13/100,
    "end" : 17/100,
    "length" : 2
}

el2 = {
    "name" : "Ni",
    "start" : 13/100,
    "end" : 17/100,
    "length" : 2
}

# in K
# just make sure that these are defined (but it's ok if your tests don't need explicitly need these)
temps = {
    "solution_temp" : 1000 + 273.15,
    "aging_temp" : 973.15,
    "start_temp" : 173.15,
    "end_temp" : 473.15
}

"""
Tests you can choose from: 

'fr' 
'csc'
'hcs' 
'meta_del_ferrite'
'laves' 
'stable_del_ferrite' 
'asp'                   # need nidata7.tdb database and MART5.TDB database
'gamma_prime'           # need nidata7.tdb database
'apbe'                  # need nidata7.tdb database
'dg_diff'               # need nidata7.tdb database
'strength'              # need nidata7.tdb database
'oxides' 

nidata7.tdb and mart5.tdb should be in the same folder as these files... 
I didn't include them in the github repo bc I am very unclear on the copyright of the database files
"""
tests = ["oxides"] # "stable_del_ferrite", "meta_del_ferrite", "hcs", ]

# If you need to overwrite exisiting data
overwrite = False 

# For final plots to allow manual contour labelling (so they don't overlap)
# The contours tend to be really ugly with where they put labels :/
manual = False

# for single point scheil calcs disp the scheil curve
disp = True 

# run your test(s, comment out as needed)!
# for a matrix of compositions:
# run(composition, tests, element1=el1, element2=el2, temps=temps, manual=manual, overwrite=overwrite)

# for a single point:
run(composition, tests, dependent_element, temps=temps, manual=manual, overwrite=overwrite, disp=disp)
