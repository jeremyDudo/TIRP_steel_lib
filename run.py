from plotting_lib import run 

composition = {
    "C": 0.0003,
    "Cr": .12,
    "Ni": .17,
    "Mo": .013,
    "Ti": 0.030,
    "V": 0.003,
    "Al": 0.002
    #"B": 0.01/100
    #"O": 0.015/100
} 

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
temps = {
    "solution_temp" : 1000 + 273.15,
    "aging_temp" : 973.15
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

nidata7.tdb and mart5.tdb should be in the same folder as these files... 
I didn't include them in the github repo bc I am very unclear on the copyright of the database files
"""
tests = ["stable_del_ferrite", "meta_del_ferrite", "hcs"]

# For if bugs are found in the code after generating data
overwrite = False 

# For final plots to allow manual contour labelling (so they don't overlap)
manual = False

# run your test!
# for a matrix of compositions:
run(composition, tests, element1=el1, element2=el2, temps=temps, manual=manual, overwrite=overwrite)

# for a single point:
# run(composition, tests, temps=temps, manual=manual, overwrite=overwrite)
