import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tc_python import *
import os
import copy 

# Generate Printability Values from Scheil Calculations
def printability(composition, scheil_calculation, system, disp=False):
    """
    Parameters
    ----------
    composition : Dict
        Ex: c_1217 = {
                       "C": 0.03/100,
                       "Cr": 12/100,
                       "Ni": 17/100,
                       "Mo": 1.3/100,
                       "Ti": 3.0/100,
                       "V": 0.3/100,
                       "Al": 0.2/100
                       #"B": 0.01/100
                       #"O": 0.015/100
                       }  # in wt-fraction 
        
    scheil_calculation : TCPython session for scheil_calculations (if taken as an input, we can evaluate over a range of values with only one instancing call)
        TCPython().set_cache_folder(os.path.basename(__file__) + "_cache")
              .select_database_and_elements(database, [dependent_element] + elements)
              .get_system_for_scheil_calculations().with_scheil_calculation()
              .set_composition_unit(CompositionUnit.MASS_FRACTION)
        
    disp : Bool, optional
        Shows intermediate plots (Would only recommend for a single composition value).
        The default is False.

    Returns
    -------
    [fr, csc, BCC_frac, laves_frac]
        fr : freezing Range [Float]
        csc : Cracking Susceptibility Coefficient [Float] 
        BCC_frac : Mole fraction of delta-Ferrite present at-print [Float]
        laves_frac : Mole fraction of laves present at-print [Float]
    """
    for element in composition:
        scheil_calculation = scheil_calculation.set_composition(element, composition[element])

    solidification = scheil_calculation.calculate()

    scheil_curve = solidification.get_values_grouped_by_stable_phases_of(
            ScheilQuantity.mole_fraction_of_all_solid_phases(),
            ScheilQuantity.temperature())

    if disp: 
        # Create plot title of composition, e.g., Fe-10Ni-7Cr-0.2C
        dependent_element = "Fe"
        plotName = dependent_element
        for key, value in sorted(composition.items(), key=lambda item: item[1], reverse=True):
            plotName = plotName+"-%s%s" % (value, key)
        # 1. Plot the solidification curve (mole fraction solid phases vs. T) including the equilibrium
        temp_min = 1e6
        temp_max = -1e6
        fig, ax_1 = plt.subplots(figsize=(10,10))
        

        for label in scheil_curve:
            section = scheil_curve[label]
            temp_min = min(np.min(section.y), temp_min)
            temp_max = max(np.max(section.y), temp_max)
            ax_1.plot(section.x, np.array(section.y) - 273.15, label=label)
            
    
        # calculate the equilibrium solidification line (starting at the liquidus temperature)
        prop_calculation = system.with_property_diagram_calculation()
        prop = (prop_calculation.
                  with_axis(CalculationAxis(ThermodynamicQuantity.temperature()).
                            set_min(temp_min).
                            set_max(temp_max)))
        for element in composition:
            prop = prop.set_condition("w("+element+")", composition[element])
    
        result = prop.calculate()
    
        temp_eq_frac, liquid_eq_frac = result.get_values_of(ThermodynamicQuantity.temperature(),
                                                            ThermodynamicQuantity.mole_fraction_of_a_phase("LIQUID"))
        solid_eq_frac = 1 - np.array(liquid_eq_frac)
        temp_eq_frac = np.array(temp_eq_frac) - 273.15
        valid_indices = np.where(solid_eq_frac < 1.0)  # cutting off all values with 100% solid

        title_font = 25
        subtitle_font = 18 
        ax_font = 25
        tick_font = 17
        ax_1.tick_params(axis='x',labelsize=tick_font)
        ax_1.tick_params(axis='y',labelsize=tick_font)
        ax_1.plot(solid_eq_frac[valid_indices], temp_eq_frac[valid_indices], '--', label="Equilibrium")
        ax_1.set_xlabel("Mole fraction of all solid phases [-]", fontsize= ax_font)
        ax_1.set_ylabel("Temperature [\N{DEGREE SIGN} C]", fontsize= ax_font)
    
        ax_1.legend(loc="lower left",fontsize=tick_font) # , prop={'size': 6}) #, bbox_to_anchor=(1,1))
        fig.suptitle("Scheil Calculation", y=0.945, fontsize=title_font)
        ax_1.set_title(plotName, fontsize=subtitle_font)
    
        # 2. Plot the mole fraction of the solid phases separated for each
        """
        groups = \
            solidification.get_values_grouped_by_stable_phases_of(ScheilQuantity.mole_fraction_of_all_solid_phases(),
                                                                  ScheilQuantity.mole_fraction_of_a_solid_phase(ALL_PHASES))
    
        for group in groups.values():
            ax_2.plot(group.x,np.array(group.y),label=group.label)
    
        ax_2.set_xlabel("Mole fraction of all solid phases [-]")
        ax_2.set_ylabel("Mole fraction of each solid phase [-]")
        ax_2.legend(loc="upper left", bbox_to_anchor=(1,1))
        """
        
        plt.show()
        sys.exit()
    
    del_ferrite = False
    laves = False
    for label in scheil_curve:        
        if "BCC_A2" in label:
            del_ferrite = True
        
        if "C14_LAVES" in label:
            laves = True
    print(composition)
    NS_T = solidification.get_values_of(
        ScheilQuantity.mole_fraction_of_all_solid_phases(),
        ScheilQuantity.temperature())
    print("Got NS_T")
    NS_NH = solidification.get_values_of(
        ScheilQuantity.mole_fraction_of_all_solid_phases(),
        ScheilQuantity.heat_per_mole())
    print("Got NS_NH")
    # check if labels contains BCC_A2 before generating this
    if del_ferrite:
        BCC_frac0 = solidification.get_values_of(
            ScheilQuantity.mole_fraction_of_all_solid_phases(),
            ScheilQuantity.mole_fraction_of_a_solid_phase("BCC_A2"))
        BCC_frac = BCC_frac0[1][-1]
        
    else:
        BCC_frac = 0
    print("Got BCC_frac")
    if laves:
        laves_frac0 = solidification.get_values_of(
            ScheilQuantity.mole_fraction_of_all_solid_phases(),
            ScheilQuantity.mole_fraction_of_a_solid_phase("C14_LAVES"))
        laves_frac = laves_frac0[1][-1]
    else:
        laves_frac = 0
    print("Got laves_frac")

    NS = NS_NH[0]
    NH = NS_NH[1]
    T = NS_T[1]
    
    fr = round(max(T) - min(T), 2)

    # Calculate Cracking Susceptibility Coeffieicent 
    def csc(NS, NH):
        """
        Parameters
        ----------
        NS : [Float]
            scheil_calculation.calculate().get_values_of(
                ScheilQuantity.mole_fraction_of_all_solid_phases(),
                ScheilQuantity.heat_per_mole()
            )[0]
            
        NH : [Float]
            scheil_calculation.calculate().get_values_of(
                ScheilQuantity.mole_fraction_of_all_solid_phases(),
                ScheilQuantity.heat_per_mole()
            )[1]

        Returns
        -------
        csc : [Float]
            Cracking Susceptibility Coefficient (csc) as defined by QuesTek.
        """
        
        NS = np.asarray(NS)
        NH = np.asarray(NH)
        t = []

        for i, e in enumerate(NH):
            t.append((e / NH[-1]) ** 2)

        df1 = (np.abs(NS - 0.4)).argmin()
        df2 = (np.abs(NS - 0.9)).argmin()
        df3 = (np.abs(NS - 0.99)).argmin()

        t1 = t[df1]
        t2 = t[df2]
        t3 = t[df3]

        t_v = t3 - t2
        t_r = t2 - t1

        csc = round((t_v / t_r), 5)

        return csc

    csc = csc(NS, NH)
    
    if disp:
        print("freezing Range = %s" % fr)
        print("Cracking susceptibility coefficient (csc) = %s" % csc)
    
    return [fr, csc, BCC_frac, laves_frac]

# Calculate how much stable del-ferrite present post solution treatment
def stable_del_ferrite(composition, solution_temp, calculation):
    """
    Parameters
    ----------
    composition : Dict
        Ex: c_1217 = {
                       "C": 0.03/100,
                       "Cr": 12/100,
                       "Ni": 17/100,
                       "Mo": 1.3/100,
                       "Ti": 3.0/100,
                       "V": 0.3/100,
                       "Al": 0.2/100
                       #"B": 0.01/100
                       #"O": 0.015/100
                       }  # in wt-fraction #
        
    solution_temp : [Float]
        Solution Temperature. For 2020 3.041 TRIP Steel we look for 1000-1200 Celsius
        !!! Value must be in Kelvin !!!
        
    calculation : TCPython session for single equilibrium calculation
        TCPython()..set_cache_folder(os.path.basename(__file__) + "_cache")
              .select_database_and_elements(database, [dependent_element] + elements)
              .get_system()
              .with_single_equilibrium_calculation()

    Returns
    -------
    del_ferrite : [Float]
        Amount of delta ferrite in matrix post solution treatment
    """

    calculation.set_condition(ThermodynamicQuantity.temperature(), solution_temp)
    
    for element in composition:
        calculation.set_condition("w("+element+")", composition[element])
    point_calculation = calculation.calculate()

    phases = point_calculation.get_stable_phases()

    del_ferrite = 0
    if "BCC_A2#1" in phases:
        del_ferrite += point_calculation.get_value_of(ThermodynamicQuantity.mole_fraction_of_a_phase('BCC_A2#1')) 
    if "BCC_B2#1" in phases:
        del_ferrite += point_calculation.get_value_of(ThermodynamicQuantity.mole_fraction_of_a_phase('BCC_B2#1'))
        
    return del_ferrite

# Calculate the Austenite Stability Parameter
def asp(composition, singlepoint, calculation):

    matrixComposition = composition 
    # matrixComposition = copy.deepcopy(composition)
    # gammaPrimeComposition = copy.deepcopy(composition)

    YS_RT = 1035

    # frictional Work
    def w_fSol(comp, t, warmWork):
        # Constants
        # Athermal work term coeffecients [J/mol]
        k_ath = {
            'c': 3807,
            'n': 3048,
            'cr': 1868,
            'mo': 1418,
            'ti': 1473,
            'v': 1618,
            'al': 280,
            'b': 0,
            'ni': 172,
            'mn': 1980,
            'si': 1879
        }

        # Thermal work term coeffecients [J/mol]
        k_th = {
            'c': 21216,
            'n': 16986,
            'cr': 3923,
            'mo': 2918,
            'ti': 3031,
            'v': 3330,
            'al': 576,
            'b': 0,
            'ni': 345,
            'mn': 4107,
            'si': 3867
        }

        w_0_Fe = 836
        A = 400
        epsilon = np.log((1 / (1 - warmWork)) ** (1 / 2))
        t_mu = 510  # kelvin

        ####################

        # Calculate frictional work term for a given composition

        # w_Mu = athermal work term
        w_Mu_i = np.sqrt(((k_ath['c'] ** 2) * comp['c']) +
                        ((k_ath['n'] ** 2) * comp['n']))

        w_Mu_j = np.sqrt(((k_ath['cr'] ** 2) * comp['cr']) +
                        ((k_ath['mn'] ** 2) * comp['mn']) +
                        ((k_ath['si'] ** 2) * comp['si']) +
                        ((k_ath['mo'] ** 2) * comp['mo']) +
                        ((k_ath['ti'] ** 2) * comp['ti']) +
                        ((k_ath['cr'] ** 2) * comp['cr']))

        w_Mu_k = np.sqrt(((k_ath['al'] ** 2) * comp['al']) +
                        ((k_ath['b'] ** 2) * comp['b']) +
                        ((k_ath['ni'] ** 2) * comp['ni']))

        w_Mu = w_Mu_i + w_Mu_j + w_Mu_k

        # w_th = thermal work term
        w_0_i = np.sqrt(((k_th['c'] * (comp['c']) ** 0.5) ** 2) +
                        ((k_th['n'] * (comp['n']) ** 0.5) ** 2))

        w_0_j = np.sqrt(((k_th['cr'] * (comp['cr']) ** 0.5) ** 2) +
                        ((k_th['mn'] * (comp['mn']) ** 0.5) ** 2) +
                        ((k_th['si'] * (comp['si']) ** 0.5) ** 2) +
                        ((k_th['mo'] * (comp['mo']) ** 0.5) ** 2) +
                        ((k_th['ti'] * (comp['ti']) ** 0.5) ** 2) +
                        ((k_th['cr'] * (comp['cr']) ** 0.5) ** 2))

        w_0_k = np.sqrt(((k_th['al'] * (comp['al']) ** 0.5) ** 2) +
                        ((k_th['ni'] * (comp['ni']) ** 0.5) ** 2) +
                        ((k_th['b'] * (comp['b']) ** 0.5) ** 2))

        w_0 = w_0_Fe + w_0_i + w_0_j + w_0_k

        w_th = w_0 * (1 - (t / t_mu) ** (2 / 3)) ** 2

        # w_Rho = forest hardening contribution
        w_Rho = A * epsilon ** (1 / 2)

        # frictional work = w_FSol
        w_fSol = w_Mu + w_th + w_Rho

        return float(w_fSol)   
    
    # Uniform Tension (?)
    def dGsigUT(sigma, Vdil):
        GsigUT = -(0.7183*sigma+6.85*Vdil*(sigma/3)-185.3*(1-np.exp(-0.003043*sigma)))
        return GsigUT

    # Crack Tip (?)
    def dGsigCT(sigma, Vdil):
        GsigCT = -(0.7183*sigma+6.85*Vdil*(2.67*sigma)-185.3*(1-np.exp(-0.003043*sigma)))
        return GsigCT

    for element in composition:
        singlepoint.set_condition("w("+element+")") 
    point_calculation = singlepoint.calculate() 

    # phases = point_calculation.get_stable_phases() 

    # matrix_mole_fraction = point_calculation.get_value_of(ThermodynamicQuantity.mole_fraction_of_a_phase("FCC_A1#1"))
    # gamma_prime_mole_fraction = point_calculation.get_value_of(ThermodynamicQuantity.mole_fraction_of_a_phase("GAMMA_PRIME"))

    for key, _ in composition.items():
        matrixComposition[key] = round(point_calculation.get_value_of(ThermodynamicQuantity.composition_of_phase_as_weight_fraction("FCC_A1#1", key)))
        # gammaPrimeComposition[key] = round(point_calculation.get_value_of(ThermodynamicQuantity.composition_of_phase_as_weight_fraction("GAMMA_PRIME", key)))
    
    for element in composition:
        calculation.set_condition("w("+element+")", matrixComposition[element])

    property_diagram = calculation.calculate() 

    dgCH_T = property_diagram.get_values_of(ThermodynamicQuantity.user_defined_function('gm(bcc)-gm(fcc)'), ThermodynamicQuantity.temperature())

    d = {"GCHEM": dgCH_T[0], "T[K]": dgCH_T[1]}

    GData = pd.Dataframe(data=d, dtype=float)
    GData.sort_values(by=["T[K]"])
    GData.reset_index(drop=True)
    # make dataframe
    for index, __ in GData.iterrows():
        GData.loc[index, 'T[C]'] = float(GData.loc[index, 'T[K]'] - 273.15)
        GData.loc[index, 'w_FSol'] = w_fSol(matrixComposition, GData.loc[index, 'T[K]'], 0.4)
        GData.loc[index, 'asp'] = GData.loc[index, 'GCHEM']+GData.loc[index, 'w_FSol']
        GData.loc[index, 'YS'] = -1.425*(GData.loc[index, 'T[C]']-25) + YS_RT
        GData.loc[index, 'GsigUT'] = dGsigUT(GData.loc[index, 'YS'], 0.04)
        GData.loc[index, 'GsigCT'] = dGsigCT(GData.loc[index, 'YS'], 0.04)
        GData.loc[index, '-Gn-GsigUT'] = -1010-GData.loc[index, 'GsigUT']
        GData.loc[index, '-Gn-GsigCT'] = -1010-GData.loc[index, 'GsigCT']
        GData.loc[index, 'GNetUT'] = GData.loc[index, 'asp'] - GData.loc[index, '-Gn-GsigUT']
        GData.loc[index, 'GNetCT'] = GData.loc[index, 'asp'] - GData.loc[index, '-Gn-GsigCT']
    
    # ASK CLAY WHAT THIS IS
    index = abs(GData['T[C]'] -25.0).idxmin() #row number that gives
    asp_RT = GData.loc[index, 'asp'] # row, col position in dataframe

    return asp_RT

# Calculate Phase frac + APBE
def phase_frac_and_apbe(composition, calculation):

    # RT = 8.314*aging_temp 
    # matrixComposition = copy.deepcopy(composition)
    # gammaPrimeComposition = copy.deepcopy(composition)

    # elements = list(composition.keys())

    def APBE(sitefrac):
        gamma_0 = (0.1135*sitefrac)+0.1834
        return gamma_0
    
    for element in composition:
        calculation = calculation.set_condition("w("+element+")", composition[element])
    point_calculation = calculation.calculate()

    # matrix_mole_fraction = point_calculation.get_value_of(ThermodynamicQuantity.
    #                                                         mole_fraction_of_a_phase("FCC_A1#1"))

    gamma_prime_mole_fraction = point_calculation.get_value_of(ThermodynamicQuantity.
                                                                mole_fraction_of_a_phase("gamma-prime"))

    gammaprime_ti = point_calculation.get_value_of("Y(GAMMA_PRIME,TI#2)")
    apbe=APBE(gammaprime_ti)

    return [gamma_prime_mole_fraction, apbe]

# Calculate Strength and Driving Force
def strength_and_df(composition, temperature, dG_calculation, strength_calculation):

    RT = 8.314*temperature

    b = 2.5*10**-10   ## Burgers Vector
    G = 76.54*10**9   ## Shear Modulus in Pa
    T_L = (G*b**2)/2  ## Line Tension
    M = 3.06          ## Taylor Factor for FCC polycrystal
    vMAustenite = 7.225379*10**-6    ## Molar volume for FCC
    vMGammaPrime = 7.4002985*10**-6  ## Molar volume for Gamma-prime

    def APBE(sitefrac):
        gamma_0 = (0.1135*sitefrac)+0.1834
        return gamma_0
    def singleDislCutting(volfrac,avgRad,APBE):
        delTauHAM = (APBE/(2*b))*((((8*APBE*avgRad*volfrac)/(np.pi*G*b**2))**0.5)-volfrac)
        delSigmaHAM = M*delTauHAM
        return delSigmaHAM
    def calcVolfrac(molefrac):
        volTot = (molefrac*vMGammaPrime)+((1-molefrac)*vMAustenite)
        volfrac = (molefrac*vMGammaPrime)/volTot
        return volfrac    
    
    for element in composition:
        dG_calculation = dG_calculation.set_condition("w("+element+")", composition[element])
    point_calculation = dG_calculation.calculate()
    dG_gammaprime_norm = point_calculation.get_value_of(ThermodynamicQuantity.
                                                        normalized_driving_force_of_a_phase('GAMMA_PRIME'))
    dG_gammaprime = dG_gammaprime_norm * RT
    dG_eta_norm = point_calculation.get_value_of(
        ThermodynamicQuantity.normalized_driving_force_of_a_phase('ETA'))
    dG_eta = dG_eta_norm * RT
    dg_diff = dG_gammaprime - dG_eta

    for element in composition:
        strength_calculation = strength_calculation.set_condition("w("+element+")", composition[element])
    point_calculation = strength_calculation.calculate()

    gamma_prime_mole_fraction = point_calculation.get_value_of(ThermodynamicQuantity.
                                                                mole_fraction_of_a_phase("GAMMA_PRIME"))

    strengthaddition=singleDislCutting(calcVolfrac(gamma_prime_mole_fraction),0.000000015,APBE(point_calculation.get_value_of("Y(GAMMA_PRIME,TI#2)")))
    strength = float(341+strengthaddition/1000000)

    return [dg_diff, strength]

# Calculate average diameter of grain boundaries based on oxides
def zener(d,f,Z):
    """
    d : average pinning precipitate diameter
    f : total volume fraction of the pinning precipitates
    Z : fitting parameter related to the size distribution of the grains
    """
    return 4/3 * (3/2 - 2/Z) * d/f

# TCPython Initializer and Caller for Single Point Calc
def single_TC_caller(calcs, composition, temperature, disp=False):
    """
    Parameters
    ----------
    calcs : List
        Ex: ["printability", "stable_del_ferrite"]
        
    composition : Dict
        Ex: c_1217 = {
                       "C": 0.03/100,
                       "Cr": 12/100,
                       "Ni": 17/100,
                       "Mo": 1.3/100,
                       "Ti": 3.0/100,
                       "V": 0.3/100,
                       "Al": 0.2/100
                       #"B": 0.01/100
                       #"O": 0.015/100
                       }  # in wt-fraction #
        
    temperature : Dict
        Ex: temps = {
            "solution_temp" : 1000 + 273.15,
            "aging_temp" = 973.15
        }

    Returns
    -------
    result : Dict 
    [May contain, depending on what you wanted]
        "printability":
        fr : freezing Range [Float]
        csc : Cracking Susceptibility Coefficient [Float] 
        BCC_frac : Mole fraction of delta-Ferrite present at-print [Float]
        laves_frac : Mole fraction of laves present at-print [Float]

        "stable_del_ferrite":
        del_ferrite : Amount of stable del_ferrite post solution treatment

        "asp":
        asp : Austenite Stability Parameter

        "phase_frac_and_apbe":
        gamma_prime_mole_fraction : How much gamma Prime is in our post-aged matrix
        apbe : (??)

        "strength_and_df":
        dg_diff : dG difference between ____ and ____ (??)
        strength : Strength of alloy
   
    """
    elements = list(composition.keys())
    results = {}
    if len(calcs) == 0:
        return results
    with TCPython() as tcpython: 
        if "printability" in calcs: 
            # system definer info
            database = "TCFE10"
            dependent_element = "Fe"
            system = (tcpython.
                    set_cache_folder(os.path.basename(__file__) + "_cache").
                    select_database_and_elements(database, [dependent_element] + elements).
                    get_system_for_scheil_calculations())
        
            scheil_calculation = (system.with_scheil_calculation().
                                set_composition_unit(CompositionUnit.MASS_FRACTION))
            [fr, csc, BCC_frac, laves_frac] = printability(composition, scheil_calculation, system, disp=disp)
            results.update( {'fr' : fr, "csc" : csc, "BCC_frac" : BCC_frac, "laves_frac" : laves_frac} )
        
        if "stable_del_ferrite" in calcs: 
            # system definer info
            database = "TCFE10"
            dependent_element = "Fe"
            system = (tcpython
                    .set_cache_folder(os.path.basename(__file__) + "_cache")
                    .select_database_and_elements(database, [dependent_element] + elements)
                    .get_system()
                    .with_single_equilibrium_calculation())
                
            del_ferrite = stable_del_ferrite(composition, temperature["solution_temp"], system)
            results.update( {'del_ferrite' : del_ferrite} )

        if "asp" in calcs:
            singlepoint = (tcpython.
                set_cache_folder(os.path.basename(__file__) + "_cache").
                select_user_database_and_elements("nidata7.tdb", ["Fe"] + elements).
                without_default_phases().
                select_phase("FCC_A1").
                select_phase("gamma_prime").
                get_system().
                with_single_equilibrium_calculation().
                set_condition(ThermodynamicQuantity.temperature(), 973.15).
                set_gibbs_energy_addition_for('gamma_prime', -1456)
                )
            calculation = (tcpython.
                set_cache_folder(os.path.basename(__file__) + "_cache").
                select_user_database_and_elements("MART5.TDB", ["Fe"] + elements).
                without_default_phases().select_phase("FCC_A1").select_phase("BCC_A2").
                get_system().
                with_property_diagram_calculation().
                with_axis(CalculationAxis(ThermodynamicQuantity.temperature()).
                            set_min(temperature["start_temp"]).
                            set_max(temperature["end_temp"]).
                            with_axis_type(Linear().set_max_step_size(1))).
                disable_global_minimization().
                enable_step_separate_phases()
                )
            asp = asp(composition, singlepoint, calculation)
            results.update( {"asp" : asp} )

        if "phase_frac_and_apbe" in calcs:
            database = "nidata7.tdb"
            dependent_element = "fe"
            calculation = (tcpython. 
                        set_cache_folder(os.path.basename(__file__) + "_cache").
                        select_user_database_and_elements(database, [dependent_element] + elements).
                        without_default_phases().
                        select_phase("FCC_A1").
                        select_phase("gamma_prime").
                        get_system().
                        with_single_equilibrium_calculation().
                        set_condition(ThermodynamicQuantity.temperature(), 973.15).
                        set_gibbs_energy_addition_for('gamma_prime', -1456)
                        )
            [gamma_prime_mole_fraction, apbe] = phase_frac_and_apbe(composition, calculation)
            results.update( {"gamma_prime_mole_fraction" : gamma_prime_mole_fraction, "apbe" : apbe} )
                        
        if "strength_and_df" in calcs:
            # create and configure a single equilibrium calculation
            tcpython.set_ges_version(5)
            dG_calculation = (tcpython
                        .set_cache_folder(os.path.basename(__file__) + "_cache")
                        .select_user_database_and_elements(database, [dependent_element] + elements)
                        .without_default_phases()
                        .select_phase('FCC_A1').select_phase('GAMMA_PRIME').select_phase('ETA')
                        .get_system()
                        .with_single_equilibrium_calculation()
                        .set_condition(ThermodynamicQuantity.temperature(), temperature)
                        .set_phase_to_dormant('GAMMA_PRIME').set_phase_to_dormant('ETA')
                        .set_gibbs_energy_addition_for('GAMMA_PRIME', -1456)
                        )
            strength_calculation = (tcpython.
                        set_cache_folder(os.path.basename(__file__) + "_cache").
                        select_user_database_and_elements(database, [dependent_element] + elements).
                        without_default_phases().
                        select_phase("FCC_A1").
                        select_phase("GAMMA_PRIME").
                        get_system().
                        with_single_equilibrium_calculation().
                        set_condition(ThermodynamicQuantity.temperature(), temperature["aging_temp"]).
                        set_gibbs_energy_addition_for('GAMMA_PRIME', -1456)
                        ) 
            [dg_diff, strength] = strength_and_df(composition, temperature["aging_temp"], dG_calculation, strength_calculation)                   
            results.update( {"dg_diff" : dg_diff, "strength" : strength} )
    
    if len(results) == 0:
        print("calcs must be a list containing at least one of these: \"printability\", \"stable_del_fettire\", \"asp\", \"phase_frac_and_apbe\", \"strength_and_df\"")
    
    return results 

# TCPython Initializer and Caller for Matrix
def matrix_TC_caller(calcs, compositions_matrix, temperature):
    """
    Parameters
    ----------
    calcs : List
        Ex: ["printability", "stable_del_ferrite"]
        
    compositions_matrix : 2D array of Dicts
        Ex: c_1217 = {
                       "C": 0.03/100,
                       "Cr": 12/100,
                       "Ni": 17/100,
                       "Mo": 1.3/100,
                       "Ti": 3.0/100,
                       "V": 0.3/100,
                       "Al": 0.2/100
                       #"B": 0.01/100
                       #"O": 0.015/100
                       }  # in wt-fraction #
        
    temperature : Dict
        Ex: temps = {
            "solution_temp" : 1000 + 273.15,
            "aging_temp" : 973.15
        }

    Returns
    -------
    result : Dict 
    [May contain, depending on what you wanted]
        "printability":
        fr : freezing Range [Float]
        csc : Cracking Susceptibility Coefficient [Float] 
        BCC_frac : Mole fraction of delta-Ferrite present at-print [Float]
        laves_frac : Mole fraction of laves present at-print [Float]

        "stable_del_ferrite":
        del_ferrite : Amount of stable del_ferrite post solution treatment

        "asp":
        asp : Austenite Stability Parameter

        "phase_frac_and_apbe":
        gamma_prime_mole_fraction : How much gamma Prime is in our post-aged matrix
        apbe : (??)

        "strength_and_df":
        dg_diff : dG difference between ____ and ____ (??)
        strength : Strength of alloy
   
    """
    elements = list(compositions_matrix[0][0].keys())
    results = {}
    if len(calcs) == 0:
        return results
    with TCPython() as tcpython: 
        if "printability" in calcs: 
            # system definer info
            database = "TCFE10"
            dependent_element = "Fe"
            system = (tcpython.
                    set_cache_folder(os.path.basename(__file__) + "_cache").
                    select_database_and_elements(database, [dependent_element] + elements).
                    get_system_for_scheil_calculations())
        
            scheil_calculation = (system.with_scheil_calculation().
                                set_composition_unit(CompositionUnit.MASS_FRACTION))

            calc_matr = [[printability(compositions_matrix[x][y], scheil_calculation, system) for y in range(len(compositions_matrix[0]))] for x in range(len(compositions_matrix))]
        
            fr_matr = [[calc_matr[i][j][0] for i in range(len(compositions_matrix))] for j in range(len(compositions_matrix[0]))]
            csc_matr = [[calc_matr[i][j][1] for i in range(len(compositions_matrix))] for j in range(len(compositions_matrix[0]))]
            bcc_matr = [[calc_matr[i][j][2] for i in range(len(compositions_matrix))] for j in range(len(compositions_matrix[0]))]
            laves_matr = [[calc_matr[i][j][3] for i in range(len(compositions_matrix))] for j in range(len(compositions_matrix[0]))]

            results.update( {'fr' : fr_matr, "csc" : csc_matr, "BCC_frac" : bcc_matr, "laves_frac" : laves_matr} )
        
        if "stable_del_ferrite" in calcs: 
            # system definer info
            database = "TCFE10"
            dependent_element = "Fe"
            system = (tcpython
                    .set_cache_folder(os.path.basename(__file__) + "_cache")
                    .select_database_and_elements(database, [dependent_element] + elements)
                    .get_system()
                    .with_single_equilibrium_calculation())
                
            del_ferrite_matr = [[stable_del_ferrite(compositions_matrix[x][y], temperature["solution_temp"], system) for y in range(len(compositions_matrix[0]))] for x in range(len(compositions_matrix))]
            results.update( {'del_ferrite' : del_ferrite_matr} )

        if "asp" in calcs:
            singlepoint = (tcpython.
                set_cache_folder(os.path.basename(__file__) + "_cache").
                select_user_database_and_elements("nidata7.tdb", ["Fe"] + elements).
                without_default_phases().
                select_phase("FCC_A1").
                select_phase("gamma_prime").
                get_system().
                with_single_equilibrium_calculation().
                set_condition(ThermodynamicQuantity.temperature(), 973.15).
                set_gibbs_energy_addition_for('gamma_prime', -1456)
                )
            calculation = (tcpython.
                set_cache_folder(os.path.basename(__file__) + "_cache").
                select_user_database_and_elements("MART5.TDB", ["Fe"] + elements).
                without_default_phases().select_phase("FCC_A1").select_phase("BCC_A2").
                get_system().
                with_property_diagram_calculation().
                with_axis(CalculationAxis(ThermodynamicQuantity.temperature()).
                            set_min(temperature["start_temp"]).
                            set_max(temperature["end_temp"]).
                            with_axis_type(Linear().set_max_step_size(1))).
                disable_global_minimization().
                enable_step_separate_phases()
                )
            asp_matr = [[asp(compositions_matrix[x][y], singlepoint, calculation) for y in range(len(compositions_matrix[0]))] for x in range(len(compositions_matrix))]
            results.update( {"asp" : asp_matr} )

        if "phase_frac_and_apbe" in calcs:
            database = "nidata7.tdb"
            dependent_element = "fe"
            calculation = (tcpython. 
                        set_cache_folder(os.path.basename(__file__) + "_cache").
                        select_user_database_and_elements(database, [dependent_element] + elements).
                        without_default_phases().
                        select_phase("FCC_A1").
                        select_phase("gamma_prime").
                        get_system().
                        with_single_equilibrium_calculation().
                        set_condition(ThermodynamicQuantity.temperature(), 973.15).
                        set_gibbs_energy_addition_for('gamma_prime', -1456)
                        )

            res = [[phase_frac_and_apbe(compositions_matrix, calculation) for y in range(len(compositions_matrix[0]))] for x in range(len(compositions_matrix))]

            gamma_prime_mole_fraction =     [[res[i][j][0] for i in range(len(compositions_matrix))] for j in range(len(compositions_matrix[0]))]
            apbe    =                       [[res[i][j][1] for i in range(len(compositions_matrix))] for j in range(len(compositions_matrix[0]))]

            results.update( {"gamma_prime_mole_fraction" : gamma_prime_mole_fraction, "apbe" : apbe} )
                        
        if "strength_and_df" in calcs:
            # create and configure a single equilibrium calculation
            tcpython.set_ges_version(5)
            dG_calculation = (tcpython
                        .set_cache_folder(os.path.basename(__file__) + "_cache")
                        .select_user_database_and_elements(database, [dependent_element] + elements)
                        .without_default_phases()
                        .select_phase('FCC_A1').select_phase('GAMMA_PRIME').select_phase('ETA')
                        .get_system()
                        .with_single_equilibrium_calculation()
                        .set_condition(ThermodynamicQuantity.temperature(), temperature)
                        .set_phase_to_dormant('GAMMA_PRIME').set_phase_to_dormant('ETA')
                        .set_gibbs_energy_addition_for('GAMMA_PRIME', -1456)
                        )
            strength_calculation = (tcpython.
                        set_cache_folder(os.path.basename(__file__) + "_cache").
                        select_user_database_and_elements(database, [dependent_element] + elements).
                        without_default_phases().
                        select_phase("FCC_A1").
                        select_phase("GAMMA_PRIME").
                        get_system().
                        with_single_equilibrium_calculation().
                        set_condition(ThermodynamicQuantity.temperature(), temperature["aging_temp"]).
                        set_gibbs_energy_addition_for('GAMMA_PRIME', -1456)
                        ) 

            res = [[strength_and_df(compositions_matrix, temperature["aging_temp"], dG_calculation, strength_calculation) for y in range(len(compositions_matrix[0]))] for x in range(len(compositions_matrix))]                   

            dg_diff =  [[res[i][j][0] for i in range(len(compositions_matrix))] for j in range(len(compositions_matrix[0]))]
            strength = [[res[i][j][1] for i in range(len(compositions_matrix))] for j in range(len(compositions_matrix[0]))]

            results.update( {"dg_diff" : dg_diff, "strength" : strength} )
    
    if len(results) == 0:
        print("calcs must be a list containing at least one of these: \"printability\", \"stable_del_fettire\", \"asp\", \"phase_frac_and_apbe\", \"strength_and_df\"")
    
    return results 