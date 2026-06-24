# GENERAL IMPORTS
import numpy as np
from scipy.optimize import root_scalar
from .constants import ft2m, kt2ms, lb2N, nm2m, gravity, gamma_air, R_air
from .geometry import control_surface_area_fraction
from .propulsion import engineTSFC
from .aerodynamics import aerodynamics
from .auxiliary import atmosphere

#========================================

def weight(W0_guess, T0_guess, airplane):

    # Unpacking dictionary
    W_payload = airplane['inputs']['W_payload']
    W_crew = airplane['inputs']['W_crew']
    range_cruise = airplane['inputs']['range_cruise']

    # # Set iterator
    # delta = 1000

    # while abs(delta) > 10:

    #     # We need to call fuel_weight first since it
    #     # calls the aerodynamics module to get Swet_f used by
    #     # the empty weight function
    #     W_fuel, W_cruise = fuel_weight(W0_guess, airplane, range_cruise=range_cruise, update_Mf_hist=True)

    #     W_empty = empty_weight(W0_guess, T0_guess, airplane)

    #     W0 = W_empty + W_fuel + W_payload + W_crew

    #     delta = W0 - W0_guess

    #     W0_guess = W0

    # return W0, W_empty, W_fuel, W_cruise
    
    # Define the residual function we want to drive to zero
    def weight_residual(W0_test):
        W_fuel, _ = fuel_weight(W0_test, airplane, range_cruise=range_cruise, update_Mf_hist=False)
        W_empty = empty_weight(W0_test, T0_guess, airplane)
        return W0_test - (W_empty + W_fuel + W_payload + W_crew)

    # Use the Secant method to rapidly find the root (the exact W0)
    # x0 is the initial guess, x1 is a slightly perturbed second guess to start it
    sol = root_scalar(weight_residual, x0=W0_guess, x1=W0_guess*0.95, method='secant', xtol=1e-7)
    
    W0_final = sol.root
    
    # Run one last time to get the final split of components
    W_fuel, W_cruise = fuel_weight(W0_final, airplane, range_cruise=range_cruise, update_Mf_hist=True)
    W_empty = empty_weight(W0_final, T0_guess, airplane)
    
    return W0_final, W_empty, W_fuel, W_cruise

#----------------------------------------

def empty_weight(W0_guess, T0_guess, airplane):

    # Unpack dictionary
    S_w = airplane['inputs']['S_w']
    AR_eff = airplane['aerodynamics']['AR_eff']
    b_w = airplane['geometry']['b_w']
    taper_w = airplane['inputs']['taper_w']
    sweep_w = airplane['inputs']['sweep_w']
    xm_w = airplane['geometry']['xm_w']
    cm_w = airplane['geometry']['cm_w']
    tcr_w = airplane['inputs']['tcr_w']

    flap_type = airplane['inputs']['flap_type']
    c_flap_c_wing = airplane['inputs']['c_flap_c_wing']
    b_flap_b_wing = airplane['inputs']['b_flap_b_wing']
    slat_type = airplane['inputs']['slat_type']
    c_slat_c_wing = airplane['inputs']['c_slat_c_wing']
    b_slat_b_wing = airplane['inputs']['b_slat_b_wing']
    c_ail_c_wing = airplane['inputs']['c_ail_c_wing']
    b_ail_b_wing = airplane['inputs']['b_ail_b_wing']
    
    S_h = airplane['geometry']['S_h']
    xm_h = airplane['geometry']['xm_h']
    cm_h = airplane['geometry']['cm_h']
    
    S_v = airplane['geometry']['S_v']
    xm_v = airplane['geometry']['xm_v']
    cm_v = airplane['geometry']['cm_v']
    
    L_f = airplane['inputs']['L_f']
    D_f = airplane['inputs']['D_f']
    Swet_f = airplane['aerodynamics']['Swet_f']
    
    n_engines = airplane['inputs']['n_engines']
    x_n = airplane['inputs']['x_n']
    L_n = airplane['inputs']['L_n']
    
    x_nlg = airplane['inputs']['x_nlg']
    x_mlg = airplane['inputs']['x_mlg']
    
    altitude_cruise = airplane['inputs']['altitude_cruise']
    Mach_cruise = airplane['inputs']['Mach_cruise']
    
    airplane_type = airplane['inputs']['type']

    engine = airplane['inputs']['engine']

    # @REMOVE

    # Select appropriate parameters for weight regression
    
    if airplane_type == 'transport':
        
        # Wing weight (Raymer Eq 15.25)
        
        # Ultimate load factor
        Nz = 1.5*max(2.5, min(3.8, 2.1 + 24000/(10000 + W0_guess/lb2N)))

        # Area of control surfaces
        Scsw_frac = 0.0

        # Flap contribution
        if flap_type is not None:

            # Area multiplier to account for multiple surfaces
            # These factors are based on Torenbeek Eq. C-10
            flap_factor = {'plain' : 1,
                           'single slotted' : 1.15*1.25,
                           'double slotted' : 1.30*1.25,
                           'triple slotted' : 1.45*1.25,}

            Scsw_frac = Scsw_frac + control_surface_area_fraction(xf=c_flap_c_wing,
                                             y1=D_f/b_w,
                                             y2=b_flap_b_wing,
                                             taper=taper_w)*flap_factor[flap_type]

        # Slat contribution
        if slat_type is not None:

            # Area multiplier to account for multiple surfaces
            # These factors are based on Torenbeek Eq. C-10
            slat_factor = {'fixed' : 0,
                           'flap' : 1,
                           'kruger' : 1,
                           'slat' : 1.25,}

            Scsw_frac = Scsw_frac + control_surface_area_fraction(xf=c_slat_c_wing,
                                             y1=D_f/b_w,
                                             y2=b_slat_b_wing,
                                             taper=taper_w)*slat_factor[slat_type]

        # Aileron contribution
        Scsw_frac = Scsw_frac + control_surface_area_fraction(xf=c_ail_c_wing,
                                         y1=1-b_ail_b_wing,
                                         y2=1.0,
                                         taper=taper_w)

        Scsw = Scsw_frac*S_w # Area of control surfaces
        
        # I increased the AR_w exponent from 0.5 to 0.55 to make it more sensitive.
        # Otherwise, the optimum would be around AR = 12, which may be too optimistic.
        W_w = 0.0051*(W0_guess*Nz/lb2N)**0.557*(S_w/ft2m**2)**0.649*AR_eff**0.55*tcr_w**(-0.4)*(1+taper_w)**0.1/np.cos(sweep_w)*(Scsw/ft2m**2)**0.1*lb2N
        xcg_w = xm_w + 0.4*cm_w
        
        # Surface densities for remaining components (kg/m2) - Raymer Tab 15.2
        W_h_dens = 27
        W_v_dens = 27
        W_f_dens = 24
        W_lg_fact = 0.043
        W_eng_fact = 1.3
        W_allelse_fact = 0.17
        
    elif airplane_type == 'fighter':

        # Wing weight (Raymer Eq 15.1)
        Nz = 1.5*6.0 # Ultimate load factor
        # Area of control surfaces
        Scsw_frac = 0.0

        # Flap contribution
        if flap_type is not None:

            # Area multiplier to account for multiple surfaces
            flap_factor = {'plain' : 1,
                           'single slotted' : 1.2,
                           'double slotted' : 2,
                           'triple slotted' : 3,}

            Scsw_frac = Scsw_frac + control_surface_area_fraction(alpha=c_flap_c_wing,
                                             beta1=D_f/b_w,
                                             beta2=b_flap_b_wing,
                                             taper=taper_w)*flap_factor[flap_type]

        # Slat contribution
        if slat_type is not None:

            # Area multiplier to account for multiple surfaces
            slat_factor = {'fixed' : 0,
                           'flap' : 1,
                           'kruger' : 1,
                           'slat' : 1.2,}

            Scsw_frac = Scsw_frac + control_surface_area_fraction(alpha=c_slat_c_wing,
                                             beta1=D_f/b_w,
                                             beta2=b_slat_b_wing,
                                             taper=taper_w)*slat_factor[slat_type]

        # Aileron contribution
        Scsw_frac = Scsw_frac + control_surface_area_fraction(alpha=c_ail_c_wing,
                                         beta1=1-b_ail_b_wing,
                                         beta2=1.0,
                                         taper=taper_w)

        Scsw = Scsw_frac*S_w # Area of control surfaces
        W_w = 0.0103*(W0_guess*Nz/lb2N)**0.5*(S_w/ft2m**2)**0.622*AR_eff**0.785*tcr_w**(-0.4)*(1+taper_w)**0.05/np.cos(sweep_w)*(Scsw/ft2m**2)**0.04*lb2N
        xcg_w = xm_w + 0.4*cm_w
        
        # Surface densities for remaining components (kg/m2) - Raymer Tab 15.2
        W_h_dens = 20
        W_v_dens = 26
        W_f_dens = 23
        W_lg_fact = 0.033
        W_eng_fact = 1.3
        W_allelse_fact = 0.17
        
    elif airplane_type == 'general':

        # Cruise dynamic pressure
        atm_data = atmosphere(altitude_cruise)
        a_cruise = atm_data['speed_of_sound']
        rho = atm_data['density']

        v_cruise = Mach_cruise*a_cruise
        q_cruise = 0.5*rho*v_cruise**2

        # Wing weight (Raymer Eq 15.46)
        Nz = 1.5*4.4 # Ultimate load factor
        W_w = 0.036*(W0_guess*Nz/lb2N)**0.49*(S_w/ft2m**2)**0.758*(AR_eff/np.cos(sweep_w)**2)**0.6*(100*tcr_w/np.cos(sweep_w))**(-0.3)*(taper_w)**0.04*q_cruise**0.006*lb2N
        xcg_w = xm_w + 0.4*cm_w
        
        # Surface densities for remaining components (kg/m2) - Raymer Tab 15.2
        W_h_dens = 10
        W_v_dens = 10
        W_f_dens = 7
        W_lg_fact = 0.057
        W_eng_fact = 1.4
        W_allelse_fact = 0.1

    # Use Raymer Tab 15.2 for the remaining components

    W_h = S_h*gravity*W_h_dens
    xcg_h = xm_h + 0.4*cm_h

    W_v = S_v*gravity*W_v_dens
    xcg_v = xm_v + 0.4*cm_v

    W_f = Swet_f*gravity*W_f_dens
    #Kdoor = 1.0
    #Klg = 1.0 #1.0 for low wing, 1.12 for high wing
    #Kws = 0.75*(1+2*taper_w)/(1+taper_w)*b_w*np.tan(sweep_w)/L_f
    #W_f = 0.3280*Kdoor*Klg*(W0_guess*Nz/lb2N)**0.5*(L_f/ft2m)**0.25*(Swet_f/ft2m**2)**0.302*(1+Kws)**0.04*(L_f/D_f)**0.1*lb2N*1.6 #+ 0.0577*(5*91*9.81/lb2N)**0.1*(107*91*9.81/lb2N)**0.393*(Swet_f/ft2m**2)**0.75*lb2N
    xcg_f = 0.45*L_f

    # Check if LG is active
    if x_nlg is not None:
        
        W_nlg = 0.15*W0_guess*W_lg_fact
        xcg_nlg = x_nlg
    
        W_mlg = 0.85*W0_guess*W_lg_fact
        xcg_mlg = x_mlg
        
    else:
        
        W_nlg = 0.0
        xcg_nlg = 0.0
    
        W_mlg = 0.0
        xcg_mlg = 0.0

    # Engine weight
    T_eng = T0_guess/n_engines
    
    if 'weight' in engine:
        
        # The user already gave the engine weight
        W_eng = engine['weight']
    
    elif 'turbofan' in engine['model']:
        BPR = engine['BPR']
        
        # Turbofan weight (Raymer Eq. 10.4)
        W_eng = gravity*14.7*(T_eng/1000.0)**1.1*np.exp(-0.045*BPR)

    W_eng_installed = n_engines*W_eng*W_eng_fact
    xcg_eng = x_n + 0.5*L_n

    # All else weight
    W_allelse = W_allelse_fact*W0_guess
    xcg_allelse = 0.45*L_f

    # Empty weight
    W_empty = W_w + W_h + W_v + W_f + W_nlg + W_mlg + W_eng_installed + W_allelse

    # Empty weight CG
    xcg_empty = (W_w*xcg_w + W_h*xcg_h + W_v*xcg_v + W_f*xcg_f +
             W_nlg*xcg_nlg + W_mlg*xcg_mlg + W_eng_installed*xcg_eng +
             W_allelse*xcg_allelse)/W_empty

    # @REMOVE

    # Update dictionary
    airplane['empty_weight'] = {}
    airplane['empty_weight']['W_w'] = W_w
    airplane['empty_weight']['W_h'] = W_h
    airplane['empty_weight']['W_v'] = W_v
    airplane['empty_weight']['W_f'] = W_f
    airplane['empty_weight']['W_nlg'] = W_nlg
    airplane['empty_weight']['W_mlg'] = W_mlg
    airplane['empty_weight']['W_eng'] = W_eng_installed
    airplane['empty_weight']['W_allelse'] = W_allelse
    airplane['empty_weight']['W_empty'] = W_empty
    airplane['empty_weight']['xcg_w'] = xcg_w
    airplane['empty_weight']['xcg_h'] = xcg_h
    airplane['empty_weight']['xcg_v'] = xcg_v
    airplane['empty_weight']['xcg_f'] = xcg_f
    airplane['empty_weight']['xcg_nlg'] = xcg_nlg
    airplane['empty_weight']['xcg_mlg'] = xcg_mlg
    airplane['empty_weight']['xcg_eng'] = xcg_eng
    airplane['empty_weight']['xcg_allelse'] = xcg_allelse
    airplane['empty_weight']['xcg_empty'] = xcg_empty

    return W_empty

#----------------------------------------

def fuel_weight(W0_guess, airplane, range_cruise, update_Mf_hist=False):
    '''
    This function estimates the fuel consumed using the weight fractions approach.

    The update_Mf_hist flag should only be activated when the user wants
    to store intermediate information of fuel estimation.

    The mission profile is:

    take-off -> climb -> cruise -> descent ->  alternate cruise -> loiter ->landing
    '''

    # Unpacking dictionary
    S_w = airplane['inputs']['S_w']
    
    altitude_cruise = airplane['inputs']['altitude_cruise']
    Mach_cruise = airplane['inputs']['Mach_cruise']
    
    time_loiter = airplane['inputs']['time_loiter']
    altitude_loiter = airplane['inputs']['altitude_loiter']
    
    altitude_altcruise = airplane['inputs']['altitude_altcruise']
    Mach_altcruise = airplane['inputs']['Mach_altcruise']
    range_altcruise = airplane['inputs']['range_altcruise']
    
    airplane_type = airplane['inputs']['type']

    # @REMOVE

    # Get engine TSFC
    C_cruise,_ = engineTSFC(Mach_cruise, altitude_cruise, airplane)
    C_altcruise,_ = engineTSFC(Mach_altcruise, altitude_altcruise, airplane)
    
    # Get fractions according to the category

    if airplane_type == 'transport':
        Mf_start   = 0.990
        Mf_taxi    = 0.990
        Mf_takeoff = 0.995
        Mf_climb   = 0.980
        Mf_descent = 0.990
        Mf_landing = 0.992
        
    if airplane_type == 'fighter':
        Mf_start   = 0.990
        Mf_taxi    = 0.990
        Mf_takeoff = 0.990
        Mf_climb   = 0.930
        Mf_descent = 0.990
        Mf_landing = 0.995
        
    if airplane_type == 'general':
        Mf_start   = 0.995
        Mf_taxi    = 0.997
        Mf_takeoff = 0.998
        Mf_climb   = 0.992
        Mf_descent = 0.993
        Mf_landing = 0.993

    ### Cruise

    # Compute weight at the beginning of the cruise
    W_cruise = W0_guess*Mf_start*Mf_taxi*Mf_takeoff*Mf_climb

    # Atmospheric conditions at cruise altitude
    atm_data = atmosphere(altitude_cruise)
    a_cruise = atm_data['speed_of_sound']
    rho = atm_data['density']

    # Cruise speed
    v_cruise = Mach_cruise*a_cruise

    # Cruise CL
    CL = 2.0*W_cruise/rho/S_w/v_cruise**2

    # Cruise Drag
    CD, _, dragDict = aerodynamics(airplane, Mach_cruise, altitude_cruise,
                                   CL, n_engines_failed=0, highlift_config='clean',
                                   lg_down=0, h_ground=0)

    # Aerodynamic efficiency
    LD_cruise = CL/CD

    # Breguet Equation
    Mf_cruise = np.exp(-range_cruise*C_cruise/v_cruise/LD_cruise)

    ### Cruise 2

    # Compute weight at the beginning of the alternative cruise
    W_altcruise = W_cruise*Mf_cruise*Mf_descent

    # Atmospheric conditions at cruise altitude
    atm_data = atmosphere(altitude_altcruise)
    a_altcruise = atm_data['speed_of_sound']
    rho = atm_data['density']

    # Cruise speed
    v_altcruise = Mach_altcruise*a_altcruise

    # Cruise CL
    CL = 2.0*W_altcruise/rho/S_w/v_altcruise**2

    # Cruise CD
    CD, _, dragDict = aerodynamics(airplane, Mach_altcruise, altitude_altcruise,
                                       CL, n_engines_failed=0, highlift_config='clean',
                                       lg_down=0, h_ground=0)

    # Aerodynamic efficiency
    LD_altcruise = CL/CD 

    # Breguet Equation
    Mf_altcruise = np.exp(-range_altcruise*C_altcruise/v_altcruise/LD_altcruise)

    ### Loiter

    # Loiter at max L/D and low altitude
    # For now, we take the cruise CD0, CDwave and K to estimate L/Dmax
    Mach_loiter = 0.5
    CD, _, dragDict = aerodynamics(airplane, Mach_loiter, altitude_loiter,
                                   CL, n_engines_failed=0, highlift_config='clean',
                                   lg_down=0, h_ground=0)
    LD_max = 0.5/np.sqrt((dragDict['CD0']+dragDict['CDwave'])*dragDict['K'])
    LD_loiter = LD_max

    # Factor to fuel consumption (Correction based on Raymer Tab. 3.3)
    C_loiter = C_cruise - 0.1/3600.0

    # Breguet Equation
    Mf_loiter = np.exp(-time_loiter*C_loiter/LD_loiter)

    ### Overall mass fraction
    Mf = Mf_start*Mf_taxi*Mf_takeoff*Mf_climb*Mf_cruise*Mf_descent*Mf_altcruise*Mf_loiter*Mf_landing

    ### Fuel weight (Raymer Eq 3.13)
    trapped_fuel_factor = 1.01
    W_fuel = trapped_fuel_factor*(1-Mf)*W0_guess

    # @REMOVE

    # Store the history history of fuel consumed only when requested
    if update_Mf_hist:

        C_hist = {} # Dictionary to hold TSFC of all mission phases
        C_hist['cruise'] = C_cruise
        C_hist['loiter'] = C_loiter
        C_hist['altcruise'] = C_altcruise

        Mf_hist = {} # Dictionary to hold mass fractions of all mission phases
        Mf_hist['start'] = Mf_start
        Mf_hist['taxi'] = Mf_taxi
        Mf_hist['takeoff'] = Mf_takeoff
        Mf_hist['climb'] = Mf_climb
        Mf_hist['cruise'] = Mf_cruise
        Mf_hist['descent'] = Mf_descent
        Mf_hist['altcruise'] = Mf_altcruise
        Mf_hist['loiter'] = Mf_loiter
        Mf_hist['landing'] = Mf_landing

        LD_hist = {} # Dictionary to hold L/D of all mission phases
        LD_hist['cruise'] = LD_cruise
        LD_hist['altcruise'] = LD_altcruise
        LD_hist['loiter'] = LD_loiter

        airplane['fuel_weight'] = {}
        airplane['fuel_weight']['Mf_hist'] = Mf_hist
        airplane['fuel_weight']['LD_hist'] = LD_hist
        airplane['fuel_weight']['C_hist'] = C_hist
        airplane['fuel_weight']['trapped_fuel_factor'] = trapped_fuel_factor

    return W_fuel, W_cruise