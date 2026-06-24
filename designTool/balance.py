# GENERAL IMPORTS
import numpy as np
from .constants import gravity
from .geometry import change_sweep

#========================================

def balance(airplane):

    # Unpack dictionary
    W0 = airplane['thrust_matching']['W0']
    W_payload = airplane['inputs']['W_payload']
    xcg_payload = airplane['inputs']['xcg_payload']
    W_crew = airplane['inputs']['W_crew']
    xcg_crew = airplane['inputs']['xcg_crew']
    W_empty = airplane['thrust_matching']['W_empty']
    xcg_empty = airplane['empty_weight']['xcg_empty']
    W_fuel = airplane['thrust_matching']['W_fuel']
    
    Mach_cruise = airplane['inputs']['Mach_cruise']
    
    S_w = airplane['inputs']['S_w']
    AR_eff = airplane['aerodynamics']['AR_eff']
    taper_w = airplane['inputs']['taper_w']
    sweep_w = airplane['inputs']['sweep_w']
    b_w = airplane['geometry']['b_w']
    xr_w = airplane['inputs']['xr_w']
    zr_w = airplane['inputs']['zr_w']
    cr_w = airplane['geometry']['cr_w']
    ct_w = airplane['geometry']['ct_w']
    xm_w = airplane['geometry']['xm_w']
    cm_w = airplane['geometry']['cm_w']
    tcr_w = airplane['inputs']['tcr_w']
    tct_w = airplane['inputs']['tct_w']
    
    S_h = airplane['geometry']['S_h']
    AR_h = airplane['inputs']['AR_h']
    sweep_h = airplane['inputs']['sweep_h']
    b_h = airplane['geometry']['b_h']
    cr_h = airplane['geometry']['cr_h']
    ct_h = airplane['geometry']['ct_h']
    xm_h = airplane['geometry']['xm_h']
    zm_h = airplane['geometry']['zm_h']
    cm_h = airplane['geometry']['cm_h']
    eta_h = airplane['inputs']['eta_h']
    Lc_h = airplane['inputs']['Lc_h']
    
    Cvt = airplane['inputs']['Cvt']
    
    L_f = airplane['inputs']['L_f']
    D_f = airplane['inputs']['D_f']
    
    y_n = airplane['inputs']['y_n']
    
    T0 = airplane['thrust_matching']['T0']
    n_engines = airplane['inputs']['n_engines']

    c_tank_c_w = airplane['inputs']['c_tank_c_w']
    x_tank_c_w = airplane['inputs']['x_tank_c_w']
    b_tank_b_w_start = airplane['inputs']['b_tank_b_w_start']
    b_tank_b_w_end = airplane['inputs']['b_tank_b_w_end']

    rho_fuel = airplane['inputs']['rho_fuel']
    
    CLmaxTO = airplane['thrust_matching']['CLmaxTO']

    # @REMOVE

    # Tank CG
    V_maxfuel, W_maxfuel, xcg_fuel, ycg_fuel = tank_properties(cr_w, ct_w, tcr_w, tct_w, b_w, sweep_w, xr_w,
                                                               x_tank_c_w, c_tank_c_w, b_tank_b_w_start, b_tank_b_w_end,
                                                               rho_fuel, gravity)

    # Check available fuel volume
    tank_excess = W_maxfuel/W_fuel - 1

    ### CG RANGE
    # Here we evaluate 5 load scenarios

    # Empty airplane
    W_e = W_empty
    xcg_e = xcg_empty

    # Crew
    W_ec = W_empty + W_crew
    xcg_ec = (W_empty*xcg_empty + W_crew*xcg_crew)/(W_empty + W_crew)

    # Payload and crew
    W_epc = W_empty + W_payload + W_crew
    xcg_epc = (W_empty*xcg_empty + W_payload*xcg_payload + W_crew*xcg_crew)/(W_empty + W_payload + W_crew)

    # Fuel and crew
    W_efc = W_empty + W_fuel + W_crew
    xcg_efc = (W_empty*xcg_empty + W_fuel*xcg_fuel + W_crew*xcg_crew)/(W_empty + W_fuel + W_crew)

    # Payload, crew, and fuel (full airplane)
    W_efpc = W0
    xcg_efpc = (W_empty*xcg_empty + W_fuel*xcg_fuel + W_payload*xcg_payload + W_crew*xcg_crew)/W0

    # Find CG range (adding a 2% margin on either side)
    xcg_list = [xcg_e, xcg_ec, xcg_epc, xcg_efc, xcg_efpc]
    xcg_fwd = min(xcg_list) - 0.02*cm_w
    xcg_aft = max(xcg_list) + 0.02*cm_w

    ### NEUTRAL POINT

    # Compressibility correction for aerodynamic center
    # Raymer, Eq. 16.12
    delta_xac = 0.0 # Removing this for now to be more conservative (the aircraft must be stable at low speeds)
    '''
    if Mach_cruise < 0.4:
        delta_xac = 0.0
    elif Mach_cruise < 1.1:
        delta_xac = 0.26*(Mach_cruise-0.4)**2.5
    else:
        delta_xac = 0.112 - 0.004*Mach_cruise
    '''

    # Wing lift slope (Raymer Eq 12.6) - ATTENTION: The Equation is wrong in Raymer
    # the main reference has a term beta**4 instead of beta**2
    # Here we assume that clalpha/2/pi = 0.95 (without beta correction, that was wrong in Raymer)
    sweep_maxt_w = change_sweep(0.25, 0.40,
                                    sweep_w, b_w/2, cr_w, ct_w) # Sweep at max. thickness
    beta2 = 1-Mach_cruise**2
    CLa_w = 2*np.pi*AR_eff/(2 + np.sqrt(4 + AR_eff**2*beta2/0.95**2*(1+np.tan(sweep_maxt_w)**2/beta2)))

    # Wing aerodynamic center at 25% mac
    xac_w = xm_w + 0.25*cm_w + delta_xac*np.sqrt(S_w)

    # Fuselage moment slope (Raymer Eq 16.25)
    K_fus = 0.1462*np.exp(4.8753*(xr_w + 0.25*cr_w)/L_f)
    CMa_f = K_fus*D_f**2*L_f/cm_w/S_w

    # Loss of lift due to fuselage
    # (Roskam: Methods for Estimating Stability and Control Derivatives of Conventional Subsonic Airplanes)
    # Eq. 3.7
    #K_wb = 1.0 - 0.25*(D_f/b_w)**2 + 0.025*(D_f/b_w)
    #CLa_wb = K_wb*CLa_w
    # Raymer
    CLa_wf = CLa_w*0.98

    # Neutral point of the wing-fuselage combination
    xac_wf = xac_w - CMa_f/CLa_wf*cm_w

    # HT lift slope (Raymer Eq 12.6)
    sweep_maxt_h = change_sweep(0.25, 0.40,
                                    sweep_h, b_h/2, cr_h, ct_h) # Sweep at max. thickness
    CLa_h = 2*np.pi*AR_h/(2 + np.sqrt(4 + AR_h**2*beta2/0.95**2*(1+np.tan(sweep_maxt_h)**2/beta2)))*0.98

    # HT aerodynamic center at 25% mac
    xac_h = xm_h + 0.25*cm_h + delta_xac*np.sqrt(S_h)

    # Downwash (Nelson Eq 2.23) - I am not using this anymore because it seems like it overestimates downwash
    #deda = 2*CLa_w/np.pi/AR_eff

    # Downwash (Roskam: Methods for Estimating Stability and Control Derivatives of Conventional Subsonic Airplanes)
    # Eq. 3.11
    K_A = 1.0/AR_eff - 1.0/(1.0 + AR_eff**1.7)
    K_lambda = (10.0 - 3.0*taper_w)/7.0
    h_H = abs(zm_h - zr_w)
    L_H = Lc_h*cm_w
    K_H = (1-h_H/b_w)/(2*L_H/b_w)**(1/3)
    CLa_w0 = 2*np.pi*AR_eff/(2 + np.sqrt(4 + AR_eff**2/0.95**2*(1+np.tan(sweep_maxt_w)**2)))
    deda = 4.44*(K_A*K_lambda*K_H*np.sqrt(np.cos(sweep_w)))**1.19*CLa_w/CLa_w0

    # Neutral point position (Raymer Eq 16.9 and Eq 16.23)
    xnp = (CLa_wf*xac_wf + eta_h*S_h/S_w*CLa_h*(1-deda)*xac_h)/(CLa_wf + eta_h*S_h/S_w*CLa_h*(1-deda))

    # Static margin
    SM_fwd = (xnp - xcg_fwd)/cm_w
    SM_aft = (xnp - xcg_aft)/cm_w

    # VERTICAL TAIL VERIFICATION FOR OEI CONDITION

    # Compute stall factor assuming that V2 = 1.1*Vmc and V2=1.2*Vs (FAR 25.107)
    ks = 1.2/1.1

    # Compute required lift for the vertical tail
    CLv = y_n/b_w*CLmaxTO/ks**2*T0/W0/n_engines/Cvt

    # @REMOVE

    # Update dictionary
    airplane['balance'] = {}
    airplane['balance']['xcg_fwd'] = xcg_fwd
    airplane['balance']['xcg_aft'] = xcg_aft
    airplane['balance']['xnp'] = xnp
    airplane['balance']['SM_fwd'] = SM_fwd
    airplane['balance']['SM_aft'] = SM_aft
    airplane['balance']['tank_excess'] = tank_excess
    airplane['balance']['V_maxfuel'] = V_maxfuel
    airplane['balance']['W_maxfuel'] = W_maxfuel
    airplane['balance']['xcg_fuel'] = xcg_fuel
    airplane['balance']['CLv'] = CLv

    airplane['balance']['CG_hist'] = {'xcg_e':xcg_e,
                                      'W_e':W_e,
                                      'xcg_ec':xcg_ec,
                                      'W_ec':W_ec,
                                      'xcg_epc':xcg_epc,
                                      'W_epc':W_epc,
                                      'xcg_efc':xcg_efc,
                                      'W_efc':W_efc,
                                      'xcg_efpc':xcg_efpc,
                                      'W_efpc':W_efpc,}

    return None

#----------------------------------------

def tank_properties(cr_w, ct_w, tcr_w, tct_w, b_w, sweep_w, xr_w,
                    x_tank_c_w, c_tank_c_w, b_tank_b_w_start, b_tank_b_w_end,
                    rho_fuel, gravity):
    '''
    This function computes the maximum fuel tank volume and center of gravity.
    We assume that the tank has a prism shape.

    c_tank_c_w: float -> fraction of the chord where tank begins (0-leading edge, 1-trailing edge)
    c_tank_c_w: float -> fraction of the chord occupied by the tank (between 0 and 1)
    bf_w_start: float -> semi-span fraction where tank begins (0-root, 1-tip)
    bf_w_end: float -> semi-span fraction where tank ends (0-root, 1-tip)
    '''

    # Compute the local chords where the tank begins and ends
    c_tank_start = cr_w + b_tank_b_w_start*(ct_w - cr_w)
    c_tank_end = cr_w + b_tank_b_w_end*(ct_w - cr_w)

    # Compute the local thickness where the tank begins and ends
    tc_tank_start = tcr_w + b_tank_b_w_start*(tct_w - tcr_w)
    tc_tank_end = tcr_w + b_tank_b_w_end*(tct_w - tcr_w)

    # Compute the prism area where the tank begins.
    # We assume that this face is rectangular, and that its height
    # is 85% of the maximum airfoil thickness (Gudmundsson, page 87).
    ll = c_tank_start*c_tank_c_w
    hh = c_tank_start*tc_tank_start*0.85
    S1 = ll*hh

    # Compute the prism area where the tank ends.
    ll = c_tank_end*c_tank_c_w
    hh = c_tank_end*tc_tank_end*0.85
    S2 = ll*hh

    # Compute distance between prism faces along the wing span
    Lprism = 0.5*b_w*(b_tank_b_w_end-b_tank_b_w_start)

    # Compute fuel volume with the prism expression (Torenbeek Fig B-4, pg 448).
    # We multiply by 2 to take into account both semi-wings.
    # The 0.91 factor is to take into account internal structures and fuel expansion,
    # as suggested by Torenbeek.
    V_maxfuel = 0.91*2*Lprism/3*(S1 + S2 + np.sqrt(S1*S2))

    # Compute corresponding fuel weight
    W_maxfuel = V_maxfuel*rho_fuel*gravity

    # Compute the span-wise distance between the first prism face and its center of gravity
    # using the expression from Jenkinson, Fig 7.13, pg 148.
    Lprism_cg = Lprism/4*(S1 + 3*S2 + 2*np.sqrt(S1*S2))/(S1 + S2 + np.sqrt(S1*S2))

    # Now find the span-wise distance between the tank CG and the aircraft centerline
    ycg_fuel = Lprism_cg + 0.5*b_w*b_tank_b_w_start

    # Find the sweep angle at the chord position located on the middle of the chord
    # fraction occupied by the fuel tank
    c_pos = x_tank_c_w + 0.5*c_tank_c_w

    # Sweep at the tank center line
    sweep_tank = change_sweep(0.25, c_pos, sweep_w, b_w/2, cr_w, ct_w)

    # Longitudinal position of the tank CG
    xcg_fuel = xr_w + cr_w*c_pos + ycg_fuel*np.tan(sweep_tank)


    return V_maxfuel, W_maxfuel, xcg_fuel, ycg_fuel