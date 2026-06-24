# GENERAL IMPORTS
import numpy as np
from .geometry import change_sweep
from .auxiliary import atmosphere

#========================================

def aerodynamics(airplane, Mach, altitude, CL,
                 n_engines_failed=0, highlift_config='clean',
                 lg_down=0, h_ground=0):
    '''
    Mach: float -> Freestream Mach number.
    
    altitude: float -> Flight altitude [meters].
    
    CL: float -> Lift coefficient
    
    n_engines_failed: integer -> number of engines failed. Windmilling drag is
                                 added here. This number should be less than the
                                 total number of engines.
    
    highlift_config: 'clean', 'takeoff', or 'landing' -> Configuration of high-lift devices
    
    lg_down: 0 or 1 -> 0 for retraced landing gear or 1 for extended landing gear
    
    h_ground: float -> Distance between wing and the ground for ground effect [m].
                       Use 0 for no ground effect.
    '''

    # Unpacking dictionary
    S_w = airplane['inputs']['S_w']
    AR_w = airplane['inputs']['AR_w']
    cr_w = airplane['geometry']['cr_w']
    ct_w = airplane['geometry']['ct_w']
    taper_w = airplane['inputs']['taper_w']
    sweep_w = airplane['inputs']['sweep_w']
    tcr_w = airplane['inputs']['tcr_w']
    tct_w = airplane['inputs']['tct_w']
    b_w = airplane['geometry']['b_w']
    cm_w = airplane['geometry']['cm_w']
    
    clmax_w = airplane['inputs']['clmax_w']
    k_korn = airplane['inputs']['k_korn']

    S_h = airplane['geometry']['S_h']
    cr_h = airplane['geometry']['cr_h']
    ct_h = airplane['geometry']['ct_h']
    taper_h = airplane['inputs']['taper_h']
    sweep_h = airplane['inputs']['sweep_h']
    tcr_h = airplane['inputs']['tcr_h']
    tct_h = airplane['inputs']['tct_h']
    b_h = airplane['geometry']['b_h']
    cm_h = airplane['geometry']['cm_h']
    
    S_v = airplane['geometry']['S_v']
    cr_v = airplane['geometry']['cr_v']
    ct_v = airplane['geometry']['ct_v']
    taper_v = airplane['inputs']['taper_v']
    sweep_v = airplane['inputs']['sweep_v']
    tcr_v = airplane['inputs']['tcr_v']
    tct_v = airplane['inputs']['tct_v']
    b_v = airplane['geometry']['b_v']
    cm_v = airplane['geometry']['cm_v']
    
    L_f = airplane['inputs']['L_f']
    D_f = airplane['inputs']['D_f']
    
    L_n = airplane['inputs']['L_n']
    D_n = airplane['inputs']['D_n']
    
    x_nlg = airplane['inputs']['x_nlg'] # This is only used to check if we have LG
    
    n_engines = airplane['inputs']['n_engines']
    n_engines_under_wing = airplane['inputs']['n_engines_under_wing']

    flap_type = airplane['inputs']['flap_type']
    c_flap_c_wing = airplane['inputs']['c_flap_c_wing']
    b_flap_b_wing = airplane['inputs']['b_flap_b_wing']
    
    slat_type = airplane['inputs']['slat_type']
    c_slat_c_wing = airplane['inputs']['c_slat_c_wing']
    b_slat_b_wing = airplane['inputs']['b_slat_b_wing']
    
    k_exc_drag = airplane['inputs']['k_exc_drag']

    has_winglet = airplane['inputs']['winglet']
    
    # Default rugosity value (smooth paint from Raymer Tab 12.5)
    rugosity = 0.634e-5

    # @REMOVE

    ### VISCOUS DRAG

    # Wetted areas from Torenbeek's Appendix B

    ### WING

    # Compute the wing planform area hidden by the fuselage
    D_f_b_wing = D_f/b_w
    S_hid_S_wing = D_f_b_wing*(2-D_f_b_wing*(1-taper_w))/(1+taper_w)

    # Exposed Area
    Sexp = S_w*(1 - S_hid_S_wing)

    # Wetted Area
    tau = tcr_w/tct_w
    Swet_w = 2*Sexp*(1 + 0.25*tcr_w*(1 + tau*taper_w)/(1 + taper_w))

    # Winglet contribution (used dimensions from Raymer 5th ed, Fig. 7.34 as reference)
    if has_winglet:
        tcr_winglet = 0.8
        tau = 1.0
        taper_winglet = 0.21
        Sexp_winglet = 2*(ct_w + taper_winglet*ct_w)*ct_w/2 # Assume height equal to wing tip chord
        Swet_w = Swet_w + 2*Sexp_winglet*(1 + 0.25*tcr_winglet*(1 + tau*taper_winglet)/(1 + taper_winglet))

        # Effective aspect ratio (Raymer Eq. 12.11)
        AR_eff = AR_w*1.2

    else:

        # Effective aspect ratio does not change
        AR_eff = AR_w

    # Friction coefficient
    Cf_w = Cf_calc(Mach, altitude,
                   length = cm_w,
                   rugosity = rugosity,
                   k_lam = 0.05)

    # Mean thickness
    tcm_w = 0.25*tcr_w + 0.75*tct_w

    # Form factor
    FF_w = FF_surface(Mach, tcm_w, sweep_w, b_w/2, cr_w, ct_w)

    # Interference factor
    Q_w = 1.0

    # Drag coefficient
    CD0_w = Cf_w*FF_w*Q_w*Swet_w/S_w

    ### HORIZONTAL TAIL

    # Exposed Area
    Sexp = S_h

    # Wetted Area
    tau = tcr_h/tct_h
    Swet_h = 2*Sexp*(1 + 0.25*tcr_h*(1 + tau*taper_h)/(1 + taper_h))

    # Friction coefficient
    Cf_h = Cf_calc(Mach, altitude,
                   length = cm_h,
                   rugosity = rugosity,
                   k_lam = 0.05)

    # Mean thickness
    tcm_h = 0.25*tcr_h + 0.75*tct_h

    # Form factor
    FF_h = FF_surface(Mach, tcm_h, sweep_h, b_h/2, cr_h, ct_h)

    # Interference factor
    Q_h = 1.04

    # Drag coefficient
    CD0_h = Cf_h*FF_h*Q_h*Swet_h/S_w

    ### VERTICAL TAIL

    # Exposed Area
    Sexp = S_v

    # Wetted Area
    tau = tcr_v/tct_v
    Swet_v = 2*Sexp*(1 + 0.25*tcr_v*(1 + tau*taper_v)/(1 + taper_v))

    # Friction coefficient
    Cf_v = Cf_calc(Mach, altitude,
                   length = cm_v,
                   rugosity = rugosity,
                   k_lam = 0.05)

    # Mean thickness
    tcm_v = 0.25*tcr_v + 0.75*tct_v

    # Form factor
    FF_v = FF_surface(Mach, tcm_v, sweep_v, b_v, cr_v, ct_v)

    # Interference factor
    Q_v = 1.04

    # Drag coefficient
    CD0_v = Cf_v*FF_v*Q_v*Swet_v/S_w

    ### FUSELAGE

    # Wetted area
    lambda_fus = L_f/D_f
    Swet_f = np.pi*D_f*L_f*(1 - 2/lambda_fus)**(2.0/3.0)*(1 + 1/lambda_fus**2)

    # Friction coefficient
    Cf_f = Cf_calc(Mach, altitude,
                   length = L_f,
                   rugosity = rugosity,
                   k_lam = 0.05)

    # Form factor
    FF_f = 1 + 60/lambda_fus**3 + lambda_fus/400

    # Interference factor
    Q_f = 1.0

    # Drag coefficient
    CD0_f = Cf_f*FF_f*Q_f*Swet_f/S_w

    ### NACELLE

    # Wetted area (where we take the number of nacelles into account)
    Swet_n = n_engines*np.pi*D_n*L_n

    # Friction coefficient
    Cf_n = Cf_calc(Mach, altitude,
                   length = L_n,
                   rugosity = rugosity,
                   k_lam = 0.05)

    # Form factor
    lambda_n = L_n/D_n
    FF_n = 1 + 0.35/lambda_n

    # Interference factor
    Q_n = 1.2

    # Drag coefficient
    CD0_n = Cf_n*FF_n*Q_n*Swet_n/S_w

    # Total wetted area
    Swet = Swet_w + Swet_h + Swet_v + Swet_f + Swet_n

    # Clean configuration parasite drag coefficient
    CD0_clean = CD0_w + CD0_h + CD0_v + CD0_f + CD0_n

    ### INDUCED

    # Nita and Scholz method from:
    # Estimating the Oswald Factor from Basic Aircraft Geometrical Parameters
    # Hamburg University of Applied Sciences
    # https://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PUB_DLRK_12-09-10.pdf
    # Note that there is a minus sign missing from the exponent in Eq. 37
    delta_taper = -0.357+0.45*np.exp(-0.0375*sweep_w*180/np.pi)
    taper_opt = taper_w - delta_taper
    f_taper = 0.0524*taper_opt**4 - 0.15*taper_opt**3 + 0.1659*taper_opt**2 - 0.0706*taper_opt + 0.0119
    e_theo = 1/(1 + f_taper*AR_eff)
    kem = 1/(1+0.12*Mach**6) # Here I used Howe's Mach correction since Scholz's method decays too fast with Mach
    kef = 1-2*(D_f/b_w)**2
    
    # Clean wing induced drag (before the effect of hogh-lift devices)
    e_clean = e_theo*kef*kem*0.873

    ### WAVE DRAG (Korn Equation)

    if Mach > 0.4:

        sweep_50 = change_sweep(0.25, 0.50, sweep_w, b_w/2, cr_w, ct_w)

        Mach_dd = k_korn/np.cos(sweep_50) - tcm_w/np.cos(sweep_50)**2 - CL/10/np.cos(sweep_50)**3
        Mach_crit = Mach_dd - (0.1/80)**(1/3)

        CDwave = 20*max(0, Mach - Mach_crit)**4
        
    else:
        CDwave = 0.0

    ### HIGH LIFT DEVICES

    ### Clean wing CLmax (Raymer Eq. 5.7)
    CLmax_clean = 0.9*clmax_w*np.cos(sweep_w)

    # Factor to adjust the effect of high-lift devices on each flight phase
    if highlift_config == 'clean':

        # Factor to turn off the increase in CLmax since the device is not deployed
        lift_factor = 0.0

    elif highlift_config == 'takeoff':

        # Factor to adjust CLmax contribution for intermediate deflections at takeoff.
        # Based on ratios of Torenbeek's Tab 7.2 and comments under Raymer Fig 5.3.
        lift_factor = 0.75

    elif highlift_config == 'approach':

        # Factor to adjust CLmax contribution for intermediate deflections at takeoff.
        # Based on FAR 25.121d, that says that approach Vs should be 110% of landing Vs
        # This leads to CLmax,approach = CLmax,landing/1.21
        lift_factor = 1/1.21

    elif highlift_config == 'landing':

        # Use full CLmax available
        lift_factor = 1.0

    ### Flaps deflection
    if flap_type is not None:

        # Compute flapped area (Raymer Fig. 12.21)
        # Here we consider a trapezoidal wing.
        # The second term is the subtraction of the wing area inside the fuselage
        S_flap_S_wing = (b_flap_b_wing*(2-b_flap_b_wing*(1-taper_w)))/(1+taper_w) - S_hid_S_wing

        # Sweep at flap hinge line
        sweep_flap = change_sweep(0.25, 1-c_flap_c_wing, sweep_w, b_w/2, cr_w, ct_w)
        
        # Take coefficients depending on the flap type
        # dclmax - Raymer Tab. 12.2
        # Fflap - Based on Howe Eq. 6.15b (I added the plain flap value to fit lower bound of Torenbeek Tab. 7.2)

        if flap_type == 'plain':
            dclmax = 0.9
            Fflap = 0.9
                                                
        elif flap_type == 'single slotted':
            dclmax = 1.3*(1+c_flap_c_wing) # Assumed that flap extends until trailing edge
            Fflap = 1.0
                        
        elif flap_type == 'double slotted':
            dclmax = 1.6*(1+c_flap_c_wing) # Assumed that flap extends until trailing edge
            Fflap = 1.2
                        
        elif flap_type == 'triple slotted':
            dclmax = 1.9*(1+c_flap_c_wing) # Assumed that flap extends until trailing edge
            Fflap = 1.5

        # Here we compute several factors that depend on the flight phase.
        if highlift_config == 'clean':

            # There is no increment in parasite nor induced drag in this case
            CD0_flap = 0.0
            delta_e_flap = 0.0

        elif highlift_config == 'takeoff':

            # Parasite drag contribution - Howe Eq. 6.15b
            CD0_flap = (0.03*Fflap - 0.004)/AR_eff**0.33

            # Roskam Part I Tab. 3.6 suggests a deterioration of the Oswald's factor due to
            # high lift devices. We apply the average offset here.
            delta_e_flap = -0.05

        elif highlift_config == 'approach':

            # Interpolation fraction between takeoff and landing
            ww = 0.3

            # Parasite drag contribution - Howe Eq. 6.15b
            CD0_flap = ((1-ww)*(0.03*Fflap - 0.004)+ww*0.12*Fflap)/AR_eff**0.33

            # Roskam Part I Tab. 3.6 suggests a deterioration of the Oswald's factor due to
            # high lift devices. We apply the average offset here.
            delta_e_flap = (1-ww)*(-0.05) + ww*(-0.10)
        
        elif highlift_config == 'landing':

            # Parasite drag contribution - Howe Eq. 6.16a
            # I adjusted the first factor so that drag values fall within the range
            # of Roskam Part I Tab. 3.6. The original values given by Howe's regression
            # seemed to overestimate the drag.
            CD0_flap = 0.12*Fflap/AR_eff**0.33

            # Roskam Part I Tab. 3.6 suggests a deterioration of the Oswald's factor due to
            # high lift devices. We apply the average offset here.
            delta_e_flap = -0.10

        # Chord factor to consider flap chord fractions different from the reference
        # value assumed as 30% based on DATCOM Fig. 6.1.1.3-12b
        chord_factor = c_flap_c_wing/0.30

        # Raymer Eq 12.21
        deltaCLmax_flap = 0.9*dclmax*S_flap_S_wing*np.cos(sweep_flap)*lift_factor*chord_factor

    else:
        CD0_flap = 0.0
        delta_e_flap = 0.0
        deltaCLmax_flap = 0.0

    ### Slats deflection
    if slat_type is not None:

        # Compute flapped area (Raymer Fig. 12.21)
        # Here we consider a trapezoidal wing.
        # The second term is the subtraction of the wing area inside the fuselage
        D_f_b_wing = D_f/b_w
        S_slat_S_wing = (b_slat_b_wing*(2-b_slat_b_wing*(1-taper_w)) - D_f_b_wing*(2-D_f_b_wing*(1-taper_w)))/(1+taper_w)

        sweep_slat = change_sweep(0.25, c_slat_c_wing, sweep_w, b_w/2, cr_w, ct_w)

        if slat_type == 'fixed':
            dclmax = 0.2
        elif slat_type == 'flap':
            dclmax = 0.3
        elif slat_type == 'kruger':
            dclmax = 0.3
        elif slat_type == 'slat':
            dclmax = 0.4*(1+c_slat_c_wing)

        # Chord factor to consider slat chord fractions different from the reference
        # value assumed as 15%.
        chord_factor = c_slat_c_wing/0.15

        # Raymer Eq 12.21 - I added the chord_factor to penalize smaller slats
        deltaCLmax_slat = 0.9*dclmax*S_slat_S_wing*np.cos(sweep_slat)*lift_factor*chord_factor # Raymer Eq 12.21
        
        # Parasite drag contribution (Roskam Part VI Eq. 4.72)
        CD0_slat = CD0_w*c_slat_c_wing*np.cos(sweep_w)*S_slat_S_wing*lift_factor

    else:
        CD0_slat = 0.0
        deltaCLmax_slat = 0.0

    # Maximum lift
    CLmax = CLmax_clean + deltaCLmax_flap + deltaCLmax_slat

    # Induced drag adjustment
    ee = e_clean + delta_e_flap

    # INDUCED DRAG

    # Induced drag term
    K = 1/np.pi/AR_eff/ee

    ### GROUND EFFECT
    if h_ground > 0:
        aux = 33*(h_ground/b_w)**1.5
        Kge = aux/(1+aux) # Raymer Eq. 12.61
        K = K*Kge

    # Induced drag
    CDind = K*CL**2

    ### Landing gear 
    if x_nlg is not None: # Check if we have a LG

        # Assume average value from Roskam Tab 3.6 and from Raymer page 142
        CD0_lg = lg_down*0.0200

    else:
        CD0_lg = 0.0

    ### Windmill engine
    CD0_wdm = n_engines_failed*0.3*np.pi/4*D_n**2/S_w # Raymer Eq 12.40

    # Add all parasite drag values found so far
    CD0 = (CD0_clean + CD0_flap + CD0_slat + CD0_lg + CD0_wdm)/(1-k_exc_drag)

    ### Excrescence
    CD0_exc = CD0*k_exc_drag

    # Total drag
    CD = CD0 + CDind + CDwave

    # @REMOVE

    # Create a drag breakdown dictionary
    dragDict = {
                'Mach': Mach,
                'altitude': altitude,
                'CL': CL,
                'n_engines_failed': n_engines_failed,
                'highlift_config': highlift_config,
                'lg_down': lg_down,
                'h_ground': h_ground,
                'CD': CD,
                'CD0_w': CD0_w,
                'CD0_h': CD0_h,
                'CD0_v': CD0_v,
                'CD0_f': CD0_f,
                'CD0_n': CD0_n,
                'CD0_flap' : CD0_flap,
                'CD0_slat' : CD0_slat,
                'CD0_lg' : CD0_lg,
                'CD0_wdm' : CD0_wdm,
                'CD0_exc' : CD0_exc,
                'CD0' : CD0,
                'CDind' : CDind,
                'CDwave' : CDwave,
                'CLmax_clean' : CLmax_clean,
                'deltaCLmax_flap' : deltaCLmax_flap,
                'deltaCLmax_slat' : deltaCLmax_slat,
                'CLmax' : CLmax,
                'K' : K,
                'e' : ee,
                'Swet' : Swet}

    # Update dictionary
    airplane['aerodynamics'] = {}
    airplane['aerodynamics']['Swet_f'] = Swet_f
    airplane['aerodynamics']['AR_eff'] = AR_eff
    airplane['aerodynamics']['dragDict'] = dragDict

    return CD, CLmax, dragDict
    
#----------------------------------------

def Cf_calc(Mach, altitude, length, rugosity, k_lam, deltaISA=0):
    '''
    This function computes the flat plate friction coefficient
    for a given Reynolds number while taking transition into account

    k_lam: float -> Fraction of the length (from 0 to 1) where
                    transition occurs
    '''
    
    # Dados atmosféricos
    atm_data = atmosphere(altitude, deltaISA)
    aa = atm_data['speed_of_sound']
    rho = atm_data['density']
    mi = atm_data['dyn_viscosity']

    # Velocidade
    v = aa*Mach

    # Reynolds na transição
    Re_conv = rho*v*k_lam*length/mi
    if Mach < (38.21/44.62)**(1/1.16):
        Re_rug = 38.21*(k_lam*length/rugosity)**1.053 # Raymer Eq. 12.28
    else:
        Re_rug = 44.62*(k_lam*length/rugosity)**1.053*Mach**1.16 # Raymer Eq. 12.29
    Re_trans = min(Re_conv, Re_rug)

    # Reynolds no fim
    Re_conv = rho*v*length/mi
    if Mach < (38.21/44.62)**(1/1.16):
        Re_rug = 38.21*(length/rugosity)**1.053 # Raymer Eq. 12.28
    else:
        Re_rug = 44.62*(length/rugosity)**1.053*Mach**1.16 # Raymer Eq. 12.29
    Re_fim = min(Re_conv, Re_rug)

    # Coeficientes de fricção
    # Laminar na transição
    Cf1 = 1.328/np.sqrt(Re_trans)

    # Turbulento na transição
    Cf2 = 0.455/(np.log10(Re_trans)**2.58*(1+0.144*Mach**2)**0.65)

    # Turbulento no fim
    Cf3 = 0.455/(np.log10(Re_fim)**2.58*(1+0.144*Mach**2)**0.65)

    # Média
    Cf = (Cf1 - Cf2)*k_lam + Cf3

    return Cf

#----------------------------------------

def FF_surface(Mach, tcm, sweep, panel_length, cr, ct, x_c_max_tc=0.4):
    '''
    This function computes the form factor for lifting surfaces

    INPUTS

    tcm: float -> Mean thickness/chord ratio of the surface (usually weighted as 25% of the root and 75% of the tip)
    sweep: float -> Quarter-chord sweep angle [rad]
    panel_length: float -> Spanwise distance between the root and tip sections.
    cr: float -> Root chord
    ct: float -> Tip chord
    x_c_max_tc: float -> Chord fraction with maximum thickness
    '''

    # Sweep at maximum thickness position
    sweep_maxtc = change_sweep(0.25, x_c_max_tc, sweep, panel_length, cr, ct)

    # Form factor
    FF = 1.34*Mach**0.18*np.cos(sweep_maxtc)**0.28*(1 + 0.6*tcm/x_c_max_tc + 100*(tcm)**4)

    return FF