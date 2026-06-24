# GENERAL IMPORTS
import numpy as np
from .constants import gravity
from .auxiliary import atmosphere
from .aerodynamics import aerodynamics
from .propulsion import engineTSFC
from .weight import weight

#========================================

def performance(W0, W_cruise, airplane,
                landing_method='Torenbeek'):

    '''
    This function computes the required thrust and wing areas
    required to meet takeoff, landing, climb, and cruise requirements.

    OUTPUTS:
    T0: real -> Total thrust required to meet all mission phases
    deltaS_wlan: real -> Wing area margin for landing. This value should be positive
                         for a feasible landing.
    '''

    # Unpacking dictionary
    S_w = airplane['inputs']['S_w']
    
    n_engines = airplane['inputs']['n_engines']
    
    h_ground = airplane['inputs']['h_ground']
    
    altitude_takeoff = airplane['inputs']['altitude_takeoff']
    distance_takeoff = airplane['inputs']['distance_takeoff']
    deltaISA_takeoff = airplane['inputs']['deltaISA_takeoff']
    
    altitude_landing = airplane['inputs']['altitude_landing']
    distance_landing = airplane['inputs']['distance_landing']
    deltaISA_landing = airplane['inputs']['deltaISA_landing']
    MLW_frac = airplane['inputs']['MLW_frac']
    
    altitude_cruise = airplane['inputs']['altitude_cruise']
    Mach_cruise = airplane['inputs']['Mach_cruise']

    altitude_maxcruise = airplane['inputs']['altitude_maxcruise']
    Mach_maxcruise = airplane['inputs']['Mach_maxcruise']

    # @REMOVE

    ### TAKEOFF

    # Compute air density at takeoff altitude
    atm_data = atmosphere(altitude_takeoff, deltaISA_takeoff)
    rho = atm_data['density']

    # density ratio
    sigma = rho/1.225

    # Takeoff aerodynamics
    _, CLmaxTO, _ = aerodynamics(airplane, Mach=0.2, altitude=altitude_takeoff,
                                 CL=0.5, n_engines_failed=0, highlift_config='takeoff',
                                 lg_down=1, h_ground=h_ground)
    
    T0W0 = 0.2387/sigma/CLmaxTO/distance_takeoff*W0/S_w

    T0_to = T0W0*W0

    ### LANDING

    # Compute air density at landing altitude
    atm_data = atmosphere(altitude_landing, deltaISA_landing)
    rho = atm_data['density']

    # Landing aerodynamics
    _, CLmaxLD, _ = aerodynamics(airplane, Mach=0.2, altitude=altitude_landing,
                                 CL=0.5, n_engines_failed=0, highlift_config='landing',
                                 lg_down=1, h_ground=h_ground)

    if landing_method == 'Roskam':
    
        # Landing Field Length (Roskam)

        # Approach speed (Roskam adapted to SI)
        Va = 1.701*np.sqrt(distance_landing)

        # Required stall speed
        Vs = Va/1.3

        # Required wing area
        S_wlan = 2*W0*MLW_frac/rho/Vs**2/CLmaxLD
    
    elif landing_method == 'Torenbeek':

        # Landing distance - Torenbeek (eq. 5-98, page 171)
    
        h_land = 15.3 # Screen height (meters)
        f_land = 5/3 # Landing distance safety factor (FAR Part 91)
        a_g = 0.5 # Mean deceleration factor

        x_land = 1.52/a_g + 1.69
        A_land = gravity/f_land/x_land
        B_land = -10.0*gravity*h_land/x_land

        S_wlan = W0*MLW_frac/rho/(A_land*distance_landing + B_land)/CLmaxLD

    # Compute wing area excess with respect to the landing requirement.
    # The aircraft should have deltaS_wlan >= 0 to satisfy landing.
    deltaS_wlan = S_w - S_wlan

    ### CRUISE

    # Compute air density at cruise altitude
    atm_data = atmosphere(altitude_cruise)
    a_cruise = atm_data['speed_of_sound']
    rho = atm_data['density']

    # Cruise speed
    v_cruise = Mach_cruise*a_cruise

    # Cruise CL
    CL = 2.0*W_cruise/rho/S_w/v_cruise**2

    # Cruise CD
    CD, _, _ = aerodynamics(airplane, Mach=Mach_cruise, altitude=altitude_cruise,
                            CL=CL, n_engines_failed=0, highlift_config='clean',
                            lg_down=0, h_ground=0)

    # Cruise traction
    T = 0.5*rho*v_cruise**2*S_w*CD

    # Cruise correction factor
    _, kT = engineTSFC(Mach_cruise, altitude_cruise, airplane)

    # Corrected thrust
    T0_cruise = T/kT

    ### HIGH-SPEED CRUISE

    # Compute air density at cruise altitude
    atm_data = atmosphere(altitude_maxcruise)
    a_maxcruise = atm_data['speed_of_sound']
    rho = atm_data['density']

    # Cruise speed
    v_maxcruise = Mach_maxcruise*a_maxcruise

    # Cruise CL
    CL = 2.0*W_cruise/rho/S_w/v_maxcruise**2

    # Cruise CD
    CD, _, _ = aerodynamics(airplane, Mach=Mach_maxcruise, altitude=altitude_maxcruise,
                            CL=CL, n_engines_failed=0, highlift_config='clean',
                            lg_down=0, h_ground=0)

    # Cruise traction
    T = 0.5*rho*v_maxcruise**2*S_w*CD

    # Cruise correction factor
    _, kT = engineTSFC(Mach_maxcruise, altitude_maxcruise, airplane)

    # Corrected thrust
    T0_maxcruise = T/kT

    ### CLIMB

    # Define standard function for climb analysis
    def climb_analysis(grad, Ks, altitude,
                       lg_down, h_ground_climb, highlift_config, n_engines_failed, Mf,
                       kT, deltaISA=0):

        '''
        kT: Thrust decay factor (e.g. use 0.94 for maximum continuous thrust)
        '''

	    # Compute air temperature and speed of sound
        atm_data = atmosphere(altitude, deltaISA)
        a = atm_data['speed_of_sound']
        rho = atm_data['density']

        # Dummy run to get CLmax
        _, CLmax, _ = aerodynamics(airplane, Mach=0.2, altitude=altitude,
                                   CL=0.5, n_engines_failed=n_engines_failed, highlift_config=highlift_config,
                                   lg_down=lg_down, h_ground=h_ground_climb)

        # Get climb CL
        CL = CLmax/Ks**2
        
        # Get climb speed
        Vclimb = np.sqrt(2*W0*Mf/rho/S_w/CL)
        
        # Compute sound speed and Mach number
        Mach = Vclimb/a

        # Get corresponding CD # ESCREVER
        CD, _, _ = aerodynamics(airplane, Mach=Mach, altitude=altitude,
                                CL=CL, n_engines_failed=n_engines_failed, highlift_config=highlift_config,
                                lg_down=lg_down, h_ground=h_ground_climb)

        # Check number of failed engines
        if n_engines_failed >= n_engines:
            print('Warning: number of failed engines is equal or greater than the number of engines')
            print('We will force n_engines_failed = n_engines-1')
            n_engines_failed = n_engines - 1

        # Compute T/W
        TW = n_engines/(n_engines-n_engines_failed)*(grad + CD/CL)

        # Compute required traction
        T0 = TW*W0*Mf/kT

        return T0

    # FAR 25.111
    if n_engines <= 2:
        grad = 0.012
    elif n_engines == 3:
        grad = 0.015
    elif n_engines == 4:
        grad = 0.017
    Ks = 1.2
    altitude = altitude_takeoff
    lg_down = 0
    h_ground_climb = h_ground
    highlift_config = 'takeoff'
    n_engines_failed = 1
    Mf = 1.0
    kT = 1.0
    deltaISA = deltaISA_takeoff
    T0_1 = climb_analysis(grad, Ks, altitude,
                          lg_down, h_ground_climb, highlift_config, n_engines_failed, Mf,
                          kT, deltaISA)

    # FAR 25.121a
    if n_engines <= 2:
        grad = 0.000
    elif n_engines == 3:
        grad = 0.003
    elif n_engines == 4:
        grad = 0.005
    Ks = 1.1
    altitude = altitude_takeoff
    lg_down = 1
    h_ground_climb = h_ground
    highlift_config = 'takeoff'
    n_engines_failed = 1
    Mf = 1.0
    kT = 1.0
    deltaISA = deltaISA_takeoff
    T0_2 = climb_analysis(grad, Ks, altitude,
                          lg_down, h_ground_climb, highlift_config, n_engines_failed, Mf,
                          kT, deltaISA)

    # FAR 25.121b
    if n_engines <= 2:
        grad = 0.024
    elif n_engines == 3:
        grad = 0.027
    elif n_engines == 4:
        grad = 0.030
    Ks = 1.2
    altitude = altitude_takeoff
    lg_down = 0
    h_ground_climb = 0
    highlift_config  = 'takeoff'
    n_engines_failed = 1
    Mf = 1.0
    kT = 1.0
    deltaISA = deltaISA_takeoff
    T0_3 = climb_analysis(grad, Ks, altitude,
                          lg_down, h_ground_climb, highlift_config, n_engines_failed, Mf,
                          kT, deltaISA)

    # FAR 25.121c
    if n_engines <= 2:
        grad = 0.012
    elif n_engines == 3:
        grad = 0.015
    elif n_engines == 4:
        grad = 0.017
    Ks = 1.25
    altitude = altitude_takeoff
    lg_down = 0
    h_ground_climb = 0
    highlift_config = 'clean'
    n_engines_failed = 1
    Mf = 1.0
    kT = 0.94
    deltaISA = deltaISA_takeoff
    T0_4 = climb_analysis(grad, Ks, altitude,
                          lg_down, h_ground_climb, highlift_config, n_engines_failed, Mf,
                          kT, deltaISA)

    # FAR 25.119
    grad = 0.032
    Ks = 1.30
    altitude = altitude_landing
    lg_down = 1
    h_ground_climb = 0
    highlift_config = 'landing'
    n_engines_failed = 0
    Mf = MLW_frac
    kT = 1.0
    deltaISA = deltaISA_landing
    T0_5 = climb_analysis(grad, Ks, altitude,
                          lg_down, h_ground_climb, highlift_config, n_engines_failed, Mf,
                          kT, deltaISA)

    # FAR 25.121d
    if n_engines <= 2:
        grad = 0.021
    elif n_engines == 3:
        grad = 0.024
    elif n_engines == 4:
        grad = 0.027
    Ks = 1.40
    altitude = altitude_landing
    lg_down = 0
    h_ground_climb = 0
    highlift_config = 'approach'
    n_engines_failed = 1
    Mf = MLW_frac
    kT = 1.0
    deltaISA = deltaISA_landing
    T0_6 = climb_analysis(grad, Ks, altitude,
                          lg_down, h_ground_climb, highlift_config, n_engines_failed, Mf,
                          kT, deltaISA)
    # @REMOVE
    # @REMOVE
    # RATO TAKEOFF

    if 'altitude_rato_takeoff' in airplane['inputs'].keys():
        '''
        RATO takeoff formulation from
        SILVA, Daniel Marques. Modelo de Desempenho para Avaliar o Perfil de Missão de Um
        AAM Aplicável ao Projeto Conceitual da Aeronave. 2021. Trabalho de Conclusão
        de Curso (Graduação). Instituto Tecnológico de Aeronáutica, São José dos Campos.

        This will be assigned as a wing loading requirement that could replace
        the landing requirement if it is more critical
        '''

        # Unpack optional arguments
        altitude_rato = airplane['altitude_rato_takeoff']
        ramp_angle_rato = airplane['ramp_angle_rato_takeoff']
        c_rato = airplane['c_rato_takeoff']
        prop_mass_ratio_rato = airplane['prop_mass_ratio_rato_takeoff']
        t_burn_rato = airplane['t_burn_rato_takeoff']
        Ks = airplane['Ks_rato_takeoff']

        # Compute air density at takeoff altitude
        atm_data = atmosphere(altitude_rato)
        rho = atm_data['density']
    
        # density ratio
        sigma = rho/1.225
    
        # Takeoff aerodynamics
        _, CLmaxTO, _ = aerodynamics(airplane, Mach=0.2, altitude=altitude_rato,
                                     CL=0.5, n_engines_failed=0, highlift_config='takeoff',
                                     lg_down=0, h_ground=0)
    
        # Compute required wing loading to reach desired stall speed
        # factor at the end of the burn
        WS_req = 0.5*rho*CLmaxTO/Ks*(c_rato*np.log(1-prop_mass_ratio_rato) + gravity*t_burn_rato*np.sin(ramp_angle_rato))**2
    
        deltaS_rato = S_w - W0/WS_req
        deltaS_wlan = min(deltaS_wlan, deltaS_rato)

        # Check climb gradient at the end of the burn
        # This will replace last climb requirement
        grad = np.tan(ramp_angle_rato)
        Ks = Ks
        altitude = altitude_rato
        CLmax_guess = CLmaxTO
        lg_down = 1
        h_ground_climb = 0
        highlift_config = 'takeoff'
        n_engines_failed = 0
        Mf = 1.0
        kT = 1.0
        T0_6 = climb_analysis(grad, Ks, altitude, CLmax_guess,
                              lg_down, h_ground_climb, highlift_config, n_engines_failed, Mf,
                              kT)
    # @REMOVE

    # Get the maximum required thrust with a 5% margin
    T0vec = [T0_to, T0_cruise, T0_maxcruise, T0_1, T0_2, T0_3, T0_4, T0_5, T0_6]
    T0 = 1.05*max(T0vec)
    T0req = {'Takeoff': T0_to,
             'Cruise': T0_cruise,
             'High speed cruise': T0_maxcruise,
             'FAR 25.111': T0_1,
             'FAR 25.121a': T0_2,
             'FAR 25.121b': T0_3,
             'FAR 25.121c': T0_4,
             'FAR 25.119': T0_5,
             'FAR 25.121d': T0_6}

    return T0, T0req, deltaS_wlan, CLmaxTO

#----------------------------------------

def thrust_matching(W0_guess, T0_guess, airplane, calcular_performance=False):

    # If the user provides the engine thrust, we only do one execution
    # to gather the required thrust
    if 'Tmax' in airplane['inputs']['engine']:

        T0 = airplane['inputs']['engine']['Tmax']*airplane['inputs']['n_engines']

        W0, W_empty, W_fuel, W_cruise = weight(W0_guess, T0, airplane)
        if calcular_performance:
            _, T0req, deltaS_wlan, CLmaxTO = performance(W0, W_cruise, airplane)
        else:
            # Takeoff aerodynamics
            _, CLmaxTO, _ = aerodynamics(airplane, Mach=0.2, altitude=airplane['inputs']['altitude_takeoff'],
                                        CL=0.5, n_engines_failed=0, highlift_config='takeoff',
                                        lg_down=1, h_ground=airplane['inputs']['h_ground'])

    else:
        # Set iterator
        delta = 1000

        # Loop to adjust T0
        while abs(delta) > 10:

            W0, W_empty, W_fuel, W_cruise = weight(W0_guess, T0_guess, airplane)

            T0, T0req, deltaS_wlan, CLmaxTO = performance(W0, W_cruise, airplane)

            # Compute change with respect to previous iteration
            delta = T0 - T0_guess

            # Update guesses for the next iteration
            T0_guess = T0
            W0_guess = W0

    # Update dictionary with converged values
    airplane['thrust_matching'] = {}
    airplane['thrust_matching']['W0'] = W0
    airplane['thrust_matching']['W_empty'] = W_empty
    airplane['thrust_matching']['W_fuel'] = W_fuel
    airplane['thrust_matching']['T0'] = T0
    if calcular_performance:
        airplane['thrust_matching']['T0req'] = T0req
        airplane['thrust_matching']['deltaS_wlan'] = deltaS_wlan
    airplane['thrust_matching']['CLmaxTO'] = CLmaxTO

    # Return
    return None