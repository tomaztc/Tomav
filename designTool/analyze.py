# coding=utf-8
'''
Conceptual Aircraft Design Tool
(for PRJ-22 and AP-701 courses)

Maj. Eng. Ney Rafael Secco (ney@ita.br)
Aircraft Design Department
Aeronautics Institute of Technology

05-2025

The code uses several historical regression from
aircraft design books to make a quick initial
sizing procedure.

Generally, the user should call only the 'analyze'
function from this module.
'''

# IMPORTS
import numpy as np
from .constants import gravity
from .standard_airplane import standard_airplane
from .plots import plot_geometry
from .geometry import geometry
from .performance import thrust_matching
from .balance import balance
from .landing_gear import landing_gear

np.set_printoptions(legacy='1.25')
#========================================
# MAIN FUNCTION

def analyze(airplane = None,
            print_log = False, # Plot results on the terminal screen
            plot = False, # Generate 3D plot of the aircraft
            ):
    '''
    This is the main function that should be used for aircraft analysis.
    '''

    # Load standard airplane if none is provided
    if airplane is None:
        airplane = standard_airplane()

    # Use an average wing loading for transports
    # to estime W0_guess and T0_guess if none are provided
    if 'W0_guess' in airplane['inputs'].keys():
        W0_guess = airplane['inputs']['W0_guess']
    else:
        W0_guess = 5e3*airplane['inputs']['S_w']
    
    if 'T0_guess' in airplane['inputs'].keys():
        T0_guess = airplane['inputs']['T0_guess']
    else:
        T0_guess = 0.3*W0_guess

    # Generate geometry
    geometry(airplane)

    if plot:
        plot_geometry(airplane)

    # Converge MTOW and Takeoff Thrust
    thrust_matching(W0_guess, T0_guess, airplane, calcular_performance=True)

    # Balance analysis
    balance(airplane)

    # Landing gear design
    landing_gear(airplane)

    if print_log:
        W_empty = airplane['thrust_matching']['W_empty']
        W_fuel = airplane['thrust_matching']['W_fuel']
        W0 = airplane['thrust_matching']['W0']
        T0 = airplane['thrust_matching']['T0']
        S_w = airplane['inputs']['S_w']
        deltaS_wlan = airplane['thrust_matching']['deltaS_wlan']
        tank_excess = airplane['balance']['tank_excess']
        V_maxfuel = airplane['balance']['V_maxfuel']
        CLv = airplane['balance']['CLv']
        SM_fwd = airplane['balance']['SM_fwd']
        SM_aft = airplane['balance']['SM_aft']
        xnp = airplane['balance']['xnp']
        xm_w = airplane['geometry']['xm_w']
        cm_w = airplane['geometry']['cm_w']
        xcg_fwd = airplane['balance']['xcg_fwd']
        xcg_aft = airplane['balance']['xcg_aft']

        print('W_empty [kgf]: %d'%(W_empty/gravity))
        print('W_fuel [kgf]: %d'%(W_fuel/gravity))
        print('W0 [kgf]: %d'%(W0/gravity))
        print('T0 [kgf]: %d'%(T0/gravity))
        print('T0/W0: %.3f'%(T0/W0))
        print('W0/S [kgf/m2]: %d'%(W0/gravity/S_w))
        print('deltaS_wlan [m2]: %.1f'%(deltaS_wlan))
        print('tank_excess [%%]: %.1f'%(tank_excess*100))
        print('V_maxfuel [L]: %d'%(V_maxfuel*1000))
        print('CLv: %.3f'%(CLv))
        print('SM_fwd [%%]: %.1f'%(SM_fwd*100))
        print('SM_aft [%%]: %.1f'%(SM_aft*100))
        print('xnp [%%MAC]: %.1f'%((xnp-xm_w)/cm_w*100))
        print('xcg_fwd [%%MAC]: %.1f'%((xcg_fwd-xm_w)/cm_w*100))
        print('xcg_aft [%%MAC]: %.1f'%((xcg_aft-xm_w)/cm_w*100))
        
        if airplane['inputs']['x_nlg'] is not None:

            x_mlg = airplane['inputs']['x_mlg']
            frac_nlg_fwd = airplane['landing_gear']['frac_nlg_fwd']
            frac_nlg_aft = airplane['landing_gear']['frac_nlg_aft']
            alpha_tipback = airplane['landing_gear']['alpha_tipback']
            alpha_tailstrike = airplane['landing_gear']['alpha_tailstrike']
            phi_overturn = airplane['landing_gear']['phi_overturn']

            print('x_mlg [%%MAC]: %.1f'%((x_mlg-xm_w)/cm_w*100))
            print('frac_nlg_fwd [%%]: %.1f'%(frac_nlg_fwd*100))
            print('frac_nlg_aft [%%]: %.1f'%(frac_nlg_aft*100))
            print('alpha_tipback [deg]: %.1f'%(alpha_tipback*180.0/np.pi))
            print('alpha_tailstrike [deg]: %.1f'%(alpha_tailstrike*180.0/np.pi))
            print('phi_overturn [deg]: %.1f'%(phi_overturn*180.0/np.pi))

    # Plot again now that we have CG and NP
    if plot:
        plot_geometry(airplane)

    return airplane