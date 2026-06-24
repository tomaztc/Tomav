# GENERAL IMPORTS
import numpy as np

#========================================

def landing_gear(airplane):

    # Unpack dictionary
    x_nlg = airplane['inputs']['x_nlg']
    x_mlg = airplane['inputs']['x_mlg']
    y_mlg = airplane['inputs']['y_mlg']
    z_lg = airplane['inputs']['z_lg']
    xcg_fwd = airplane['balance']['xcg_fwd']
    xcg_aft = airplane['balance']['xcg_aft']
    x_tailstrike = airplane['inputs']['x_tailstrike']
    z_tailstrike = airplane['inputs']['z_tailstrike']

    # @REMOVE

    # Check if there is a landing gear
    if x_nlg is not None:

        # Weight fractions on NLG for both load cases
        frac_nlg_fwd = (x_mlg-xcg_fwd)/(x_mlg-x_nlg)
        frac_nlg_aft = (x_mlg-xcg_aft)/(x_mlg-x_nlg)
    
        # Tipback angle (for now assume that CG is along fuselage axis)
        alpha_tipback = np.arctan((x_mlg - xcg_aft)/(-z_lg))
    
        # Tailstrike angle
        alpha_tailstrike = np.arctan((z_tailstrike - z_lg)/(x_tailstrike - x_mlg))
    
        # Overturn angle
        sgl = (xcg_fwd - x_nlg)*y_mlg/np.sqrt((x_mlg - x_nlg)**2 + y_mlg**2)
        phi_overturn = np.arctan(-z_lg/sgl)
        
    else:
        
        # Add dummy data
        frac_nlg_fwd = None
        frac_nlg_aft = None
        alpha_tipback = None
        alpha_tailstrike = None
        phi_overturn = None

    # @REMOVE

    # Update dictionary
    airplane['landing_gear'] = {}
    airplane['landing_gear']['frac_nlg_fwd'] = frac_nlg_fwd
    airplane['landing_gear']['frac_nlg_aft'] = frac_nlg_aft
    airplane['landing_gear']['alpha_tipback'] = alpha_tipback
    airplane['landing_gear']['alpha_tailstrike'] = alpha_tailstrike
    airplane['landing_gear']['phi_overturn'] = phi_overturn

    return None