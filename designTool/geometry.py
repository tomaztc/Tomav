# GENERAL IMPORTS
import numpy as np

#========================================

def geometry(airplane):

    # Unpack dictionary
    S_w = airplane['inputs']['S_w']
    AR_w = airplane['inputs']['AR_w']
    taper_w = airplane['inputs']['taper_w']
    sweep_w = airplane['inputs']['sweep_w']
    dihedral_w = airplane['inputs']['dihedral_w']
    xr_w = airplane['inputs']['xr_w']
    zr_w = airplane['inputs']['zr_w']
    Cht = airplane['inputs']['Cht']
    AR_h = airplane['inputs']['AR_h']
    taper_h = airplane['inputs']['taper_h']
    sweep_h = airplane['inputs']['sweep_h']
    dihedral_h = airplane['inputs']['dihedral_h']
    Lc_h = airplane['inputs']['Lc_h']
    zr_h = airplane['inputs']['zr_h']
    Cvt = airplane['inputs']['Cvt']
    AR_v = airplane['inputs']['AR_v']
    taper_v = airplane['inputs']['taper_v']
    sweep_v = airplane['inputs']['sweep_v']
    Lb_v = airplane['inputs']['Lb_v']
    zr_v = airplane['inputs']['zr_v']

    ### ADD CODE FROM SECTION 3.1 HERE ###

    # @REMOVE

    ### WING

    # Wingspan
    b_w = np.sqrt(S_w*AR_w)

    # Root chord
    cr_w = 2*S_w/(1+taper_w)/b_w

    # Tip chord
    ct_w = taper_w*cr_w

    # Tip leading edge position
    yt_w = b_w/2
    xt_w = xr_w + (cr_w - ct_w)/4 + yt_w*np.tan(sweep_w)
    zt_w = zr_w + yt_w*np.tan(dihedral_w)
    
    # Mean aerodynamic chord
    cm_w = 2*cr_w/3*(1+taper_w + taper_w**2)/(1+taper_w)

    # Mean aerodynamic chord leading edge position
    ym_w = b_w/6*(1+2*taper_w)/(1+taper_w)
    xm_w = xr_w + (cr_w - cm_w)/4 + ym_w*np.tan(sweep_w)
    zm_w = zr_w + ym_w*np.tan(dihedral_w)
    
    ### HORIZONTAL TAIL

    # Lever arm
    L_h = Lc_h*cm_w

    # Area
    S_h = Cht*S_w*cm_w/L_h

    # Wingspan
    b_h = np.sqrt(S_h*AR_h)

    # Root chord
    cr_h = 2*S_h/(1+taper_h)/b_h

    # Tip chord
    ct_h = taper_h*cr_h

    # Mean aerodynamic chord
    cm_h = 2*cr_h/3*(1+taper_h + taper_h**2)/(1+taper_h)

    # Mean aerodynamic chord leading edge position
    ym_h = b_h/6*(1+2*taper_h)/(1+taper_h)
    xm_h = xm_w + cm_w/4 + L_h - cm_h/4
    zm_h = zr_h + ym_h*np.tan(dihedral_h)

    # Root leading edge position
    xr_h = xm_h + (cm_h - cr_h)/4 - ym_h*np.tan(sweep_h)

    # Tip leading edge position
    yt_h = b_h/2
    xt_h = xr_h + (cr_h - ct_h)/4 + yt_h*np.tan(sweep_h)
    zt_h = zr_h + yt_h*np.tan(dihedral_h)

    ### VERTICAL TAIL

    # Lever arm
    L_v = Lb_v*b_w

    # Area
    S_v = Cvt*S_w*b_w/L_v

    # Wingspan
    b_v = np.sqrt(S_v*AR_v)

    # Root chord
    cr_v = 2*S_v/(1+taper_v)/b_v

    # Tip chord
    ct_v = taper_v*cr_v

    # Mean aerodynamic chord
    cm_v = 2*cr_v/3*(1+taper_v + taper_v**2)/(1+taper_v)

    # Mean aerodynamic chord leading edge position
    zm_v = zr_v + b_v/3*(1+2*taper_v)/(1+taper_v)
    xm_v = xm_w + cm_w/4 + L_v - cm_v/4

    # Root leading edge position
    xr_v = xm_v + (cm_v - cr_v)/4 - (zm_v-zr_v)*np.tan(sweep_v)

    # Tip leading edge position
    zt_v = zr_v + b_v
    xt_v = xr_v + (cr_v - ct_v)/4 + (zt_v-zr_v)*np.tan(sweep_v)

    # @REMOVE

    # Update dictionary with new results
    airplane['geometry'] = {}
    airplane['geometry']['b_w'] = b_w
    airplane['geometry']['cr_w'] = cr_w
    airplane['geometry']['xt_w'] = xt_w
    airplane['geometry']['yt_w'] = yt_w
    airplane['geometry']['zt_w'] = zt_w
    airplane['geometry']['ct_w'] = ct_w
    airplane['geometry']['xm_w'] = xm_w
    airplane['geometry']['ym_w'] = ym_w
    airplane['geometry']['zm_w'] = zm_w
    airplane['geometry']['cm_w'] = cm_w
    airplane['geometry']['S_h'] = S_h
    airplane['geometry']['b_h'] = b_h
    airplane['geometry']['xr_h'] = xr_h
    airplane['geometry']['cr_h'] = cr_h
    airplane['geometry']['xt_h'] = xt_h
    airplane['geometry']['yt_h'] = yt_h
    airplane['geometry']['zt_h'] = zt_h
    airplane['geometry']['ct_h'] = ct_h
    airplane['geometry']['xm_h'] = xm_h
    airplane['geometry']['ym_h'] = ym_h
    airplane['geometry']['zm_h'] = zm_h
    airplane['geometry']['cm_h'] = cm_h
    airplane['geometry']['S_v'] = S_v
    airplane['geometry']['b_v'] = b_v
    airplane['geometry']['xr_v'] = xr_v
    airplane['geometry']['cr_v'] = cr_v
    airplane['geometry']['xt_v'] = xt_v
    airplane['geometry']['zt_v'] = zt_v
    airplane['geometry']['ct_v'] = ct_v
    airplane['geometry']['xm_v'] = xm_v
    airplane['geometry']['zm_v'] = zm_v
    airplane['geometry']['cm_v'] = cm_v

    # All variables are stored in the dictionary.
    # There is no need to return anything
    return None

#----------------------------------------

def change_sweep(x,y,sweep_x,panel_length,chord_root,chord_tip):

    '''
    This function converts sweep computed at chord fraction x into
    sweep measured at chord fraction y
    (x and y should be between 0 (leading edge) and 1 (trailing edge).
    panel_length is the spanwise distance between the root and tip section.
    '''

    sweep_y = np.arctan(np.tan(sweep_x)+(x-y)*(chord_root-chord_tip)/panel_length)

    return sweep_y

#----------------------------------------

def control_surface_area_fraction(xf, y1, y2, taper):
    '''
    xf: flap_chord/wing_chord
    y1: spanwise fraction where flap starts
    y2: spanwise fraction where flap ends
    taper: taper ratio of the wing
    '''

    S_flap_S_wing = xf/(1+taper)*(y2*(2-y2*(1-taper)) - y1*(2-y1*(1-taper)))

    return S_flap_S_wing