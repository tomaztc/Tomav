# GENERAL IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#========================================

def plot_geometry(airplane, figname='3dview.png', az1=45, az2=-135):
    '''
    az1 and az2: degrees of azimuth and elevation for the 3d plot view
    '''

    from matplotlib.patches import Ellipse
    import mpl_toolkits.mplot3d.art3d as art3d
    from .auxiliary import lin_interp

    xr_w = airplane['inputs']['xr_w']
    zr_w = airplane['inputs']['zr_w']
    b_w = airplane['geometry']['b_w']

    tct_w = airplane['inputs']['tct_w']
    tcr_w = airplane['inputs']['tcr_w']

    cr_w = airplane['geometry']['cr_w']
    xt_w = airplane['geometry']['xt_w']
    yt_w = airplane['geometry']['yt_w']
    zt_w = airplane['geometry']['zt_w']
    ct_w = airplane['geometry']['ct_w']

    ym_w = airplane['geometry']['ym_w']
    zm_w = airplane['geometry']['zm_w']

    xr_h = airplane['geometry']['xr_h']
    zr_h = airplane['inputs']['zr_h']

    tcr_h = airplane['inputs']['tcr_h']
    tct_h = airplane['inputs']['tct_h']

    cr_h = airplane['geometry']['cr_h']
    xt_h = airplane['geometry']['xt_h']
    yt_h = airplane['geometry']['yt_h']
    zt_h = airplane['geometry']['zt_h']
    ct_h = airplane['geometry']['ct_h']
    b_h  = airplane['geometry']['b_h']

    cm_h = airplane['geometry']['cm_h']
    xm_h = airplane['geometry']['xm_h']
    ym_h = airplane['geometry']['ym_h']
    zm_h = airplane['geometry']['zm_h']

    xr_v = airplane['geometry']['xr_v']
    zr_v = airplane['inputs']['zr_v']

    tcr_v = airplane['inputs']['tcr_v']
    tct_v = airplane['inputs']['tct_v']

    cr_v = airplane['geometry']['cr_v']
    xt_v = airplane['geometry']['xt_v']
    zt_v = airplane['geometry']['zt_v']
    ct_v = airplane['geometry']['ct_v']
    b_v  = airplane['geometry']['b_v']

    L_f = airplane['inputs']['L_f']
    D_f = airplane['inputs']['D_f']
    x_n = airplane['inputs']['x_n']
    y_n = airplane['inputs']['y_n']
    z_n = airplane['inputs']['z_n']
    L_n = airplane['inputs']['L_n']
    D_n = airplane['inputs']['D_n']

    has_winglet = airplane['inputs']['winglet']

    if 'balance' in airplane:
        xcg_fwd = airplane['balance']['xcg_fwd']
        xcg_aft = airplane['balance']['xcg_aft']
        xnp = airplane['balance']['xnp']
    else:
        xcg_fwd = None
        xcg_aft = None
        xnp = None

    x_nlg = airplane['inputs']['x_nlg']
    y_nlg = 0
    z_nlg = airplane['inputs']['z_lg']
    x_mlg = airplane['inputs']['x_mlg']
    y_mlg = airplane['inputs']['y_mlg']
    z_mlg = airplane['inputs']['z_lg']
    x_tailstrike = airplane['inputs']['x_tailstrike']
    z_tailstrike = airplane['inputs']['z_tailstrike']

    flap_type = airplane['inputs']['flap_type']
    b_flap_b_wing = airplane['inputs']['b_flap_b_wing']
    c_flap_c_wing = airplane['inputs']['c_flap_c_wing']

    slat_type = airplane['inputs']['slat_type']
    b_slat_b_wing = airplane['inputs']['b_slat_b_wing']
    c_slat_c_wing = airplane['inputs']['c_slat_c_wing']

    b_ail_b_wing = airplane['inputs']['b_ail_b_wing']
    c_ail_c_wing = airplane['inputs']['c_ail_c_wing']

    c_tank_c_w = airplane['inputs']['c_tank_c_w']
    x_tank_c_w = airplane['inputs']['x_tank_c_w']
    b_tank_b_w_start = airplane['inputs']['b_tank_b_w_start']
    b_tank_b_w_end = airplane['inputs']['b_tank_b_w_end']

    ### PLOT

    #fig = plt.figure(fignum,figsize=(20, 10))
    fig = plt.figure(num=figname, constrained_layout=True)
    ax = fig.add_subplot(projection='3d', proj_type='ortho')
    # ax.set_aspect('equal')

    ax.plot([xr_w, xt_w, xt_w+ct_w, xr_w+cr_w, xt_w+ct_w, xt_w, xr_w],
            [0.0, yt_w, yt_w, 0.0, -yt_w, -yt_w, 0.0],
            [zr_w+cr_w*tcr_w/2, zt_w+ct_w*tct_w/2, zt_w+ct_w*tct_w/2, zr_w+cr_w*tcr_w/2, zt_w+ct_w*tct_w/2, zt_w+ct_w*tct_w/2, zr_w+cr_w*tcr_w/2],color='blue')

    if has_winglet:
        ttw = 0.21 # Winglet taper ratio
        ax.plot([xt_w, xt_w + (1-ttw)*ct_w, xt_w+ct_w, xt_w+ct_w, xt_w],
            [yt_w, yt_w, yt_w, yt_w, yt_w],
            [zt_w, zt_w+ct_w, zt_w+ct_w, zt_w, zt_w],color='blue')
        ax.plot([xt_w, xt_w + (1-ttw)*ct_w, xt_w+ct_w, xt_w+ct_w, xt_w],
            [-yt_w, -yt_w, -yt_w, -yt_w, -yt_w],
            [zt_w, zt_w+ct_w, zt_w+ct_w, zt_w, zt_w],color='blue')

    ax.plot([xr_h, xt_h, xt_h+ct_h, xr_h+cr_h, xt_h+ct_h, xt_h, xr_h],
            [0.0, yt_h, yt_h, 0.0, -yt_h, -yt_h, 0.0],
            [zr_h+cr_h*tcr_h/2, zt_h+ct_h*tct_h/2, zt_h+ct_h*tct_h/2, zr_h+cr_h*tcr_h/2, zt_h+ct_h*tct_h/2, zt_h+ct_h*tct_h/2, zr_h+cr_h*tcr_h/2],color='green')


    ax.plot([xr_v        , xt_v        , xt_v+ct_v   , xr_v+cr_v   , xr_v        ],
            [tcr_v*cr_v/2, tct_v*ct_v/2, tct_v*ct_v/2, tcr_v*cr_v/2, tcr_v*cr_v/2],
            [zr_v        , zt_v        , zt_v        , zr_v        , zr_v        ],\
            color='orange')

    ax.plot([ xr_v        ,  xt_v        ,  xt_v+ct_v   ,  xr_v+cr_v   ,  xr_v        ],
            [-tcr_v*cr_v/2, -tct_v*ct_v/2, -tct_v*ct_v/2, -tcr_v*cr_v/2, -tcr_v*cr_v/2],
            [ zr_v        ,  zt_v        ,  zt_v        ,  zr_v        ,  zr_v     ],\
            color='orange')



    ax.plot([0.0, L_f],
            [0.0, 0.0],
            [0.0, 0.0])
    ax.plot([x_n, x_n+L_n],
            [y_n, y_n],
            [z_n, z_n])
    ax.plot([x_n, x_n+L_n],
            [-y_n, -y_n],
            [z_n, z_n])

    # Forward CG point
    if xcg_fwd is not None:
        ax.plot([xcg_fwd], [-ym_w], [zm_w],'ko')
    
    # Rear CG point
    if xcg_aft is not None:
        ax.plot([xcg_aft], [-ym_w], [zm_w],'ko')
    
    # Neutral point
    if xnp is not None:
        ax.plot([xnp], [-ym_w], [zm_w],'x')

    # Define a parametrized fuselage by setting height and width
    # values along its axis
    # xx is non-dimensionalized by fuselage length
    # hh and ww are non-dimensionalized by fuselage diameter
    # There are 6 stations where we define the arrays:
    # nose1; nose2; nose3; cabin start; tailstrike; tail
    xx = [0.0, 1.24/41.72, 3.54/41.72, 7.55/41.72, x_tailstrike/L_f, 1.0]
    hh = [0.0, 2.27/4.0, 3.56/4.0, 1.0, 1.0, 1.07/4.0]
    ww = [0.0, 1.83/4.0, 3.49/4.0, 1.0, 1.0, 0.284/4]
    num_tot_ell = 70 # Total number of ellipses
    
    # Loop over every section
    for ii in range(len(xx)-1):
        
        # Define number of ellipses based on the section length
        num_ell = int((xx[ii+1]-xx[ii])*num_tot_ell)+1
        
        # Define arrays of dimensional positions, heights and widths
        # for the current section
        xdim = np.linspace(xx[ii], xx[ii+1], num_ell)*L_f
        hdim = np.linspace(hh[ii], hh[ii+1], num_ell)*D_f
        wdim = np.linspace(ww[ii], ww[ii+1], num_ell)*D_f
        
        # Loop over every ellipse
        for xc, hc, wc in zip(xdim, hdim, wdim):

            # Define ellipse center to make flat top at the fuselage tail
            if xc > x_tailstrike:
                yye = (D_f-hc)/2
            else:
                yye = 0

            p = Ellipse((0, yye), wc, hc, angle=0,
                        facecolor = 'none', edgecolor = 'k', lw=1.0)
            ax.add_patch(p)
            art3d.pathpatch_2d_to_3d(p, z=xc, zdir="x")


    #____________________________________________________________
    #                                                            \
    # MLG / NLG
    
    # Check if LG is activated
    d_lg = 0
    if x_nlg is not None:
    
        # Make landing gear dimensions based on the fuselage
        w_lg = 0.05*D_f
        d_lg = 4*w_lg
        
        mlg_len = np.linspace(y_mlg-w_lg/2, y_mlg+w_lg/2, 2)
        nlg_len = np.linspace(y_nlg-w_lg/2, y_nlg+w_lg/2, 2)
        
        for i in range(len(mlg_len)):
            p = Ellipse((x_mlg, z_mlg), d_lg, d_lg, angle=0,\
            facecolor = 'gray', edgecolor = 'k', lw=2)
            ax.add_patch(p)
            art3d.pathpatch_2d_to_3d(p, z=mlg_len[i], zdir="y")
            
            p = Ellipse((x_mlg, z_mlg), d_lg, d_lg, angle=0,\
            facecolor = 'gray', edgecolor = 'k', lw=2)
            ax.add_patch(p)
            art3d.pathpatch_2d_to_3d(p, z=-mlg_len[i], zdir="y")
    
            # NLG
            p = Ellipse((x_nlg, z_nlg), d_lg, d_lg, angle=0,\
            facecolor = 'gray', edgecolor = 'k', lw=1.5)
            ax.add_patch(p)
            art3d.pathpatch_2d_to_3d(p, z=nlg_len[i], zdir="y")

    # Nacelle
    nc_len = np.linspace(x_n,x_n+L_n,11)
    for i in range(len(nc_len)):
        p = Ellipse((y_n, z_n), D_n, D_n, angle=0,\
        facecolor = 'none', edgecolor = 'orange', lw=1.0)
        ax.add_patch(p)
        art3d.pathpatch_2d_to_3d(p, z=nc_len[i], zdir="x")

        # Inner wall
        #p = Ellipse((y_n, z_n), D_n*0.8, D_n*0.8, angle=0,\
        #facecolor = 'none', edgecolor = 'k', lw=.1)
        #ax.add_patch(p)
        #art3d.pathpatch_2d_to_3d(p, z=nc_len[i], zdir="x")


        p = Ellipse((-y_n, z_n), D_n, D_n, angle=0, \
        facecolor = 'none', edgecolor = 'orange', lw=1.0)
        ax.add_patch(p)
        art3d.pathpatch_2d_to_3d(p, z=nc_len[i], zdir="x")

        # Inner wall
        #p = Ellipse((-y_n, z_n), D_n*0.8, D_n*0.8, angle=0, \
        #facecolor = 'none', edgecolor = 'k', lw=.1)
        #ax.add_patch(p)
        #art3d.pathpatch_2d_to_3d(p, z=nc_len[i], zdir="x")

    # Aileron
    ail_tip_margin = 0.02 # Margem entre flap e aileron em % de b_w

    # Spanwise positions (root and tip)
    yr_a = (1.0 - (ail_tip_margin + b_ail_b_wing))*b_w/2
    yt_a = (1.0 - (ail_tip_margin))*b_w/2

    cr_a = lin_interp(0, b_w/2, cr_w, ct_w, yr_a)*c_ail_c_wing
    ct_a = lin_interp(0, b_w/2, cr_w, ct_w, yt_a)*c_ail_c_wing

    # To find the longitudinal position of the aileron LE, we find the TE position first
    # then we subtract the chord
    xr_a = lin_interp(0, b_w/2, xr_w+cr_w, xt_w+ct_w, yr_a) - cr_a
    xt_a = lin_interp(0, b_w/2, xr_w+cr_w, xt_w+ct_w, yt_a) - ct_a

    zr_a = lin_interp(0, b_w/2, zr_w, zt_w, yr_a)
    zt_a = lin_interp(0, b_w/2, zr_w, zt_w, yt_a)

    # Airfoil thickness at aileron location
    tcr_a = lin_interp(0, b_w/2, tcr_w, tct_w, yr_a)
    tct_a = lin_interp(0, b_w/2, tcr_w, tct_w, yt_a)

    ax.plot([xr_a, xt_a, xt_a+ct_a, xr_a+cr_a, xr_a],
            [yr_a, yt_a, yt_a     , yr_a     , yr_a],
            [zr_a+cr_a*tcr_a/2/c_ail_c_wing, zt_a+ct_a*tct_a/2/c_ail_c_wing, zt_a+ct_a*tct_a/2/c_ail_c_wing     , zr_a+cr_a*tcr_a/2/c_ail_c_wing, zr_a+cr_a*tcr_a/2/c_ail_c_wing],lw=1,color='green')

    ax.plot([ xr_a,  xt_a,  xt_a+ct_a,  xr_a+cr_a,  xr_a],
            [-yr_a, -yt_a, -yt_a     , -yr_a     , -yr_a],
            [ zr_a+cr_a*tcr_a/2/c_ail_c_wing,  zt_a+ct_a*tct_a/2/c_ail_c_wing,  zt_a+ct_a*tct_a/2/c_ail_c_wing,  zr_a+cr_a*tcr_a/2/c_ail_c_wing     ,  zr_a+cr_a*tcr_a/2/c_ail_c_wing],lw=1,color='green')

    # Fuel tank

    # Spanwise positions (root and tip)
    yr_tk = b_tank_b_w_start*b_w/2
    yt_tk = b_tank_b_w_end*b_w/2

    cr_tk = lin_interp(0, b_w/2, cr_w, ct_w, yr_tk)*c_tank_c_w
    ct_tk = lin_interp(0, b_w/2, cr_w, ct_w, yt_tk)*c_tank_c_w

    # To find the longitudinal position of the tank LE
    xr_tk = lin_interp(0, b_w/2, xr_w, xt_w, yr_tk) + cr_tk*x_tank_c_w/c_tank_c_w
    xt_tk = lin_interp(0, b_w/2, xr_w, xt_w, yt_tk) + ct_tk*x_tank_c_w/c_tank_c_w

    zr_tk = lin_interp(0, b_w/2, zr_w, zt_w, yr_tk)
    zt_tk = lin_interp(0, b_w/2, zr_w, zt_w, yt_tk)

    # Airfoil thickness at tank location
    tcr_tk = lin_interp(0, b_w/2, tcr_w, tct_w, yr_tk)
    tct_tk = lin_interp(0, b_w/2, tcr_w, tct_w, yt_tk)

    ax.plot([xr_tk, xt_tk, xt_tk+ct_tk, xr_tk+cr_tk, xr_tk],
            [yr_tk, yt_tk, yt_tk     , yr_tk     , yr_tk],
            [zr_tk+cr_tk*tcr_tk/2, zt_tk+ct_tk*tct_tk/2, zt_tk+ct_tk*tct_tk/2     , zr_tk+cr_tk*tcr_tk/2, zr_tk+cr_tk*tcr_tk/2],lw=1,color='magenta')

    ax.plot([ xr_tk,  xt_tk,  xt_tk+ct_tk,  xr_tk+cr_tk,  xr_tk],
            [-yr_tk, -yt_tk, -yt_tk     , -yr_tk     , -yr_tk],
            [ zr_tk+cr_tk*tcr_tk/2,  zt_tk+ct_tk*tct_tk/2,  zt_tk+ct_tk*tct_tk/2,  zr_tk+cr_tk*tcr_tk/2     ,  zr_tk+cr_tk*tcr_tk/2],lw=1,color='magenta')

    # Slat
    if slat_type is not None:
        
        #slat_tip_margin = 0.02  # Margem da ponta como % da b_w
        #slat_root_margin = 0.12 # Margem da raiz como % da b_w
        #hist_c_s = 0.25        # Corda do Flap
        #hist_b_s = 1 - slat_root_margin - slat_tip_margin

        # Spanwise positions (root and tip)
        yr_s = D_f/2
        yt_s = b_slat_b_wing*b_w/2

        cr_s = lin_interp(0, b_w/2, cr_w, ct_w, yr_s)*c_slat_c_wing
        ct_s = lin_interp(0, b_w/2, cr_w, ct_w, yt_s)*c_slat_c_wing

        # Find the longitudinal position of the slat LE
        xr_s = lin_interp(0, b_w/2, xr_w, xt_w, yr_s)
        xt_s = lin_interp(0, b_w/2, xr_w, xt_w, yt_s)

        zr_s = lin_interp(0, b_w/2, zr_w, zt_w, yr_s)
        zt_s = lin_interp(0, b_w/2, zr_w, zt_w, yt_s)

        # Airfoil thickness at slat location
        tcr_s = lin_interp(0, b_w/2, tcr_w, tct_w, yr_s)
        tct_s = lin_interp(0, b_w/2, tcr_w, tct_w, yt_s)


        ax.plot([xr_s, xt_s, xt_s+ct_s, xr_s+cr_s, xr_s],
                [yr_s, yt_s, yt_s     , yr_s     , yr_s],
                [zr_s+cr_s*tcr_s/2/c_slat_c_wing, zt_s+ct_s*tct_s/2/c_slat_c_wing, zt_s+ct_s*tct_s/2/c_slat_c_wing     , zr_s+cr_s*tcr_s/2/c_slat_c_wing, zr_s+cr_s*tcr_s/2/c_slat_c_wing],lw=1,color='m')

        ax.plot([ xr_s,  xt_s,  xt_s+ct_s,  xr_s+cr_s,  xr_s],
                [-yr_s, -yt_s, -yt_s     , -yr_s     , -yr_s],
                [ zr_s+cr_s*tcr_s/2/c_slat_c_wing,  zt_s+ct_s*tct_s/2/c_slat_c_wing,  zt_s+ct_s*tct_s/2/c_slat_c_wing,  zr_s+cr_s*tcr_s/2/c_slat_c_wing     ,  zr_s+cr_s*tcr_s/2/c_slat_c_wing],lw=1,color='m')

    # Flap outboard
    if flap_type is not None:

        # Spanwise positions (root and tip)
        yr_f = D_f/2
        yt_f = b_flap_b_wing*b_w/2

        cr_f = lin_interp(0, b_w/2, cr_w, ct_w, yr_f)*c_flap_c_wing
        ct_f = lin_interp(0, b_w/2, cr_w, ct_w, yt_f)*c_flap_c_wing

        # To find the longitudinal position of the flap LE, we find the TE position first
        # then we subtract the chord
        xr_f = lin_interp(0, b_w/2, xr_w+cr_w, xt_w+ct_w, yr_f) - cr_f
        xt_f = lin_interp(0, b_w/2, xr_w+cr_w, xt_w+ct_w, yt_f) - ct_f

        zr_f = lin_interp(0, b_w/2, zr_w, zt_w, yr_f)
        zt_f = lin_interp(0, b_w/2, zr_w, zt_w, yt_f)

        # Airfoil thickness at flap location
        tcr_f = lin_interp(0, b_w/2, tcr_w, tct_w, yr_f)
        tct_f = lin_interp(0, b_w/2, tcr_w, tct_w, yt_f)


        ax.plot([xr_f, xt_f, xt_f+ct_f, xr_f+cr_f, xr_f],
                [yr_f, yt_f, yt_f     , yr_f     , yr_f],
                [zr_f+cr_f*tcr_f/2/c_flap_c_wing, zt_f+ct_f*tct_f/2/c_flap_c_wing, zt_f+ct_f*tct_f/2/c_flap_c_wing     , zr_f+cr_f*tcr_f/2/c_flap_c_wing, zr_f+cr_f*tcr_f/2/c_flap_c_wing],lw=1,color='r')

        ax.plot([ xr_f,  xt_f,  xt_f+ct_f,  xr_f+cr_f,  xr_f],
                [-yr_f, -yt_f, -yt_f     , -yr_f     , -yr_f],
                [ zr_f+cr_f*tcr_f/2/c_flap_c_wing,  zt_f+ct_f*tct_f/2/c_flap_c_wing,  zt_f+ct_f*tct_f/2/c_flap_c_wing,  zr_f+cr_f*tcr_f/2/c_flap_c_wing     ,  zr_f+cr_f*tcr_f/2/c_flap_c_wing],lw=1,color='r')

    # Elevator
    ele_tip_margin = 0.1  # Margem do profundor para a ponta
    ele_root_margin = 0.1 # Margem do profundor para a raiz
    hist_b_e = 1-ele_root_margin-ele_tip_margin
    hist_c_e = 0.25


    ct_e_loc = (1-ele_tip_margin)*(ct_h - cr_h)+cr_h
    cr_e_loc = (1-hist_b_e-ele_tip_margin)*(ct_h - cr_h)+cr_h

    ct_e = ct_e_loc*hist_c_e
    cr_e = cr_e_loc*hist_c_e

    xr_e = (1-hist_b_e-ele_tip_margin)*(xt_h - xr_h)+xr_h + cr_e_loc*(1-hist_c_e)
    xt_e = (1-ele_tip_margin)*(xt_h - xr_h)+xr_h + ct_e_loc*(1-hist_c_e)

    yr_e = (1-hist_b_e-ele_tip_margin)*b_h/2
    yt_e = (1-ele_tip_margin)*b_h/2

    zr_e = (1-hist_b_e-ele_tip_margin)*(zt_h - zr_h)+zr_h
    zt_e = (1-ele_tip_margin)*(zt_h - zr_h)+zr_h



    ax.plot([xr_e, xt_e, xt_e+ct_e, xr_e+cr_e, xr_e],
            [yr_e, yt_e, yt_e     , yr_e     , yr_e],
            [zr_e, zt_e, zt_e     , zr_e     , zr_e],lw=1,color='g')

    ax.plot([ xr_e,  xt_e,  xt_e+ct_e,  xr_e+cr_e,  xr_e],
            [-yr_e, -yt_e, -yt_e     , -yr_e     , -yr_e],
            [ zr_e,  zt_e,  zt_e     ,  zr_e     ,  zr_e],lw=1,color='g')

    # Rudder
    ver_base_margin = 0.1               # Local da base % de b_v
    ver_tip_margin1 = 0.1               # Local da base % de b_v
    ver_tip_margin = 1-ver_tip_margin1  # Local do topo % de b_v
    hist_c_v = 0.32

    cr_v_loc = ver_base_margin*(ct_v - cr_v)+cr_v
    ct_v_loc = ver_tip_margin*(ct_v - cr_v)+cr_v


    cr_v2 = cr_v_loc*hist_c_v
    ct_v2 = ct_v_loc*hist_c_v


    xr_v2 = ver_base_margin*(xt_v - xr_v)+xr_v+cr_v_loc*(1-hist_c_v)
    xt_v2 = ver_tip_margin*(xt_v - xr_v)+xr_v+ct_v_loc*(1-hist_c_v)


    zr_v2 = ver_base_margin*(zt_v - zr_v)+zr_v
    zt_v2 = ver_tip_margin*(zt_v - zr_v)+zr_v



    ax.plot([xr_v2  , xt_v2  , xt_v2+ct_v2   , xr_v2+cr_v2   , xr_v2        ],
            [tcr_v*cr_v_loc/2, tct_v*ct_v_loc/2, tct_v*ct_v_loc/2, \
            tcr_v*cr_v_loc/2, tcr_v*cr_v_loc/2],
            [zr_v2  , zt_v2   , zt_v2       , zr_v2        , zr_v2        ],\
            color='orange')


    ax.plot([xr_v2  , xt_v2  , xt_v2+ct_v2   , xr_v2+cr_v2   , xr_v2        ],
            [-tcr_v*cr_v_loc/2, -tct_v*ct_v_loc/2, -tct_v*ct_v_loc/2, \
            -tcr_v*cr_v_loc/2, -tcr_v*cr_v_loc/2],
            [zr_v2  , zt_v2   , zt_v2       , zr_v2        , zr_v2        ],\
            color='orange')

    # _______ONLY FRONT VIEW_______

    # Wing Lower
    #------------------------------
    ax.plot([xr_w    , xt_w, xt_w+ct_w, xr_w+cr_w, xt_w+ct_w, xt_w, xr_w],
            [0.0     , yt_w, yt_w, 0.0, -yt_w, -yt_w, 0.0],
            [zr_w-tcr_w*cr_w/2, zt_w-tct_w*ct_w/2, zt_w-tct_w*ct_w/2, zr_w-tcr_w*cr_w/2, \
             zt_w-tct_w*ct_w/2, zt_w-tct_w*ct_w/2, zr_w-tcr_w*cr_w/2],color='blue')

    ax.plot([xr_w         , xr_w],
            [0.0          , 0.0 ],
            [zr_w-tcr_w*cr_w/2, zr_w+tcr_w*cr_w/2],color='blue')
    ax.plot([xr_w+cr_w         , xr_w+cr_w],
            [0.0          , 0.0 ],
            [zr_w-tcr_w*cr_w/2, zr_w+tcr_w*cr_w/2],color='blue')

    ax.plot([xt_w         , xt_w],
            [yt_w         , yt_w ],
            [zt_w-tct_w*ct_w/2, zt_w+tct_w*ct_w/2],color='blue')
    ax.plot([xt_w+ct_w    , xt_w+ct_w],
            [yt_w         , yt_w ],
            [zt_w-tct_w*ct_w/2, zt_w+tct_w*ct_w/2],color='blue')

    ax.plot([xt_w         , xt_w],
            [-yt_w         , -yt_w ],
            [zt_w-tct_w*ct_w/2, zt_w+tct_w*ct_w/2],color='blue')
    ax.plot([xt_w+ct_w    , xt_w+ct_w],
            [-yt_w         , -yt_w ],
            [zt_w-tct_w*ct_w/2, zt_w+tct_w*ct_w/2],color='blue')

    #------------------------------



    # HT Lower
    #------------------------------
    ax.plot([xr_h    , xt_h, xt_h+ct_h, xr_h+cr_h, xt_h+ct_h, xt_h, xr_h],
            [0.0     , yt_h, yt_h, 0.0, -yt_h, -yt_h, 0.0],
            [zr_h-tcr_h*cr_h/2, zt_h-tct_h*ct_h/2, zt_h-tct_h*ct_h/2, zr_h-tcr_h*cr_h/2, \
             zt_h-tct_h*ct_h/2, zt_h-tct_h*ct_h/2, zr_h-tcr_h*cr_h/2],color='green')

    ax.plot([xr_h         , xr_h],
            [0.0          , 0.0 ],
            [zr_h-tcr_h*cr_h/2, zr_h+tcr_h*cr_h/2],color='green')
    ax.plot([xr_h+cr_h         , xr_h+cr_h],
            [0.0          , 0.0 ],
            [zr_h-tcr_h*cr_h/2, zr_h+tcr_h*cr_h/2],color='green')

    ax.plot([xt_h         , xt_h],
            [yt_h         , yt_h ],
            [zt_h-tct_h*ct_h/2, zt_h+tct_h*ct_h/2],color='green')
    ax.plot([xt_h+ct_h    , xt_h+ct_h],
            [yt_h         , yt_h ],
            [zt_h-tct_h*ct_h/2, zt_h+tct_h*ct_h/2],color='green')

    ax.plot([ xt_h         ,  xt_h],
            [-yt_h         , -yt_h ],
            [ zt_h-tct_h*ct_h/2, zt_h+tct_h*ct_h/2],color='green')
    ax.plot([ xt_h+ct_h    ,  xt_h+ct_h],
            [-yt_h         , -yt_h ],
            [ zt_h-tct_h*ct_h/2, zt_h+tct_h*ct_h/2],color='green')


    # Slat Lower
    #------------------------------
    if slat_type is not None:
        ax.plot([xr_s, xt_s, xt_s+ct_s, xr_s+cr_s, xr_s],
                [yr_s, yt_s, yt_s     , yr_s     , yr_s],
                [zr_s-tcr_s*cr_s/2/c_slat_c_wing ,\
                 zt_s-tct_s*ct_s/2/c_slat_c_wing ,\
                 zt_s-tct_s*ct_s/2/c_slat_c_wing ,\
                 zr_s-tcr_s*cr_s/2/c_slat_c_wing ,\
                 zr_s-tcr_s*cr_s/2/c_slat_c_wing],\
                 lw=1,color='m')

        ax.plot([ xr_s,  xt_s,  xt_s+ct_s,  xr_s+cr_s,  xr_s],
                [-yr_s, -yt_s, -yt_s     , -yr_s     , -yr_s],
                [ zr_s-tcr_s*cr_s/2/c_slat_c_wing,\
                  zt_s-tct_s*ct_s/2/c_slat_c_wing,\
                  zt_s-tct_s*ct_s/2/c_slat_c_wing,\
                  zr_s-tcr_s*cr_s/2/c_slat_c_wing,\
                  zr_s-tcr_s*cr_s/2/c_slat_c_wing],\
                  lw=1,color='m')
    #------------------------------



    # Flap Lower
    #------------------------------
    if flap_type is not None:
        ax.plot([xr_f, xt_f, xt_f+ct_f, xr_f+cr_f, xr_f],
                [yr_f, yt_f, yt_f     , yr_f     , yr_f],
                [zr_f-tcr_f*cr_f/2/c_flap_c_wing ,\
                 zt_f-tct_f*ct_f/2/c_flap_c_wing ,\
                 zt_f-tct_f*ct_f/2/c_flap_c_wing ,\
                 zr_f-tcr_f*cr_f/2/c_flap_c_wing ,\
                 zr_f-tcr_f*cr_f/2/c_flap_c_wing],\
                 lw=1,color='r')

        ax.plot([ xr_f,  xt_f,  xt_f+ct_f, xr_f+cr_f, xr_f],
                [-yr_f, -yt_f, -yt_f     ,-yr_f     ,-yr_f],
                [zr_f-tcr_f*cr_f/2/c_flap_c_wing ,\
                 zt_f-tct_f*ct_f/2/c_flap_c_wing ,\
                 zt_f-tct_f*ct_f/2/c_flap_c_wing ,\
                 zr_f-tcr_f*cr_f/2/c_flap_c_wing ,\
                 zr_f-tcr_f*cr_f/2/c_flap_c_wing],\
                 lw=1,color='r')
    #------------------------------



    # Aleron Lower
    #------------------------------
    ax.plot([xr_a, xt_a, xt_a+ct_a, xr_a+cr_a, xr_a],
            [yr_a, yt_a, yt_a     , yr_a     , yr_a],
            [zr_a-tcr_a*cr_a/2/c_ail_c_wing ,\
             zt_a-tct_a*ct_a/2/c_ail_c_wing ,\
             zt_a-tct_a*ct_a/2/c_ail_c_wing ,\
             zr_a-tcr_a*cr_a/2/c_ail_c_wing ,\
             zr_a-tcr_a*cr_a/2/c_ail_c_wing],\
             lw=1,color='green')

    ax.plot([ xr_a,  xt_a,  xt_a+ct_a, xr_a+cr_a, xr_a],
            [-yr_a, -yt_a, -yt_a     ,-yr_a     ,-yr_a],
            [zr_a-tcr_a*cr_a/2/c_ail_c_wing ,\
             zt_a-tct_a*ct_a/2/c_ail_c_wing ,\
             zt_a-tct_a*ct_a/2/c_ail_c_wing ,\
             zr_a-tcr_a*cr_a/2/c_ail_c_wing ,\
             zr_a-tcr_a*cr_a/2/c_ail_c_wing],\
             lw=1,color='green')
    #------------------------------

    # Avoiding blanketing the rudder
    ax.plot([xr_h         , xr_h+b_v/np.tan(60*np.pi/180)],
            [0.0          , 0.0 ],
            [zr_h, zr_h+b_v],'k--')


    ax.plot([xr_h+cr_h         , xr_h+0.6*b_v/np.tan(30*np.pi/180)+cr_h],
            [0.0          , 0.0 ],
            [zr_h, zr_h+0.6*b_v],'k--')

    # Auxiliary landing gear lines
    if x_nlg is not None:

        # Water Spray
        ax.plot([x_nlg         , x_nlg+0.25*b_w/np.tan(22*np.pi/180)],
                [0.0          , 0.25*b_w ],
                [z_nlg, z_nlg],'k--')
    
        ax.plot([x_nlg         , x_nlg+0.25*b_w/np.tan(22*np.pi/180)],
                [0.0          , -0.25*b_w ],
                [z_nlg, z_nlg],'k--')

        # Tailstrike
        #tailstrike_angle = np.arctan((-D_f/2-z_mlg)/(x_tailstrike-x_mlg)) # This one uses the fuselage diameter as reference
        tailstrike_angle = np.arctan((z_tailstrike-z_mlg)/(x_tailstrike-x_mlg)) # This one uses the fuselage diameter as reference
        ax.plot([x_mlg         , L_f],
                [0.0          , 0.0 ],
                [z_mlg, z_mlg+(L_f-x_mlg)*np.tan(tailstrike_angle)],'k--')
    
        ax.plot([x_mlg         , L_f],
                [0.0          , 0.0 ],
                [z_mlg, z_mlg],'k--')


    # Equal aspect ratio: set all axes to the same range centered on the data
    X = np.array([0, xr_w, xt_h+ct_h, xt_v+ct_v, L_f, xr_h+b_v/np.tan(60*np.pi/180), xr_h+0.6*b_v/np.tan(30*np.pi/180)+cr_h])
    Y = np.array([-yt_w, yt_w])
    Z = np.array([-D_f/2, zt_w, zt_h, zt_v, z_mlg-d_lg/2, zr_h+b_v])
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2
    mid_x = 0.5*(X.max()+X.min())
    mid_y = 0.5*(Y.max()+Y.min())
    mid_z = 0.5*(Z.max()+Z.min())
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(az1, az2)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.zaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(False)
    fig.tight_layout()
    fig.subplots_adjust(left=-0.15, right=1.15, top=1.15, bottom=-0.15)
    
   
#----------------------------------------