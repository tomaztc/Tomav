# GENERAL IMPORTS
import numpy as np
from .constants import gamma_air, R_air
from .auxiliary import atmosphere

#======================================

def engineTSFC(Mach, altitude, airplane):
    '''
    This function computes the engine thrust-specific fuel
    consumption and thrust correction factor compared to
    static sea-level conditions. The user has to define the
    engine parameters in a 'engine' dictionary within
    the airplane dictionary. The engine model must be
    identified by the 'model' field of the engine dictionary.
    The following engine models are available:
    
    Howe TSFC turbofan model:
    requires the bypass ratio. An optional sea-level TSFC
    could also be provided. Otherwise, standard parameters
    are used.
    airplane['engine'] = {'model': 'howe turbofan',
                          'BPR': 3.04,
                          'Cbase': 0.58/3600} # Could also be None
                          
    Thermodynamic cycle turbojet:
    This model uses a simplified thermodynamic model of
    turbofans to estimate maximum thrust and TSFC
    
    airplane['engine'] = {'model': 'thermo turbojet'
                          'data': dictionary (check turbojet_model function)}
    
    The user can also leave a 'weight' field in the dictionary
    to replace the weight estimation.
    '''

    # Get a reference to the engine dictionary
    engine = airplane['inputs']['engine']

    # @REMOVE

    # Check which model was given
    if engine['model'].lower() == 'howe turbofan':

        # Unpack dictionary
        BPR = engine['BPR']
    
        # Determine base TSFC of Howe's regression
        if 'C_ref' not in engine.keys(): # User did not provide a base TSFC

            if BPR < 4.0:
                C_base = 0.85/3600
            else:
                C_base = 0.70/3600
        
        else:

            # Get reference point information from the dictionary
            C_ref = engine['C_ref']
            altitude_ref = engine['altitude_ref']
            Mach_ref = engine['Mach_ref']

            # Atmospheric conditions at cruise altitude
            atm_data = atmosphere(altitude_ref)
            rho_ref = atm_data['density']
    
            # Density ratio
            sigma_ref = rho_ref/1.225
    
            # Correct Cbase so that the equation gives the desired static TSFC at sea-level
            # Howe Eq 3.12a
            C_base = C_ref/(1-0.15*BPR**0.65)/(1+0.28*(1+0.063*BPR**2)*Mach_ref)/sigma_ref**0.08


        # Atmospheric conditions at cruise altitude
        atm_data = atmosphere(altitude)
        rho = atm_data['density']
    
        # Density ratio
        sigma = rho/1.225
        
        # Howe Eq 3.12a
        C = C_base*(1-0.15*BPR**0.65)*(1+0.28*(1+0.063*BPR**2)*Mach)*sigma**0.08

        # Cruise traction correction for takeoff conditions by Scholz
        # (https://www.fzt.haw-hamburg.de/pers/Scholz/HOOU/AircraftDesign_5_PreliminarySizing.pdf)
        kT = (0.0013*BPR-0.0397)*altitude/1000.0 - 0.0248*BPR + 0.7125
        
        '''
        # Cruise traction correction for takeoff conditions I created
     tsfc   # by fitting curves to Fig. 5.1 from Raymer, 5ed
        # Factor to migrate from low bypass to high bypass curves
        ff = (BPR-0.4)/7.6
        
        # Curve fit coefficients
        aa = -0.0012 + 0.0002*ff
        bb = -0.0317 + 0.0271*ff
        cc =  0.6936 - 0.3315*ff
        
        # Compute thrust correction
        kT = aa*(altitude/1000)**2 + bb*(altitude/1000) + cc

        # Add another factor due to Mach since Raymer's corretion is
        # for Mach 0.8. We developed this correction assuming.
        # KM is the ratio between T@Mach/T@static at the given altitude.
        # We assume that the correction is not necessary at high altitude since
        # The thrust curves are flat with respect to Mach at high altitudes
        Mref = 0.8
        href = 40000*ft2m
        KMref = cc + (1-cc)*hh/href
        kT = kT/KMref*(1 + (KMref-1)*Mach/Mref)
        '''
    
    elif engine['model'].lower() == 'raymer turbofan':

        # Check if the interpolants where already loaded
        if not 'thrust_interp' in engine.keys():
            load_raymer_turbofan(airplane)

        # Get sea-level static TSFC
        C_base = engine['Cbase']

        # Use the interpolant
        kT = engine['thrust_interp'](Mach, altitude)
        C = C_base*engine['tsfc_interp'](Mach, altitude)

        print('')
        print('Mach:',Mach)
        print('altitude:',altitude)
        print('TSFC:',C*3600)
        print('kT:',kT)

    elif engine['model'].lower() == 'thermo turbojet':
        
        C, F = turbojet_model(Mach, altitude, engine['data'])
        
        # Check if maximum sea-level thrust was already computed
        if 'T0_eng' not in engine:
            
            _, F0 = turbojet_model(0.0, 0.0, engine['data'])
            airplane['engine']['T0_eng'] = F0
            
        # Compute thrust correction factor
        kT = F/engine['T0_eng']

    # @REMOVE

    return C, kT

#----------------------------------------

def load_raymer_turbofan(airplane):
    '''
    This function loads max thrust and TSFC data from the high bypass turbofan
    from Raymer, 5th Ed, Appendix E.
    This information isstored in a pickle file generated after digitalizing the plot.
    This function stores interpolation functions in the airplane dictionary:
    T = thrust_interp(Mach, altitude)
    TSFC = tsfc_interp(Mach, altitude)
    '''

    import pickle

    # Load pickle data
    with open(f'max_thrust_data.pickle','rb') as fid:
        [Mach,altitude,thrust] = pickle.load(fid)

    # Normalize thrust by sea-level conditions
    thrust = thrust/thrust[0]

    # Create interpolator
    from scipy.interpolate import CloughTocher2DInterpolator as cubic
    thrust_interp = cubic(list(zip(Mach, altitude)), thrust)

    # Load pickle data
    with open(f'tsfc_data.pickle','rb') as fid:
        [Mach,altitude,tsfc] = pickle.load(fid)

    # Normalize tsfc by sea-level conditions
    tsfc = tsfc/tsfc[0]

    # Create interpolator
    from scipy.interpolate import CloughTocher2DInterpolator as cubic
    tsfc_interp = cubic(list(zip(Mach, altitude)), tsfc)

    # Store interpolants
    airplane['inputs']['engine']['thrust_interp'] = thrust_interp
    airplane['inputs']['engine']['tsfc_interp'] = tsfc_interp

#----------------------------------------
# @REMOVE
def turbojet_model(Mach, altitude, data):
    
    '''
    Model implemented by Marcelo Yuri Sampaio de Freitas
    during his undergraduate thesis:
    "Projeto conceitual de alvo aéreo manobrável baseado na propulsão de turbojato"
    Instituto Tecnológico de Aeronáutica, 2020.
    http://www.bdita.bibl.ita.br/TGsDigitais/lista_resumo.php?num_tg=77532
    
    This function returns the TSFC and maximum thrust of the engine

    data should be a dictionary with the following fields:
    data['p02_p01']: compression ratio
    data['T03']: turbine inlet temperature
    data['n_comp']: compressor efficiency
    data['n_turbine']: turbine efficiency
    data['n_intake']: intake efficiency
    data['n_isen_noz']: insentropic nozzle efficiency
    data['n_mec']: mechanical efficiency
    data['n_comb']: combustion efficiency
    data['delta_pb']: pressure loss at the combustion chamber
    data['A5']: nozzle area
    '''
        
    # Thermodynamic parameters
    Cp_air = 1005 # J/kg*K
    gamma_comb = 1.333
    Cp_comb = 1148 # J/kg*K
    LHV = 43.1e6 #J/kg fuel calorific power
        
    # Air properties
    atm_data = atmosphere(altitude)
    sound_speed = atm_data['speed_of_sound']
    Ta = atm_data['temperature']
    pa = atm_data['pressure']
        
    # Generic engine parameters
    p02_p01 = data['p02_p01']
    T03 = data['T03']
    n_comp = data['n_comp']
    n_turbine = data['n_turbine']
    n_intake = data['n_intake']
    n_isen_noz = data['n_isen_noz']
    n_mec = data['n_mec']
    n_comb = data['n_comb']
    delta_pb = data['delta_pb'] # Pressure loss in combustion chamber

    #O requisito pode ser alterado para fluxo de massa
    A5 = data['A5']
    #m_dot = 7.823
        
    # Intake
    T_estag = ((Mach*sound_speed)**2)/(2*Cp_air)
    T01 = Ta + T_estag
    p01_pa = (1 + n_intake*T_estag/Ta)**(gamma_air/(gamma_air-1))
    p01 = p01_pa*pa
        
    # Compressor
    p02 = p02_p01*p01
    T02 = T01 + (T01/n_comp)*(p02_p01**((gamma_air-1)/gamma_air) - 1)
        
    # Combustion
    p03 = p02*(1-(delta_pb/100))
    T03 = T03
        
    # Turbine
    deltaT = Cp_air*(T02 - T01)/Cp_comb/n_mec
    T04 = T03 - deltaT
    T04_isen = T03 - deltaT/n_turbine
    p04 = p03*(T04_isen/T03)**(gamma_comb/(gamma_comb-1))
    
    # Exhaust nozzle
    p04_pa = p04/pa
    p04_pc = 1/((1 - (1/n_isen_noz)*((gamma_comb - 1)/(gamma_comb + 1)))**(gamma_comb/(gamma_comb - 1)))
        
    if p04_pa >= p04_pc:
            
        # Choked exhaust
        T05 = (2/(gamma_comb + 1))*T04
        p05 = p04/p04_pc
        rho5 = p05/(R_air*T05)
        C5 = np.sqrt(gamma_comb*R_air*T05)
        A5s = 1/(rho5*C5)
        Fs = (C5 - Mach*sound_speed) + A5s*(p05 - pa)
        m_dot = A5/A5s
        F = Fs*m_dot # Thrust
        f_comb = (Cp_comb*T03 - Cp_air*T02)/(LHV*n_comb - Cp_comb*T03)
        m_dot_fuel = f_comb*m_dot
        SFC = m_dot_fuel/F #kg/s.N
            
    else:
            
        # Exhaust is not choked
        T05 = T04 - n_isen_noz*T04*(1 - (1/p04_pa)**((gamma_comb-1)/gamma_comb))
        C5 = np.sqrt(2*Cp_comb*(T04 - T05))
        p05 = pa
        rho5 = p05/(R_air*T05)
        A5s = 1/(rho5*C5)
        Fs = (C5 - Mach*sound_speed) + A5s*(p05 - pa)
        m_dot = A5/A5s
        F = Fs*m_dot # Thrust
        f_comb = (Cp_comb*T03 - Cp_air*T02)/(LHV*n_comb - Cp_comb*T03)
        m_dot_fuel = f_comb*m_dot
        SFC = m_dot_fuel/F #kg/s.N
        #sigma = rho/1.225
        
    # Convert to 1/s
    C = SFC*9.81
        
    return C, F