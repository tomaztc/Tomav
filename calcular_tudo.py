
# coding: utf-8
# %% IMPORTS
from importlib.util import find_spec
if find_spec("matplotlib") is None or find_spec("scipy") is None or find_spec("numpy") is None:
    import sys
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "matplotlib", "scipy", "numpy"])

import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import minimize, NonlinearConstraint, minimize_scalar, brentq
from scipy import integrate
import numpy as np

import os
import csv
import json
import warnings

from designTool.aerodynamics import aerodynamics
from designTool.auxiliary import atmosphere
from designTool.balance import balance
from designTool.constants import gravity, ft2m
from designTool.geometry import geometry
from designTool.performance import thrust_matching
from designTool.plots import plot_geometry
from designTool.standard_airplane import standard_airplane
from designTool.weight import weight
from designTool.landing_gear import landing_gear
from designTool.propulsion import engineTSFC
from designTool.analyze import analyze

# %% CONFIGURAÇÕES DE PLOT
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['legend.fontsize'] = 7
mpl.rcParams['figure.max_open_warning'] = 0
plt.rcParams['axes.grid'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0

# %% FUNÇÕES
def flatten_dict(d, parent_key=""):
    flattened = {}
    def flatten(current_dict, current_parent_key):
        for k, v in current_dict.items():
            new_key = f"{current_parent_key}_{k}" if current_parent_key else str(k)
            if isinstance(v, dict) and v:
                flatten(v, new_key)
            else:
                if new_key in flattened:
                    raise KeyError
                flattened[new_key] = v

    flatten(d, parent_key)
    return flattened

def aerodynamics_condicao_e_CL(airplane, condicao, CL):
    CD, CLmax, dragDict = aerodynamics(airplane, condicao["Mach"], condicao["altitude"], CL,
                                       n_engines_failed=condicao["n_engines_failed"],
                                       highlift_config=condicao["highlift_config"], lg_down=condicao["lg_down"],
                                       h_ground=condicao["h_ground"])
    return CD, CLmax, dragDict

def calcular_CL_voo(condicao, peso):
    atm = atmosphere(condicao["altitude"])
    rho = atm["density"]
    Mach = condicao["Mach"]
    V = Mach*atm["speed_of_sound"]
    CL = 2*peso / (rho*V**2*S_w)
    return CL

def plotar_curva_polar_condicao(condicao, titulo, LD=False):
    CL = condicao["CL"]
    CD, CLmax, _ = aerodynamics_condicao_e_CL(airplane, condicao, CL)

    if LD == True:
        plt.scatter(CL, CL/CD, s=30, zorder=5, label=f"Voo no {titulo}: CL={CL:.2f}, L/D={(CL / CD):.2f}")
    else:
        plt.scatter(CD, CL, s=30, zorder=5, label=f"Voo no {titulo}: CL={CL:.2f}, CD={CD:.4f}")

    CL_space = np.linspace(0, CLmax, 200)
    CD_space = []
    for CL in CL_space:
        CD, _, _ = aerodynamics_condicao_e_CL(airplane, condicao, CL)
        CD_space.append(CD)

    for CL in (0, CLmax):
        CD, _, _ = aerodynamics_condicao_e_CL(airplane, condicao, CL)
        if LD == True:
            plt.scatter(CL, CL / CD, color="black", s=20, zorder=5, label="_")
            plt.annotate(f"CL={CL:.2f}\nL/D={(CL/CD):.2f}", (CL, CL/CD), fontsize=7)
        else:
            plt.scatter(CD, CL, color="black", s=20, zorder=5, label="_")
            plt.annotate(f"CL={CL:.2f}\nCD={CD:.4f}", (CD, CL), fontsize=7)

    if LD == True:
        LD_space = CL_space / CD_space
        LD_max, CL_max_LD = max(zip(LD_space, CL_space))
        plt.scatter(CL_max_LD, LD_max, color="black", s=20, zorder=5, label="_")
        plt.annotate(f"Máximo:\nL/D={LD_max:.2f}\nCL={CL_max_LD:.2f}", (CL_max_LD, LD_max), ha="center", fontsize=7)
        plt.plot(CL_space, LD_space, label=f"Polar do {titulo}")
    else:
        plt.plot(CD_space, CL_space, label=f"Polar do {titulo}")

def calcular_CD_Mach_CL_variavel(Mach):
    condicao_Mach = CondicoesPolar["Cruzeiro"].copy()
    condicao_Mach["Mach"] = Mach
    CL_Mach = calcular_CL_voo(condicao_Mach, PesosIniciais["cruise"])
    CD_Mach, _, _ = aerodynamics_condicao_e_CL(airplane, condicao_Mach, CL_Mach)
    return CD_Mach

def calcular_CD_Mach_CL_fixo(Mach, CL):
    condicao_Mach = CondicoesPolar["Cruzeiro"].copy()
    condicao_Mach["Mach"] = Mach
    CD_Mach, _, _ = aerodynamics_condicao_e_CL(airplane, condicao_Mach, CL)
    return CD_Mach

def calcular_ld_max(airplane, condicao):
    _, CLmax, _ = aerodynamics_condicao_e_CL(airplane, condicao, 0)
    CL_space = np.linspace(0.01, CLmax, 500)
    LD_space = []
    for CL in CL_space:
        CD, _, _ = aerodynamics_condicao_e_CL(airplane, condicao, CL)
        LD_space.append(CL / CD)

    indice_ld_max = int(np.argmax(LD_space))
    return LD_space[indice_ld_max], CL_space[indice_ld_max], CLmax

def calcular_ld_operacional(airplane, condicao, peso):
    atm_data = atmosphere(condicao["altitude"])
    velocidade = condicao["Mach"] * atm_data["speed_of_sound"]
    CL = 2 * peso / (atm_data["density"] * velocidade ** 2 * airplane["inputs"]["S_w"])
    CD, CLmax, _ = aerodynamics_condicao_e_CL(airplane, condicao, CL)
    return CL / CD, CL, CLmax

def velocidade_em_ld_max(airplane, condicao, peso, CL_ld_max):
    atm_data = atmosphere(condicao["altitude"])
    return np.sqrt(2 * peso / (atm_data["density"] * airplane["inputs"]["S_w"] * CL_ld_max))

def calcular_dados_otimizacao(airplane_parametro):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        try:
            geometry(airplane_parametro)
            thrust_matching(airplane_parametro["inputs"]["W0_guess"], airplane_parametro["inputs"]["n_engines"]*airplane_parametro["inputs"]["engine"]["Tmax"], airplane_parametro)
            balance(airplane_parametro)
            landing_gear(airplane_parametro)

            W0_parametro = airplane_parametro["thrust_matching"]["W0"]

            sigma_parametro =  atmosphere(airplane_parametro["inputs"]["altitude_takeoff"], airplane_parametro["inputs"]["deltaISA_takeoff"])["density"] / 1.225
            _, CLmax_TO_parametro, _ = aerodynamics_condicao_e_CL(airplane_parametro, CondicoesPolar["2º Segmento"], 0)
            d_TO_parametro = 0.2387 / (sigma_parametro * CLmax_TO_parametro * airplane_parametro["inputs"]["S_w"]) * W0_parametro**2 / T0

            _, CLmax_LD_parametro, _ = aerodynamics_condicao_e_CL(airplane_parametro, CondicoesPolar["Pouso"], 0)
            W_LD_parametro = W0_parametro * airplane_parametro["inputs"]["MLW_frac"]
            d_LD_parametro = ((W_LD_parametro / airplane_parametro["inputs"]["S_w"]) / (rho_LD * CLmax_LD_parametro) - B_LD) / A_LD

            return {
                "W0": W0_parametro,
                "W_empty": airplane_parametro["thrust_matching"]["W_empty"],
                "W_fuel": airplane_parametro["thrust_matching"]["W_fuel"],
                "d_TO": d_TO_parametro,
                "d_LD": d_LD_parametro,
                "SM_fwd": airplane_parametro["balance"]["SM_fwd"],
                "SM_aft": airplane_parametro["balance"]["SM_aft"],
                "CLv": airplane_parametro["balance"]["CLv"],
                "tank_excess": airplane_parametro["balance"]["tank_excess"],
                "mlg_xcg_margin": airplane_parametro["inputs"]["x_mlg"] - airplane_parametro["balance"]["xcg_aft"],
                "alpha_tailstrike": np.degrees(airplane_parametro["landing_gear"]["alpha_tailstrike"]),
                "alpha_tipback": np.degrees(airplane_parametro["landing_gear"]["alpha_tipback"]),
                "phi_overturn": np.degrees(airplane_parametro["landing_gear"]["phi_overturn"]),
                "frac_nlg_aft": airplane_parametro["landing_gear"]["frac_nlg_aft"],
                "frac_nlg_fwd": airplane_parametro["landing_gear"]["frac_nlg_fwd"],
            }
        except Exception:
            return None

def sombrear_mascara(ax, xs, mascara, color, alpha, label):
    dx = (xs[1] - xs[0]) / 2
    em_bloco = False
    inicio = None
    label_usado = False
    for i, x in enumerate(xs):
        if mascara[i] and not em_bloco:
            em_bloco = True
            inicio = x - dx if i > 0 else x
        elif not mascara[i] and em_bloco:
            ax.axvspan(inicio, x - dx, color=color, alpha=alpha, label=(label if not label_usado else "_"))
            label_usado = True
            em_bloco = False
    if em_bloco:
        ax.axvspan(inicio, xs[-1], color=color, alpha=alpha, label=(label if not label_usado else "_"))
        
def construir_airplane_parametros(vetor_parametros):
    kwargs = {key: vetor_parametros[i] for i, (_, _, key, _, _) in enumerate(PARAMETROS_VARREDURA) if key in KWARGS_PADRAO}
    airplane_parametros = standard_airplane(airplane_name, **kwargs)
    for i, (_, _, key, _, _) in enumerate(PARAMETROS_VARREDURA):
        if key not in KWARGS_PADRAO:
            airplane_parametros["inputs"][key] = vetor_parametros[i]
    return airplane_parametros

def normalizar_vetor_parametros(vetor_parametros):
    return (np.array(vetor_parametros) - LIMITES_INFERIORES_PARAMETROS) / (LIMITES_SUPERIORES_PARAMETROS - LIMITES_INFERIORES_PARAMETROS)

def desnormalizar_vetor_parametros(vetor_parametros_normalizados):
    vetor_parametros = LIMITES_INFERIORES_PARAMETROS + np.array(vetor_parametros_normalizados) * (LIMITES_SUPERIORES_PARAMETROS - LIMITES_INFERIORES_PARAMETROS)
    return np.clip(vetor_parametros, LIMITES_INFERIORES_PARAMETROS, LIMITES_SUPERIORES_PARAMETROS)

def dados_otimizacao_cache(vetor_parametros_normalizados):
    key = tuple(np.round(vetor_parametros_normalizados, 11))
    if key not in cache_otimizacao:
        vetor_parametros = desnormalizar_vetor_parametros(vetor_parametros_normalizados)
        cache_otimizacao[key] = calcular_dados_otimizacao(construir_airplane_parametros(vetor_parametros))
    return cache_otimizacao[key]

def funcao_objetivo_otimizacao(vetor_parametros_normalizados):
    dados = dados_otimizacao_cache(vetor_parametros_normalizados)
    return 1e10 if dados is None else dados["W0"] / W0

def avaliar_restricoes(vetor_parametros_normalizados):
    dados = dados_otimizacao_cache(vetor_parametros_normalizados)
    if dados is None:
        # Se quebrou, retorna valores muito fora de 1.0
        return np.array([100.0 if op == ">" else -100.0 for _, op, *_ in RESTRICOES_VARREDURA])
        
    # NORMALIZAÇÃO DAS RESTRIÇÕES: Divide pelo limite para igualar a escala
    return np.array([dados[key] / (abs(limite) if limite != 0 else 1.0) for key, _, limite, _, _ in RESTRICOES_VARREDURA])

def calcular_T0_req_cruzeiro(Mach, h, W):
    atm = atmosphere(h)
    rho = atm["density"]
    V = Mach * atm["speed_of_sound"]
    CL = 2 * W / (rho * V**2 * S_w)
    condicao = CondicoesPolar["Cruzeiro"].copy()
    condicao["Mach"] = Mach
    condicao["altitude"] = h
    CD, _, _ = aerodynamics_condicao_e_CL(airplane, condicao, CL)
    T_req = 0.5 * rho * V**2 * S_w * CD
    _, kT = engineTSFC(Mach, h, airplane)
    T0_req = T_req / kT
    return T0_req

def calcular_altitude_tracao_limite(Mach, W):
    try:
        return float(brentq(lambda h: T0 - calcular_T0_req_cruzeiro(Mach, h, W), ALTITUDE_CRUZEIRO_MINIMA, ALTITUDE_CRUZEIRO_MAXIMA))# type: ignore
    except ValueError:
        return np.nan

def calcular_h_max_alcance_altitude_constante(W_i, W_f, Mach): # R = integral de W_f a W_i de SR dW = integral de (V/c)*(L/D)/W dW
    def alcance_breguet(h):
        atm = atmosphere(h)
        rho = atm["density"]
        V = Mach * atm["speed_of_sound"]
        C, _ = engineTSFC(Mach, h, airplane)
        condicao = CondicoesPolar["Cruzeiro"].copy()
        condicao["Mach"] = Mach
        condicao["altitude"] = h
        def SR(W):
            CL = 2 * W / (rho * V**2 * S_w)
            CD, _, _ = aerodynamics_condicao_e_CL(airplane, condicao, CL)
            SR = V * (CL / CD) / (C * W)
            return SR

        R, _ = integrate.quad(SR, W_f, W_i)
        return R

    res = minimize_scalar(lambda h: -alcance_breguet(h), bounds=(ALTITUDE_CRUZEIRO_MINIMA, ALTITUDE_CRUZEIRO_MAXIMA), method="bounded")
    return np.array(res.x)

def calcular_h_max_alcance_LD_constante(W_i, Mach):
    def SR(h):
        atm = atmosphere(h)
        rho = atm["density"]
        V = Mach * atm["speed_of_sound"]
        C, _ = engineTSFC(Mach, h, airplane)
        condicao = CondicoesPolar["Cruzeiro"].copy()
        condicao["Mach"] = Mach
        condicao["altitude"] = h
        # O weight.py avalia o L/D APENAS para o peso de entrada (W_i)
        CL = 2 * W_i / (rho * V**2 * S_w)
        CD, _, _ = aerodynamics_condicao_e_CL(airplane, condicao, CL)
        SR = V * (CL / CD) / C
        return SR

    res = minimize_scalar(lambda h: -SR(h), bounds=(ALTITUDE_CRUZEIRO_MINIMA, ALTITUDE_CRUZEIRO_MAXIMA), method="bounded")
    return np.array(res.x)

# %% PREPARAR O AVIÃO
airplane_name = "Tomav"
airplane = standard_airplane(airplane_name)

# Execute the geometry module to compute all dimensions.
# This updates the airplane dictionary with new entries.
geometry(airplane)

# Execute the weight and thrust estimation
thrust_matching(airplane["inputs"]["W0_guess"], airplane["inputs"]["n_engines"]*airplane["inputs"]["engine"]["Tmax"], airplane, calcular_performance=True)
T0 = airplane["thrust_matching"]["T0"]
W0 = airplane["thrust_matching"]["W0"]
W_empty = airplane["thrust_matching"]["W_empty"]
W_fuel = airplane["thrust_matching"]["W_fuel"]
airplane["inputs"]["W0_guess"] = W0

# Execute the balance analysis
balance(airplane)

# Execute the landing gear analysis
landing_gear(airplane)

# %% PLOTAR A GEOMETRIA
plot_geometry(airplane, az1=0, az2=180, figname="Vista Frontal")

plot_geometry(airplane, az1=90, az2=0, figname="Vista Superior")

plot_geometry(airplane, az1=0, az2=90, figname="Vista Lateral")

# %% CONSTANTES
FRACAO_CL_CLMAX_DECOLAGEM = 1 / 1.41
FRACAO_CL_CLMAX_POUSO = 1 / 1.69
KWARGS_PADRAO = {"delta_xr_w", "x_tank_c_w"}

# RESTRIÇÕES DE OTIMIZAÇÃO
D_TO_MAX = airplane["inputs"]["distance_takeoff"]
D_LD_MAX = airplane["inputs"]["distance_landing"]
SM_FWD_MAX = 0.40
SM_AFT_MIN = 0.05
CLV_MAX = 0.75
TANK_EXC_MIN = 0
MLG_XCG_MARGIN_MIN = 0
ALPHA_TAILSTRIKE_MIN = 10
ALPHA_TIPBACK_MIN = 10
PHI_OVERTURN_MAX = 63
FRAC_NLG_AFT_MIN = 0.04
FRAC_NLG_FWD_MAX = 0.15
ALTITUDE_CRUZEIRO_MINIMA = 10000*ft2m
ALTITUDE_CRUZEIRO_MAXIMA = 43000*ft2m

# Inputs
MLW_frac = airplane["inputs"]["MLW_frac"]
delta_xr_w = airplane['inputs']['delta_xr_w']
AR_w = airplane["inputs"]["AR_w"]
k_korn = airplane["inputs"]["k_korn"]
sweep_w = airplane["inputs"]["sweep_w"]
S_w = airplane["inputs"]["S_w"]  # m^2
W_crew = airplane["inputs"]["W_crew"]
W_payload = airplane["inputs"]["W_payload"]
n_engines = airplane["inputs"]["n_engines"]
trapped_fuel_factor = airplane["fuel_weight"]["trapped_fuel_factor"]
tcr_w = airplane["inputs"]["tcr_w"]
tct_w = airplane["inputs"]["tct_w"]
c_flap_c_wing = airplane["inputs"]["c_flap_c_wing"]
b_flap_b_wing = airplane["inputs"]["b_flap_b_wing"]
c_slat_c_wing = airplane["inputs"]["c_slat_c_wing"]
b_slat_b_wing = airplane["inputs"]["b_slat_b_wing"]
Cht = airplane["inputs"]["Cht"]
Cvt = airplane["inputs"]["Cvt"]
AR_h = airplane["inputs"]["AR_h"]
AR_v = airplane["inputs"]["AR_v"]
taper_w = airplane["inputs"]["taper_w"]
taper_h = airplane["inputs"]["taper_h"]
taper_v = airplane["inputs"]["taper_v"]
sweep_h = airplane["inputs"]["sweep_h"]
sweep_v = airplane["inputs"]["sweep_v"]
tcr_h = airplane["inputs"]["tcr_h"]
tct_h = airplane["inputs"]["tct_h"]
tcr_v = airplane["inputs"]["tcr_v"]
tct_v = airplane["inputs"]["tct_v"]
z_n = airplane["inputs"]["z_n"]
D_n = airplane["inputs"]["D_n"]
D_f = airplane["inputs"]["D_f"]
x_tank_c_w = airplane["inputs"]["x_tank_c_w"]
x_mlg = airplane['inputs']['x_mlg']
y_mlg = airplane['inputs']['y_mlg']
z_mlg = airplane['inputs']['z_lg']
h_TO = airplane["inputs"]["altitude_takeoff"]
deltaISA_TO = airplane["inputs"]["deltaISA_takeoff"]
h_CR = airplane["inputs"]["altitude_cruise"]
Mach_CR = airplane["inputs"]["Mach_cruise"]
h_LD = airplane["inputs"]["altitude_landing"]
deltaISA_LD = airplane["inputs"]["deltaISA_landing"]
Mach_ALT = airplane["inputs"]["Mach_altcruise"]
h_ALT = airplane["inputs"]["altitude_altcruise"]
Mach_MAXCR = airplane["inputs"]["Mach_maxcruise"]
h_MAXCR = airplane["inputs"]["altitude_maxcruise"]
Mach_LOITER = 0.50
h_LOITER = airplane["inputs"]["altitude_loiter"]
SM_fwd = airplane["balance"]["SM_fwd"]
SM_aft = airplane["balance"]["SM_aft"]
CLv = airplane["balance"]["CLv"]
tank_excess = airplane["balance"]["tank_excess"]
alpha_tailstrike = np.degrees(airplane["landing_gear"]["alpha_tailstrike"])
alpha_tipback = np.degrees(airplane["landing_gear"]["alpha_tipback"])
phi_overturn = np.degrees(airplane["landing_gear"]["phi_overturn"])
frac_nlg_aft = airplane["landing_gear"]["frac_nlg_aft"]
frac_nlg_fwd = airplane["landing_gear"]["frac_nlg_fwd"]

# Decolagem
rho_TO = atmosphere(h_TO, deltaISA_TO)["density"]
sigma_TO = rho_TO / 1.225

# Cruzeiro
rho_CR = atmosphere(h_CR)["density"]

# Pouso
rho_LD = atmosphere(h_LD, deltaISA_LD)["density"]
W_LD = W0 * MLW_frac

# Geometria
tc_w_medio = 0.25 * tcr_w + 0.75 * tct_w
h_nac = z_n - D_n/2 - z_mlg
h_fus = -D_f/2 - z_mlg

# Distância vertical do MLG até a parte de baixo da asa (na estação do MLG)
zr_w = airplane["inputs"]["zr_w"]
zt_w = airplane["geometry"]["zt_w"]
yt_w = airplane["geometry"]["yt_w"]
cr_w = airplane["geometry"]["cr_w"]
ct_w = airplane["geometry"]["ct_w"]
frac_mlg = y_mlg / yt_w                                   # fração da semi-envergadura
z_w_mlg = zr_w + (zt_w - zr_w) * frac_mlg                 # linha média da asa em y_mlg (com diedro)
c_w_mlg = cr_w + (ct_w - cr_w) * frac_mlg                 # corda em y_mlg
tc_w_mlg = tcr_w + (tct_w - tcr_w) * frac_mlg             # espessura relativa em y_mlg
z_w_low_mlg = z_w_mlg - c_w_mlg * tc_w_mlg / 2            # intradorso (parte de baixo) da asa em y_mlg
h_mlg_asa = z_w_low_mlg - z_mlg                           # distância vertical MLG -> parte de baixo da asa

# %% CONSTANTES - OTIMIZAÇÃO
RESTRICOES_VARREDURA = [ # (key, operador, limite+MARGEM, label_legenda, cor)
    ("d_TO", ">", D_TO_MAX-50, f"d_TO > {D_TO_MAX-50:.0f} m", "C0"),
    ("d_LD", ">", D_LD_MAX-50, f"d_LD > {D_LD_MAX-50:.0f} m", "C1"),
    ("CLv", ">", CLV_MAX-0.03, f"CLv > {CLV_MAX-0.03:.2f}", "C2"),
    ("tank_excess", "<", TANK_EXC_MIN+0.01, f"tank_excess < {TANK_EXC_MIN+0.01:.2f}", "C3"),
    ("mlg_xcg_margin", "<", MLG_XCG_MARGIN_MIN+0.5, f"mlg_xcg_margin < {MLG_XCG_MARGIN_MIN+0.5:.2f}", "C4"),
    ("SM_fwd", ">", SM_FWD_MAX-0.02, f"SM_fwd > {SM_FWD_MAX-0.03:.2%}", "C5"),
    ("SM_aft", "<", SM_AFT_MIN+0.02, f"SM_aft < {SM_AFT_MIN+0.03:.2%}", "C6"),
    ("alpha_tailstrike", "<", ALPHA_TAILSTRIKE_MIN+0.1, f"alpha_tailstrike < {ALPHA_TAILSTRIKE_MIN+0.1:.1f}°", "C7"),
    ("alpha_tipback", "<", ALPHA_TIPBACK_MIN+0.1, f"alpha_tipback < {ALPHA_TIPBACK_MIN+0.1:.1f}°", "C8"),
    ("phi_overturn", ">", PHI_OVERTURN_MAX-0.1, f"phi_overturn > {PHI_OVERTURN_MAX-0.1:.1f}°", "C9"),
    ("frac_nlg_aft", "<", FRAC_NLG_AFT_MIN+0.01, f"frac_nlg_aft < {FRAC_NLG_AFT_MIN+0.01:.2%}", "C10"),
    ("frac_nlg_fwd", ">", FRAC_NLG_FWD_MAX-0.01, f"frac_nlg_fwd > {FRAC_NLG_FWD_MAX-0.01:.2%}", "C11"),
]

LIMITES_INFERIORES_RESTRICOES = np.array([-np.inf if op == ">" else (limite / (abs(limite) if limite != 0 else 1.0)) for _, op, limite, _, _ in RESTRICOES_VARREDURA])
LIMITES_SUPERIORES_RESTRICOES = np.array([(limite / (abs(limite) if limite != 0 else 1.0)) if op == ">" else np.inf for _, op, limite, _, _ in RESTRICOES_VARREDURA])

PARAMETROS_VARREDURA = [ # (titulo, xlabel, key_input, valores, padrao)
    # Asa
    ("Variação na Posição da Asa", "delta_xr_w (m)", "delta_xr_w", np.linspace(-8, 0, 200), delta_xr_w),
    ("Área da Asa", "S_w (m²)", "S_w", np.linspace(300, 500, 200), S_w),
    ("Enflechamento da Asa", "sweep_w (°)", "sweep_w", np.radians(np.linspace(25, 45, 200)), sweep_w),
    # ("Alongamento da Asa", "AR_w", "AR_w", np.linspace(6, 10.5, 200), AR_w),
    # ("Afilamento da Asa", "taper_w", "taper_w", np.linspace(0.17, 1, 200), taper_w),
    # ("Espessura Relativa da Asa (raiz)", "tcr_w", "tcr_w", np.linspace(0.11, 0.175, 200), tcr_w),
    # ("Espessura Relativa da Asa (ponta)","tct_w", "tct_w", np.linspace(0.08, 0.12, 200), tct_w),
    # ("Corda do Flap",  "c_flap_c_wing", "c_flap_c_wing", np.linspace(0.15, 0.33, 200), c_flap_c_wing),
    ("Envergadura do Flap",  "b_flap_b_wing", "b_flap_b_wing", np.linspace(0.40, 0.65, 200), b_flap_b_wing),
    # ("Corda do Slat",  "c_slat_c_wing", "c_slat_c_wing", np.linspace(0, 0.17, 200), c_slat_c_wing),
    # ("Envergadura do Slat",  "b_slat_b_wing", "b_slat_b_wing", np.linspace(0, 0.90, 200), b_slat_b_wing),
    ("Posição do tanque", "x_tank_c_w", "x_tank_c_w", np.linspace(0.1, 0.7, 200), x_tank_c_w),
    
    # Empenagem Horizontal
    # ("Vol. Estabilizador Horiz.","Cht", "Cht", np.linspace(0.85, 1.2, 200), Cht),
    # ("Alongamento da EH", "AR_h", "AR_h", np.linspace(3.4, 6, 200), AR_h),
    # ("Enflechamento da EH", "sweep_h (°)", "sweep_h", np.radians(np.linspace(20, 45, 200)), sweep_h),
    # ("Afilamento da EH", "taper_h", "taper_h", np.linspace(0.35, 0.60, 200), taper_h),
    # ("Espessura Relativa da EH (raiz)", "tcr_h", "tcr_h", np.linspace(0.09, 0.12, 200), tcr_h),
    # ("Espessura Relativa da EH (ponta)","tct_h", "tct_h", np.linspace(0.08, 0.12, 200), tct_h),
    
    # Empenagem Vertical
    # ("Vol. Estabilizador Vert.", "Cvt", "Cvt", np.linspace(0.075, 0.11, 200), Cvt),
    # ("Alongamento da EV", "AR_v", "AR_v", np.linspace(1.2, 2, 200), AR_v),
    # ("Enflechamento da EV", "sweep_v (°)", "sweep_v", np.radians(np.linspace(30, 50, 200)), sweep_v),
    # ("Afilamento da EV", "taper_v", "taper_v", np.linspace(0.35, 0.60, 200), taper_v),
    # ("Espessura Relativa da EV (raiz)", "tcr_v", "tcr_v", np.linspace(0.095, 0.12, 200), tcr_v),
    # ("Espessura Relativa da EV (ponta)","tct_v", "tct_v", np.linspace(0.08, 0.12, 200), tct_v),
    
    # Voo
    ("Altitude de Cruzeiro", "altitude_cruise (m)", "altitude_cruise", np.linspace(ALTITUDE_CRUZEIRO_MINIMA, ALTITUDE_CRUZEIRO_MAXIMA, 200), h_CR),
    # ("Altitude de Cruzeiro Alternativo", "altitude_altcruise (m)", "altitude_altcruise", np.linspace(ALTITUDE_CRUZEIRO_MINIMA, ALTITUDE_CRUZEIRO_MAXIMA, 200), h_ALT),
    ("Mach de Cruzeiro Alternativo", "Mach_altcruise", "Mach_altcruise", np.linspace(0.5, 0.9, 200), Mach_ALT),
    # ("Altitude de Loiter", "altitude_loiter (m)", "altitude_loiter", np.linspace(1500*ft2m, ALTITUDE_CRUZEIRO_MAXIMA, 200), h_LOITER),
]

LIMITES_PARAMETROS = [(xs[0], xs[-1]) for _, _, _, xs, _ in PARAMETROS_VARREDURA]
VETOR_PARAMETROS_PADRAO = np.array([x_padrao for _, _, _, _, x_padrao in PARAMETROS_VARREDURA])
LIMITES_INFERIORES_PARAMETROS = np.array([b[0] for b in LIMITES_PARAMETROS])
LIMITES_SUPERIORES_PARAMETROS = np.array([b[1] for b in LIMITES_PARAMETROS])

# %% PESOS EM CADA FASE
fases = ["start", "taxi", "takeoff", "climb", "cruise", "descent", "altcruise", "loiter", "landing"]
fracoes = [airplane["fuel_weight"]["Mf_hist"][fase] for fase in fases]
pesos = [0 for _ in range(len(fracoes) + 1)]
pesos[0] = W0
for i, f in enumerate(fracoes):
    pesos[i + 1] = pesos[i] * fracoes[i]
PesosFinais = {fase: peso for fase, peso in zip(fases, pesos[1:])}
PesosIniciais = {fase: peso for fase, peso in zip(fases, pesos[:-1])}

# %% CONDIÇÕES EM CADA FASE DA MISSÃO
CondicoesFase = {
    "start": {
        "Mach": 0.20,
        "altitude": 0,
        "n_engines_failed": 0,
        "highlift_config": "clean",
        "lg_down": 0,
        "h_ground": 0,
    },
    "taxi": {
        "Mach": 0.20,
        "altitude": 0,
        "n_engines_failed": 0,
        "highlift_config": "clean",
        "lg_down": 1,
        "h_ground": 10.668,
    },
    "takeoff": {
        "Mach": 0.30,
        "altitude": 0,
        "n_engines_failed": 0,
        "highlift_config": "takeoff",
        "lg_down": 1,
        "h_ground": 10.668,
    },
    "climb": {
        "Mach": 0.30,
        "altitude": 0,
        "n_engines_failed": 0,
        "highlift_config": "takeoff",
        "lg_down": 0,
        "h_ground": 0,
    },
    "cruise": {
        "Mach": airplane["inputs"]["Mach_cruise"],
        "altitude": airplane["inputs"]["altitude_cruise"],
        "n_engines_failed": 0,
        "highlift_config": "clean",
        "lg_down": 0,
        "h_ground": 0,
    },
    "descent": {
        "Mach": airplane["inputs"]["Mach_cruise"],
        "altitude": airplane["inputs"]["altitude_cruise"],
        "n_engines_failed": 0,
        "highlift_config": "clean",
        "lg_down": 0,
        "h_ground": 0,
    },
    "altcruise": {
        "Mach": airplane["inputs"]["Mach_altcruise"],
        "altitude": airplane["inputs"]["altitude_altcruise"],
        "n_engines_failed": 0,
        "highlift_config": "clean",
        "lg_down": 0,
        "h_ground": 0,
    },
    "loiter": {
        "Mach": 0.50,
        "altitude": airplane["inputs"]["altitude_loiter"],
        "n_engines_failed": 0,
        "highlift_config": "clean",
        "lg_down": 0,
        "h_ground": 0,
    },
    "landing": {
        "Mach": 0.30,
        "altitude": 0,
        "n_engines_failed": 0,
        "highlift_config": "landing",
        "lg_down": 1,
        "h_ground": 10.668,
    },
}

# %% CONDIÇÕES EM CADA FASE DAS POLARES
CondicoesPolar = dict()

# CRUZEIRO
CondicoesPolar["Cruzeiro"] = CondicoesFase["cruise"].copy()
CondicoesPolar["Cruzeiro"]["CL"] = calcular_CL_voo(CondicoesPolar["Cruzeiro"], PesosIniciais["cruise"])

# 2º SEGMENTO
CondicoesPolar["2º Segmento"] = CondicoesFase["takeoff"].copy()
CondicoesPolar["2º Segmento"]["n_engines_failed"] = 1
_, CLmax_TO, _ = aerodynamics_condicao_e_CL(airplane, CondicoesPolar["2º Segmento"], 0)
CL_TO = CLmax_TO * FRACAO_CL_CLMAX_DECOLAGEM
CondicoesPolar["2º Segmento"]["CL"] = CL_TO

# POUSO
CondicoesPolar["Pouso"] = CondicoesFase["landing"].copy()
_, CLmax_LD, _ = aerodynamics_condicao_e_CL(airplane, CondicoesPolar["Pouso"], 0)
CL_LD = CLmax_LD * FRACAO_CL_CLMAX_POUSO
CondicoesPolar["Pouso"]["CL"] = CL_LD

# %% POLARES CD x CL
plt.figure(num="Polares de Arrasto")
plt.title("Polares de Arrasto")
plt.xlabel("CD")
plt.ylabel("CL")

for titulo, condicao in CondicoesPolar.items():
    plotar_curva_polar_condicao(condicao, titulo, LD=False)

plt.legend()

# %% POLARES CL x L/D
plt.figure(num="Eficiências Aerodinâmicas LD")
plt.title("Eficiências Aerodinâmicas L/D")
plt.xlabel("CL")
plt.ylabel("L/D")

for titulo, condicao in CondicoesPolar.items():
    plotar_curva_polar_condicao(condicao, titulo, LD=True)

plt.legend()

# %% BREAKDOWN DO ARRASTO (DRAGDICT + PIE CHART)
_, ax_pie = plt.subplot_mosaic([["Cruzeiro", "2º Segmento", "Pouso"]], constrained_layout=False, num="Breakdown do Arrasto")

for titulo, condicao in CondicoesPolar.items():
    # DRAGDICT
    print(f"\n=== DRAGDICT NO {titulo.upper()} ===")
    CL = condicao["CL"]
    CD, CLmax, dragDict = aerodynamics_condicao_e_CL(airplane, condicao, CL)
    CD0 = dragDict["CD0"]
    dragDictCD = {k: v for k, v in dragDict.items() if isinstance(v, float) and k.startswith("CD")}
    
    for k, v in dragDict.items():
        if k not in dragDictCD:
            print(f"{k} = {v}")
    
    for k, v in dragDictCD.items():
        print(f"{k} = {v * 10000:.0f} counts = {v} ({(v / CD):.2%} do CD, {v / CD0:.2%} do CD0)")

    # PIE CHART
    labels = []
    counts = []
    for k, v in sorted(dragDictCD.items(), key=lambda item: item[1], reverse=True):
        if isinstance(v, float) and k.startswith("CD") and (k not in ("CD", "CD0")) and v != 0:
            labels.append(k)
            counts.append(v)

    wedges, *_ = ax_pie[titulo].pie(counts, autopct=lambda pct: f"{pct:.1f}%" if pct >= 3 else "", startangle=145)
    ax_pie[titulo].axis("equal")
    ax_pie[titulo].legend(wedges, [f"{label}: {count * 10000:.0f} counts ({(count / CD):.2%} do CD)" for label, count in zip(labels, counts)],
                        loc="upper center", bbox_to_anchor=(0.5, -0.03), frameon=False, handlelength=1.2)
    ax_pie[titulo].set_title(f"CD no {titulo}", pad=12)

plt.subplots_adjust(wspace=0.38, bottom=0.26, top=0.88)

# %% MACH x CD NO CRUZEIRO
# Pontos Mach, Mc e Mdd
Mach_CR = CondicoesPolar["Cruzeiro"]["Mach"]
CL_CR = CondicoesPolar["Cruzeiro"]["CL"]
Mdd = k_korn / np.cos(sweep_w) - tc_w_medio / np.cos(sweep_w) ** 2 - CL_CR / (10 * np.cos(sweep_w) ** 3)
Mc = Mdd - (0.1 / 80) ** (1 / 3)

print("\n=== MACH NO CRUZEIRO ===")
print(f"Mach no Cruzeiro = {Mach_CR:.2f}")
print(f"Mach de Divergência = {Mdd:.2f}")
print(f"Mach Crítico = {Mc:.2f}")

Mach_space = np.linspace(0.4, 1, 200)

plt.figure(num="Cruzeiro - Mach x CD, com CL = CL do Cruzeiro")
plt.title(f"Cruzeiro - Mach x CD, com CL = CL do Cruzeiro")
plt.xlabel("Mach")
plt.ylabel("CD")

CD_space = [calcular_CD_Mach_CL_fixo(Mach, CL_CR) for Mach in Mach_space]
plt.plot(Mach_space, CD_space)

for Mach, titulo in zip((Mach_CR, Mdd, Mc), ("Cruzeiro", "Mdd", "Mc")):
    CD_Mach = calcular_CD_Mach_CL_fixo(Mach, CL_CR)
    plt.scatter(Mach, CD_Mach, color=("red" if titulo == "Cruzeiro" else "black"), s=20, zorder=5, label="_")
    plt.annotate(f"{titulo}\nMach={Mach:.2f}\nCD={CD_Mach:.4f}", (Mach, CD_Mach), fontsize=7)

plt.figure(num="Cruzeiro - Mach x CD, com CL variável")
plt.title(f"Cruzeiro - Mach x CD, com CL variável")
plt.xlabel("Mach")
plt.ylabel("CD")

CD_space = [calcular_CD_Mach_CL_variavel(Mach) for Mach in Mach_space]
plt.plot(Mach_space, CD_space)

for Mach, titulo in zip((Mach_CR, Mdd, Mc), ("Cruzeiro", "Mdd", "Mc")):
    CD_Mach = calcular_CD_Mach_CL_variavel(Mach)
    plt.scatter(Mach, CD_Mach, color=("red" if titulo == "Cruzeiro" else "black"), s=20, zorder=5, label="_")
    plt.annotate(f"{titulo}\nMach={Mach:.2f}\nCD={CD_Mach:.4f}", (Mach, CD_Mach), fontsize=7)

# %% DISTÂNCIAS DE DECOLAGEM E POUSO
print("\n=== DISTÂNCIA DE DECOLAGEM (2º SEGMENTO) ===")
d_TO = 0.2387 / (sigma_TO * CLmax_TO) * (W0 / S_w) * (W0 / T0)
print(f"CLmax_TO = {CLmax_TO:.2f}")
print(f"d_TO = {d_TO:.0f} m")
print(f"Limite: {D_TO_MAX} m")

print("\n=== DISTÂNCIA DE POUSO ===")
a_g = 0.5 # Fator de desaceleração média
x_LD = 1.52 / a_g + 1.69
f_LD = 5/3 # Fator de segurança para distância de pouso: (5/3) para FAR 91
A_LD = gravity / (f_LD * x_LD)
h_obstaculo_LD = 15.3 # Altura de obstáculo
B_LD = -10 * gravity * h_obstaculo_LD / x_LD
# W_S_LD = rho_LD * CLmax_LD * (A_LD*d_LD + B_LD)
d_LD = ((W_LD/S_w) / (rho_LD * CLmax_LD) - B_LD) / A_LD

print(f"CLmax_LD = {CLmax_LD:.2f}")
print(f"W_LD = {W_LD/gravity:.1f} kgf")
print(f"d_LD = {d_LD:.0f} m")
print(f"Limite: {D_LD_MAX:.0f} m")

# %% ENFLECHAMENTO x CLMAX, d_TO NO 2º SEGMENTO
enflechamento_space = np.linspace(25, 60, 200)
CLmax_space = []
d_TO_space = []

for enflechamento in enflechamento_space:
    airplane_enflechamento = standard_airplane(airplane_name)
    airplane_enflechamento["inputs"]["sweep_w"] = np.radians(enflechamento)
    geometry(airplane_enflechamento)
    W0_enflechamento, _, _, _ = weight(W0, T0, airplane_enflechamento)
    _, CLmax_enflechamento, _ = aerodynamics_condicao_e_CL(airplane_enflechamento, CondicoesPolar["2º Segmento"], 0)
    d_TO_enflechamento = 0.2387 / (sigma_TO * CLmax_enflechamento) * (W0_enflechamento / S_w) * (W0_enflechamento / T0)
    CLmax_space.append(CLmax_enflechamento)
    d_TO_space.append(d_TO_enflechamento)

_, ax1 = plt.subplots(num="2º Segmento - Enflechamento x CLmax e d_TO")
ax2 = ax1.twinx()

ax1.set_title("2º Segmento - Enflechamento x CLmax e d_TO")
ax1.set_xlabel("Enflechamento (°)")
ax1.set_ylabel("CLmax")
ax2.set_ylabel("d_TO (m)")

ax1.plot(enflechamento_space, CLmax_space, color="#1F77B4", label="CLmax")
ax2.plot(enflechamento_space, d_TO_space, color="#FF7F0E", label="d_TO")

ax1.scatter(np.degrees(sweep_w), CLmax_TO, color="red", s=20, zorder=5)
ax1.annotate(f"Padrão: {np.degrees(sweep_w):.1f}°\nCLmax={CLmax_TO:.2f}", (np.degrees(sweep_w), CLmax_TO), fontsize=7)

ax2.scatter(np.degrees(sweep_w), d_TO, color="red", s=20, zorder=5)
ax2.annotate(f"Padrão: {np.degrees(sweep_w):.1f}°\nd_TO={d_TO:.0f}m", (np.degrees(sweep_w), d_TO), fontsize=7)


# %% BREAKDOWN DO PESO
print("\n=== BREAKDOWN DO PESO ===")
print(f"W0 = {W0 / gravity:.1f} kgf  ({W0 / W0:.2%} do MTOW)")
print(f"W_empty = {W_empty / gravity:.1f} kgf ({W_empty / W0:.2%} do MTOW)")
print(f"W_fuel = {W_fuel / gravity:.1f} kgf ({W_fuel / W0:.2%} do MTOW)")
print(f"W_payload = {W_payload / gravity:.1f} kgf ({W_payload / W0:.2%} do MTOW)")
print(f"W_crew = {W_crew / gravity:.1f} kgf ({W_crew / W0:.2%} do MTOW)")

print("\n=== BREAKDOWN DO PESO VAZIO ===")
for k, v in sorted(airplane["empty_weight"].items(), key=lambda item: item[1], reverse=True):
    if k.startswith("xcg"):
        pass
    else:
        print(f"{k}: {v / gravity:.1f} kgf ({v / W0:.2%} do MTOW, {v / W_empty:.2%} de W_empty)")
        k_xcg = k.replace("W", "xcg", 1)
        print(f"{k_xcg} = {airplane['empty_weight'][k_xcg]:.2f} m")

print("\n=== DESEMPENHO EM CADA FASE ===")
print(f"L/D em cruzeiro: {airplane['fuel_weight']['LD_hist']['cruise']:.2f}")
print(f"TSFC em cruzeiro: {airplane['fuel_weight']['C_hist']['cruise'] * 3600:.2f}")
print(f"L/D em altcruise: {airplane['fuel_weight']['LD_hist']['altcruise']:.2f}")
print(f"TSFC em altcruise: {airplane['fuel_weight']['C_hist']['altcruise'] * 3600:.2f}")
print(f"L/D em loiter: {airplane['fuel_weight']['LD_hist']['loiter']:.2f}")
print(f"TSFC em loiter: {airplane['fuel_weight']['C_hist']['loiter'] * 3600:.2f}")

print("\n=== PESOS EM CADA FASE ===")
print(f"Peso Inicial (MTOW): {W0 / gravity:.2f} kgf")
for i, fase in enumerate(PesosFinais.keys()):
    print(f"{i}. {fase}:")
    print(f"    -> Peso Final: {PesosFinais[fase]/gravity:.2f} kgf ({PesosFinais[fase]/W0:.2%} do MTOW)")
    diferenca = PesosFinais[fase] - PesosIniciais[fase]
    print(f"    -> Diferença: {diferenca/gravity:.2f} kgf ({diferenca/W0:.2%} do MTOW, {diferenca/W_fuel:.2%} do W_fuel)")

W_trapped_fuel = (W0 - PesosFinais["landing"])*(trapped_fuel_factor - 1)
print(f"* Trapped fuel: {W_trapped_fuel/gravity:.2f} kgf ({W_trapped_fuel/W0:.2%} do MTOW, {W_trapped_fuel/W_fuel:.2%} do W_fuel)")

print("\n=== L/D MAX EM CADA FASE ===")
for i, fase in enumerate(PesosFinais.keys()):
    condicao = CondicoesFase[fase]
    peso_inicio = PesosIniciais[fase]
    peso_fim = PesosFinais[fase]
    LD_max, CL_ld_max, CLmax = calcular_ld_max(airplane, condicao)
    velocidade_inicio = velocidade_em_ld_max(airplane, condicao, peso_inicio, CL_ld_max)
    velocidade_fim = velocidade_em_ld_max(airplane, condicao, peso_fim, CL_ld_max)
    LD_inicio, CL_inicio, _ = calcular_ld_operacional(airplane, condicao, peso_inicio)
    LD_fim, CL_fim, _ = calcular_ld_operacional(airplane, condicao, peso_fim)
    aviso_inicio = " (fora do CLmax!)" if CL_inicio > CLmax else ""
    aviso_fim = " (fora do CLmax!)" if CL_fim > CLmax else ""
    print(f"{i}. {fase}:")
    print(f"    -> Início: CL={CL_inicio:.2f}{aviso_inicio}, V@L/Dmax={velocidade_inicio:.1f} m/s, L/D={LD_inicio:.2f}")
    print(f"    -> Fim:    CL={CL_fim:.2f}{aviso_fim}, V@L/Dmax={velocidade_fim:.1f} m/s, L/D={LD_fim:.2f}")
    print(f"    -> V operacional={condicao['Mach']*atmosphere(condicao['altitude'])['speed_of_sound']:.1f} m/s, CL_max={CLmax:.2f}, L/D max={LD_max:.2f} em CL={CL_ld_max:.2f}")

# %% ALONGAMENTO x PESO
AR_space = np.linspace(0.8*AR_w, 1.2*AR_w, 200)
W0_space = []
We_space = []
Wf_space = []
for AR in AR_space:
    airplane_alongamento = standard_airplane(airplane_name)
    airplane_alongamento["inputs"]["AR_w"] = AR
    geometry(airplane_alongamento)
    W0_AR, W_empty_AR, W_fuel_AR, _ = weight(W0, T0, airplane_alongamento)
    W0_space.append(W0_AR / gravity)
    We_space.append(W_empty_AR / gravity)
    Wf_space.append(W_fuel_AR / gravity)

AR_min = AR_space[np.argmin(W0_space)]
W0_min = min(W0_space)

_, ax1 = plt.subplots(num="Alongamento x Peso")
plt.title("Alongamento x Peso")

ax2 = ax1.twinx()
ax1.autoscale(enable=True, axis="both", tight=True)
ax2.autoscale(enable=True, axis="both", tight=True)

ax1.set_xlabel("Alongamento")
ax1.set_ylabel("Peso (kgf)")
ax2.set_ylabel("Peso (kgf)")

ax1.plot(AR_space, W0_space, label="MTOW", color="#1F77B4")
ax2.plot(AR_space, Wf_space, label="Fuel Weight", color="#FF7F0E")
ax2.plot(AR_space, We_space, label="Empty Weight", color="#2CA02C")

ax1.scatter(AR_min, W0_min, color="black", s=20, zorder=5, label="_")
ax1.annotate(f"Peso mínimo:\nAR={AR_min:.2f}\nW0={W0_min:.2f} kgf", (AR_min, W0_min), fontsize=7)

ax1.scatter(AR_w, W0 / gravity, color="red", s=20, zorder=5, label="_")
ax1.annotate(f"Alongamento padrão:\nAR={AR_w:.2f}\nW0={W0/gravity:.2f} kgf", (AR_w, W0 / gravity), fontsize=7)

plt.legend(["MTOW", "Empty Weight", "Fuel Weight"])
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")

# %% AREA x PESO E TRAÇÃO REQUERIDA
area_space = np.linspace(int(0.8 * S_w), int(1.2 * S_w), 200)
W0_space = []
T0reqs_space = []
deltaS_wlan_space = []
for S in area_space:
    airplane_area = standard_airplane(airplane_name)
    airplane_area["inputs"]["S_w"] = S
    geometry(airplane_area)
    thrust_matching(W0, T0, airplane_area, calcular_performance=True)
    W0_space.append(airplane_area["thrust_matching"]["W0"] / gravity)
    T0reqs_space.append(airplane_area["thrust_matching"]["T0req"])
    deltaS_wlan_space.append(airplane_area["thrust_matching"]["deltaS_wlan"])

_, ax1 = plt.subplots(num="Área da Asa x W0 e T0_req")
ax2 = ax1.twinx()

ax1.set_title("Área da Asa x W0 e T0_req")
ax1.set_xlabel("S (m^2)")
ax1.set_ylabel("W0 (kgf)")
ax2.set_ylabel("T0_req (kN)")

ax1.plot(area_space, W0_space, label="W0")
ax1.axvline(area_space[np.argmin(np.abs(deltaS_wlan_space))], color="red", linestyle="--", label="deltaS_wlan = 0")
ax1.scatter(S_w, W0/gravity, color="red", s=20, zorder=5)
ax1.annotate(f"Padrão: S={S_w:.1f} m^2\nW0={W0/gravity:.1f} kgf", (S_w, W0/gravity), fontsize=7)

for key in T0reqs_space[0].keys():
    T0req_key_space = [T0reqs[key]/1000 for T0reqs in T0reqs_space]
    ax2.plot(area_space, T0req_key_space, label=f"T0req {key}", linestyle="--")

ax2.axhline(T0/1000, color="green", label="T0 máximo")
ax2.scatter(S_w, T0/1000, color="red", s=20, zorder=5, label="_")
ax2.annotate(f"Padrão: S={S_w:.1f} m^2\nT0={T0/1000:.1f} kN", (S_w, T0/1000), fontsize=7)

handles1, labels1 = ax1.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(handles1 + handles2, labels1 + labels2, loc="upper right")


# %% MACH X ALTITUDE IDEAL DE CRUZEIRO
Mach_space = np.linspace(0.50, 1, 100)

casos_cruzeiro = {
    "cruise":    (PesosIniciais["cruise"],    PesosFinais["cruise"],    "C0"),
    "altcruise": (PesosIniciais["altcruise"], PesosFinais["altcruise"], "C1"),
    "maxcruise": (PesosIniciais["cruise"], PesosFinais["cruise"], "C3"),
}

condicoes_cruzeiro = {
    "cruise":    (Mach_CR,     h_CR,     "C0"),
    "altcruise": (Mach_ALT,    h_ALT,    "C1"),
    "maxcruise": (Mach_MAXCR,  h_MAXCR,  "C3"),
    "loiter":    (Mach_LOITER, h_LOITER, "C2"),
}

plt.figure(num="Mach x Altitude Ideal de Cruzeiro")
plt.title("Mach x Altitude Ideal de Cruzeiro")
plt.xlabel("Mach")
plt.ylabel("Altitude (ft)")

print("\n=== ALTITUDES IDEAIS DE CRUZEIRO ===")
for label, (W_i, W_f, cor) in casos_cruzeiro.items():
    if label != "maxcruise":
        h_max_alcance   = [calcular_h_max_alcance_altitude_constante(W_i, W_f, Mach) / ft2m for Mach in Mach_space]
        h_weight_cruise = [calcular_h_max_alcance_LD_constante(W_i, Mach) / ft2m for Mach in Mach_space]
        h_tracao_limite = [calcular_altitude_tracao_limite(Mach, W_i) / ft2m for Mach in Mach_space]
        plt.plot(Mach_space, h_max_alcance, color=cor, linestyle=":", linewidth=1, label=f"Alcance máx. com altitude constante – {label}")
        plt.plot(Mach_space, h_weight_cruise, color=cor,  linewidth=1, label=f"Alcance máx. com L/D constante (Cruise Climb) - {label}")
        plt.plot(Mach_space, h_tracao_limite, color=cor, linestyle="--", label=f"Tração limite – {label}")
        
    Mach_atual, h_atual, _ = condicoes_cruzeiro[label]
    print(f"{label}: Mach={Mach_atual:.2f}\n"
          f"    -> Altitude atual = {h_atual/ft2m:.0f} ft\n"
          f"    -> Altitude ideal (Cruise Climb com L/D constante) = {calcular_h_max_alcance_LD_constante(W_i, Mach_atual)/ft2m:.0f} ft\n"
          f"    -> Altitude ideal (altitude constante) = {calcular_h_max_alcance_altitude_constante(W_i, W_f, Mach_atual)/ft2m:.0f} ft\n"
          f"    -> Altitude limite de tração = {calcular_altitude_tracao_limite(Mach_atual, W_i)/ft2m:.0f} ft")



for label, (Mach, h, cor) in condicoes_cruzeiro.items():
    plt.scatter(Mach, h / ft2m, color=cor, s=40, zorder=5,
               label=f"{label}: Mach={Mach:.2f}, {h/ft2m:.0f} ft")

plt.axhline(ALTITUDE_CRUZEIRO_MAXIMA / ft2m, color="gray",  linewidth=2, label="Teto de Serviço = 43000 ft")
plt.axhline(ALTITUDE_CRUZEIRO_MINIMA / ft2m, color="gray", linestyle="-.", linewidth=2, label="Altitude Mínima em Rota = 10000 ft")
plt.legend(loc="lower right")


# %% T0 REQUERIDO X ALTITUDE
h_space = np.linspace(1000*ft2m, 50000*ft2m, 200)

plt.figure(num="Altitude x T0 requerido")
plt.title("Altitude x T0 requerido")
plt.xlabel("T0 (kN)")
plt.ylabel("Altitude (ft)")

for label, (Mach, h_atual, cor) in condicoes_cruzeiro.items():
    T0_req_space = np.array([calcular_T0_req_cruzeiro(Mach, h, PesosIniciais[label] if label != "maxcruise" else PesosIniciais["cruise"]) for h in h_space])
    plt.plot(T0_req_space/1000, h_space/ft2m, color=cor, linestyle="--", label=f"{label} (Mach={Mach:.2f})")
    plt.scatter(calcular_T0_req_cruzeiro(Mach, h_atual, PesosIniciais[label] if label != "maxcruise" else PesosIniciais["cruise"])/1000, h_atual / ft2m, color=cor, s=40, zorder=5, label=f"{label} atual ({h_atual/ft2m:.0f} ft)")

plt.axvline(T0/1000, color="black", linewidth=2, label=f"T0 atual = {T0/1000:.0f} kN")
plt.axhline(ALTITUDE_CRUZEIRO_MAXIMA / ft2m, color="gray", linewidth=2, label="Teto de Serviço = 43000 ft")
plt.axhline(ALTITUDE_CRUZEIRO_MINIMA / ft2m, color="gray", linestyle="-.", linewidth=2, label="Altitude Mínima em Rota = 10000 ft")
plt.legend(loc="lower right")

# %% DIAGRAMA DE PROJETO
W0_space = np.linspace(0.3*W0, 2*W0, 500)

# DECOLAGEM -> T/W >= 0.2387/(sigma_TO*CLmax_TO*d_TO_max)*W/Sw
T0_W0_TO_space = 0.2387 / (sigma_TO * CLmax_TO * D_TO_MAX) * W0_space / S_w

# SUBIDA -> T/W >= ks^2/CLmax*CD0 + CLmax/ks^2*K + gamma
CLmax_CL = CLmax_TO
ks_CL = 1.2
if n_engines == 2:
    G_CL = 2.4 / 100
elif n_engines == 3:
    G_CL = 2.7 / 100
elif n_engines == 4:
    G_CL = 3.0 / 100
else:
    raise ValueError
W_CL = PesosIniciais["climb"]
condicao_climb = CondicoesFase["climb"].copy()
condicao_climb["n_engines_failed"] = 1
condicao_climb["CL"] = CLmax_CL / ks_CL**2
_, _, dragDict_CL = aerodynamics_condicao_e_CL(airplane, condicao_climb, condicao_climb["CL"])
CD0_CL = dragDict_CL["CD0"]
K_CL = dragDict_CL["K"]
T_W_CL = ks_CL ** 2 / CLmax_CL * CD0_CL + CLmax_CL / ks_CL ** 2 * K_CL + G_CL
T0_W0_CL = (W_CL / W0) * (n_engines / (n_engines - 1)) * T_W_CL

# CRUZEIRO
_, _, dragDict_CR = aerodynamics_condicao_e_CL(airplane, CondicoesPolar["Cruzeiro"], CondicoesPolar["Cruzeiro"]["CL"])
CDwave_CR = dragDict_CR["CDwave"]
CD0_CR = dragDict_CR["CD0"]
K_CR = dragDict_CR["K"]
_, kT = engineTSFC(Mach_CR, h_CR, airplane)
T_CR = kT * T0
V_CR = Mach_CR * atmosphere(h_CR)["speed_of_sound"]
q_CR = rho_CR * V_CR ** 2 / 2
W_CR_inicial = PesosIniciais["cruise"]
W_S_CR_space = (W_CR_inicial / W0) * W0_space / S_w
T_W_CR_space = q_CR / W_S_CR_space * (CD0_CR+CDwave_CR) + W_S_CR_space / q_CR * K_CR
T0_W0_CR_space = (W_CR_inicial / W0) / (T_CR / T0) * T_W_CR_space

# POUSO
W_S_LD = rho_LD * CLmax_LD * (A_LD * D_LD_MAX + B_LD)
W0_S_LD = W_S_LD / (W_LD / W0)

# PLOT W0/S x T0/W0
plt.figure(num="Diagrama de Projeto")
plt.title("Diagrama de Projeto")
plt.xlabel("W0/S (kN/m^2)")
plt.ylabel("T0/W0")

plt.plot(W0_space / S_w / 1000, T0_W0_TO_space, color="#1F77B4", label=f"Decolagem (d_TO={D_TO_MAX:.0f} m)")
plt.plot(W0_space / S_w / 1000, T0_W0_CR_space, color="#FF7F0E", label=f"Cruzeiro (Mach={Mach_CR:.2f}, altitude={h_CR/ft2m:.0f} ft)")
plt.axhline(T0_W0_CL, color="#9467BD", label=f"Subida (OEI)")
plt.axvline(W0_S_LD / 1000, color="#2CA02C", label=f"Pouso (d_LD={D_LD_MAX:.0f} m)")
plt.xlim(min(W0_space / S_w / 1000), max(W0_space / S_w / 1000))

# PONTO DE PROJETO
plt.scatter(W0 / S_w / 1000, T0 / W0, color="red", marker="x", label=f"{airplane_name}")
plt.annotate(f"Ponto de projeto:\nW0/S = {W0/S_w/1000:.1f} kN/m^2\nT0/W0 = {T0/W0:.2f}", (W0/S_w/1000, T0/W0))
plt.legend()

# PLOT S_req x T0_req
plt.figure(num="Diagrama de Projeto (Tração e Área requeridas)")
plt.title("Diagrama de Projeto (Tração e Área requeridas)")
plt.xlabel("S_req (m^2)")
plt.ylabel("T0_req (kN)")

plt.plot(W0/(W0_space/S_w), T0_W0_TO_space * W0 / 1000, color="#1F77B4", label=f"Decolagem (d_TO={D_TO_MAX:.0f} m)")
plt.plot(W0/(W0_space/S_w), T0_W0_CR_space * W0 / 1000, color="#FF7F0E", label=f"Cruzeiro (Mach={Mach_CR:.2f}), altitude={h_CR/ft2m:.0f} ft")
plt.axhline(T0_W0_CL * W0 / 1000, color="#9467BD", label=f"Subida (OEI)")
plt.axvline(W0 / W0_S_LD, color="#2CA02C", label=f"Pouso (d_LD={D_LD_MAX:.0f} m)")
plt.xlim(min(W0/(W0_space/S_w)), max(W0/(W0_space/S_w)))

# PONTO DE PROJETO
plt.scatter(S_w, T0/1000, color="red", marker="x", label=f"{airplane_name}")
plt.annotate(f"Ponto de projeto:\nS = {S_w:.1f} m^2\nT0 = {T0/1000:.1f} kN", (S_w, T0/1000))
plt.legend()

# %% RESTRIÇÕES
print("\n=== RESTRIÇÕES ===")
print(f"d_TO = {d_TO:.0f} m (máximo: {D_TO_MAX:.0f} m)")
print(f"d_LD = {d_LD:.0f} m (máximo: {D_LD_MAX:.0f} m)")
print(f"CLv = {CLv:.3f} (máximo: {CLV_MAX:.3f})")
print(f"Tank excess = {tank_excess:.3f} (mínimo: {TANK_EXC_MIN:.3f})")
print(f"x_mlg - xcg_aft = {x_mlg - airplane['balance']['xcg_aft']:.3f} m (mínimo: {MLG_XCG_MARGIN_MIN:.3f} m)")
print(f"SM_fwd = {SM_fwd:.2%} (máximo: {SM_FWD_MAX:.2%})")
print(f"SM_aft = {SM_aft:.2%} (mínimo: {SM_AFT_MIN:.2%})")
print(f"Ground Clearance da nacele: {h_nac:.2f} m (mínimo: 0.75 m)")
print(f"Ground Clearance da fuselagem: {h_fus:.2f} m (mínimo: 0.75 m)")
print(f"Distância vertical do MLG até a parte de baixo da asa: {h_mlg_asa:.2f} m")
print(f"Ângulo de tailstrike = {alpha_tailstrike:.2f}° (mínimo: {ALPHA_TAILSTRIKE_MIN:.2f}°)")
print(f"Ângulo de tipback = {alpha_tipback:.2f}° (mínimo: {ALPHA_TIPBACK_MIN:.2f}°)")
print(f"Ângulo de overturn = {phi_overturn:.2f}° (máximo: {PHI_OVERTURN_MAX:.2f}°)")
print(f"Fração do peso no NLG (aft) = {frac_nlg_aft:.2%} (mínimo: {FRAC_NLG_AFT_MIN:.2%})")
print(f"Fração do peso no NLG (fwd) = {frac_nlg_fwd:.2%} (máximo: {FRAC_NLG_FWD_MAX:.2%})")

# %% ANALYZE
print("\n=== ANALYZE ===")
analyze(airplane, print_log=True, plot=True)

# %% OTIMIZAÇÃO
print("\n=== OTIMIZAÇÃO ===")
cache_otimizacao = {}
resultados_parametros_otimizados = minimize(
    fun=funcao_objetivo_otimizacao,
    x0=normalizar_vetor_parametros(VETOR_PARAMETROS_PADRAO),
    method="SLSQP",
    bounds=[(0.0, 1.0)] * len(LIMITES_PARAMETROS),
    constraints=NonlinearConstraint(avaliar_restricoes, lb=LIMITES_INFERIORES_RESTRICOES, ub=LIMITES_SUPERIORES_RESTRICOES),
    options={"ftol": 1e-7, "disp": True, "maxiter": 500, "eps": 1e-3}, 
)

resultados_parametros_otimizados.x = desnormalizar_vetor_parametros(resultados_parametros_otimizados.x)
dados_airplane_otimizado = calcular_dados_otimizacao(construir_airplane_parametros(resultados_parametros_otimizados.x))

if dados_airplane_otimizado is not None:
    print(f"\nResultados do avião otimizado (não é o avião atual!):")
    print(f"  W0 = {dados_airplane_otimizado['W0'] / gravity:.1f} kgf  | Atual: {W0 / gravity:.1f} kgf | Redução = {(W0 - dados_airplane_otimizado['W0']) / W0:.2%})")
    print(f"  W_empty = {dados_airplane_otimizado['W_empty'] / gravity:.1f} kgf  | Atual: {W_empty / gravity:.1f} kgf | Redução = {(W_empty - dados_airplane_otimizado['W_empty']) / W_empty:.2%})")
    print(f"  W_fuel = {dados_airplane_otimizado['W_fuel'] / gravity:.1f} kgf  | Atual: {W_fuel / gravity:.1f} kgf | Redução = {(W_fuel - dados_airplane_otimizado['W_fuel']) / W_fuel:.2%})")
    print(f"  d_TO = {dados_airplane_otimizado['d_TO']:.0f} m  | Atual: {d_TO:.0f} m | Aumento = {-(d_TO - dados_airplane_otimizado['d_TO']) / d_TO:.2%})")
    print(f"  d_LD = {dados_airplane_otimizado['d_LD']:.0f} m  | Atual: {d_LD:.0f} m | Aumento = {-(d_LD - dados_airplane_otimizado['d_LD']) / d_LD:.2%})")
    print(f"Parâmetros do avião otimizado (não é o avião atual!):")
    for i, (titulo, xlabel, key, xs, x_padrao) in enumerate(PARAMETROS_VARREDURA):
        parametro = resultados_parametros_otimizados.x[i]
        print(f"  {key} = {np.degrees(parametro) if key.startswith('sweep') else parametro:.4f} | Atual = {np.degrees(x_padrao) if key.startswith('sweep') else x_padrao:.4f} | Diferença = {(parametro - x_padrao) / x_padrao:.2%}")
else:
    print("Erro na otimização.")

# %% VARREDURA DE PARÂMETROS
dados_airplane_padrao = calcular_dados_otimizacao(airplane)
if dados_airplane_padrao is None:
    raise ValueError

print("\n=== VARREDURA DE PARÂMETROS ===")
for titulo, xlabel, key, xs, x_padrao in PARAMETROS_VARREDURA:
    print(f"Varrendo {titulo}...")
    _, ax1 = plt.subplots(num=f"Otimização - {titulo}")
    ax2 = ax1.twinx()
    ax1.set_title(f"Otimização - {titulo}")
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("W0 (kgf)", color="#1F77B4")
    ax2.set_ylabel("Distância (m)")
    ax1.tick_params(axis="y", labelcolor="#1F77B4")
    
    dados_airplanes = []
    for x in xs:
        if key in KWARGS_PADRAO:
            airplane_parametro = standard_airplane(airplane_name, **{key: x})
        else:
            airplane_parametro = standard_airplane(airplane_name)
            airplane_parametro["inputs"][key] = x
        dados_airplane_parametro = calcular_dados_otimizacao(airplane_parametro)
        dados_airplanes.append(dados_airplane_parametro)

    sombrear_mascara(ax1, xs, [d is None for d in dados_airplanes], color="gray", alpha=0.30, label="Não converge")
    for key_requisito, op, limite, label, cor in RESTRICOES_VARREDURA:
        mascara = np.array([(dados_airplanes[i] is not None) and (dados_airplanes[i][key_requisito] > limite if op == ">" else dados_airplanes[i][key_requisito] < limite) for i in range(len(xs))], dtype=bool)
        sombrear_mascara(ax1, xs if not key.startswith("sweep") else np.degrees(xs), mascara, color=cor, alpha=0.20, label=label)

    ax1.plot(xs if not key.startswith("sweep") else np.degrees(xs), [d['W0'] / gravity if d else np.nan for d in dados_airplanes], color="#1F77B4", label="W0 (kgf)")
    ax2.plot(xs if not key.startswith("sweep") else np.degrees(xs), [d['d_TO'] if d else np.nan for d in dados_airplanes], color="#FF7F0E", label="d_TO (m)")
    ax2.plot(xs if not key.startswith("sweep") else np.degrees(xs), [d['d_LD'] if d else np.nan for d in dados_airplanes], color="#BCBD22", label="d_LD (m)")
    ax2.axhline(D_TO_MAX, color="#FF7F0E", linestyle="--", linewidth=1, label=f"Limite d_TO = {D_TO_MAX:.0f} m")
    ax2.axhline(D_LD_MAX, color="#BCBD22", linestyle="--", linewidth=1, label=f"Limite d_LD = {D_LD_MAX:.0f} m")

    ax1.scatter(x_padrao if not key.startswith("sweep") else np.degrees(x_padrao), dados_airplane_padrao['W0'] / gravity, color="black", s=30, zorder=6)
    ax1.annotate(
        f"Atual: W0={dados_airplane_padrao['W0'] / gravity:.0f} kgf, {xlabel}={x_padrao if not key.startswith('sweep') else np.degrees(x_padrao):.2f}\n"
        f"d_TO={dados_airplane_padrao['d_TO']:.0f} m, d_LD={dados_airplane_padrao['d_LD']:.0f} m\n"
        f"SM_fwd={dados_airplane_padrao['SM_fwd']:.2%}, SM_aft={dados_airplane_padrao['SM_aft']:.2%}\n"
        f"CLv={dados_airplane_padrao['CLv']:.3f}, tank_excess={dados_airplane_padrao['tank_excess']:.2f}\n"
        f"alpha_tailstrike={dados_airplane_padrao['alpha_tailstrike']:.2f}°, alpha_tipback={dados_airplane_padrao['alpha_tipback']:.2f}°, phi_overturn={dados_airplane_padrao['phi_overturn']:.2f}°\n"
        f"frac_nlg_aft={dados_airplane_padrao['frac_nlg_aft']:.2%}, frac_nlg_fwd={dados_airplane_padrao['frac_nlg_fwd']:.2%}, mlg_xcg_margin={dados_airplane_padrao['mlg_xcg_margin']:.3f} m\n",
        xy=(x_padrao if not key.startswith("sweep") else np.degrees(x_padrao), dados_airplane_padrao['W0'] / gravity), fontsize=6, ha="right", zorder=7)
    ax2.scatter(x_padrao if not key.startswith("sweep") else np.degrees(x_padrao), dados_airplane_padrao['d_TO'], color="black", s=30, zorder=6)
    ax2.scatter(x_padrao if not key.startswith("sweep") else np.degrees(x_padrao), dados_airplane_padrao['d_LD'], color="black", s=30, zorder=6, marker="^")

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc="upper left")

# %% DICIONÁRIO DO AVIÃO
with open(f"{airplane_name.lower()}.json", "w") as jsonfile:
    json.dump(airplane, jsonfile, indent=4)
    
flattened_airplane = flatten_dict(airplane)
with open(f"{airplane_name.lower()}.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    for key, value in flattened_airplane.items():
        writer.writerow([key, value])
        
print(f"\nDicionário do avião salvo em '{airplane_name.lower()}.csv' e '{airplane_name.lower()}.json'.")

# %% SALVAR E MOSTRAR FIGURAS
print("Salvando figuras em 'Figuras/'...")
os.makedirs("Figuras", exist_ok=True)
for fig_num in plt.get_fignums():
    fig = plt.figure(fig_num)
    title = fig.canvas.manager.get_window_title()  # type: ignore
    fig.savefig(f"Figuras/{title.replace('/','')}.png", dpi=300, bbox_inches="tight")

print("Todas as figuras salvas.")
print("Abrindo figuras.")
plt.show()