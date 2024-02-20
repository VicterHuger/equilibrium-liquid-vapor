""" Mixing Rules functions"""
import math
from typing import List
from service.helpers.pr_parameters import solve_pr_cubic_equation
from service.helpers.uniquac import calculate_gibbs_excess_activity_coefficients
from utils.constants import GAS_CONSTANT, THRESHOLD, GAS_CONSTANT_SI
from utils.phase import Phase
from utils.uniquac_properties import UniquacProperties


class MixParams():
    """Object Class with mixing params
    """

    def __init__(self, av: float, amm: float, clamb: float, alf: List[float]) -> None:
        self.av = av
        self.amm = amm
        self.clamb = clamb
        self.alf = alf


def calculate_heidemann_kokal_elv(number_components: int, number_molecular_components: int, system_phase: Phase, temperature: float, pressure: float, liquid_compositions: List[float], atractive_params: List[float], repulsive_params: List[float], activity_coefficients: List[float]):
    """calculate Heidemann kokal liquid vapor equilibrium values

    Args:
        number_components (int): _description_
        number_molecular_components (int): _description_
        system_phase (Phase): _description_
        temperature (float): _description_
        pressure (float): _description_
        liquid_compositions (List[float]): _description_
        atractive_params (List[float]): _description_
        repulsive_params (List[float]): _description_
        activity_coefficients (List[float]): _description_

    Returns:
        (volume, fugacities): tuple(float, List[float])
    """
    # Cálculo de ELV pela equação PR via metodologia de Heidemann-Kokal (1990) (PR + UNIQUAC + Debye-Hückel)
    # SUBROUTINE HKOKAL

    # VARIÁVEIS PRINCIPAIS --------------------------------------------------
    # NC    ------ Número de componentes (espécies) no sistema
    # NM    ------ Número de espécies moleculares no sistema
    # GAMA(I) ---- Coeficiente de atividade do componente I
    # RG    ------ Constante universal do gases
    # T     ------ Temperatura (K)
    # P     ------ Pressão (bar)
    # X(I)  ------ Fração molar do componente I
    # A(I)  ------ Parâmetro de energia do componente I puro
    # AM    ------ Parâmetro de energia da mistura
    # B(I)  ------ Co-volume do componente I puro
    # BM    ------ Co-volume da mistura
    # ALFA(I) ---- a/bRT para o componente I puro
    # ALFAM ------ a/bRT para a mistura
    # QSI(I) ----- b/v para o componente I puro
    # QSIM  ------ b/v para a mistura
    # FI(I) ------ Coeficiente de fugacidade do componente I na mistura
    # V_p   ------ Volume molar da fase
    # DGE   ------ Variação da energia livre de Gibbs na mistura
    # Fug(I) ----- Fugacidade do componente I na mistura

    # VARIÁVEIS AUXILIARES --------------------------------------------------
    # ALFAMIN
    # C1
    # C2
    # D1(I)
    # D2(I)
    # BETA
    # DELTA
    # F1
    # F2
    # G(I)
    # SOMA1
    # H(I)
    # GM
    # HM
    # FMIN
    # QSIMIN
    # DELTAF
    # D3
    # D4(I)
    # D5(I)
    # ALNFI(I)
    # SOMA0
    # GX, HX e BX
    # CCC

    # CÁLCULO DE DGE
    delta_mixture_gibbs_energy = 0.0
    for i in range(number_components):
        delta_mixture_gibbs_energy = GAS_CONSTANT * \
            temperature*liquid_compositions[i] * \
            math.log(activity_coefficients[i])

    # NORMALIZAÇÃO DAS FRAÇÕES MOLARES --------------------------------------
    sum0 = 0.0
    for i in range(number_molecular_components):
        sum0 = sum0 + liquid_compositions[i]
    for i in range(number_molecular_components):
        liquid_compositions[i] = liquid_compositions[i]/sum0

    # MÉTODO DE HEIDEMANN-KOKAL ---------------------------------------------

    # CÁLCULO DE QSI(I)
    c1 = 1.0 + math.sqrt(2)
    c2 = 1.0 - math.sqrt(2)
    alpha_min = 1.15*((2.0 + c1 + c2) + 2.0*math.sqrt((c1+1.0)*(c2+1.0)))
    beta = - 0.030082157186
    delta = 0.003041741701
    alpha: List[float] = []
    d1: List[float] = []
    d2: List[float] = []
    qsi: List[float] = []

    for i in range(number_molecular_components):
        alpha[i] = atractive_params[i] / \
            (repulsive_params[i]*GAS_CONSTANT*temperature)
        d1[i] = (c1 + c2) - alpha[i]
        d2[i] = c1*c2 + alpha[i]
        if alpha[i] >= alpha_min:
            qsi[i] = (- d1[i] + math.sqrt(d1[i]**2.0-4.0*d2[i]))/(2.0*d2[i])
        else:
            qsi[i] = 1.0 + beta*alpha[i]**2.0 + delta*alpha[i]**3.0

    # CÁLCULO DE BM
    mixture_repulsive_param = 0.0
    for i in range(number_molecular_components):
        mixture_repulsive_param = mixture_repulsive_param + \
            liquid_compositions[i]*repulsive_params[i]

    # CÁLCULO DE F
    sum_1 = 0.0
    g: List[float] = []
    for i in range(number_molecular_components):
        g[i] = (1.0/(c2-c1))*math.log((1+c2*qsi[i])/(1+c1*qsi[i]))/qsi[i]
        sum_1 = sum_1 + liquid_compositions[i]*(alpha[i]*g[i]*qsi[i]+math.log(
            (1.0-qsi[i])/qsi[i]) + math.log(repulsive_params[i]/mixture_repulsive_param))
    f1 = - delta_mixture_gibbs_energy + GAS_CONSTANT*temperature * \
        liquid_compositions[i]*math.log(activity_coefficients[i])

    # CÁLCULO DE QSIM
    hm = 0.0
    h: List[float] = []
    fmin = 0.0
    qsim = 0.0
    qsimin = 0.0
    for i in range(number_molecular_components):
        h[i] = 1.0/(1.0-qsi[i])-alpha[i]*qsi[i] / \
            ((1.0+c1*qsi[i])*(1.0+c2*qsi[i]))
        hm = hm + liquid_compositions[i]*h[i]

    if abs(hm) < 1.0e-5:
        fmin = 2.4744
        qsimin = math.sqrt(2.0)-1.0
    else:
        fmin = 2.47-1.25*hm
        qsimin = 0.414+0.015*hm

    if f1 <= fmin:
        aas: List[float] = []
        mixture_atractive_param: float = 0.0
        aij: List[List[float]] = [[]]
        for i in range(number_molecular_components):
            for j in range(number_molecular_components):
                aij[i][j] = math.sqrt(atractive_params[i]*atractive_params[j])
                mixture_atractive_param = mixture_atractive_param + \
                    liquid_compositions[i]*liquid_compositions[j]*aij[i][j]
                aas[i] = aas[i] + liquid_compositions[j]*aij[i][j]

        r1 = -c1
        r2 = -c2
        a10 = - (mixture_repulsive_param*(r1+r2+1.0) +
                 GAS_CONSTANT*temperature/pressure)
        a11 = (mixture_repulsive_param**2)*(r1*r2+r1+r2) + GAS_CONSTANT * \
            temperature*mixture_repulsive_param * \
            (r1+r2)/pressure + mixture_atractive_param/pressure
        a12 = - mixture_repulsive_param*(r1*r2*(mixture_repulsive_param**2) +
                                         r1*r2*mixture_repulsive_param*GAS_CONSTANT*temperature/pressure + mixture_atractive_param/pressure)
        volume: float = 0.0

        compressability_factor_list, nr = solve_pr_cubic_equation(
            a10, a11, a12)
        if nr <= 0:
            if system_phase == Phase.LIQUID:
                volume = compressability_factor_list[2]
            if system_phase == Phase.VAPOR:
                volume = compressability_factor_list[0]
        else:
            volume = compressability_factor_list[0]

        compressability_factor = volume/(GAS_CONSTANT*temperature/pressure)
        a2 = pressure*mixture_atractive_param/((GAS_CONSTANT*temperature)**2.0)
        y2 = pressure*mixture_repulsive_param/(GAS_CONSTANT*temperature)
        fugacities: List[float] = []
        fugacities_i: List[float] = []
        alnfi: List[float] = []
        for i in range(number_molecular_components):
            fugacities = liquid_compositions[i]*pressure*math.exp(repulsive_params[i]/mixture_repulsive_param)-math.log(compressability_factor-y2)-(a2/y2)*(1.0/(c1-c2))*(
                2.0*aas[i]/mixture_atractive_param - repulsive_params[i]/mixture_repulsive_param)*math.log((compressability_factor+c1*y2)/(compressability_factor+c2*y2))
            fugacities_i[i] = fugacities[i]/liquid_compositions[i]/pressure
            alnfi[i] = math.log(fugacities_i[i])
        return volume, fugacities

    qsim = 0.9999
    delta_f = 1.0
    ccc = 1
    while abs(delta_f) > THRESHOLD and ccc <= 100:
        gm = (1.0/(c2-c1))*math.log((1.0+c2*qsim)/(1.0+c1*qsim))/qsim
        f2 = (1.0+c1*qsim)*(1.0+c2*qsim)*gm * \
            (1.0/(1.0-qsim)-hm) + math.log((1.0-qsim)/qsim)
        delta_f = f2 - f1
        d3 = -hm*(1.0+(c1*c2*qsim**2.0-1.0)*gm)+gm*(-1.0+2.0 *
                                                    qsim+(c1+c2+c1*c2)*qsim**2.0)/((1.0-qsim)**2)
        qsim = qsim - delta_f/(d3/qsim)
        ccc = ccc + 1
    if qsim < qsimin or qsim > 1.0:
        print('ERRO 2 (HK)')
        print('qsim = ', qsim)
        print('QSImin = ', qsimin)
    # CÁLCULO DE ALFAM
    alpha_m = (1.0/(1.0-qsim)-hm)*(1.0+c1*qsim)*(1.0+c2*qsim)/qsim
    mixture_atractive_param = alpha_m*mixture_repulsive_param*GAS_CONSTANT*temperature

    # CÁLCULO DO COEFICIENTE DE FUGACIDADE
    r1 = - c1
    r2 = - c2
    a10 = - (mixture_repulsive_param*(r1+r2+1.0) +
             GAS_CONSTANT*temperature/pressure)
    a11 = (mixture_repulsive_param**2)*(r1*r2+r1+r2) + \
        GAS_CONSTANT*temperature*mixture_repulsive_param * \
        (r1+r2)/pressure + mixture_atractive_param/pressure
    a12 = - mixture_repulsive_param * \
        (r1*r2*(mixture_repulsive_param**2) + r1 *
         r2*mixture_repulsive_param*GAS_CONSTANT*temperature/pressure + mixture_atractive_param/pressure)
    compressability_factor_list, nr = solve_pr_cubic_equation(a10, a11, a12)
    volume = 0.0
    if nr <= 0:
        if system_phase == Phase.LIQUID:
            volume = compressability_factor_list[2]
        if system_phase == Phase.VAPOR:
            volume = compressability_factor_list[0]
    else:
        volume = compressability_factor_list[0]

    compressability_factor = volume/(GAS_CONSTANT*temperature/pressure)
    b1: List[float] = []
    b2: List[float] = []
    b3: List[float] = []
    b4: List[float] = []
    g1: List[float] = []
    g2: List[float] = []
    g3: List[float] = []
    g4: List[float] = []
    g5: List[float] = []
    d4: List[float] = []
    d5: List[float] = []
    h10: List[float] = []
    fugacities: List[float] = []
    fugacities_i: List[float] = []
    alnfi: List[float] = []
    for i in range(number_molecular_components):
        # Cálculo baseado em parâmetros de componente puro e BM -----
        g1[i] = -math.log(activity_coefficients[i])
        g2[i] = qsi[i]*alpha[i]*g[i]
        g3[i] = math.log((1.0-qsi[i])/qsi[i])
        g4[i] = math.log(repulsive_params[i]/mixture_repulsive_param)
        g5[i] = -(repulsive_params[i]-mixture_repulsive_param) / \
            mixture_repulsive_param
        d4[i] = g1[i] + g2[i] + g3[i] + g4[i] + g5[i]
        # Cálculo baseado no estado de referência P=0 ---------------
        h1 = (c1 + c2 + c1*c2)*qsim*qsim + 2.0*qsim - 1.0
        h2 = (1.0-qsim)**2.0
        h3 = (1.0+c2*qsim)/(1.0+c1*qsim)
        h4 = (1.0+c2*qsim)*(1.0+c1*qsim)
        h5 = math.log(h3)/(c2-c1)/qsim
        h6 = h5*h1/h2 - ((c1*c2*qsim*qsim-1.0)*h5+1.0)*hm
        h7 = hm/h6
        h8 = 1.0/(1.0-qsim)
        h9 = math.log((1.0-qsim)/qsim)
        h10[i] = ((d4[i]-h9)*(1.0+h7)+h5*h4*(h[i]-h8)*h7)/(qsim*h5)
        # Cálculo baseado nas condições de mistura ------------------
        h11 = mixture_repulsive_param / volume
        h12 = (1.0+c2*h11)/(1.0+c1*h11)
        h13 = math.log(h12)/(c2-c1)
        d5[i] = h13*h10[i]
        b1[i] = math.log(volume/(volume-mixture_repulsive_param))
        b2[i] = mixture_repulsive_param[i]/(volume-mixture_repulsive_param)
        b3[i] = -alpha_m*mixture_repulsive_param[i]*volume / \
            ((volume+c1*mixture_repulsive_param)
             * (volume+c2*mixture_repulsive_param))
        b4[i] = -math.log(compressability_factor)
        alnfi[i] = b1[i] + b2[i] + b3[i] + b4[i] - d5[i]
        fugacities_i[i] = math.exp(alnfi[i])
        fugacities[i] = liquid_compositions[i]*pressure*fugacities_i[i]
    return volume, fugacities


def calculate_activity_coefficient_fugchk(system_phase: Phase, number_components: int, temperature: float, pressure: float, compositions: List[float], atractive_params: List[float], repulsive_params: List[float], uniquac_params: UniquacProperties):
    """Calculate activity coefficient with fugchk

    Args:
        system_phase (Phase): _description_
        number_components (int): _description_
        temperature (float): _description_
        pressure (float): _description_
        compositions (List[float]): _description_
        atractive_params (List[float]): _description_
        repulsive_params (List[float]): _description_
        uniquac_params (UniquacProperties): _description_

    Returns:
       volume, fugacities (tuple(float, list[float])): _description_
    """
    # ********************************************************************
    #
    #                   PROGRAMA FUGCHK
    #
    #   CALCULO DO COEFICIENTE DE FUGACIDADE
    #
    #    VARIAVEIS:
    #       IFASE    = ESCOLHA DO TIPO DE FASE (1=VAP, 2=LIQ)
    #       NC       = NUM. DE ESPEC.
    #       T        = TEMP.
    #       P        = PRES.
    #       COMP(NC) = COMPOSICAO
    #       FUG(NC)  = COEF. DE FUGACIDADE
    #
    # ********************************************************************
    itdm = 1

    for _ in range(15):
        itdm = itdm + number_components

    alff, activity_coefficients, mixture_atractive_param, mixture_repulsive_param, mix_params = calculate_cubic_params_fugacity_terms(
        number_components, temperature, compositions, atractive_params, repulsive_params, uniquac_params)

    compressability_factor = calculate_compressability_factor(
        system_phase, temperature, pressure, mixture_atractive_param, mixture_repulsive_param)

    fugacities = calculate_coefficient_fugacity(number_components, temperature, pressure, repulsive_params,
                                                mixture_repulsive_param, alff, activity_coefficients, compressability_factor, mix_params)

    volume = compressability_factor*GAS_CONSTANT_SI*temperature/pressure
    return volume, fugacities


def calculate_cubic_params_fugacity_terms(number_components: int, temperature: float, compositions: List[float], atractive_params: List[float], repulsive_params: List[float], uniquac_params: UniquacProperties):
    """calculate cubic params fugacity terms

    Args:
        number_components (int): _description_
        temperature (float): _description_
        compositions (List[float]): _description_
        atractive_params (List[float]): _description_
        repulsive_params (List[float]): _description_
        uniquac_params (UniquacProperties): _description_

    Return
        alff, activity_coefficients, mixture_atractive_param, mixture_repulsive_param, mix_params: tuple(float, List[float], float, float,  MixParams)

    """
    # --------------------------------------------------------------------
    #
    #     SUBROTINA: MIX
    #
    #     FINALIDADE: Calcular os coeficientes da equação cúbica em Z
    #                 e partes que vão interessar ao cálculo do
    #                 coeficiente de fugacidade.
    #
    #     REFERÊNCIAS: Boukouvalas et al., FPE, 92, (1994), 75-106
    #
    #     AUTORA: G.M.N. Costa
    #
    #     MODIFICAÇÕES: S.A.B. Vieira de Melo (MAI.97)
    #                   F.L.P.PESSOA (NOV.97)
    #
    # ----------------------------------------------------------------------
    # Cálculo de B da mistura
    mixture_repulsive_param = 0.0
    for i in range(number_components):
        mixture_repulsive_param = mixture_repulsive_param + \
            compositions[i]*repulsive_params[i]

    alf1 = 0.0
    bblog = 0.0
    alff: List[float] = []
    for i in range(number_components):
        #       Cálculo do somatório de X(I) e ALFA(I)
        alff[i] = atractive_params[i] / \
            repulsive_params[i]/GAS_CONSTANT_SI/temperature
        alf1 = alf1+compositions[i]*alff[i]
#       CáLculo do somatório de X[i]LNB[i]
        bblog = bblog + \
            compositions[i] * \
            math.log(mixture_repulsive_param/repulsive_params[i])
    # Valores dos parâmetros LAMBDA,AV,AM da EOS e da regra de mistura
    av = -0.623
    amm = -0.52
    clamb = 0.36

#    Cálculo do valor de ALFA dado pela equação (5) do artigo do
#    BOUKOUVALAS ET AL.(1994).
    alf2 = (1.0-clamb)/amm
    alf3 = (clamb/av)+alf2
    alf4 = bblog*alf2

    # Chamada do UNIQUAC
    activity_coefficients = calculate_gibbs_excess_activity_coefficients(
        number_components, temperature, compositions, uniquac_params)

    gert = 0.0
    for i in range(number_components):
        gert = gert + compositions[i]*math.log(activity_coefficients[i])

    # GE=GE/RT e GAM é o LN GAM(I)

    alf = alf3*gert+alf1
    alf = alf+alf4

    mix_params = MixParams(av, amm, clamb, alf)

    # Cálculo de A da mistura
    mixture_atractive_param = alf*mixture_repulsive_param*GAS_CONSTANT_SI*temperature

    return alff, activity_coefficients, mixture_atractive_param, mixture_repulsive_param, mix_params


def calculate_compressability_factor(system_phase: Phase, temperature: float, pressure: float, mixture_atractive_param: float, mixture_repulsive_param: float):
    """calculate compressebality factor

    Args:
        system_phase (Phase): _description_
        temperature (float): _description_
        pressure (float): _description_
        mixture_atractive_param (float): _description_
        mixture_repulsive_param (float): _description_

    Returns:
        compressability_factor: float
    """
    # *******************************************************************
    #
    #                   PROGRAMA FATORZ
    #
    #   CALCULO DO FATOR DE COMP.
    #
    #    VARIAVEIS:
    #       ITYP  = TIPO DE FASE (-1=VAP  1=LIQ)
    #       P     = PRESSAO      .
    #       T     = TEMPERATURA
    #       BM    = PAR. B DE MIST.
    #       AM    = PARAM. A DE MIST.
    #       Z     = FATOR DE COMP.
    #
    # *******************************************************************
    r1 = - (1.0 + math.sqrt(2.0))
    r2 = - (1.0 - math.sqrt(2.0))
    a10 = - (mixture_repulsive_param*(r1+r2+1.0) +
             GAS_CONSTANT_SI*temperature/pressure)
    a11 = (mixture_repulsive_param**2)*(r1*r2+r1+r2) + \
        GAS_CONSTANT_SI*temperature*mixture_repulsive_param * \
        (r1+r2)/pressure + mixture_atractive_param/pressure
    a12 = - mixture_repulsive_param * \
        (r1*r2*(mixture_repulsive_param**2) + r1 *
         r2*mixture_repulsive_param*GAS_CONSTANT_SI*temperature/pressure + mixture_atractive_param/pressure)

    compressability_factor_list, nr = solve_pr_cubic_equation(a10, a11, a12)

    volume = 0.0
    if nr <= 0:
        if system_phase == Phase.LIQUID:
            volume = compressability_factor_list[2]
        if system_phase == Phase.VAPOR:
            volume = compressability_factor_list[0]
    else:
        volume = compressability_factor_list[0]

    compressability_factor = volume/(GAS_CONSTANT_SI*temperature/pressure)

    return compressability_factor


def calculate_coefficient_fugacity(number_components: int, temperature: float, pressure: float, repulsive_params: List[float], mixture_repulsive_param: float, alff: List[float], activity_coefficients: List[float], compressability_factor: float, mix_params: MixParams):
    """calculate coefficient fugacity

    Args:
        number_components (int): _description_
        temperature (float): _description_
        pressure (float): _description_
        repulsive_params (List[float]): _description_
        mixture_repulsive_param (float): _description_
        alff (List[float]): _description_
        activity_coefficients (List[float]): _description_
        compressability_factor (float): _description_
        mix_params (MixParams): _description_

    Returns:
        fugacities: List[float]
    """
    # *******************************************************************
    # *                                                                 *
    # *   TITULO:                                                       *
    # *     SUBROTINA FUGCOF
    # *     BOUKOUVALAS ET AL. (1994).
    # *                                                                 *
    # *******************************************************************

    d1 = 1.0-math.sqrt(2.0)
    d2 = 1.0+math.sqrt(2.0)
    prt = pressure/GAS_CONSTANT_SI/temperature
    bmix = mixture_repulsive_param*prt
    t1 = compressability_factor+bmix*d2
    t2 = compressability_factor+bmix*d1
    ta = -math.log(t1/t2)
    d3 = 2.0*math.sqrt(2.0)
    tta = ta/d3

    ttb = -math.log(compressability_factor-bmix)

    ttc = (compressability_factor-1.0)/mixture_repulsive_param

    alf2 = (1.0-mix_params.clamb)/mix_params.amm
    alf3 = (mix_params.clamb/mix_params.av)+alf2

    fugacities: List[float] = []
    for i in range(number_components):
        alf5 = repulsive_params[i]/mixture_repulsive_param
        alf6 = -math.log(alf5)+alf5-1.0
        alf7 = alf2*alf6+alff[i]
        add = alf3*math.log(activity_coefficients[i])+alf7
        fugacities[i] = math.exp(ttc*repulsive_params[i]+ttb+add*tta)

    return fugacities
