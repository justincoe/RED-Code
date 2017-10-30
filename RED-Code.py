import math
import scipy.optimize as op
import numpy


class tube:
    def __init__(self, name, section, material, cost, l, d, thickness, e, style, kb, ba, br):
        # type: (str, str, str, float, float, float, float, float, str, float, float, float) -> tube
        self.name = name
        self.section = section
        self.material = material
        self.cost = cost
        self.l = l
        self.d = d
        self.thickness = thickness
        self.e = e
        self.style = style
        self.kb = kb
        self.ba = ba
        self.br = br


class valve:
    def __init__(self, name, material, weight, cost, d, style, k):
        self.name = name
        self.material = material
        self.weight = weight
        self.cost = cost
        self.d = d
        self.style = style
        self.k = k


def velocity(m_dot, current_d, rho):
    # Function computes the velocity for changing dimensions
    # m_dot is in slug / in ^ 3, d is in inches, rho is in slug / in ^ 3
    # output velocity is in [ in / s]
    a = (math.pi / 4) * (math.pow(current_d, 2))
    v = m_dot / (rho * a)
    return v


def reynolds(v, d, rho, mu):
    # Function computes reynolds number when called
    # input: velocity, diameter, density, dynamic viscosity
    # density and dynamic viscosity could be changing...
    re = rho * v * d / mu
    return re


def relative_roughness(e, d):
    # sort function to calculate the relative roughness inside a tube
    # value is used to calculate the friction factor
    ed = e / d
    return ed


def moody(ed, re):
    # moody Find friction factor b solving the Colebrook equation(Moody Chart)
    #
    # Synopsis: f = moody(ed, Re)
    #
    # Input: ed = relative roughness = epsilon / diameter
    # Re = Reynolds number
    # ]
    # Output: f = friction factor
    #
    # Note: Laminar and turbulent flow are correctly accounted for

    if re < 0:
        print 'Reynolds number = ' + re + ' cannot be negative'
        return 0

    if re < 2000:
        f = 64 / re
        return f
    # laminar flow end
    if ed > 0.05:
        print 'epsilon/diameter ratio = ' + ' is not on Moody chart'
    if re < 4000:
        print 'Re = ' + re + ' in transition range'

    # --- Use fzero to find f from Colebrook equation.
    # coleFun is an inline function object to evaluateF(f, e / d, Re)
    # fzero returns the value of f such that F(f, e / d / Re) = 0(approximately)
    # fi = initial guess from Haaland equation, see White, equation 6.64 a
    # Iterations of fzero are terminated when f is known to whithin + / - dfTol

    def cole_fun(f_in):
        1.0 / math.sqrt(f_in) + 2.0 * math.log10(ed / 3.7 + 2.51 / (re * math.sqrt(f)))

    fi = numpy.empty()
    fi[1] = [1 / (1.8 * math.pow(math.log10(6.9 / re + math.pow((ed / 3.7), 1.11)), 2))]
    #  initial guess at f
    # dfTol = 5e-6; --- TOOK OUT ---
    f = op.fsolve(cole_fun, fi)
    # --- sanity check:

    if f < 0:
        print 'Friction factor = ', f, ' , but cannot be negative'
        return f
    return f


def tube_dp(l, d, v, rho, f, style, br, ba, kb):
    # type: (float, float, float, float, float, str, float, float, float, float) -> float
    # re removed because not used
    # Function computes tube pressure drop
    # To calculate pressure accross the tube
    # 1. calculate the Reynolds number(Re function)
    # 2. determine the roughness
    # 3. determine the frictionfactor(moodyfunction)
    # 4. calculate the head loss
    # 5. convert to pressure drop

    if style == 'straight':
        head = f * l / d * (math.pow(v, 2)) / 2
        pdrop = head * rho
        return pdrop
    elif style == 'curved':
        pdrop = (0.5 * f * rho * (math.pow(v, 2)) * math.pi * br / d * ba / 180) + (0.5 * kb * rho * (math.pow(v, 2)))
        return pdrop


def valvedp(k, v, rho):
    # function computes pressure drop across different valves
    # things that are needed, K factor, dimensions, flow characteristics, etc.
    # depending on type, pressure drop is calculated accordingly
    # The velocity is the average velocity
    head = k * (math.pow(v, 2)) / 2
    pdrop = head * rho
    return pdrop


def flowcoef(q, sg, param):
    # type: (float, float, float) -> float
    # Short function computes the flow coefficient
    # Function uses the volumetric flow rate, specific gravity, and pressure
    # drop through the element...
    cv = q * math.sqrt(sg / param)
    return cv


def pressure2gg(pi, m_dot, rho, mu, sg, q, t_tube, v_valve):
    # Function computes the pressure to the RP1 gas generator inlet
    # Detailed explanation goes here
    dp = [0, 0]
    cv = [0, 0]
    # Flow paths: RP1 pump outlet to gas generator inlet:
    # tube section 1:
    v = velocity(m_dot, t_tube.d, rho)
    re = reynolds(v, t_tube.d, rho, mu)
    ed = relative_roughness(t_tube.e, t_tube.d)
    f = moody(ed, re)
    dp[0] = tube_dp(t_tube.l, t_tube.d, v, rho, f, t_tube.style, t_tube.br, t_tube.ba, t_tube.kb)
    # Cv(1) = flowcoef(Q,SG,dp(1));
    po = pi - dp[0]
    # through check valve 2:
    p = po
    v = velocity(m_dot, v_valve.d, rho)
    dp[1] = valvedp(v_valve.k, v, rho)
    cv[1] = flowcoef(q, sg, dp[1])
    p_gg = p - dp[1]
    print "The result of this function is p_gg =", p_gg
    # tube section 2: DEALING WITH DIVERGING SECTIONS>....
    # outlet of pump
    # passes through some tube
    # passes a tee with a poppet valve a the end
    # pass through a dividing section
    # passes a check valve
    # passes a tee section leading to another path
    # throught tube
    # passes through poppet valve
    # tube
    # into GG.
    return


check_valve2 = valve('FIT6-23', 'Stainless_Steel', 123, 1200, 3.5, 'Check', 2)
rp1_tube1 = tube('name1', 'section1', 'Stainless-Steel', 1500, 12*2.5, 3.5, 1, 12*0.000049, 'straight', 0, 0, 0)

pressure2gg(100, 100, 10, .5, 200, 100, rp1_tube1, check_valve2)  # Test case numbers are Sam's (Dumb) guesses
