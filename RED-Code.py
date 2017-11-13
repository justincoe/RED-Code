import math
import scipy.optimize as op
import numpy


class tube:
    def __init__(self, name, section, material, cost, length, diameter, thickness, roughness_e, style, kb, ba, br):
        # type: (str, str, str, float, float, float, float, float, str, float, float, float) -> tube
        # English Units
        self.name = name
        self.section = section
        self.material = material
        self.cost = cost
        self.length = length
        self.diameter = diameter
        self.thickness = thickness
        self.roughness_e = roughness_e
        self.style = style
        self.kb = kb
        self.ba = ba
        self.br = br


class valve:
    def __init__(self, name, material, weight, cost, diameter, style, k):
        self.name = name
        self.material = material
        self.weight = weight
        self.cost = cost
        self.diameter = diameter
        self.style = style
        self.k = k


def velocity_calc(m_dot, current_d, rho):
    # Function computes the velocity for changing dimensions
    # m_dot is in slug / s, d is in inches, rho is in slug / in ^ 3
    # output velocity is in [ in / s]
    a = (math.pi / 4) * (math.pow(current_d, 2))
    v = m_dot / (rho * a)
    return v


def reynolds(velocity, diameter, rho, dyn_visc_mu):
    # Function computes reynolds number when called
    # input: velocity, diameter, density, dynamic viscosity
    # density and dynamic viscosity could be changing...
    reynolds_num = rho * velocity * diameter / dyn_visc_mu
    return reynolds_num


def relative_roughness_calc(roughness_e, diameter):
    # sort function to calculate the relative roughness inside a tube
    # value is used to calculate the friction factor
    realitve_roughness = roughness_e / diameter
    return realitve_roughness


def moody(relative_roughness, reynolds_num):
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

    if reynolds_num < 0:
        print('Reynolds number = ' + reynolds_num + ' cannot be negative')
        return 0

    if reynolds_num < 2000:
        friction_factor = 64 / reynolds_num
        return friction_factor
    # laminar flow end
    if relative_roughness > 0.05:
        print('epsilon/diameter ratio = ' + ' is not on Moody chart')
    if reynolds_num < 4000:
        print('Re = ' + reynolds_num + ' in transition range')

    # --- Use fzero to find f from Colebrook equation.
    # coleFun is an inline function object to evaluateF(f, e / d, Re)
    # fzero returns the value of f such that F(f, e / d / Re) = 0(approximately)
    # fi = initial guess from Haaland equation, see White, equation 6.64 a
    # Iterations of fzero are terminated when f is known to whithin + / - dfTol

    def cole_fun(f_in):
        1.0 / math.sqrt(f_in) + 2.0 * math.log10(relative_roughness / 3.7 + 2.51 / (relative_roughness * math.sqrt(friction_factor)))

    fi = numpy.empty([1,2])
    fi[1] = [1 / (1.8 * math.pow(math.log10(6.9 / reynolds_num + math.pow((relative_roughness / 3.7), 1.11)), 2))]
    #  initial guess at f
    # dfTol = 5e-6; --- TOOK OUT ---
    friction_factor = op.fsolve(cole_fun, fi)
    # --- sanity check:

    if friction_factor < 0:
        print('Friction factor = ', friction_factor, ' , but cannot be negative')
        return friction_factor
    return friction_factor


def tube_dp(length, diameter, velocity, rho, friction_factor, style, br, ba, kb):
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
        head = friction_factor * length / diameter * (math.pow(velocity, 2)) / 2
        pdrop = head * rho
        return pdrop
    elif style == 'curved':
        pdrop = (0.5 * friction_factor * rho * (math.pow(velocity, 2)) * math.pi * br / diameter * ba / 180) + (0.5 * kb * rho * (math.pow(velocity, 2)))
        return pdrop


def valvedp(k, velocity, rho):
    # function computes pressure drop across different valves
    # things that are needed, K factor, dimensions, flow characteristics, etc.
    # depending on type, pressure drop is calculated accordingly
    # The velocity is the average velocity
    head = k * (math.pow(velocity, 2)) / 2
    pdrop = head * rho
    return pdrop


def flowcoef(vol_flow, spec_grav, pressure_diff):
    # type: (float, float, float) -> float
    # Short function computes the flow coefficient
    # Function uses the volumetric flow rate, specific gravity, and pressure
    # drop through the element...
    cv = vol_flow * math.sqrt(spec_grav / pressure_diff)
    return cv


def pressure_section(pressure_in, m_dot, rho, dyn_visc_mu, spec_grav, vol_flow, t_tube, v_valve, section_num):
    # Pressure psi, m_dot slugs/s, rho slugs/in^3, dyn_visc_mu slugs/(in*s), vol_flow gallon/min
    # Function computes the pressure to the RP1 gas generator inlet
    # Detailed explanation goes here
    dp = [0, 0]
    cv = [0, 0]
    # Flow paths: RP1 pump outlet to gas generator inlet:
    # tube section 1:
    velocity = velocity_calc(m_dot, t_tube.diameter, rho)
    reynolds_num = reynolds(velocity, t_tube.diameter, rho, dyn_visc_mu)
    relative_roughness = relative_roughness_calc(t_tube.roughness_e, t_tube.diameter)
    friction_factor = moody(relative_roughness, reynolds_num)
    dp[0] = tube_dp(t_tube.length, t_tube.diameter, velocity, rho, friction_factor, t_tube.style, t_tube.br, t_tube.ba, t_tube.kb)
    # Cv(1) = flowcoef(Q,SG,dp(1));
    pressure_out = pressure_in - dp[0]
    # through check valve 2:
    pressure = pressure_out
    velocition = velocity_calc(m_dot, v_valve.diameter, rho)
    dp[1] = valvedp(v_valve.k, velocity, rho)
    cv[1] = flowcoef(vol_flow, spec_grav, dp[1])
    pressure_end = pressure - dp[1]
    print("The pressure at the end of section ", section_num, " is ", pressure_end, " psi")
    # tube section 2: DEALING WITH DIVERGING SECTIONS>....
    # outlet of pump
    # passes through some tube
    # passes a tee with a poppet valve a the end
    # pass through a dividing section
    # passes a check valve
    # passes a tee section leading to another path
    # throughout tube
    # passes through poppet valve
    # tube
    # into GG.
    return pressure_end


check_valve2 = valve('FIT6-23', 'Stainless_Steel', 123, 1200, 3.5, 'Check', 2)
rp1_tube1 = tube('name1', 'section1', 'Stainless-Steel', 1500, 12*2.5, 3.5, 1, 12*0.000049, 'straight', 0, 0, 0)
RHO_RP1 = 0.9095    # Slugs / in^3
# MU is in slug / in seconds
pressure_in_0 = 100
pressure_in_1 = pressure_section(pressure_in_0, 3.85, RHO_RP1, .0141, RHO_RP1, 1100, rp1_tube1, check_valve2, 1)  # Test case numbers
pressure_in_2 = pressure_section(pressure_in_1, 3.85, RHO_RP1, .0141, RHO_RP1, 1100, rp1_tube1, check_valve2, 2)
pressure_in_3 = pressure_section(pressure_in_2, 3.85, RHO_RP1, .0141, RHO_RP1, 1100, rp1_tube1, check_valve2, 3)