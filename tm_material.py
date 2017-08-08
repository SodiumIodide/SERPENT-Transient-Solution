'''
Import only module. Material definition for transient model.

Dependencies:
CoolProp for material properties
numpy for mathematical constants (should be included as part of main file import already)
'''

# External
from CoolProp.CoolProp import PropsSI  # Assuming that the solution is approximated by water
from numpy import pi  # Bring pi into namespace for cylindrical volume calculations

# Shared
import tm_constants as c  # Required for the included tag dictionaries

class Material():
    '''Generic class for the material in question, and associated geometry'''
    def __init__(self, matnum, elems, ndens, dens, height, base_height, radius,
                 inner_radius, pres=101325, temp=300):
        self.matnum = matnum
        self.elems = elems
        self.ndens = ndens  # a/b-cm
        self.temp = temp  # K
        self.height = height  # cm
        self.base_height = base_height  # cm
        self.radius = radius  # cm
        self.inner_radius = inner_radius  # cm
        self.dens = dens  # g/cm^3
        self.mass = 0.0  # g, placeholder until calculated
        self.volume = 0.0  # cm^3, placeholder until calculated
        self.base = 0.0  # cm^2, placeholder until calculated
        self.kappa = 0.0  # 1/Pa, isothermal compressibility, placeholder until needed
        self.beta = 0.0  # 1/K, isobaric compressibility, placeholder until needed
        self.delta_temp = 0.0  # K, placeholder until required for expansion
        self.vol_accel = 0.0  # cm^3/s^2, initially zero upon definition
        self.vol_vel = 0.0  # cm^3/s, initially zero upon definition
        self.delta_vol = 0.0  # cm^3, initially zero upon definition
        self.delta_temp = 0.0  # K, initially zero upon definition
        self.delta_pres = 0.0  # Pa, initially zero upon definition
        self.atoms = [0.0] * len(ndens)  # _, placeholder until calculated
        self.volfrac_gas = 0.0  # _, initially no radiolytic gas in solution
        self.mass_h2 = 0.0  # g, initially no radiolytic gas in solution
        # Pressure: Assume atmospheric conditions (allow for expansion of fluid)
        self.av_pressure = pres  # Pa
        self.xs_tag = "03c"
        self.sab_tag = "00t"
        if temp != 300:  # K
            self.__update_xs_tag()
            self.__update_sab_tag()
        self.__calc_init()

    # Setter method used here for clarity
    def set_vol_accel(self, newva):
        '''Adjust volume acceleration'''
        self.vol_accel = newva  # cm^3/s^2s

    def update_height(self, delta_t):
        '''Calculate new height based on volume acceleration and time'''
        # Reliant on Newtonian kinematics equations
        vol_vel_i = self.vol_vel  # cm^3/s
        vol_vel_f = vol_vel_i + self.vol_accel * delta_t  # cm^3/s
        self.vol_vel = vol_vel_f  # cm^3/s, update information for next snapshot
        self.delta_vol = vol_vel_i * delta_t + 1 / 2 * self.vol_accel * delta_t**2  # cm^3
        height_diff = self.delta_vol / self.base  # cm
        self.append_height(height_diff, add=True)  # Number density bookkeeping

    def append_height(self, newheight, add=False):
        '''Append a new volume after adjusting height'''
        self.height = self.height + newheight if add else newheight
        self.volume = self.base * self.height  # cm
        self.ndens = [atom * 1e-24 / self.volume for atom in self.atoms]  # a/b-cm

    def old_calc_temp(self, fissions):
        '''Adjust material temperature'''
        energydep = fissions * 180 * 1.6022e-13  # J
        # Based on intensive properties (specific heat) of water
        # This is an aqueous solution of material
        spec_heat = PropsSI('C', 'T', self.temp, 'P', self.av_pressure, 'WATER')  # J/kg-K
        newtemp = self.temp + energydep / (self.mass / 1e3) / spec_heat  # K
        self.delta_temp = newtemp - self.temp  # K
        self.temp = newtemp  # K
        self.__update_xs_tag()
        self.__update_sab_tag()

    def __update_xs_tag(self):
        '''Update the materials cross-section tag; called when new temperature calculated'''
        round_value = 300  # technically K
        rounded_temp = int(round_value * round(self.temp / round_value))
        self.xs_tag = c.XS_TEMP.get(rounded_temp, "18c")

    def __update_sab_tag(self):
        '''Update the thermal scattering correction tag; called when new temperature calculated'''
        round_value = 50  # technically K
        rounded_temp = int(round_value * round(self.temp / round_value))
        self.sab_tag = c.SAB_TEMP.get(rounded_temp, "18t")

    def __produce_gas(self, fissions):
        '''Produce a mass of radiolytic gas in the solution; called from self.update_pressure()'''
        energydep = fissions * 180 * 1.6022e-19  # MJ
        self.mass_h2 = c.RADIOLYTIC_G * energydep * 1000  # g
        self.__update_volfrac_gas()

    def __update_volfrac_gas(self):
        '''Update the volume fraction of gas, f_e; called from self.produce_gas()'''
        self.volfrac_gas = self.mass_h2 / 1000 * c.RH2 * self.temp \
                           / (self.volume / 100**3 * (self.av_pressure + 2 * 
                                             PropsSI('I', 'T', self.temp, 'Q', 0.0, 'WATER')
                                             / c.RAD_GAS_BUBBLE))

    def expand_old(self):
        '''Self-regulate the expansion of the material'''
        # The isobaric expansion coefficient of water is calculated
        expansion = PropsSI('ISOBARIC_EXPANSION_COEFFICIENT', 'T', self.temp,
                            'P', self.av_pressure, 'WATER')  # 1/K
        # Coefficient is defined as: alpha = 1/V * (dV/dT) over constant P
        # Where V is volume, T is temperature, P is pressure, and alpha is the coefficient
        # In this case, changing differential into delta, and solving for a delta height
        # while maintaining base area
        delta_height = self.delta_temp * self.volume * expansion / self.base  # cm
        newheight = self.height + delta_height  # cm
        self.append_height(newheight)  # Book-keeping for number densities and such

    def update_state(self, fissions):
        '''Self-regulate the expansion state of the material'''
        # Material constants are assumed to closely match those of water for an
        # aqueous solution
        beta_0 = PropsSI('ISOBARIC_EXPANSION_COEFFICIENT', 'T', self.temp,
                         'Q', 0.0, 'WATER')  # 1/K
        kappa_0 = PropsSI('ISOTHERMAL_EXPANSION_COEFFICIENT', 'T', self.temp,
                          'Q', 0.0, 'WATER')  # 1/Pa
        surf_tens = PropsSI('I', 'T', self.temp, 'Q', 0.0, 'WATER')  # N/m
        self.__produce_gas(fissions)
        self.beta = beta_0 * (1 - self.volfrac_gas) + self.volfrac_gas / self.temp \
                    * (self.av_pressure + 2 * surf_tens / c.RAD_GAS_BUBBLE) \
                    / (self.av_pressure + 4 * surf_tens / 3 / c.RAD_GAS_BUBBLE)
        self.kappa = kappa_0 * (1 - self.volfrac_gas) + self.volfrac_gas \
                     / (self.av_pressure + 4 * surf_tens / 3 / c.RAD_GAS_BUBBLE)
        # Constant volume specific heat
        spec_heat = PropsSI('O', 'T', self.temp, 'P', self.av_pressure, 'WATER')  # J/kg-K
        self.delta_temp = 1 / spec_heat * (fissions * 180 * 1.6022e-13 - self.beta \
                          / self.kappa * self.delta_vol) / self.mass  # K
        self.temp += self.delta_temp  # K
        self.__update_xs_tag()
        self.__update_sab_tag()
        self.delta_pres = self.beta / self.kappa * self.delta_temp - 1 \
                          / (self.kappa * self.volume / 100**3) * self.delta_vol  # Pa
        self.av_pressure += self.delta_pres  # Pa

    def __calc_init(self):
        '''Self-called method to calculate some constants after initial file run'''
        # Ideally only called ONE time after each material definition
        self.base = pi * (self.radius**2 - self.inner_radius**2)  # cm^2
        self.volume = self.base * self.height  # cm^3
        self.atoms = [nden * 1e24 * self.volume for nden in self.ndens]  # a
        self.mass = self.dens * self.volume  # g

    def __str__(self):
        # Material representation, overload built-in string definition
        ret = "mat solution{0} sum moder lwtr{0} 1001 tmp {1}\n".format(self.matnum, self.temp)
        template = "{0}.{1} {2}\n"
        for ind, elem in enumerate(self.elems):
            ret += template.format(elem, self.xs_tag, self.ndens[ind])
        return ret + "\n"  # str
