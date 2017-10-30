# -*- coding: utf-8 -*-
# tm_material.py

'''
Import only module. Material definition for transient model, defined as a discretized
segment of overall solution volume.

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
                 inner_radius, com_height, pres=c.ATM, temp=300):
        self.matnum = matnum
        self.elems = elems
        self.ndens = ndens  # a/b-cm
        self.temp = temp  # K
        self.height = height  # cm, note that this is the ABSOLUTE HEIGHT of the material from 0
        self.base_height = base_height  # cm
        self.radius = radius  # cm
        self.inner_radius = inner_radius  # cm
        self.dens = dens  # g/cm^3
        self.mass = 0.0  # g, placeholder until calculated
        self.volume = 0.0  # cm^3, placeholder until calculated
        self.base = 0.0  # cm^2, placeholder until calculated
        self.kappa = 0.0  # 1/Pa, isothermal compressibility, placeholder until needed
        self.beta = 0.0  # 1/K, isobaric compressibility, placeholder until needed
        self.delta_temp = 0.0  # K, initially zero upon definition
        self.delta_pres = 0.0  # Pa, initially zero upon definition
        self.atoms = [0.0] * len(ndens)  # _, placeholder until calculated
        self.volfrac_gas = 0.0  # _, initially no radiolytic gas in solution
        self.mass_h2 = 0.0  # g, initially no radiolytic gas in solution
        # Pressure: Assume atmospheric conditions (allow for expansion of fluid)
        self.av_pressure = pres  # MPa, gauge value
        self.__calc_bottom_pressure()  # MPa, calculate from average pressure input
        self.com_height = com_height  # cm
        self.com_accel = 0.0  # cm/s^2, initially not in motion
        self.com_vel = 0.0  # cm/s, initially not in motion
        self.delta_com = 0.0  # cm, placeholder, no initial delta
        self.delta_vol = 0.0  # cm^3, placeholder, no initial delta
        self.xs_tag = "03c"
        self.sab_tag = "00t"
        self.gas_production_flag = False  # Once radiolytic gas is produced, keep producing it
        if temp != 300:  # K
            self.__update_xs_tag()
            self.__update_sab_tag()
        self.__calc_init()

    # Setter method used here for clarity
    def set_com_accel(self, newcoma):
        '''Adjust center of mass acceleration'''
        self.com_accel = newcoma  # cm/s^2
        self.__update_com_height()
        self.update_mat_height()

    def update_mat_height(self):
        '''Change the material height'''
        half_height = self.com_height - self.base_height  # cm
        new_height = self.base_height + 2 * half_height  # cm
        self.append_height(new_height)

    def height_shift(self, shift):
        '''Change the positioning of each material layer'''
        self.base_height += shift  # cm
        self.height += shift  # cm

    # Public function for refactoring purposes
    def append_height(self, newheight, add=False):
        '''Append a new volume after adjusting height'''
        oldvol = self.volume  # cm^3
        self.height = self.height + newheight if add else newheight  # cm
        self.volume = self.base * (self.height - self.base_height)  # cm^3
        self.delta_vol = self.volume - oldvol  # cm^3
        self.dens = self.mass / self.volume  # g/cm^3
        self.ndens = [atom * 1e-24 / self.volume for atom in self.atoms]  # a/b-cm

    def update_state(self, fissions, top_pres, initial=False):
        '''Self-regulate the expansion state of the material'''
        # Material constants are assumed to closely match those of water for an
        # aqueous solution
        beta_0 = PropsSI('ISOBARIC_EXPANSION_COEFFICIENT', 'T', self.temp,
                         'Q', 0.0, 'WATER')  # 1/K
        kappa_0 = PropsSI('ISOTHERMAL_COMPRESSIBILITY', 'T', self.temp,
                          'Q', 0.0, 'WATER') * 1e6  # 1/MPa
        surf_tens = PropsSI('I', 'T', self.temp, 'Q', 0.0, 'WATER') * 1e-6  # MN/m
        if fissions / self.volume / 0.001 > c.THRESHOLD or self.gas_production_flag:  # fis/liter
            self.__produce_gas(fissions)
            if not self.gas_production_flag:
                self.gas_production_flag = True
        self.beta = beta_0 * (1 - self.volfrac_gas) + self.volfrac_gas / self.temp \
                    * (self.av_pressure + 2 * surf_tens / c.RAD_GAS_BUBBLE) \
                    / (self.av_pressure + 4 * surf_tens / 3 / c.RAD_GAS_BUBBLE)  # 1/K
        self.kappa = kappa_0 * (1 - self.volfrac_gas) + self.volfrac_gas \
                     / (self.av_pressure + 4 * surf_tens / 3 / c.RAD_GAS_BUBBLE)  # 1/MPa
        # Constant volume specific heat
        spec_heat = PropsSI('O', 'T', self.temp, 'P', self.av_pressure + c.ATM_ABS, 'WATER') \
                    * 1e-6  # MJ/kg-K
        if c.TEMPERATURE:
            self.delta_temp = 1 / spec_heat * (fissions * 180 * 1.6022e-19 - self.beta \
                              / self.kappa * self.delta_vol / 100**3 * self.temp) / (self.mass / 1000)  # K
            self.temp += self.delta_temp  # K
            self.__update_xs_tag()
            self.__update_sab_tag()
        # Don't account delta pressures into the back-end excursion energy calculation
        # This can have negative impacts on the COM acceleration
        if not initial:
            self.delta_pres = self.beta / self.kappa * self.delta_temp - 1 \
                              / (self.kappa * self.volume / 100**3) * self.delta_vol / 100**3  # MPa
            self.av_pressure += self.delta_pres  # MPa
            self.__calc_bottom_pressure(top_pres)

    def __update_com_height(self):
        '''Calculate new height based on volume acceleration and time'''
        # Reliant on Newtonian kinematics equations
        com_vel_i = self.com_vel  # cm/s
        com_vel_f = com_vel_i + self.com_accel * c.DELTA_T  # cm/s
        self.com_vel = com_vel_f  # cm/s, update information for next snapshot
        self.delta_com = com_vel_i * c.DELTA_T + 1 / 2 * self.com_accel * c.DELTA_T**2  # cm
        self.com_height += self.delta_com  # cm

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
        self.mass_h2 += c.RADIOLYTIC_G * energydep * 1000  # g
        self.__update_volfrac_gas()

    def __update_volfrac_gas(self):
        '''Update the volume fraction of gas, f_e; called from self.produce_gas()'''
        # Surface tension is re-called for program clarity (prevents multiple function passes)
        self.volfrac_gas = self.mass_h2 / 1000 * c.RH2 * self.temp \
                           / (self.volume / 100**3
                              * (self.av_pressure + 2 \
                                 * PropsSI('I', 'T', self.temp, 'Q', 0.0, 'WATER')
                                 / c.RAD_GAS_BUBBLE))

    def __calc_bottom_pressure(self, top_pres=None):
        '''Called to calculate bottom pressure, overloaded with top pressure if known'''
        if top_pres is None:
            self.bot_pressure = self.av_pressure + self.dens / 2 * c.GRAV * self.height \
                                * 100**2 / 1000 / 1e6  # MPa
        else:
            self.bot_pressure = 2 * self.av_pressure - top_pres

    def __calc_init(self):
        '''Self-called method to calculate some constants after initial file run'''
        # Ideally only called ONE time after each material definition
        self.base = pi * (self.radius**2 - self.inner_radius**2)  # cm^2
        self.volume = self.base * (self.height - self.base_height)  # cm^3
        self.atoms = [nden * 1e24 * self.volume for nden in self.ndens]  # a
        self.mass = self.dens * self.volume  # g

    def __str__(self):
        # Material representation, overload built-in string definition
        #round_value = 300  # K
        #rounded_temp = int(round_value * round(self.temp / round_value))  # K
        ret = "mat solution{0} sum moder lwtr{0} 1001 tmp {1}\n".format(self.matnum, self.temp)
        template = "{0}.{1} {2}\n"
        for ind, elem in enumerate(self.elems):
            ret += template.format(elem, self.xs_tag, self.ndens[ind])
        return ret + "\n"  # str
