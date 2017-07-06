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
                 inner_radius, temp=300):
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
        self.base = 0.0  # cm, placeholder until calculated
        self.delta_temp = 0.0  # K, placeholder until required for expansion
        self.atoms = [0.0] * len(ndens)  # _, placeholder until calculated
        # Pressure: Assume atmospheric conditions (allow for expansion of fluid)
        self.pressure = 101325  # Pa
        self.xs_tag = "03c"
        self.sab_tag = "00t"
        if temp != 300:  # K
            self.update_xs_tag()
            self.update_sab_tag()
        self.__calc_init()

    def append_height(self, newheight):
        '''Append a new volume after adjusting height'''
        self.height = newheight  # cm
        self.volume = self.base * self.height  # cm
        self.ndens = [atom * 1e-24 / self.volume for atom in self.atoms]  # a/b-cm

    def calc_temp(self, fissions):
        '''Adjust material temperature'''
        heatgen = fissions * 180 * 1.6022e-13  # J
        # Based on intensive properties (specific heat) of water
        # This is an aqueous solution of material
        spec_heat = PropsSI('C', 'T', self.temp, 'P', self.pressure, 'WATER')  # J/kg-K
        newtemp = self.temp + heatgen / (self.mass / 1e3) / spec_heat  # K
        self.delta_temp = newtemp - self.temp  # K
        self.temp = newtemp  # K
        self.update_xs_tag()
        self.update_sab_tag()

    def update_xs_tag(self):
        '''Update the materials cross-section tag'''
        round_value = 300
        rounded_temp = int(round_value * round(self.temp / round_value))
        self.xs_tag = c.XS_TEMP.get(rounded_temp, "18c")

    def update_sab_tag(self):
        '''Update the thermal scattering correction tag'''
        round_value = 50
        rounded_temp = int(round_value * round(self.temp / round_value))
        self.sab_tag = c.SAB_TEMP.get(rounded_temp, "18t")

    def expand(self):
        '''Self-regulate the expansion of the material'''
        # The isobaric expansion coefficient of water is calculated
        expansion = PropsSI('ISOBARIC_EXPANSION_COEFFICIENT', 'T', self.temp,
                            'P', self.pressure, 'WATER')  # 1/K
        # Coefficient is defined as: alpha = 1/V * (dV/dT) over constant P
        # Where V is volume, T is temperature, P is pressure, and alpha is the coefficient
        # In this case, changing differential into delta, and solving for a delta height
        # while maintaining base area
        delta_height = self.delta_temp * self.volume * expansion / self.base  # cm
        newheight = self.height + delta_height  # cm
        self.append_height(newheight)  # Book-keeping for number densities and such

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
