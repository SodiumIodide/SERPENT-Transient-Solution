'''
Import only module. Contains constants for the Python calculations.

Dependencies:
numpy for mathematical constants (if elected for use)
'''

from numpy import pi as PI

# Calculation parameters
TIMESTEP_MAGNITUDE = -3
DELTA_T = 10**TIMESTEP_MAGNITUDE  # s
INIT_NEUTRONS = 1e16  # Small initiating fission accident source -> Flux build-up
INIT_HEIGHT = 53  # cm
RAD = 15  # cm
NUM_AXIAL = 5  # Currently limited to < 10 by Serpent definitions
NUM_RADIAL = 3  # Currently limited to < 10 by Serpent definitions
NUM_MATERIALS = NUM_AXIAL * NUM_RADIAL  # Equivalent to the number of regions in the model
# Calculated as a threshold value, 1e15 fissions per liter before radiolytic nucleation
THRESHOLD = 1e15 * PI * RAD**2 * INIT_HEIGHT  # fissions/liter

# Equation of state
GRAV = 9.80665  # m/s^2
ATM = 101325  # Pa, atmospheric value
# Radius of radiolytic gas bubbles, from Forehand dissertation, 1981
# -> Independent of pressure, temperature, surface tension, gas and fissile concentration
RAD_GAS_BUBBLE = 5e-6 / 100  # m

# Uranyl nitrate
ELEMS = ["1001", "7014", "8016", "92234", "92235", "92236", "92238"]
NDENS = [6.258e-2, 1.569e-3, 3.576e-2, 1.060e-6, 1.686e-4, 4.350e-7, 1.170e-5]  # a/b-cm
AWEIGHT = [1.0078250321, 14.003074005, 15.9949146196, 234.040952, 235.043930,
           236.045568, 238.050788]  # g/mol

# Radiolytic gas constants
RADIOLYTIC_G = 2.0e-3  # g/J

# Convergence and entropy constants
NUM_PARTICLES = 10000  # Neutrons per cycle
NUM_CYCLES = 500  # Calculation cycles
NUM_DROPS = 25  # Initial cycles to drop from k-convergence

# Thermal scattering cross section for light water
SAB_TEMP = {
    300: "00t",
    350: "02t",
    400: "04t",
    450: "06t",
    500: "08t",
    550: "10t",
    600: "12t",
    650: "14t",
    800: "18t"
}

# Cross section temperatures
XS_TEMP = {
    300: "03c",
    600: "06c",
    900: "09c",
    1200: "12c",
    1500: "15c",
    1800: "18c"
}
