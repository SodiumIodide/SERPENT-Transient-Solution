#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# transientmodel.py

'''
A transient model for critical solution systems.

By Corey Skinner.

Dependencies:
re for regular expression pattern matching in output files
os for calling the sss and bash processes via command interface, file existence for debug
numpy for arrays and mathematical methods
CoolProp for water properties (assuming aqueous mixture is approximated by water)
    - Imported in module file "tm_material"
'''

# External
from os import system, path
import re
import numpy as np

# Shared
from tm_material import Material  # Requires CoolProp
import tm_constants as c
import tm_fileops as fo

def set_materials(elems, ndens, tot_height, tot_radius, **kwargs):
    '''Return a list of materials and geometries, radial columns dominant'''
    mat_counter = 1  # Materials begin with number 1
    materials = []  # Two dimensional list of materials for return
    inner_radius = 0  # cm
    # All materials should start out with the same initial density
    dens = [0.0] * len(elems)  # g/cm^3
    for ind, nden in enumerate(ndens):
        dens[ind] = nden / 6.022e23 * 1e24 * c.AWEIGHT[ind]  # g/cm^3
    den = sum(dens)  # g/cm^3
    # Initial pressures are assumed to be linear, based on rho-g-h model
    # All pressures are absolute, not gauge
    for ind, radius in enumerate(calc_radii(tot_radius)):
        base_height = 0.0  # cm, start at planar origin
        r_list = []  # Second dimension empty list for appending
        # From bottom to top, heights being added
        half_height = 0.0  # m, placeholder
        for height in calc_heights(tot_height):
            # The half_height value is calculated for the initial value from base
            if half_height == 0.0:
                half_height = height / 2
            av_height = (tot_height - height + half_height) / 100  # m
            av_pres = den * c.GRAV * av_height / 1000 * 100**3 / 1e6 + c.ATM  # MPa
            # Equivalent to the height of the center of mass:
            com_height = tot_height - av_height * 100  # cm
            if (ind == 0 and c.INNER_RAD > 0.0) or (ind == c.NUM_RADIAL - 1):
                dissipate = 1
            else:
                dissipate = 0
            if 'temp' in kwargs:
                temperature = kwargs['temp'][mat_counter - 1]  # K
                r_list.append(Material(mat_counter, elems, ndens, den, height,
                                       base_height, radius, inner_radius, com_height,
                                       dissipate, av_pres, temperature))
            else:
                r_list.append(Material(mat_counter, elems, ndens, den, height,
                                       base_height, radius, inner_radius, com_height,
                                       dissipate, av_pres))
            mat_counter += 1
            base_height = height  # cm
        materials.append(r_list)
        inner_radius = radius  # cm
    return materials

def propagate_power(k_eff, lifetime, beta_eff, power):
    '''Propagate the number of neutrons over delta-t'''
    reactivity = (k_eff - 1) / k_eff
    prompt_gen_time = lifetime / k_eff  # s
    return power * np.exp((reactivity - beta_eff) / prompt_gen_time * c.DELTA_T)  # fis/s

def energy_dep_init(k_eff, lifetime, beta_eff):
    '''Calculate the initial energy deposition at the start of the reaction'''
    reactivity = (k_eff - 1) / k_eff
    prompt_gen_time = lifetime / k_eff  # s
    period = prompt_gen_time / (reactivity - beta_eff)  # s
    # Assuming an initiating accident of 1 fission per second at time t=0
    # NOTE: np.log() is the natural logarithm
    time = np.log(c.INIT_POWER / 1) * period  # s
    return period * (np.exp(time / period) - 1) * 1  # fissions

def calc_heights(tot_height):
    '''Returns a list of the ranges of height based on total'''
    height_diff = tot_height / c.NUM_AXIAL  # cm
    heights = list(map(lambda ind: ind * height_diff, range(1, c.NUM_AXIAL + 1)))  # cm
    return heights  # cm

def calc_radii(tot_rad):
    '''Returns a list of the ranges of radii based on total'''
    rad_diff = (tot_rad - c.INNER_RAD) / c.NUM_RADIAL  # cm
    radii = list(map(lambda ind: ind * rad_diff + c.INNER_RAD, range(1, c.NUM_RADIAL + 1)))  # cm
    return radii  # cm

def update_material_states(materials, fissions, tot_height, initial=False):
    '''Updates the state of all materials in profile'''
    counter = 0  # Two inner loops prevent use of enumerate()
    temperatures = []  # K
    pressures = []  # MPa
    heights = []  # cm
    top_pressure = c.ATM  # MPa
    for material_layer in materials:
        bot_height = 0.0  # cm
        for material in material_layer:
            material.update_temp(fissions[counter])
            counter += 1
        if not initial:
            for material in reversed(material_layer):
                material.update_pres()
                material.update_vol(top_pressure)
                top_pressure = material.bot_pressure  # MPa
        for material in material_layer:
            material.update_temp(0.0, decrease=True)
            material.update_pres(decrease=True)
            temperatures.append(material.temp)  # K
            pressures.append(material.av_pressure)  # MPa
            material.shift_height(bot_height)
            bot_height = material.height  # cm
            heights.append(material.height)  # cm
    tot_height = max(heights)  # cm
    return tot_height, temperatures, pressures  # cm, [K], [MPa]

def main():
    '''Main wrapper'''
    print("\nWelcome to the Transient Solution Modeling software.")
    print("Developed by Corey Skinner for the purposes of a revision of")
    print("DOE-HDBK-3010, Chapter 6: Accidental Criticality using the Python 3.6")
    print("programming language in 2017. Supplementary work for a Master's Thesis,")
    print("\"Evaluation of Energy Released in Nuclear Criticality Excursions in")
    print("Process Solutions\"")
    print("\nPlease enter a filename (no extension necessary):")
    filename = input(">>> ")
    if not filename.endswith(".inp"):
        filename += ".inp"
    if not filename.strip():
        print("Must include a filename...")
        raise ValueError
    tot_height = c.INIT_HEIGHT  # cm
    tot_radius = c.RAD  # cm
    # Set materials
    materials = set_materials(c.ELEMS, c.NDENS, tot_height, tot_radius)
    print("Running preliminary file, to determine masses, volumes, etc...")
    timer = 0 #s
    fo.write_file(filename, materials, tot_height)
    outfilename = filename + "_res.m"
    detfilename = filename + "_det0.m"
    # Execute SERPENT calculation
    if not path.isfile(outfilename):
        system("bash -c \"sss {}\"".format(filename))
    temperatures = []  # K
    pressures = []  # MPa
    for material_layer in materials:
        for material in material_layer:
            temperatures.append(material.temp)  # K
            pressures.append(material.av_pressure)  # MPa
    maxtemp = max(temperatures)  # K
    maxpres = max(pressures)  # MPa
    power = c.INIT_POWER  # Start of the flux
    lifetime, keff, keffmax, nubar, beff = fo.get_transient(outfilename)  # s, _, _, n/fis, _
    timer = 0  # s
    integrated_fissions = energy_dep_init(keff, lifetime, beff)  # fissions
    total_fissions = integrated_fissions  # fissions
    number_fissions = power * c.DELTA_T  # fissions
    # Start results file
    with open('results.txt', 'w') as resfile:
        resfile.write("Time (s), Num Fissions, Total Fissions, Max Temperature (K), " + \
                      "Max Pressure (bar), Neutron Lifetime (s), nu-bar, b-eff, k-eff, " + \
                      "k-eff+2sigma, Max Height (cm)\n")
    fo.record(timer, number_fissions, total_fissions, maxtemp, maxpres * 10, lifetime,
              nubar, beff, keff, keffmax, tot_height)
    initial = True  # Flag to calculate integrated energy deposition with initial profile
    # Material addition loop
    print("Beginning main calculation...")
    # # Material expansion loop
    print("Now expanding system by temperature...")
    while keff > c.SUBCRITICAL_LIMIT:
        # Proceed in time
        timer += c.DELTA_T  # s
        # Read previous output file for information and calculate new changes
        fission_profile = fo.count_fissions(detfilename)
        power = propagate_power(keff, lifetime, beff, power)  # fissions/s
        # Correlation between flux profile and fission density
        number_fissions = power * c.DELTA_T  # fissions
        fissions = [frac * number_fissions for frac in fission_profile]  # fissions
        total_fissions += number_fissions  # fissions
        # Assume that initiating fission profile does not change from the initial calculation
        if initial:
            initial = False  # Only calculate one time
            integrated_dist = [frac * integrated_fissions for frac in fission_profile]  # fissions
            _, _, _ = update_material_states(materials, integrated_dist, tot_height,
                                             initial)  # cm, [K], [MPa]
        # Begin total expansion of material
        #if c.EXPANSION:
        #    update_pressures(materials)
        #    update_heights(materials)
        tot_height, temperatures, pressures = update_material_states(materials, fissions,
                                                                     tot_height)  # cm, [K], [MPa]
        maxtemp = max(temperatures)  # K
        maxpres = max(pressures)  # MPa
        timer_string = f"{round(timer, abs(c.TIMESTEP_MAGNITUDE)):.6f}"
        filename = re.sub(r'\d', r'', filename[:filename.rfind(".inp")]).replace('.', '') \
                   + timer_string + ".inp"
        outfilename = filename + "_res.m"
        detfilename = filename + "_det0.m"
        # Do not need to recalculate masses (thus volumes) for materials at this stage
        fo.write_file(filename, materials, tot_height)
        # NOTE: Conditional should be "outfilename", but can be "filename" for debugging purposes
        if not path.isfile(outfilename):
            system("bash -c \"sss {}\"".format(filename))
        lifetime, keff, keffmax, nubar, beff = fo.get_transient(outfilename)
        fo.record(timer, number_fissions, total_fissions, maxtemp, maxpres * 10, lifetime,
                  nubar, beff, keff, keffmax, tot_height)
        print("Current time: {} s".format(round(timer, abs(c.TIMESTEP_MAGNITUDE) + 1)))
        print("Current k-eff: {}".format(keff))
        print("Maximum k-eff: {}".format(keffmax))
        print("Number of fissions: {0:E}".format(sum(fissions)))
        print("Maximum temperature: {}".format(maxtemp))

if __name__ == '__main__':
    try:
        main()
    finally:
        print("\nProgram terminated\n")
