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
import re
from os import system, path
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
    for radius in calc_radii(tot_radius):
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
            if 'temp' in kwargs:
                temperature = kwargs['temp'][mat_counter - 1]  # K
                r_list.append(Material(mat_counter, elems, ndens, den, height,
                                       base_height, radius, inner_radius, com_height,
                                       av_pres, temperature))
            else:
                r_list.append(Material(mat_counter, elems, ndens, den, height,
                                       base_height, radius, inner_radius, com_height,
                                       av_pres))
            mat_counter += 1
            base_height = height  # cm
        materials.append(r_list)
        inner_radius = radius  # cm
    return materials

def propagate_neutrons(k_eff, lifetime, neutrons):
    '''Propagate the number of neutrons over delta-t'''
    return neutrons * np.exp((k_eff - 1)/lifetime * c.DELTA_T)

# This function is largely present for refactoring purposes
def increase_height(height, incr):
    '''Return an incremented height'''
    return height + incr  # cm

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

def update_com_accel(materials):
    '''
    Updates center of mass acceleration in each material
    Called by the update_heights2() function, not explicitly calculated independently
    '''
    # Materials are arranged by radius, then by height
    # Annular materials are subject to the same pressure fluctuations
    # -> i.e. assume no inter-region mixing, most expansion and transfer will occur
    # -> axially, along the axis of possible expansion (not into walls)
    for ind, material_layer in enumerate(materials):
        # Dissipation effects on the fluid regions next to the wall
        dissipate = 1 if ind == (c.NUM_RADIAL - 1) else 0
        # All materials in each layer should be of the same mass and base area
        mat_mass = material_layer[0].mass / 1000  # kg
        mat_area = material_layer[0].base / 100**2  # m^2
        top_pressure = c.ATM  # MPa, gauge, revised value at each stage
        # Materials/pressures assigned bottom to top, need top to bottom
        for material in reversed(material_layer):
            com_accel = mat_area / mat_mass * (material.bot_pressure - top_pressure) \
                        * 1e6 * 100 - c.DISSIPATION * material.com_vel * dissipate \
                        - c.GRAV * 100  # cm/s^2
            material.set_com_accel(com_accel)
            top_pressure = material.bot_pressure  # MPa, gauge

def update_heights3(materials):
    '''Ensures that the materials are appropriately layered'''
    update_com_accel(materials)
    # Each material needs to have its height and relative base height adjusted
    # -> from the bottom, up; shift heights for no overlap
    # Working on the assumption that pressure acceleration is a relative
    # acceleration value, not an absolute acceleration
    for material_layer in materials:
        base_height = 0.0  # cm, adjusted value by iteration
        for material in material_layer:
            shift = base_height - material.base_height  # cm
            material.height_shift(shift)
            base_height = material.height  # cm, update for next axial layer

def main():
    '''Main wrapper'''
    print("\nWelcome to the Transient Solution Modeling software.")
    print("Developed by Corey Skinner for the purposes of a revision of")
    print("DOE-HDBK-3010, Chapter 6: Accidental Criticality using the Python 3.6")
    print("programming language in 2017.")
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
    if not path.isfile(outfilename):
        system("bash -c \"sss {}\"".format(filename))
    temperatures = []  # K
    for material_layer in materials:
        for material in material_layer:
            temperatures.append(material.temp)
    maxtemp = max(temperatures)  # K
    # update_heights2(materials)
    number_neutrons = c.INIT_NEUTRONS  # Start of the flux
    lifetime, keff, keffmax, nubar = fo.get_transient(outfilename)  # s, _, _, n/fis
    timer = 0  # s
    number_fissions = number_neutrons / nubar
    total_fissions = number_fissions
    # Start results file
    with open('results.txt', 'w') as resfile:
        resfile.write("Time (s), Num Fissions, Total Fissions, Max Temperature, " + \
                      "Neutron Lifetime (s), k-eff, k-eff+2sigma, Max Height (cm)\n")
    fo.record(timer, number_fissions, total_fissions, maxtemp, lifetime, keff, keffmax, tot_height)
    # Material addition loop
    print("Beginning main calculation...")
    # # Material expansion loop
    print("Now expanding system by temperature...")
    # with open("results.txt", 'a') as appfile:
        # appfile.write("# Expanding material #\n")
    while keff > c.SUBCRITICAL_LIMIT:
        # Proceed in time
        timer += c.DELTA_T  # s
        # Read previous output file for information and calculate new changes
        fission_profile = fo.count_fissions(detfilename)
        number_neutrons = propagate_neutrons(keff, lifetime, number_neutrons)
        # Correlation between flux profile and fission density
        fissions = [frac * number_neutrons / nubar for frac in fission_profile]
        number_fissions = sum(fissions)
        total_fissions += number_fissions
        # Begin total expansion of material
        update_heights3(materials)
        counter = 0  # Two inner loops prevent use of enumerate()
        temperatures = []  # K, reset of list
        for material_layer in materials:
            # Pressure goes from top down, so materials reversed
            top_pres = c.ATM  # MPa, gauge, changes by iteration
            for material in reversed(material_layer):
                material.update_state([fis for fis in reversed(fissions)][counter], top_pres)
                top_pres = material.bot_pressure  # MPa, gauge
                # Keep checks on total height such that void data doesn't get overwritten
                if tot_height < material.height:
                    tot_height = material.height
                temperatures.append(material.temp)  # K
                counter += 1
        maxtemp = max(temperatures)  # K
        filename = re.sub(r'\d', r'', filename.rstrip(".inp")) + \
                   re.sub(r'\.', r'', str(round(timer, abs(c.TIMESTEP_MAGNITUDE) + 1))) + \
                   ".inp"
        outfilename = filename + "_res.m"
        detfilename = filename + "_det0.m"
        # Do not need to recalculate masses (thus volumes) for materials at this stage
        fo.write_file(filename, materials, tot_height)
        # NOTE: This should be "outfilename", but can be "filename" for debugging
        if not path.isfile(outfilename):
            system("bash -c \"sss {}\"".format(filename))
        lifetime, keff, keffmax, nubar = fo.get_transient(outfilename)
        fo.record(timer, number_fissions, total_fissions, maxtemp, lifetime, keff, keffmax, tot_height)
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
