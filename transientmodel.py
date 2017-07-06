#!/usr/bin/env python3

'''
A transient model for critical solution systems.

By Corey Skinner.

Dependencies:
re for regular expression pattern matching in output files
os for calling the sss and bash processes via command interface, file existence for debug
numpy for arrays and mathematical methods
coolprop for water properties (assuming aqueous mixture is approximated by water)
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
    for radius in calc_radii(tot_radius):
        base_height = 0  # cm, start at planar origin
        r_list = []  # Second dimension empty list for appending
        for height in calc_heights(tot_height):
            if 'temp' in kwargs:
                temperature = kwargs['temp'][mat_counter - 1]  # K
                r_list.append(Material(mat_counter, elems, ndens, den, height,
                                       base_height, radius, inner_radius, temperature))
            else:
                r_list.append(Material(mat_counter, elems, ndens, den, height,
                                       base_height, radius, inner_radius))
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
    rad_diff = tot_rad / c.NUM_RADIAL  # cm
    radii = list(map(lambda ind: ind * rad_diff, range(1, c.NUM_RADIAL + 1)))  # cm
    return radii  # cm

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
    number_neutrons = c.INIT_NEUTRONS  # Start of the flux
    lifetime, keff, keffmax, nubar = fo.get_transient(outfilename)  # s, _, _, n/fis
    timer = 0  # s
    number_fissions = number_neutrons / nubar
    total_fissions = number_fissions
    # Start results file
    with open('results.txt', 'w') as resfile:
        resfile.write("Time (s), Num Fissions, Total Fissions, Max Temperature, " + \
                      "Neutron Lifetime(s), k-eff, k-eff+2sigma\n")
    fo.record(timer, number_fissions, total_fissions, maxtemp, lifetime, keff, keffmax)
    # Material addition loop
    print("Beginning main calculation...")
    while keff < 1.01:
        # Proceed in time
        timer += c.DELTA_T  # s
        # Read previous output file for information and calculate new changes
        fission_profile = fo.count_fissions(detfilename)
        number_neutrons = propagate_neutrons(keff, lifetime, number_neutrons)
        # Correlation between flux profile and fission density
        fissions = [frac * number_neutrons / nubar for frac in fission_profile]
        number_fissions = sum(fissions)
        total_fissions += number_fissions
        # Increase total height, requires creating new material information
        tot_height = increase_height(tot_height, 0.2)  # cm
        # Re-apply materials with incresed height
        # (Adding material to system, maintaining even dimension split)
        materials = set_materials(c.ELEMS, c.NDENS, tot_height, tot_radius, temp=temperatures)
        filename = re.sub(r'\d', r'', filename.strip(".inp")) + \
                   re.sub(r'\.', r'', str(round(timer, abs(c.TIMESTEP_MAGNITUDE) + 1))) + \
                   ".inp"
        fo.write_file(filename, materials, tot_height)
        outfilename = filename + "_res.m"
        detfilename = filename + "_det0.m"
        if not path.isfile(outfilename):
            system("bash -c \"sss {}\"".format(filename))
        lifetime, keff, keffmax, nubar = fo.get_transient(outfilename)
        temperatures = []  # K, reset of list
        counter = 0  # Two dimensional loops prevent use of enumerate()
        for material_layer in materials:
            for material in material_layer:
                material.calc_temp(fissions[counter])
                temperatures.append(material.temp)  # K
                counter += 1
        maxtemp = max(temperatures)
        fo.record(timer, number_fissions, total_fissions, maxtemp, lifetime, keff, keffmax)
        print("Current time: {} s".format(round(timer, abs(c.TIMESTEP_MAGNITUDE) + 1)))
        print("Current k-eff: {}".format(keff))
        print("Maximum k-eff: {}".format(keffmax))
        print("Number of fissions: {0:E}".format(sum(fissions)))
        print("Maximum temperature: {}".format(maxtemp))
    # Material expansion loop
    print("Finished adding material")
    print("Now expanding system by temperature...")
    # Store heights in two-dimensional matrix
    heights = np.zeros([c.NUM_RADIAL, c.NUM_AXIAL])
    for rad_ind, material_layer in enumerate(materials):
        for ax_ind, material in enumerate(material_layer):
            heights[rad_ind, ax_ind] = material.height
    with open("results.txt", 'a') as appfile:
        appfile.write("# Expanding material #\n")
    while keff > 1.0:
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
        for rad_ind, material_layer in enumerate(materials):
            for ax_ind, material in enumerate(material_layer):
                if ax_ind != 0:
                    height_shift = heights[rad_ind, ax_ind - 1] - material.base_height
                    material.base_height = heights[rad_ind, ax_ind - 1]
                    material.height += height_shift
                material.expand()
                heights[rad_ind, ax_ind] = material.height
                # Keep checks on total height such that void data doesn't get overwritten
                if tot_height < material.height:
                    tot_height = material.height
        filename = re.sub(r'\d', r'', filename.strip(".inp")) + \
                   re.sub(r'\.', r'', str(round(timer, abs(c.TIMESTEP_MAGNITUDE) + 1))) + \
                   ".inp"
        # Do not need to recalculate masses (thus volumes) for materials at this stage
        fo.write_file(filename, materials, tot_height)
        system("bash -c \"sss {}\"".format(filename))
        outfilename = filename + "_res.m"
        detfilename = filename + "_det0.m"
        lifetime, keff, keffmax, nubar = fo.get_transient(outfilename)
        temperatures = []  # K, reset of list
        counter = 0  # Two dimensional loops prevent use of enumerate()
        for material_layer in materials:
            for material in material_layer:
                material.calc_temp(fissions[counter])
                temperatures.append(material.temp)  # K
                counter += 1
        maxtemp = max(temperatures)
        fo.record(timer, number_fissions, total_fissions, maxtemp, lifetime, keff, keffmax)
        print("Current time: {} s".format(round(timer, abs(c.TIMESTEP_MAGNITUDE) + 1)))
        print("Current k-eff: {}".format(keff))
        print("Maximum k-eff: {}".format(keffmax))
        print("Number of fissions: {0:E}".format(sum(fissions)))
        print("Maximum temperature: {}".format(maxtemp))

if __name__ == '__main__':
    try:
        main()
    except ValueError:
        pass
    finally:
        print("\nProgram terminated\n")
