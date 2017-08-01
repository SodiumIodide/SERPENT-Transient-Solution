#!/usr/bin/env python3

'''
A representation of the volume acceleration developed by Adrienne Bobbette Smith
in her 1989 thesis to the University of Arizona: "Nuclear Excursions in Aqueous
Solutions of Fissile Material".

Dependencies:
numpy for matrix computation
'''

import numpy as np

def pressure_matrix(num_axial):
    '''
    Return a matrix of the pressure substitutions
    Atmospheric pressure is always present in the equation system.
    All other pressures are averaged, in order, based on the number of axial
    regions in the system.
    '''
    num_units = num_axial + 1
    p_sub_matrix = np.zeros([num_units, num_units]).astype(int)
    for row_counter in range(num_units):
        inverse_switch = 1
        for col_counter in range(num_units):
            if col_counter == num_axial:
                p_sub_matrix[row_counter, col_counter] = 1 * inverse_switch
                inverse_switch *= -1
            elif col_counter >= num_axial - row_counter:
                p_sub_matrix[row_counter, col_counter] = 2 * inverse_switch
                inverse_switch *= -1
    return np.flipud(p_sub_matrix)

def volume_matrix(num_axial, p_matrix):
    '''
    Perform substitutions to fill in volume matrix
    '''
    num_cols = num_axial + 1
    v_sub_matrix = np.zeros([num_axial, num_cols]).astype(int)
    for row_counter in range(num_axial):
        bar_vector = np.zeros([num_axial]).astype(int)
        vol_num = num_axial - (row_counter + 1)
        inverse_switch = 1
        for col_counter in range(num_axial):
            if col_counter == vol_num:
                bar_vector[col_counter] = 1
                inverse_switch *= -1
            if col_counter > vol_num:
                bar_vector[col_counter] = 2 * inverse_switch
                inverse_switch *= -1
        bar_vector = bar_vector[::-1]
        unit_vector = bar_vector * -1
        atmosphere_correction = 0
        for num, unit in enumerate(unit_vector):
            sub_vector = unit * p_matrix[num + 1][:-1]
            bar_vector += sub_vector
            atmosphere_correction += unit * p_matrix[num + 1][-1]
        for col_counter in range(num_cols):
            if col_counter < num_axial:
                v_sub_matrix[row_counter, col_counter] = bar_vector[col_counter]
            elif col_counter == num_axial:
                v_sub_matrix[row_counter, col_counter] = atmosphere_correction
    return v_sub_matrix

def test():
    '''
    Test of the calculation to ensure that the constants match published and calculated values
    '''
    num_axial = 17  # Can be changed to anything for testing purposes
    p_mat = pressure_matrix(num_axial)
    print("\np_mat:\n{}\n".format(p_mat))
    v_mat = volume_matrix(num_axial, p_mat)
    print("\nv_mat:\n{}\n".format(v_mat))

if __name__ == '__main__':
    try:
        test()
    finally:
        print("\nProgram terminated\n")
