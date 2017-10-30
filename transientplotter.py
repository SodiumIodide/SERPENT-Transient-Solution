#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# transientplotter.py

'''
Python script to plot the output files based on selection.
'''

from os import listdir
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DEFAULT_NAME = 'results.txt'

def main():
    '''Main driver'''
    filename = getfile()
    results_frame = pd.read_csv(filename)
    data = selectdata(results_frame)
    print(results_frame.iloc[:, [0]])
    plt.scatter(results_frame.iloc[:, [0]], results_frame.iloc[:, [1]])
    plt.yscale('log')
    plt.show()

def selectdata(frame):
    '''Return data option'''
    print("Select data to plot:")
    for ind, column in enumerate(frame):
        if ind > 0:
            print(f"{ind}] {column.title()}")
    response = input(">>> ")
    try:
        value = int(response)
    except ValueError:
        exit("Must be a valid option")
    result = value
    return result

def getfile():
    '''Return file name to read'''
    # A default name is assumed, which is the standard name of output files
    chosen_name = DEFAULT_NAME
    # txt files and csv files are prioritized
    filenames = [f for f in listdir('.') if f.endswith('.txt') or f.endswith('.csv')]
    if filenames:
        print("Select filename to use")
    else:
        sys.exit("No files ending with suffix '.txt' or '.csv' located\n")
    # Default name flags
    if DEFAULT_NAME in filenames:
        filenames.remove(DEFAULT_NAME)
        print(f"Press enter to use\n*] {DEFAULT_NAME}\n\nOtherwise press relevant option")
    else:
        chosen_name = ""
    # List applicable files found
    for ind, filename in enumerate(filenames):
        print(f"{ind + 1}] {filename}")
    # Read user input, simply exit if nonsensical
    response = input(">>> ")
    try:
        value = int(response) - 1
        if value > len(filenames) or value < 0:
            sys.exit("No such file\n")
        chosen_name = filenames[value]
    except ValueError:
        if response is not '':
            sys.exit("Not recognized as a valid input\n")
    return chosen_name

if __name__ == '__main__':
    try:
        main()
    finally:
        print("\nProgram terminated\n")
