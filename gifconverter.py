#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Converts all .png files in a directory into one sorted-sequence .gif file.

Dependencies:
os for directory reading
imageio for image processing and .gif file creation
'''

from os import listdir
import imageio

def main():
    '''Uses a writer object to project a series of images into a gif file'''
    filenames = [f for f in listdir('.') if f.endswith('.png')]
    with imageio.get_writer('results.gif', mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

if __name__ == '__main__':
    try:
        main()
    finally:
        print("\nProgram terminated\n")
