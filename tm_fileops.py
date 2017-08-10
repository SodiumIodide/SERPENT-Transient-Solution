'''Import only module. Contains operations on input and output files.'''

import re
import tm_constants as c

def count_fissions(filename):
    '''Read an output file to determine the fission (flux) distribution'''
    pat = re.compile(r'(\s+1)\1+\s+(\S+)')
    profile = []
    counter = 0
    with open(filename, mode='r') as detfile:
        for line in detfile:
            if re.match(pat, line):
                matches = re.findall(pat, line)
                if counter != 0:
                    profile.append(float(matches[0][1]))
                counter += 1
    return profile

def get_transient(filename):
    '''
    Read an output file to determine the k-eff and neutron lifetime
    Returns tuple of lifetime, k-eff, and maximum k-eff
    '''
    ltpat = re.compile(r'IMPL_PROMPT_LIFETIME\s+\(idx,\s\[1:\s+2\]\)\s=\s\[\s+(\S+)\s\S+\s\];')
    kepat = re.compile(r'IMP_KEFF\s+\(idx,\s\[1:\s+2\]\)\s=\s\[\s+(\S+)\s(\S+)\s\];')
    nbpat = re.compile(r'NUBAR\s+\(idx,\s\S+\s+\S+\)\s=\s\[\s+(\S+).+;')
    with open(filename, mode='r') as ofile:
        for line in ofile:
            if re.match(ltpat, line):
                matches = re.findall(ltpat, line)
                lifetime = float(matches[0])  # s
            if re.match(kepat, line):
                matches = re.findall(kepat, line)
                keff = float(matches[0][0])
                maxkeff = round(keff + 2 * float(matches[0][1]), 5)
            if re.match(nbpat, line):
                matches = re.findall(nbpat, line)
                nubar = float(matches[0])  # n/fis
    return (lifetime, keff, maxkeff, nubar)  # s, _, _, n/fis

def write_file(filename, materials, tot_height):
    '''Function to create the series of input files'''
    with open(filename, mode='w', newline='\n') as fhan:
        fhan.write("% Serpent Input File\n")
        fhan.write("set title \"solutionmodel\"\n\n")
        # Material data
        fhan.write("\n% Materials\n")
        for material_level in materials:
            for material in material_level:
                fhan.write(str(material))
        # Surface data
        fhan.write("\n% Surfaces\n")
        # Axial distributions occur on order nx, where n is the material number
        # and x is the plane number (1 for bottom and 2 for top)
        # Radial distributions read as n, where n is increasing to NUM_RADIAL
        for rad_ind, material_level in enumerate(materials):
            fhan.write("surf {0} cyl 0 0 ".format(rad_ind + 1))
            for mat_ind, material in enumerate(material_level):
                if mat_ind == 0:
                    fhan.write("{}\n".format(material.radius))
                    fhan.write("surf {0}0 pz {1}\n".format(rad_ind + 1, material.base_height))
                fhan.write("surf {0}{1} pz {2}\n".format(rad_ind + 1, mat_ind + 1, material.height))
        # Overall boundaries
        fhan.write("surf 1001 pz 0\n")
        fhan.write("surf 1002 pz {}\n\n".format(tot_height))
        # Cell data
        fhan.write("\n% Cells\n")
        # Cells are numbered by n, where n is the material number
        for rad_ind, material_level in enumerate(materials):
            # Void super-radials are marked with 10n, where n is the radial number
            # Void sub-radials are marked with 11n, where n is the radial number
            fhan.write("cell 10{0} {0} void {0}\n".format(rad_ind + 1))
            if rad_ind > 0:
                fhan.write("cell 11{0} {0} void -{1}\n".format(rad_ind + 1, rad_ind))
            # Void sub-axials are marked with 12n, where n is the radial number
            for mat_ind, material in enumerate(material_level):
                if mat_ind == 0:
                    fhan.write("cell 12{0} {0} void -{0}0\n".format(rad_ind + 1))
                fhan.write("cell {0} {1} solution{0} {1}{2} -{1}{3} -{1}"
                           .format(material.matnum, rad_ind + 1, mat_ind, mat_ind + 1))
                if rad_ind > 0:
                    fhan.write(" {}".format(rad_ind))
                fhan.write("\n")
            # Void super-axials are marked with 13n, where n is the radial number
            fhan.write("cell 13{0} {0} void {0}{1}\n".format(rad_ind + 1, c.NUM_AXIAL))
            fhan.write("\n")
        # Global universe cells, marked with 100n, where n is the radial number
        for radial_num in range(1, c.NUM_RADIAL + 1):
            fhan.write("cell 100{0} 0 fill {0} 1001 -1002 -{0}".format(radial_num))
            if radial_num > 1:
                fhan.write(" {}".format(radial_num - 1))
            fhan.write("\n")
        # Now mark global outside, still of form 100n, but with n > radial number
        fhan.write("cell 100{} 0 outside -1001\n".format(c.NUM_RADIAL + 1))
        fhan.write("cell 100{} 0 outside 1002\n".format(c.NUM_RADIAL + 2))
        fhan.write("cell 100{0} 0 outside {1}\n".format(c.NUM_RADIAL + 3, c.NUM_RADIAL))
        fhan.write("\n")
        # Thermal scattering S(a,b)
        fhan.write("\n% Thermal scattering library\n")
        for material_level in materials:
            for material in material_level:
                fhan.write("therm lwtr{0} lwe7.{1}\n".format(material.matnum, material.sab_tag))
        fhan.write("\n")
        # Cross-section library
        fhan.write("\n% Cross-section library\n")
        fhan.write("set acelib \"/xs/sss_endfb7u.xsdata\"\n\n")
        # Criticality parameters
        fhan.write("\n% Criticality parameters\n")
        # Leave a guess of 1.0 for k-eff
        fhan.write("set pop {0} {1} {2} 1.0\n\n"
                   .format(c.NUM_PARTICLES, c.NUM_CYCLES, c.NUM_DROPS))
        # Include unresolved resonance probability tables
        fhan.write("\n% Unresolved resonance probability calculations\n")
        fhan.write("set ures 1\n\n")
        # Geometry plot
        fhan.write("\n% Geometry plot\n")
        fhan.write("plot 1 1000 1000 0 -{0} {0} -1 {1}\n\n".format(c.RAD + 1, tot_height + 1))
        # Flux detectors
        fhan.write("\n% Flux detectors\n")
        # Total detector of fissions
        fhan.write("det 100 du 0 dr -6 void\n")
        for material_level in materials:
            for material in material_level:
                fhan.write("det {0} dc {0} dr -6 void dt 3 100\n".format(material.matnum))

def record(time, numfissions, totfissions, maxt, lifetime, keff, keffmax):
    '''Record current step to a results file'''
    with open("results.txt", 'a') as appfile:
        appfile.write("{0}, {1:E}, {2:E}, {3}, {4}, {5}, {6}\n"
                      .format(round(time, abs(c.TIMESTEP_MAGNITUDE) + 1), numfissions,
                              totfissions, maxt, lifetime, keff, keffmax))
