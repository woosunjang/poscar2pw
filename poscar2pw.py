#!/usr/local/bin/python

import argparse
import numpy as np
from ase.data import atomic_masses, atomic_numbers
from collections import OrderedDict


def strucparser(infile):
    print("Reading " + str(infile) + "...")

    with open(infile, 'r') as f:
        fileindex = f.readline().strip(' \n')
        scalfac = float(f.readline().strip(' \n'))

        unitvec = []
        for i in range(3):
            unitvec.append(f.readline().split())
        unitvec = np.array(unitvec, dtype='d') * scalfac

        atom_species = f.readline().split()
        atom_counts = f.readline().split()
        atominfo = OrderedDict()
        for i in range(len(atom_species)):
            atominfo[atom_species[i]] = int(atom_counts[i])
        atominfo = atominfo
        count = 0
        for x in atominfo:
            count += atominfo[x]

        determ_line = f.readline()
        if str(determ_line[0]) == "S":
            seldyn = True
            coordtype = f.readline().strip(' \n')
        else:
            seldyn = False
            coordtype = determ_line.strip(' \n')

        cart = None
        if coordtype[0] == 'D':
            cart = False
        elif coordtype[0] == 'd':
            cart = False
        elif coordtype[0] == 'C':
            cart = True
        elif coordtype[0] == 'c':
            cart = True

        tmp = f.read()

        coord = []
        dyn = []
        vel = []

        if not seldyn:
            tmp = np.array(tmp.split(), dtype='d')
            tmp = np.reshape(tmp, (int(len(tmp) / 3), 3))
            dyn = np.array([], dtype='str')
            for i in range(count):
                coord.append(tmp[i])
                if len(tmp) == count * 2:
                    vel.append(tmp[i + count])
            coord = np.array(coord)
            vel = np.array(vel)

        elif seldyn:
            tmp = np.array(tmp.split())
            tmpc = tmp[:(count * 6)]
            if len(tmp) == count * 9:
                vel = np.reshape(tmp[(count * 6):], (count, 3))
            vel = np.array(vel, dtype='d')
            for i in range(count):
                for j in range(3):
                    coord.append(tmpc[i * 6 + j])
                    dyn.append(tmpc[i * 6 + j + 3])
            coord = np.reshape(np.array(coord, dtype='d'), (count, 3))
            dyn = np.reshape(np.array(dyn, dtype='str'), (count, 3))

        atomlist = []
        for i in atominfo:
            for j in range(atominfo[i]):
                atomlist.append(i)
        atomlist = np.vstack(atomlist)

        dic = {"index": fileindex,
               "unitvec": unitvec,
               "atominfo": atominfo,
               "cart": cart,
               "seldyn": seldyn,
               "coord": coord,
               "dyn": dyn,
               "vel": vel,
               "atomlist": atomlist
               }

    return dic


def converter(cellinfo, outfile, writetopw):
    print("Converting to PWscf input lines...")
    print("----------------------------------------------")

    lines = []
    lines.append("\nATOMIC_SPECIES\n")

    for x in cellinfo['atominfo'].keys():
        lines.append(str(x) + " " + str(atomic_masses[atomic_numbers[x]]) + "\n")

    lines.append("\nCELL_PARAMETERS {angstrom}\n")

    for x in np.array(cellinfo['unitvec'], dtype=str):
        lines.append(x[0] + " " + x[1] + " " + x[2] + "\n")

    if cellinfo['cart'] is False:
        lines.append("\nATOMIC_POSITIONS " + "{crystal}\n")
    elif cellinfo['cart'] is True:
        lines.append("\nATOMIC_POSITIONS " + "{angstrom}\n")

    for x in enumerate(cellinfo['coord']):
        lines.append(str(cellinfo['atomlist'][x[0]][0]) + "   " +
                     str(x[1][0]) + " " + str(x[1][1]) + " " + str(x[1][2]) + "\n")

    if writetopw is False:
        for x in lines:
            print(x.strip("\n"))
    else:
        with open(outfile, 'a') as out:
            for x in lines:
                out.write(x)

    return


def executescript(args):
    p = strucparser(args.input)
    converter(p, args.output, args.write)
    return


def main():
    # Main parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-i", dest="input", type=str, default="POSCAR")
    parser.add_argument("-o", dest="output", type=str, default="pw.in")
    parser.add_argument("-w", dest="write", action="store_true")
    parser.set_defaults(func=executescript)

    args = parser.parse_args()

    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        raise AttributeError("Error!")
    args.func(args)


if __name__ == "__main__":
    main()
