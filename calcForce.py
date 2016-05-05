#!/usr/bin/env python3

import sys
import os.path
import math
import enum
import copy
import subprocess

bohr = 0.52917706


class AutoNumber(enum.Enum):
    def __new__(cls):
        value = len(cls.__members__) + 1
        obj = object.__new__(cls)
        obj._value_ = value
        return obj


class MolType(AutoNumber):
    A = ()
    B = ()


class Coordinate:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __sub__(self, other):
        return Coordinate(self.x - other.x, self.y - other.y, self.z - other.z)

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self


class Force:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __iadd__(self, other):
        """

        :type other: Force
        """
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __add__(self, other):
        """

        :type other: Force
        """
        return Force(self.x + other.x, self.y + other.y, self.z + other.z)

    def __isub__(self, other):
        """

        :type other: Force
        """
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def __sub__(self, other):
        """

        :type other: Force
        """
        return Force(self.x - other.x, self.y - other.y, self.z - other.z)

    def __abs__(self):
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def __mul__(self, other):
        """

        :type other: Force
        """
        return self.x * other.x + self.y * other.y + self.z * other.z

    @staticmethod
    def angle(force1, force2):
        """

        :type force1: Force
        :type force2: Force
        """
        return math.acos(force1 * force2 / (abs(force1) * abs(force2)))


def multipy_cord_force(cord, force):
    """

    :type cord: Coordinate
    :type force: Force
    """
    return Force(cord.z * force.y - cord.y * force.z,
                 cord.x * force.z - cord.z * force.x,
                 cord.y * force.x - cord.x * force.y)


class Atom:
    """
    :type atom_typ : int
    :type mol : Molecule
    :type i12_lst : list[Atom]
    """

    def __init__(self):
        self.qm_force = Force()
        self.mm_force = Force()
        self.cord = Coordinate()
        self.seq = 0
        self.name = "Unknown"
        self.atom_typ = None
        self.mol = None
        self.i12_lst = []


def checkMoltype(mol):
    if len(mol.atom_list) == 1:
        return MolType.B
    else:
        return MolType.A


class Molecule:
    """
    :type atom_list : list[Atom]
    """

    def __init__(self):
        self.cord = Coordinate()
        self.atom_list = []
        self.qm_force = Force()
        self.qm_moment_of_force = Force()
        self.mm_force = Force()
        self.mm_moment_of_force = Force()

    def calcCenter(self):
        cord = Coordinate()
        for atom in self.atom_list:
            cord += atom.cord
        leng = len(self.atom_list)
        cord.x /= leng
        cord.y /= leng
        cord.z /= leng
        self.cord = cord


def analyze(xyzfile, keyfile):
    while len(xyzfile) == 0:
        xyzfile = input("Enter xyz file:")
    while len(keyfile) == 0:
        keyfile = input("Enter key file:")
    analyz = subprocess.Popen(["analyze.exe", xyzfile, "-k", keyfile, "F"], stdout=subprocess.PIPE)
    out, err = analyz.communicate()
    lines = out.decode().splitlines()
    it = iter(lines)
    next(it)
    force_map = {}
    try:
        while True:
            line = next(it)
            if len(line) == 0: continue
            if "Forces" in line:
                next(it)
                next(it)
                next(it)
                while True:
                    line = next(it)
                    if len(line) == 0: raise StopIteration
                    keywords = line.split()
                    force_map[int(keywords[0])] = Force(float(keywords[3]),
                                                        float(keywords[4]),
                                                        float(keywords[5]))

    except StopIteration:
        pass

    return force_map


def readGaussianlog(logfile):
    force_map = {}
    with open(logfile) as log:
        it = iter(log)
        next(it)
        try:
            while True:
                line = next(it)
                if len(line) == 0: continue
                if "Forces" in line:
                    next(it)
                    next(it)
                    while True:
                        line = next(it)
                        keywords = line.split()
                        if len(keywords) != 5: raise StopIteration

                        force_map[int(keywords[0])] = Force(float(keywords[2]) * 627.5,
                                                            float(keywords[3]) * 627.5,
                                                            float(keywords[4]) * 627.5)
        except StopIteration:
            pass

        return force_map


def getQMForceMap(logfile):
    return readGaussianlog(logfile)


def getMMForceMap(xyzfile, keyfile):
    return analyze(xyzfile, keyfile)


def getMMForceMapFormFile(file):
    force_map = {}
    with open(sys.argv[2]) as forcefile:
        lines = forcefile.readlines()
        for line in lines:
            if len(line) == 0: continue
            keywords = line.split()
            if len(keywords) == 0: continue
            force_map[int(keywords[0])] = Force(float(keywords[3]),
                                                float(keywords[4]),
                                                float(keywords[5]))
    return force_map


def getForceMap():
    if len(sys.argv) == 3:
        if sys.argv[2][-4:] == ".log":
            return readGaussianlog(sys.argv[2])
        else:
            return getMMForceMapFormFile(sys.argv[2])
    else:
        return analyze(sys.argv[2], sys.argv[3])


def readxyz(xyzfile):
    atom_map = {}
    mol_list = []
    with open(xyzfile) as xyz:
        lines = xyz.readlines()
        line = lines.pop(0)
        keywords = line.split()
        atom_num = int(keywords[0])
        atom_map = {}
        atom_map_list = {}
        for line in lines:
            if len(line) == 0: continue
            keywords = line.split()
            if len(keywords) == 0: continue
            atom = Atom()
            atom.seq = int(keywords[0])
            atom.name = keywords[1]
            atom.cord.x = float(keywords[2])
            atom.cord.y = float(keywords[3])
            atom.cord.z = float(keywords[4])
            atom.atom_typ = int(keywords[5])
            atom_map[atom.seq] = atom
            atom_map_list[atom] = list()
            for atom_num in keywords[6:]:
                atom_map_list[atom].append(int(atom_num))

        for seq, atom in atom_map.items():
            for atom_num in atom_map_list[atom]:
                atom.i12_lst.append(atom_map[atom_num])

        def add_to_mol(atom, mol):
            if atom.mol is not None: return
            mol.atom_list.append(atom)
            atom.mol = mol
            for a in atom.i12_lst:
                add_to_mol(a, mol)

        for seq, atom in atom_map.items():
            if atom.mol is not None: continue
            mol = Molecule()
            add_to_mol(atom, mol)
            mol_list.append(mol)

        for mol in mol_list:
            mol.moltype = checkMoltype(mol)

    return atom_map, mol_list


def main():
    if len(sys.argv) != 5:
        print("Wrong number arguments")
        exit(1)

    atom_map, mol_list = readxyz(sys.argv[1])
    mm_force_map = getMMForceMap(sys.argv[2], sys.argv[3])
    qm_force_map = getQMForceMap(sys.argv[4])

    for seq, atom in atom_map.items():
        atom.qm_force = qm_force_map[seq]
        atom.mm_force = mm_force_map[seq]

    # start to calculate force on molecules

    i = 1
    for mol in mol_list:
        qm_force = Force()
        qm_moment_of_force = Force()
        mm_force = Force()
        mm_moment_of_force = Force()
        mol.calcCenter()
        mol.mol_no = i
        i += 1
        for atom in mol.atom_list:
            qm_force += atom.qm_force
            qm_moment_of_force += multipy_cord_force(atom.cord - mol.cord, atom.qm_force)
            mm_force += atom.mm_force
            mm_moment_of_force += multipy_cord_force(atom.cord - mol.cord, atom.mm_force)

        mol.qm_force = qm_force
        mol.qm_moment_of_force = qm_moment_of_force
        mol.mm_force = mm_force
        mol.mm_moment_of_force = mm_moment_of_force


    # print force to screen
    print("QM     Forces on Molecules (kcal/bohr)")
    print("NO.      type                 X              Y               Z          leng")
    for mol in mol_list:
        print("%d     %s   %17.8f%15.8f%15.8f%15.8f"
              % (mol.mol_no, mol.moltype, mol.qm_force.x, mol.qm_force.y, mol.qm_force.z, abs(mol.qm_force)))

    print("MM     Forces on Molecules (kcal/bohr)")
    print("NO.      type                 X              Y               Z          leng")
    for mol in mol_list:
        print("%d     %s   %17.8f%15.8f%15.8f%15.8f"
              % (mol.mol_no, mol.moltype, mol.mm_force.x, mol.mm_force.y, mol.mm_force.z, abs(mol.mm_force)))

    print("\nQM    Moment of forces on Molecules (kcal)")
    print("NO.      type                 X              Y               Z")
    for mol in mol_list:
        print("%d     %s   %17.8f%15.8f%15.8f" %
              (mol.mol_no, mol.moltype, mol.qm_moment_of_force.x / bohr,
               mol.qm_moment_of_force.y / bohr,
               mol.qm_moment_of_force.z / bohr))

    print("\nMM    Moment of forces on Molecules (kcal)")
    print("NO.      type                 X              Y               Z")
    for mol in mol_list:
        print("%d     %s   %17.8f%15.8f%15.8f" %
              (mol.mol_no, mol.moltype, mol.mm_moment_of_force.x / bohr,
               mol.mm_moment_of_force.y / bohr,
               mol.mm_moment_of_force.z / bohr))




    print("\nQM force angle (degree)")
    print(" mol1 NO.     mol1 NO.      angle(degree)")

    it1 = iter(mol_list)
    try:
        while True:
            mol1 = next(it1)
            it2 = copy.deepcopy(it1)
            try:
                while True:
                    mol2 = next(it2)
                    angle = math.degrees(Force.angle(mol1.qm_force, mol2.qm_force))
                    print("%4d%14d%19.3f" % (mol1.mol_no, mol2.mol_no, angle))
            except StopIteration:
                pass
    except  StopIteration:
        pass

    print("\nMM force angle (degree)")
    print(" mol1 NO.     mol1 NO.      angle(degree)")

    it1 = iter(mol_list)
    try:
        while True:
            mol1 = next(it1)
            it2 = copy.deepcopy(it1)
            try:
                while True:
                    mol2 = next(it2)
                    angle = math.degrees(Force.angle(mol1.mm_force, mol2.mm_force))
                    print("%4d%14d%19.3f" % (mol1.mol_no, mol2.mol_no, angle))
            except StopIteration:
                pass
    except  StopIteration:
        pass
    #
    # mol1 = mol_list.pop(0)
    # mol2 = mol_list.pop(0)
    # cord = mol1.cord - mol2.cord
    # vector12 = Force(cord.x, cord.y, cord.z)
    # angle1 = math.degrees(Force.angle(mol1.force, vector12))
    # angle2 = math.degrees(Force.angle(mol2.force, vector12))
    #
    # print(angle1, angle2)
    #
    # qmforce1 = Force(0.83151531, -1.859342112, -0.911324525)
    # qmforce2 = Force(0.828494525, 0.52586006, 0.265257428)
    # qmangle1 = math.degrees(Force.angle(qmforce1, vector12))
    # qmangle2 = math.degrees(Force.angle(qmforce2, vector12))
    #
    # print(qmangle1, qmangle2)

    print("\n Missing Complete")


if __name__ == "__main__":
    main()
