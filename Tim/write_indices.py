import numpy as np
import optparse
def getPhiPsiOmega(file, cyclic, gro):
    fi= open(file)
    index = {}
    index["N"] = []
    index["CA"] = []
    index["C"] = []
    for line in fi:
        if not gro and line[0:6].strip() != "ATOM":
            continue
        if len(line) < 20:
            continue

        if gro:
            atomName = line[10:15].strip()
            atomID = line[15:20].strip()
        else:
            atomID = line[6:11].strip()
            atomName = line[12:16].strip()

        if atomName in ["C", "N", "CA"]:
            if index[atomName]:
                index[atomName].append(atomID)
            else:
                index[atomName] = [atomID]
    fi.close()

    assert len(index["N"]) == len(index["CA"])
    assert len(index["CA"]) == len(index["C"])

    index_Psi = []
    for i in range(len(index["N"]) -  int(not(cyclic))):
        index_Psi.append(index["N"][i])
        index_Psi.append(index["CA"][i])
        index_Psi.append(index["C"][i])
        index_Psi.append(index["N"][(i + 1) % len(index["N"])])

    index_Phi = []
    for i in range(int(not(cyclic)), len(index["N"])):
        index_Phi.append(index["C"][(i - 1) % len(index["C"])])
        index_Phi.append(index["N"][i])
        index_Phi.append(index["CA"][i])
        index_Phi.append(index["C"][i])

    index_Omega = []
    for i in range(len(index["N"]) - int(not(cyclic))):
        index_Omega.append(index["CA"][i])
        index_Omega.append(index["C"][i])
        index_Omega.append(index["N"][(i + 1) % len(index["N"])])
        index_Omega.append(index["CA"][(i + 1) % len(index["CA"])])

    assert len(index_Psi) == len(index_Phi)
    assert len(index_Omega) == len(index_Psi)

    return (np.array(index_Phi).reshape((-1, 4)).tolist(), \
            np.array(index_Psi).reshape((-1, 4)).tolist(), \
            np.array(index_Omega).reshape((-1, 4)).tolist())

def write_Omega(Omega):
    with open("Omega.ndx", "w+") as fo:
        fo.write("[ Omega ]\n")
        for i in range(len(Omega)):
            fo.write(" ".join(Omega[i]))
            fo.write("\n")

def write_PhiPsi(Phi, Psi):
    with open("PhiPsi.ndx", "w+") as fo:
        fo.write("[ PhiPsi ]\n")
        for i in range(len(Phi)):
            fo.write(" ".join(Phi[i]))
            fo.write("\n")
        for j in range(len(Psi)):
            fo.write(" ".join(Psi[j]))
            fo.write("\n")

def getArgument(arg):
    if arg[0].upper() == "T":
        return True
    return False

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('--PhiPsi', dest = 'PhiPsi',
                      default = "False")
    parser.add_option('--Omega', dest = 'Omega',
                      default = "False")
    parser.add_option('--gro', dest = 'gro',
                      default = '')
    parser.add_option('--pdb', dest = 'pdb',
                      default = '')
    parser.add_option('--cyclic', dest = 'cyclic',
                      default = 'True')
    (options, args) = parser.parse_args()
    pdb = options.pdb
    gro = options.gro
    PhiPsi = getArgument(options.PhiPsi)
    Omega = getArgument(options.Omega)
    cyclic = getArgument(options.cyclic)

    if pdb and gro:
        sys.exit("\nExiting...Both pdb and gro files declared...\n")
    elif not pdb and not gro:
        sys.exit("\nExiting...No pdb or gro file declared...\n")

    if gro:
        inputFile = gro
    else:
        inputFile = pdb

    Phi, Psi, Omega = getPhiPsiOmega(inputFile, cyclic, gro)
    if PhiPsi:
        write_PhiPsi(Phi, Psi)
    if Omega:
        write_Omega(Omega)
