import optparse

def get_improper_atom(gro):
    fi = open(gro)
    lines = fi.readlines()
    fi.close()
    length = lines[-2][0:5].strip()
    target = ["CA", "C", "N", "H", "O"]
    mem = {}
    for line in lines:
        resID = line[0:5].strip()
        if resID != "1" and resID != length:
            continue

        atomName = line[10:15].strip()
        if atomName not in target:
            continue

        atomID = line[15:20].strip()
        entryID = atomName + resID
        mem[entryID] = atomID
    return length, mem

def add_improper(topOri, topOut, improper1, improper2):
    fi = open(topOri)
    lines = fi.readlines()
    fi.close()
    with open(topOut, "w+") as fo:
        for i in range(1, len(lines)):
            if lines[i] == "; Include Position restraint file\n":
                fo.write(improper1)
                fo.write("\n")
                fo.write(improper2)
                fo.write("\n")
            fo.write(lines[i - 1])
        fo.write(lines[-1])


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('--gro', dest = 'gro', default = 'prot.gro')
    parser.add_option('--ori', dest = 'topOri', default = 'cx_amber99sbMod_tip3p_temp.top')
    parser.add_option('--out', dest = 'topOut', default = 'cx_amber99sbMod_tip3p.top')

    (options, args) = parser.parse_args()
    gro = options.gro
    topOri = options.topOri
    topOut = options.topOut


    length, improper_dic = get_improper_atom(gro)
    improper1 = "%5s %5s %5s %5s %5s     ; added by TL" % (improper_dic["CA" + length],
                                                           improper_dic["N1"],
                                                           improper_dic["C" + length],
                                                           improper_dic["O" + length],
                                                           4)
    improper2 = "%5s %5s %5s %5s %5s     ; added by TL" % (improper_dic["C" + length],
                                                           improper_dic["CA1"],
                                                           improper_dic["N1"],
                                                           improper_dic["H1"],
                                                           4)
    add_improper(topOri, topOut, improper1, improper2)
