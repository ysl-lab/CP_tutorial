import re
import MDAnalysis as mda
import os
import sys
import optparse

import numpy as np
from MDAnalysis.lib.distances import calc_dihedrals
from MDAnalysis.analysis.dihedrals import Dihedral
def get_trajectory_files(seq):
    neutral_start_ind = len(seq) * 2
    neutral_end_ind = neutral_start_ind + 5
    neutral_range = [i for i in range(neutral_start_ind, neutral_end_ind)]
    all_files = os.listdir()
    trajectory_files = [i for i in all_files if i[-3:] == "xtc" and int(re.search(r"\d{1,2}", i).group(0)) in neutral_range]
    topology_files = [i for i in all_files if i[-3:] == "gro"]

    trajectory_files.sort()
    return topology_files, trajectory_files

def calculate_omega(u, start, stop, cyclic):
    u_protein = u.select_atoms('protein')
    num_res = len(u_protein.residues)
    ags_omega = [res.omega_selection() for res in u_protein.residues[:-1]]
    if cyclic:
        ags_omega += [u_protein.select_atoms("resid %d and name CA" % num_res, \
                                             "resid %d and name C" % num_res, \
                                             "resid 1 and name N", \
                                             "resid 1 and name CA")]
    R_omega = Dihedral(ags_omega).run(start=start, stop=stop)
    return R_omega.angles

def pick_out_trans(omega_angles, cutoff):
    res = np.ones(len(omega_angles), dtype=bool)
    for i, omega_angle in enumerate(omega_angles):
        if np.all(np.absolute(omega_angle) > cutoff):
            continue
        else:
            res[i] = False

    return res

def pick_out_chirality(chirality_set, seq):
    assert len(seq) == len(chirality_set[0])
    res = np.ones(len(chirality_set), dtype=bool)
    for j, frame in enumerate(chirality_set):
        for i, chirality in enumerate(frame):
            if chirality == "L" and seq[i].isupper():
                continue
            elif chirality == "D" and seq[i].islower():
                continue
            else:
                res[j] = False
    return res

def calculate_chirality(u, start, stop):
    u_protein = u.select_atoms('protein')
    num_frame = len(u.trajectory)
    sequence = [str(i)[9:12] for i in u_protein.residues]
    chirality = np.ones((num_frame, 1))
    for i, aa in enumerate(sequence, start=1):
        if aa == "GLY":
            chirality_temp = np.empty((num_frame, 1), dtype=str)
            chirality_temp.fill("L")
            chirality = np.hstack((chirality, chirality_temp))
        else:
            ags_improper = u_protein.select_atoms("resid %d and name CA" % i, \
                                                  "resid %d and name N" % i, \
                                                  "resid %d and name C" % i, \
                                                  "resid %d and name CB" % i)
            improper_angles = Dihedral([ags_improper]).run(start=start, stop=stop).angles
            chirality_temp = np.empty((num_frame, 1), dtype=str)
            for i, angle in enumerate(improper_angles):
                if angle == 0:
                    sys.exit("\nExiting...Ambigious Improper Dihedral Angle...\n")
                elif angle > 0:
                    chirality_temp[i] = "L"
                else:
                    chirality_temp[i] = "D"
            chirality = np.hstack((chirality, chirality_temp))
    return chirality[:, 1:]

def write_frame_ndx(omega_mask, chirality_mask, file_name):
    total_mask = omega_mask & chirality_mask
    file_name = ".".join([file_name.split(".")[0], "ndx"])
    counter = 0
    with open(file_name, "w+") as fo:
        fo.write("[ frame ]\n")
        for i in total_mask:
            counter += 1
            if i:
                fo.write("%d\n" % counter)


def calculate_phi_psi(u, start, stop, cyclic):
    u_protein = u.select_atoms("protein")
    num_res = len(u_protein.residues)
    num_dihedral = num_res * 2
    ags_phi = [res.phi_selection() for res in u_protein.residues[1:]]
    ags_psi = [res.psi_selection() for res in u_protein.residues[0:-1]]

    if cyclic:
        ags_phi = [u_protein.select_atoms("resid %d and name C" % num_res, \
                                          "resid 1 and name N", \
                                          "resid 1 and name CA", \
                                          "resid 1 and name C")] + ags_phi
        ags_psi += [u_protein.select_atoms("resid %d and name N" % num_res, \
                                           "resid %d and name CA" % num_res, \
                                           "resid 1 and name C", \
                                           "resid 1 and name N")]

    R_phi = Dihedral(ags_phi).run(start=start, stop=stop)
    R_psi = Dihedral(ags_psi).run(start=start, stop=stop)

    phi = R_phi.angles
    psi = R_psi.angles

    dihedral_angle = np.hstack((phi, psi))
    return dihedral_angle

if __name__ == "__main__":
    print("=" * 80)
    print("\nRunning check_trajectory.py......\n")
    print("\n!!!NOTE:Test Version Use With Causion!!!\n")
    print("\n!!!If multiple gro files are found, one will be chosen randomly!!!\n")
    parser = optparse.OptionParser()
    parser.add_option('--seq', dest = 'seq', default = '')
    parser.add_option('--gro', dest = 'gro', default = '')
    parser.add_option('--cutoff', dest = 'cutoff', default = '150')
    parser.add_option('--cyclic', dest = 'cyclic', default = 'True')
    (options, args) = parser.parse_args()
    seq = options.seq
    gro = options.gro

    if options.cyclic.upper()[0] == "T":
        cyclic = True
    else:
        cyclic = False

    cutoff = float(options.cutoff)
    while not seq:
        seq = input("Enter sequence:\n")
    if gro:
        u = mda.Universe(gro)
        omega_angle = calculate_omega(u, 0, 1, cyclic)
        chirality = calculate_chirality(u, 0, 1)
        mask1 = pick_out_trans(omega_angle, cutoff)
        mask2 = pick_out_chirality(chirality, seq)
        if not mask1.flatten()[0]:
            print("\n!!!Cis Bond Found!!!\n")
            input("\nPress Enter to continue...\n")
        else:
            print("\nNo Cis Bond Found\n")
        if not mask2.flatten()[0]:
            print("\n!!!Chirality Not the Same As Declared!!!\n")
            input("\nPress Enter to continue...\n")
        else:
            print("\nChirality Is the Same As Declared\n")

    else:
        topology_files, trajectory_files = get_trajectory_files(seq)
        if not trajectory_files:
            sys.exit("\nExiting...\nNo Trajectory Files Found...\n")
        if len(topology_files) > 1:
            print("\n!!!Multiple Topology Files (.gro) Found!!!\n")
            print("!!!A Random One Will Be Used!!!\n")
        for trajectory_file in trajectory_files:
            print("\nWorking on %s\n" % trajectory_file)
            u = mda.Universe(topology_files[0], trajectory_file)
            num_frame = len(u.trajectory)
            omega_angles = calculate_omega(u, 0, num_frame, cyclic)
            chirality = calculate_chirality(u, 0, num_frame)
            mask1 = pick_out_trans(omega_angles, cutoff)
            mask2 = pick_out_chirality(chirality, seq)
            write_frame_ndx(mask1, mask2, trajectory_file)
            print("\n=> %s <=\n" % trajectory_file)
            print("%.1f%% of the frames have cis bond\n" % (len(np.argwhere(~mask1).flatten()) / num_frame * 100))
            print("%.1f%% of the frames have wrong chirality\n\n" % (len(np.argwhere(~mask2).flatten()) / num_frame * 100))
    print("\nFinish Running check_trajectory.py......\n")
    print("=" * 80)
