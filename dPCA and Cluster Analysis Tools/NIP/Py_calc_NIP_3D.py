import numpy as np
import sys
def main():

    if len (sys.argv) != 3:
       print ("Usage: python Py_calc_NIP_3D.py s1*.txt s2*.txt")
       sys.exit (1)

    global BINS
    BINS = 50

    PC_s1 = np.loadtxt(sys.argv[1])
    PC_s2 = np.loadtxt(sys.argv[2])
    PC_bound = np.concatenate((PC_s1, PC_s2))

    bound_h = np.histogramdd(PC_bound, bins = BINS, density = True)
    s1_h = np.histogramdd(PC_s1, bins = (bound_h[1][0], bound_h[1][1], bound_h[1][2]), density = True)
    s2_h = np.histogramdd(PC_s2, bins = (bound_h[1][0], bound_h[1][1], bound_h[1][2]), density = True)
    ref_h = s1_h

    print("NIP, ref S1: ")
    print("s1:" + str(calc_NIP(s1_h[0],ref_h[0])))
    print("s2:" + str(calc_NIP(s2_h[0],ref_h[0])))

    ref_h = s2_h

    print("NIP, ref S2: ")
    print("s1:" + str(calc_NIP(s1_h[0],ref_h[0])))
    print("s2:" + str(calc_NIP(s2_h[0],ref_h[0])))

def calc_NIP(s_pop,ref_pop):
    numerator = 0.0
    s_denom = 0.0
    ref_denom = 0.0
    for i in range(len(s_pop)):
        for j in range(len(s_pop[i])):
            for k in range(len(s_pop[i][j])):
                numerator += (s_pop[i][j][k] * ref_pop[i][j][k])
                s_denom += (s_pop[i][j][k])**2
                ref_denom += (ref_pop[i][j][k])**2

    NIP = (2*numerator)/(s_denom + ref_denom)

    return NIP

if __name__ == '__main__':
    main()
