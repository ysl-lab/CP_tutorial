# MODIFIED BY SMM TO TREAT CYCLIC PEPTIDES and Aidan to prevent HIS/unknown residue issue
#Commented by Marc, August 2021
#This program takes an input topology using the modified Amber forcefield and outputs a converted topology using the RSFF2 forcefield. This program was modified to handle cyclic peptides but should also handle uncapped linear peptides and capped linear peptides.
#Example usage: python g_mod_top_RSFF2_JD_0726Comment.py [input forcefield] [output forcefield]
#          e.g. python g_mod_top_RSFF2_JD_0726Comment.py modifiedAmber.top RSFF2.top

#If the script cannot find and replace the lines in the topology for forcefield.itp it will exit with an error.
#If the script cannot find and replace the lines in the topology for tip3p.itp with the RSFF2 folder's version of this file, it will ask for your approval to continue (in case you are using a different water model, thus the tip3p.itp line would not be in the topology.
#If the script cannot find a residue in its database, it will ask for approval to continue (in case the residues not found are caps), you must enter 'y' to get an output in this case.
#The script will print the residues in the sequence, which residues were converted to RSFF2, and the linear/cyclic status of the input.
#The script also will remind users to make sure their structure generation for D-Valine has correctly swapped CG1 and CG2 when switching the chirality from L- to D-
#The program expects residues to begin their numbering at 1 in the input topology.
#Always check that the conversion to RSFF2 did not miss any interactions or dihedral parameters that are specified by RSFF2. Be especially vigilent when dealing with caps and termini, ensuring that no dihedrals or interactions between the cap/termini and its neighbor residues were forgotten.
#The termini supported for uncapped peptides are NH3+ COOH, COO- and NH2 (Make sure that NH2 works correctly, it was added in Aug 2021 by Marc)
import sys

#Define the parameters for RSFF2's special 1-5 and 1-6 vdw interactions
vdw_H_N       = "   1    0.270   5.0  ;  !  H ... +N  destab. (-140,20), alpha shape "
vdw_H_O       = "   1    0.200   2.0  ;  !  H ... O   stabilize C5 conformer "
vdw_O_C       = "   1    0.260   2.0  ;  ! -O ... C   lower the phi = 0 barrier"

vdw_chiS_OO   = "   1    0.360   0.2  ;  ! OG ... O for Ser, Thr "
vdw_chiS_OC   = "   1    0.338   0.2  ;  ! OG ... -C for Ser, Thr "
vdw_chiS_ON   = "   1    0.370   0.2  ;  ! OG ... +N for Ser, Thr "

vdw_chiX_CO   = "   1    0.250   4.0  ;  ! CG ... O  for Asn "
vdw_chiX_CC   = "   1    0.360   0.2  ;  ! CG ... -C for Asn "
vdw_chiX_CN   = "   1    0.320   0.2  ;  ! CG ... +N for Asn "
vdw_chiX_ODC  = "   1    0.260   4.0  ;  ! OD ... C  for Asn "
vdw_chiN_NDO  = "   1    0.380   1.0  ;  ! ND ... O  for Asn "

vdw_chiD_HCG  = "   1    0.290   0.2  ;  ! CG ... H  for Asp "
vdw_chiD_CGO  = "   1    0.250   4.0  ;  ! CG ... O  for Asp "
vdw_chiD_CG_O = "   1    0.280   2.0  ;  ! CG ...-O  for Asp "
vdw_chiD_CG_C = "   1    0.370   0.4  ;  ! CG ...-C  for Asp "
vdw_chiD_CGN  = "   1    0.340   0.2  ;  ! CG ...+N  for Asp "
vdw_chiD_ODO  = "   1    0.340   0.2  ;  ! OD ... O  for Asp "

vdw_chiQ_CGO  = "   1    0.320   4.0  ;  ! CG ... O  for type E,Q,K,R,M "

vdw_chiT_CGH  = "   1    0.320   1.0  ;  ! CG ... H  for Thr "
vdw_chiT_CGO  = "   1    0.330   4.0  ;  ! CG ... O  for Thr "

vdw_chiL_CG_O = "   1    0.340   2.0  ;  ! CG ... O  for Leu "

vdw_chiI_CG_O = "   1    0.320   0.2  ;  ! CG ... -O  for Ile Val "
vdw_chiI_CGO  = "   1    0.320   0.2  ;  ! CG ... O  for Ile Val"

vdw_chiF_CG_O = "   1    0.310   2.0  ;  ! CG ... -O  for Phe and Tyr and His? "
vdw_chiF_CGO  = "   1    0.310   8.0  ;  ! CG ... O  for Phe and Tyr and His? "

vdw_chiW_CG_O = "   1    0.330   2.0  ;  ! CG ... -O  for Trp "
vdw_chiW_CGO  = "   1    0.320   8.0  ;  ! CG ... O  for Trp "

#Define the parameters for the phi, psi, phi prime (phi_), and psi prime (psi_) dihedral angles for each of the residues that RSFF2 adjusts.
dih_G_phi  = "  -2.54    1.92    5.54   -8.86    1.72    3.28    ; ! phi  for Gly "
dih_G_psi  = "   2.25   -4.75  -15.84    7.66    2.52   -1.02    ; ! psi  for Gly "

dih_A_phi  = "  -2.43   -3.84   13.25    8.40   -9.37   -8.27    ; ! phi  for Ala "
dih_A_phi_ = "  -2.46    0.83    8.63    8.89   -4.93  -13.98    ; ! phi_ for Ala "
dih_A_psi  = "   2.94   -1.39  -12.68   -1.56    2.72    3.14    ; ! psi  for Ala "
dih_A_psi_ = "   2.46   -0.25    4.25   -2.99   -7.01    4.03    ; ! psi_ for Ala " 

dih_C_phi  = "  -2.02   -2.44   17.81    4.69  -12.77   -4.53    ; ! phi  for Cys "
dih_C_phi_ = "  -2.94    1.42   12.38    6.50   -4.33  -12.46    ; ! phi_ for Cys "
dih_C_psi  = "   2.02   -0.77  -11.40   -0.50    1.16   -0.13    ; ! psi  for Cys "
dih_C_psi_ = "   2.50   -2.37    2.87   -5.60   -5.92   10.79    ; ! psi_ for Cys " 

dih_P_phi  = "  -1.07    2.56   -0.98  -16.05    0.00    0.00    ; ! phi for Pro "
dih_P_phi_ = "   0.00    0.00    0.00    0.00    0.00    0.00    ; ! phi_ for Pro "
dih_P_psi  = "   2.15   -0.78  -11.03   -4.12   -0.60    9.06    ; ! psi for Pro  "
dih_P_psi_ = "   0.68    2.70   11.35   -0.28  -18.45   -0.64    ; ! psi_ for Pro "

dih_D_phi  = "  -2.93   -1.96   10.87    3.50   -5.61   -4.35    ; ! phi  for Asp "
dih_D_phi_ = "  -2.60    2.15    5.74    2.34   -2.52   -6.48    ; ! phi_ for Asp "
dih_D_psi  = "   2.17    0.90   -5.56   -4.22   -3.54    3.58    ; ! psi  for Asp "
dih_D_psi_ = "   1.48    1.26    4.56   -8.60   -8.07   11.02    ; ! psi_ for Asp "

dih_N_phi  = "  -2.22    0.66   11.97   -2.71   -5.26   -0.56    ; ! phi  for Asn "
dih_N_phi_ = "  -2.86    2.97    2.06    2.35    2.84   -6.56    ; ! phi_ for Asn "
dih_N_psi  = "   2.28   -0.58  -12.16   -2.83    1.82    4.16    ; ! psi  for Asn "
dih_N_psi_ = "   2.08    0.92    0.71   -8.10   -4.24    8.71    ; ! psi_ for Asn "

dih_E_phi  = "  -2.73   -1.00   13.86    4.57   -8.04   -6.31    ; ! phi  for Glu "
dih_E_phi_ = "  -2.79    1.00   13.28    3.32  -10.20   -4.99    ; ! phi_ for Glu "
dih_E_psi  = "   2.20   -1.00   -8.37   -0.90   -0.88    0.19    ; ! psi  for Glu "
dih_E_psi_ = "   1.71   -1.14    5.69   -5.92   -8.46   10.27    ; ! psi_ for Glu "

dih_Q_phi  = "  -9.05    1.17   16.92    3.26   -9.58   -7.60    ; ! phi  for Gln "
dih_Q_phi_ = "  -8.70   -0.57   15.22    8.82   -7.91  -10.30    ; ! phi_ for Gln "
dih_Q_psi  = "   2.79   -0.47   -8.98   -0.11   -1.30    0.00    ; ! psi  for Gln "
dih_Q_psi_ = "   2.15   -1.91    4.00   -3.17   -6.82    8.20    ; ! psi_ for Gln "

dih_K_phi  = "  -2.18    0.17   15.08    2.28   -6.05   -7.58    ; ! phi  for Lys "
dih_K_phi_ = "  -2.54    1.59   14.06    2.90  -10.08   -4.60    ; ! phi_ for Lys "
dih_K_psi  = "   2.57   -0.17  -11.41   -0.41    1.21   -0.33    ; ! psi  for Lys "
dih_K_psi_ = "   2.27   -2.88    3.83   -4.10   -5.53   10.31    ; ! psi_ for Lys " 

dih_R_phi  = "  -2.90   -1.83   20.24    1.76  -13.08   -3.12    ; ! phi  for Arg "
dih_R_phi_ = "  -2.47    1.93   14.25    5.88   -6.12  -14.40    ; ! phi_ for Arg "
dih_R_psi  = "   2.22   -1.09  -12.39   -1.48    1.41    1.46    ; ! psi  for Arg "
dih_R_psi_ = "   2.63   -1.88    3.65   -5.69   -7.08   10.65    ; ! psi_ for Arg " 

dih_M_phi  = "  -2.24   -2.86   19.91    2.37  -15.10   -0.31    ; ! phi  for Met "
dih_M_phi_ = "  -2.83    1.86   10.20    6.65   -4.45  -13.78    ; ! phi_ for Met "
dih_M_psi  = "   2.28   -1.03  -11.39   -5.89    0.60    7.51    ; ! psi  for Met "
dih_M_psi_ = "   2.38    2.70    3.92  -12.41   -8.60   12.63    ; ! psi_ for Met " 

dih_L_phi  = "  -2.47   -1.02   12.75    3.57   -6.88   -8.22    ; ! phi  for Leu "
dih_L_phi_ = "  -2.17   -2.98   16.17    2.39  -13.62   -1.94    ; ! phi_ for Leu "
dih_L_psi  = "   1.91   -0.75  -14.93    0.19    5.66    1.63    ; ! psi  for Leu "
dih_L_psi_ = "   0.84   -0.31   -5.18   -2.42    1.60    5.34    ; ! psi_ for Leu "

dih_F_phi  = "  -7.82   -2.42   15.98    4.88   -9.71   -5.81    ; ! phi  for Phe "
dih_F_phi_ = "  -9.40    1.20    9.63    5.03   -2.32  -12.28    ; ! phi_ for Phe "
dih_F_psi  = "   2.92   -1.34  -10.79   -4.51    1.11    5.61    ; ! psi  for Phe "
dih_F_psi_ = "   2.79    0.69   -1.80   -9.23   -1.84   10.71    ; ! psi_ for Phe " 

dih_Y_phi  = "  -8.19   -2.32   16.29    3.85   -9.84   -4.62    ; ! phi  for Tyr "
dih_Y_phi_ = "  -9.05    1.12    9.43    5.27   -1.95  -12.58    ; ! phi_ for Tyr "
dih_Y_psi  = "   2.30   -0.87  -10.39   -4.54    1.48    5.76    ; ! psi  for Tyr "
dih_Y_psi_ = "   2.57    1.98   -1.51  -10.07   -2.04   11.02    ; ! psi_ for Tyr " 

dih_W_phi  = "  -9.50   -1.28   20.17    2.00  -15.19   -2.55    ; ! phi  for Trp "
dih_W_phi_ = "  -9.95    0.39   16.00    3.66  -13.68   -5.26    ; ! phi_ for Trp "
dih_W_psi  = "   0.67   -1.69   -5.54   -2.24   -2.42    4.10    ; ! psi  for Trp "
dih_W_psi_ = "   0.65    0.40    1.95   -4.62   -6.55    8.07    ; ! psi_ for Trp " 

dih_V_phi  = "  -2.21   -6.14   22.96   11.38  -15.99   -5.77    ; ! phi  for Val "
dih_V_phi_ = "  -2.35    0.10   19.58   10.94  -12.29  -16.63    ; ! phi_ for Val "
dih_V_psi  = "   2.07   -0.18  -18.32   -3.39    5.67    2.08    ; ! psi  for Val "
dih_V_psi_ = "   2.68    0.00   -6.90  -10.60    1.66   16.28    ; ! psi_ for Val " 

dih_I_phi  = "  -2.53   -7.82   26.77   12.70  -21.00   -6.01    ; ! phi  for Ile "
dih_I_phi_ = "  -2.62    0.27   20.57   12.05  -12.06  -19.68    ; ! phi_ for Ile "
dih_I_psi  = "   2.73    0.41  -18.35   -2.56    6.09    0.69    ; ! psi  for Ile "
dih_I_psi_ = "   2.59   -0.75   -6.97  -12.12    3.65   18.45    ; ! psi_ for Ile " 

dih_S_phi  = "  -2.45   -2.34   20.68    1.22  -12.89    0.30    ; ! phi  for Ser "
dih_S_phi_ = "  -2.27    0.57   17.10    7.08   -7.75  -14.11    ; ! phi_ for Ser "
dih_S_psi  = "   2.58   -0.60   -8.40   -2.23   -1.70    3.15    ; ! psi  for Ser "
dih_S_psi_ = "   1.74    2.03    5.95   -5.43  -11.81    9.57    ; ! psi_ for Ser " 

dih_T_phi  = "  -2.88   -4.43   21.90    9.14  -16.02   -5.42    ; ! phi  for Thr "
dih_T_phi_ = "  -2.70   -2.07   22.08    8.27  -14.68  -12.31    ; ! phi_ for Thr "
dih_T_psi  = "   2.53   -1.19   -8.10    1.34   -3.31    0.65    ; ! psi  for Thr "
dih_T_psi_ = "   2.23    1.35   -0.60    7.88   -7.19   -3.93    ; ! psi_ for Thr " 

dih_Hd_phi  = "  -9.80  -0.67   18.40   -3.25  -13.49    2.72    ; ! phi  for HisD "
dih_Hd_phi_ = "  -9.65   2.33   11.33    4.23   -2.21  -13.13    ; ! phi_ for HisD "
dih_Hd_psi  = "   1.64  -0.80   -8.80   -0.93   -0.77    1.56    ; ! psi  for HisD "
dih_Hd_psi_ = "   2.42   0.56   -0.78   -3.79   -2.88    6.06    ; ! psi_ for HisD " 

dih_He_phi  = "  -9.88  -1.43   19.34    6.96  -10.82   -6.14    ; ! phi  for HisE "
dih_He_phi_ = "  -9.79   3.70   13.59    8.58   -5.67  -14.95    ; ! phi_ for HisE "
dih_He_psi  = "   2.74  -1.20  -10.55   -2.24    0.74    3.36    ; ! psi  for HisE "
dih_He_psi_ = "   2.56   0.17    2.33   -4.21   -7.06    6.04    ; ! psi_ for HisE " 

#Define the parameters for the chi1, chi1 prime (chi1_), chi2, chi2 prime (chi2_), chi3, and chi4 dihedral angles for each of the residues that RSFF2 adjusts.
dih_Zeroes  = "   0.0 0.0 0.0 0.0 0.0 0.0 ;  ! zeroes for chi "

dih_E_chi1  = "   8.6    1.7   -2.7   -1.7   0.0   0.0      ;  ! chi1 for Glu "
dih_E_chi1_ = "   8.7    0.7   -3.1    0.0   0.0   0.0      ;  ! chi1_ for Glu "
dih_E_chi2  = "   8.5   -2.7   -0.8    4.5   0.0   0.0      ;  ! chi2 for Glu "
dih_E_chi3  = "  13.1    4.8   -8.5   -2.6   0.0   0.0      ;  ! chi3 for Glu "

dih_Q_chi1  = "   7.3    3.5   -2.7   -3.4   0.0   0.0      ;  ! chi1 for Gln "
dih_Q_chi1_ = "   5.3    0.0   -3.8    0.0   0.0   0.0      ;  ! chi1_ for Gln "
dih_Q_chi2  = "   5.0    1.2    0.1   -1.5   0.0   0.0      ;  ! chi2 for Gln "
dih_Q_chi3  = "   6.1    2.3   -6.4   -8.7   0.0   0.0      ;  ! chi3 for Gln "

dih_D_chi1  = "   3.4    4.2    1.5   -4.7   0.0   0.0      ;  ! chi1 for Asp "
dih_D_chi1_ = "   3.4   -0.2    1.9    0.0   0.0   0.0      ;  ! chi1_ for Asp "
dih_D_chi2  = "   5.0   -2.2  -13.1    2.3   0.0   0.0      ;  ! chi2 for Asp "

dih_N_chi1  = "   4.8    3.9    2.0   -5.8   0.0   0.0      ;  ! chi1 for Asn "
dih_N_chi1_ = "   6.4   -5.8    4.7    0.0   0.0   0.0      ;  ! chi1_ for Asn "
dih_N_chi2  = "   8.4    0.5   -8.4   -9.0   0.0   0.0      ;  ! chi2 for Asn "

dih_K_chi1  = "   8.3    4.4   -2.2   -4.8   0.0   0.0      ;  ! chi1 for Lys "
dih_K_chi1_ = "   4.9   -0.6   -3.6    0.0   0.0   0.0      ;  ! chi1_ for Lys "
dih_K_chi2  = "  10.0    1.5   -1.1    0.4   0.0   0.0      ;  ! chi2 for Lys "
dih_K_chi3  = "  10.1   -3.2   -0.6    4.0   0.0   0.0      ;  ! chi3 for Lys "
dih_K_chi4  = "  10.9    0.2   -4.3    3.9   0.0   0.0      ;  ! chi4 for Lys "

dih_R_chi1  = "   8.7    4.1   -1.8   -5.0   0.0   0.0      ;  ! chi1 for Arg "
dih_R_chi1_ = "   7.0   -0.0   -3.1    0.0   0.0   0.0      ;  ! chi1_ for Arg "
dih_R_chi2  = "   8.1   -0.9   -0.5    2.0   0.0   0.0      ;  ! chi2 for Arg "
dih_R_chi3  = "   8.5    1.1   -1.5    0.2   0.0   0.0      ;  ! chi3 for Arg "
dih_R_chi4  = "   8.4    2.9   -1.6   -0.9   0.0   0.0      ;  ! chi4 for Arg "

dih_M_chi1  = "   7.7    3.2   -0.4   -3.7   0.0   0.0      ;  ! chi1 for Met "
dih_M_chi1_ = "   8.5   -0.9   -2.4    0.0   0.0   0.0      ;  ! chi1_ for Met "
dih_M_chi2  = "   8.7    4.6   -0.8   -7.1   0.0   0.0      ;  ! chi2 for Met "
dih_M_chi3  = "   0.0    3.3    0.6   -1.3   0.0   0.0      ;  ! chi3 for Met "

dih_L_chi1  = "   8.9    5.0    0.4   -5.7   0.0   0.0      ;  ! chi1 for Leu "
dih_L_chi1_ = "   5.1    1.2   -0.1    0.0   0.0   0.0      ;  ! chi1_ for Leu "
dih_L_chi2  = "   2.7    1.3   -9.0    1.1   0.0   0.0      ;  ! chi2 for Leu "

dih_F_chi1  = "   6.9    8.0    4.5  -13.5   0.0   0.0      ;  ! chi1 for Phe and Tyr "
dih_F_chi1_ = "   6.5   -4.3    1.3    0.0   0.0   0.0      ;  ! chi1_ for Phe and Tyr"
dih_F_chi2  = "   4.2    4.5   -2.4   -2.6   0.0   0.0      ;  ! chi2 for Phe and Tyr"

dih_W_chi1  = "   1.2    4.8    2.1   -6.9   0.0   0.0      ;  ! chi1 for Trp "
dih_W_chi1_ = "   1.1   -4.3   -0.3    0.0   0.0   0.0      ;  ! chi1_ for Trp "
dih_W_chi2  = "   0.6   -3.6    0.6    9.4   0.0   0.0      ;  ! chi2 for Trp "

dih_C_chi1  = "   1.7   11.0   -1.5   -9.4   0.0   0.0      ;  ! chi1 for Cys "
dih_C_chi1_ = "   1.2    1.0   -3.7    0.0   0.0   0.0      ;  ! chi1_ for Cys "

dih_S_chi1  = "   1.1    4.4   -1.7   -1.7   0.0   0.0      ;  ! chi1 for Ser "
dih_S_chi1_ = "   1.1    1.8   -3.4    0.0   0.0   0.0      ;  ! chi1_ for Ser "

dih_T_chi1  = "   8.6    7.9    1.5  -12.4   0.0   0.0      ;  ! chi1 for Thr "
dih_T_chi1_ = "   8.4    2.1   -4.8    0.0   0.0   0.0      ;  ! chi1_ for Thr "

dih_V_chi1  = "   4.3    2.4    0.5   -2.9   0.0   0.0      ;  ! chi1 for Val "
dih_V_chi1_ = "   4.3    1.8   -6.9    0.0   0.0   0.0      ;  ! chi1_ for Val "

dih_I_chi1  = "   4.5    9.0    0.4  -11.8   0.0   0.0      ;  ! chi1 for Ile "
dih_I_chi1_ = "   5.1   -2.3   -0.6    0.0   0.0   0.0      ;  ! chi1_ for Ile "
dih_I_chi2  = "   0.0   -0.0   -3.2   -1.3   0.0   0.0      ;  ! chi2 for Ile "

dih_Hd_chi1  = "  5.6    6.7    2.1   -9.7   0.0   0.0      ;  ! chi1 for HisD, HisE "
dih_Hd_chi1_ = "  8.6   -4.5    0.8    0.0   0.0   0.0      ;  ! chi1_ for HisD, HisE "
dih_Hd_chi2  = "  6.8    3.3    4.9   -9.2   0.0   0.0      ;  ! chi2 for HisD, HisE "
dih_Hd_chi2_ = " 11.4   -2.1   -2.5    0.0   0.0   0.0      ;  ! chi2_ for HisD, HisE "

#These were commented out prior to Aug 2021
#dih_He_chi1  = "  0.0    0.0    0.0     0.0    0.0     0.0          ;  ! chi1 for HisE "
#dih_He_chi1_ = "  0.0    0.0    0.0     0.0    0.0     0.0          ;  ! chi1_ for HisE "
#dih_He_chi2  = "  0.0    0.0    0.0     0.0    0.0     0.0          ;  ! chi2 for HisE "
#dih_He_chi2_ = "  0.0    0.0    0.0     0.0    0.0     0.0          ;  ! chi2_ for HisE "

# dih_Hp_chi1  = "  0.0    0.0    0.0     0.0    0.0     0.0          ;  ! chi1 for HisP "
# dih_Hp_chi1_ = "  0.0    0.0    0.0     0.0    0.0     0.0          ;  ! chi1_ for HisP "
# dih_Hp_chi2  = "  0.0    0.0    0.0     0.0    0.0     0.0          ;  ! chi2 for HisP "
# dih_Hp_chi2_ = "  0.0    0.0    0.0     0.0    0.0     0.0          ;  ! chi2_ for HisP "

dih_P_chi1  = "   0.0    0.0    0.0    0.0    0.0     0.0      ;  ! chi1 for Pro "
dih_P_chi1_ = "   0.0    0.0    0.0    0.0    0.0     0.0      ;  ! chi1_ for Pro "
dih_P_chi2  = "   0.0    0.0    0.0    0.0    0.0     0.0      ;  ! chi2 for Pro "
dih_P_chi2_ = "   0.0    0.0    0.0    0.0    0.0     0.0      ;  ! chi2_ for Pro "


from sys import argv
ifile_name = argv[1]
ofile_name = argv[2]
#ifile_name = 'test.old.top'
#ofile_name = 'test.new.top'

#This list is how we record if the input residues are valid for conversion to RSFF2
availRes = ['GLY', 'ALA', 'PRO', 'ASP', 'ASN', 'GLU', 'GLN', 'LYS', 'NLE', 'ARG', 'MET', 'LEU', 'PHE', 'TYR', 'TRP', 'VAL', 'ILE', 'THR', 'SER', 'CYS', 'HID', 'HIE', 'HIP']
#Read the input topology
ifile = file( ifile_name, 'r' )
Lines = ifile.readlines()
ifile.close

i_res_old = -999
#A list to store the amino acid residues in the Protein
Protein = [ ]
#A list that stores the names of all atoms in the topology
Anames = [ '*' ]
#Define a residue with a number, an amino acid identity, and a way to determine if it is a terminal residue.
class Residue :
    def __init__(self, i_res, aa) :
        self.i  = i_res
        self.aa = aa
    def Get_Ter_Type(self) :
        #Search for the terminus NH3+
        if hasattr( self, 'H3' ) :
            self.ter = 'NH3+'
            print self.i, self.aa, "Terminus NH3+"
            return
        if hasattr( self, 'H2') :
            self.ter = 'NH2'
            print self.i, self.aa, "Terminus NH2"
            return
        #Search for the terminus COOH
        if hasattr( self, 'OT' ) :
            self.ter = 'COOH'
            self.O1 = self.O
            self.O2 = self.OT
            print self.i, self.aa, "Terminus COOH"
            return
        #Search for the terminus COO-
        if hasattr( self, 'OC2' ) :
            self.ter = 'COO-'
            print self.i, self.aa, "Terminus COO-"
            return
        #Otherwise there is no program-supported terminus in this residue.
        self.ter = 'None'
        print self.i, self.aa, "No terminus"
#ReadibleProtein = []

aa = "NULL"
#Scan through all lines to collect the amino acids contained in the topology.
for line in Lines :
    #If we have already seen the beginning of the amino acid section, check here for Histidine.
    if i_res_old > -999:
    #All that follows this is indented until if '; nr type etc.'
        words = line.split()
        if len(words) == 0 :
        #Something very strange is happening here Somehow it does not look indented on mine, but it is indented on the cluster
	    if Res.aa == "HIS":
		while Res.aa != "HID" and Res.aa != "HIE" and Res.aa != "HIP":
            	#print "We made it in the strange loop"
			Res.aa = raw_input("Please enter which protonation state of HIS you are using at the C terminal (HID, HIE, or HIP):")
            Protein.append( Res )
            break
    	#This triggers once per each amino acid described in the topology, before it reads the atoms for that amino acid.
	if words[0] == ";" and words[1] == "residue":
        	#If this residue is a Histidine, ask the user what type of histidine you are using.
		if words[2] == "1" and words[3] == "HIS":
			while aa != "HID" and aa != "HIE" and aa != "HIP":
				aa = raw_input("Please enter which protonation state of HIS you are using at the N terminal (HID, HIE, or HIP):")
        	#If Histidine status is not ambiguous in the topology, accept the input topology
		elif words[5] == "HID" or words[5] == "HIE" or words[5] == "HIP":
			aa = words[5]
        	#Read the amino acid type
        	else:
			aa = words[3]
    	#Here it reads the atoms for an amino acid that is not
    	else:
        	#Read the atom number and residue number, as well as the atom name.
		i_atom = int(words[0])
		i_res = int(words[2])
		atom = words[4]
        	# If we are not on the same residue we were last iteration, try to add to the Residue the Protein.
		if i_res != i_res_old :
		    try :
			Protein.append( Res )
		    except :
			pass
            	    #Create the residue based on the defined residue name, update i_res/i_res_old
		    Res = Residue( i_res, aa )
            	    #ReadibleProtein.append([aa])
		    i_res_old = i_res
		#Add the newly read atom to the Residue
	        #ReadibleProtein[i_res].append(atom,i_atom)
		setattr( Res, atom, i_atom )
		
        	#Check the atom count matches up between the list and file.
		if len(Anames) == i_atom :
		    Anames.append( atom )
		else :
		    print 'Fatal Error: wrong atom numbers !'
    
    #This triggers only in the line directly before any amino acid is listed in the topology.
    if ';   nr       type  resnr residue' in line :
        i_res_old = 0
#Print the current status of the peptide.
Len = len(Protein)
print Len, 'residues  and ', len(Anames)-1, 'atoms'
print('The input protein is:')
#This is the structure of the input protein
#readibleCounter =
#for i in ReadibleProtein:

#Here we check if the residues we found are in the list of amino acids that RSFF2 is built for. If not, get user confirmation before continuing.
for res in Protein:
	if res.aa not in availRes:
		print "Error: RSFF2 info on ", res.aa, "not available"
        	isok = raw_input("If this is expected, type 'y'. All other inputs will exit the program: ")
        	if isok != 'y':
            		sys.exit()
    	if res.aa == 'VAL':
        	print "\nNote: Valine detected, if using D-Valine, make sure that your structure correctly swaps Cg1 and Cg2 as discussed in Google Drive 0.Discrepancies...RSFF2_v4.\n"


#Determine what terminal residues exist using the Residue Class's function Get_Ter_Type()
for i in range( Len ) :
    Protein[i].Get_Ter_Type()
print("End of protein \n")
#Tell the user what the protein sequence is and if there are termini.
resCounter = 1
terExist = False
twoTerExist = False
for res in Protein:
    #print("Residue ", resCounter, "has been found to be ", res.aa)
    if Res.ter not in ['None']:
        #print("\nResidue ", resCounter, "has terminus type ",Res.ter,"\n")
    	if terExist == True:
        	twoTerExist = True
    	terExist = True
    resCounter += 1
    #print("\n")
if twoTerExist == False:
    if Protein[0].aa in availRes and Protein[Len-1].aa in availRes:
        print("There are not two termini in this protein, and neither is an RSFF2 parametrized residue, thus the input appears cyclic.\n ")
    elif Protein[0].aa not in availRes and Protein[Len-1].aa not in availRes:
    	print("The input appears to be linear and capped, as both the first and last residues are not parametrized by RSFF2.\n")
elif Protein[0].aa in availRes and Protein[Len-1].aa in availRes:
    print("There are two termini in this protein, both are RSFF2 parametrized residues, thus the input appears to be linear and uncapped.\n ")   
else:
    print("The protein was not successfully identified as either cyclic, linear and uncapped, or linear and capped.")
    


Pair_15 = [ ]    # special 1-5 & 1-6 vdw
Dih = [ ]        # torsions

#Iterate through the built protein, and treat each residue in turn.
for i in range( Len ) :
    #Initialize
    Res = Protein[i]
    has_Res_prev, has_Res_next = False, False
    #If residue already has a terminus, print it. This means that both has_Res_prev or has_Res_next will remain False for terminal residues.
    if Res.ter not in ['None']:
        print Res.ter
        
        #if Res.ter in ( 'NH3+' ) or Res.ter in ( 'NH2' ):
            #Res_next = Protein[i+1]
            #if Res_next.aa not in ( 'NA', 'CL' ) :   # true residues
                #has_Res_next = True
        #if Res.ter in ( 'COOH' ) or Res.ter in ( 'COO-' ) :
            #Res_prev = Protein[i-1]
            #if Res_prev.aa not in ( 'NA', 'CL' ) :   # true residues
                #has_Res_prev = True
    #Otherwise, we must find its neighboring amino acid partners.
    else :
      #Commented out before Aug 2021
      #if i > 0 : 

      #All residues:
      if i >= 0 :
        #If this is not the first residue, the previous residue is the residue before the current residue in Protein.
        if i > 0 :
          Res_prev = Protein[i-1]
        #If this is the first residue, the previous residue is the final residue in Protein
        if i == 0 :
          Res_prev = Protein[Len-1]
        #If the residue is not an ion, it has a previous residue.
        if Res_prev.aa not in ( 'NA', 'CL' ) :   # true residues
            has_Res_prev = True
      #Commented out before Aug 2021
      #if i < Len - 1 :

      #All residues:
      if i <= Len - 1 :
        #If this not the last residue
        if i < Len -1 :
          Res_next = Protein[i+1]
        #If this is the last residue, its next residue is the first residue in Protein
        if i == Len - 1 :
          Res_next = Protein[0]
        #If the residue is not an ion, it has a next residue.
        if Res_next.aa not in ( 'NA', 'CL' ) :   # true residues
            has_Res_next = True
    
    #Here we add interactions between a residue and its neighbors. Create the correct lines in the output topology for dihedrals and 1-5,1-6 interactions by picking the correct atoms in this residue and its neighbors, and adding the predefined dihedral parameters or 1-5 1-6 interactions.
    if Res.aa == 'GLY' :
    
        #As long as Gly is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi of Glycine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_G_phi) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            
        #As long as Gly is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi of Glycine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_G_psi) )
            
            #Only add these interactions if the residue does not have a NH3+ or NH2 terminus
            if Res.ter not in ( 'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )

    if Res.aa == 'ALA' :
        #As long as Ala is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Alanine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_A_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_A_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        #As long as Ala is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Alanine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_A_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_A_psi_) )
            
            #Only add these interactions if the residue does not have a NH3+ or NH2 terminus
            if Res.ter not in ( 'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )

    if Res.aa == 'PRO' :
    
        #As long as Pro is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime and the additional dihedral in Proline
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_P_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_P_phi_) )
            Dih.append( (Res_prev.C, Res.N, Res.CD, Res.CG, dih_Zeroes) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            
        #As long as Pro is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime and the additional dihedral in Proline
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_P_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_P_psi_) )
        
    if Res.aa == 'ASP' :
    
        #As long as Asp is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Aspartic Acid
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_D_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_D_phi_) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiD_CG_O) )
            Pair_15.append( (Res_prev.C, Res.CG, vdw_chiD_CG_C) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            #Only add these interactions if the residue does not have an COO- terminus
            if Res.ter not in ('COO-'):
                Pair_15.append( (Res.O, Res.CG, vdw_chiD_CGO) )
                Pair_15.append( (Res.OD1, Res.O, vdw_chiD_ODO) )
                Pair_15.append( (Res.OD2, Res.O, vdw_chiD_ODO) )
                
        #As long as Asp is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Aspartic Acid
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_D_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_D_psi_) )
            Pair_15.append( (Res.CG, Res_next.N, vdw_chiD_CGN) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in ( 'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
                Pair_15.append( (Res.H, Res.CG, vdw_chiD_HCG) )
    
    if Res.aa == 'ASN' :
    
        #As long as Asn is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Asparagine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_N_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_N_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.C, Res.CG, vdw_chiX_CC) )
            #Only add these interactions if the residue does not have an COO- terminus
            if Res.ter not in ('COO-'):
                Pair_15.append( (Res.CG, Res.O, vdw_chiX_CO) )
                Pair_15.append( (Res.ND2, Res.O, vdw_chiN_NDO) )
                
        #As long as Asn is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Asparagine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_N_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_N_psi_) )
            Pair_15.append( (Res.CG, Res_next.N, vdw_chiX_CN) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.C, Res.OD1, vdw_chiX_ODC) )
    
    if Res.aa == 'GLU' :
        
	#As long as Glu is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Glutamate
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_E_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_E_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            
        #As long as Glu is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Glutamate
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_E_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_E_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+', 'NH2'):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )

    if Res.aa == 'GLN' :
        
	#As long as Gln is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Glutamine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_Q_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_Q_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            
        #As long as Gln is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Glutamine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_Q_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_Q_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in ( 'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )
    
    if Res.aa == 'LYS' :
        
	#As long as Lys is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Lysine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_K_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_K_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            
        #As long as Gln is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Lysine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_K_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_K_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )
    
    if Res.aa == 'NLE' :
        
	#As long as Gln is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Norleucine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_K_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_K_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        
        #As long as Nle is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Lysine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_K_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_K_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )            
    
    if Res.aa == 'ARG' :
    
        #As long as Arg is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Arginine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_R_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_R_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            
        #As long as Arg is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Arginine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_R_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_R_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )
 
    if Res.aa == 'MET' :
    
        #As long as Met is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Methionine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_M_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_M_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            
        #As long as Met is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Methionine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_M_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_M_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )

    if Res.aa == 'LEU' :
    
        #As long as Leu is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Leucine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_L_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_L_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiL_CG_O) )
            
        #As long as Leu is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Leucine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_L_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_L_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiL_CG_O) )

    if Res.aa == 'PHE' :
    
        #As long as Phe is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Phenylalanine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_F_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_F_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiF_CG_O) )
            
        #As long as Phe is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Phenylalanine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N,  dih_F_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_F_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiF_CGO) )
    
    if Res.aa == 'TYR' :
    
        #As long as Tyr is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Tyrosine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_Y_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_Y_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiF_CG_O) )
            
        #As long as Tyr is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Tyrosine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N,  dih_Y_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_Y_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiF_CGO) )

    if Res.aa == 'TRP' :
        
        #As long as Trp is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Tryptophan
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_W_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_W_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiW_CG_O) )
        
        #As long as Trp is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Tryptophan
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_W_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_W_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiW_CGO) )

    if Res.aa == 'VAL' :
        
	#As long as Val is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Valine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C, dih_V_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_V_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG1, vdw_chiI_CG_O) )
            Pair_15.append( (Res_prev.O, Res.CG2, vdw_chiI_CG_O) )
            #These were commented out prior to Aug 2021
#            Pair_15.append( (Res.H, Res.CG1, vdw_chi_HG) )
#            Pair_15.append( (Res.H, Res.CG2, vdw_chi_HG) )

        #As long as Val is not at a terminus
        if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Valine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_V_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_V_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG1, Res.O, vdw_chiI_CGO) )
            Pair_15.append( (Res.CG2, Res.O, vdw_chiI_CGO) )

    if Res.aa == 'ILE' :
	
	#As long as Ile is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Isoleucine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C, dih_I_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_I_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG1, vdw_chiI_CG_O) )
            Pair_15.append( (Res_prev.O, Res.CG2, vdw_chiI_CG_O) )
        
	#As long as Ile is not at a terminus
	if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Isoleucine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_I_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_I_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG1, Res.O, vdw_chiI_CGO) )
            Pair_15.append( (Res.CG2, Res.O, vdw_chiI_CGO) )

    if Res.aa == 'THR' :
        
	#As long as Thr is not at a terminus
	if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Threonine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_T_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_T_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.C, Res.OG1, vdw_chiS_OC) )
            Pair_15.append( (Res.H, Res.CG2, vdw_chiT_CGH) )
        
	#As long as Thr is not at a terminus
	if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Threonine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_T_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_T_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2'):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.OG1, Res_next.N, vdw_chiS_ON) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.OG1, Res.O, vdw_chiS_OO) )
            Pair_15.append( (Res.CG2, Res.O, vdw_chiT_CGO) )

    if Res.aa == 'SER' :
	
	#As long as Ser is not at a terminus
        if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Serine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_S_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_S_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.C, Res.OG, vdw_chiS_OC) )
        
	#As long as Ser is not at a terminus
	if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Serine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_S_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_S_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.OG, Res_next.N, vdw_chiS_ON) )
            Pair_15.append( (Res.OG, Res.O, vdw_chiS_OO) )

    if Res.aa in ( 'CYS',) :
        
	#As long as Cys is not at a terminus
	if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Cysteine
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_C_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_C_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        
	#As long as Cys is not at a terminus
	if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Cysteine
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_C_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_C_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )

    if Res.aa in ( 'HID',) :
        
	#As long as Hid is not at a terminus
	if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Histidine with hydrogen on the delta nitrogen
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_Hd_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_Hd_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiF_CG_O) )
        
	#As long as Hid is not at a terminus
	if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Histidine with hydrogen on the delta nitrogen
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_Hd_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_Hd_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiF_CGO) )

    if Res.aa in ( 'HIE',) :
        
	#As long as Hie is not at a terminus
	if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Histidine with hydrogen on the epsilon nitrogen
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_He_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_He_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiF_CG_O) )
        
	#As long as Hie is not at a terminus
	if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Histidine with hydrogen on the epsilon nitrogen
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_He_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_He_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiF_CGO) )
 
    if Res.aa in ( 'HIP',) :
        
	#As long as Hip is not at a terminus
	if has_Res_prev :
            #Build RSFF2 dihedral for phi and phi prime of Histidine with hydrogen on both nitrogen atoms of the sidechain
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_He_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_He_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiF_CG_O) )
        
	#As long as Hip is not at a terminus
	if has_Res_next :
            #Build RSFF2 dihedral for psi and psi prime of Histidine with hydrogen on both nitrogen atoms of the sidechain
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_He_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_He_psi_) )
            #Only add these interactions if the residue does not have an NH3+ or NH2 terminus
            if Res.ter not in (  'NH3+','NH2' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiF_CGO) )
            
#These were commented out prior to Aug 2021
#    if hasattr( Res, 'CG' ) or hasattr( Res, 'CG1' ) or hasattr( Res, 'SG' )\
#    or hasattr( Res, 'OG1' ) or hasattr( Res, 'OG' ):

    #Here we add dihedrals for chi dihedrals of each residue that RSFF2 has parametrized, within the residue itself.
    if Res.aa == 'PRO' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_P_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_P_chi1_) )
        Dih.append( (Res.CA, Res.N, Res.CD, Res.CG, dih_Zeroes) )
    
    if Res.aa == 'GLU' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_E_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_E_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD, dih_E_chi2) )
        Dih.append( (Res.CB, Res.CG, Res.CD, Res.OE1, dih_E_chi3) )
        Dih.append( (Res.CB, Res.CG, Res.CD, Res.OE2, dih_Zeroes) )
        
    if Res.aa == 'GLN' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_Q_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_Q_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD, dih_Q_chi2) )
        Dih.append( (Res.CB, Res.CG, Res.CD, Res.OE1, dih_Q_chi3) )
        Dih.append( (Res.CB, Res.CG, Res.CD, Res.NE2, dih_Zeroes) )
        
    if Res.aa in ('HID',) :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_Hd_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_Hd_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.ND1, dih_Hd_chi2) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD2, dih_Hd_chi2_) )
    if Res.aa in ('HIE',) :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_Hd_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_Hd_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.ND1, dih_Hd_chi2) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD2, dih_Hd_chi2_) )
    if Res.aa in ('HIP',) :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_Hd_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_Hd_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.ND1, dih_Hd_chi2) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD2, dih_Hd_chi2_) )
    
    if Res.aa == 'VAL' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG2, dih_V_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG2, dih_V_chi1_) )
        
    if Res.aa == 'ILE' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG1, dih_I_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG1, dih_I_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG1, Res.CD, dih_I_chi2) )
        Dih.append( (Res.CG2, Res.CB, Res.CG1, Res.CD, dih_Zeroes) )
    
    if Res.aa == 'THR' :       
        Dih.append( (Res.N, Res.CA, Res.CB, Res.OG1, dih_T_chi1) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD2, dih_Hd_chi2_) )
    
    if Res.aa == 'VAL' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG2, dih_V_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG2, dih_V_chi1_) )
        
    if Res.aa == 'ILE' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG1, dih_I_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG1, dih_I_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG1, Res.CD, dih_I_chi2) )
        Dih.append( (Res.CG2, Res.CB, Res.CG1, Res.CD, dih_Zeroes) )
    
    if Res.aa == 'THR' :       
        Dih.append( (Res.N, Res.CA, Res.CB, Res.OG1, dih_T_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.OG1, dih_T_chi1_) )
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG2, dih_Zeroes) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG2, dih_Zeroes) )
        Dih.append( (Res.CA, Res.CB, Res.OG1, Res.HG1, dih_Zeroes) )
        
    if Res.aa == 'SER' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.OG, dih_S_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.OG, dih_S_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.OG, Res.HG, dih_Zeroes) )
        
    if Res.aa == 'CYS' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.SG, dih_C_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.SG, dih_C_chi1_) )
    
    if Res.aa == 'LEU' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_L_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_L_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD1, dih_L_chi2) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD2, dih_L_chi2) )
        
    if Res.aa in ('PHE', 'TYR') :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_F_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_F_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD1, dih_F_chi2) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD2, dih_Zeroes) )
        
    if Res.aa == 'TRP' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_W_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_W_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD1, dih_W_chi2) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD2, dih_Zeroes) )
        
    if Res.aa == 'MET' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_M_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_M_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.SD, dih_M_chi2) )
        Dih.append( (Res.CB, Res.CG, Res.SD, Res.CE, dih_M_chi3) )
        
    if Res.aa == 'ARG' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_R_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_R_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD, dih_R_chi2) )
        Dih.append( (Res.CB, Res.CG, Res.CD, Res.NE, dih_R_chi3) )
        Dih.append( (Res.CG, Res.CD, Res.NE, Res.CZ, dih_R_chi4) )
        
    if Res.aa in ('LYS', 'NLE') :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_K_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_K_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.CD, dih_K_chi2) )
        Dih.append( (Res.CB, Res.CG, Res.CD, Res.CE, dih_K_chi3) )
        if Res.aa == 'LYS' :
            Dih.append( (Res.CG, Res.CD, Res.CE, Res.NZ, dih_K_chi4) )

    if Res.aa in ('ASP', 'ASPH') :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_D_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_D_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.OD1, dih_D_chi2) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.OD2, dih_Zeroes) )
    
    if Res.aa == 'ASN' :
        Dih.append( (Res.N, Res.CA, Res.CB, Res.CG, dih_N_chi1) )
        Dih.append( (Res.C, Res.CA, Res.CB, Res.CG, dih_N_chi1_) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.OD1, dih_N_chi2) )
        Dih.append( (Res.CA, Res.CB, Res.CG, Res.ND2, dih_Zeroes) )

#Begin to write the new topology file, in_angle seems to be without use.
NewLines = [ ]
in_pairs, in_angle, in_dih = False, False, False
addedITP = False
addedTIP3P = False
#Iterate through our topology file
for line in Lines :
    words = line.split()
    writen = False
    #We now look for matching dihedrals that we constructed for RSFF2 to the dihedral in the current line being read, if the current line is identified as being in the dihedral section of the input topology.
    if in_dih and line[0] != ';' and len(words) > 0 :
        writen = True
        #Initialize the target atoms
        i = int( words[0] )
        j = int( words[1] )
        k = int( words[2] )
        l = int( words[3] )
        found = False
        #Look for matching atoms to replace the existing dihedral with an RSFF2 parametrized dihedral.
        for each in Dih :
            if ( i == each[0] and j == each[1] and k == each[2] and l == each[3] ) \
            or ( l == each[0] and k == each[1] and j == each[2] and i == each[3] ) :
                NewLines.append( line[:28] + '3' + each[4] + '\n' )     # will have problem if i_atom > 99999
                found = True
                break
        #If we did not find a dihedral in our list of RSFF2 dihedrals, keep the line as it was originally in the input.
        if not found :
            NewLines.append( line )
            
    #We now look for the first blank line in the pairs section, so that we may add all the new 1-5 and 1-6 interactions that RSFF2 includes.
    if len(words) == 0 and in_pairs :
        writen = True
        in_pairs = False
        #Cycle through our list of new pairs, and add the atom numbers and the comment.
        for each in Pair_15 :
            print
            NewLines.append( '%5i %5i' %(each[0],each[1]) )
            NewLines.append( each[2] + '\n' )
            
        #Now we have to prevent double counting the interactions as we trainsition from the standard interaction to the bonded vdw interaction in Gromacs. Thus, we add ot the exclusions section the pair interactions just added above.
        NewLines.append( '\n[ exclusions ]\n' )
        for each in Pair_15 :
            NewLines.append( '%5i %5i\n' %(each[0],each[1]) )
        NewLines.append( '\n' )
        
    #If we did not already find something to write for this line iteration, search for the beginning of the next section or other things to replace.
    if not writen :
        #Can we replace the forcefield.itp reference with the RSFF2 forcefield file during this line?
        if 'amber99sbHY.ff/forcefield.itp' in line :
            NewLines.append( '#include "RSFF2/RSFF2.itp" \n' )
            addedITP = True
        #Can we replace the modified amber tip3p.itp file with the file in the RSFF2 folder during this line?
        elif 'amber99sbHY.ff/tip3p.itp' in line :
            NewLines.append( '#include "RSFF2/tip3p.itp" \n' )
            addedTIP3P = True
            #Here track that tip3p.itp was corrected!!!!
        #Keep the original topology line.
        else :
            NewLines.append( line )
        #If this line is the start of the dihedral section track that with the boolean.
        if '[ dihedrals ]' in line :
            in_dih = True
        #If this line is the start of the pairs section track that with the boolean.
        if '[ pairs ]' in line :
            in_pairs = True
        #If the current line is blank, reset all parameters
        if len(words) == 0 :
	    in_pairs, in_angle, in_dih = False, False, False

#Before writing the new topology, make sure everything was correctly completed in constructing the new topology.
if not addedITP:
    print("Fatal error: Modified AMBER force field not found! RSFF2.itp not added.")
    sys.exit()
if not addedTIP3P:
    userinput = raw_input("Error: Modified AMBER tip3p.itp not found! RSFF2/tip3p.itp not added. Is this expected? Type 'y' to continue, all other inputs will exit the program: ")
    if userinput != 'y':
        sys.exit()
ofile = file( ofile_name, 'w' )
for line in NewLines :
        ofile.write( line )
ofile.close

print("Topology conversion complete. Remember to always check that the dihedral parameters match RSFF2 and that there are no extraneous or missing pairs or exclusions that do not make sense for your peptide.")
