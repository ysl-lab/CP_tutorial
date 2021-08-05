# MODIFIED BY SMM TO TREAT CYCLIC PEPTIDES and Aidan to prevent HIS/unknown residue issue

import sys
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

availRes = ['GLY', 'ALA', 'PRO', 'ASP', 'ASN', 'GLU', 'GLN', 'LYS', 'NLE', 'ARG', 'MET', 'LEU', 'PHE', 'TYR', 'TRP', 'VAL', 'ILE', 'THR', 'SER', 'CYS', 'HID', 'HIE', 'HIP']

ifile = file( ifile_name, 'r' )
Lines = ifile.readlines()
ifile.close

i_res_old = -999
Protein = [ ]
Anames = [ '*' ]

class Residue :
    def __init__(self, i_res, aa) :
        self.i  = i_res
        self.aa = aa
    def Get_Ter_Type(self) :
        if hasattr( self, 'H3' ) :
            self.ter = 'NH3+'
            return
        if hasattr( self, 'OT' ) :
            self.ter = 'COOH'
            self.O1 = self.O
            self.O2 = self.OT
            return
        if hasattr( self, 'OC2' ) :
            self.ter = 'COO-'
            return
        
        self.ter = 'None'

aa = "NULL"
for line in Lines :
    if i_res_old > -999:
        words = line.split()
        if len(words) == 0 :
	    if Res.aa == "HIS":
		while Res.aa != "HID" and Res.aa != "HIE" and Res.aa != "HIP":
			Res.aa = raw_input("Please enter which protonation state of HIS you are using at the C terminal (HID, HIE, or HIP):")
            Protein.append( Res )
            break

	if words[0] == ";" and words[1] == "residue":
		if words[2] == "1" and words[3] == "HIS":
			while aa != "HID" and aa != "HIE" and aa != "HIP":
				aa = raw_input("Please enter which protonation state of HIS you are using at the N terminal (HID, HIE, or HIP):")
		elif words[5] == "HID" or words[5] == "HIE" or words[5] == "HIP":
			aa = words[5]
		else: 
			aa = words[3]
	else:
		i_atom = int(words[0])
		i_res = int(words[2])
		atom = words[4]
		if i_res != i_res_old :
		    try :
			Protein.append( Res )
		    except :
			pass
		    Res = Residue( i_res, aa )
		    i_res_old = i_res
		
		setattr( Res, atom, i_atom )
		
		if len(Anames) == i_atom :
		    Anames.append( atom )
		else :
		    print 'Fatal Error: wrong atom numbers !'
            
    if ';   nr       type  resnr residue' in line :
        i_res_old = 0

Len = len(Protein)
print Len, 'residues  and ', len(Anames)-1, 'atoms'

for res in Protein:
	if res.aa not in availRes:
		print "ERROR: RSFF2 info on ", res.aa, "not available"
		sys.exit()

for i in range( Len ) :
    Protein[i].Get_Ter_Type()

Pair_15 = [ ]    # special 1-5 & 1-6 vdw
Dih = [ ]        # torsions

for i in range( Len ) :
    Res = Protein[i]
    has_Res_prev, has_Res_next = False, False
    if Res.ter not in ['None']:
        print Res.ter
    else :
      #if i > 0 :
      if i >= 0 :
        if i > 0 :
          Res_prev = Protein[i-1]
        if i == 0 :
          Res_prev = Protein[Len-1]
        if Res_prev.aa not in ( 'NA', 'CL' ) :   # true residues
            has_Res_prev = True

      #if i < Len - 1 :
      if i <= Len - 1 :
        if i < Len -1 :
          Res_next = Protein[i+1]
        if i == Len - 1 :
          Res_next = Protein[0]
        if Res_next.aa not in ( 'NA', 'CL' ) :   # true residues
            has_Res_next = True
    
    if Res.aa == 'GLY' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_G_phi) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_G_psi) )
            if Res.ter not in ( 'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )

    if Res.aa == 'ALA' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_A_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_A_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_A_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_A_psi_) )
            if Res.ter not in ( 'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )

    if Res.aa == 'PRO' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_P_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_P_phi_) )
            Dih.append( (Res_prev.C, Res.N, Res.CD, Res.CG, dih_Zeroes) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_P_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_P_psi_) )
        
    if Res.aa == 'ASP' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_D_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_D_phi_) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiD_CG_O) )
            Pair_15.append( (Res_prev.C, Res.CG, vdw_chiD_CG_C) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            if Res.ter not in ('COO-'):
                Pair_15.append( (Res.O, Res.CG, vdw_chiD_CGO) )
                Pair_15.append( (Res.OD1, Res.O, vdw_chiD_ODO) )
                Pair_15.append( (Res.OD2, Res.O, vdw_chiD_ODO) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_D_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_D_psi_) )
            Pair_15.append( (Res.CG, Res_next.N, vdw_chiD_CGN) )
            if Res.ter not in ( 'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
                Pair_15.append( (Res.H, Res.CG, vdw_chiD_HCG) )
    
    if Res.aa == 'ASN' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_N_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_N_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.C, Res.CG, vdw_chiX_CC) )
            if Res.ter not in ('COO-'):
                Pair_15.append( (Res.CG, Res.O, vdw_chiX_CO) )
                Pair_15.append( (Res.ND2, Res.O, vdw_chiN_NDO) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_N_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_N_psi_) )
            Pair_15.append( (Res.CG, Res_next.N, vdw_chiX_CN) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.C, Res.OD1, vdw_chiX_ODC) )
    
    if Res.aa == 'GLU' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_E_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_E_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_E_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_E_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )
    
    if Res.aa == 'GLN' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_Q_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_Q_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_Q_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_Q_psi_) )
            if Res.ter not in ( 'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )
    
    if Res.aa == 'LYS' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_K_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_K_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_K_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_K_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )
    
    if Res.aa == 'NLE' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_K_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_K_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_K_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_K_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )            
    
    if Res.aa == 'ARG' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_R_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_R_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_R_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_R_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )
 
    if Res.aa == 'MET' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_M_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_M_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_M_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_M_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiQ_CGO) )

    if Res.aa == 'LEU' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_L_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_L_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiL_CG_O) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_L_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_L_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiL_CG_O) )

    if Res.aa == 'PHE' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_F_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_F_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiF_CG_O) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N,  dih_F_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_F_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiF_CGO) )
    
    if Res.aa == 'TYR' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_Y_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_Y_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiF_CG_O) )            
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N,  dih_Y_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_Y_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiF_CGO) )

    if Res.aa == 'TRP' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_W_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_W_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiW_CG_O) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_W_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_W_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiW_CGO) )

    if Res.aa == 'VAL' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C, dih_V_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_V_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG1, vdw_chiI_CG_O) )
            Pair_15.append( (Res_prev.O, Res.CG2, vdw_chiI_CG_O) )
#            Pair_15.append( (Res.H, Res.CG1, vdw_chi_HG) )
#            Pair_15.append( (Res.H, Res.CG2, vdw_chi_HG) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_V_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_V_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG1, Res.O, vdw_chiI_CGO) )
            Pair_15.append( (Res.CG2, Res.O, vdw_chiI_CGO) )

    if Res.aa == 'ILE' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C, dih_I_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_I_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG1, vdw_chiI_CG_O) )
            Pair_15.append( (Res_prev.O, Res.CG2, vdw_chiI_CG_O) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_I_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_I_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG1, Res.O, vdw_chiI_CGO) )
            Pair_15.append( (Res.CG2, Res.O, vdw_chiI_CGO) )

    if Res.aa == 'THR' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_T_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_T_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.C, Res.OG1, vdw_chiS_OC) )
            Pair_15.append( (Res.H, Res.CG2, vdw_chiT_CGH) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_T_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_T_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.OG1, Res_next.N, vdw_chiS_ON) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.OG1, Res.O, vdw_chiS_OO) )
            Pair_15.append( (Res.CG2, Res.O, vdw_chiT_CGO) )

    if Res.aa == 'SER' :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_S_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_S_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.C, Res.OG, vdw_chiS_OC) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_S_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_S_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.OG, Res_next.N, vdw_chiS_ON) )
            Pair_15.append( (Res.OG, Res.O, vdw_chiS_OO) )

    if Res.aa in ( 'CYS',) :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_C_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_C_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_C_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_C_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )

    if Res.aa in ( 'HID',) :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_Hd_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_Hd_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiF_CG_O) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_Hd_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_Hd_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiF_CGO) )

    if Res.aa in ( 'HIE',) :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_He_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_He_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiF_CG_O) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_He_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_He_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiF_CGO) )
 
    if Res.aa in ( 'HIP',) :
        if has_Res_prev :
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.C,  dih_He_phi) )
            Dih.append( (Res_prev.C, Res.N, Res.CA, Res.CB, dih_He_phi_) )
            Pair_15.append( (Res_prev.O, Res.C, vdw_O_C) )
            Pair_15.append( (Res_prev.O, Res.CG, vdw_chiF_CG_O) )
        if has_Res_next :
            Dih.append( (Res.N,  Res.CA, Res.C, Res_next.N, dih_He_psi) )
            Dih.append( (Res.CB, Res.CA, Res.C, Res_next.N, dih_He_psi_) )
            if Res.ter not in (  'NH3+' ):
                Pair_15.append( (Res.H, Res_next.N, vdw_H_N) )
                Pair_15.append( (Res.H, Res.O, vdw_H_O) )
            Pair_15.append( (Res.CG, Res.O, vdw_chiF_CGO) )

#    if hasattr( Res, 'CG' ) or hasattr( Res, 'CG1' ) or hasattr( Res, 'SG' )\
#    or hasattr( Res, 'OG1' ) or hasattr( Res, 'OG' ):

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
            
NewLines = [ ]
in_pairs, in_angle, in_dih = False, False, False

for line in Lines :
    words = line.split()
    writen = False
    if in_dih and line[0] != ';' and len(words) > 0 :
        writen = True
        i = int( words[0] )
        j = int( words[1] )
        k = int( words[2] )
        l = int( words[3] )
        found = False
        for each in Dih :
            if ( i == each[0] and j == each[1] and k == each[2] and l == each[3] ) \
            or ( l == each[0] and k == each[1] and j == each[2] and i == each[3] ) :
                NewLines.append( line[:28] + '3' + each[4] + '\n' )     # will have problem if i_atom > 99999
                found = True
                break
        if not found :
            NewLines.append( line )
    
    if len(words) == 0 and in_pairs :
        writen = True
        in_pairs = False
        for each in Pair_15 :
            NewLines.append( '%5i %5i' %(each[0],each[1]) )
            NewLines.append( each[2] + '\n' )
        NewLines.append( '\n[ exclusions ]\n' )
        for each in Pair_15 :
            NewLines.append( '%5i %5i\n' %(each[0],each[1]) )
        NewLines.append( '\n' )
    
    if not writen :
        if 'amber99sb.mod4CPs.ff/forcefield.itp' in line :
            NewLines.append( '#include "RSFF2/RSFF2.itp" \n' )
            
        elif 'amber99sb.mod4CPs.ff/tip3p.itp' in line :
            NewLines.append( '#include "RSFF2/tip3p.itp" \n' )
                 
        else :
            NewLines.append( line )

        if '[ dihedrals ]' in line :
            in_dih = True
        if '[ pairs ]' in line :
            in_pairs = True
        if len(words) == 0 :
            in_pairs, in_angle, in_dih = False, False, False

ofile = file( ofile_name, 'w' )
for line in NewLines :
        ofile.write( line )
ofile.close
