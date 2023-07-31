# POD-FDM-FIV
The code and data of 'Reduced-order modeling of nonlinear flow-induced vibration of loosely supported elastic tube via proper orthogonal decomposition'
The remainder codes as follow:
1. beam_vibration_fdm_ddeNSD_FOM.m
   The full order modeling (FOM) finite differences method (FDM) for loosely supported tube bundle flow-induced vibration (FIV) system with various reduced velocity.
2. U_340.mat
   snapshots of flexible tube transverse deflection from 0 to 250 at U=3.4, the first line is the corresponding time.
3. U_400.mat
   snapshots of flexible tube transverse deflection from 0 to 250 at U=4.0, the first line is the corresponding time.
4. U_480.mat
   snapshots of flexible tube transverse deflection from 0 to 250 at U=4.8, the first line is the corresponding time.
5. beam_vibration_fdm_dde_POD_34.m
   The reduced order modeling (ROM) finite differences method (FDM) for loosely supported tube bundle flow-induced vibration (FIV) system at U=3.4.
6. beam_vibration_fdm_dde_POD_40.m
   The reduced order modeling (ROM) finite differences method (FDM) for loosely supported tube bundle flow-induced vibration (FIV) system at U=4.0.
7. beam_vibration_fdm_dde_POD_48.m
   The reduced order modeling (ROM) finite differences method (FDM) for loosely supported tube bundle flow-induced vibration (FIV) system at U=4.8.
8. fixed_end_mode.m
   Dominant modes for linear system and fluid induced vibration systems with different reduced velocity.
9. Relative error_POD_mode.m
   Relative error of truncated POD bases.
