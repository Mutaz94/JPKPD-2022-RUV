Sat Sep 18 15:54:52 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/All/dat51.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK

ET1 = EXP(ETA(1)*THETA(4))
ET2 = EXP(ETA(2)*THETA(5))
ET3 = EXP(ETA(3)*THETA(6))


CL = 5.0 * THETA(1) * ET1
V = 85  * THETA(2) * ET2
KA = 0.7 * THETA(3) * ET3

SC = V
$ERROR
CVERR 	= 0.05
W  	= THETA(7)*F*CVERR
Y  	= F + W * ERR(1)
$THETA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvKA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvK
(0,1) ; RUV
$OMEGA
0.9 FIX ;     IIV CL
0.9 FIX  ;     IIV V
0.9 FIX ;      IIV KA
$SIGMA  1  FIX;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E+00
 0.0000E+00   0.9000E+00
 0.0000E+00   0.0000E+00   0.9000E+00
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
   ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3

 TRANSLATOR WILL CONVERT PARAMETERS
 CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            4           *           *           *           *
    3            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m51.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:   10526.9288803176        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   1.0208E+02  2.3451E+01 -1.1617E+02 -7.3889E+01  2.5269E+01 -4.9714E+01 -2.2929E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:   1436.48864911310        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:       61
 NPARAMETR:  1.8480E+00  1.4013E+00  2.3847E+00  8.4245E-01  7.0008E-01  1.1157E+00  2.2573E+00
 PARAMETER:  7.1408E-01  4.3739E-01  9.6906E-01 -7.1444E-02 -2.5657E-01  2.0951E-01  9.1416E-01
 GRADIENT:   2.6554E+02  1.4187E+02  2.3912E+01 -1.0699E+02 -1.7012E+02  1.0736E+01 -4.4988E+03

0ITERATION NO.:   10    OBJECTIVE VALUE:  -522.329581848729        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      115
 NPARAMETR:  1.3596E+00  1.4917E+00  2.3834E+00  8.8077E-01  1.2709E+00  3.5621E-01  7.5535E+00
 PARAMETER:  4.0716E-01  4.9989E-01  9.6854E-01 -2.6962E-02  3.3973E-01 -9.3223E-01  2.1220E+00
 GRADIENT:   5.5642E+01  1.2308E+01 -3.3760E+00  6.2694E+01  1.1745E+02  1.0385E+00 -1.6506E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -610.132985286958        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      165
 NPARAMETR:  1.0260E+00  1.2866E+00  4.8361E+00  5.1569E-01  3.8601E-01  3.9698E-01  1.0848E+01
 PARAMETER:  1.2564E-01  3.5197E-01  1.6761E+00 -5.6224E-01 -8.5190E-01 -8.2386E-01  2.4840E+00
 GRADIENT:  -8.6288E+01  1.9595E+01 -1.9328E+00  7.2927E+00  2.7258E+01 -1.8228E-02  5.5160E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -625.085788818109        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      215
 NPARAMETR:  1.0721E+00  1.1135E+00  7.1040E+00  3.9177E-01  1.7110E-01  4.7713E-01  1.0931E+01
 PARAMETER:  1.6965E-01  2.0754E-01  2.0607E+00 -8.3708E-01 -1.6655E+00 -6.3996E-01  2.4916E+00
 GRADIENT:   4.6557E+00 -1.0505E+01 -1.4322E+00 -9.3833E+00  3.2632E+00 -3.0671E-02 -8.4938E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -625.864828068635        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      265
 NPARAMETR:  1.0588E+00  1.1008E+00  9.9159E+00  3.9445E-01  1.1594E-01  5.7042E-01  1.1204E+01
 PARAMETER:  1.5709E-01  1.9607E-01  2.3941E+00 -8.3026E-01 -2.0547E+00 -4.6139E-01  2.5163E+00
 GRADIENT:  -9.0448E-01  1.0329E-01 -1.0439E+00  1.9780E-01 -8.3250E-01 -2.2467E-02  8.6397E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -626.913829780104        NO. OF FUNC. EVALS.:  56
 CUMULATIVE NO. OF FUNC. EVALS.:      321
 NPARAMETR:  1.0644E+00  1.1070E+00  2.8552E+02  4.0285E-01  1.4524E-01  3.0119E+00  1.1172E+01
 PARAMETER:  1.6242E-01  2.0161E-01  5.7543E+00 -8.0918E-01 -1.8294E+00  1.2026E+00  2.5134E+00
 GRADIENT:   5.3955E+00 -6.2214E+00 -3.3388E-02  2.7732E+00  6.4231E-01 -6.9736E-04  4.9432E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -626.961876462726        NO. OF FUNC. EVALS.:  53
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0554E+00  1.1013E+00  1.8708E+03  3.9692E-01  1.3457E-01  7.7443E+00  1.1141E+01
 PARAMETER:  1.5389E-01  1.9652E-01  7.6341E+00 -8.2402E-01 -1.9057E+00  2.1470E+00  2.5106E+00
 GRADIENT:  -2.5913E-01 -1.2832E-01 -5.2255E-03 -1.1448E-01  2.0493E-03 -1.1205E-04  2.8915E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -626.968222912649        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      425
 NPARAMETR:  1.0558E+00  1.1019E+00  1.9285E+04  3.9736E-01  1.3535E-01  2.4782E+01  1.1136E+01
 PARAMETER:  1.5431E-01  1.9705E-01  9.9671E+00 -8.2291E-01 -1.8999E+00  3.3101E+00  2.5102E+00
 GRADIENT:   3.5663E-02  2.6390E-02 -5.0392E-04  2.1816E-02 -1.1378E-03 -1.0768E-05 -5.4136E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -626.968492422120        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      477
 NPARAMETR:  1.0557E+00  1.1018E+00  3.9589E+05  3.9729E-01  1.3523E-01  1.1179E+02  1.1137E+01
 PARAMETER:  1.5424E-01  1.9696E-01  1.2989E+01 -8.2310E-01 -1.9008E+00  4.8166E+00  2.5103E+00
 GRADIENT:  -3.3692E-03 -2.2937E-03 -2.4528E-05 -1.7686E-03  8.4811E-05 -5.1844E-07  4.4966E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -626.972387741512        NO. OF FUNC. EVALS.: 113
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  1.0566E+00  1.1029E+00  1.5302E+13  3.9786E-01  1.3620E-01  6.7714E+05  1.1143E+01
 PARAMETER:  1.5509E-01  1.9798E-01  3.0459E+01 -8.2166E-01 -1.8936E+00  1.3526E+01  2.5108E+00
 GRADIENT:  -6.7885E-01 -4.5291E-01 -1.8662E-12 -6.6833E-01 -7.1657E-02  0.0000E+00 -1.5221E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -626.977064560250        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:      699
 NPARAMETR:  1.0594E+00  1.1061E+00  5.6118E+12  3.9939E-01  1.3813E-01  6.7257E+05  1.1163E+01
 PARAMETER:  1.5774E-01  2.0081E-01  2.9456E+01 -8.1783E-01 -1.8795E+00  1.3519E+01  2.5126E+00
 GRADIENT:   4.1907E-01  6.6956E-02 -1.9298E-12 -2.5830E-02 -5.8343E-02  0.0000E+00 -1.4543E-01

0ITERATION NO.:   57    OBJECTIVE VALUE:  -626.977259242748        NO. OF FUNC. EVALS.:  37
 CUMULATIVE NO. OF FUNC. EVALS.:      736
 NPARAMETR:  1.0592E+00  1.1061E+00  5.6118E+12  3.9946E-01  1.3869E-01  6.7257E+05  1.1162E+01
 PARAMETER:  1.5748E-01  2.0088E-01  2.9456E+01 -8.1764E-01 -1.8755E+00  1.3519E+01  2.5125E+00
 GRADIENT:   3.3034E-02  2.1320E-02 -1.9298E-12 -4.0861E-02 -9.9620E-03  0.0000E+00  2.3504E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      736
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9135E-02 -5.2129E-02  0.0000E+00
 SE:             8.6152E-02  3.6943E-02  0.0000E+00
 N:                     100         100         100

 P VAL.:         7.3523E-01  1.5822E-01  1.0000E+00

 ETASHRINKSD(%)  8.7299E+00  6.0863E+01  1.0000E+02
 ETASHRINKVR(%)  1.6698E+01  8.4683E+01  1.0000E+02
 EBVSHRINKSD(%)  8.1508E+00  6.2042E+01  1.0000E+02
 EBVSHRINKVR(%)  1.5637E+01  8.5592E+01  1.0000E+02
 RELATIVEINF(%)  7.3366E+01  1.2530E+01  1.0000E-10
 EPSSHRINKSD(%)  6.1765E+00
 EPSSHRINKVR(%)  1.1972E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -626.97725924274755     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       108.17356732099063     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           200
  
 #TERE:
 Elapsed estimation  time in seconds:     5.29
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     1.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -626.977       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.06E+00  1.11E+00  5.61E+12  3.99E-01  1.39E-01  6.73E+05  1.12E+01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-01
 
 ETA2
+        0.00E+00  9.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.49E-01
 
 ETA2
+        0.00E+00  9.49E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.49E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        5.25E+02
 
 TH 2
+       -2.18E+02  6.14E+02
 
 TH 3
+        1.73E-16 -3.98E-16  1.22E-31
 
 TH 4
+        3.52E+01 -2.07E+02 -3.35E-16  8.69E+02
 
 TH 5
+       -6.09E+01 -2.39E+02 -2.98E-16 -8.53E+01  3.73E+02
 
 TH 6
+        1.18E-01 -5.89E-02  6.33E-25  2.98E-06  2.00E-06  9.71E-07
 
 TH 7
+       -1.03E+01 -1.43E+01  1.35E-18  8.96E+00  1.86E+01 -4.21E-03  4.15E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,        6.796
Stop Time:
Sat Sep 18 15:55:00 CDT 2021
