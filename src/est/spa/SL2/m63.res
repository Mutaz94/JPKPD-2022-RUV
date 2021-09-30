Wed Sep 29 16:01:44 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/SL2/dat63.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
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
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
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

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
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
 NO. OF SIG. FIGURES REQUIRED:            2
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
 RAW OUTPUT FILE (FILE): m63.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1687.78362849635        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2558E+02 -8.3974E+01 -4.7115E+01 -6.7891E+01  4.6444E+01  1.8704E+01 -1.0345E+01  1.4152E+01 -1.0505E+01  2.2234E+01
             9.3034E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1701.99521897359        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0352E+00  1.1095E+00  1.2281E+00  1.0333E+00  1.0816E+00  1.1132E+00  1.0790E+00  8.9958E-01  1.0742E+00  8.5172E-01
             9.9271E-01
 PARAMETER:  1.3458E-01  2.0388E-01  3.0549E-01  1.3279E-01  1.7843E-01  2.0727E-01  1.7604E-01 -5.8265E-03  1.7156E-01 -6.0499E-02
             9.2681E-02
 GRADIENT:  -1.7574E+01  4.4671E+00  2.4182E+01 -1.3832E+01 -3.9313E+00  1.7894E+01 -7.9988E-01 -5.6839E+00  2.0509E+00 -1.6131E+01
            -5.3760E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1703.65301983132        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0354E+00  9.8324E-01  1.0378E+00  1.1192E+00  9.7428E-01  1.1050E+00  1.3808E+00  5.3783E-01  9.3281E-01  8.4246E-01
             1.0009E+00
 PARAMETER:  1.3479E-01  8.3096E-02  1.3708E-01  2.1259E-01  7.3946E-02  1.9982E-01  4.2266E-01 -5.2021E-01  3.0446E-02 -7.1432E-02
             1.0086E-01
 GRADIENT:  -1.8621E+01  8.2910E+00 -1.6308E+00  1.4367E+01  2.6549E+01  1.4770E+01  7.1266E+00 -1.8730E+00 -2.9686E+00 -4.0850E+00
             7.9370E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1705.15716193080        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0465E+00  9.7030E-01  9.3162E-01  1.1107E+00  9.0501E-01  1.0590E+00  1.3047E+00  4.9998E-01  9.4373E-01  7.9168E-01
             9.9230E-01
 PARAMETER:  1.4541E-01  6.9854E-02  2.9171E-02  2.0501E-01  1.8711E-04  1.5729E-01  3.6596E-01 -5.9319E-01  4.2086E-02 -1.3360E-01
             9.2269E-02
 GRADIENT:  -5.8413E-01 -5.6632E-02 -4.1993E+00  4.4900E+00  7.6322E+00 -1.4325E+00 -7.6684E-01  5.8060E-01 -2.1122E-03  2.2794E-02
            -2.3925E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1705.49773817056        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      728
 NPARAMETR:  1.0466E+00  7.8423E-01  8.4945E-01  1.2064E+00  7.8550E-01  1.0624E+00  1.5985E+00  3.0903E-01  8.6296E-01  7.0628E-01
             9.9516E-01
 PARAMETER:  1.4552E-01 -1.4305E-01 -6.3165E-02  2.8763E-01 -1.4143E-01  1.6049E-01  5.6907E-01 -1.0743E+00 -4.7383E-02 -2.4774E-01
             9.5152E-02
 GRADIENT:   4.0758E-01  8.1122E-01 -2.6589E-02 -6.2504E-03 -1.1661E+00 -1.0807E-01  7.1833E-01  1.9332E-01  3.2291E-01  5.7779E-01
             2.1275E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1705.49944373981        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  1.0464E+00  7.7025E-01  8.4222E-01  1.2125E+00  7.7738E-01  1.0624E+00  1.6164E+00  2.8470E-01  8.5742E-01  7.0017E-01
             9.9554E-01
 PARAMETER:  1.4535E-01 -1.6104E-01 -7.1712E-02  2.9272E-01 -1.5182E-01  1.6052E-01  5.8021E-01 -1.1563E+00 -5.3826E-02 -2.5643E-01
             9.5530E-02
 GRADIENT:   1.1418E-01 -3.9530E-01 -5.1696E-01 -2.5912E-01  4.3172E-01 -8.9037E-02  1.1469E-01  1.0250E-01  7.3408E-02  2.1400E-01
             9.3612E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1705.57275017463        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1082
 NPARAMETR:  1.0472E+00  8.3351E-01  7.9793E-01  1.1719E+00  7.7448E-01  1.0641E+00  1.5186E+00  1.1257E-01  8.6947E-01  6.8823E-01
             9.9578E-01
 PARAMETER:  1.4612E-01 -8.2109E-02 -1.2573E-01  2.5863E-01 -1.5556E-01  1.6216E-01  5.1781E-01 -2.0841E+00 -3.9873E-02 -2.7363E-01
             9.5771E-02
 GRADIENT:  -1.4677E-01  1.7509E+00  3.7795E+00  1.9068E-01 -5.0567E+00  2.7695E-01 -8.8941E-01 -5.0374E-03 -6.4864E-01 -4.0615E-01
            -2.8416E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1705.58807344608        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1241
 NPARAMETR:  1.0483E+00  8.2723E-01  7.9020E-01  1.1727E+00  7.6964E-01  1.0636E+00  1.5372E+00  7.2339E-02  8.6839E-01  6.8297E-01
             9.9575E-01
 PARAMETER:  1.4720E-01 -8.9678E-02 -1.3547E-01  2.5929E-01 -1.6184E-01  1.6168E-01  5.2996E-01 -2.5264E+00 -4.1111E-02 -2.8130E-01
             9.5744E-02
 GRADIENT:   1.9352E+00 -2.4316E-01  4.7140E-01 -9.1779E-01 -3.3033E-01  6.2841E-02  4.5729E-02  2.3313E-03  1.5579E-02 -1.3433E-01
            -1.3191E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1705.59230736963        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1431             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0491E+00  8.3086E-01  7.8793E-01  1.1704E+00  7.6999E-01  1.0641E+00  1.5324E+00  1.0000E-02  8.6901E-01  6.8475E-01
             9.9595E-01
 PARAMETER:  1.4792E-01 -8.5288E-02 -1.3835E-01  2.5735E-01 -1.6138E-01  1.6212E-01  5.2686E-01 -4.7276E+00 -4.0404E-02 -2.7870E-01
             9.5938E-02
 GRADIENT:   6.1600E+02  2.2120E+01  4.7379E+00  2.3222E+02  1.6939E+01  6.5228E+01  2.4385E+01  0.0000E+00  6.7082E+00  1.3798E+00
             8.1079E-01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1705.59230736963        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1488
 NPARAMETR:  1.0491E+00  8.3086E-01  7.8793E-01  1.1704E+00  7.6999E-01  1.0641E+00  1.5324E+00  1.0000E-02  8.6901E-01  6.8475E-01
             9.9595E-01
 PARAMETER:  1.4792E-01 -8.5288E-02 -1.3835E-01  2.5735E-01 -1.6138E-01  1.6212E-01  5.2686E-01 -4.7276E+00 -4.0404E-02 -2.7870E-01
             9.5938E-02
 GRADIENT:   3.2797E+00 -3.7664E-01 -5.6434E-02 -7.2530E-01  2.5597E-01  2.2564E-01  1.3839E-01  0.0000E+00 -1.5818E-03  1.1963E-01
             1.0805E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1488
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1972E-04  1.1895E-02 -5.2660E-04 -1.1842E-02 -2.6709E-03
 SE:             2.9853E-02  2.2134E-02  2.2492E-04  2.5010E-02  2.1188E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9680E-01  5.9099E-01  1.9219E-02  6.3585E-01  8.9969E-01

 ETASHRINKSD(%)  1.0000E-10  2.5848E+01  9.9246E+01  1.6215E+01  2.9019E+01
 ETASHRINKVR(%)  1.0000E-10  4.5015E+01  9.9994E+01  2.9801E+01  4.9616E+01
 EBVSHRINKSD(%)  3.8223E-01  2.5739E+01  9.9300E+01  1.5992E+01  2.8363E+01
 EBVSHRINKVR(%)  7.6299E-01  4.4853E+01  9.9995E+01  2.9427E+01  4.8682E+01
 RELATIVEINF(%)  9.8951E+01  5.7373E+00  5.3180E-04  9.6030E+00  3.9052E+00
 EPSSHRINKSD(%)  4.3334E+01
 EPSSHRINKVR(%)  6.7889E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1705.5923073696340     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -970.44148080589582     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1705.592       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  8.31E-01  7.88E-01  1.17E+00  7.70E-01  1.06E+00  1.53E+00  1.00E-02  8.69E-01  6.85E-01  9.96E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


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
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        8.85E+02
 
 TH 2
+       -6.33E+00  3.94E+02
 
 TH 3
+        1.61E+01  2.88E+02  9.42E+02
 
 TH 4
+       -7.23E+00  2.49E+02 -3.41E+02  7.71E+02
 
 TH 5
+       -5.87E+00 -5.15E+02 -1.24E+03  4.23E+02  2.04E+03
 
 TH 6
+        3.17E-01 -1.31E+00  2.81E+00 -2.10E+00 -3.14E+00  1.73E+02
 
 TH 7
+        9.68E-01  3.34E+01 -9.27E+00 -1.09E+01  1.71E+00 -3.42E-01  3.13E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.11E+00 -2.03E+01 -3.91E+01  2.50E+01  4.19E+00 -8.03E-01  1.15E+01  0.00E+00  1.42E+02
 
 TH10
+       -1.84E+00 -8.42E+00 -8.87E+01 -3.52E+01 -3.49E+01  1.16E+00  1.43E+01  0.00E+00  1.29E+01  1.19E+02
 
 TH11
+       -5.88E+00 -9.81E+00 -4.53E+01 -8.63E+00  1.42E+01  2.35E+00  3.82E+00  0.00E+00  1.12E+01  3.00E+01  2.20E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       25.672
Stop Time:
Wed Sep 29 16:02:11 CDT 2021
