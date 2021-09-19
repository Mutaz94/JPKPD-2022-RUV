Sat Sep 18 09:17:32 CDT 2021
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
$DATA ../../../../data/spa/A1/dat52.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
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
 RAW OUTPUT FILE (FILE): m52.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1517.82076084265        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0843E+01 -6.4710E+01  2.1974E+01 -1.4215E+02 -3.3793E+01  7.5585E+00 -1.8048E+01  5.0002E+00 -2.8832E+01 -5.4599E+00
            -2.4542E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1556.70332774990        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0397E+00  1.0083E+00  1.0496E+00  9.7063E-01  1.0600E+00  9.9345E-01  1.0805E+00  9.5848E-01  1.1024E+00  9.9572E-01
             1.4953E+00
 PARAMETER:  1.3889E-01  1.0829E-01  1.4838E-01  7.0194E-02  1.5830E-01  9.3432E-02  1.7741E-01  5.7589E-02  1.9751E-01  9.5711E-02
             5.0232E-01
 GRADIENT:   1.1634E+02 -9.6813E+01  5.8972E-01 -1.4601E+02 -4.5780E+00  3.4878E+00 -6.7676E+00  5.6922E+00 -2.2078E-01  2.1752E+00
             9.7312E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1561.17530349056        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0263E+00  1.4345E+00  1.2577E+00  8.0921E-01  1.4356E+00  9.4026E-01  1.0386E+00  5.2746E-01  1.1256E+00  1.5244E+00
             1.4586E+00
 PARAMETER:  1.2601E-01  4.6085E-01  3.2932E-01 -1.1170E-01  4.6158E-01  3.8404E-02  1.3784E-01 -5.3968E-01  2.1828E-01  5.2158E-01
             4.7746E-01
 GRADIENT:   9.0772E+01  1.1474E+00  1.7089E+01 -3.0513E+01 -1.3524E+01 -1.6963E+01  2.4172E+00 -7.8162E-01 -3.3059E+00  1.4534E+01
            -2.6037E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1569.77610742187        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.9067E-01  1.1047E+00  8.4512E-01  1.0208E+00  9.9024E-01  9.8955E-01  1.1863E+00  5.0826E-01  9.4954E-01  8.8281E-01
             1.4739E+00
 PARAMETER:  9.0631E-02  1.9958E-01 -6.8271E-02  1.2062E-01  9.0190E-02  8.9496E-02  2.7084E-01 -5.7677E-01  4.8226E-02 -2.4649E-02
             4.8794E-01
 GRADIENT:  -6.1596E+00 -1.2872E+00 -7.7980E+00  4.3980E+00  1.4178E+01  3.9121E+00 -7.1416E+00  1.0479E+00 -3.5356E+00 -3.5299E+00
            -5.1653E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1571.45531542342        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.9008E-01  9.0683E-01  8.6216E-01  1.1446E+00  8.9114E-01  9.7656E-01  1.4394E+00  2.4453E-01  8.6389E-01  8.5910E-01
             1.4795E+00
 PARAMETER:  9.0030E-02  2.1997E-03 -4.8313E-02  2.3505E-01 -1.5249E-02  7.6283E-02  4.6426E-01 -1.3084E+00 -4.6315E-02 -5.1873E-02
             4.9171E-01
 GRADIENT:  -6.2753E+00  1.0239E+01  5.0248E+00  1.1996E+01 -5.6279E+00 -8.5331E-01 -2.3364E+00  5.9387E-02 -4.0134E+00 -1.0808E+00
             1.0007E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1572.33060943150        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.9261E-01  7.7259E-01  7.9835E-01  1.2073E+00  7.9388E-01  9.7856E-01  1.6340E+00  8.9805E-02  8.4135E-01  7.8604E-01
             1.4670E+00
 PARAMETER:  9.2585E-02 -1.5801E-01 -1.2521E-01  2.8835E-01 -1.3082E-01  7.8325E-02  5.9106E-01 -2.3101E+00 -7.2746E-02 -1.4075E-01
             4.8324E-01
 GRADIENT:   9.6280E-01  5.5250E+00  3.6347E+00  6.1891E+00 -6.3174E+00 -1.2523E-03  2.5604E-01  2.5003E-02 -4.3258E-02 -5.1243E-01
            -4.1883E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1572.74212443094        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  9.9142E-01  6.8431E-01  7.9877E-01  1.2522E+00  7.6602E-01  9.7750E-01  1.7852E+00  3.9339E-02  8.1411E-01  7.8425E-01
             1.4670E+00
 PARAMETER:  9.1386E-02 -2.7934E-01 -1.2469E-01  3.2489E-01 -1.6654E-01  7.7240E-02  6.7954E-01 -3.1355E+00 -1.0565E-01 -1.4303E-01
             4.8324E-01
 GRADIENT:   1.6298E-02  1.0297E+00 -2.6378E+00  4.9912E+00  2.3992E+00 -1.5905E-01  2.2596E-01  6.5946E-03 -1.0365E+00  2.9201E-01
             5.8733E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1573.15299299939        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:      567
 NPARAMETR:  9.9270E-01  5.7552E-01  8.2259E-01  1.3151E+00  7.3722E-01  9.7912E-01  2.0176E+00  1.2760E-02  8.0009E-01  7.8507E-01
             1.4649E+00
 PARAMETER:  9.2670E-02 -4.5248E-01 -9.5299E-02  3.7392E-01 -2.0487E-01  7.8895E-02  8.0189E-01 -4.2615E+00 -1.2303E-01 -1.4199E-01
             4.8177E-01
 GRADIENT:  -1.0776E+01  1.8847E+00  1.1670E+01 -1.1853E+01 -1.5756E+01 -3.3675E-01 -1.1164E+00  1.7473E-04  1.7131E-02 -1.1935E+00
            -1.6493E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1573.48212170991        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      742
 NPARAMETR:  9.9512E-01  4.7327E-01  8.9057E-01  1.3862E+00  7.5493E-01  9.7836E-01  2.3299E+00  1.0000E-02  7.7862E-01  8.4879E-01
             1.4704E+00
 PARAMETER:  9.5107E-02 -6.4809E-01 -1.5894E-02  4.2654E-01 -1.8113E-01  7.8124E-02  9.4584E-01 -5.6727E+00 -1.5023E-01 -6.3942E-02
             4.8552E-01
 GRADIENT:  -2.0635E-01  7.5550E-01  6.7779E-01  2.4747E+00 -5.2585E-01  1.9209E-01  4.4014E-02  0.0000E+00 -1.5007E-01 -3.0791E-01
            -3.6011E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1573.48631952356        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      904
 NPARAMETR:  9.9508E-01  4.6395E-01  8.8609E-01  1.3891E+00  7.4990E-01  9.7773E-01  2.3575E+00  1.0000E-02  7.7746E-01  8.4759E-01
             1.4703E+00
 PARAMETER:  9.5070E-02 -6.6799E-01 -2.0931E-02  4.2867E-01 -1.8782E-01  7.7481E-02  9.5759E-01 -5.8689E+00 -1.5173E-01 -6.5356E-02
             4.8549E-01
 GRADIENT:  -5.4772E-04  1.1353E-03  5.6442E-04  5.5594E-03  1.8162E-03 -3.0653E-03 -1.8000E-04  0.0000E+00  5.9699E-04 -5.6197E-04
            -1.7445E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      904
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.7101E-04  3.0818E-02 -3.3499E-04 -2.9424E-02 -4.9779E-03
 SE:             2.9653E-02  1.9543E-02  1.9052E-04  2.3890E-02  2.0981E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7926E-01  1.1482E-01  7.8691E-02  2.1808E-01  8.1246E-01

 ETASHRINKSD(%)  6.5872E-01  3.4528E+01  9.9362E+01  1.9966E+01  2.9710E+01
 ETASHRINKVR(%)  1.3131E+00  5.7134E+01  9.9996E+01  3.5945E+01  5.0593E+01
 EBVSHRINKSD(%)  9.5469E-01  3.8617E+01  9.9344E+01  1.7530E+01  2.6392E+01
 EBVSHRINKVR(%)  1.9003E+00  6.2321E+01  9.9996E+01  3.1987E+01  4.5819E+01
 RELATIVEINF(%)  9.7296E+01  6.1065E+00  3.7173E-04  1.3901E+01  4.1527E+00
 EPSSHRINKSD(%)  3.8994E+01
 EPSSHRINKVR(%)  6.2782E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1573.4863195235553     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -838.33549295981709     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1573.486       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  4.64E-01  8.86E-01  1.39E+00  7.50E-01  9.78E-01  2.36E+00  1.00E-02  7.77E-01  8.48E-01  1.47E+00
 


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
+        1.16E+03
 
 TH 2
+       -2.28E+01  3.84E+02
 
 TH 3
+        1.08E+01  1.81E+02  6.50E+02
 
 TH 4
+       -1.23E+01  2.89E+02 -1.74E+02  6.68E+02
 
 TH 5
+       -6.68E+00 -3.60E+02 -9.15E+02  1.60E+02  1.49E+03
 
 TH 6
+        3.78E-01 -4.46E+00 -9.96E-01 -5.53E+00 -2.59E-01  2.04E+02
 
 TH 7
+        1.42E+00  3.26E+01  9.05E-01 -6.50E+00 -3.24E+00 -1.03E-01  1.12E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.12E-01 -1.90E+01 -1.43E+01 -1.09E+01  2.15E+01  2.07E+00  7.26E+00  0.00E+00  1.77E+02
 
 TH10
+        9.46E-01  8.15E+00 -3.54E+01 -2.70E+01 -4.15E+01  2.10E+00  2.72E+00  0.00E+00  2.68E+00  9.12E+01
 
 TH11
+       -1.08E+01 -6.46E+00 -3.60E+01 -1.18E+01  1.27E+01  4.27E+00  1.90E+00  0.00E+00  1.05E+01  2.44E+01  1.12E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.093
Stop Time:
Sat Sep 18 09:17:49 CDT 2021
