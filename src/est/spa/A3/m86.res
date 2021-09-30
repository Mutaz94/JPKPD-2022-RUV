Wed Sep 29 13:52:12 CDT 2021
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
$DATA ../../../../data/spa/A3/dat86.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -126.660851266696        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2143E+02  6.9751E+01  6.8275E+01  2.8864E+01  2.2519E+02  5.8116E+01 -4.1078E+01 -5.4893E+01 -1.4200E+02 -1.3993E+02
            -2.6909E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1251.09081349242        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0225E+00  9.1639E-01  9.0603E-01  1.0680E+00  8.5364E-01  8.4498E-01  9.9966E-01  1.0110E+00  1.2017E+00  1.0160E+00
             2.9078E+00
 PARAMETER:  1.2227E-01  1.2689E-02  1.3154E-03  1.6577E-01 -5.8249E-02 -6.8439E-02  9.9657E-02  1.1095E-01  2.8376E-01  1.1583E-01
             1.1674E+00
 GRADIENT:   5.6834E+01 -3.3712E+00  5.4996E+00 -1.2033E+01  4.2382E+01 -1.9001E+01  1.1356E-01 -5.3646E-01 -1.3423E+01 -8.0388E+00
            -1.3154E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1262.47033059532        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:      166
 NPARAMETR:  1.0271E+00  5.4231E-01  4.3693E-01  1.3074E+00  4.2603E-01  9.5565E-01  6.6368E-01  7.3480E-01  1.3843E+00  4.9968E-01
             2.7582E+00
 PARAMETER:  1.2673E-01 -5.1192E-01 -7.2798E-01  3.6801E-01 -7.5324E-01  5.4632E-02 -3.0995E-01 -2.0815E-01  4.2519E-01 -5.9379E-01
             1.1146E+00
 GRADIENT:   5.4146E+01  3.8860E+01  3.8612E+01  1.3015E+02  6.2559E+00  1.5052E+01 -1.1547E+00 -2.7042E+00  3.8457E+01 -1.0008E+01
            -1.4421E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1281.44372654958        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      326
 NPARAMETR:  1.0338E+00  8.0204E-01  8.3902E-01  1.1630E+00  7.5193E-01  8.9333E-01  2.2960E-01  7.6358E-01  1.2181E+00  3.8668E-01
             3.6978E+00
 PARAMETER:  1.3325E-01 -1.2060E-01 -7.5517E-02  2.5102E-01 -1.8512E-01 -1.2800E-02 -1.3714E+00 -1.6973E-01  2.9729E-01 -8.5015E-01
             1.4077E+00
 GRADIENT:  -2.9786E+00  1.6752E+00  1.8706E+01 -4.8506E+00 -1.0559E+01  2.1057E+00 -8.1151E-02  4.3012E-01  9.2733E+00  2.2940E+00
            -2.2829E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1286.60445298061        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  1.0440E+00  6.5808E-01  5.0034E-01  1.2135E+00  5.1471E-01  9.0216E-01  9.4086E-01  3.8133E-01  1.0834E+00  1.8277E-01
             3.7204E+00
 PARAMETER:  1.4308E-01 -3.1842E-01 -5.9247E-01  2.9348E-01 -5.6414E-01 -2.9581E-03  3.9042E-02 -8.6410E-01  1.8007E-01 -1.5995E+00
             1.4138E+00
 GRADIENT:   2.9172E+00  7.1253E+00  2.5284E+00  1.8029E+01 -5.5565E+00 -9.7261E-01  1.1550E+00 -1.3712E-01  1.8357E+00  3.4227E-01
            -3.0383E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1286.95275507695        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      678
 NPARAMETR:  1.0423E+00  7.3187E-01  4.4921E-01  1.1418E+00  5.0725E-01  9.0741E-01  7.9197E-01  2.5833E-01  1.0992E+00  1.3636E-01
             3.7413E+00
 PARAMETER:  1.4142E-01 -2.1216E-01 -7.0026E-01  2.3261E-01 -5.7875E-01  2.8420E-03 -1.3323E-01 -1.2535E+00  1.9459E-01 -1.8924E+00
             1.4194E+00
 GRADIENT:  -2.5621E-01  5.1043E+00  5.3903E+00  2.5142E+00 -9.9940E+00 -1.2342E-01 -5.2269E-01 -1.8274E-01 -2.0118E+00  1.1311E-01
            -1.3209E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1287.13430560412        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      855
 NPARAMETR:  1.0420E+00  8.8434E-01  3.6699E-01  1.0260E+00  5.0672E-01  9.0732E-01  7.5328E-01  1.0561E-01  1.1901E+00  7.2829E-02
             3.6858E+00
 PARAMETER:  1.4111E-01 -2.2914E-02 -9.0242E-01  1.2562E-01 -5.7979E-01  2.7402E-03 -1.8332E-01 -2.1480E+00  2.7401E-01 -2.5196E+00
             1.4045E+00
 GRADIENT:   4.2959E+00  8.3931E+00  8.9850E+00 -1.9411E+00 -1.4707E+01 -1.7245E+00 -4.9554E-01 -4.4186E-02  5.6257E-01  4.5099E-02
            -1.6721E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1288.50708554174        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1033
 NPARAMETR:  1.0322E+00  7.9048E-01  2.0030E-01  9.4371E-01  3.6258E-01  9.3531E-01  9.4768E-01  2.4440E-02  1.1678E+00  1.5390E-02
             3.5172E+00
 PARAMETER:  1.3167E-01 -1.3511E-01 -1.5080E+00  4.2067E-02 -9.1451E-01  3.3127E-02  4.6264E-02 -3.6115E+00  2.5511E-01 -4.0740E+00
             1.3577E+00
 GRADIENT:   7.9459E+00  1.5607E+01  3.4457E+00 -3.5541E+00 -1.4554E+01  2.2826E+00  8.2231E+00 -4.5098E-03 -1.6986E+01  1.4013E-03
             2.0760E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1291.24251297092        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1212
 NPARAMETR:  1.0173E+00  7.6810E-01  1.7082E-01  9.1699E-01  3.4000E-01  9.3604E-01  8.5599E-01  1.4153E-02  1.3638E+00  1.1204E-02
             3.3196E+00
 PARAMETER:  1.1714E-01 -1.6383E-01 -1.6672E+00  1.3342E-02 -9.7882E-01  3.3899E-02 -5.5498E-02 -4.1578E+00  4.1030E-01 -4.3915E+00
             1.2998E+00
 GRADIENT:  -5.7999E-01  2.6609E+00  6.0588E-04  4.1870E-01 -2.2016E+00 -1.2268E-01  7.5530E-02 -2.1501E-03 -6.4668E-01 -4.3318E-03
             3.6676E-01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1291.26042624521        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1339
 NPARAMETR:  1.0181E+00  7.4971E-01  1.7183E-01  9.2413E-01  3.3516E-01  9.3711E-01  8.7061E-01  1.5074E-02  1.3629E+00  1.1624E-02
             3.3131E+00
 PARAMETER:  1.1793E-01 -1.8807E-01 -1.6612E+00  2.1093E-02 -9.9314E-01  3.5043E-02 -3.8566E-02 -4.0948E+00  4.0958E-01 -4.3547E+00
             1.2979E+00
 GRADIENT:   2.6768E-02 -6.9978E-02 -4.6370E-02 -1.1207E-01 -2.5327E-02  1.2178E-02  1.5096E-02 -2.5701E-03  5.4643E-04 -5.1555E-03
             4.1569E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1339
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.3717E-04 -1.5650E-03 -9.1768E-05 -1.1234E-02  2.3564E-04
 SE:             2.8683E-02  2.0231E-02  2.3217E-04  2.5638E-02  4.2357E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7394E-01  9.3834E-01  6.9265E-01  6.6125E-01  5.7799E-01

 ETASHRINKSD(%)  3.9077E+00  3.2222E+01  9.9222E+01  1.4110E+01  9.8581E+01
 ETASHRINKVR(%)  7.6627E+00  5.4061E+01  9.9994E+01  2.6229E+01  9.9980E+01
 EBVSHRINKSD(%)  3.7811E+00  3.1584E+01  9.9315E+01  1.2296E+01  9.8714E+01
 EBVSHRINKVR(%)  7.4193E+00  5.3192E+01  9.9995E+01  2.3080E+01  9.9983E+01
 RELATIVEINF(%)  8.6424E+01  1.9004E+00  9.9303E-04  2.8339E+01  5.4364E-04
 EPSSHRINKSD(%)  2.9029E+01
 EPSSHRINKVR(%)  4.9631E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1291.2604262452130     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -556.10959968147483     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1291.260       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  7.50E-01  1.72E-01  9.24E-01  3.35E-01  9.37E-01  8.71E-01  1.51E-02  1.36E+00  1.16E-02  3.31E+00
 


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
+        1.13E+03
 
 TH 2
+       -8.42E+01  1.20E+03
 
 TH 3
+       -4.18E+02  2.19E+03  8.23E+03
 
 TH 4
+       -3.76E+01  1.02E+02 -8.70E+02  5.40E+02
 
 TH 5
+        3.20E+02 -3.65E+03 -8.36E+03  2.91E+02  1.25E+04
 
 TH 6
+       -1.32E+00 -3.41E+00  2.91E+01 -1.51E+01  3.90E+01  1.91E+02
 
 TH 7
+       -3.84E+00 -1.15E+01 -1.46E+02  2.56E+00  1.64E+02  5.96E+00  5.02E+01
 
 TH 8
+       -2.44E-02 -3.37E-01 -1.97E+00 -3.68E-02  2.32E+00 -1.95E-02  2.10E-01 -1.14E+01
 
 TH 9
+        1.02E+01 -1.72E+01  1.05E+02 -5.62E+00  2.87E+01  1.72E+00  5.14E+00  6.97E-02  5.78E+01
 
 TH10
+       -7.34E-02 -1.44E+00 -8.97E+00 -7.96E-02  1.12E+01  1.31E-01  1.31E+00 -1.04E-02 -6.38E-02 -3.88E+01
 
 TH11
+       -1.73E+01 -1.11E+01 -6.02E+01 -6.11E+00  3.13E+01  2.35E+00  1.36E+01  1.87E-01  5.84E+00  7.74E-01  2.96E+01
 
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
 #CPUT: Total CPU Time in Seconds,       23.182
Stop Time:
Wed Sep 29 13:52:37 CDT 2021
