Sat Sep 25 12:03:34 CDT 2021
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
$DATA ../../../../data/spa/S2/dat3.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1690.99162470345        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -8.0088E+00 -7.3715E+01 -2.6163E+01 -7.3581E+01  3.4240E+01  4.4760E+00  3.4948E+00  2.9164E+00  1.5562E+00  2.4429E+01
            -2.5972E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1697.75059398196        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0198E+00  1.1039E+00  1.0539E+00  1.0196E+00  1.0154E+00  9.8804E-01  9.1119E-01  1.0372E+00  9.8875E-01  7.0269E-01
             1.1263E+00
 PARAMETER:  1.1956E-01  1.9885E-01  1.5250E-01  1.1936E-01  1.1531E-01  8.7968E-02  6.9989E-03  1.3652E-01  8.8688E-02 -2.5284E-01
             2.1891E-01
 GRADIENT:   2.8145E+01  2.5936E+01  7.8415E+00  3.9738E+01  2.1435E+01 -1.2944E+00 -3.4937E+00 -1.0968E+01 -6.3017E+00 -1.1893E+01
             1.0968E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1698.76742875111        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0232E+00  1.1595E+00  9.9858E-01  9.8297E-01  9.9718E-01  9.8452E-01  9.5022E-01  1.2668E+00  1.0161E+00  6.3942E-01
             1.1296E+00
 PARAMETER:  1.2298E-01  2.4801E-01  9.8584E-02  8.2823E-02  9.7173E-02  8.4401E-02  4.8943E-02  3.3647E-01  1.1593E-01 -3.4720E-01
             2.2188E-01
 GRADIENT:   3.5625E+01  3.3492E+01  3.5678E+00  3.6696E+01 -2.8688E+00 -3.2529E+00  2.5513E+00 -1.0147E+00 -6.0714E-01 -7.3807E+00
             1.5251E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.04253628296        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  1.0081E+00  1.1429E+00  1.0783E+00  9.7234E-01  1.0352E+00  9.9225E-01  8.5291E-01  1.3343E+00  1.0444E+00  7.6438E-01
             1.0790E+00
 PARAMETER:  1.0807E-01  2.3355E-01  1.7538E-01  7.1954E-02  1.3457E-01  9.2218E-02 -5.9099E-02  3.8843E-01  1.4342E-01 -1.6869E-01
             1.7600E-01
 GRADIENT:  -3.8577E+01 -5.4766E+00  6.4969E-01 -1.5199E+00 -2.0352E-01 -6.0265E+00 -8.0301E-01 -3.0124E-01 -1.0739E+00 -1.3274E+00
             1.2214E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1702.46740796572        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      431
 NPARAMETR:  1.0289E+00  1.6701E+00  6.4097E-01  6.4396E-01  1.0941E+00  1.0132E+00  7.1318E-01  1.1164E+00  1.4019E+00  7.8543E-01
             1.0782E+00
 PARAMETER:  1.2848E-01  6.1287E-01 -3.4478E-01 -3.4012E-01  1.8990E-01  1.1316E-01 -2.3802E-01  2.1007E-01  4.3785E-01 -1.4152E-01
             1.7525E-01
 GRADIENT:   5.9694E-01  5.5481E+01  9.5564E+00  2.0762E+01 -2.5959E+01  1.2938E+00 -1.9096E+00 -8.8561E-01 -1.8189E+00  1.1411E+00
             3.9056E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1704.49298150757        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  1.0263E+00  1.8346E+00  4.3259E-01  5.0389E-01  1.1092E+00  1.0063E+00  6.8953E-01  8.3678E-01  1.5960E+00  7.4974E-01
             1.0753E+00
 PARAMETER:  1.2598E-01  7.0683E-01 -7.3797E-01 -5.8540E-01  2.0368E-01  1.0630E-01 -2.7174E-01 -7.8188E-02  5.6753E-01 -1.8803E-01
             1.7262E-01
 GRADIENT:  -4.8707E+00  3.0162E+00 -1.4589E+00  6.8093E+00  1.5721E+00 -1.4956E+00 -5.8141E-01  8.4220E-01 -5.8549E-01 -8.8709E-02
             3.8858E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1704.62162274163        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  1.0267E+00  1.8558E+00  4.1557E-01  4.8947E-01  1.1134E+00  1.0068E+00  6.8530E-01  8.0931E-01  1.6230E+00  7.4951E-01
             1.0754E+00
 PARAMETER:  1.2632E-01  7.1832E-01 -7.7810E-01 -6.1443E-01  2.0740E-01  1.0679E-01 -2.7790E-01 -1.1157E-01  5.8426E-01 -1.8833E-01
             1.7268E-01
 GRADIENT:  -4.1008E+00  3.8914E+00 -1.5867E+00  6.9835E+00  1.7394E+00 -1.2922E+00 -6.5342E-01  8.3163E-01 -6.1609E-01 -1.1374E-01
             1.2385E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1704.62879430370        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:      971
 NPARAMETR:  1.0286E+00  1.8552E+00  4.1534E-01  4.8926E-01  1.1133E+00  1.0069E+00  6.8601E-01  8.0923E-01  1.6226E+00  7.4735E-01
             1.0753E+00
 PARAMETER:  1.2817E-01  7.1800E-01 -7.7865E-01 -6.1486E-01  2.0732E-01  1.0684E-01 -2.7686E-01 -1.1167E-01  5.8405E-01 -1.9122E-01
             1.7259E-01
 GRADIENT:  -5.5956E-02  2.6281E+00 -1.6345E+00  6.6182E+00  2.3683E+00 -1.2622E+00 -5.3756E-01  8.3168E-01 -6.2736E-01 -3.4282E-01
            -1.1358E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1704.63226214972        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1144
 NPARAMETR:  1.0285E+00  1.8551E+00  4.1538E-01  4.8916E-01  1.1133E+00  1.0069E+00  6.8774E-01  8.0866E-01  1.6226E+00  7.4940E-01
             1.0753E+00
 PARAMETER:  1.2812E-01  7.1792E-01 -7.7856E-01 -6.1507E-01  2.0730E-01  1.0683E-01 -2.7435E-01 -1.1238E-01  5.8403E-01 -1.8849E-01
             1.7259E-01
 GRADIENT:  -1.6422E-01  2.2688E+00 -1.6609E+00  6.3796E+00  1.9124E+00 -1.2651E+00 -1.5883E-03  8.5618E-01 -5.3241E-01 -6.3442E-03
             2.7336E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1704.63226214972        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1173
 NPARAMETR:  1.0285E+00  1.8551E+00  4.1538E-01  4.8916E-01  1.1133E+00  1.0069E+00  6.8774E-01  8.0866E-01  1.6226E+00  7.4940E-01
             1.0753E+00
 PARAMETER:  1.2812E-01  7.1792E-01 -7.7856E-01 -6.1507E-01  2.0730E-01  1.0683E-01 -2.7435E-01 -1.1238E-01  5.8403E-01 -1.8849E-01
             1.7259E-01
 GRADIENT:  -7.2344E+05  1.2910E+05  2.3812E+05  1.5070E+05  4.4712E+05 -8.6764E+05 -1.6853E-03 -1.6495E+06  3.1738E+05 -6.4940E-03
             1.0738E+06

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1173
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.6277E-04 -3.1791E-02 -1.5175E-02  2.3968E-02 -4.1166E-02
 SE:             2.9953E-02  2.4124E-02  6.8732E-03  2.3166E-02  1.9679E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8501E-01  1.8756E-01  2.7251E-02  3.0085E-01  3.6448E-02

 ETASHRINKSD(%)  1.0000E-10  1.9183E+01  7.6974E+01  2.2391E+01  3.4074E+01
 ETASHRINKVR(%)  1.0000E-10  3.4686E+01  9.4698E+01  3.9769E+01  5.6537E+01
 EBVSHRINKSD(%)  4.7490E-01  1.9233E+01  7.7991E+01  2.2498E+01  3.3954E+01
 EBVSHRINKVR(%)  9.4755E-01  3.4767E+01  9.5156E+01  3.9934E+01  5.6379E+01
 RELATIVEINF(%)  9.8992E+01  4.8380E+00  4.2361E-01  4.5798E+00  9.1482E+00
 EPSSHRINKSD(%)  4.3831E+01
 EPSSHRINKVR(%)  6.8451E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1704.6322621497238     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -969.48143558598565     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.86
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.72
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1704.632       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.86E+00  4.15E-01  4.89E-01  1.11E+00  1.01E+00  6.88E-01  8.09E-01  1.62E+00  7.49E-01  1.08E+00
 


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
+        2.67E+09
 
 TH 2
+        1.16E+03  2.61E+07
 
 TH 3
+        3.87E+03 -6.13E+03  4.43E+08
 
 TH 4
+        4.26E+03 -6.56E+03 -2.51E+04  5.12E+08
 
 TH 5
+        7.94E+03  3.24E+02  3.26E+03  4.16E+03  8.70E+08
 
 TH 6
+        3.27E+09  2.32E+03  7.76E+03  8.54E+03 -1.87E+09  4.01E+09
 
 TH 7
+       -2.38E+03  5.01E+02  2.16E+03  2.30E+03  2.69E+03 -5.48E+03  1.97E+02
 
 TH 8
+       -1.44E+04  2.28E+04 -8.82E+04 -9.18E+04 -1.25E+04 -2.87E+04 -7.62E+03  5.61E+09
 
 TH 9
+        1.28E+03 -2.17E+03  3.22E+04 -8.07E+03  1.30E+03  2.57E+03  7.64E+02 -3.08E+04  5.16E+07
 
 TH10
+       -3.48E+03  4.04E+02  1.76E+03  1.88E+03  2.25E+03 -4.87E+03  2.66E+01 -6.27E+03  6.16E+02  8.90E+01
 
 TH11
+        8.58E+03  1.87E+08 -4.39E+05 -4.72E+05  2.84E+03  1.72E+04  3.50E+03  1.56E+06 -1.50E+05  2.98E+03  1.34E+09
 
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
 #CPUT: Total CPU Time in Seconds,       21.637
Stop Time:
Sat Sep 25 12:03:57 CDT 2021
