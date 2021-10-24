Sun Oct 24 03:27:26 CDT 2021
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
$DATA ../../../../data/SD4/SL3/dat15.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1667.15998485480        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2213E+02 -9.1314E+01 -5.6161E+01 -3.1495E+01  1.6519E+02  6.8338E+01 -1.2847E+01 -8.5898E+00 -2.5822E+01 -1.9632E+01
             3.0117E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1679.25066657251        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  1.0916E+00  1.0964E+00  1.0158E+00  1.0830E+00  8.5223E-01  8.7587E-01  1.0608E+00  1.0684E+00  1.1499E+00  1.0374E+00
             9.9443E-01
 PARAMETER:  1.8765E-01  1.9200E-01  1.1570E-01  1.7978E-01 -5.9897E-02 -3.2540E-02  1.5903E-01  1.6614E-01  2.3970E-01  1.3672E-01
             9.4412E-02
 GRADIENT:   1.1769E+02  7.1379E+01  3.9232E+01  5.7993E+01 -9.7792E+01 -1.4505E+01  1.9356E+00 -8.0183E+00  9.0412E+00  1.2119E+01
             9.2187E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1682.09372705833        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      347
 NPARAMETR:  1.0770E+00  1.2182E+00  1.3903E+00  1.0263E+00  1.0185E+00  9.0719E-01  9.9232E-01  2.3415E+00  1.2497E+00  9.7010E-01
             1.0072E+00
 PARAMETER:  1.7422E-01  2.9735E-01  4.2949E-01  1.2592E-01  1.1838E-01  2.5994E-03  9.2290E-02  9.5078E-01  3.2290E-01  6.9644E-02
             1.0720E-01
 GRADIENT:   7.5902E+01  5.1615E+01  2.9994E+00  6.9287E+01 -3.9218E+01  2.3770E+00  9.1999E+00  1.5732E+01  1.3772E+01  1.5081E+00
             7.0775E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1689.83080058247        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      524
 NPARAMETR:  1.0485E+00  1.2894E+00  1.1541E+00  9.0232E-01  1.0286E+00  8.9644E-01  7.4533E-01  1.7269E+00  1.2664E+00  1.0238E+00
             9.9080E-01
 PARAMETER:  1.4736E-01  3.5414E-01  2.4334E-01 -2.7859E-03  1.2817E-01 -9.3229E-03 -1.9392E-01  6.4634E-01  3.3615E-01  1.2351E-01
             9.0762E-02
 GRADIENT:   1.6034E+00  4.3166E+00  4.6386E+00  5.8390E-01 -8.1987E+00 -5.3081E-01 -3.1283E+00 -1.4628E+00 -5.5624E+00 -1.4779E+00
            -2.7979E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1691.46729858633        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  1.0497E+00  1.6872E+00  6.5557E-01  6.4127E-01  1.0338E+00  9.0075E-01  7.8689E-01  1.4641E+00  1.5203E+00  9.7856E-01
             9.9163E-01
 PARAMETER:  1.4848E-01  6.2310E-01 -3.2225E-01 -3.4431E-01  1.3323E-01 -4.5222E-03 -1.3967E-01  4.8121E-01  5.1891E-01  7.8325E-02
             9.1599E-02
 GRADIENT:  -4.8411E-01  1.4202E+01  4.2211E-01  1.1095E+01 -6.2405E+00  3.8670E-01 -1.0860E+00  1.6576E-02 -3.5678E+00 -8.8929E-02
             9.0968E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1692.01663404677        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  1.0483E+00  1.9342E+00  4.6355E-01  4.7662E-01  1.0872E+00  9.0111E-01  7.5397E-01  1.2450E+00  1.8305E+00  1.0094E+00
             9.9009E-01
 PARAMETER:  1.4715E-01  7.5967E-01 -6.6884E-01 -6.4103E-01  1.8359E-01 -4.1235E-03 -1.8240E-01  3.1914E-01  7.0461E-01  1.0940E-01
             9.0045E-02
 GRADIENT:  -4.1125E+00  2.3884E+01  1.5663E+00  1.0078E+01 -6.6283E+00  2.7172E-01 -3.9378E-02 -9.4870E-01 -2.2405E+00 -2.7021E-01
            -7.5975E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1692.39105449011        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  1.0502E+00  1.9603E+00  4.3939E-01  4.3975E-01  1.1037E+00  9.0056E-01  7.4538E-01  1.3460E+00  1.9182E+00  1.0160E+00
             9.9131E-01
 PARAMETER:  1.4895E-01  7.7312E-01 -7.2236E-01 -7.2155E-01  1.9864E-01 -4.7374E-03 -1.9387E-01  3.9713E-01  7.5138E-01  1.1589E-01
             9.1275E-02
 GRADIENT:   1.6118E+00 -1.1519E+01 -2.1556E-01  3.6569E-01 -2.4009E+00  1.5011E-01  1.1369E+00 -1.2899E-01 -7.0249E-01  2.1215E-01
             1.3879E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1692.41392032154        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1241
 NPARAMETR:  1.0504E+00  1.9651E+00  4.4159E-01  4.3573E-01  1.1090E+00  9.0031E-01  7.3739E-01  1.3927E+00  1.9387E+00  1.0176E+00
             9.9154E-01
 PARAMETER:  1.4916E-01  7.7556E-01 -7.1738E-01 -7.3074E-01  2.0342E-01 -5.0206E-03 -2.0464E-01  4.3122E-01  7.6202E-01  1.1748E-01
             9.1507E-02
 GRADIENT:   2.2816E+00 -1.4302E+01 -3.2583E-01  2.5234E-01 -1.1711E+00  7.0399E-02 -1.9658E-01 -8.1130E-02 -1.0970E-01 -2.6418E-01
            -5.4910E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1692.42055235852        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1426             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0510E+00  1.9643E+00  4.4512E-01  4.3571E-01  1.1118E+00  9.0027E-01  7.3788E-01  1.4289E+00  1.9454E+00  1.0217E+00
             9.9166E-01
 PARAMETER:  1.4970E-01  7.7513E-01 -7.0941E-01 -7.3078E-01  2.0598E-01 -5.0590E-03 -2.0398E-01  4.5688E-01  7.6548E-01  1.2144E-01
             9.1629E-02
 GRADIENT:   6.6962E+02  9.9273E+02  3.0136E+00  9.7047E+01  1.4909E+01  3.1411E+01  1.5138E+01  5.9289E-01  3.4783E+01  1.1444E+00
             9.0245E-01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1692.42170472665        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1518
 NPARAMETR:  1.0503E+00  1.9652E+00  4.4566E-01  4.3545E-01  1.1118E+00  9.0014E-01  7.3759E-01  1.4223E+00  1.9439E+00  1.0213E+00
             9.9162E-01
 PARAMETER:  1.4912E-01  7.7562E-01 -7.0820E-01 -7.3137E-01  2.0594E-01 -5.2024E-03 -2.0436E-01  4.5231E-01  7.6470E-01  1.2112E-01
             9.1583E-02
 GRADIENT:  -1.5826E+00  3.0988E-01 -2.2979E-01  2.4172E-01 -7.2064E-01 -5.5375E-02  9.4387E-02  1.6932E-02 -1.2094E-01 -3.2500E-02
             2.4960E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1518
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0935E-03 -4.2788E-02 -2.4370E-02  3.8421E-02 -4.9624E-02
 SE:             2.9838E-02  2.3365E-02  9.3090E-03  2.2704E-02  2.1760E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7077E-01  6.7058E-02  8.8471E-03  9.0591E-02  2.2578E-02

 ETASHRINKSD(%)  4.0373E-02  2.1724E+01  6.8814E+01  2.3940E+01  2.7101E+01
 ETASHRINKVR(%)  8.0729E-02  3.8729E+01  9.0274E+01  4.2149E+01  4.6857E+01
 EBVSHRINKSD(%)  5.4696E-01  2.1043E+01  7.1606E+01  2.5056E+01  2.4280E+01
 EBVSHRINKVR(%)  1.0909E+00  3.7657E+01  9.1938E+01  4.3834E+01  4.2665E+01
 RELATIVEINF(%)  9.8832E+01  6.3681E+00  1.1054E+00  5.6978E+00  1.8118E+01
 EPSSHRINKSD(%)  4.6246E+01
 EPSSHRINKVR(%)  7.1105E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1692.4217047266461     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -957.27087816290793     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1692.422       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.97E+00  4.46E-01  4.35E-01  1.11E+00  9.00E-01  7.38E-01  1.42E+00  1.94E+00  1.02E+00  9.92E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       44.660
Stop Time:
Sun Oct 24 03:27:35 CDT 2021
