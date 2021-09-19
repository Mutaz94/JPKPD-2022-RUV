Sat Sep 18 12:48:26 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1632.33651718191        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.7715E+01 -4.7575E+01 -6.1965E+01  1.6566E+01  1.3082E+02  1.4024E+01 -5.6526E+00  9.0268E+00 -5.9966E+00 -3.9731E+01
            -9.1643E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1651.11625364518        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0426E+00  1.0124E+00  1.0556E+00  9.9365E-01  9.6848E-01  9.4947E-01  1.0045E+00  9.5431E-01  1.0204E+00  1.1009E+00
             1.2172E+00
 PARAMETER:  1.4171E-01  1.1236E-01  1.5415E-01  9.3633E-02  6.7970E-02  4.8148E-02  1.0451E-01  5.3236E-02  1.2023E-01  1.9613E-01
             2.9653E-01
 GRADIENT:   6.8527E+00 -2.0091E+01 -2.2226E+01 -6.3695E+00  2.9668E+01 -2.7374E-01  1.3098E+00  8.0559E+00  5.7201E-01 -6.2162E+00
             1.0012E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1653.58297293165        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      180
 NPARAMETR:  1.0433E+00  9.1518E-01  1.2473E+00  1.0768E+00  1.0109E+00  9.4958E-01  8.0536E-01  6.7971E-01  1.0376E+00  1.2951E+00
             1.1849E+00
 PARAMETER:  1.4242E-01  1.1363E-02  3.2102E-01  1.7399E-01  1.1087E-01  4.8260E-02 -1.1647E-01 -2.8609E-01  1.3691E-01  3.5856E-01
             2.6968E-01
 GRADIENT:  -3.2162E+01 -1.2726E+00 -8.3123E+00  1.0685E+01  1.5960E+01 -2.8303E+00 -1.1814E+00  2.4548E-01  1.2549E+00  1.0148E+00
            -7.4360E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1654.15272401480        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      357
 NPARAMETR:  1.0571E+00  8.7896E-01  1.1801E+00  1.0913E+00  9.5396E-01  9.5747E-01  9.9990E-01  6.0846E-01  9.8320E-01  1.2169E+00
             1.1892E+00
 PARAMETER:  1.5549E-01 -2.9019E-02  2.6556E-01  1.8736E-01  5.2865E-02  5.6538E-02  9.9897E-02 -3.9683E-01  8.3052E-02  2.9627E-01
             2.7326E-01
 GRADIENT:  -5.6511E-01 -2.3728E+00 -1.0102E+00 -2.0205E+00  1.6974E+00  5.3759E-01  1.3000E-01  3.6629E-01  2.9813E-01 -1.4830E-01
             3.0499E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1654.32486930042        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      536
 NPARAMETR:  1.0589E+00  1.0527E+00  1.0326E+00  9.8138E-01  9.6129E-01  9.5771E-01  9.5146E-01  3.8336E-01  1.0521E+00  1.1953E+00
             1.1889E+00
 PARAMETER:  1.5720E-01  1.5136E-01  1.3206E-01  8.1208E-02  6.0520E-02  5.6795E-02  5.0238E-02 -8.5877E-01  1.5074E-01  2.7844E-01
             2.7305E-01
 GRADIENT:  -7.5141E-01  2.9936E+00  4.1041E-01  4.0791E+00 -1.5961E+00 -1.8407E-01  5.3758E-01  4.6183E-02 -6.1899E-02  4.9825E-01
            -2.5596E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1654.44477500714        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      713
 NPARAMETR:  1.0613E+00  1.2230E+00  9.3157E-01  8.6733E-01  9.9396E-01  9.5966E-01  8.6161E-01  2.3311E-01  1.1520E+00  1.1957E+00
             1.1890E+00
 PARAMETER:  1.5950E-01  3.0130E-01  2.9115E-02 -4.2331E-02  9.3943E-02  5.8820E-02 -4.8951E-02 -1.3563E+00  2.4146E-01  2.7874E-01
             2.7315E-01
 GRADIENT:   2.3830E+00  8.0372E-01  4.7654E-01  3.4373E-01 -1.8778E+00  2.0041E-01 -1.2912E-01  2.6183E-02 -2.1161E-01  2.4438E-02
            -2.1671E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1654.48144006377        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      888
 NPARAMETR:  1.0606E+00  1.3449E+00  8.7691E-01  7.8781E-01  1.0315E+00  9.5946E-01  8.1353E-01  1.5368E-01  1.2392E+00  1.2132E+00
             1.1883E+00
 PARAMETER:  1.5886E-01  3.9633E-01 -3.1354E-02 -1.3849E-01  1.3097E-01  5.8618E-02 -1.0637E-01 -1.7729E+00  3.1447E-01  2.9328E-01
             2.7250E-01
 GRADIENT:  -7.9746E-02  5.8470E-01  4.4705E-02  4.4545E-01 -7.4150E-02 -8.9882E-03 -5.1819E-02  1.3254E-02 -6.3938E-02 -1.1490E-01
            -3.9926E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1654.48246971400        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1066
 NPARAMETR:  1.0607E+00  1.3679E+00  8.6362E-01  7.7208E-01  1.0379E+00  9.5959E-01  8.0632E-01  1.2931E-01  1.2564E+00  1.2153E+00
             1.1880E+00
 PARAMETER:  1.5897E-01  4.1330E-01 -4.6623E-02 -1.5867E-01  1.3717E-01  5.8753E-02 -1.1528E-01 -1.9456E+00  3.2828E-01  2.9499E-01
             2.7229E-01
 GRADIENT:   1.1598E-02 -3.0105E-01 -2.1334E-01 -1.5147E-01  4.2701E-01  1.6684E-02 -3.2511E-02  1.1536E-02 -1.0646E-01 -1.1685E-01
            -3.3531E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1654.48759482614        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1243
 NPARAMETR:  1.0608E+00  1.3593E+00  8.6340E-01  7.7768E-01  1.0330E+00  9.5961E-01  8.1064E-01  2.0333E-02  1.2493E+00  1.2125E+00
             1.1884E+00
 PARAMETER:  1.5899E-01  4.0698E-01 -4.6880E-02 -1.5144E-01  1.3245E-01  5.8771E-02 -1.0993E-01 -3.7955E+00  3.2255E-01  2.9271E-01
             2.7258E-01
 GRADIENT:   5.0578E-02 -4.3204E-02 -4.2587E-02 -2.8798E-02 -3.5775E-02  1.4054E-02 -1.6339E-02  2.1206E-04 -8.4734E-03  2.5127E-02
            -9.7080E-04

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1654.48770643681        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1370
 NPARAMETR:  1.0607E+00  1.3576E+00  8.6455E-01  7.7882E-01  1.0327E+00  9.5956E-01  8.1141E-01  1.0000E-02  1.2480E+00  1.2124E+00
             1.1884E+00
 PARAMETER:  1.5896E-01  4.0575E-01 -4.5551E-02 -1.4998E-01  1.3217E-01  5.8723E-02 -1.0898E-01 -4.6163E+00  3.2155E-01  2.9264E-01
             2.7264E-01
 GRADIENT:  -5.3174E-03 -7.9669E-04 -7.1225E-03  5.1652E-03  5.3550E-03 -1.8267E-03  2.9354E-03  0.0000E+00  1.0448E-03  2.4568E-03
             6.3877E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1370
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.2469E-04 -2.8354E-02 -2.6014E-04  1.2577E-02 -3.0566E-02
 SE:             2.9792E-02  1.8418E-02  1.1191E-04  2.4352E-02  2.4267E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8863E-01  1.2370E-01  2.0094E-02  6.0553E-01  2.0782E-01

 ETASHRINKSD(%)  1.9222E-01  3.8297E+01  9.9625E+01  1.8418E+01  1.8704E+01
 ETASHRINKVR(%)  3.8408E-01  6.1927E+01  9.9999E+01  3.3443E+01  3.3909E+01
 EBVSHRINKSD(%)  6.2920E-01  3.7385E+01  9.9671E+01  1.9317E+01  1.5998E+01
 EBVSHRINKVR(%)  1.2544E+00  6.0793E+01  9.9999E+01  3.4903E+01  2.9437E+01
 RELATIVEINF(%)  9.8457E+01  1.6540E+00  1.7528E-04  3.2166E+00  1.3144E+01
 EPSSHRINKSD(%)  4.2439E+01
 EPSSHRINKVR(%)  6.6868E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1654.4877064368122     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -919.33687987307405     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.86
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1654.488       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.36E+00  8.65E-01  7.79E-01  1.03E+00  9.60E-01  8.11E-01  1.00E-02  1.25E+00  1.21E+00  1.19E+00
 


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
+        1.06E+03
 
 TH 2
+       -9.00E+00  3.83E+02
 
 TH 3
+        1.12E+01  1.03E+02  1.75E+02
 
 TH 4
+       -2.09E+01  3.94E+02 -9.64E+01  7.66E+02
 
 TH 5
+       -5.08E+00 -1.86E+02 -2.18E+02  1.02E+02  4.99E+02
 
 TH 6
+        6.42E-01 -5.25E-01 -1.12E+00 -1.62E+00 -1.40E+00  2.11E+02
 
 TH 7
+        1.01E+00 -3.22E-01  1.27E+01 -1.72E+01 -1.07E+01 -1.26E+00  5.35E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.59E-01 -2.34E+01 -1.43E+01  3.83E+01  1.37E+00  3.27E-01  2.30E+01  0.00E+00  6.23E+01
 
 TH10
+       -3.53E-01 -9.47E+00 -2.63E+01 -4.61E+00 -4.44E+01 -1.03E+00  1.18E+01  0.00E+00  2.64E+00  6.59E+01
 
 TH11
+       -7.33E+00 -1.81E+01 -2.24E+01 -4.76E+00  8.33E+00  3.75E+00  6.47E+00  0.00E+00  6.19E+00  1.41E+01  1.58E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.836
Stop Time:
Sat Sep 18 12:48:50 CDT 2021
