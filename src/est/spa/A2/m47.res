Wed Sep 29 12:50:35 CDT 2021
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
$DATA ../../../../data/spa/A2/dat47.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m47.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1099.06883928369        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1112E+02 -7.4426E+00  5.2395E+01 -5.8915E+01  7.4286E+01  5.8130E+01 -1.5913E+01 -4.1229E+01 -5.3057E+01 -5.5481E+01
            -9.2531E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1335.65287363224        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.2426E+00  9.8771E-01  8.9252E-01  1.1189E+00  8.6979E-01  1.2442E+00  9.8104E-01  1.0884E+00  1.0542E+00  1.1344E+00
             3.0891E+00
 PARAMETER:  3.1718E-01  8.7632E-02 -1.3703E-02  2.1233E-01 -3.9506E-02  3.1850E-01  8.0856E-02  1.8469E-01  1.5283E-01  2.2611E-01
             1.2279E+00
 GRADIENT:   4.1956E+02  3.0369E+00  6.7911E+00 -1.2455E-01 -3.1877E+01  4.2956E+01  7.9396E+00  6.9828E+00  1.6501E+01  1.5924E+01
             1.2666E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1372.98269617617        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0611E+00  5.7031E-01  3.3108E-01  1.3254E+00  3.6234E-01  9.4380E-01  3.5393E-01  8.5613E-01  1.1227E+00  4.3322E-01
             2.3531E+00
 PARAMETER:  1.5932E-01 -4.6157E-01 -1.0054E+00  3.8172E-01 -9.1518E-01  4.2157E-02 -9.3866E-01 -5.5335E-02  2.1571E-01 -7.3652E-01
             9.5574E-01
 GRADIENT:   2.2352E+02  1.4188E+02  8.1051E+01  2.2857E+02 -1.2975E+02  8.3254E+00 -5.3661E+00 -6.7007E+00  1.4249E+01 -2.3987E+01
            -8.9186E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1404.69370602068        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      317
 NPARAMETR:  1.0092E+00  6.8104E-01  5.5976E-01  1.2133E+00  5.7960E-01  9.1882E-01  8.4772E-01  6.8765E-01  9.5217E-01  8.0247E-01
             2.2808E+00
 PARAMETER:  1.0913E-01 -2.8414E-01 -4.8025E-01  2.9334E-01 -4.4541E-01  1.5340E-02 -6.5210E-02 -2.7448E-01  5.0992E-02 -1.2007E-01
             9.2454E-01
 GRADIENT:   3.0307E+01 -1.2402E+01 -2.5530E+01 -1.2071E+01  3.4944E+01  6.4149E+00 -2.8014E+00  4.3083E+00 -5.0331E+00 -2.2359E+00
            -2.1723E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1407.70531181946        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      492
 NPARAMETR:  9.9129E-01  5.7016E-01  7.1401E-01  1.3125E+00  6.3433E-01  8.9630E-01  1.0381E+00  4.1331E-01  9.1874E-01  9.0833E-01
             2.3103E+00
 PARAMETER:  9.1255E-02 -4.6184E-01 -2.3686E-01  3.7191E-01 -3.5519E-01 -9.4750E-03  1.3739E-01 -7.8356E-01  1.5247E-02  3.8577E-03
             9.3736E-01
 GRADIENT:  -7.6372E+00  6.0898E+00  3.1956E+00 -8.0648E-04 -2.1646E+00  3.4152E-01 -5.2785E-01 -3.6369E-01 -3.7529E+00 -3.5289E+00
            -5.2502E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1409.48977223954        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      674
 NPARAMETR:  9.8550E-01  2.8297E-01  8.9870E-01  1.5033E+00  6.7018E-01  8.8260E-01  3.6193E-01  7.0858E-02  8.5212E-01  1.0143E+00
             2.3951E+00
 PARAMETER:  8.5389E-02 -1.1624E+00 -6.8109E-03  5.0763E-01 -3.0021E-01 -2.4878E-02 -9.1632E-01 -2.5471E+00 -6.0024E-02  1.1420E-01
             9.7344E-01
 GRADIENT:  -4.1509E+00  4.2196E+00  8.5863E+00  1.5093E+01 -9.3725E+00 -1.4809E+00  3.3153E-02 -9.2448E-03 -1.1625E+00 -7.6166E-01
            -5.1420E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1410.85274564104        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      850
 NPARAMETR:  9.8477E-01  1.2644E-01  6.2240E-01  1.5202E+00  4.8284E-01  8.8787E-01  5.6251E-02  1.2085E-02  8.2056E-01  8.7781E-01
             2.3549E+00
 PARAMETER:  8.4652E-02 -1.9680E+00 -3.7418E-01  5.1887E-01 -6.2807E-01 -1.8933E-02 -2.7779E+00 -4.3158E+00 -9.7774E-02 -3.0328E-02
             9.5649E-01
 GRADIENT:  -1.1356E+00  1.3096E+00  4.6266E+00  1.0260E+01 -8.9582E+00 -5.1345E-01  1.3471E-04 -6.2018E-04 -1.3464E+00  7.5386E-01
            -7.7829E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1411.04516341450        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  9.8042E-01  5.5290E-02  6.3981E-01  1.5569E+00  4.8374E-01  8.8588E-01  1.0000E-02  1.0000E-02  7.9984E-01  8.7103E-01
             2.3808E+00
 PARAMETER:  8.0230E-02 -2.7952E+00 -3.4658E-01  5.4270E-01 -6.2622E-01 -2.1168E-02 -4.6565E+00 -6.5355E+00 -1.2334E-01 -3.8079E-02
             9.6743E-01
 GRADIENT:  -5.1697E+00  2.7406E-01  1.5284E+00  3.9561E+00 -2.8724E+00 -2.0076E-01  0.0000E+00  0.0000E+00 -7.5064E-01 -1.3803E-01
             2.1987E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1411.10316366608        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1203
 NPARAMETR:  9.8157E-01  2.0652E-02  6.3586E-01  1.5703E+00  4.7713E-01  8.8608E-01  1.0000E-02  1.0000E-02  7.9502E-01  8.7097E-01
             2.3686E+00
 PARAMETER:  8.1401E-02 -3.7799E+00 -3.5277E-01  5.5129E-01 -6.3997E-01 -2.0943E-02 -7.0183E+00 -9.1874E+00 -1.2939E-01 -3.8151E-02
             9.6228E-01
 GRADIENT:   2.3506E+00  2.9234E-02 -1.6476E-01  4.5727E-01 -2.9359E-02  1.3138E-01  0.0000E+00  0.0000E+00  3.4014E-01  5.8167E-02
            -7.9768E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1411.10795822686        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1380
 NPARAMETR:  9.8005E-01  1.0270E-02  6.4204E-01  1.5778E+00  4.7889E-01  8.8519E-01  1.0000E-02  1.0000E-02  7.8999E-01  8.7113E-01
             2.3765E+00
 PARAMETER:  7.9851E-02 -4.4785E+00 -3.4310E-01  5.5601E-01 -6.3629E-01 -2.1955E-02 -8.7037E+00 -1.1093E+01 -1.3574E-01 -3.7967E-02
             9.6561E-01
 GRADIENT:  -8.7695E-01  2.4569E-02  4.3097E-01  1.3338E+00 -8.7253E-01 -2.8831E-02  0.0000E+00  0.0000E+00 -1.5146E-01 -9.3495E-03
             3.0925E-01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1411.11045582790        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1480
 NPARAMETR:  9.8044E-01  1.0000E-02  6.4176E-01  1.5759E+00  4.7912E-01  8.8528E-01  1.0000E-02  1.0000E-02  7.9045E-01  8.7119E-01
             2.3749E+00
 PARAMETER:  8.0242E-02 -4.5984E+00 -3.4354E-01  5.5484E-01 -6.3580E-01 -2.1853E-02 -8.7037E+00 -1.1093E+01 -1.3515E-01 -3.7896E-02
             9.6497E-01
 GRADIENT:   4.8091E-01  0.0000E+00 -3.2729E-01 -2.4054E+00  8.7081E-01  3.1838E-02  0.0000E+00  0.0000E+00  3.5937E-02 -6.3115E-02
             7.3312E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1480
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.5352E-04 -2.2109E-06 -5.8700E-06 -1.1273E-02 -1.5834E-02
 SE:             2.9116E-02  1.4468E-06  2.0460E-04  2.7015E-02  2.2200E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9579E-01  1.2648E-01  9.7711E-01  6.7648E-01  4.7569E-01

 ETASHRINKSD(%)  2.4589E+00  9.9995E+01  9.9315E+01  9.4960E+00  2.5627E+01
 ETASHRINKVR(%)  4.8574E+00  1.0000E+02  9.9995E+01  1.8090E+01  4.4687E+01
 EBVSHRINKSD(%)  2.4967E+00  9.9995E+01  9.9334E+01  9.0786E+00  2.4774E+01
 EBVSHRINKVR(%)  4.9311E+00  1.0000E+02  9.9996E+01  1.7333E+01  4.3410E+01
 RELATIVEINF(%)  8.1340E+01  1.1638E-08  1.7507E-04  7.5870E+00  1.5304E+00
 EPSSHRINKSD(%)  3.3629E+01
 EPSSHRINKVR(%)  5.5949E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1411.1104558278969     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -675.95962926415871     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.29
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1411.110       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.00E-02  6.42E-01  1.58E+00  4.79E-01  8.85E-01  1.00E-02  1.00E-02  7.90E-01  8.71E-01  2.37E+00
 


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
+        1.40E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -4.38E+00  0.00E+00  1.19E+03
 
 TH 4
+       -4.73E+01  0.00E+00 -1.29E+02  5.93E+02
 
 TH 5
+        5.70E+01  0.00E+00 -1.96E+03 -1.35E+02  3.58E+03
 
 TH 6
+       -2.08E+00  0.00E+00  7.27E+00 -1.01E+01  1.37E+00  2.29E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.86E+00  0.00E+00  2.14E+01 -1.18E+01  1.94E+00  1.49E+00  0.00E+00  0.00E+00  2.14E+02
 
 TH10
+       -5.67E+00  0.00E+00 -1.69E+01 -1.54E+00 -5.19E+01 -1.66E+00  0.00E+00  0.00E+00  5.60E-01  8.91E+01
 
 TH11
+       -1.53E+01  0.00E+00 -1.05E+01 -9.47E+00 -5.39E-01  4.34E+00  0.00E+00  0.00E+00  1.36E+01  2.00E+01  4.85E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.278
Stop Time:
Wed Sep 29 12:51:01 CDT 2021
