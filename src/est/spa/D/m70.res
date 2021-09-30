Wed Sep 29 20:16:43 CDT 2021
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
$DATA ../../../../data/spa/D/dat70.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m70.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13555.4256585459        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.0814E+02  3.1583E+02 -7.5946E+01  2.0782E+01  3.8831E+02 -2.4438E+03 -7.5365E+02 -6.1377E+01 -1.6645E+03 -8.8152E+02
            -2.4328E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -622.819790025951        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4954E+00  9.7302E-01  7.9980E-01  1.8373E+00  1.5795E+00  2.3035E+00  1.2136E+00  9.6315E-01  1.2098E+00  1.1591E+00
             1.4305E+01
 PARAMETER:  5.0242E-01  7.2648E-02 -1.2339E-01  7.0830E-01  5.5713E-01  9.3444E-01  2.9357E-01  6.2454E-02  2.9047E-01  2.4768E-01
             2.7606E+00
 GRADIENT:   3.0627E+01  2.1691E+01 -1.1260E+01  4.2584E+01 -4.7695E+00  6.1244E+01  2.7573E-01  7.3948E+00  3.8244E+00  3.2253E-01
             8.2203E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -649.482242525262        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.4029E+00  7.4163E-01  8.3841E-01  2.0342E+00  4.5738E+00  2.2262E+00  1.8607E+00  3.3332E-01  2.6935E+00  1.0602E+01
             1.1003E+01
 PARAMETER:  4.3851E-01 -1.9890E-01 -7.6250E-02  8.1012E-01  1.6203E+00  9.0029E-01  7.2098E-01 -9.9865E-01  1.0908E+00  2.4611E+00
             2.4982E+00
 GRADIENT:   3.3486E+01  1.0905E+01 -7.3924E+00  5.9140E+01 -1.2732E+01 -2.6923E+01  8.6410E+00  4.2008E-02  5.1157E+01  1.9278E+01
             5.3526E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -674.571692567776        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.1461E+00  3.6654E-01  5.9518E-01  1.6314E+00  6.3429E+00  2.0992E+00  1.3947E+00  4.3705E-02  1.3911E+00  1.0379E+01
             1.1063E+01
 PARAMETER:  2.3634E-01 -9.0366E-01 -4.1889E-01  5.8943E-01  1.9473E+00  8.4155E-01  4.3269E-01 -3.0303E+00  4.3011E-01  2.4398E+00
             2.5036E+00
 GRADIENT:  -2.6897E+01  9.0064E+00  1.0441E+01  1.2192E+01 -1.1853E+00  1.5018E+01  1.8332E+00 -1.6493E-02 -7.0452E+00  3.8396E+00
             7.4729E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -697.465081527459        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.1276E+00  9.7088E-02  3.9317E-01  1.4020E+00  1.1862E+01  2.0458E+00  1.1736E+00  4.4389E-02  1.5836E+00  5.0245E+00
             8.0670E+00
 PARAMETER:  2.2010E-01 -2.2321E+00 -8.3351E-01  4.3791E-01  2.5733E+00  8.1578E-01  2.6005E-01 -3.0148E+00  5.5969E-01  1.7143E+00
             2.1878E+00
 GRADIENT:   5.1450E+01 -3.1074E-01  6.1978E+01  1.5755E+00  3.2623E-01 -4.3091E+01  9.2161E-02 -4.3201E-02  4.9528E+00 -3.2457E-02
            -1.2931E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -759.998672312406        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      381
 NPARAMETR:  5.2240E-01  1.0000E-02  3.2649E-02  3.8633E-01  8.4027E+02  1.2783E+00  1.0000E-02  1.0000E-02  1.7486E-01  2.8915E+00
             1.1024E+01
 PARAMETER: -5.4932E-01 -5.6962E+00 -3.3219E+00 -8.5107E-01  6.8337E+00  3.4551E-01 -7.0400E+00 -1.7117E+01 -1.6438E+00  1.1618E+00
             2.5001E+00
 GRADIENT:  -1.4324E+01  0.0000E+00 -9.9628E+01  2.5712E+02  1.1994E-03 -1.5865E+02  0.0000E+00  0.0000E+00  3.1090E-01 -3.9200E-06
             1.5914E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -798.167406447179        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      559
 NPARAMETR:  5.0447E-01  1.0000E-02  2.4244E-02  2.8192E-01  1.1235E+03  1.8815E+00  1.0000E-02  1.0000E-02  1.8999E-01  2.4498E+00
             1.0918E+01
 PARAMETER: -5.8425E-01 -6.0746E+00 -3.6196E+00 -1.1661E+00  7.1242E+00  7.3209E-01 -1.0010E+01 -2.1212E+01 -1.5608E+00  9.9602E-01
             2.4904E+00
 GRADIENT:   3.0743E+01  0.0000E+00 -1.3731E+01  4.7545E+00  2.1702E-03  1.9506E+01  0.0000E+00  0.0000E+00 -8.3863E-01 -1.2276E-07
             3.0841E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -804.103749678316        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      716
 NPARAMETR:  4.4900E-01  1.0000E-02  2.1056E-02  2.5193E-01  1.7086E+03  1.7682E+00  1.0000E-02  1.0000E-02  6.5744E-01  2.4724E+00
             1.0320E+01
 PARAMETER: -7.0072E-01 -6.3013E+00 -3.7605E+00 -1.2786E+00  7.5434E+00  6.6997E-01 -1.0061E+01 -2.1412E+01 -3.1940E-01  1.0052E+00
             2.4341E+00
 GRADIENT:   3.5422E+01  0.0000E+00  4.8858E+01  5.5179E+01  2.4911E-05  2.7351E+00  0.0000E+00  0.0000E+00  1.1502E+00  2.0312E-08
             7.4765E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -804.615887789656        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      788
 NPARAMETR:  4.3697E-01  1.0000E-02  1.9766E-02  2.4788E-01  1.7023E+03  1.7878E+00  1.0000E-02  1.0000E-02  7.2715E-01  3.1304E+00
             9.7026E+00
 PARAMETER: -7.2789E-01 -6.3013E+00 -3.8238E+00 -1.2948E+00  7.5398E+00  6.8098E-01 -1.0061E+01 -2.1412E+01 -2.1863E-01  1.2412E+00
             2.3724E+00
 GRADIENT:   4.5789E+01  0.0000E+00 -1.5251E+00  1.2591E+02  7.0376E-04  4.3155E+00  0.0000E+00  0.0000E+00 -2.7790E+00  6.3898E-08
             4.0662E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -808.220062890602        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      969
 NPARAMETR:  4.5078E-01  1.0000E-02  2.0813E-02  2.4546E-01  1.6988E+03  1.8624E+00  1.0000E-02  1.0000E-02  7.8605E-01  3.2184E+00
             9.1847E+00
 PARAMETER: -6.9679E-01 -6.3013E+00 -3.7722E+00 -1.3046E+00  7.5377E+00  7.2185E-01 -1.0061E+01 -2.1412E+01 -1.4073E-01  1.2689E+00
             2.3175E+00
 GRADIENT:   5.1336E+00  0.0000E+00 -5.6115E+00  2.3505E+00  4.8738E-04  1.1506E+00  0.0000E+00  0.0000E+00  1.3368E+00  6.2753E-08
             1.9071E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -808.234693977350        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  4.4962E-01  1.0000E-02  2.0989E-02  2.4643E-01  1.6983E+03  1.8570E+00  1.0000E-02  1.0000E-02  7.7333E-01  3.1625E+00
             9.2062E+00
 PARAMETER: -6.9935E-01 -6.3013E+00 -3.7638E+00 -1.3007E+00  7.5374E+00  7.1896E-01 -1.0061E+01 -2.1412E+01 -1.5705E-01  1.2514E+00
             2.3199E+00
 GRADIENT:   1.5175E+00  0.0000E+00 -5.9315E-01 -1.8518E+00  3.5723E-04  2.5942E-01  0.0000E+00  0.0000E+00  9.0718E-02  5.5365E-08
             1.3087E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1026
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0595E-04  2.4073E-06  1.3657E-04 -1.9903E-02 -1.0542E-06
 SE:             2.9154E-02  1.1035E-06  2.5682E-04  2.2531E-02  1.6906E-06
 N:                     100         100         100         100         100

 P VAL.:         9.9163E-01  2.9144E-02  5.9488E-01  3.7704E-01  5.3292E-01

 ETASHRINKSD(%)  2.3314E+00  9.9996E+01  9.9140E+01  2.4520E+01  9.9994E+01
 ETASHRINKVR(%)  4.6085E+00  1.0000E+02  9.9993E+01  4.3027E+01  1.0000E+02
 EBVSHRINKSD(%)  2.1236E+00  9.9995E+01  9.9205E+01  2.5638E+01  9.9994E+01
 EBVSHRINKVR(%)  4.2022E+00  1.0000E+02  9.9994E+01  4.4703E+01  1.0000E+02
 RELATIVEINF(%)  2.4259E+00  6.3294E-08  4.3711E-05  3.2869E-01  1.5225E-08
 EPSSHRINKSD(%)  1.3706E+01
 EPSSHRINKVR(%)  2.5533E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -808.23469397735005     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -73.083867413611870     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.01
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -808.235       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.50E-01  1.00E-02  2.10E-02  2.46E-01  1.70E+03  1.86E+00  1.00E-02  1.00E-02  7.73E-01  3.16E+00  9.21E+00
 


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
+        1.52E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.06E+04  0.00E+00  1.48E+06
 
 TH 4
+       -2.28E+02  0.00E+00 -1.46E+05  1.57E+04
 
 TH 5
+        1.13E-05  0.00E+00 -2.63E-04  1.92E-05 -6.02E-11
 
 TH 6
+        2.82E+00  0.00E+00  3.42E+01 -3.38E+01  2.76E-07  5.16E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.03E+01  0.00E+00  1.31E+03 -1.35E+02 -4.63E-07 -2.48E+00  0.00E+00  0.00E+00  8.10E+01
 
 TH10
+       -1.41E-04  0.00E+00 -2.49E-04  4.59E-05 -1.40E-09  1.09E-05  0.00E+00  0.00E+00 -3.76E-05  4.18E-05
 
 TH11
+       -1.49E+01  0.00E+00  3.20E+02 -2.43E+01 -1.79E-07  1.52E+00  0.00E+00  0.00E+00  5.61E+00  4.08E-06  4.13E+00
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       22.308
Stop Time:
Wed Sep 29 20:17:07 CDT 2021
