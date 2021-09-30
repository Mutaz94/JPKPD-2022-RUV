Wed Sep 29 15:21:45 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat80.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1713.44509180243        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0551E+02  2.2245E+01  1.7154E+01  4.9931E+01 -2.8027E+01  8.2508E+01  6.9979E+00 -9.1330E-01  4.5246E+01  4.6002E+00
            -1.2099E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1718.31480163529        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0126E+00  1.0906E+00  9.8929E-01  9.7279E-01  1.0659E+00  8.1758E-01  9.8436E-01  1.0081E+00  7.5254E-01  9.8836E-01
             1.0779E+00
 PARAMETER:  1.1248E-01  1.8677E-01  8.9235E-02  7.2416E-02  1.6381E-01 -1.0140E-01  8.4236E-02  1.0806E-01 -1.8430E-01  8.8289E-02
             1.7506E-01
 GRADIENT:  -3.1457E+01  1.3773E+01  1.2302E+00  2.0106E+01  6.8058E+00 -4.1402E+01 -7.8404E+00 -3.2582E+00 -1.1969E+01 -2.9698E+00
             1.2107E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1719.26074960757        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0160E+00  1.1505E+00  8.7759E-01  9.2674E-01  1.0265E+00  8.3781E-01  9.8896E-01  1.1338E+00  7.4108E-01  9.0046E-01
             1.0553E+00
 PARAMETER:  1.1591E-01  2.4021E-01 -3.0575E-02  2.3916E-02  1.2620E-01 -7.6959E-02  8.8899E-02  2.2562E-01 -1.9965E-01 -4.8493E-03
             1.5379E-01
 GRADIENT:  -2.1957E+01  9.3445E+00 -4.2966E+00  1.4739E+01  1.8669E+00 -3.0687E+01 -5.2647E+00  2.4992E+00 -1.2859E+01 -2.8818E+00
             5.2910E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1721.53876029443        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0235E+00  1.2131E+00  8.0825E-01  8.7847E-01  1.0285E+00  9.0331E-01  9.5675E-01  9.1674E-01  8.4269E-01  9.4791E-01
             1.0358E+00
 PARAMETER:  1.2326E-01  2.9321E-01 -1.1289E-01 -2.9570E-02  1.2809E-01 -1.6846E-03  5.5783E-02  1.3072E-02 -7.1154E-02  4.6502E-02
             1.3518E-01
 GRADIENT:   4.5955E-01  8.9542E-01  3.2595E-01  1.1188E+00 -9.6014E-01  5.2572E-01  1.0293E-01  3.6179E-02 -1.6051E-01  1.3436E-01
            -5.5253E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1721.56847283909        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0244E+00  1.3382E+00  6.7053E-01  7.9390E-01  1.0140E+00  8.9877E-01  8.9932E-01  7.7869E-01  8.8700E-01  9.1649E-01
             1.0363E+00
 PARAMETER:  1.2409E-01  3.9135E-01 -2.9968E-01 -1.3080E-01  1.1395E-01 -6.7234E-03 -6.1178E-03 -1.5014E-01 -1.9907E-02  1.2797E-02
             1.3564E-01
 GRADIENT:  -6.3083E-01  4.8302E+00  1.7202E+00  1.5793E+00 -4.8909E+00 -2.0568E+00  2.7481E-01  2.4428E-01 -3.7548E-01  7.0057E-03
             1.7697E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1721.58495683867        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  1.0252E+00  1.4083E+00  6.0609E-01  7.4879E-01  1.0155E+00  9.0288E-01  8.6730E-01  6.8121E-01  9.1908E-01  9.1459E-01
             1.0354E+00
 PARAMETER:  1.2486E-01  4.4240E-01 -4.0073E-01 -1.8930E-01  1.1534E-01 -2.1664E-03 -4.2365E-02 -2.8388E-01  1.5619E-02  1.0721E-02
             1.3476E-01
 GRADIENT:   2.2287E-01  8.9780E+00  1.2947E+00  5.2146E+00 -4.6681E+00 -5.2765E-01  6.1724E-02  2.1020E-01 -2.2316E-01  3.8482E-01
            -1.6401E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1721.61860013236        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1082
 NPARAMETR:  1.0258E+00  1.5108E+00  5.1451E-01  6.7805E-01  1.0239E+00  9.0906E-01  8.2271E-01  4.8187E-01  9.7567E-01  9.1140E-01
             1.0348E+00
 PARAMETER:  1.2543E-01  5.1262E-01 -5.6454E-01 -2.8853E-01  1.2366E-01  4.6562E-03 -9.5146E-02 -6.3008E-01  7.5371E-02  7.2295E-03
             1.3425E-01
 GRADIENT:   6.8148E-01  5.6482E+00 -5.1604E-01  5.4655E+00 -5.0284E-02  1.8614E+00 -3.9843E-01  1.5296E-01  1.0972E-02  3.7994E-01
            -3.5632E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1721.65245171411        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1263
 NPARAMETR:  1.0255E+00  1.5224E+00  5.0183E-01  6.6697E-01  1.0254E+00  9.0470E-01  8.2086E-01  3.9753E-01  9.8563E-01  9.1048E-01
             1.0363E+00
 PARAMETER:  1.2513E-01  5.2028E-01 -5.8949E-01 -3.0501E-01  1.2511E-01 -1.5470E-04 -9.7403E-02 -8.2248E-01  8.5526E-02  6.2123E-03
             1.3570E-01
 GRADIENT:  -1.6062E-01  5.5088E-01 -1.8491E-02  1.5790E+00  1.1201E+00 -1.7108E-02  1.7158E-01  6.0881E-02  6.6225E-02  8.4152E-02
             2.1704E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1721.70317409774        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1443
 NPARAMETR:  1.0252E+00  1.5397E+00  4.6332E-01  6.5351E-01  1.0065E+00  9.0476E-01  8.1440E-01  6.4328E-02  9.9715E-01  8.9268E-01
             1.0346E+00
 PARAMETER:  1.2489E-01  5.3160E-01 -6.6934E-01 -3.2540E-01  1.0653E-01 -8.7435E-05 -1.0531E-01 -2.6438E+00  9.7143E-02 -1.3528E-02
             1.3398E-01
 GRADIENT:  -1.1179E+00  2.6629E+00 -9.1690E-01  4.2174E+00  6.2190E-01 -1.0474E-01  2.5012E-01  4.4275E-03  7.3641E-01  5.0374E-01
            -6.1170E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1721.71558306494        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1602
 NPARAMETR:  1.0264E+00  1.5424E+00  4.5476E-01  6.4741E-01  1.0026E+00  9.0518E-01  8.1182E-01  3.0318E-02  9.9504E-01  8.8227E-01
             1.0351E+00
 PARAMETER:  1.2602E-01  5.3336E-01 -6.8798E-01 -3.3478E-01  1.0255E-01  3.8081E-04 -1.0848E-01 -3.3960E+00  9.5023E-02 -2.5256E-02
             1.3447E-01
 GRADIENT:   2.0806E+00 -3.2686E+00 -2.2193E-01 -1.1619E+00 -5.7582E-02  9.3012E-02  7.5239E-03  1.2003E-03  1.8568E-02  2.7136E-02
             5.6909E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1721.71664590151        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:     1723
 NPARAMETR:  1.0267E+00  1.5422E+00  4.5515E-01  6.4790E-01  1.0027E+00  9.0527E-01  8.1180E-01  1.0000E-02  9.9492E-01  8.8208E-01
             1.0348E+00
 PARAMETER:  1.2638E-01  5.3324E-01 -6.8714E-01 -3.3402E-01  1.0269E-01  4.7346E-04 -1.0850E-01 -5.2628E+00  9.4912E-02 -2.5475E-02
             1.3425E-01
 GRADIENT:   5.3254E+02  4.2337E+02  7.6176E+00  1.0563E+02  7.9192E+00  4.0836E+01  6.3200E+00  0.0000E+00  3.2051E+00  4.9964E-01
             1.0050E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1721.71700023536        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1900
 NPARAMETR:  1.0266E+00  1.5415E+00  4.5566E-01  6.4879E-01  1.0031E+00  9.0522E-01  8.1213E-01  1.0000E-02  9.9449E-01  8.8212E-01
             1.0348E+00
 PARAMETER:  1.2617E-01  5.3297E-01 -6.8606E-01 -3.3386E-01  1.0255E-01  4.6368E-04 -1.0802E-01 -5.3156E+00  9.4820E-02 -2.4423E-02
             1.3450E-01
 GRADIENT:  -1.1033E-01  2.6449E-01 -1.1312E-02 -7.6287E-01 -4.3886E-01  8.1626E-03  7.0255E-03  0.0000E+00  2.5334E-02  3.3572E-02
             6.3947E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1900
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.0359E-05 -1.8928E-02 -3.2146E-04  1.6369E-02 -2.7043E-02
 SE:             2.9819E-02  2.4546E-02  1.2864E-04  2.2377E-02  2.1800E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9865E-01  4.4062E-01  1.2455E-02  4.6447E-01  2.1479E-01

 ETASHRINKSD(%)  1.0147E-01  1.7769E+01  9.9569E+01  2.5035E+01  2.6967E+01
 ETASHRINKVR(%)  2.0284E-01  3.2381E+01  9.9998E+01  4.3802E+01  4.6662E+01
 EBVSHRINKSD(%)  5.3146E-01  1.7799E+01  9.9618E+01  2.6372E+01  2.5776E+01
 EBVSHRINKVR(%)  1.0601E+00  3.2430E+01  9.9999E+01  4.5789E+01  4.4908E+01
 RELATIVEINF(%)  9.8903E+01  3.2404E+00  9.0382E-05  2.1591E+00  6.8818E+00
 EPSSHRINKSD(%)  4.3892E+01
 EPSSHRINKVR(%)  6.8519E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1721.7170002353578     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -986.56617367161959     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.68
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
 





 #OBJV:********************************************    -1721.717       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.54E+00  4.56E-01  6.48E-01  1.00E+00  9.05E-01  8.12E-01  1.00E-02  9.95E-01  8.83E-01  1.04E+00
 


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
+        1.27E+03
 
 TH 2
+       -9.98E+00  4.93E+02
 
 TH 3
+        9.79E+00  2.30E+02  8.65E+02
 
 TH 4
+       -2.52E+01  4.16E+02 -6.63E+02  1.46E+03
 
 TH 5
+       -6.41E+00 -2.71E+02 -6.75E+02  4.83E+02  7.85E+02
 
 TH 6
+       -4.66E-01 -1.59E+00  2.89E+00 -4.96E+00 -1.07E+00  2.38E+02
 
 TH 7
+        5.26E-01  1.60E+01 -4.22E+01 -1.13E+01  1.82E+00  6.77E-02  1.46E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.75E+00 -1.91E+01 -5.04E+01  6.55E+01 -1.47E+00  5.53E-02  2.04E+01  0.00E+00  7.33E+01
 
 TH10
+       -2.83E-01 -1.58E+01 -4.84E+01 -8.36E+00 -6.94E+01 -4.76E-02  1.84E+01  0.00E+00  1.33E+01  8.43E+01
 
 TH11
+       -8.14E+00 -1.73E+01 -3.25E+01  5.28E+00 -3.28E+00  2.28E+00  1.04E+01  0.00E+00  1.07E+01  1.83E+01  1.98E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.711
Stop Time:
Wed Sep 29 15:22:28 CDT 2021
