Sat Sep 25 11:28:49 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat4.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m4.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1524.00071483321        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.5131E+01 -4.9777E+01 -2.2601E+01 -3.2611E+00  1.4517E+02 -3.2622E+01 -1.2477E+01 -1.7969E+01 -2.9763E+00 -6.5657E+01
            -1.3012E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1552.09662243658        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8979E-01  1.0209E+00  9.3094E-01  9.9749E-01  8.7458E-01  1.0753E+00  1.0422E+00  1.0799E+00  9.8329E-01  1.1703E+00
             1.2682E+00
 PARAMETER:  8.9737E-02  1.2065E-01  2.8442E-02  9.7490E-02 -3.4015E-02  1.7258E-01  1.4134E-01  1.7689E-01  8.3147E-02  2.5722E-01
             3.3758E-01
 GRADIENT:   3.8220E+01  8.7195E-01  1.0026E+01 -1.0232E+01 -2.0601E+01  1.5921E+00 -1.7659E+00 -1.5945E+00  1.0412E+00  4.9103E+00
             3.7707E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1552.22289531507        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.8528E-01  1.0831E+00  9.0749E-01  9.6013E-01  8.8779E-01  1.0897E+00  1.1524E+00  1.2258E+00  9.1738E-01  1.0996E+00
             1.2667E+00
 PARAMETER:  8.5169E-02  1.7980E-01  2.9311E-03  5.9313E-02 -1.9018E-02  1.8587E-01  2.4188E-01  3.0356E-01  1.3766E-02  1.9493E-01
             3.3639E-01
 GRADIENT:   2.8729E+01  8.5354E+00  9.8384E+00 -1.0076E+01 -1.9095E+01  7.3456E+00  2.6598E+00  8.4455E-01 -4.2122E+00 -1.2861E+00
             2.3609E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1552.70255825195        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  9.8396E-01  1.0966E+00  8.6615E-01  9.5586E-01  8.9062E-01  1.0888E+00  1.0824E+00  1.1113E+00  9.7705E-01  1.1212E+00
             1.2577E+00
 PARAMETER:  8.3828E-02  1.9219E-01 -4.3700E-02  5.4859E-02 -1.5842E-02  1.8504E-01  1.7921E-01  2.0551E-01  7.6783E-02  2.1444E-01
             3.2932E-01
 GRADIENT:  -1.1481E+00  1.0743E+00  8.4203E-01  5.9651E-01 -2.4479E+00 -5.5794E-01  7.0663E-04  1.0933E-01  3.8271E-01 -6.5075E-02
             1.2711E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1552.75329509569        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      474
 NPARAMETR:  9.8568E-01  1.2481E+00  7.2571E-01  8.5458E-01  8.9364E-01  1.0923E+00  9.9483E-01  1.0054E+00  1.0417E+00  1.1052E+00
             1.2572E+00
 PARAMETER:  8.5579E-02  3.2159E-01 -2.2061E-01 -5.7149E-02 -1.2450E-02  1.8832E-01  9.4812E-02  1.0537E-01  1.4084E-01  2.0007E-01
             3.2890E-01
 GRADIENT:  -1.9338E-01  2.3746E+00  4.2552E-01  1.8030E+00 -1.6627E+00  4.1333E-02  1.5715E-01  1.0138E-01 -1.3002E-01  1.7970E-01
            -1.6962E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1552.82118469875        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      651
 NPARAMETR:  9.8659E-01  1.4408E+00  5.4416E-01  7.2132E-01  8.9553E-01  1.0933E+00  9.0215E-01  7.7147E-01  1.1499E+00  1.0809E+00
             1.2594E+00
 PARAMETER:  8.6503E-02  4.6520E-01 -5.0852E-01 -2.2667E-01 -1.0336E-02  1.8924E-01 -2.9729E-03 -1.5945E-01  2.3966E-01  1.7779E-01
             3.3062E-01
 GRADIENT:  -1.1703E-01  3.1325E+00  6.3224E-01  1.2649E+00 -2.2036E+00 -1.3533E-01  1.0782E-01  5.9776E-02  8.5393E-02  5.8637E-02
             1.5983E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1552.84472652278        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      827
 NPARAMETR:  9.8670E-01  1.5368E+00  4.6990E-01  6.5367E-01  9.0952E-01  1.0938E+00  8.6589E-01  6.3411E-01  1.2134E+00  1.0817E+00
             1.2602E+00
 PARAMETER:  8.6607E-02  5.2967E-01 -6.5523E-01 -3.2516E-01  5.1628E-03  1.8962E-01 -4.3997E-02 -3.5554E-01  2.9343E-01  1.7855E-01
             3.3128E-01
 GRADIENT:   8.1003E-02 -5.8012E-01 -3.5590E-01 -1.3842E-01  2.6846E-01 -4.1521E-02  2.2481E-01  1.4144E-01 -1.9183E-02  1.1234E-01
             1.6090E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1552.85674266869        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1004
 NPARAMETR:  9.8677E-01  1.5864E+00  4.2527E-01  6.1730E-01  9.1325E-01  1.0941E+00  8.4747E-01  4.5103E-01  1.2515E+00  1.0803E+00
             1.2603E+00
 PARAMETER:  8.6685E-02  5.6147E-01 -7.5502E-01 -3.8240E-01  9.2516E-03  1.8995E-01 -6.5501E-02 -6.9622E-01  3.2434E-01  1.7722E-01
             3.3133E-01
 GRADIENT:   4.0019E-01 -3.1307E+00 -7.2140E-01 -8.3375E-01  1.5085E+00  9.3824E-02  4.7024E-02  8.3883E-02  1.6599E-02  1.8814E-01
             4.4836E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1552.87854475152        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1179
 NPARAMETR:  9.8653E-01  1.5889E+00  4.0837E-01  6.1489E-01  9.0136E-01  1.0939E+00  8.4762E-01  1.7600E-01  1.2518E+00  1.0697E+00
             1.2601E+00
 PARAMETER:  8.6443E-02  5.6305E-01 -7.9557E-01 -3.8631E-01 -3.8549E-03  1.8972E-01 -6.5320E-02 -1.6373E+00  3.2460E-01  1.6739E-01
             3.3117E-01
 GRADIENT:  -2.5150E-02  8.2658E-01  3.1730E-01  1.0847E-01 -4.7033E-01 -7.2534E-03  2.5990E-03  2.9807E-03 -3.8243E-02 -6.6921E-02
            -7.9018E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1552.88312152903        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1354
 NPARAMETR:  9.8651E-01  1.6027E+00  3.9859E-01  6.0517E-01  9.0389E-01  1.0939E+00  8.4264E-01  2.7440E-02  1.2626E+00  1.0701E+00
             1.2604E+00
 PARAMETER:  8.6421E-02  5.7170E-01 -8.1981E-01 -4.0224E-01 -1.0528E-03  1.8972E-01 -7.1211E-02 -3.4958E+00  3.3319E-01  1.6777E-01
             3.3141E-01
 GRADIENT:   2.4647E-02  1.5856E-01 -2.1915E-02  1.0787E-01 -3.1936E-02  5.9215E-03 -9.8120E-03  2.0199E-04 -9.2879E-03 -5.7263E-03
            -3.2303E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1552.88322173490        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1529
 NPARAMETR:  9.8649E-01  1.6042E+00  3.9835E-01  6.0417E-01  9.0478E-01  1.0938E+00  8.4224E-01  1.0000E-02  1.2642E+00  1.0710E+00
             1.2604E+00
 PARAMETER:  8.6403E-02  5.7262E-01 -8.2043E-01 -4.0390E-01 -5.8148E-05  1.8970E-01 -7.1686E-02 -4.7524E+00  3.3440E-01  1.6856E-01
             3.3144E-01
 GRADIENT:  -1.9614E-03 -1.9245E-02 -5.2605E-03 -2.7854E-03  6.4166E-03 -1.4363E-04  2.7435E-03  0.0000E+00  2.4964E-03  2.9694E-03
             1.1353E-03

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1552.88322173490        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1557
 NPARAMETR:  9.8652E-01  1.6042E+00  3.9839E-01  6.0416E-01  9.0476E-01  1.0939E+00  8.4216E-01  1.0000E-02  1.2641E+00  1.0709E+00
             1.2604E+00
 PARAMETER:  8.6403E-02  5.7262E-01 -8.2043E-01 -4.0390E-01 -5.8148E-05  1.8970E-01 -7.1686E-02 -4.7524E+00  3.3440E-01  1.6856E-01
             3.3144E-01
 GRADIENT:  -1.2187E-02  7.3543E-03 -5.1109E-03  1.0003E-03  6.8174E-03 -4.8539E-03  3.2355E-03  0.0000E+00  2.2831E-03  2.5507E-03
             1.0898E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1557
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8608E-04 -2.7389E-02 -2.7030E-04  2.1341E-02 -3.1377E-02
 SE:             2.9793E-02  2.3368E-02  1.0375E-04  2.2927E-02  2.2628E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9502E-01  2.4116E-01  9.1808E-03  3.5194E-01  1.6556E-01

 ETASHRINKSD(%)  1.9120E-01  2.1715E+01  9.9652E+01  2.3191E+01  2.4192E+01
 ETASHRINKVR(%)  3.8203E-01  3.8714E+01  9.9999E+01  4.1004E+01  4.2532E+01
 EBVSHRINKSD(%)  5.6591E-01  2.1328E+01  9.9710E+01  2.4891E+01  2.2492E+01
 EBVSHRINKVR(%)  1.1286E+00  3.8107E+01  9.9999E+01  4.3587E+01  3.9925E+01
 RELATIVEINF(%)  9.8768E+01  3.6316E+00  7.4454E-05  2.9706E+00  9.1032E+00
 EPSSHRINKSD(%)  4.3982E+01
 EPSSHRINKVR(%)  6.8620E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1552.8832217348956     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -817.73239517115746     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1552.883       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.60E+00  3.98E-01  6.04E-01  9.05E-01  1.09E+00  8.42E-01  1.00E-02  1.26E+00  1.07E+00  1.26E+00
 


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
+        9.43E+02
 
 TH 2
+       -8.41E+00  3.92E+02
 
 TH 3
+        4.13E-01  1.76E+02  6.64E+02
 
 TH 4
+       -1.53E+01  3.36E+02 -4.73E+02  1.08E+03
 
 TH 5
+       -4.13E+00 -2.41E+02 -5.34E+02  3.47E+02  7.08E+02
 
 TH 6
+       -2.64E-01 -1.74E+00  1.22E+00 -4.56E+00 -6.35E-01  1.62E+02
 
 TH 7
+       -1.71E+00  4.57E+00 -3.24E+01 -2.37E+00 -1.16E+01  1.75E+00  1.03E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.51E+00 -1.52E+01 -3.91E+01  4.92E+01  7.65E-01 -7.88E-01  2.36E+01  0.00E+00  4.70E+01
 
 TH10
+       -1.03E+00 -1.51E+01 -4.51E+01 -1.07E+00 -5.39E+01  4.21E-01  1.12E+01  0.00E+00  6.39E+00  6.90E+01
 
 TH11
+       -7.22E+00 -1.47E+01 -2.25E+01  1.71E-01 -2.42E+00  2.54E+00  1.11E+01  0.00E+00  5.93E+00  1.30E+01  1.32E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.459
Stop Time:
Sat Sep 25 11:29:15 CDT 2021
