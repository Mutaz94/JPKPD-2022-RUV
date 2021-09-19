Sat Sep 18 09:11:36 CDT 2021
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
$DATA ../../../../data/spa/A1/dat33.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m33.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1045.09761204200        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6752E+02  9.5043E+00  2.0289E+01  1.1171E+01  2.4909E+01  1.8484E+01 -1.5892E+01 -7.0308E+00  2.8049E+00 -2.3631E+01
            -1.1250E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1416.10135101864        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.4177E-01  1.1288E+00  1.1836E+00  9.3124E-01  1.1530E+00  8.3804E-01  1.0139E+00  9.3041E-01  8.2078E-01  7.9156E-01
             2.2917E+00
 PARAMETER:  4.0006E-02  2.2113E-01  2.6859E-01  2.8764E-02  2.4236E-01 -7.6693E-02  1.1385E-01  2.7866E-02 -9.7502E-02 -1.3374E-01
             9.2928E-01
 GRADIENT:  -2.7505E+01 -1.1177E+01 -9.2723E+00  1.5926E+00  2.2815E+01 -3.9415E+01 -7.1632E+00  1.0397E-01 -5.4297E+00  1.8565E-01
            -4.8899E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1420.98657013880        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.4754E-01  1.3480E+00  1.4759E+00  8.2224E-01  1.2747E+00  9.3431E-01  8.9081E-01  1.2337E+00  9.7428E-01  4.5416E-01
             2.5503E+00
 PARAMETER:  4.6116E-02  3.9865E-01  4.8924E-01 -9.5721E-02  3.4267E-01  3.2053E-02 -1.5624E-02  3.0999E-01  7.3942E-02 -6.8931E-01
             1.0362E+00
 GRADIENT:  -1.3257E+01  3.4507E+01  3.7606E+00  2.4690E+01 -1.7673E+01  3.7703E+00  1.4992E-01 -1.2644E-01  3.1754E+00 -4.2387E-03
            -4.1036E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1421.64208346450        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.5177E-01  1.1395E+00  1.4785E+00  9.3900E-01  1.2167E+00  9.2508E-01  1.0425E+00  1.3067E+00  8.2255E-01  4.2201E-01
             2.5509E+00
 PARAMETER:  5.0572E-02  2.3061E-01  4.9103E-01  3.7061E-02  2.9618E-01  2.2125E-02  1.4159E-01  3.6747E-01 -9.5341E-02 -7.6272E-01
             1.0364E+00
 GRADIENT:   7.0243E-01  3.3178E+00  5.1214E-02  3.8000E+00 -1.3549E+00  1.0463E+00  1.6046E-01  1.1679E-01  3.4876E-01  3.9365E-01
             9.2682E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1421.70371750238        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.5132E-01  1.1155E+00  1.3736E+00  9.4538E-01  1.1814E+00  9.2092E-01  1.0838E+00  1.2862E+00  7.9520E-01  2.2336E-01
             2.5592E+00
 PARAMETER:  5.0091E-02  2.0927E-01  4.1746E-01  4.3831E-02  2.6672E-01  1.7621E-02  1.8043E-01  3.5167E-01 -1.2917E-01 -1.3990E+00
             1.0397E+00
 GRADIENT:   2.8460E-01 -3.2457E+00 -3.8702E-01 -3.0542E+00  1.0572E+00 -3.5643E-01 -6.4477E-02  8.1271E-02  2.5091E-03  6.0751E-02
             5.8148E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1421.80116918888        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  9.5213E-01  1.2307E+00  1.1068E+00  8.6782E-01  1.1398E+00  9.2234E-01  1.0395E+00  1.2018E+00  8.1417E-01  2.2353E-02
             2.5508E+00
 PARAMETER:  5.0948E-02  3.0761E-01  2.0151E-01 -4.1772E-02  2.3087E-01  1.9163E-02  1.3876E-01  2.8385E-01 -1.0559E-01 -3.7008E+00
             1.0364E+00
 GRADIENT:  -8.4616E-01  3.0744E+00  6.1745E-01  1.7198E+00 -2.5093E+00 -1.6708E-01  1.0243E-01  1.4690E-01 -4.7435E-02  7.8899E-04
             3.5168E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1421.89203664290        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      553
 NPARAMETR:  9.5514E-01  1.3350E+00  1.0266E+00  8.0151E-01  1.1622E+00  9.2498E-01  9.7744E-01  1.1387E+00  8.5212E-01  1.0000E-02
             2.5575E+00
 PARAMETER:  5.4099E-02  3.8893E-01  1.2624E-01 -1.2126E-01  2.5030E-01  2.2016E-02  7.7186E-02  2.2991E-01 -6.0025E-02 -4.5226E+00
             1.0390E+00
 GRADIENT:  -7.9612E-01  3.4566E+00  1.1382E+00  2.4299E+00 -3.1031E+00  9.4151E-02 -2.9965E-01 -1.6382E-01  4.5017E-02  0.0000E+00
            -1.1200E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1421.99989612771        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      731
 NPARAMETR:  9.5799E-01  1.5262E+00  8.1834E-01  6.7257E-01  1.1833E+00  9.2613E-01  8.9375E-01  1.1410E+00  9.1190E-01  1.0000E-02
             2.5630E+00
 PARAMETER:  5.7079E-02  5.2277E-01 -1.0048E-01 -2.9666E-01  2.6827E-01  2.3261E-02 -1.2329E-02  2.3195E-01  7.7770E-03 -5.6944E+00
             1.0412E+00
 GRADIENT:   3.5693E+00  1.1800E+00  2.8954E-01  6.0331E-01 -1.1304E+00  2.0482E-01 -2.5770E-02  5.9527E-02 -1.7145E-01  0.0000E+00
             5.0621E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1422.03618279127        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      906
 NPARAMETR:  9.5747E-01  1.6525E+00  6.6881E-01  5.8583E-01  1.1947E+00  9.2713E-01  8.4297E-01  9.7236E-01  9.8734E-01  1.0000E-02
             2.5678E+00
 PARAMETER:  5.6539E-02  6.0232E-01 -3.0225E-01 -4.3473E-01  2.7787E-01  2.4334E-02 -7.0824E-02  7.1966E-02  8.7255E-02 -5.8187E+00
             1.0430E+00
 GRADIENT:   7.3514E-01 -2.0045E+00 -2.0116E-01 -5.2112E-01  7.4517E-01  3.5799E-01  5.1225E-02 -1.7766E-02  1.2913E-01  0.0000E+00
             7.7163E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1422.04411613194        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1082
 NPARAMETR:  9.5741E-01  1.7285E+00  6.0642E-01  5.3639E-01  1.2128E+00  9.2651E-01  8.1598E-01  9.7989E-01  1.0264E+00  1.0000E-02
             2.5690E+00
 PARAMETER:  5.6474E-02  6.4724E-01 -4.0019E-01 -5.2290E-01  2.9294E-01  2.3667E-02 -1.0336E-01  7.9688E-02  1.2610E-01 -5.7273E+00
             1.0435E+00
 GRADIENT:  -9.4228E-02 -7.2447E-02 -4.7769E-02  7.2220E-02  2.0340E-01 -4.5962E-03 -9.4363E-02 -1.5672E-02 -5.3337E-02  0.0000E+00
             3.4535E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1422.04441914414        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1257
 NPARAMETR:  9.5743E-01  1.7346E+00  6.0604E-01  5.3251E-01  1.2165E+00  9.2648E-01  8.1408E-01  1.0113E+00  1.0310E+00  1.0000E-02
             2.5689E+00
 PARAMETER:  5.6498E-02  6.5075E-01 -4.0081E-01 -5.3015E-01  2.9595E-01  2.3638E-02 -1.0569E-01  1.1123E-01  1.3048E-01 -5.7121E+00
             1.0435E+00
 GRADIENT:  -6.3628E-04  4.8223E-03  7.0355E-03 -8.8597E-03 -1.2744E-02 -2.7263E-04  1.8092E-03 -1.2220E-03 -3.7729E-04  0.0000E+00
            -2.4135E-04

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1422.04441995196        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1422
 NPARAMETR:  9.5743E-01  1.7346E+00  6.0602E-01  5.3252E-01  1.2165E+00  9.2648E-01  8.1407E-01  1.0125E+00  1.0309E+00  1.0000E-02
             2.5689E+00
 PARAMETER:  5.6499E-02  6.5075E-01 -4.0084E-01 -5.3013E-01  2.9595E-01  2.3638E-02 -1.0571E-01  1.1243E-01  1.3043E-01 -5.7121E+00
             1.0435E+00
 GRADIENT:  -6.0559E-04  5.9163E-03 -9.5386E-05  3.2396E-03  2.4807E-03 -4.7015E-04 -1.9944E-03 -1.0457E-05 -4.3184E-04  0.0000E+00
            -9.0040E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1422
 NO. OF SIG. DIGITS IN FINAL EST.:  4.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1366E-04 -4.3317E-03 -6.2759E-03 -4.3355E-03 -9.5832E-05
 SE:             2.9111E-02  2.4144E-02  5.8294E-03  1.3738E-02  1.8924E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7770E-01  8.5761E-01  2.8165E-01  7.5232E-01  6.1256E-01

 ETASHRINKSD(%)  2.4743E+00  1.9114E+01  8.0471E+01  5.3975E+01  9.9366E+01
 ETASHRINKVR(%)  4.8875E+00  3.4574E+01  9.6186E+01  7.8817E+01  9.9996E+01
 EBVSHRINKSD(%)  2.5915E+00  1.9201E+01  8.1259E+01  5.4521E+01  9.9360E+01
 EBVSHRINKVR(%)  5.1159E+00  3.4715E+01  9.6488E+01  7.9317E+01  9.9996E+01
 RELATIVEINF(%)  9.3287E+01  1.2412E+00  9.0188E-02  4.2054E-01  1.8648E-04
 EPSSHRINKSD(%)  2.5395E+01
 EPSSHRINKVR(%)  4.4341E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1422.0444199519609     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -686.89359338822271     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1422.044       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.57E-01  1.73E+00  6.06E-01  5.33E-01  1.22E+00  9.26E-01  8.14E-01  1.01E+00  1.03E+00  1.00E-02  2.57E+00
 


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
+        1.34E+03
 
 TH 2
+       -4.35E+01  3.66E+02
 
 TH 3
+        2.43E+01  7.63E+01  1.13E+02
 
 TH 4
+       -9.96E+01  3.85E+02 -9.57E+01  7.57E+02
 
 TH 5
+       -1.71E+01 -1.53E+02 -1.72E+02  9.65E+01  2.94E+02
 
 TH 6
+       -1.00E+00 -9.53E+00  4.64E+00 -2.20E+01 -6.54E+00  2.07E+02
 
 TH 7
+        5.66E+00  7.96E+00 -1.98E+00 -3.38E+01  8.47E+00  1.67E+00  1.28E+02
 
 TH 8
+       -2.14E-01 -2.35E+00 -3.56E+00  2.36E+00  2.63E+00  8.62E-01  1.32E+00  1.01E-01
 
 TH 9
+        8.24E-01 -5.78E+00 -2.48E+00  8.12E-01  1.20E+00  1.44E+00  1.54E+01  8.83E-01  1.32E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.46E+01 -8.88E+00 -9.30E+00 -9.76E-01 -3.87E+00  3.34E+00  1.33E+01  1.98E+00  4.86E+00  0.00E+00  5.82E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.682
Stop Time:
Sat Sep 18 09:12:00 CDT 2021
