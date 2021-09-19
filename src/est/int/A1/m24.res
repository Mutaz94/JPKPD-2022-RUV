Fri Sep 17 23:51:04 CDT 2021
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
$DATA ../../../../data/int/A1/dat24.csv ignore=@
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
Current Date:       17 SEP 2021
Days until program expires : 212
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m24.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3106.00681693799        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1428E+02  6.8334E+01  5.4961E+01 -1.2829E+02  3.3411E+01  4.6970E+00 -5.6583E+01 -1.9454E+01 -8.1647E+00 -2.7824E+01
            -1.3574E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3371.98668034564        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7736E-01  7.8224E-01  8.5865E-01  1.1822E+00  8.0568E-01  9.4997E-01  1.0005E+00  8.5136E-01  8.9809E-01  8.3096E-01
             1.4102E+00
 PARAMETER:  7.7099E-02 -1.4559E-01 -5.2395E-02  2.6739E-01 -1.1606E-01  4.8676E-02  1.0054E-01 -6.0922E-02 -7.4893E-03 -8.5173E-02
             4.4372E-01
 GRADIENT:   3.7426E+01 -1.2696E+01  2.2877E+01  3.5549E+01  1.4515E+01 -1.3759E+01 -1.3316E+01  1.0866E+01 -2.2517E+01 -8.6105E+00
            -1.8095E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3376.39866833929        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.8049E-01  7.1245E-01  6.8200E-01  1.2262E+00  6.8698E-01  9.5079E-01  1.0816E+00  3.9049E-01  9.1204E-01  8.4910E-01
             1.4060E+00
 PARAMETER:  8.0299E-02 -2.3905E-01 -2.8272E-01  3.0393E-01 -2.7545E-01  4.9541E-02  1.7845E-01 -8.4035E-01  7.9316E-03 -6.3576E-02
             4.4077E-01
 GRADIENT:   4.5828E+01  4.6483E+01 -5.5002E+01  7.4182E+01  3.7185E+01 -1.3688E+01  3.0935E+00  1.5515E+00 -2.7580E+01 -7.3900E+00
            -1.8083E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3386.59367084636        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.6315E-01  6.6350E-01  6.6197E-01  1.2148E+00  6.4700E-01  9.8351E-01  9.8685E-01  3.0021E-01  9.7179E-01  8.8766E-01
             1.5158E+00
 PARAMETER:  6.2454E-02 -3.1022E-01 -3.1254E-01  2.9457E-01 -3.3541E-01  8.3372E-02  8.6761E-02 -1.1033E+00  7.1385E-02 -1.9164E-02
             5.1591E-01
 GRADIENT:  -7.7076E-01  9.1803E+00 -9.3356E+00  1.4046E+01  4.5025E+00  7.4036E-01 -1.8336E+00  2.2551E+00 -1.2033E+00 -6.0115E-01
            -9.4169E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3386.60787319043        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.6307E-01  6.4754E-01  6.4819E-01  1.2161E+00  6.3181E-01  9.8280E-01  9.9191E-01  2.4857E-01  9.7203E-01  8.8726E-01
             1.5169E+00
 PARAMETER:  6.2368E-02 -3.3458E-01 -3.3357E-01  2.9564E-01 -3.5917E-01  8.2653E-02  9.1873E-02 -1.2920E+00  7.1636E-02 -1.9622E-02
             5.1669E-01
 GRADIENT:  -9.4345E-01  5.1097E+00 -5.6906E+00  7.8515E+00  2.2838E+00  5.2578E-01 -1.1946E+00  1.6193E+00 -9.4260E-01 -3.2504E-01
             1.2869E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3386.61010524546        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.6313E-01  6.4005E-01  6.4138E-01  1.2167E+00  6.2459E-01  9.8232E-01  9.9425E-01  2.1454E-01  9.7254E-01  8.8781E-01
             1.5169E+00
 PARAMETER:  6.2435E-02 -3.4621E-01 -3.4413E-01  2.9618E-01 -3.7066E-01  8.2165E-02  9.4233E-02 -1.4392E+00  7.2155E-02 -1.8999E-02
             5.1666E-01
 GRADIENT:  -7.6010E-01  3.4289E+00 -4.0716E+00  5.2359E+00  1.4694E+00  3.6821E-01 -8.5633E-01  1.2204E+00 -7.2461E-01 -1.7913E-01
             1.2959E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3386.61109960980        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  9.6319E-01  6.3549E-01  6.3716E-01  1.2172E+00  6.2018E-01  9.8202E-01  9.9562E-01  1.8887E-01  9.7293E-01  8.8833E-01
             1.5168E+00
 PARAMETER:  6.2492E-02 -3.5335E-01 -3.5073E-01  2.9652E-01 -3.7775E-01  8.1861E-02  9.5610E-02 -1.5667E+00  7.2556E-02 -1.8411E-02
             5.1660E-01
 GRADIENT:  -6.0793E-01  2.4651E+00 -3.0669E+00  3.7428E+00  1.0234E+00  2.7008E-01 -6.4446E-01  9.4899E-01 -5.6934E-01 -1.0442E-01
             1.1469E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3386.61179072762        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  9.6323E-01  6.3264E-01  6.3448E-01  1.2174E+00  6.1741E-01  9.8184E-01  9.9644E-01  1.6992E-01  9.7321E-01  8.8874E-01
             1.5167E+00
 PARAMETER:  6.2533E-02 -3.5785E-01 -3.5495E-01  2.9674E-01 -3.8223E-01  8.1674E-02  9.6437E-02 -1.6724E+00  7.2841E-02 -1.7950E-02
             5.1655E-01
 GRADIENT:  -4.9978E-01  1.8907E+00 -2.4366E+00  2.8584E+00  7.6693E-01  2.0969E-01 -5.1035E-01  7.6870E-01 -4.6296E-01 -6.3037E-02
             9.9713E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3386.61244675012        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  9.6326E-01  6.3037E-01  6.3233E-01  1.2177E+00  6.1520E-01  9.8170E-01  9.9708E-01  1.5262E-01  9.7344E-01  8.8911E-01
             1.5166E+00
 PARAMETER:  6.2568E-02 -3.6144E-01 -3.5834E-01  2.9692E-01 -3.8581E-01  8.1528E-02  9.7072E-02 -1.7798E+00  7.3082E-02 -1.7538E-02
             5.1650E-01
 GRADIENT:  -4.0819E-01  1.4508E+00 -1.9342E+00  2.1840E+00  5.7603E-01  1.6266E-01 -4.0344E-01  6.2004E-01 -3.7446E-01 -3.5011E-02
             8.5171E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3386.61288189947        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      658
 NPARAMETR:  9.6328E-01  6.2882E-01  6.3085E-01  1.2178E+00  6.1368E-01  9.8160E-01  9.9750E-01  1.3924E-01  9.7361E-01  8.8938E-01
             1.5166E+00
 PARAMETER:  6.2593E-02 -3.6391E-01 -3.6069E-01  2.9705E-01 -3.8828E-01  8.1429E-02  9.7496E-02 -1.8716E+00  7.3252E-02 -1.7233E-02
             5.1647E-01
 GRADIENT:  -3.4267E-01  1.1583E+00 -1.5892E+00  1.7363E+00  4.5131E-01  1.3108E-01 -3.3027E-01  5.1595E-01 -3.1239E-01 -1.8680E-02
             7.3965E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3386.61336218092        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      731
 NPARAMETR:  9.6330E-01  6.2765E-01  6.2972E-01  1.2179E+00  6.1253E-01  9.8153E-01  9.9781E-01  1.2796E-01  9.7374E-01  8.8960E-01
             1.5166E+00
 PARAMETER:  6.2613E-02 -3.6578E-01 -3.6248E-01  2.9715E-01 -3.9016E-01  8.1355E-02  9.7808E-02 -1.9561E+00  7.3384E-02 -1.6988E-02
             5.1644E-01
 GRADIENT:  -2.9131E-01  9.4285E-01 -1.3274E+00  1.4075E+00  3.6078E-01  1.0746E-01 -2.7493E-01  4.3564E-01 -2.6433E-01 -8.0230E-03
             6.4648E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3386.61365840572        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      802
 NPARAMETR:  9.6332E-01  6.2643E-01  6.2854E-01  1.2181E+00  6.1133E-01  9.8145E-01  9.9813E-01  1.1496E-01  9.7387E-01  8.8983E-01
             1.5165E+00
 PARAMETER:  6.2634E-02 -3.6772E-01 -3.6436E-01  2.9726E-01 -3.9211E-01  8.1280E-02  9.8126E-02 -2.0631E+00  7.3524E-02 -1.6722E-02
             5.1642E-01
 GRADIENT:  -2.3669E-01  7.2354E-01 -1.0577E+00  1.0740E+00  2.7003E-01  8.3396E-02 -2.1768E-01  3.5165E-01 -2.1374E-01  2.5402E-03
             5.4485E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3386.61389132220        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      874
 NPARAMETR:  9.6334E-01  6.2555E-01  6.2769E-01  1.2182E+00  6.1047E-01  9.8140E-01  9.9835E-01  1.0448E-01  9.7397E-01  8.9001E-01
             1.5165E+00
 PARAMETER:  6.2650E-02 -3.6913E-01 -3.6571E-01  2.9734E-01 -3.9352E-01  8.1226E-02  9.8351E-02 -2.1588E+00  7.3626E-02 -1.6524E-02
             5.1640E-01
 GRADIENT:  -1.9692E-01  5.6856E-01 -8.6247E-01  8.3837E-01  2.0626E-01  6.6220E-02 -1.7653E-01  2.9048E-01 -1.7701E-01  8.8779E-03
             4.6796E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3386.61421056892        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      944
 NPARAMETR:  9.6335E-01  6.2466E-01  6.2682E-01  1.2183E+00  6.0960E-01  9.8135E-01  9.9858E-01  9.2516E-02  9.7407E-01  8.9019E-01
             1.5165E+00
 PARAMETER:  6.2666E-02 -3.7055E-01 -3.6709E-01  2.9742E-01 -3.9495E-01  8.1172E-02  9.8574E-02 -2.2804E+00  7.3730E-02 -1.6316E-02
             5.1638E-01
 GRADIENT:  -1.5583E-01  4.1506E-01 -6.6522E-01  6.0583E-01  1.4416E-01  4.9165E-02 -1.3499E-01  2.2790E-01 -1.3930E-01  1.4382E-02
             3.8626E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3386.61453563972        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1014
 NPARAMETR:  9.6337E-01  6.2385E-01  6.2603E-01  1.2183E+00  6.0881E-01  9.8130E-01  9.9877E-01  7.9954E-02  9.7417E-01  8.9037E-01
             1.5164E+00
 PARAMETER:  6.2680E-02 -3.7184E-01 -3.6835E-01  2.9749E-01 -3.9626E-01  8.1123E-02  9.8772E-02 -2.4263E+00  7.3828E-02 -1.6118E-02
             5.1636E-01
 GRADIENT:  -1.1778E-01  2.7860E-01 -4.8702E-01  3.9996E-01  8.9983E-02  3.3811E-02 -9.7310E-02  1.7041E-01 -1.0445E-01  1.8962E-02
             3.0845E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -3386.61479519859        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1085
 NPARAMETR:  9.6338E-01  6.2321E-01  6.2541E-01  1.2184E+00  6.0818E-01  9.8126E-01  9.9893E-01  6.8319E-02  9.7424E-01  8.9051E-01
             1.5164E+00
 PARAMETER:  6.2692E-02 -3.7287E-01 -3.6935E-01  2.9755E-01 -3.9729E-01  8.1085E-02  9.8928E-02 -2.5836E+00  7.3904E-02 -1.5961E-02
             5.1634E-01
 GRADIENT:  -8.7631E-02  1.7193E-01 -3.4478E-01  2.3861E-01  4.7542E-02  2.1736E-02 -6.7519E-02  1.2464E-01 -7.6892E-02  2.1608E-02
             2.4523E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -3386.61534438277        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1155
 NPARAMETR:  9.6339E-01  6.2257E-01  6.2478E-01  1.2185E+00  6.0754E-01  9.8123E-01  9.9908E-01  5.3758E-02  9.7432E-01  8.9066E-01
             1.5164E+00
 PARAMETER:  6.2704E-02 -3.7390E-01 -3.7036E-01  2.9761E-01 -3.9834E-01  8.1047E-02  9.9081E-02 -2.8233E+00  7.3983E-02 -1.5794E-02
             5.1633E-01
 GRADIENT:  -5.6575E-02  6.6834E-02 -1.9974E-01  8.0910E-02  6.6098E-03  9.7591E-03 -3.7423E-02  7.7437E-02 -4.8509E-02  2.3134E-02
             1.7713E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -3386.61667993146        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1225
 NPARAMETR:  9.6340E-01  6.2204E-01  6.2425E-01  1.2186E+00  6.0702E-01  9.8120E-01  9.9920E-01  3.6982E-02  9.7439E-01  8.9080E-01
             1.5164E+00
 PARAMETER:  6.2715E-02 -3.7474E-01 -3.7120E-01  2.9766E-01 -3.9919E-01  8.1017E-02  9.9196E-02 -3.1973E+00  7.4054E-02 -1.5641E-02
             5.1631E-01
 GRADIENT:  -2.8494E-02 -1.2021E-02 -8.0774E-02 -3.5699E-02 -2.2052E-02  4.6717E-04 -1.2733E-02  3.6884E-02 -2.3672E-02  2.2616E-02
             1.1087E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -3386.62054187121        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1298
 NPARAMETR:  9.6341E-01  6.2175E-01  6.2394E-01  1.2186E+00  6.0672E-01  9.8118E-01  9.9924E-01  1.5955E-02  9.7444E-01  8.9091E-01
             1.5163E+00
 PARAMETER:  6.2724E-02 -3.7522E-01 -3.7170E-01  2.9771E-01 -3.9969E-01  8.1004E-02  9.9235E-02 -4.0380E+00  7.4108E-02 -1.5511E-02
             5.1630E-01
 GRADIENT:  -6.8221E-03 -4.1202E-02 -3.0855E-03 -7.2906E-02 -2.5894E-02 -3.5351E-03  2.1375E-03  6.9792E-03 -4.9584E-03  1.4357E-02
             4.5038E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -3387.46300815629        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:     1453
 NPARAMETR:  9.7104E-01  6.6704E-01  6.6748E-01  1.2136E+00  6.5382E-01  9.8544E-01  9.8727E-01  1.0000E-02  9.7646E-01  9.1578E-01
             1.5236E+00
 PARAMETER:  7.0610E-02 -3.0490E-01 -3.0424E-01  2.9363E-01 -3.2493E-01  8.5333E-02  8.7184E-02 -4.7843E+00  7.6182E-02  1.2018E-02
             5.2109E-01
 GRADIENT:   7.4520E-01 -1.3072E+00  1.5487E+00 -5.0441E-02 -1.6363E-01  1.9698E-02  2.7662E-01  0.0000E+00  1.1801E-01  1.9026E-01
             3.2571E+00

0ITERATION NO.:   99    OBJECTIVE VALUE:  -3387.46719786240        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1580
 NPARAMETR:  9.7069E-01  6.6736E-01  6.6658E-01  1.2136E+00  6.5344E-01  9.8544E-01  9.8504E-01  1.0000E-02  9.7650E-01  9.1539E-01
             1.5214E+00
 PARAMETER:  7.0251E-02 -3.0443E-01 -3.0560E-01  2.9356E-01 -3.2551E-01  8.5331E-02  8.4925E-02 -4.7945E+00  7.6220E-02  1.1599E-02
             5.1963E-01
 GRADIENT:  -1.5061E-02  4.1720E-03  8.8388E-03  8.5396E-03 -9.8423E-03  5.1397E-03 -4.5239E-03  0.0000E+00  6.8212E-06 -3.6448E-03
             1.8307E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1580
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4103E-04 -1.4896E-02 -1.7570E-04  2.0810E-03 -1.0892E-02
 SE:             2.9734E-02  2.1063E-02  2.7185E-04  2.8731E-02  2.6232E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7743E-01  4.7942E-01  5.1808E-01  9.4226E-01  6.7797E-01

 ETASHRINKSD(%)  3.8738E-01  2.9437E+01  9.9089E+01  3.7463E+00  1.2121E+01
 ETASHRINKVR(%)  7.7325E-01  5.0209E+01  9.9992E+01  7.3522E+00  2.2772E+01
 EBVSHRINKSD(%)  6.4836E-01  2.8870E+01  9.9084E+01  3.7881E+00  1.2843E+01
 EBVSHRINKVR(%)  1.2925E+00  4.9406E+01  9.9992E+01  7.4328E+00  2.4037E+01
 RELATIVEINF(%)  9.8698E+01  1.4821E+01  2.0094E-03  7.7435E+01  1.3479E+01
 EPSSHRINKSD(%)  1.9177E+01
 EPSSHRINKVR(%)  3.4677E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3387.4671978623978     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1733.3778380939871     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.86
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3387.467       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.71E-01  6.67E-01  6.67E-01  1.21E+00  6.53E-01  9.85E-01  9.85E-01  1.00E-02  9.77E-01  9.15E-01  1.52E+00
 


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
+        1.20E+03
 
 TH 2
+       -1.96E+00  1.31E+03
 
 TH 3
+       -7.10E-01 -2.73E+02  1.97E+03
 
 TH 4
+       -4.70E+00  1.76E+02 -5.61E+01  7.33E+02
 
 TH 5
+       -4.81E-01 -8.73E+02 -1.44E+03  1.11E+02  2.35E+03
 
 TH 6
+       -4.39E-02 -2.68E+00  2.04E+00 -2.07E+00 -2.35E+00  2.00E+02
 
 TH 7
+        3.16E-01  2.04E+01 -1.13E+01  2.16E+00 -6.69E+00  1.82E+00  5.36E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.77E+00 -2.84E+01  4.15E+01  8.66E+00 -1.56E+01 -3.43E-01 -4.51E-01  0.00E+00  1.78E+02
 
 TH10
+       -3.15E+00 -1.80E+01 -2.96E+01  3.95E+00 -1.64E+01  1.95E+00  2.62E+01  0.00E+00 -5.77E-01  1.39E+02
 
 TH11
+       -1.13E+01 -2.18E+01 -4.96E+01 -1.29E+01  1.02E+01  1.87E+00  1.46E+01  0.00E+00  9.29E+00  1.09E+01  4.80E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.391
Stop Time:
Fri Sep 17 23:51:43 CDT 2021
