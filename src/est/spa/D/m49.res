Sat Sep 25 14:23:42 CDT 2021
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
$DATA ../../../../data/spa/D/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12174.2924504366        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7668E+02  1.9501E+02  9.2335E-01  1.9961E+02  1.9266E+02 -1.6962E+03 -6.4040E+02 -1.0083E+02 -8.3419E+02 -6.2712E+02
            -2.3716E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -572.805712445508        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2541E+00  1.1879E+00  9.5721E-01  1.8989E+00  1.5636E+00  2.6019E+00  1.8041E+00  1.0304E+00  1.9086E+00  1.1671E+00
             1.3092E+01
 PARAMETER:  3.2643E-01  2.7218E-01  5.6272E-02  7.4127E-01  5.4698E-01  1.0562E+00  6.9003E-01  1.2992E-01  7.4637E-01  2.5454E-01
             2.6720E+00
 GRADIENT:  -1.8258E+01  2.7004E+01 -2.0793E+01  5.6874E+01 -6.4027E+00  6.5012E+01 -5.1733E+00  5.8070E+00 -1.5865E+00  8.1503E-01
             1.0988E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -599.619917156424        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.2429E+00  1.0861E+00  1.7327E+00  1.9522E+00  4.8688E+00  1.9816E+00  3.7390E+00  1.8830E+00  1.7657E+00  1.0285E+01
             1.1600E+01
 PARAMETER:  3.1745E-01  1.8258E-01  6.4969E-01  7.6898E-01  1.6829E+00  7.8393E-01  1.4188E+00  7.3287E-01  6.6855E-01  2.4307E+00
             2.5510E+00
 GRADIENT:   4.5580E+00  2.9494E+01 -4.7617E+00  7.0545E+01 -8.1941E+00 -6.5727E-01  4.4017E+00  1.6362E+00  1.1419E+01  1.9889E+01
             4.5760E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -628.164354093528        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      224
 NPARAMETR:  1.0978E+00  5.8088E-01  1.6680E+00  1.6068E+00  6.9399E+00  1.7564E+00  3.7031E+00  1.1898E+00  1.1720E+00  7.5180E+00
             1.0687E+01
 PARAMETER:  1.9332E-01 -4.4321E-01  6.1164E-01  5.7426E-01  2.0373E+00  6.6329E-01  1.4092E+00  2.7380E-01  2.5867E-01  2.1173E+00
             2.4690E+00
 GRADIENT:  -1.3539E+01  1.0748E+01  1.2523E+01  1.9537E+00  6.4461E-01  8.3343E+00  8.7074E+00 -1.1482E-01 -9.7054E+00 -2.4550E+00
             2.3972E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -688.687301269268        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  8.1680E-01  7.9688E-02  1.4226E-01  9.9436E-01  2.8720E+01  2.1379E+00  8.2466E-01  1.0000E-02  1.2255E+00  5.9781E+00
             7.4573E+00
 PARAMETER: -1.0236E-01 -2.4296E+00 -1.8501E+00  9.4347E-02  3.4576E+00  8.5981E-01 -9.2782E-02 -8.7912E+00  3.0334E-01  1.8881E+00
             2.1092E+00
 GRADIENT:   2.8097E+01  4.2222E+00 -1.6320E+01  7.7825E+01  1.2452E-01  1.2334E+01  5.3698E-01  0.0000E+00 -1.6137E+01  2.5511E-02
            -1.6396E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -727.776561461719        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  4.3229E-01  1.1243E-02  2.4495E-02  2.9167E-01  1.4371E+02  1.8258E+00  1.1505E-02  1.0000E-02  1.0196E+00  2.7664E+00
             8.8749E+00
 PARAMETER: -7.3866E-01 -4.3880E+00 -3.6093E+00 -1.1321E+00  5.0678E+00  7.0200E-01 -4.3650E+00 -2.0180E+01  1.1946E-01  1.1175E+00
             2.2832E+00
 GRADIENT:  -1.1562E+01  2.3020E-01  2.4342E+01 -1.0396E+01 -5.9298E-03  8.7585E+00  4.7277E-07  0.0000E+00  4.4139E+00 -1.7814E-05
             1.7036E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -728.148415246721        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  4.2497E-01  1.0239E-02  2.2144E-02  2.7006E-01  1.6607E+02  1.7786E+00  1.0000E-02  1.0000E-02  9.9753E-01  2.7648E+00
             8.8633E+00
 PARAMETER: -7.5573E-01 -4.4816E+00 -3.7102E+00 -1.2091E+00  5.2124E+00  6.7581E-01 -4.5674E+00 -2.0736E+01  9.7527E-02  1.1170E+00
             2.2819E+00
 GRADIENT:  -3.6381E+00  1.0742E-01  8.5649E+00 -1.0211E+01 -1.5920E-03  5.9004E-01  0.0000E+00  0.0000E+00  2.9822E-01 -1.6810E-05
             1.0641E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -728.292425321831        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  4.3924E-01  1.0000E-02  2.3889E-02  2.8866E-01  1.5578E+02  1.7846E+00  1.2000E-02  1.0000E-02  1.0022E+00  2.8944E+00
             8.7470E+00
 PARAMETER: -7.2270E-01 -4.5327E+00 -3.6343E+00 -1.1425E+00  5.1484E+00  6.7921E-01 -4.3229E+00 -2.0130E+01  1.0217E-01  1.1628E+00
             2.2687E+00
 GRADIENT:   6.5196E-01  0.0000E+00 -6.0574E-01  4.1091E-01  1.7119E-03 -1.9735E-01  2.4177E-07  0.0000E+00 -7.5630E-02 -2.2449E-05
            -2.6266E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -728.292901251898        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  4.3889E-01  1.0000E-02  2.3900E-02  2.8866E-01  1.5507E+02  1.7856E+00  1.1999E-02  1.0000E-02  1.0025E+00  2.8945E+00
             8.7482E+00
 PARAMETER: -7.2351E-01 -4.5451E+00 -3.6339E+00 -1.1425E+00  5.1439E+00  6.7977E-01 -4.3229E+00 -2.0130E+01  1.0245E-01  1.1628E+00
             2.2689E+00
 GRADIENT:  -7.2373E-03  0.0000E+00  1.7461E-02 -8.0802E-02  1.4766E-03  7.6951E-04  2.4124E-07  0.0000E+00  1.0353E-02 -2.2543E-05
             1.0499E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -728.303607898835        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1067
 NPARAMETR:  4.4010E-01  1.0000E-02  2.3892E-02  2.8897E-01  1.4832E+01  1.7867E+00  1.2008E-02  1.0000E-02  1.0017E+00  2.9846E+00
             8.7457E+00
 PARAMETER: -7.2075E-01 -4.5451E+00 -3.6342E+00 -1.1414E+00  2.7968E+00  6.8036E-01 -4.3221E+00 -2.0130E+01  1.0169E-01  1.1935E+00
             2.2686E+00
 GRADIENT:   8.4135E-03  0.0000E+00 -4.5817E-01  3.6218E-01 -1.6472E-03 -2.0159E-02  3.5091E-07  0.0000E+00 -3.9319E-02 -2.0793E-03
            -1.5829E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -728.307634198941        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1246
 NPARAMETR:  4.4182E-01  1.0000E-02  2.4187E-02  2.9139E-01  1.9855E+01  1.7874E+00  1.2008E-02  1.0000E-02  1.0012E+00  1.1243E+01
             8.7536E+00
 PARAMETER: -7.1684E-01 -4.5451E+00 -3.6219E+00 -1.1331E+00  3.0885E+00  6.8077E-01 -4.3222E+00 -2.0130E+01  1.0124E-01  2.5197E+00
             2.2695E+00
 GRADIENT:  -5.5520E-02  0.0000E+00  7.0624E-01 -9.2832E-01  1.2536E-02  1.6643E-02  2.8565E-07  0.0000E+00 -3.9644E-02 -3.4255E-03
             4.5051E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -728.309210020695        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1424
 NPARAMETR:  4.4187E-01  1.0000E-02  2.4173E-02  2.9144E-01  1.7727E+01  1.7874E+00  1.2008E-02  1.0000E-02  1.0017E+00  1.1284E+01
             8.7492E+00
 PARAMETER: -7.1673E-01 -4.5451E+00 -3.6225E+00 -1.1329E+00  2.9751E+00  6.8078E-01 -4.3222E+00 -2.0130E+01  1.0174E-01  2.5234E+00
             2.2690E+00
 GRADIENT:   2.7020E-02  0.0000E+00  1.1016E-01 -2.3158E-01  3.9517E-03  1.1555E-03  3.2182E-07  0.0000E+00  1.1988E-03  2.6100E-03
             7.5222E-02

0ITERATION NO.:   56    OBJECTIVE VALUE:  -728.309210020695        NO. OF FUNC. EVALS.:  35
 CUMULATIVE NO. OF FUNC. EVALS.:     1459
 NPARAMETR:  4.4195E-01  1.0000E-02  2.4194E-02  2.9136E-01  1.7740E+01  1.7871E+00  1.2060E-02  1.0000E-02  1.0017E+00  1.1278E+01
             8.7541E+00
 PARAMETER: -7.1673E-01 -4.5451E+00 -3.6225E+00 -1.1329E+00  2.9751E+00  6.8078E-01 -4.3222E+00 -2.0130E+01  1.0174E-01  2.5234E+00
             2.2690E+00
 GRADIENT:  -1.3806E+03  0.0000E+00 -2.7693E+02  8.8170E+02 -3.3340E+02  1.4525E+03 -2.1469E-05  0.0000E+00  3.0966E-03  1.9732E+02
            -4.3683E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1459
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0026E-03  1.6645E-06  1.0437E-04 -2.2728E-02 -8.6621E-04
 SE:             2.9136E-02  1.0561E-06  2.0944E-04  2.4280E-02  6.8106E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7255E-01  1.1502E-01  6.1826E-01  3.4924E-01  2.0342E-01

 ETASHRINKSD(%)  2.3902E+00  9.9996E+01  9.9298E+01  1.8659E+01  9.7718E+01
 ETASHRINKVR(%)  4.7233E+00  1.0000E+02  9.9995E+01  3.3837E+01  9.9948E+01
 EBVSHRINKSD(%)  2.2924E+00  9.9995E+01  9.9340E+01  1.8884E+01  9.7926E+01
 EBVSHRINKVR(%)  4.5323E+00  1.0000E+02  9.9996E+01  3.4202E+01  9.9957E+01
 RELATIVEINF(%)  3.1813E+00  6.0115E-08  3.2142E-05  4.5782E-01  3.1414E-03
 EPSSHRINKSD(%)  1.5040E+01
 EPSSHRINKVR(%)  2.7818E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -728.30921002069454     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       6.8416165430436422     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.91
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -728.309       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.42E-01  1.00E-02  2.42E-02  2.91E-01  1.77E+01  1.79E+00  1.20E-02  1.00E-02  1.00E+00  1.13E+01  8.75E+00
 


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
+        2.47E+06
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.12E+04  0.00E+00  3.39E+07
 
 TH 4
+        2.75E+03  0.00E+00 -2.26E+05  2.32E+06
 
 TH 5
+       -1.82E+01  0.00E+00  1.79E+02 -1.44E+04  8.93E+01
 
 TH 6
+       -6.41E+05  0.00E+00  5.88E+03 -1.61E+03  1.01E+01  5.41E+01
 
 TH 7
+        2.98E+00  0.00E+00  1.28E+01  2.09E+00  4.07E-03 -1.42E+00 -4.95E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.78E+03  0.00E+00  7.55E+03 -1.82E+03  1.06E+01 -4.66E+02  1.83E+00  0.00E+00  8.86E+01
 
 TH10
+        3.37E+01  0.00E+00 -1.11E+03  2.67E+04 -1.89E+00 -1.87E+01 -6.33E-03  0.00E+00 -1.96E+01  3.09E+02
 
 TH11
+       -6.58E+01  0.00E+00  5.62E+02 -1.03E+02  5.61E-01  2.75E+01 -2.42E-02  0.00E+00  3.13E+01 -1.04E+00  6.34E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.064
Stop Time:
Sat Sep 25 14:24:16 CDT 2021
