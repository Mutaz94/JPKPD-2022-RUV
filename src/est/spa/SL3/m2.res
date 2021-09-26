Sat Sep 25 11:28:01 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat2.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1631.39626426571        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.1536E+01 -7.6759E+01 -4.3062E+01 -6.0314E+01  2.8501E+01  7.7482E+00  8.4064E+00  1.4505E+01  2.8753E+01  1.2406E+01
            -1.0628E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1645.94461682893        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.7784E-01  1.1685E+00  1.3371E+00  9.5346E-01  1.1612E+00  9.6218E-01  8.8071E-01  8.4837E-01  7.3868E-01  8.7625E-01
             1.3808E+00
 PARAMETER:  7.7591E-02  2.5575E-01  3.9053E-01  5.2344E-02  2.4947E-01  6.1441E-02 -2.7028E-02 -6.4438E-02 -2.0289E-01 -3.2104E-02
             4.2265E-01
 GRADIENT:   7.0149E+00  1.9491E+00  1.2205E+01 -1.7276E+01  2.5784E+00 -5.3309E+00 -1.0194E+01 -2.0338E+00 -1.9178E+01 -2.2155E+01
             2.1810E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1648.96847986857        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.9574E-01  1.0541E+00  1.7425E+00  1.0482E+00  1.2464E+00  9.5974E-01  8.6774E-01  5.7786E-01  8.4618E-01  1.0793E+00
             1.3559E+00
 PARAMETER:  9.5727E-02  1.5265E-01  6.5534E-01  1.4712E-01  3.2026E-01  5.8902E-02 -4.1857E-02 -4.4842E-01 -6.7029E-02  1.7635E-01
             4.0449E-01
 GRADIENT:   5.4061E+01  1.2022E+01  8.0575E+00  1.0331E+01  4.9574E+00 -6.5177E+00  2.2466E+00 -6.0370E-01 -1.6013E+00 -1.0565E+01
             1.7125E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1652.24534029980        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.7318E-01  1.0771E+00  1.0976E+00  1.0098E+00  1.0545E+00  9.7707E-01  9.0485E-01  2.5432E-01  8.7012E-01  1.0091E+00
             1.2502E+00
 PARAMETER:  7.2810E-02  1.7428E-01  1.9310E-01  1.0971E-01  1.5302E-01  7.6808E-02  1.4669E-05 -1.2692E+00 -3.9119E-02  1.0905E-01
             3.2332E-01
 GRADIENT:   1.7785E+00 -5.9710E-01 -2.0031E+00  2.0028E+00  5.4474E+00 -2.8098E-01 -1.3529E-01 -7.7165E-02  8.0810E-01  3.2249E-01
             2.0769E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1652.25206353127        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.7265E-01  1.0689E+00  1.0821E+00  1.0138E+00  1.0413E+00  9.7753E-01  9.2264E-01  2.6273E-01  8.5909E-01  9.9586E-01
             1.2462E+00
 PARAMETER:  7.2272E-02  1.6664E-01  1.7894E-01  1.1372E-01  1.4046E-01  7.7278E-02  1.9482E-02 -1.2366E+00 -5.1876E-02  9.5849E-02
             3.2008E-01
 GRADIENT:   6.6342E-01 -1.1542E-01 -9.3357E-01  1.1107E+00  2.7790E+00 -1.4616E-01 -1.5602E-01 -7.3886E-02  2.4396E-01  6.5826E-02
             9.8624E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1652.41730146852        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      460
 NPARAMETR:  9.8424E-01  1.0939E+00  1.0638E+00  1.0000E+00  1.0407E+00  9.8345E-01  9.1645E-01  2.8469E-01  8.6366E-01  9.8714E-01
             1.2447E+00
 PARAMETER:  8.4113E-02  1.8975E-01  1.6189E-01  1.0002E-01  1.3992E-01  8.3315E-02  1.2755E-02 -1.1564E+00 -4.6577E-02  8.7058E-02
             3.1886E-01
 GRADIENT:   1.5805E+00  6.2720E-01  1.0498E+00 -7.7338E-02 -1.8070E+00 -7.5757E-02 -1.8761E-01 -8.9534E-02 -3.8362E-01 -4.7506E-01
            -3.6001E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1652.45948365248        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      637
 NPARAMETR:  9.8403E-01  1.2179E+00  1.0219E+00  9.1967E-01  1.0835E+00  9.8477E-01  8.4274E-01  3.8902E-01  9.2086E-01  1.0014E+00
             1.2415E+00
 PARAMETER:  8.3899E-02  2.9717E-01  1.2167E-01  1.6259E-02  1.8024E-01  8.4651E-02 -7.1101E-02 -8.4413E-01  1.7553E-02  1.0141E-01
             3.1632E-01
 GRADIENT:  -2.8118E-01 -1.2282E+00 -1.0638E-01 -1.6358E+00  5.0500E-01  1.4688E-01 -3.1226E-01 -5.8897E-02 -7.5502E-01 -3.3254E-01
            -5.2082E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1652.48541730422        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      813
 NPARAMETR:  9.8341E-01  1.1637E+00  1.1111E+00  9.6052E-01  1.0915E+00  9.8329E-01  8.4974E-01  5.8135E-01  9.1058E-01  1.0151E+00
             1.2437E+00
 PARAMETER:  8.3267E-02  2.5161E-01  2.0539E-01  5.9722E-02  1.8751E-01  8.3144E-02 -6.2829E-02 -4.4240E-01  6.3223E-03  1.1495E-01
             3.1812E-01
 GRADIENT:  -4.5147E-01  4.2084E+00  9.2462E-01  3.5060E+00 -3.1252E+00 -2.1965E-01  1.0499E-01 -2.1111E-03  2.9069E-01  1.5008E-01
             7.7680E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1652.50940979587        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  9.8271E-01  1.0814E+00  1.2154E+00  1.0147E+00  1.0973E+00  9.8307E-01  8.7101E-01  6.9205E-01  8.8223E-01  1.0283E+00
             1.2410E+00
 PARAMETER:  8.2559E-02  1.7829E-01  2.9507E-01  1.1454E-01  1.9287E-01  8.2923E-02 -3.8105E-02 -2.6809E-01 -2.5300E-02  1.2793E-01
             3.1592E-01
 GRADIENT:  -1.6828E-01  1.9050E+00  9.0606E-01  1.4229E+00 -1.2321E+00  1.0286E-02  7.4719E-02 -9.5684E-02  4.6326E-02 -1.4076E-01
            -3.6862E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1652.51993144684        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  9.8208E-01  9.9140E-01  1.2991E+00  1.0733E+00  1.0914E+00  9.8241E-01  8.9933E-01  7.6615E-01  8.5390E-01  1.0352E+00
             1.2406E+00
 PARAMETER:  8.1915E-02  9.1364E-02  3.6168E-01  1.7070E-01  1.8748E-01  8.2249E-02 -6.1020E-03 -1.6638E-01 -5.7942E-02  1.3461E-01
             3.1563E-01
 GRADIENT:   1.6663E-01 -4.4602E-02 -1.0405E-01  4.7586E-02  1.4837E-01  5.3549E-02  2.2181E-02  2.6336E-03  1.6202E-01  1.5421E-01
            -6.6283E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1652.52169809532        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1340
 NPARAMETR:  9.8156E-01  9.3925E-01  1.3418E+00  1.1076E+00  1.0852E+00  9.8182E-01  9.2699E-01  7.9483E-01  8.3475E-01  1.0348E+00
             1.2408E+00
 PARAMETER:  8.1392E-02  3.7327E-02  3.9398E-01  2.0218E-01  1.8174E-01  8.1652E-02  2.4185E-02 -1.2963E-01 -8.0627E-02  1.3419E-01
             3.1576E-01
 GRADIENT:   2.2595E-02  1.1886E-01  1.4136E-02  1.6965E-01 -4.6601E-02  1.4406E-03  3.4049E-02 -1.1760E-04  2.3334E-02  1.3168E-02
             8.0647E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1652.52177360340        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1515
 NPARAMETR:  9.8145E-01  9.2922E-01  1.3549E+00  1.1142E+00  1.0859E+00  9.8171E-01  9.2634E-01  8.0685E-01  8.3294E-01  1.0364E+00
             1.2408E+00
 PARAMETER:  8.1275E-02  2.6588E-02  4.0376E-01  2.0814E-01  1.8237E-01  8.1543E-02  2.3488E-02 -1.1462E-01 -8.2792E-02  1.3579E-01
             3.1573E-01
 GRADIENT:   4.1261E-04 -3.6681E-03 -2.2021E-03 -3.5601E-03  2.4063E-03  6.6805E-05 -5.3372E-04  1.4909E-04 -7.2207E-05  6.4430E-04
             1.1960E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1652.52177449210        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1685
 NPARAMETR:  9.8145E-01  9.2964E-01  1.3546E+00  1.1139E+00  1.0859E+00  9.8172E-01  9.2631E-01  8.0660E-01  8.3304E-01  1.0364E+00
             1.2408E+00
 PARAMETER:  8.1279E-02  2.7047E-02  4.0348E-01  2.0789E-01  1.8240E-01  8.1547E-02  2.3451E-02 -1.1493E-01 -8.2672E-02  1.3575E-01
             3.1573E-01
 GRADIENT:   8.6323E-04  1.3441E-03  1.3118E-03  3.9781E-04 -9.2954E-04 -7.8373E-05  1.8892E-03  1.2400E-04  2.3654E-03  8.0753E-05
             1.6543E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1685
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.8512E-04 -1.2753E-02 -1.5829E-02 -1.4911E-03 -3.2981E-02
 SE:             2.9730E-02  1.5578E-02  9.5653E-03  2.4791E-02  2.2290E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9503E-01  4.1299E-01  9.7962E-02  9.5204E-01  1.3897E-01

 ETASHRINKSD(%)  4.0042E-01  4.7812E+01  6.7955E+01  1.6948E+01  2.5325E+01
 ETASHRINKVR(%)  7.9923E-01  7.2764E+01  8.9731E+01  3.1023E+01  4.4237E+01
 EBVSHRINKSD(%)  6.5713E-01  4.7552E+01  7.0124E+01  1.7470E+01  2.3608E+01
 EBVSHRINKVR(%)  1.3099E+00  7.2493E+01  9.1074E+01  3.1888E+01  4.1643E+01
 RELATIVEINF(%)  9.7132E+01  3.5968E-01  1.0173E+00  1.0209E+00  6.7310E+00
 EPSSHRINKSD(%)  4.0464E+01
 EPSSHRINKVR(%)  6.4555E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1652.5217744920960     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -917.37094792835785     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.21
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1652.522       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  9.30E-01  1.35E+00  1.11E+00  1.09E+00  9.82E-01  9.26E-01  8.07E-01  8.33E-01  1.04E+00  1.24E+00
 


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
+        1.18E+03
 
 TH 2
+       -2.02E+01  4.21E+02
 
 TH 3
+        5.34E+00  6.37E+01  8.91E+01
 
 TH 4
+       -1.66E+01  5.35E+02 -3.84E+01  8.77E+02
 
 TH 5
+       -1.79E+00 -1.85E+02 -1.75E+02  2.40E+01  5.15E+02
 
 TH 6
+       -2.90E-01 -3.87E+00  2.05E+00 -5.25E+00 -2.26E+00  1.99E+02
 
 TH 7
+        3.45E+00 -1.45E+00  5.26E+00 -4.49E+00 -6.66E+00  1.14E+00  1.52E+01
 
 TH 8
+       -3.03E+00 -4.13E+00 -1.33E+01 -2.36E-01  2.28E+00 -1.28E+00 -3.54E-01  5.68E+00
 
 TH 9
+        3.34E+00 -1.35E+01 -5.02E-01  1.31E+01  6.73E-01  4.53E-01  4.01E+01  3.44E+00  1.31E+02
 
 TH10
+        2.57E+00  3.50E+00 -7.25E+00 -6.24E+00 -5.62E+01  7.22E-02  4.08E+00  8.68E+00  1.92E+00  7.30E+01
 
 TH11
+       -1.06E+01 -1.55E+01 -1.71E+01 -6.28E+00 -6.93E-02  1.95E+00  3.90E+00  8.19E+00  1.28E+01  1.68E+01  1.50E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.069
Stop Time:
Sat Sep 25 11:28:28 CDT 2021
