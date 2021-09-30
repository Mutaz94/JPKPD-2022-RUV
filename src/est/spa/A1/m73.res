Wed Sep 29 12:19:18 CDT 2021
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
$DATA ../../../../data/spa/A1/dat73.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m73.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1204.15108461815        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8180E+02  1.0685E+01  3.7855E+01 -6.7008E+00 -2.5553E+01  7.4649E+01 -3.7010E+00 -1.0052E+01 -2.5001E+00  9.3391E+00
            -9.5101E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1479.18677876230        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1842E+00  1.1041E+00  1.0862E+00  1.0134E+00  1.1178E+00  1.0758E+00  9.7646E-01  9.7423E-01  8.9985E-01  7.9333E-01
             2.3067E+00
 PARAMETER:  2.6910E-01  1.9903E-01  1.8267E-01  1.1331E-01  2.1137E-01  1.7306E-01  7.6181E-02  7.3889E-02 -5.5250E-03 -1.3152E-01
             9.3583E-01
 GRADIENT:   4.7946E+02  1.1274E+01  3.6583E+00  1.1385E+01 -7.8057E+00  4.0731E+01 -7.6881E-04 -1.0399E+00  6.1505E-01  1.0688E+01
             9.3519E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1485.39466838224        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.1447E+00  9.9741E-01  8.1621E-01  1.0552E+00  8.9235E-01  9.7229E-01  1.3872E+00  1.3698E+00  9.1626E-01  3.5764E-01
             2.2301E+00
 PARAMETER:  2.3513E-01  9.7403E-02 -1.0308E-01  1.5369E-01 -1.3898E-02  7.1899E-02  4.2730E-01  4.1466E-01  1.2544E-02 -9.2821E-01
             9.0207E-01
 GRADIENT:   4.4196E+02  8.8592E+00 -1.2559E+01  2.4983E+01 -2.1036E+01  2.0145E+01  2.7232E+01  2.1892E+01  1.9733E+01  4.3872E+00
             1.0921E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1508.87684750542        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0050E+00  1.0121E+00  9.6109E-01  1.0328E+00  9.9124E-01  8.3334E-01  1.1993E+00  1.0894E+00  8.2525E-01  3.6621E-01
             2.2096E+00
 PARAMETER:  1.0501E-01  1.1201E-01  6.0313E-02  1.3229E-01  9.1205E-02 -8.2310E-02  2.8171E-01  1.8561E-01 -9.2066E-02 -9.0455E-01
             8.9279E-01
 GRADIENT:   3.2062E+01  2.1775E+00  2.1431E+00  3.4709E+00 -9.5097E+00 -3.8317E+00  5.0531E+00  1.9199E+00  9.1237E-01  2.2637E+00
            -9.1035E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1510.67804274188        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      365
 NPARAMETR:  1.0239E+00  1.1192E+00  1.0637E+00  9.8815E-01  1.1050E+00  8.6685E-01  1.0254E+00  1.1925E+00  8.8933E-01  3.3280E-01
             2.3221E+00
 PARAMETER:  1.2359E-01  2.1258E-01  1.6179E-01  8.8081E-02  1.9981E-01 -4.2887E-02  1.2511E-01  2.7603E-01 -1.7285E-02 -1.0002E+00
             9.4246E-01
 GRADIENT:  -4.0601E+00 -1.1309E+00 -5.1249E-02  3.0326E-01 -2.5628E+00  6.3930E+00 -4.6105E-01 -4.1542E-01  6.4726E-01  1.2231E+00
             3.5003E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1511.38776623754        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      540
 NPARAMETR:  1.0255E+00  1.2732E+00  1.2095E+00  8.9750E-01  1.2407E+00  8.5111E-01  8.9599E-01  1.6337E+00  9.6769E-01  1.4573E-01
             2.3413E+00
 PARAMETER:  1.2519E-01  3.4157E-01  2.9024E-01 -8.1473E-03  3.1564E-01 -6.1217E-02 -9.8207E-03  5.9082E-01  6.7158E-02 -1.8260E+00
             9.5071E-01
 GRADIENT:   2.3708E-01  1.3943E+00  1.4221E+00  8.7161E-01 -1.7583E+00 -1.4716E-02 -2.0092E-01 -3.6354E-01  3.6894E-02  1.4503E-01
             8.0531E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1511.78533434506        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      719
 NPARAMETR:  1.0273E+00  1.6486E+00  7.7807E-01  6.5429E-01  1.2544E+00  8.5497E-01  7.7705E-01  1.6544E+00  1.1712E+00  1.0000E-02
             2.3222E+00
 PARAMETER:  1.2696E-01  5.9991E-01 -1.5094E-01 -3.2420E-01  3.2669E-01 -5.6685E-02 -1.5225E-01  6.0341E-01  2.5805E-01 -5.3477E+00
             9.4250E-01
 GRADIENT:  -2.3356E+00  2.2249E+01  3.4715E+00  8.4707E+00 -1.0062E+01  5.3606E-01 -9.4636E-01 -4.2090E-02 -2.3935E+00  0.0000E+00
            -1.5286E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1512.63732629454        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      896
 NPARAMETR:  1.0274E+00  1.9560E+00  3.0856E-01  4.3213E-01  1.1772E+00  8.6122E-01  7.1588E-01  8.2454E-01  1.5041E+00  1.0000E-02
             2.2794E+00
 PARAMETER:  1.2702E-01  7.7093E-01 -1.0759E+00 -7.3902E-01  2.6311E-01 -4.9407E-02 -2.3424E-01 -9.2924E-02  5.0823E-01 -1.1129E+01
             9.2390E-01
 GRADIENT:  -4.4561E+00  2.8259E+01  1.2019E+00  8.0393E+00 -1.1691E+01  2.1374E+00 -4.8042E-01  1.6678E-01 -2.6586E-01  0.0000E+00
            -3.9765E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1513.23307713314        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1081
 NPARAMETR:  1.0290E+00  2.0302E+00  2.7492E-01  3.6939E-01  1.2396E+00  8.5503E-01  6.9548E-01  5.4307E-01  1.7018E+00  1.0000E-02
             2.3063E+00
 PARAMETER:  1.2856E-01  8.0814E-01 -1.1913E+00 -8.9591E-01  3.1483E-01 -5.6615E-02 -2.6316E-01 -5.1052E-01  6.3171E-01 -1.3012E+01
             9.3564E-01
 GRADIENT:   1.5352E+00 -5.9458E+00 -1.5782E-01  8.2356E-01 -2.3691E+00  1.3157E-01  7.6635E-01  4.6299E-02  5.0238E-01  0.0000E+00
             1.5988E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1513.46498295216        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1259
 NPARAMETR:  1.0264E+00  2.1547E+00  2.4071E-01  3.0130E-01  1.3456E+00  8.5455E-01  6.6791E-01  1.7694E-02  2.0248E+00  1.0000E-02
             2.3104E+00
 PARAMETER:  1.2609E-01  8.6767E-01 -1.3242E+00 -1.0997E+00  3.9686E-01 -5.7184E-02 -3.0359E-01 -3.9345E+00  8.0548E-01 -1.3012E+01
             9.3744E-01
 GRADIENT:  -7.4400E+00  5.9445E+00 -2.0689E+00  7.8923E+00  7.6975E+00 -4.8489E-01  1.3316E-01  2.0504E-05  1.6863E+00  0.0000E+00
            -2.9971E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1513.75244388974        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1436
 NPARAMETR:  1.0291E+00  2.2915E+00  1.9103E-01  2.0766E-01  1.4513E+00  8.5456E-01  6.3758E-01  1.0000E-02  2.4977E+00  1.0000E-02
             2.3446E+00
 PARAMETER:  1.2866E-01  9.2921E-01 -1.5553E+00 -1.4718E+00  4.7249E-01 -5.7172E-02 -3.5007E-01 -7.7441E+00  1.0154E+00 -1.3012E+01
             9.5212E-01
 GRADIENT:   8.8003E-01  3.9940E+00  6.8783E-01  1.0991E+00  2.4166E+00  4.6783E-02 -4.6608E-01  0.0000E+00 -8.2063E-02  0.0000E+00
             4.3792E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1513.80973376193        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1621
 NPARAMETR:  1.0287E+00  2.3233E+00  1.5729E-01  1.8229E-01  1.4641E+00  8.5417E-01  6.3312E-01  1.0000E-02  2.6562E+00  1.0000E-02
             2.3448E+00
 PARAMETER:  1.2830E-01  9.4301E-01 -1.7497E+00 -1.6022E+00  4.8122E-01 -5.7626E-02 -3.5710E-01 -8.9383E+00  1.0769E+00 -1.3012E+01
             9.5218E-01
 GRADIENT:   5.6000E-01 -7.4023E+00 -1.9464E-01  1.3419E+00  6.0572E+00 -1.1926E-01 -1.5374E-01  0.0000E+00  1.3029E+00  0.0000E+00
             1.3272E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1515.04152240365        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1803
 NPARAMETR:  1.0256E+00  2.3452E+00  5.0430E-02  1.4730E-01  1.3136E+00  8.5572E-01  6.4325E-01  1.0000E-02  2.3192E+00  1.0000E-02
             2.3011E+00
 PARAMETER:  1.2526E-01  9.5237E-01 -2.8872E+00 -1.8153E+00  3.7278E-01 -5.5810E-02 -3.4121E-01 -8.9383E+00  9.4123E-01 -1.3012E+01
             9.3339E-01
 GRADIENT:   3.2772E+00 -1.5832E+00  1.4776E-01  1.9849E+00 -4.3096E+00  9.1695E-01  1.7656E+00  0.0000E+00 -1.1394E+00  0.0000E+00
            -4.3008E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1516.09145812996        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1985
 NPARAMETR:  1.0238E+00  2.4256E+00  3.3288E-02  9.9034E-02  1.4199E+00  8.4645E-01  6.1715E-01  1.0000E-02  2.8113E+00  1.0000E-02
             2.3775E+00
 PARAMETER:  1.2351E-01  9.8606E-01 -3.3026E+00 -2.2123E+00  4.5060E-01 -6.6702E-02 -3.8265E-01 -8.9383E+00  1.1337E+00 -1.3012E+01
             9.6605E-01
 GRADIENT:  -6.2356E+00 -1.5621E+01 -2.4620E-01  1.3545E+00  2.2268E+00 -1.9456E+00 -2.4684E+00  0.0000E+00 -2.9765E+00  0.0000E+00
             9.4474E+00

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1516.09145812996        NO. OF FUNC. EVALS.:  36
 CUMULATIVE NO. OF FUNC. EVALS.:     2021
 NPARAMETR:  1.0238E+00  2.4244E+00  3.3341E-02  9.8929E-02  1.4199E+00  8.4730E-01  6.1726E-01  1.0000E-02  2.8098E+00  1.0000E-02
             2.3775E+00
 PARAMETER:  1.2351E-01  9.8606E-01 -3.3026E+00 -2.2123E+00  4.5060E-01 -6.6702E-02 -3.8265E-01 -8.9383E+00  1.1337E+00 -1.3012E+01
             9.6605E-01
 GRADIENT:  -8.8373E+03  2.2038E+03 -6.4668E+02  4.8996E+02  6.6655E+00 -1.8667E+00 -2.8530E+03  0.0000E+00  1.8848E+03  0.0000E+00
            -6.5727E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2021
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3427E-03 -1.2377E-02  3.5434E-06  2.7008E-02 -1.8360E-04
 SE:             2.9291E-02  2.6560E-02  2.0494E-05  1.4768E-02  1.6276E-04
 N:                     100         100         100         100         100

 P VAL.:         9.0914E-01  6.4121E-01  8.6273E-01  6.7440E-02  2.5932E-01

 ETASHRINKSD(%)  1.8704E+00  1.1021E+01  9.9931E+01  5.0524E+01  9.9455E+01
 ETASHRINKVR(%)  3.7057E+00  2.0828E+01  1.0000E+02  7.5521E+01  9.9997E+01
 EBVSHRINKSD(%)  2.5840E+00  1.2632E+01  9.9919E+01  5.3499E+01  9.9400E+01
 EBVSHRINKVR(%)  5.1013E+00  2.3668E+01  1.0000E+02  7.8376E+01  9.9996E+01
 RELATIVEINF(%)  9.2473E+01  2.2874E+01  3.7612E-05  4.8400E+00  9.2996E-04
 EPSSHRINKSD(%)  2.8268E+01
 EPSSHRINKVR(%)  4.8545E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1516.0914581299621     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -780.94063156622394     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.90
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.33
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1516.091       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.43E+00  3.33E-02  9.90E-02  1.42E+00  8.46E-01  6.17E-01  1.00E-02  2.81E+00  1.00E-02  2.38E+00
 


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
+        6.82E+06
 
 TH 2
+        8.13E-01  9.91E+03
 
 TH 3
+       -1.24E+03  2.12E+05  4.33E+06
 
 TH 4
+       -1.98E+06  2.83E+02  3.34E+04  1.12E+06
 
 TH 5
+        6.74E+05  3.56E+04  7.91E+05 -3.92E+05  1.34E+05
 
 TH 6
+       -7.52E+02  3.00E+01 -8.69E+02  4.09E+02 -4.67E+00  2.50E+02
 
 TH 7
+       -2.42E+02  9.65E+04 -3.05E+02  1.35E+02  3.61E+05 -4.06E+02  9.78E+05
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.38E+05 -6.81E+00  5.83E+03  7.91E+04 -2.73E+04  3.02E+01  1.17E+01  0.00E+00  1.08E+04
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.87E+05  9.88E+03 -2.20E+05  1.09E+05 -2.70E+02  5.34E+00 -1.00E+05  0.00E+00  7.57E+03  0.00E+00  2.07E+04
 
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
 #CPUT: Total CPU Time in Seconds,       32.290
Stop Time:
Wed Sep 29 12:19:52 CDT 2021
