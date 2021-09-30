Wed Sep 29 21:12:24 CDT 2021
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
$DATA ../../../../data/spa1/B/dat54.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1667.52869715815        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8150E+02  5.3637E+01  3.9170E+01  1.0739E+02  1.2188E+01  3.3826E+01  1.3869E+00 -1.8446E+02 -2.9844E+01  1.0182E+01
            -7.2876E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2078.16002167164        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  1.0123E+00  1.0053E+00  9.9440E-01  8.7690E-01  1.0321E+00  8.7526E-01  1.0338E+00  1.0878E+00  1.0803E+00  9.9274E-01
             1.6232E+00
 PARAMETER:  1.1219E-01  1.0524E-01  9.4379E-02 -3.1358E-02  1.3158E-01 -3.3235E-02  1.3322E-01  1.8418E-01  1.7728E-01  9.2714E-02
             5.8439E-01
 GRADIENT:   9.3428E+01 -8.7546E+01 -3.0848E+01 -8.4719E+01  4.3816E+00 -7.5845E+01  1.4139E+00  1.5923E+01  1.0638E+01  2.4866E+01
             2.9419E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2124.77421285911        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  9.3821E-01  1.1730E+00  1.0399E+00  6.8309E-01  1.1789E+00  9.4704E-01  1.2717E+00  1.7823E-02  1.5486E+00  1.0216E+00
             1.2027E+00
 PARAMETER:  3.6216E-02  2.5954E-01  1.3909E-01 -2.8112E-01  2.6455E-01  4.5589E-02  3.4036E-01 -3.9273E+00  5.3737E-01  1.2137E-01
             2.8456E-01
 GRADIENT:   8.0165E+01  3.4599E+01  1.3868E+01  1.5892E+01  8.8400E+00 -2.0917E+01  4.9263E+01  9.7894E-04  7.0304E+01  8.0616E+00
             1.7120E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2181.96251460430        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      277
 NPARAMETR:  1.0145E+00  1.4289E+00  8.5755E-01  6.7446E-01  1.1660E+00  9.4078E-01  9.1730E-01  1.2872E-01  1.4880E+00  1.0555E+00
             1.0721E+00
 PARAMETER:  1.1440E-01  4.5693E-01 -5.3673E-02 -2.9384E-01  2.5358E-01  3.8951E-02  1.3681E-02 -1.9501E+00  4.9746E-01  1.5405E-01
             1.6960E-01
 GRADIENT:  -1.1025E+02 -2.4494E+01  5.1171E+00 -1.7156E+01 -2.3215E+01 -7.6357E+01  1.5985E+01  9.3431E-02  2.3297E+01  6.4070E+00
             1.0130E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2199.79323081822        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  1.0596E+00  1.5005E+00  7.8333E-01  6.5095E-01  1.1870E+00  1.1039E+00  8.0721E-01  8.0090E-02  1.3799E+00  1.0659E+00
             9.2791E-01
 PARAMETER:  1.5792E-01  5.0579E-01 -1.4420E-01 -3.2933E-01  2.7145E-01  1.9884E-01 -1.1418E-01 -2.4246E+00  4.2203E-01  1.6379E-01
             2.5175E-02
 GRADIENT:   4.1001E+00  5.9786E-01  1.8093E-01  2.8737E+00  9.1539E-01  2.9223E+00  4.6262E-02  3.7670E-02 -8.8764E-01 -1.6503E-01
            -3.9635E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2199.83966339710        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      634             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0613E+00  1.5019E+00  7.8042E-01  6.4598E-01  1.1866E+00  1.1011E+00  8.0240E-01  1.1099E-02  1.3953E+00  1.0663E+00
             9.3224E-01
 PARAMETER:  1.5949E-01  5.0670E-01 -1.4792E-01 -3.3699E-01  2.7111E-01  1.9628E-01 -1.2015E-01 -4.4009E+00  4.3308E-01  1.6420E-01
             2.9831E-02
 GRADIENT:   1.0088E+03  6.8653E+02  1.5447E+00  1.7795E+02  1.8890E+01  1.9480E+02  9.2820E+00  1.2148E-03  3.6557E+01  1.5579E+00
             1.0203E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2199.83985399166        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  1.0613E+00  1.5012E+00  7.8066E-01  6.4630E-01  1.1861E+00  1.1011E+00  8.0294E-01  1.0000E-02  1.3950E+00  1.0663E+00
             9.3222E-01
 PARAMETER:  1.5948E-01  5.0626E-01 -1.4761E-01 -3.3650E-01  2.7070E-01  1.9628E-01 -1.1947E-01 -5.3557E+00  4.3293E-01  1.6418E-01
             2.9816E-02
 GRADIENT:   6.9710E+00 -5.6451E+00  2.1539E-01 -8.1408E-01 -5.9855E-01  1.8478E+00  2.9741E-01  0.0000E+00  4.3624E-01  1.8996E-01
             2.3209E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2199.84107391171        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1009             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0613E+00  1.4989E+00  7.8071E-01  6.4772E-01  1.1856E+00  1.1011E+00  8.0038E-01  1.0000E-02  1.3917E+00  1.0638E+00
             9.3215E-01
 PARAMETER:  1.5947E-01  5.0475E-01 -1.4755E-01 -3.3430E-01  2.7028E-01  1.9627E-01 -1.2267E-01 -5.3557E+00  4.3050E-01  1.6187E-01
             2.9733E-02
 GRADIENT:   1.0087E+03  6.8107E+02  1.2096E+00  1.7748E+02  1.9928E+01  1.9478E+02  8.9087E+00  0.0000E+00  3.6118E+01  1.1486E+00
             9.1676E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2199.84199671059        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     1204
 NPARAMETR:  1.0613E+00  1.4979E+00  7.8118E-01  6.4832E-01  1.1850E+00  1.1011E+00  8.0154E-01  1.0000E-02  1.3899E+00  1.0641E+00
             9.3216E-01
 PARAMETER:  1.5947E-01  5.0408E-01 -1.4695E-01 -3.3337E-01  2.6978E-01  1.9627E-01 -1.2122E-01 -5.3557E+00  4.2920E-01  1.6218E-01
             2.9749E-02
 GRADIENT:   6.9544E+00 -6.0558E+00 -9.5927E-02 -7.1353E-01  3.9673E-01  1.8485E+00 -1.9263E-01  0.0000E+00  4.1050E-02 -1.1937E-01
            -4.5808E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2199.84265312940        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1394             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0613E+00  1.4972E+00  7.8150E-01  6.4878E-01  1.1845E+00  1.1011E+00  8.0332E-01  1.0000E-02  1.3899E+00  1.0644E+00
             9.3218E-01
 PARAMETER:  1.5947E-01  5.0362E-01 -1.4654E-01 -3.3266E-01  2.6934E-01  1.9627E-01 -1.1900E-01 -5.3557E+00  4.2926E-01  1.6243E-01
             2.9774E-02
 GRADIENT:   1.0085E+03  6.7794E+02  1.3648E+00  1.7712E+02  1.9169E+01  1.9472E+02  8.9711E+00  0.0000E+00  3.6203E+01  1.4588E+00
             1.0111E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2199.84307019714        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1583
 NPARAMETR:  1.0613E+00  1.4957E+00  7.8215E-01  6.4976E-01  1.1837E+00  1.1010E+00  8.0572E-01  1.0000E-02  1.3887E+00  1.0644E+00
             9.3221E-01
 PARAMETER:  1.5947E-01  5.0260E-01 -1.4571E-01 -3.3116E-01  2.6862E-01  1.9626E-01 -1.1602E-01 -5.3557E+00  4.2838E-01  1.6243E-01
             2.9799E-02
 GRADIENT:   6.9648E+00 -5.8390E+00  4.2733E-02 -8.1581E-01 -3.3356E-01  1.8453E+00  2.8304E-01  0.0000E+00  4.4877E-01  2.1125E-01
             7.0674E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2199.84445833941        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1777             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0613E+00  1.4940E+00  7.8267E-01  6.5092E-01  1.1830E+00  1.1010E+00  8.0434E-01  1.0000E-02  1.3863E+00  1.0632E+00
             9.3215E-01
 PARAMETER:  1.5946E-01  5.0142E-01 -1.4504E-01 -3.2937E-01  2.6806E-01  1.9626E-01 -1.1774E-01 -5.3557E+00  4.2666E-01  1.6131E-01
             2.9742E-02
 GRADIENT:   1.0084E+03  6.7208E+02  1.3498E+00  1.7654E+02  1.9137E+01  1.9465E+02  8.7784E+00  0.0000E+00  3.6015E+01  1.4010E+00
             9.9531E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2199.84574277311        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1969
 NPARAMETR:  1.0612E+00  1.4916E+00  7.8364E-01  6.5249E-01  1.1821E+00  1.1010E+00  8.0598E-01  1.0000E-02  1.3828E+00  1.0623E+00
             9.3214E-01
 PARAMETER:  1.5945E-01  4.9983E-01 -1.4381E-01 -3.2697E-01  2.6730E-01  1.9625E-01 -1.1569E-01 -5.3557E+00  4.2413E-01  1.6048E-01
             2.9727E-02
 GRADIENT:   6.9480E+00 -5.7756E+00 -3.4147E-02 -7.1390E-01  1.8075E-01  1.8440E+00 -2.9998E-02  0.0000E+00  1.5323E-01 -3.0125E-02
            -1.3144E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2199.84622920212        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2163             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0612E+00  1.4903E+00  7.8407E-01  6.5327E-01  1.1814E+00  1.1010E+00  8.0773E-01  1.0000E-02  1.3825E+00  1.0625E+00
             9.3216E-01
 PARAMETER:  1.5944E-01  4.9901E-01 -1.4325E-01 -3.2576E-01  2.6667E-01  1.9624E-01 -1.1353E-01 -5.3557E+00  4.2390E-01  1.6060E-01
             2.9751E-02
 GRADIENT:   1.0081E+03  6.6550E+02  1.3704E+00  1.7591E+02  1.8887E+01  1.9456E+02  8.6864E+00  0.0000E+00  3.5980E+01  1.5259E+00
             1.0498E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2199.84641479379        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2354
 NPARAMETR:  1.0612E+00  1.4898E+00  7.8428E-01  6.5360E-01  1.1811E+00  1.1010E+00  8.0830E-01  1.0000E-02  1.3820E+00  1.0623E+00
             9.3216E-01
 PARAMETER:  1.5944E-01  4.9866E-01 -1.4299E-01 -3.2526E-01  2.6645E-01  1.9624E-01 -1.1282E-01 -5.3557E+00  4.2350E-01  1.6048E-01
             2.9748E-02
 GRADIENT:   6.9566E+00 -5.7125E+00  2.7357E-03 -7.4435E-01 -1.8316E-01  1.8421E+00  2.1321E-01  0.0000E+00  4.0982E-01  1.5176E-01
             5.3345E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2199.84735956878        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2548             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0612E+00  1.4876E+00  7.8488E-01  6.5510E-01  1.1804E+00  1.1010E+00  8.0638E-01  1.0000E-02  1.3780E+00  1.0607E+00
             9.3208E-01
 PARAMETER:  1.5943E-01  4.9714E-01 -1.4223E-01 -3.2297E-01  2.6583E-01  1.9623E-01 -1.1520E-01 -5.3557E+00  4.2061E-01  1.5894E-01
             2.9664E-02
 GRADIENT:   1.0080E+03  6.6061E+02  1.2424E+00  1.7539E+02  1.9475E+01  1.9453E+02  8.3344E+00  0.0000E+00  3.5343E+01  1.2248E+00
             9.4059E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2199.84818728303        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2737
 NPARAMETR:  1.0612E+00  1.4863E+00  7.8537E-01  6.5588E-01  1.1798E+00  1.1010E+00  8.0753E-01  1.0000E-02  1.3770E+00  1.0603E+00
             9.3209E-01
 PARAMETER:  1.5943E-01  4.9631E-01 -1.4160E-01 -3.2178E-01  2.6532E-01  1.9623E-01 -1.1377E-01 -5.3557E+00  4.1992E-01  1.5852E-01
             2.9671E-02
 GRADIENT:   6.9427E+00 -5.7578E+00 -1.0706E-01 -6.3974E-01  3.8955E-01  1.8429E+00 -1.8431E-01  0.0000E+00  9.3330E-02 -1.3691E-01
            -4.2656E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2199.84871220062        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2927             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0612E+00  1.4858E+00  7.8568E-01  6.5626E-01  1.1793E+00  1.1010E+00  8.0928E-01  1.0000E-02  1.3768E+00  1.0607E+00
             9.3211E-01
 PARAMETER:  1.5943E-01  4.9592E-01 -1.4121E-01 -3.2120E-01  2.6494E-01  1.9622E-01 -1.1161E-01 -5.3557E+00  4.1973E-01  1.5889E-01
             2.9691E-02
 GRADIENT:   1.0079E+03  6.5735E+02  1.3526E+00  1.7505E+02  1.8947E+01  1.9448E+02  8.3833E+00  0.0000E+00  3.5552E+01  1.4132E+00
             1.0067E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2199.84940651638        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3114
 NPARAMETR:  1.0612E+00  1.4843E+00  7.8620E-01  6.5719E-01  1.1787E+00  1.1010E+00  8.1008E-01  1.0000E-02  1.3750E+00  1.0601E+00
             9.3209E-01
 PARAMETER:  1.5942E-01  4.9496E-01 -1.4054E-01 -3.1978E-01  2.6442E-01  1.9622E-01 -1.1062E-01 -5.3557E+00  4.1848E-01  1.5840E-01
             2.9678E-02
 GRADIENT:   6.9460E+00 -5.5696E+00 -2.4564E-02 -6.9328E-01  3.3675E-02  1.8407E+00  3.1846E-02  0.0000E+00  2.3430E-01  2.4139E-02
             1.3170E-03

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2199.84986714496        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     3304             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0612E+00  1.4824E+00  7.8667E-01  6.5841E-01  1.1781E+00  1.1010E+00  8.0910E-01  1.0000E-02  1.3724E+00  1.0591E+00
             9.3206E-01
 PARAMETER:  1.5941E-01  4.9368E-01 -1.3995E-01 -3.1792E-01  2.6388E-01  1.9621E-01 -1.1184E-01 -5.3557E+00  4.1659E-01  1.5742E-01
             2.9644E-02
 GRADIENT:   1.0077E+03  6.5137E+02  1.1930E+00  1.7443E+02  1.9427E+01  1.9442E+02  8.0851E+00  0.0000E+00  3.5119E+01  1.2232E+00
             9.6199E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2199.85058689777        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3493
 NPARAMETR:  1.0612E+00  1.4811E+00  7.8719E-01  6.5931E-01  1.1775E+00  1.1010E+00  8.1017E-01  1.0000E-02  1.3708E+00  1.0586E+00
             9.3205E-01
 PARAMETER:  1.5941E-01  4.9276E-01 -1.3929E-01 -3.1656E-01  2.6337E-01  1.9621E-01 -1.1051E-01 -5.3557E+00  4.1537E-01  1.5696E-01
             2.9629E-02
 GRADIENT:   6.9369E+00 -5.6603E+00 -1.6532E-01 -5.8693E-01  5.2154E-01  1.8401E+00 -2.1624E-01  0.0000E+00  4.4188E-02 -1.3953E-01
            -4.4041E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -2199.85140190395        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3685             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0612E+00  1.4795E+00  7.8814E-01  6.6037E-01  1.1763E+00  1.1010E+00  8.1384E-01  1.0000E-02  1.3700E+00  1.0589E+00
             9.3207E-01
 PARAMETER:  1.5940E-01  4.9168E-01 -1.3807E-01 -3.1495E-01  2.6235E-01  1.9620E-01 -1.0600E-01 -5.3557E+00  4.1484E-01  1.5723E-01
             2.9653E-02
 GRADIENT:   1.0075E+03  6.4628E+02  1.4851E+00  1.7384E+02  1.8385E+01  1.9433E+02  8.1449E+00  0.0000E+00  3.5362E+01  1.5017E+00
             1.0218E+00

0ITERATION NO.:  110    OBJECTIVE VALUE:  -2199.85233971489        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3872
 NPARAMETR:  1.0612E+00  1.4776E+00  7.8890E-01  6.6158E-01  1.1754E+00  1.1010E+00  8.1413E-01  1.0000E-02  1.3678E+00  1.0584E+00
             9.3209E-01
 PARAMETER:  1.5940E-01  4.9042E-01 -1.3711E-01 -3.1312E-01  2.6162E-01  1.9619E-01 -1.0564E-01 -5.3557E+00  4.1321E-01  1.5678E-01
             2.9669E-02
 GRADIENT:   6.9452E+00 -5.1548E+00  1.4825E-01 -7.6266E-01 -4.8230E-01  1.8379E+00  1.2937E-01  0.0000E+00  2.9801E-01  1.2643E-01
             3.6019E-02

0ITERATION NO.:  115    OBJECTIVE VALUE:  -2199.85275560963        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     4062             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0612E+00  1.4769E+00  7.8886E-01  6.6199E-01  1.1754E+00  1.1010E+00  8.1366E-01  1.0000E-02  1.3668E+00  1.0578E+00
             9.3205E-01
 PARAMETER:  1.5939E-01  4.8998E-01 -1.3716E-01 -3.1251E-01  2.6165E-01  1.9619E-01 -1.0622E-01 -5.3557E+00  4.1250E-01  1.5622E-01
             2.9634E-02
 GRADIENT:   1.0073E+03  6.4172E+02  1.3036E+00  1.7337E+02  1.8890E+01  1.9428E+02  7.9191E+00  0.0000E+00  3.5039E+01  1.3654E+00
             1.0059E+00

0ITERATION NO.:  119    OBJECTIVE VALUE:  -2199.85308256037        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     4204
 NPARAMETR:  1.0612E+00  1.4752E+00  7.8870E-01  6.6304E-01  1.1750E+00  1.1010E+00  8.1435E-01  1.0000E-02  1.3654E+00  1.0573E+00
             9.3208E-01
 PARAMETER:  1.5939E-01  4.8923E-01 -1.3601E-01 -3.1131E-01  2.6083E-01  1.9619E-01 -1.0515E-01 -5.3557E+00  4.1135E-01  1.5588E-01
             2.9588E-02
 GRADIENT:   7.6703E-03  7.7999E-01  2.1346E-01 -3.1034E-01 -6.2111E-01  2.1812E-03  2.1388E-02  0.0000E+00 -1.9187E-02  2.4865E-02
            -5.9442E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4204
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3381E-04 -3.1389E-02 -3.1713E-04  2.2036E-02 -3.2250E-02
 SE:             2.9913E-02  2.0848E-02  1.2019E-04  2.4475E-02  2.3813E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9376E-01  1.3216E-01  8.3273E-03  3.6792E-01  1.7563E-01

 ETASHRINKSD(%)  1.0000E-10  3.0157E+01  9.9597E+01  1.8007E+01  2.0224E+01
 ETASHRINKVR(%)  1.0000E-10  5.1220E+01  9.9998E+01  3.2772E+01  3.6358E+01
 EBVSHRINKSD(%)  2.3837E-01  2.9252E+01  9.9637E+01  1.9178E+01  1.8140E+01
 EBVSHRINKVR(%)  4.7616E-01  4.9947E+01  9.9999E+01  3.4679E+01  3.2990E+01
 RELATIVEINF(%)  9.9379E+01  3.3948E+00  2.3840E-04  5.0357E+00  1.5733E+01
 EPSSHRINKSD(%)  3.3455E+01
 EPSSHRINKVR(%)  5.5718E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2199.8530825603693     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1280.9145493556966     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    75.02
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     7.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2199.853       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.48E+00  7.90E-01  6.63E-01  1.17E+00  1.10E+00  8.15E-01  1.00E-02  1.37E+00  1.06E+00  9.32E-01
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        8.40E+02
 
 TH 2
+       -7.21E+01  4.47E+02
 
 TH 3
+       -2.39E+01  1.85E+02  2.70E+02
 
 TH 4
+       -8.28E+01  4.12E+02 -1.45E+02  9.26E+02
 
 TH 5
+       -1.27E+01 -1.63E+02 -2.56E+02  2.37E+02  4.63E+02
 
 TH 6
+        7.49E+00 -1.44E+01 -2.11E+01  1.24E+01  2.20E+01  1.65E+02
 
 TH 7
+        1.43E+01 -2.28E+01 -1.47E+01 -3.01E+01 -3.67E+01 -6.67E+00  1.31E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.88E+00 -2.57E+01 -4.85E+01  2.62E+01  1.19E+01 -7.68E-01  1.05E+01  0.00E+00  1.46E+01
 
 TH10
+        3.20E+01 -6.99E+01 -6.48E+01 -4.48E+01 -5.05E+01  3.89E+00  3.00E+01  0.00E+00  2.86E+01  7.41E+01
 
 TH11
+        9.00E+01 -5.90E+00 -8.84E+00 -1.82E+01 -3.69E+01  2.93E+01  6.32E+00  0.00E+00  1.10E-01  2.26E+01  4.92E+02
 
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
+        8.10E+02
 
 TH 2
+       -4.33E+00  4.11E+02
 
 TH 3
+        4.30E+00  1.33E+02  2.49E+02
 
 TH 4
+       -3.01E+00  3.77E+02 -1.64E+02  8.92E+02
 
 TH 5
+        8.69E-01 -1.86E+02 -2.42E+02  1.76E+02  4.80E+02
 
 TH 6
+        2.32E-01 -7.84E-01  1.35E+00 -1.03E+00  6.17E-02  1.63E+02
 
 TH 7
+        5.89E-01  4.66E+00  9.55E+00 -1.76E+01 -1.79E+01 -1.67E-01  7.39E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.02E-01 -2.39E+01 -2.25E+01  4.88E+01  3.53E+00 -3.09E-01  2.54E+01  0.00E+00  5.28E+01
 
 TH10
+        8.61E-01 -1.04E+01 -3.58E+01 -2.51E+00 -4.65E+01  2.91E-01  2.05E+01  0.00E+00  2.26E+00  7.88E+01
 
 TH11
+       -5.67E+00 -2.29E+01 -3.33E+01  1.28E+00  2.96E+00  1.37E+00  5.12E+00  0.00E+00  3.69E+00  1.93E+01  4.75E+02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        8.10E+02
 
 TH 2
+        5.73E+01  3.94E+02
 
 TH 3
+        2.28E+01  1.24E+02  2.09E+02
 
 TH 4
+        5.36E+01  3.72E+02 -1.53E+02  8.94E+02
 
 TH 5
+        1.26E+01 -2.28E+02 -2.38E+02  1.31E+02  4.92E+02
 
 TH 6
+        1.56E+00  1.57E+01  1.95E+01 -9.93E+00 -2.88E+01  1.63E+02
 
 TH 7
+        2.00E+01  4.07E+01  1.61E+01  3.47E+01 -3.68E+01 -8.10E+00  5.86E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.41E+01  8.39E+00 -2.19E+01  9.36E+01 -1.55E+00 -3.36E+00  2.39E+01  0.00E+00  5.83E+01
 
 TH10
+       -1.83E+01  2.35E+01 -2.28E+01  2.22E+01 -4.50E+01  6.51E-01  2.16E+01  0.00E+00 -7.01E+00  7.74E+01
 
 TH11
+       -9.88E+01 -4.74E+01 -6.05E+01  1.23E+01  4.77E+01 -2.80E+01  3.66E+00  0.00E+00  8.17E+00  2.00E+01  4.77E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       82.633
Stop Time:
Wed Sep 29 21:13:49 CDT 2021
