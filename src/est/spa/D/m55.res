Wed Sep 29 20:08:23 CDT 2021
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
$DATA ../../../../data/spa/D/dat55.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m55.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17964.6279329186        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.9894E+02  3.2401E+02 -1.8817E+01  3.2827E+02  2.7913E+02 -1.9831E+03 -7.2709E+02 -5.8135E+01 -1.0923E+03 -6.4433E+02
            -3.4304E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -561.121659942230        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4185E+00  1.3095E+00  9.2795E-01  1.7391E+00  1.1740E+00  1.7703E+00  1.0821E+00  9.7770E-01  9.5114E-01  1.0799E+00
             1.5060E+01
 PARAMETER:  4.4957E-01  3.6964E-01  2.5222E-02  6.5337E-01  2.6042E-01  6.7117E-01  1.7893E-01  7.7447E-02  4.9906E-02  1.7687E-01
             2.8121E+00
 GRADIENT:   2.7759E+01  4.6136E+01  2.2353E+00  8.2145E+01 -1.2728E+01  3.0353E+01 -1.4956E+00  3.9100E+00  3.1729E+00  3.1221E+00
             4.1307E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -584.971371094703        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2730E+00  9.9525E-01  1.2762E+00  1.5096E+00  1.8587E+00  1.4330E+00  1.4288E+00  2.7512E-01  6.7208E-01  2.4498E+00
             1.4549E+01
 PARAMETER:  3.4137E-01  9.5242E-02  3.4389E-01  5.1187E-01  7.1987E-01  4.5979E-01  4.5682E-01 -1.1905E+00 -2.9737E-01  9.9602E-01
             2.7775E+00
 GRADIENT:   1.2799E+01  1.4470E+01  3.7901E+00  1.3700E+01 -8.9014E+00  6.1748E+00  4.9912E+00  1.3934E-01  6.2397E+00  1.2893E+00
             8.7213E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -595.465600629076        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0946E+00  6.3620E-01  1.9903E+00  1.4302E+00  4.5602E+00  1.3359E+00  1.0761E+00  2.5534E-01  2.7965E-01  6.2524E+00
             1.3377E+01
 PARAMETER:  1.9034E-01 -3.5223E-01  7.8829E-01  4.5783E-01  1.6174E+00  3.8964E-01  1.7334E-01 -1.2651E+00 -1.1742E+00  1.9330E+00
             2.6935E+00
 GRADIENT:  -3.0199E+01  3.9321E-02  1.0816E+01 -1.9005E+01 -3.8748E+00  1.2878E+01  1.9388E+00  1.0061E-02  1.7241E+00 -1.0060E+00
             7.2907E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -618.962653900391        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0352E+00  2.6436E-01  5.4836E-01  1.3946E+00  6.6538E+00  1.3510E+00  4.3252E-01  1.0000E-02  2.8010E-02  8.0735E+00
             1.1808E+01
 PARAMETER:  1.3461E-01 -1.2304E+00 -5.0083E-01  4.3258E-01  1.9952E+00  4.0086E-01 -7.3813E-01 -5.5119E+00 -3.4752E+00  2.1886E+00
             2.5688E+00
 GRADIENT:   2.2237E+01  6.5075E+00 -2.6371E+00  6.2852E+00 -9.9152E+00 -1.1981E-01  6.8032E-02  0.0000E+00  4.2688E-02  1.4599E+01
            -2.0190E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -637.012914647746        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  8.7504E-01  6.6018E-02  2.0813E-01  1.1933E+00  1.4979E+01  1.3440E+00  4.1175E-02  1.0000E-02  1.0000E-02  2.1610E+00
             1.1341E+01
 PARAMETER: -3.3486E-02 -2.6178E+00 -1.4696E+00  2.7674E-01  2.8066E+00  3.9566E-01 -3.0899E+00 -9.5534E+00 -5.6361E+00  8.7057E-01
             2.5284E+00
 GRADIENT:   6.9914E+00  2.7806E+00 -5.0627E+01  2.3044E+02  5.5879E-01 -6.2062E+01  1.6657E-05  0.0000E+00  0.0000E+00 -1.5110E-02
            -1.1122E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -674.611609258807        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  6.4532E-01  5.7484E-02  7.2337E-02  7.1786E-01  2.2904E+01  1.7375E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.2887E-01
             1.2247E+01
 PARAMETER: -3.3801E-01 -2.7563E+00 -2.5264E+00 -2.3149E-01  3.2313E+00  6.5247E-01 -6.4283E+00 -1.1338E+01 -6.1787E+00 -7.4659E-01
             2.6053E+00
 GRADIENT:  -3.4617E+01 -3.9754E+00 -1.1222E+01  9.3485E+01  4.5795E-01  2.8330E+01  0.0000E+00  0.0000E+00  0.0000E+00  3.6093E-03
             5.3910E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -680.899962707815        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      564
 NPARAMETR:  6.3521E-01  5.3490E-02  6.5641E-02  6.4893E-01  2.5953E+01  1.5653E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.8815E-01
             1.1550E+01
 PARAMETER: -3.5380E-01 -2.8283E+00 -2.6236E+00 -3.3242E-01  3.3563E+00  5.4806E-01 -6.3422E+00 -1.1225E+01 -6.2474E+00 -6.1714E-01
             2.5467E+00
 GRADIENT:  -3.2454E+01  1.1765E+00 -2.1024E+01  6.1842E+01  8.0510E-03 -4.6856E+00  0.0000E+00  0.0000E+00  0.0000E+00  3.6851E-03
            -5.8847E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -687.204069508093        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      739
 NPARAMETR:  5.1516E-01  2.2225E-02  3.2170E-02  3.6829E-01  4.0671E+01  1.5380E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.4704E-01
             1.1469E+01
 PARAMETER: -5.6327E-01 -3.7065E+00 -3.3367E+00 -8.9889E-01  3.8055E+00  5.3050E-01 -9.3919E+00 -1.2629E+01 -7.7542E+00 -1.8170E+00
             2.5397E+00
 GRADIENT:   1.6798E+01  2.4897E+00  2.3508E+00 -1.4147E+01 -4.0352E-02  2.7341E+00  0.0000E+00  0.0000E+00  0.0000E+00  7.7367E-06
            -4.7757E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -687.916190210441        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      916
 NPARAMETR:  4.8133E-01  1.3416E-02  2.8615E-02  3.3808E-01  2.0164E+02  1.5221E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.1773E-01
             1.1482E+01
 PARAMETER: -6.3120E-01 -4.2113E+00 -3.4538E+00 -9.8448E-01  5.4065E+00  5.2010E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -2.0393E+00
             2.5408E+00
 GRADIENT:   8.6642E-01  1.1535E-01 -2.8564E+00  2.2340E+00  1.7774E-03 -7.3525E-01  0.0000E+00  0.0000E+00  0.0000E+00 -8.2116E-09
             7.6370E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -687.926211825567        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1104             RESET HESSIAN, TYPE I
 NPARAMETR:  4.8018E-01  1.1518E-02  2.8605E-02  3.3758E-01  2.4223E+02  1.5264E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.4712E-01
             1.1475E+01
 PARAMETER: -6.3359E-01 -4.3638E+00 -3.4542E+00 -9.8595E-01  5.5899E+00  5.2293E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -1.8165E+00
             2.5402E+00
 GRADIENT:   3.6701E+01  2.4602E-02  4.4713E+01  1.6680E+01  1.5756E-03  6.1971E+00  0.0000E+00  0.0000E+00  0.0000E+00 -2.8163E-09
             2.0416E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -687.927562542741        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1302             RESET HESSIAN, TYPE I
 NPARAMETR:  4.8022E-01  1.1394E-02  2.8572E-02  3.3754E-01  1.8398E+02  1.5261E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.4041E-01
             1.1472E+01
 PARAMETER: -6.3351E-01 -4.3747E+00 -3.4553E+00 -9.8608E-01  5.3148E+00  5.2272E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -1.8632E+00
             2.5399E+00
 GRADIENT:   3.7135E+01  1.8483E-02  4.3533E+01  1.8140E+01  2.3285E-03  6.0502E+00  0.0000E+00  0.0000E+00  0.0000E+00 -4.3445E-09
             1.9958E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -687.928833931956        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1501
 NPARAMETR:  4.8011E-01  1.1398E-02  2.8549E-02  3.3735E-01  1.4906E+02  1.5261E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6451E-01
             1.1472E+01
 PARAMETER: -6.3373E-01 -4.3743E+00 -3.4561E+00 -9.8663E-01  5.1044E+00  5.2269E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -1.7048E+00
             2.5399E+00
 GRADIENT:   6.5094E-01  8.3408E-05 -1.8515E+00  1.0721E+00  2.9281E-03 -3.4558E-03  0.0000E+00  0.0000E+00  0.0000E+00 -3.2470E-08
            -3.5080E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -687.930189692198        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1697
 NPARAMETR:  4.7985E-01  1.1377E-02  2.8545E-02  3.3709E-01  1.1596E+02  1.5261E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.6100E-01
             1.1474E+01
 PARAMETER: -6.3429E-01 -4.3762E+00 -3.4563E+00 -9.8740E-01  4.8532E+00  5.2269E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -1.7264E+00
             2.5401E+00
 GRADIENT:   1.9118E-01  2.0701E-03 -6.1762E-01 -3.5023E-01  3.3637E-03  5.7226E-02  0.0000E+00  0.0000E+00  0.0000E+00 -5.0980E-08
             2.9317E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -687.931820414763        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     1894
 NPARAMETR:  4.7994E-01  1.1311E-02  2.8512E-02  3.3705E-01  9.1324E+01  1.5261E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.8314E-01
             1.1471E+01
 PARAMETER: -6.3410E-01 -4.3820E+00 -3.4574E+00 -9.8753E-01  4.6144E+00  5.2273E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -1.5975E+00
             2.5398E+00
 GRADIENT:   6.4324E-01 -7.4405E-04 -1.8620E+00  1.0627E+00  4.7618E-03 -4.3683E-03  0.0000E+00  0.0000E+00  0.0000E+00 -1.0718E-07
            -3.9153E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -687.933691857671        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     2093             RESET HESSIAN, TYPE I
 NPARAMETR:  4.7990E-01  1.1268E-02  2.8497E-02  3.3694E-01  6.9651E+01  1.5262E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7700E-01
             1.1471E+01
 PARAMETER: -6.3417E-01 -4.3858E+00 -3.4580E+00 -9.8786E-01  4.3435E+00  5.2276E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -1.6316E+00
             2.5398E+00
 GRADIENT:   3.7169E+01  1.8132E-02  4.3580E+01  1.8220E+01  6.1906E-03  6.0420E+00  0.0000E+00  0.0000E+00  0.0000E+00 -6.3155E-08
             1.9928E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -687.934832838572        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2278
 NPARAMETR:  4.7974E-01  1.1253E-02  2.8508E-02  3.3681E-01  5.8669E+01  1.5262E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7761E-01
             1.1473E+01
 PARAMETER: -6.3450E-01 -4.3872E+00 -3.4576E+00 -9.8823E-01  4.1719E+00  5.2279E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -1.6282E+00
             2.5400E+00
 GRADIENT:   1.6330E-01  1.6611E-03 -6.2094E-01 -3.6393E-01  6.5739E-03  5.6718E-02  0.0000E+00  0.0000E+00  0.0000E+00 -2.4206E-07
             2.2716E-03

0ITERATION NO.:   85    OBJECTIVE VALUE:  -687.936307229908        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2464
 NPARAMETR:  4.7962E-01  1.1021E-02  2.8549E-02  3.3701E-01  3.6182E+01  1.5258E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.8038E-01
             1.1472E+01
 PARAMETER: -6.3476E-01 -4.4079E+00 -3.4561E+00 -9.8763E-01  3.6886E+00  5.2252E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -1.6127E+00
             2.5399E+00
 GRADIENT:  -7.3784E-01  2.7582E-03  6.8853E-01 -1.6347E+00  8.8505E-03 -5.1065E-02  0.0000E+00  0.0000E+00  0.0000E+00 -6.5988E-07
             2.9088E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -687.937992932465        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2639
 NPARAMETR:  4.7988E-01  1.0119E-02  2.8689E-02  3.3799E-01  1.1176E+01  1.5254E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.9259E-01
             1.1469E+01
 PARAMETER: -6.3421E-01 -4.4934E+00 -3.4512E+00 -9.8475E-01  2.5138E+00  5.2223E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -1.5472E+00
             2.5396E+00
 GRADIENT:  -3.6033E+00  2.3309E-02  4.5953E+00 -5.3868E+00 -3.3838E-03 -3.6355E-01  0.0000E+00  0.0000E+00  0.0000E+00 -2.5644E-07
             1.1829E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -687.953137250254        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2819
 NPARAMETR:  4.8387E-01  1.0000E-02  2.8862E-02  3.4088E-01  9.1260E+00  1.5293E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.5026E-02
             1.1470E+01
 PARAMETER: -6.2595E-01 -4.5543E+00 -3.4452E+00 -9.7622E-01  2.3111E+00  5.2480E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -4.0979E+00
             2.5397E+00
 GRADIENT:   2.9729E-01  0.0000E+00 -1.1412E+00  3.6785E-01 -2.1525E-03  3.1713E-02  0.0000E+00  0.0000E+00  0.0000E+00  2.2190E-06
            -3.2902E-01

0ITERATION NO.:   98    OBJECTIVE VALUE:  -687.954261475342        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     2924
 NPARAMETR:  4.8409E-01  1.0000E-02  2.8808E-02  3.4049E-01  9.2089E+00  1.5298E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.1479E+01
 PARAMETER: -6.2548E-01 -4.5543E+00 -3.4471E+00 -9.7738E-01  2.3202E+00  5.2517E-01 -9.8855E+00 -1.2868E+01 -7.9798E+00 -4.7659E+00
             2.5405E+00
 GRADIENT:   1.1463E+00  0.0000E+00 -1.7972E+00  7.7146E-01  1.1854E-02  1.7101E-01  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -5.3807E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2924
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8622E-03  2.0565E-06  1.0623E-04 -2.2526E-04  3.7227E-07
 SE:             2.8725E-02  9.7106E-07  2.4898E-04  3.0682E-04  3.0192E-06
 N:                     100         100         100         100         100

 P VAL.:         9.4831E-01  3.4192E-02  6.6964E-01  4.6284E-01  9.0187E-01

 ETASHRINKSD(%)  3.7667E+00  9.9997E+01  9.9166E+01  9.8972E+01  9.9990E+01
 ETASHRINKVR(%)  7.3915E+00  1.0000E+02  9.9993E+01  9.9989E+01  1.0000E+02
 EBVSHRINKSD(%)  3.9609E+00  9.9996E+01  9.9069E+01  9.8822E+01  9.9988E+01
 EBVSHRINKVR(%)  7.7649E+00  1.0000E+02  9.9991E+01  9.9986E+01  1.0000E+02
 RELATIVEINF(%)  4.7252E+00  1.3145E-08  7.3733E-05  1.1732E-04  2.5241E-07
 EPSSHRINKSD(%)  6.5412E+00
 EPSSHRINKVR(%)  1.2655E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -687.95426147534249     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       47.196565088395687     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.41
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -687.954       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.84E-01  1.00E-02  2.88E-02  3.40E-01  9.21E+00  1.53E+00  1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.15E+01
 


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
+        8.04E+01
 
 TH 2
+        1.13E-16 -9.12E-34
 
 TH 3
+       -7.11E+03  5.72E-14  6.28E+05
 
 TH 4
+        7.40E+02 -5.96E-15 -6.55E+04  6.82E+03
 
 TH 5
+        2.41E-01 -1.94E-18 -2.13E+01  2.22E+00  7.24E-04
 
 TH 6
+       -6.32E+00  5.09E-17  5.59E+02 -5.82E+01 -1.90E-02  4.97E-01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.70E+00  2.18E-17  2.39E+02 -2.49E+01 -8.12E-03  2.13E-01  0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.10E-02
 
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
+        1.84E+03
 
 TH 2
+        0.00E+00  3.27E+02
 
 TH 3
+       -8.57E+03  0.00E+00  7.50E+05
 
 TH 4
+       -3.22E+02  0.00E+00 -7.80E+04  9.06E+03
 
 TH 5
+        6.43E-01  0.00E+00 -2.55E+01  2.48E+00  3.64E-03
 
 TH 6
+       -1.25E+00  0.00E+00  6.65E+02 -8.79E+01  6.01E-03  7.12E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.11E+01  0.00E+00  2.86E+02 -1.96E+01 -1.50E-02  7.76E-01  0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.08E+00
 
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
+        1.87E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.29E+04  0.00E+00  8.96E+05
 
 TH 4
+       -7.23E+01  0.00E+00 -8.77E+04  9.73E+03
 
 TH 5
+        8.08E-01  0.00E+00 -3.10E+01  2.80E+00  1.29E-03
 
 TH 6
+        2.92E+02  0.00E+00 -3.18E+03  8.59E+01  1.86E-01  1.32E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.67E+01  0.00E+00  5.21E+02 -4.69E+01 -2.72E-02 -5.92E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  2.26E+01
 
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
 #CPUT: Total CPU Time in Seconds,       46.916
Stop Time:
Wed Sep 29 20:09:11 CDT 2021
