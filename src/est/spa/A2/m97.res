Wed Sep 29 13:12:10 CDT 2021
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
$DATA ../../../../data/spa/A2/dat97.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1112.96098494114        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1261E+02 -4.1227E+01  2.0273E+01 -6.0358E+01  2.8199E+01  1.4952E+01  3.9877E+00 -7.5124E+00  1.2874E+01 -4.3398E+00
            -9.5172E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1416.51698407902        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9033E-01  1.1848E+00  1.1119E+00  1.0343E+00  1.1352E+00  1.0906E+00  8.5532E-01  9.4743E-01  6.5956E-01  7.1981E-01
             2.5873E+00
 PARAMETER:  9.0286E-02  2.6957E-01  2.0604E-01  1.3373E-01  2.2681E-01  1.8675E-01 -5.6278E-02  4.5997E-02 -3.1619E-01 -2.2876E-01
             1.0506E+00
 GRADIENT:   1.4521E+02  5.6996E+01 -2.2851E+00  7.5199E+01 -4.9549E+00  4.1601E+01 -8.1943E+00  5.7806E-01 -5.8540E+00  6.7414E+00
             4.3321E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1422.85354774354        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.6395E-01  1.0970E+00  7.8241E-01  1.0427E+00  8.6471E-01  1.0624E+00  1.0625E+00  9.5916E-01  6.8581E-01  3.1457E-01
             2.5135E+00
 PARAMETER:  6.3289E-02  1.9260E-01 -1.4537E-01  1.4182E-01 -4.5360E-02  1.6055E-01  1.6060E-01  5.8304E-02 -2.7715E-01 -1.0565E+00
             1.0217E+00
 GRADIENT:   1.0157E+02  4.6537E+01  1.8736E+01  3.4175E+01 -4.4389E+01  3.4526E+01  3.2446E+00  4.1464E+00  5.0600E+00  2.1940E+00
             5.2618E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1428.26838805569        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  9.2009E-01  1.0514E+00  4.5081E-01  1.0022E+00  6.5581E-01  9.8728E-01  1.0551E+00  4.8282E-01  6.5314E-01  2.4089E-01
             2.3316E+00
 PARAMETER:  1.6721E-02  1.5009E-01 -6.9671E-01  1.0216E-01 -3.2188E-01  8.7196E-02  1.5364E-01 -6.2812E-01 -3.2596E-01 -1.3234E+00
             9.4657E-01
 GRADIENT:  -4.5459E+01  3.2809E+00 -1.2872E+01  2.9083E+01  1.6827E+01 -2.8288E+00  1.7512E+00  1.6780E+00  2.6529E+00  2.0009E+00
             1.4257E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1430.75003802742        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      471
 NPARAMETR:  9.3730E-01  1.0271E+00  2.9440E-01  9.4028E-01  5.3023E-01  9.9982E-01  9.8767E-01  2.1771E-01  6.6137E-01  1.0844E-01
             2.2233E+00
 PARAMETER:  3.5252E-02  1.2672E-01 -1.1228E+00  3.8419E-02 -5.3445E-01  9.9825E-02  8.7591E-02 -1.4246E+00 -3.1344E-01 -2.1216E+00
             8.9898E-01
 GRADIENT:   8.3255E+00 -3.0762E+00 -6.7322E-01 -7.6716E-01  2.1046E-01  1.4120E+00  3.6181E-02  9.8174E-02  7.7269E-01  6.0827E-01
            -4.3068E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1431.80245010360        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      649
 NPARAMETR:  9.3466E-01  1.3308E+00  2.3463E-01  7.6674E-01  6.2073E-01  9.8019E-01  8.1547E-01  2.4071E-01  7.3830E-01  1.0000E-02
             2.2268E+00
 PARAMETER:  3.2428E-02  3.8578E-01 -1.3497E+00 -1.6561E-01 -3.7687E-01  7.9992E-02 -1.0400E-01 -1.3242E+00 -2.0340E-01 -6.1565E+00
             9.0055E-01
 GRADIENT:   1.1970E+01  2.5067E+01  8.3736E+00  6.1255E+00 -2.1085E+01 -4.1217E+00 -2.2027E+00 -1.7113E-01 -1.8277E+00  0.0000E+00
            -5.0646E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1443.83781339719        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      830
 NPARAMETR:  9.3652E-01  1.9683E+00  7.6049E-02  3.7764E-01  8.3396E-01  9.9080E-01  6.5193E-01  5.1871E-02  1.3528E+00  1.0000E-02
             2.0087E+00
 PARAMETER:  3.4420E-02  7.7715E-01 -2.4764E+00 -8.7383E-01 -8.1576E-02  9.0759E-02 -3.2781E-01 -2.8590E+00  4.0219E-01 -2.4379E+01
             7.9750E-01
 GRADIENT:   4.4502E+01  1.2709E+02  1.9968E+01  4.7677E+00 -1.0449E+02 -3.2568E+00 -9.4964E+00 -6.4925E-02 -5.8154E+00  0.0000E+00
            -2.7012E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1448.00689254165        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     1002
 NPARAMETR:  8.9609E-01  2.0118E+00  6.7598E-02  3.5193E-01  9.2676E-01  9.9267E-01  6.7246E-01  4.5460E-02  1.4584E+00  1.0000E-02
             2.0069E+00
 PARAMETER: -9.7149E-03  7.9902E-01 -2.5942E+00 -9.4432E-01  2.3935E-02  9.2638E-02 -2.9681E-01 -2.9909E+00  4.7735E-01 -2.6525E+01
             7.9658E-01
 GRADIENT:   2.3127E+01  1.7930E+02  1.1889E+01  6.2750E+01  6.5644E+01  4.9946E+00  1.0748E+01 -4.4924E-02 -8.6598E+00  0.0000E+00
            -1.5042E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1449.85146949496        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1160
 NPARAMETR:  9.1900E-01  2.0115E+00  6.7582E-02  3.5190E-01  8.9820E-01  9.9951E-01  6.5436E-01  5.7968E-02  1.4583E+00  1.0000E-02
             2.0074E+00
 PARAMETER:  1.5532E-02  7.9889E-01 -2.5944E+00 -9.4442E-01 -7.3646E-03  9.9513E-02 -3.2410E-01 -2.7479E+00  4.7729E-01 -2.6525E+01
             7.9682E-01
 GRADIENT:   6.5386E-01  3.9843E+01  1.0907E+01  3.0340E+01 -1.2919E+00  1.5006E-01 -4.6058E-01 -7.5364E-02 -7.8494E+00  0.0000E+00
            -2.0045E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1470.94708658100        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1341
 NPARAMETR:  9.2476E-01  2.0108E+00  6.7512E-02  3.5174E-01  9.0739E-01  9.8736E-01  6.4570E-01  2.2098E+00  1.4584E+00  1.0000E-02
             2.0088E+00
 PARAMETER:  2.1781E-02  7.9852E-01 -2.5954E+00 -9.4486E-01  2.8130E-03  8.7283E-02 -3.3742E-01  8.9290E-01  4.7736E-01 -2.6525E+01
             7.9753E-01
 GRADIENT:   1.6105E+01  1.8722E+01 -2.2443E+01  2.9673E+01 -2.1194E+01  7.3925E-01 -3.9447E-01  3.2272E-01 -1.6745E+00  0.0000E+00
             4.2547E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1474.31503509174        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1523
 NPARAMETR:  9.0938E-01  2.0112E+00  7.5061E-02  3.4516E-01  9.2293E-01  1.0798E+00  6.5374E-01  2.3777E+00  1.4621E+00  1.0000E-02
             1.8645E+00
 PARAMETER:  5.0061E-03  7.9873E-01 -2.4895E+00 -9.6375E-01  1.9795E-02  1.7676E-01 -3.2505E-01  9.6613E-01  4.7988E-01 -2.6525E+01
             7.2302E-01
 GRADIENT:  -1.7290E+01  8.0904E+00 -1.5458E+01  1.2446E+01 -8.7303E+00  3.1198E+01 -8.2609E-01  6.2071E-01 -5.6785E+00  0.0000E+00
             1.7293E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1478.47271041560        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1699
 NPARAMETR:  9.2069E-01  2.0142E+00  9.0512E-02  3.4404E-01  9.4048E-01  9.8827E-01  6.5795E-01  2.6535E+00  1.5596E+00  1.0000E-02
             1.7363E+00
 PARAMETER:  1.7370E-02  8.0024E-01 -2.3023E+00 -9.6698E-01  3.8631E-02  8.8200E-02 -3.1862E-01  1.0759E+00  5.4443E-01 -2.6525E+01
             6.5175E-01
 GRADIENT:   5.4103E+00  1.3797E+01 -3.8237E+00 -1.3503E+00 -3.5520E+00  3.5566E-01 -7.5183E-01 -6.9169E-01 -2.5373E+00  0.0000E+00
            -2.7836E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1478.91512023604        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1880
 NPARAMETR:  9.1886E-01  1.9941E+00  1.0400E-01  3.5497E-01  9.4224E-01  9.8639E-01  6.6074E-01  2.9079E+00  1.5489E+00  1.0000E-02
             1.7159E+00
 PARAMETER:  1.5377E-02  7.9018E-01 -2.1634E+00 -9.3573E-01  4.0507E-02  8.6292E-02 -3.1439E-01  1.1674E+00  5.3753E-01 -2.6525E+01
             6.3992E-01
 GRADIENT:  -1.5877E-01  1.7167E+00  6.4668E-02 -6.3318E+00  3.5208E-01  1.5144E-03 -1.9074E-01  7.5742E-01  2.8608E-01  0.0000E+00
            -6.5159E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1479.06831327447        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2058
 NPARAMETR:  9.1896E-01  1.9851E+00  1.0678E-01  3.6352E-01  9.3480E-01  9.8679E-01  6.6251E-01  2.8807E+00  1.5127E+00  1.0000E-02
             1.7164E+00
 PARAMETER:  1.5489E-02  7.8566E-01 -2.1370E+00 -9.1192E-01  3.2577E-02  8.6697E-02 -3.1173E-01  1.1580E+00  5.1389E-01 -2.6525E+01
             6.4026E-01
 GRADIENT:  -5.4614E-02  9.6845E+00  6.2666E-01 -4.6000E+00 -1.5399E+00  1.4931E-01 -3.1253E-01 -2.9422E-01  2.0441E-01  0.0000E+00
             4.7874E-01

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1479.06831327447        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     2087
 NPARAMETR:  9.1914E-01  1.9774E+00  1.0788E-01  3.6200E-01  9.3574E-01  9.8645E-01  6.6339E-01  2.8638E+00  1.5096E+00  1.0000E-02
             1.7220E+00
 PARAMETER:  1.5489E-02  7.8566E-01 -2.1370E+00 -9.1192E-01  3.2577E-02  8.6697E-02 -3.1173E-01  1.1580E+00  5.1389E-01 -2.6525E+01
             6.4026E-01
 GRADIENT:  -4.0451E-01  6.8439E+01 -1.9408E+01  4.9236E+01 -1.0071E+00  1.2110E-01 -3.2722E-01  3.7324E+01  1.9486E-01  0.0000E+00
            -6.8553E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2087
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.9724E-04 -2.1756E-02 -1.7008E-03  2.9630E-02 -5.9574E-04
 SE:             2.9541E-02  2.6728E-02  1.6415E-02  2.0981E-02  2.9577E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7307E-01  4.1566E-01  9.1748E-01  1.5788E-01  4.3988E-02

 ETASHRINKSD(%)  1.0340E+00  1.0458E+01  4.5007E+01  2.9712E+01  9.9009E+01
 ETASHRINKVR(%)  2.0573E+00  1.9822E+01  6.9757E+01  5.0596E+01  9.9990E+01
 EBVSHRINKSD(%)  1.2432E+00  1.1093E+01  4.3140E+01  3.0507E+01  9.9058E+01
 EBVSHRINKVR(%)  2.4709E+00  2.0955E+01  6.7669E+01  5.1708E+01  9.9991E+01
 RELATIVEINF(%)  9.6527E+01  1.9231E+01  1.9391E+01  1.3015E+01  2.0253E-03
 EPSSHRINKSD(%)  3.9337E+01
 EPSSHRINKVR(%)  6.3200E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1479.0683132744684     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -743.91748671073026     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.97
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.64
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1479.068       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.19E-01  1.99E+00  1.07E-01  3.64E-01  9.35E-01  9.87E-01  6.63E-01  2.88E+00  1.51E+00  1.00E-02  1.72E+00
 


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
+        1.32E+03
 
 TH 2
+       -9.66E+01  9.82E+02
 
 TH 3
+        1.43E+02  2.02E+02  2.42E+04
 
 TH 4
+       -7.38E+01  2.34E+02  6.22E+01  1.34E+04
 
 TH 5
+       -4.19E+01  7.32E+03 -2.15E+03  1.03E+03  1.27E+05
 
 TH 6
+        1.58E+00 -4.90E+01  1.32E+02 -3.68E+01 -6.72E+00  1.94E+02
 
 TH 7
+        2.89E+00 -4.07E+01  4.57E+01  9.18E+00  3.40E+01 -2.05E+00  2.92E+02
 
 TH 8
+       -5.46E+00  2.05E+00  4.84E+01 -1.02E+02  2.62E+01 -3.04E+00  1.02E+00  9.88E+01
 
 TH 9
+        2.17E+00 -9.39E+02  1.37E+02  3.75E+01 -5.40E+00 -1.49E-01  9.12E+00 -3.17E+00  2.18E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        5.65E+01 -6.75E+02  4.38E+03  1.23E+02 -4.62E+02  4.34E+01  4.49E+01  8.23E+00  4.95E+01  0.00E+00  9.95E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.663
Stop Time:
Wed Sep 29 13:12:53 CDT 2021
