Wed Sep 29 19:00:46 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat42.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1703.09497373419        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9238E+02  5.6164E+01 -2.1939E+01  1.3364E+02  6.9729E+01  5.1782E+01  7.2740E+00  7.5822E-01  6.8621E+00 -6.6499E-01
            -3.0970E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1708.45099987685        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0543E+00  8.8971E-01  9.0240E-01  1.0284E+00  8.6601E-01  9.9281E-01  8.9506E-01  9.9883E-01  1.0231E+00  8.8410E-01
             1.0246E+00
 PARAMETER:  1.5284E-01 -1.6865E-02 -2.7018E-03  1.2802E-01 -4.3855E-02  9.2784E-02 -1.0867E-02  9.8831E-02  1.2282E-01 -2.3188E-02
             1.2428E-01
 GRADIENT:   5.5092E+00 -1.2721E+01 -1.2525E+01 -7.2709E+00  6.6216E+00 -2.8550E+00  1.0536E+00  5.5895E+00  7.6810E+00  1.7806E+00
             8.8512E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1708.79713848419        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0542E+00  8.1487E-01  9.3337E-01  1.0884E+00  8.4728E-01  9.9166E-01  8.2991E-01  9.0235E-01  9.7416E-01  8.9303E-01
             1.0296E+00
 PARAMETER:  1.5280E-01 -1.0473E-01  3.1047E-02  1.8468E-01 -6.5722E-02  9.1629E-02 -8.6433E-02 -2.7483E-03  7.3819E-02 -1.3130E-02
             1.2919E-01
 GRADIENT:   6.4658E+00  5.0431E+00 -5.9332E+00  1.3298E+01 -9.0097E-01 -3.0049E+00 -1.5773E+00  1.5487E+00 -7.9531E-01 -1.0794E+00
             8.7900E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1709.18835025586        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0485E+00  7.2614E-01  1.0401E+00  1.1439E+00  8.6420E-01  9.9847E-01  9.0716E-01  9.4413E-01  9.3622E-01  9.2806E-01
             1.0046E+00
 PARAMETER:  1.4734E-01 -2.2001E-01  1.3934E-01  2.3447E-01 -4.5955E-02  9.8466E-02  2.5664E-03  4.2514E-02  3.4095E-02  2.5340E-02
             1.0454E-01
 GRADIENT:  -1.4405E+00  3.1870E+00  1.5683E+00  3.6226E+00 -1.5423E+00  3.2869E-01  5.2891E-03 -4.9987E-01 -6.9388E-01 -2.4402E-01
            -1.5820E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1709.20386985369        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  1.0477E+00  6.4666E-01  1.0817E+00  1.1944E+00  8.5374E-01  9.9611E-01  9.3008E-01  9.6673E-01  9.1144E-01  9.3719E-01
             1.0075E+00
 PARAMETER:  1.4662E-01 -3.3594E-01  1.7854E-01  2.7768E-01 -5.8132E-02  9.6101E-02  2.7518E-02  6.6166E-02  7.2731E-03  3.5133E-02
             1.0751E-01
 GRADIENT:  -6.5209E-01  3.3217E+00  2.4276E+00  3.9375E+00 -4.5577E+00 -1.2852E-01  4.3555E-01 -2.7561E-01  3.8851E-01  1.2174E+00
            -1.6680E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1709.20706554749        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  1.0472E+00  5.9603E-01  1.1086E+00  1.2261E+00  8.4826E-01  9.9482E-01  9.2816E-01  9.8807E-01  8.9696E-01  9.4032E-01
             1.0092E+00
 PARAMETER:  1.4611E-01 -4.1747E-01  2.0307E-01  3.0386E-01 -6.4563E-02  9.4805E-02  2.5451E-02  8.7994E-02 -8.7430E-03  3.8470E-02
             1.0911E-01
 GRADIENT:  -1.0603E-01  2.9204E+00  2.3974E+00  3.5677E+00 -5.2033E+00 -3.2836E-01  5.3688E-01 -9.1967E-02  8.7705E-01  1.6656E+00
             5.6934E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1709.32386278133        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1083
 NPARAMETR:  1.0489E+00  5.8022E-01  1.1168E+00  1.2294E+00  8.5034E-01  9.9543E-01  5.7742E-01  1.0130E+00  9.1320E-01  9.2999E-01
             1.0062E+00
 PARAMETER:  1.4773E-01 -4.4436E-01  2.1049E-01  3.0649E-01 -6.2115E-02  9.5424E-02 -4.4918E-01  1.1292E-01  9.1944E-03  2.7422E-02
             1.0620E-01
 GRADIENT:   4.1490E+00 -8.0911E-01 -1.8167E+00 -4.4016E+00  2.5818E+00 -7.3578E-02 -6.0683E-03 -6.2082E-01 -7.1499E-01 -1.7761E+00
            -1.0630E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1709.35614468926        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1259
 NPARAMETR:  1.0448E+00  5.5889E-01  1.1540E+00  1.2453E+00  8.5765E-01  9.9556E-01  3.4814E-01  1.0609E+00  9.2044E-01  9.4839E-01
             1.0090E+00
 PARAMETER:  1.4380E-01 -4.8180E-01  2.4327E-01  3.1941E-01 -5.3561E-02  9.5552E-02 -9.5516E-01  1.5912E-01  1.7101E-02  4.7012E-02
             1.0892E-01
 GRADIENT:  -3.7678E+00  5.5995E-01 -9.6015E-01  6.8361E-02 -1.0170E+00  9.8703E-02  2.2034E-02  1.4038E-01  1.3010E+00  5.9756E-01
             4.5562E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1709.36027046022        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1434
 NPARAMETR:  1.0442E+00  5.3925E-01  1.1867E+00  1.2589E+00  8.6404E-01  9.9525E-01  2.4057E-01  1.0953E+00  9.1438E-01  9.5109E-01
             1.0088E+00
 PARAMETER:  1.4330E-01 -5.1758E-01  2.7117E-01  3.3027E-01 -4.6138E-02  9.5235E-02 -1.3248E+00  1.9102E-01  1.0494E-02  4.9851E-02
             1.0876E-01
 GRADIENT:  -3.9391E+00  6.2804E-01  7.9731E-02  9.4210E-02 -7.6522E-01  1.2229E-01  1.8157E-02  5.6551E-02  1.2874E+00  5.0338E-01
             3.0358E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1709.37843032041        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     1577
 NPARAMETR:  1.0483E+00  5.3841E-01  1.1868E+00  1.2574E+00  8.6370E-01  9.9487E-01  2.0195E-01  1.0951E+00  9.0952E-01  9.4721E-01
             1.0082E+00
 PARAMETER:  1.4716E-01 -5.1913E-01  2.7128E-01  3.2905E-01 -4.6535E-02  9.4857E-02 -1.4998E+00  1.9084E-01  5.1628E-03  4.5765E-02
             1.0820E-01
 GRADIENT:   4.9264E+00  1.1974E-01  4.0691E-01 -3.8302E+00 -2.9526E-01  2.0726E-03  3.8225E-03 -1.5764E-01 -8.5472E-01 -1.4904E-01
            -1.5343E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1709.38036005445        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1756             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0480E+00  5.3805E-01  1.1868E+00  1.2578E+00  8.6373E-01  9.9515E-01  2.4178E-01  1.0964E+00  9.1024E-01  9.4771E-01
             1.0083E+00
 PARAMETER:  1.4684E-01 -5.1981E-01  2.7126E-01  3.2936E-01 -4.6491E-02  9.5139E-02 -1.3197E+00  1.9204E-01  5.9503E-03  4.6297E-02
             1.0828E-01
 GRADIENT:   7.2430E+02  8.0594E+01  6.2794E+00  4.4754E+02  8.2279E+00  4.6365E+01  1.0740E+00  5.2423E-01  1.2298E+01  8.8573E-01
             8.4402E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1709.38069333268        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1941             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0480E+00  5.3845E-01  1.1863E+00  1.2576E+00  8.6373E-01  9.9513E-01  2.7123E-01  1.0963E+00  9.0995E-01  9.4714E-01
             1.0083E+00
 PARAMETER:  1.4687E-01 -5.1906E-01  2.7087E-01  3.2917E-01 -4.6499E-02  9.5123E-02 -1.2048E+00  1.9191E-01  5.6396E-03  4.5691E-02
             1.0827E-01
 GRADIENT:   7.2449E+02  8.0331E+01  6.1984E+00  4.4726E+02  8.4252E+00  4.6354E+01  1.2316E+00  5.4055E-01  1.2428E+01  8.4003E-01
             8.5453E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1709.38078702446        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     2101
 NPARAMETR:  1.0479E+00  5.3886E-01  1.1861E+00  1.2575E+00  8.6362E-01  9.9515E-01  2.7877E-01  1.0964E+00  9.1006E-01  9.4739E-01
             1.0084E+00
 PARAMETER:  1.4683E-01 -5.1831E-01  2.7066E-01  3.2915E-01 -4.6620E-02  9.5138E-02 -1.1774E+00  1.9204E-01  5.7548E-03  4.5952E-02
             1.0833E-01
 GRADIENT:   4.1712E+00 -1.5835E-01  1.1190E-01 -3.5288E+00 -1.4900E-01  1.1627E-01  1.5756E-02  3.4475E-02  1.0021E-01  6.4016E-02
             3.1408E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1709.38089944771        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2284
 NPARAMETR:  1.0479E+00  5.3918E-01  1.1857E+00  1.2574E+00  8.6361E-01  9.9516E-01  2.7947E-01  1.0959E+00  9.1012E-01  9.4725E-01
             1.0084E+00
 PARAMETER:  1.4682E-01 -5.1770E-01  2.7037E-01  3.2905E-01 -4.6629E-02  9.5151E-02 -1.1749E+00  1.9159E-01  5.8203E-03  4.5805E-02
             1.0832E-01
 GRADIENT:   4.1287E+00 -1.1047E-01  9.5547E-02 -3.3600E+00 -1.0409E-01  1.1908E-01  1.5153E-02  1.8219E-02  7.0601E-02  3.2908E-02
             1.6631E-02

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1709.38089944771        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     2308
 NPARAMETR:  1.0479E+00  5.3918E-01  1.1857E+00  1.2574E+00  8.6361E-01  9.9516E-01  2.7947E-01  1.0959E+00  9.1012E-01  9.4725E-01
             1.0084E+00
 PARAMETER:  1.4682E-01 -5.1770E-01  2.7037E-01  3.2905E-01 -4.6629E-02  9.5151E-02 -1.1749E+00  1.9159E-01  5.8203E-03  4.5805E-02
             1.0832E-01
 GRADIENT:  -8.9982E-02 -1.0287E-02  1.0886E-01  4.8742E-01 -9.3125E-02  1.5630E-03  2.3843E-04  1.7012E-02  4.6577E-02  3.1752E-02
             1.7143E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2308
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8188E-04 -7.3716E-03 -3.0732E-02 -3.3279E-03 -3.0639E-02
 SE:             2.9868E-02  3.0916E-03  1.6736E-02  2.9108E-02  2.2189E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8446E-01  1.7107E-02  6.6319E-02  9.0898E-01  1.6733E-01

 ETASHRINKSD(%)  1.0000E-10  8.9643E+01  4.3932E+01  2.4831E+00  2.5664E+01
 ETASHRINKVR(%)  1.0000E-10  9.8927E+01  6.8564E+01  4.9045E+00  4.4742E+01
 EBVSHRINKSD(%)  4.2734E-01  9.0617E+01  4.7398E+01  2.7318E+00  2.3252E+01
 EBVSHRINKVR(%)  8.5285E-01  9.9120E+01  7.2331E+01  5.3889E+00  4.1097E+01
 RELATIVEINF(%)  9.8116E+01  5.7517E-02  6.0663E+00  7.9800E+00  9.1855E+00
 EPSSHRINKSD(%)  4.5370E+01
 EPSSHRINKVR(%)  7.0156E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1709.3808994477081     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -974.23007288396991     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.59
 Elapsed covariance  time in seconds:     5.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1709.381       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  5.39E-01  1.19E+00  1.26E+00  8.64E-01  9.95E-01  2.79E-01  1.10E+00  9.10E-01  9.47E-01  1.01E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.07E-02  3.90E-01  2.24E-01  2.39E-01  1.31E-01  6.28E-02  2.13E+00  3.47E-01  1.57E-01  1.49E-01  6.64E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.43E-04
 
 TH 2
+       -7.45E-05  1.52E-01
 
 TH 3
+        7.13E-04 -1.23E-02  5.02E-02
 
 TH 4
+        1.79E-05 -9.20E-02  1.08E-02  5.73E-02
 
 TH 5
+        2.66E-04  4.00E-02  1.33E-02 -2.30E-02  1.71E-02
 
 TH 6
+        4.63E-06  1.95E-03 -1.28E-03 -9.90E-04  1.54E-04  3.95E-03
 
 TH 7
+       -4.19E-03  3.49E-01 -6.48E-02 -2.02E-01  7.70E-02 -6.05E-03  4.52E+00
 
 TH 8
+        1.13E-03  6.66E-02  5.27E-02 -3.62E-02  3.65E-02 -2.78E-04  8.25E-02  1.21E-01
 
 TH 9
+        2.50E-04  5.12E-02 -5.28E-03 -3.15E-02  1.31E-02  7.71E-04  1.63E-02  2.05E-02  2.47E-02
 
 TH10
+       -3.42E-04 -8.57E-03 -2.25E-03  4.66E-03 -2.30E-03 -1.37E-03 -1.34E-02 -1.94E-02 -2.49E-03  2.23E-02
 
 TH11
+        1.55E-04 -3.50E-03  2.28E-03  2.46E-03 -3.27E-05  3.53E-04 -2.77E-02 -2.26E-04 -9.05E-04  3.55E-04  4.41E-03
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.07E-02
 
 TH 2
+       -6.22E-03  3.90E-01
 
 TH 3
+        1.04E-01 -1.40E-01  2.24E-01
 
 TH 4
+        2.44E-03 -9.87E-01  2.02E-01  2.39E-01
 
 TH 5
+        6.63E-02  7.85E-01  4.54E-01 -7.36E-01  1.31E-01
 
 TH 6
+        2.40E-03  7.95E-02 -9.09E-02 -6.58E-02  1.88E-02  6.28E-02
 
 TH 7
+       -6.42E-02  4.20E-01 -1.36E-01 -3.96E-01  2.77E-01 -4.53E-02  2.13E+00
 
 TH 8
+        1.06E-01  4.92E-01  6.78E-01 -4.36E-01  8.03E-01 -1.27E-02  1.12E-01  3.47E-01
 
 TH 9
+        5.19E-02  8.36E-01 -1.50E-01 -8.39E-01  6.38E-01  7.81E-02  4.87E-02  3.77E-01  1.57E-01
 
 TH10
+       -7.45E-02 -1.47E-01 -6.72E-02  1.30E-01 -1.18E-01 -1.46E-01 -4.24E-02 -3.74E-01 -1.06E-01  1.49E-01
 
 TH11
+        7.62E-02 -1.35E-01  1.53E-01  1.55E-01 -3.77E-03  8.47E-02 -1.96E-01 -9.78E-03 -8.68E-02  3.58E-02  6.64E-02
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.11E+03
 
 TH 2
+        9.29E+01  5.03E+02
 
 TH 3
+        2.26E+01  1.32E+02  2.31E+02
 
 TH 4
+        7.20E+01  5.32E+02 -3.34E+01  8.42E+02
 
 TH 5
+       -7.41E+01 -3.68E+02 -3.57E+02 -6.98E+01  1.06E+03
 
 TH 6
+       -3.38E+00 -4.62E+01  4.14E+00 -5.67E+01 -7.67E+00  2.73E+02
 
 TH 7
+       -1.05E+00 -5.89E+00 -9.47E-01 -2.97E+00  1.38E+00  1.30E+00  4.77E-01
 
 TH 8
+       -1.57E+01 -4.38E+01 -7.65E+01 -2.93E+00  5.16E+00  1.02E+01  1.00E+00  6.61E+01
 
 TH 9
+       -5.41E+01 -1.01E+02 -1.67E+01  7.24E+00  2.25E+01  1.21E+01  6.41E+00  1.65E+01  2.26E+02
 
 TH10
+        1.24E+01  3.86E+00 -2.40E+01  1.28E+01 -4.62E+01  2.18E+01  3.16E-01  3.66E+01  3.58E+00  7.04E+01
 
 TH11
+       -3.70E+01 -2.82E+01 -1.09E+01 -4.71E+01 -3.99E+01 -1.96E+01  1.75E+00  1.63E+01  1.32E+01  4.69E+00  2.53E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.365
Stop Time:
Wed Sep 29 19:01:23 CDT 2021
