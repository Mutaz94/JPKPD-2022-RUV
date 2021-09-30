Wed Sep 29 12:44:40 CDT 2021
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
$DATA ../../../../data/spa/A2/dat33.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -343.172744859204        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1771E+02  6.5562E+01  4.9224E+01  8.0381E+01  7.8030E+01  4.2921E+01 -2.3169E+01 -2.3581E+01 -8.7166E+00 -5.8546E+01
            -2.4284E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1228.74331253825        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0298E+00  1.0433E+00  9.7379E-01  1.0547E+00  1.0151E+00  8.1053E-01  9.6622E-01  9.9273E-01  8.6162E-01  9.3892E-01
             5.3309E+00
 PARAMETER:  1.2938E-01  1.4238E-01  7.3436E-02  1.5327E-01  1.1499E-01 -1.1007E-01  6.5639E-02  9.2704E-02 -4.8939E-02  3.6980E-02
             1.7735E+00
 GRADIENT:   9.3932E+01 -1.0390E+01 -1.4256E+01 -1.1792E+00 -5.0861E+00 -2.4841E+01  1.1785E+01  4.8398E+00  1.7383E+01  1.6673E+01
             2.6890E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1303.30241314486        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.5078E-01  7.8719E-01  1.5324E+00  1.1575E+00  1.0836E+00  8.8305E-01  5.8157E-01  5.1706E-01  8.1064E-01  4.7469E-01
             3.8768E+00
 PARAMETER:  4.9523E-02 -1.3929E-01  5.2684E-01  2.4629E-01  1.8027E-01 -2.4370E-02 -4.4202E-01 -5.5960E-01 -1.0993E-01 -6.4509E-01
             1.4550E+00
 GRADIENT:  -1.5230E+01 -5.6853E+00  1.1966E+00 -2.8793E+00 -1.1407E+01 -2.0576E+00 -5.5591E-02  7.5720E-01  7.1089E+00  3.3381E+00
             1.1064E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1311.64172888281        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.4993E-01  9.3481E-01  1.7087E+00  1.0763E+00  1.2204E+00  8.9703E-01  9.5580E-01  5.1325E-01  7.3629E-01  3.5315E-01
             3.3279E+00
 PARAMETER:  4.8636E-02  3.2585E-02  6.3571E-01  1.7355E-01  2.9919E-01 -8.6617E-03  5.4788E-02 -5.6700E-01 -2.0613E-01 -9.4087E-01
             1.3024E+00
 GRADIENT:   7.9562E+00 -1.6780E+00 -2.9415E+00  7.5337E+00  6.2621E+00  7.2201E-01  8.8164E-02  2.4908E-01  4.6525E-03  1.0489E+00
            -1.0945E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1312.38974000148        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  9.4754E-01  1.0820E+00  1.8702E+00  9.7535E-01  1.2817E+00  8.9368E-01  8.3947E-01  3.3100E-01  7.8412E-01  1.7452E-01
             3.3741E+00
 PARAMETER:  4.6116E-02  1.7882E-01  7.2603E-01  7.5041E-02  3.4815E-01 -1.2404E-02 -7.4987E-02 -1.0056E+00 -1.4320E-01 -1.6457E+00
             1.3161E+00
 GRADIENT:  -9.9383E-02 -2.0126E+00 -1.3968E-01 -4.6716E+00  6.7353E-01 -9.9009E-02  3.7413E-02  6.3923E-02  4.4912E-01  2.3434E-01
            -1.2122E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1313.03514518616        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      436
 NPARAMETR:  9.6547E-01  1.1029E+00  2.1627E+00  9.8048E-01  1.3500E+00  9.0070E-01  8.1210E-01  5.5773E-02  7.7373E-01  5.1791E-02
             3.4459E+00
 PARAMETER:  6.4864E-02  1.9794E-01  8.7137E-01  8.0289E-02  4.0012E-01 -4.5864E-03 -1.0814E-01 -2.7865E+00 -1.5654E-01 -2.8605E+00
             1.3372E+00
 GRADIENT:   6.3468E+00 -8.3809E-01 -7.5660E-01  3.9729E-01  1.7726E+00 -3.9031E-01 -1.4255E-01  1.1810E-03  3.2193E-02  1.8140E-02
             1.6383E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1313.21461299395        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      612
 NPARAMETR:  9.6530E-01  1.3495E+00  3.4527E+00  8.2394E-01  1.4868E+00  9.0279E-01  7.3881E-01  1.0000E-02  8.1814E-01  1.0000E-02
             3.4444E+00
 PARAMETER:  6.4683E-02  3.9974E-01  1.3392E+00 -9.3658E-02  4.9664E-01 -2.2631E-03 -2.0271E-01 -4.5179E+00 -1.0073E-01 -5.1771E+00
             1.3368E+00
 GRADIENT:   9.9691E-01 -2.0714E+00 -1.9291E-01 -1.3534E+00  1.8508E+00  1.9352E-02  1.6568E-01  1.6367E-06  1.2016E-01  0.0000E+00
             7.4368E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1313.26851450155        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      787
 NPARAMETR:  9.6559E-01  1.5336E+00  5.4573E+00  7.0572E-01  1.5347E+00  9.0290E-01  6.8825E-01  1.0000E-02  8.8130E-01  1.0000E-02
             3.4411E+00
 PARAMETER:  6.4983E-02  5.2761E-01  1.7970E+00 -2.4853E-01  5.2832E-01 -2.1410E-03 -2.7361E-01 -6.7089E+00 -2.6361E-02 -7.8827E+00
             1.3358E+00
 GRADIENT:  -1.9531E+00  3.2969E+00 -4.8864E-02  2.4584E+00 -5.8912E-01 -3.4680E-01 -1.8013E-01  0.0000E+00  7.1750E-03  0.0000E+00
            -8.9938E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1313.29644755545        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      962
 NPARAMETR:  9.6744E-01  1.6636E+00  9.0208E+00  6.1904E-01  1.5634E+00  9.0424E-01  6.7139E-01  1.0000E-02  9.1131E-01  1.0000E-02
             3.4475E+00
 PARAMETER:  6.6896E-02  6.0896E-01  2.2995E+00 -3.7958E-01  5.4689E-01 -6.6055E-04 -2.9841E-01 -9.5003E+00  7.1229E-03 -1.0996E+01
             1.3377E+00
 GRADIENT:   7.2305E-01  2.0779E+00 -4.0262E-02  1.2919E+00  5.4305E-01 -1.3451E-03  7.5069E-02  0.0000E+00 -7.8372E-02  0.0000E+00
             3.9204E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1313.31259002814        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1140
 NPARAMETR:  9.6714E-01  1.7904E+00  1.7002E+01  5.3133E-01  1.5678E+00  9.0411E-01  6.4744E-01  1.0000E-02  9.6737E-01  1.0000E-02
             3.4444E+00
 PARAMETER:  6.6583E-02  6.8245E-01  2.9333E+00 -5.3237E-01  5.4965E-01 -8.0699E-04 -3.3473E-01 -1.3353E+01  6.6821E-02 -1.5145E+01
             1.3367E+00
 GRADIENT:  -6.0506E-01  2.5480E+00 -1.1031E-02  1.2540E+00 -8.1249E-01 -9.4959E-02 -2.4576E-01  0.0000E+00 -2.5272E-01  0.0000E+00
            -6.1213E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1313.31615566895        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1315
 NPARAMETR:  9.6708E-01  1.8687E+00  2.7662E+01  4.7808E-01  1.5712E+00  9.0399E-01  6.3372E-01  1.0000E-02  1.0216E+00  1.0000E-02
             3.4438E+00
 PARAMETER:  6.6526E-02  7.2522E-01  3.4201E+00 -6.3797E-01  5.5184E-01 -9.3640E-04 -3.5615E-01 -1.6393E+01  1.2141E-01 -1.8343E+01
             1.3366E+00
 GRADIENT:  -1.2548E+00  4.1297E+00 -5.7246E-03  1.7195E+00 -1.0950E+00 -1.9859E-01 -1.6919E-01  0.0000E+00 -1.8184E-01  0.0000E+00
            -7.4628E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1313.32161690240        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1490
 NPARAMETR:  9.6716E-01  1.9589E+00  5.5869E+01  4.1535E-01  1.5741E+00  9.0401E-01  6.1591E-01  1.0000E-02  1.1156E+00  1.0000E-02
             3.4441E+00
 PARAMETER:  6.6613E-02  7.7240E-01  4.1230E+00 -7.7864E-01  5.5372E-01 -9.1699E-04 -3.8466E-01 -2.0890E+01  2.0937E-01 -2.2989E+01
             1.3367E+00
 GRADIENT:  -1.2615E+00  3.8184E+00 -2.2095E-03  1.4443E+00 -8.9952E-01 -2.0220E-01 -1.0044E-01  0.0000E+00 -4.7663E-02  0.0000E+00
            -5.7599E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1313.32704636818        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1673             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6749E-01  2.0026E+00  1.6701E+02  3.8241E-01  1.5769E+00  9.0433E-01  6.0720E-01  1.0000E-02  1.1746E+00  1.0000E-02
             3.4453E+00
 PARAMETER:  6.6950E-02  7.9445E-01  5.2181E+00 -8.6126E-01  5.5547E-01 -5.5977E-04 -3.9889E-01 -2.3606E+01  2.6091E-01 -2.5764E+01
             1.3370E+00
 GRADIENT:   3.3653E+01  7.3220E+01 -5.9961E-04  6.4296E+00  1.8826E+00  2.7209E+00  1.4235E+00  0.0000E+00  2.8331E-01  0.0000E+00
             1.0703E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1313.32797903066        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1849
 NPARAMETR:  9.6756E-01  2.0058E+00  1.4908E+03  3.8069E-01  1.5764E+00  9.0435E-01  6.0822E-01  1.0000E-02  1.1650E+00  1.0000E-02
             3.4460E+00
 PARAMETER:  6.7027E-02  7.9602E-01  7.4071E+00 -8.6577E-01  5.5514E-01 -5.4204E-04 -3.9721E-01 -2.3606E+01  2.5271E-01 -2.5764E+01
             1.3372E+00
 GRADIENT:   1.2034E-01 -1.0895E+00 -6.7970E-05  6.2252E-03  3.3953E-02  2.5250E-02  1.5939E-02  0.0000E+00 -6.3337E-03  0.0000E+00
             1.0443E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1313.32803068876        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2034            RESET HESSIAN, TYPE II
 NPARAMETR:  9.6755E-01  2.0063E+00  2.7349E+03  3.8016E-01  1.5762E+00  9.0431E-01  6.0814E-01  1.0000E-02  1.1675E+00  1.0000E-02
             3.4453E+00
 PARAMETER:  6.7017E-02  7.9629E-01  8.0139E+00 -8.6718E-01  5.5502E-01 -5.8687E-04 -3.9735E-01 -2.3606E+01  2.5490E-01 -2.5764E+01
             1.3370E+00
 GRADIENT:   3.3698E+01  7.4474E+01 -3.4640E-05  6.6171E+00  1.5577E+00  2.7018E+00  1.4388E+00  0.0000E+00  2.3440E-01  0.0000E+00
             1.0653E+01

0ITERATION NO.:   72    OBJECTIVE VALUE:  -1313.32803068876        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     2091
 NPARAMETR:  9.6755E-01  2.0063E+00  2.7349E+03  3.8016E-01  1.5762E+00  9.0431E-01  6.0814E-01  1.0000E-02  1.1675E+00  1.0000E-02
             3.4453E+00
 PARAMETER:  6.7017E-02  7.9629E-01  8.0139E+00 -8.6718E-01  5.5502E-01 -5.8687E-04 -3.9735E-01 -2.3606E+01  2.5490E-01 -2.5764E+01
             1.3370E+00
 GRADIENT:   1.1580E-01 -1.3166E+00 -3.6146E-05 -5.2273E-02  4.5540E-02  1.5911E-02  4.0151E-02  0.0000E+00  4.9769E-03  0.0000E+00
            -4.3244E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2091
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1740E-03 -9.4210E-03  1.9101E-09 -4.9272E-03 -4.3823E-05
 SE:             2.8631E-02  2.0763E-02  1.4907E-09  1.0791E-02  1.0980E-04
 N:                     100         100         100         100         100

 P VAL.:         9.6729E-01  6.5002E-01  2.0007E-01  6.4795E-01  6.8979E-01

 ETASHRINKSD(%)  4.0817E+00  3.0440E+01  1.0000E+02  6.3850E+01  9.9632E+01
 ETASHRINKVR(%)  7.9967E+00  5.1614E+01  1.0000E+02  8.6932E+01  9.9999E+01
 EBVSHRINKSD(%)  4.0090E+00  3.0390E+01  1.0000E+02  6.3826E+01  9.9589E+01
 EBVSHRINKVR(%)  7.8572E+00  5.1544E+01  1.0000E+02  8.6914E+01  9.9998E+01
 RELATIVEINF(%)  9.8618E+01  4.3658E-06  1.0000E-10  1.1792E-06  1.0000E-10
 EPSSHRINKSD(%)  1.9742E+01
 EPSSHRINKVR(%)  3.5586E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1313.3280306887623     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -578.17720412502410     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.75
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1313.328       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  2.01E+00  2.73E+03  3.80E-01  1.58E+00  9.04E-01  6.08E-01  1.00E-02  1.17E+00  1.00E-02  3.45E+00
 


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
+       -9.06E+01  3.56E+02
 
 TH 3
+       -2.40E-06 -6.10E-07 -1.08E-11
 
 TH 4
+       -1.29E+02  5.07E+02 -8.21E-07  7.28E+02
 
 TH 5
+       -7.54E+00 -7.40E+01 -4.27E-08 -1.06E+02  7.21E+01
 
 TH 6
+        2.71E+00 -1.82E+01  3.98E-07 -2.58E+01 -2.72E+00  2.05E+02
 
 TH 7
+        4.77E+00 -1.36E+01  1.52E-07 -2.34E+01  7.67E+00  6.28E+00  1.22E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.31E-01 -2.32E+00  2.48E-07 -1.24E+00  1.08E+00  5.95E-02  1.46E+01  0.00E+00  3.29E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.80E+01 -9.92E+00  1.14E-07 -1.42E+01  8.92E-01  3.35E+00  1.58E+01  0.00E+00  2.23E+00  0.00E+00  3.79E+01
 
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
 #CPUT: Total CPU Time in Seconds,       32.137
Stop Time:
Wed Sep 29 12:45:15 CDT 2021
