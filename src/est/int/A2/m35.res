Fri Sep 24 21:22:23 CDT 2021
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
$DATA ../../../../data/int/A2/dat35.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 (2E4.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2953.85020121856        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -4.7424E+01  3.6343E+01  1.1898E+02  9.8701E+00  1.3192E+02  3.0590E+01 -9.5374E+01 -7.9812E+01 -1.0018E+01 -7.8987E+01
            -1.6378E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3331.92048974941        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0689E+00  1.0003E+00  8.9527E-01  1.0712E+00  8.9466E-01  8.6346E-01  1.0778E+00  8.0890E-01  1.0167E+00  1.0255E+00
             1.5337E+00
 PARAMETER:  1.6661E-01  1.0028E-01 -1.0630E-02  1.6881E-01 -1.1311E-02 -4.6810E-02  1.7493E-01 -1.1209E-01  1.1661E-01  1.2517E-01
             5.2765E-01
 GRADIENT:   9.4130E+01  6.6572E+01  1.1454E+01  8.9158E+01 -1.4473E+01 -2.2116E+01 -4.5841E+00  9.9982E+00 -3.4253E-01  2.5316E+00
            -1.7938E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3339.23203982093        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0697E+00  8.0160E-01  6.6640E-01  1.1321E+00  6.8014E-01  9.1195E-01  1.3090E+00  4.0986E-01  1.0241E+00  6.9075E-01
             1.5143E+00
 PARAMETER:  1.6742E-01 -1.2114E-01 -3.0587E-01  2.2405E-01 -2.8546E-01  7.8334E-03  3.6924E-01 -7.9194E-01  1.2379E-01 -2.6998E-01
             5.1498E-01
 GRADIENT:   9.2088E+01  8.2298E+01  3.9679E+00  5.5099E+01 -1.3478E+01 -8.7015E-01 -1.1055E+01 -1.2005E+00 -8.3793E-01 -1.9184E+01
            -1.7510E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3352.14657353634        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0284E+00  5.7798E-01  5.0555E-01  1.2043E+00  5.1168E-01  9.0886E-01  1.4427E+00  4.0074E-01  1.0292E+00  5.8533E-01
             1.6220E+00
 PARAMETER:  1.2801E-01 -4.4821E-01 -5.8211E-01  2.8592E-01 -5.7006E-01  4.4315E-03  4.6652E-01 -8.1445E-01  1.2874E-01 -4.3559E-01
             5.8366E-01
 GRADIENT:  -2.4984E+01  2.3542E+01  4.2701E+01  8.3872E+01 -2.1973E+01 -2.7729E+00 -6.4455E+00 -8.2541E-01  1.0582E+01 -9.4381E+00
             2.9207E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3352.51481876919        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0261E+00  4.6323E-01  3.9317E-01  1.2216E+00  4.0941E-01  9.1781E-01  1.5636E+00  3.9726E-01  1.0411E+00  5.2930E-01
             1.6125E+00
 PARAMETER:  1.2572E-01 -6.6952E-01 -8.3351E-01  3.0018E-01 -7.9304E-01  1.4237E-02  5.4697E-01 -8.2317E-01  1.4030E-01 -5.3619E-01
             5.7780E-01
 GRADIENT:  -3.0612E+01  1.4486E+01  6.5825E+01  1.0417E+02 -4.3900E+01 -9.0226E-01  2.6997E+00 -4.7109E+00  1.4797E+01 -1.2226E+01
             5.0075E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3356.95766846952        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  1.0442E+00  4.1468E-01  3.4858E-01  1.1606E+00  3.7683E-01  9.3545E-01  1.6716E+00  3.2867E-01  9.9962E-01  6.1143E-01
             1.5903E+00
 PARAMETER:  1.4322E-01 -7.8025E-01 -9.5387E-01  2.4890E-01 -8.7597E-01  3.3267E-02  6.1375E-01 -1.0127E+00  9.9616E-02 -3.9195E-01
             5.6394E-01
 GRADIENT:  -3.5466E+00 -1.7474E+01  2.3989E+01  4.1445E+00 -8.3332E+00  5.7165E+00  1.3508E+01 -2.1412E+00  2.8204E+00  8.1713E-01
             4.4530E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3363.54034036155        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      637
 NPARAMETR:  1.0346E+00  2.6223E-01  1.6618E-01  1.0031E+00  2.1868E-01  9.3952E-01  1.6861E+00  2.7522E-01  1.0856E+00  7.0608E-01
             1.5065E+00
 PARAMETER:  1.3403E-01 -1.2385E+00 -1.6947E+00  1.0306E-01 -1.4201E+00  3.7616E-02  6.2239E-01 -1.1902E+00  1.8209E-01 -2.4803E-01
             5.0978E-01
 GRADIENT:  -2.8898E+01  2.1300E+01  3.1486E+01  3.6327E+00 -4.6664E+01  2.9488E+00 -1.7187E+00 -3.7881E+00  1.0743E+01  4.8565E+00
             3.4932E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3371.09150024306        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      812
 NPARAMETR:  1.0441E+00  1.7844E-01  8.8283E-02  7.7549E-01  1.5399E-01  9.3650E-01  1.8232E+00  2.5260E-01  1.2410E+00  7.7181E-01
             1.4502E+00
 PARAMETER:  1.4311E-01 -1.6235E+00 -2.3272E+00 -1.5426E-01 -1.7709E+00  3.4391E-02  7.0060E-01 -1.2760E+00  3.1592E-01 -1.5902E-01
             4.7167E-01
 GRADIENT:  -8.6708E+00 -3.3255E+00  5.8940E+00  1.3088E+01 -2.3538E+01  1.8215E+00  7.1276E+00 -3.6984E+00 -1.6374E-01 -1.1354E+00
             1.1485E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3377.22305963071        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      991
 NPARAMETR:  1.0479E+00  1.6551E-01  7.3120E-02  6.8761E-01  1.4344E-01  9.3099E-01  1.7682E+00  7.0915E-01  1.4309E+00  7.3726E-01
             1.4277E+00
 PARAMETER:  1.4682E-01 -1.6987E+00 -2.5156E+00 -2.7454E-01 -1.8418E+00  2.8497E-02  6.6996E-01 -2.4369E-01  4.5827E-01 -2.0482E-01
             4.5606E-01
 GRADIENT:   3.0103E+00  3.4385E+00 -8.0204E+00 -2.0823E+01  3.1164E+01  2.5114E-01  1.0077E+00 -8.5076E+00  1.1574E+01 -3.1345E+00
             1.8977E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3377.34745127172        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1171
 NPARAMETR:  1.0480E+00  1.6533E-01  7.3186E-02  6.8771E-01  1.4327E-01  9.2809E-01  1.7683E+00  7.1445E-01  1.4316E+00  7.4832E-01
             1.4270E+00
 PARAMETER:  1.4686E-01 -1.6998E+00 -2.5148E+00 -2.7438E-01 -1.8430E+00  2.5377E-02  6.7004E-01 -2.3624E-01  4.5876E-01 -1.8993E-01
             4.5555E-01
 GRADIENT:   2.4249E+01  2.7157E+01  1.5497E+01 -1.6181E+01  2.2311E+02  3.4543E-01  3.7338E+00 -7.8194E+00  1.3376E+01  2.0646E+00
             1.9356E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3377.35923202272        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1306
 NPARAMETR:  1.0480E+00  1.6533E-01  7.3186E-02  6.8772E-01  1.4322E-01  9.3256E-01  1.7683E+00  7.1445E-01  1.4316E+00  7.4488E-01
             1.4270E+00
 PARAMETER:  1.4686E-01 -1.6998E+00 -2.5147E+00 -2.7438E-01 -1.8434E+00  3.0179E-02  6.7004E-01 -2.3624E-01  4.5876E-01 -1.9453E-01
             4.5554E-01
 GRADIENT:   3.4402E+00  2.9092E+00 -4.3639E+00 -2.0530E+01  2.4979E+01  7.9428E-01  9.2076E-01 -8.0388E+00  1.1954E+01 -1.0359E-01
             1.8465E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3377.67654410732        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1492             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0480E+00  1.6528E-01  7.3191E-02  6.8773E-01  1.4282E-01  9.3653E-01  1.7680E+00  7.4022E-01  1.4315E+00  7.4864E-01
             1.4269E+00
 PARAMETER:  1.4686E-01 -1.7001E+00 -2.5147E+00 -2.7436E-01 -1.8462E+00  3.4424E-02  6.6987E-01 -2.0080E-01  4.5875E-01 -1.8949E-01
             4.5549E-01
 GRADIENT:   2.4797E+01  2.7563E+01  2.0082E+01 -1.5910E+01  2.1324E+02  3.6852E+00  4.0273E+00 -7.0861E+00  1.2951E+01  3.2984E+00
             2.1011E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3377.96153286244        NO. OF FUNC. EVALS.:  84
 CUMULATIVE NO. OF FUNC. EVALS.:     1576
 NPARAMETR:  1.0480E+00  1.6528E-01  7.3191E-02  6.8773E-01  1.4282E-01  9.3570E-01  1.7644E+00  7.7042E-01  1.4315E+00  7.4577E-01
             1.4269E+00
 PARAMETER:  1.4686E-01 -1.7001E+00 -2.5147E+00 -2.7436E-01 -1.8462E+00  3.3540E-02  6.6780E-01 -1.6082E-01  4.5875E-01 -1.9333E-01
             4.5549E-01
 GRADIENT:   2.4864E+01  2.7616E+01  2.0773E+01 -1.6475E+01  2.1312E+02  3.3551E+00  3.8627E+00 -6.3004E+00  1.2459E+01  3.6271E+00
             2.2978E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3378.00481188894        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1712
 NPARAMETR:  1.0480E+00  1.6513E-01  7.3280E-02  6.8783E-01  1.4263E-01  9.3121E-01  1.7650E+00  7.7048E-01  1.4319E+00  7.3793E-01
             1.4265E+00
 PARAMETER:  1.4693E-01 -1.7010E+00 -2.5135E+00 -2.7422E-01 -1.8475E+00  2.8726E-02  6.6812E-01 -1.6074E-01  4.5898E-01 -2.0390E-01
             4.5526E-01
 GRADIENT:   2.5023E+01  2.7037E+01  2.2871E+01 -1.6465E+01  2.0777E+02  1.5218E+00  4.1579E+00 -6.6255E+00  1.2133E+01  1.0852E+00
             2.2303E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3378.00824187329        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:     1792
 NPARAMETR:  1.0481E+00  1.6509E-01  7.3269E-02  6.8783E-01  1.4235E-01  9.3464E-01  1.7650E+00  7.7048E-01  1.4319E+00  7.4364E-01
             1.4265E+00
 PARAMETER:  1.4693E-01 -1.7013E+00 -2.5136E+00 -2.7421E-01 -1.8495E+00  3.2405E-02  6.6813E-01 -1.6075E-01  4.5898E-01 -1.9620E-01
             4.5524E-01
 GRADIENT:   2.5360E+01  2.7307E+01  2.5403E+01 -1.5847E+01  2.0223E+02  2.8380E+00  4.0155E+00 -6.3597E+00  1.2331E+01  2.9379E+00
             2.2310E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -3378.51898459374        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1991
 NPARAMETR:  1.0477E+00  1.6492E-01  7.3258E-02  6.9153E-01  1.4237E-01  9.3377E-01  1.7634E+00  8.0121E-01  1.4150E+00  7.4155E-01
             1.4227E+00
 PARAMETER:  1.4663E-01 -1.7023E+00 -2.5138E+00 -2.6885E-01 -1.8493E+00  3.1475E-02  6.6724E-01 -1.2163E-01  4.4711E-01 -1.9901E-01
             4.5252E-01
 GRADIENT:   2.7415E+00  2.0514E+00  3.2519E+00 -1.6103E+01  2.3516E+00  1.4384E+00  1.0358E+00 -5.5780E+00  7.9251E+00  2.0789E+00
             1.8131E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -3378.53290488752        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2182
 NPARAMETR:  1.0478E+00  1.6484E-01  7.3289E-02  6.9159E-01  1.4231E-01  9.2865E-01  1.7637E+00  8.0119E-01  1.4151E+00  7.3431E-01
             1.4224E+00
 PARAMETER:  1.4667E-01 -1.7028E+00 -2.5133E+00 -2.6877E-01 -1.8498E+00  2.5971E-02  6.6741E-01 -1.2166E-01  4.4721E-01 -2.0882E-01
             4.5235E-01
 GRADIENT:   2.8638E+00  1.7053E+00  3.7091E+00 -1.6097E+01  1.5511E-01 -7.0267E-01  1.2543E+00 -5.9046E+00  7.5808E+00 -3.3162E-01
             1.7614E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -3378.55899232213        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     2341
 NPARAMETR:  1.0478E+00  1.6472E-01  7.3181E-02  6.9179E-01  1.4231E-01  9.2576E-01  1.7635E+00  8.0226E-01  1.4145E+00  7.3442E-01
             1.4217E+00
 PARAMETER:  1.4666E-01 -1.7035E+00 -2.5148E+00 -2.6847E-01 -1.8498E+00  2.2857E-02  6.6732E-01 -1.2033E-01  4.4680E-01 -2.0868E-01
             4.5184E-01
 GRADIENT:   2.4575E+01  2.5954E+01  2.2719E+01 -1.1071E+01  2.0145E+02 -7.8391E-01  4.0710E+00 -5.8728E+00  8.6454E+00  7.1215E-01
             1.7473E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -3378.56294776724        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:     2492
 NPARAMETR:  1.0478E+00  1.6458E-01  7.3273E-02  6.9189E-01  1.4217E-01  9.3087E-01  1.7641E+00  8.0221E-01  1.4148E+00  7.3429E-01
             1.4214E+00
 PARAMETER:  1.4674E-01 -1.7043E+00 -2.5136E+00 -2.6834E-01 -1.8507E+00  2.8366E-02  6.6766E-01 -1.2039E-01  4.4702E-01 -2.0886E-01
             4.5161E-01
 GRADIENT:   3.0723E+00  9.1828E-01  4.7352E+00 -1.5311E+01 -2.6007E+00  2.2550E-01  1.2859E+00 -5.9151E+00  7.4331E+00 -3.6570E-01
             1.6225E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -3378.68347857754        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:     2638
 NPARAMETR:  1.0478E+00  1.6458E-01  7.3262E-02  6.9209E-01  1.4217E-01  9.3024E-01  1.7639E+00  8.1801E-01  1.4143E+00  7.3557E-01
             1.4211E+00
 PARAMETER:  1.4673E-01 -1.7043E+00 -2.5137E+00 -2.6804E-01 -1.8507E+00  2.7687E-02  6.6752E-01 -1.0088E-01  4.4663E-01 -2.0711E-01
             4.5141E-01
 GRADIENT:   2.4804E+01  2.5472E+01  2.5306E+01 -1.1025E+01  1.9746E+02  1.1074E+00  4.3114E+00 -5.2283E+00  8.4785E+00  1.7499E+00
             1.7696E+01

0ITERATION NO.:   99    OBJECTIVE VALUE:  -3378.68591659172        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     2757
 NPARAMETR:  1.0478E+00  1.6458E-01  7.3265E-02  6.9211E-01  1.4217E-01  9.3025E-01  1.7639E+00  8.1826E-01  1.4143E+00  7.3539E-01
             1.4210E+00
 PARAMETER:  1.4673E-01 -1.7043E+00 -2.5137E+00 -2.6802E-01 -1.8507E+00  2.7688E-02  6.6750E-01 -1.0058E-01  4.4661E-01 -2.0715E-01
             4.5140E-01
 GRADIENT:  -2.1256E+04  9.1609E+02 -1.2281E+03 -5.8376E+03  1.6864E+03 -2.3754E-02 -4.6708E+03  3.0998E+04 -6.9829E+03  6.8788E-01
             3.4744E+03
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         2.9         3.3         3.3         3.3         1.6
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2757
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.8771E-04  1.1839E-02  2.7019E-02  1.5583E-02  2.0482E-02
 SE:             2.9772E-02  2.7830E-02  1.6240E-02  2.5731E-02  2.6448E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8425E-01  6.7055E-01  9.6173E-02  5.4476E-01  4.3868E-01

 ETASHRINKSD(%)  2.5846E-01  6.7657E+00  4.5593E+01  1.3799E+01  1.1394E+01
 ETASHRINKVR(%)  5.1626E-01  1.3074E+01  7.0399E+01  2.5694E+01  2.1490E+01
 EBVSHRINKSD(%)  6.1007E-01  5.3349E+00  4.8201E+01  9.1215E+00  1.1875E+01
 EBVSHRINKVR(%)  1.2164E+00  1.0385E+01  7.3169E+01  1.7411E+01  2.2340E+01
 RELATIVEINF(%)  9.8767E+01  4.2529E+01  3.4872E+00  2.4201E+01  1.2215E+01
 EPSSHRINKSD(%)  2.3562E+01
 EPSSHRINKVR(%)  4.1573E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3378.6859165917235     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1724.5965568233128     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    68.29
 Elapsed covariance  time in seconds:    12.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3378.686       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.65E-01  7.33E-02  6.92E-01  1.42E-01  9.30E-01  1.76E+00  8.18E-01  1.41E+00  7.36E-01  1.42E+00
 


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
 
         6.76E-04  9.74E-03  7.87E-04  1.42E-01  5.89E-04  7.41E-02  1.23E-01  8.17E-02  2.08E-02  2.27E-01  1.34E-01
 


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
+        4.56E-07
 
 TH 2
+       -3.11E-06  9.49E-05
 
 TH 3
+        5.31E-07 -3.64E-06  6.19E-07
 
 TH 4
+       -9.37E-05  5.87E-04 -1.09E-04  2.02E-02
 
 TH 5
+       -3.98E-07  2.73E-06 -4.64E-07  8.16E-05  3.47E-07
 
 TH 6
+        3.40E-05 -2.39E-04  3.96E-05 -6.49E-03 -2.96E-05  5.49E-03
 
 TH 7
+       -2.60E-05  4.70E-05 -3.03E-05  5.73E-03  2.26E-05 -1.82E-03  1.52E-02
 
 TH 8
+       -5.35E-05  3.31E-04 -6.23E-05  1.15E-02  4.66E-05 -3.73E-03  4.14E-03  6.67E-03
 
 TH 9
+       -5.13E-06  2.01E-04 -6.02E-06  9.41E-04  4.52E-06 -3.99E-04 -4.47E-07  5.27E-04  4.32E-04
 
 TH10
+        1.50E-04 -1.03E-03  1.74E-04 -3.03E-02 -1.31E-04  1.19E-02 -7.74E-03 -1.73E-02 -1.70E-03  5.14E-02
 
 TH11
+        7.75E-05 -5.14E-04  9.05E-05 -1.63E-02 -6.76E-05  5.57E-03 -4.11E-03 -9.85E-03 -8.40E-04  2.55E-02  1.80E-02
 
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
+        6.76E-04
 
 TH 2
+       -4.72E-01  9.74E-03
 
 TH 3
+        1.00E+00 -4.75E-01  7.87E-04
 
 TH 4
+       -9.75E-01  4.24E-01 -9.75E-01  1.42E-01
 
 TH 5
+       -1.00E+00  4.76E-01 -1.00E+00  9.74E-01  5.89E-04
 
 TH 6
+        6.79E-01 -3.31E-01  6.79E-01 -6.16E-01 -6.78E-01  7.41E-02
 
 TH 7
+       -3.13E-01  3.92E-02 -3.13E-01  3.27E-01  3.11E-01 -1.99E-01  1.23E-01
 
 TH 8
+       -9.70E-01  4.16E-01 -9.70E-01  9.90E-01  9.68E-01 -6.17E-01  4.12E-01  8.17E-02
 
 TH 9
+       -3.66E-01  9.93E-01 -3.68E-01  3.19E-01  3.69E-01 -2.59E-01 -1.75E-04  3.11E-01  2.08E-02
 
 TH10
+        9.77E-01 -4.66E-01  9.77E-01 -9.40E-01 -9.78E-01  7.09E-01 -2.77E-01 -9.35E-01 -3.62E-01  2.27E-01
 
 TH11
+        8.56E-01 -3.94E-01  8.58E-01 -8.55E-01 -8.55E-01  5.61E-01 -2.49E-01 -8.99E-01 -3.02E-01  8.39E-01  1.34E-01
 
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
+        1.60E+12
 
 TH 2
+       -2.01E+10  7.46E+11
 
 TH 3
+       -1.53E+10  1.35E+11  4.99E+10
 
 TH 4
+       -1.41E+12  1.10E+12  2.02E+11  2.83E+12
 
 TH 5
+       -9.89E+10  1.32E+11  4.88E+10  2.78E+11  6.22E+10
 
 TH 6
+        9.33E+06 -2.78E+06 -1.47E+06 -1.23E+07 -2.48E+06  5.19E+02
 
 TH 7
+       -2.23E+11  1.74E+11  3.18E+10  4.46E+11  4.38E+10 -1.94E+06  7.03E+10
 
 TH 8
+        3.19E+12 -2.48E+12 -4.56E+11 -6.38E+12 -6.27E+11  2.78E+07 -1.01E+12  1.44E+13
 
 TH 9
+        9.33E+09 -3.32E+11 -5.98E+10 -4.90E+11 -5.86E+10  1.24E+06 -7.72E+10  1.10E+12  1.47E+11
 
 TH10
+       -1.79E+06  6.81E+06  2.19E+06  1.13E+07  2.64E+06 -1.94E+02  1.79E+06 -2.56E+07 -3.03E+06  6.24E+02
 
 TH11
+        4.09E+11 -3.18E+11 -5.84E+10 -8.18E+11 -8.03E+10  3.56E+06 -1.29E+11  1.84E+12  1.42E+11 -3.28E+06  2.36E+11
 
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
 #CPUT: Total CPU Time in Seconds,       80.905
Stop Time:
Fri Sep 24 21:23:46 CDT 2021
