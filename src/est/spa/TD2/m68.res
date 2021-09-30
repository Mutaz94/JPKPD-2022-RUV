Wed Sep 29 19:13:16 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1708.00252217131        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0844E+02  4.6943E+01 -2.4724E+01  1.3104E+02  5.2471E+01  4.8798E+01  2.1172E+01  4.0586E+00  4.5505E+01  5.2094E-01
             2.0994E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1713.43052367759        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0097E+00  1.0188E+00  1.0077E+00  9.6765E-01  9.5876E-01  9.9120E-01  8.4357E-01  9.8144E-01  7.2465E-01  9.8679E-01
             9.1486E-01
 PARAMETER:  1.0961E-01  1.1858E-01  1.0767E-01  6.7117E-02  5.7883E-02  9.1161E-02 -7.0117E-02  8.1264E-02 -2.2207E-01  8.6704E-02
             1.1018E-02
 GRADIENT:  -1.7213E+01  3.5631E+01  2.1026E+01  9.6788E+00 -2.9965E+01 -1.6820E+00 -7.6426E+00 -2.8951E+00 -2.6257E+01 -2.0150E+00
            -2.4100E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1714.77733230638        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0108E+00  9.0228E-01  8.7923E-01  1.0274E+00  8.5862E-01  9.9570E-01  6.9331E-01  6.7852E-01  8.0344E-01  9.3539E-01
             9.0817E-01
 PARAMETER:  1.1069E-01 -2.8292E-03 -2.8714E-02  1.2707E-01 -5.2430E-02  9.5694E-02 -2.6628E-01 -2.8784E-01 -1.1885E-01  3.3205E-02
             3.6735E-03
 GRADIENT:  -1.6004E+01  2.0425E+01  1.2256E+01  1.5584E+01 -2.0611E+01 -3.4865E-01 -7.5315E+00 -2.6329E+00 -1.0735E+01 -5.0943E+00
            -2.6420E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1716.94346238079        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  1.0182E+00  8.8336E-01  8.7292E-01  1.0333E+00  8.5877E-01  9.9655E-01  9.2859E-01  7.0804E-01  7.7445E-01  9.1666E-01
             9.6246E-01
 PARAMETER:  1.1803E-01 -2.4025E-02 -3.5908E-02  1.3277E-01 -5.2256E-02  9.6542E-02  2.5913E-02 -2.4525E-01 -1.5560E-01  1.2979E-02
             6.1735E-02
 GRADIENT:   1.2812E-01  2.8490E-01 -1.3818E+00  7.9309E-01  1.3345E+00  4.0505E-01 -1.4670E-01  2.4418E-01  1.9424E-01  4.8527E-01
             1.0592E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1716.99816577278        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0157E+00  7.6035E-01  9.3622E-01  1.1123E+00  8.3527E-01  9.9362E-01  1.0497E+00  7.3316E-01  7.2896E-01  9.0503E-01
             9.5852E-01
 PARAMETER:  1.1561E-01 -1.7397E-01  3.4094E-02  2.0639E-01 -8.0006E-02  9.3602E-02  1.4846E-01 -2.1039E-01 -2.1614E-01  2.1104E-04
             5.7640E-02
 GRADIENT:  -2.1726E+00  5.2713E+00  6.2401E+00  1.6951E+00 -8.8611E+00 -2.2452E-01  2.7986E-01 -1.0750E+00 -6.0105E-02 -1.2937E+00
            -1.5751E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1717.45701854868        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      906
 NPARAMETR:  1.0090E+00  4.9118E-01  1.2413E+00  1.2897E+00  8.6678E-01  9.8884E-01  1.2567E+00  1.0963E+00  6.6395E-01  9.5801E-01
             9.6127E-01
 PARAMETER:  1.0897E-01 -6.1094E-01  3.1620E-01  3.5444E-01 -4.2968E-02  8.8778E-02  3.2848E-01  1.9195E-01 -3.0955E-01  5.7103E-02
             6.0498E-02
 GRADIENT:  -5.5369E+00  3.6020E+00  6.2732E+00 -6.7785E+00 -1.7375E+01  8.3700E-02 -3.5798E-01  1.7204E+00 -1.5242E+00  1.1915E+00
             4.8563E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1717.96992972412        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1084
 NPARAMETR:  1.0100E+00  2.9965E-01  1.5377E+00  1.4304E+00  9.2560E-01  9.8496E-01  1.4656E+00  1.3054E+00  6.2766E-01  1.0267E+00
             9.6115E-01
 PARAMETER:  1.0992E-01 -1.1051E+00  5.3032E-01  4.5796E-01  2.2691E-02  8.4849E-02  4.8227E-01  3.6647E-01 -3.6575E-01  1.2635E-01
             6.0377E-02
 GRADIENT:   4.6991E+00  5.1661E+00 -7.6950E-01  2.3542E+01  5.2144E+00  2.2371E-01 -3.9896E-01 -2.8021E+00 -1.5594E+00 -6.3359E-02
            -5.0833E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1718.08828602253        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1259
 NPARAMETR:  1.0090E+00  2.1476E-01  1.7180E+00  1.4964E+00  9.5511E-01  9.8134E-01  1.6229E+00  1.5174E+00  6.1371E-01  1.0484E+00
             9.5882E-01
 PARAMETER:  1.0900E-01 -1.4383E+00  6.4117E-01  5.0306E-01  5.4072E-02  8.1168E-02  5.8419E-01  5.1697E-01 -3.8824E-01  1.4726E-01
             5.7953E-02
 GRADIENT:   5.5811E+00  5.7915E+00 -5.8476E+00  4.3293E+01  8.1784E+00 -5.6295E-01 -2.5759E-01  1.4881E-01 -5.1828E-01  3.6215E-01
            -5.6005E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1718.65762682918        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1436
 NPARAMETR:  1.0049E+00  9.7435E-02  1.8529E+00  1.5728E+00  9.4924E-01  9.7824E-01  2.2258E+00  1.6738E+00  5.9100E-01  1.0407E+00
             9.5747E-01
 PARAMETER:  1.0491E-01 -2.2286E+00  7.1673E-01  5.5287E-01  4.7907E-02  7.7997E-02  9.0013E-01  6.1511E-01 -4.2593E-01  1.3986E-01
             5.6542E-02
 GRADIENT:   2.6059E-01  3.0402E+00 -4.5234E+00  4.2592E+01 -1.0207E+00 -9.5011E-01 -1.4249E-01  2.3999E+00 -5.4789E-01  8.1790E-01
            -7.4009E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1719.33668050527        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1619             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0045E+00  2.0130E-02  2.0539E+00  1.6074E+00  9.7432E-01  9.7942E-01  3.6738E+00  1.8047E+00  5.7623E-01  1.0518E+00
             9.5925E-01
 PARAMETER:  1.0445E-01 -3.8055E+00  8.1974E-01  5.7463E-01  7.3989E-02  7.9202E-02  1.4012E+00  6.9038E-01 -4.5124E-01  1.5048E-01
             5.8399E-02
 GRADIENT:   5.0831E+02  1.4577E+00  8.9267E+00  1.2441E+03  1.2363E+01  4.7705E+01  2.6139E-01  3.8945E+00  3.0669E+01  9.7522E-01
             1.2744E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1719.40284832159        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1802             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0045E+00  1.0000E-02  2.0532E+00  1.6141E+00  9.7123E-01  9.7911E-01  3.7992E+00  1.8027E+00  5.7436E-01  1.0518E+00
             9.5866E-01
 PARAMETER:  1.0449E-01 -4.9948E+00  8.1940E-01  5.7876E-01  7.0805E-02  7.8889E-02  1.4348E+00  6.8928E-01 -4.5450E-01  1.5046E-01
             5.7785E-02
 GRADIENT:   5.0935E+02  0.0000E+00  9.1441E+00  1.2669E+03  1.1478E+01  4.7628E+01  7.5609E-02  3.8188E+00  3.1045E+01  1.2686E+00
             1.0548E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1719.42244661931        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     1947
 NPARAMETR:  1.0041E+00  1.0000E-02  2.0523E+00  1.6164E+00  9.6680E-01  9.7884E-01  5.8835E+00  1.7970E+00  5.7332E-01  1.0515E+00
             9.5815E-01
 PARAMETER:  1.0408E-01 -4.9948E+00  8.1896E-01  5.8019E-01  6.6241E-02  7.8608E-02  1.8722E+00  6.8614E-01 -4.5632E-01  1.5018E-01
             5.7245E-02
 GRADIENT:   5.0681E+02  0.0000E+00  1.0892E+01  1.2816E+03  6.3639E+00  4.7610E+01  2.1542E-01  3.5915E+00  3.0860E+01  1.8022E+00
             7.2794E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1719.42769854333        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2137             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0042E+00  1.0000E-02  2.0420E+00  1.6162E+00  9.6441E-01  9.7898E-01  6.3523E+00  1.7891E+00  5.7327E-01  1.0498E+00
             9.5833E-01
 PARAMETER:  1.0424E-01 -4.9948E+00  8.1394E-01  5.8010E-01  6.3766E-02  7.8761E-02  1.9488E+00  6.8170E-01 -4.5639E-01  1.4856E-01
             5.7437E-02
 GRADIENT:   5.0778E+02  0.0000E+00  1.0965E+01  1.2814E+03  5.9731E+00  4.7612E+01  2.5669E-01  3.5453E+00  3.0822E+01  1.8028E+00
             7.9609E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1719.43240669683        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2328             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0042E+00  1.0000E-02  2.0340E+00  1.6159E+00  9.6321E-01  9.7901E-01  6.8995E+00  1.7834E+00  5.7328E-01  1.0488E+00
             9.5845E-01
 PARAMETER:  1.0424E-01 -4.9948E+00  8.1002E-01  5.7989E-01  6.2513E-02  7.8787E-02  2.0314E+00  6.7855E-01 -4.5637E-01  1.4769E-01
             5.7563E-02
 GRADIENT:   5.0747E+02  0.0000E+00  1.0745E+01  1.2800E+03  6.3883E+00  4.7606E+01  3.1052E-01  3.5391E+00  3.0817E+01  1.7763E+00
             8.5905E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1719.43613197407        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2521
 NPARAMETR:  1.0042E+00  1.0000E-02  2.0273E+00  1.6157E+00  9.6237E-01  9.7899E-01  7.3360E+00  1.7801E+00  5.7323E-01  1.0475E+00
             9.5835E-01
 PARAMETER:  1.0421E-01 -4.9948E+00  8.0670E-01  5.7978E-01  6.1640E-02  7.8771E-02  2.0928E+00  6.7668E-01 -4.5647E-01  1.4639E-01
             5.7463E-02
 GRADIENT:   1.9635E+00  0.0000E+00  5.7652E-01 -2.3716E+01 -1.1003E-01  1.0678E-01 -5.6524E-03  2.2607E-01  2.3434E-01  1.7296E-01
             1.8346E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1719.44059695821        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2719             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0042E+00  1.0000E-02  2.0213E+00  1.6155E+00  9.6077E-01  9.7900E-01  8.6069E+00  1.7746E+00  5.7319E-01  1.0468E+00
             9.5837E-01
 PARAMETER:  1.0420E-01 -4.9948E+00  8.0375E-01  5.7961E-01  5.9984E-02  7.8778E-02  2.2526E+00  6.7355E-01 -4.5653E-01  1.4574E-01
             5.7474E-02
 GRADIENT:   5.0726E+02  0.0000E+00  1.0597E+01  1.2792E+03  6.5660E+00  4.7576E+01  5.2913E-01  3.5369E+00  3.0890E+01  1.7290E+00
             8.3253E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1719.44702128155        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     2789
 NPARAMETR:  1.0045E+00  1.0000E-02  1.9849E+00  1.6153E+00  9.5470E-01  9.7959E-01  6.9194E+00  1.7504E+00  5.7339E-01  1.0406E+00
             9.5831E-01
 PARAMETER:  1.0446E-01 -4.9948E+00  7.8558E-01  5.7950E-01  5.3642E-02  7.9376E-02  2.0343E+00  6.5987E-01 -4.5618E-01  1.3980E-01
             5.7420E-02
 GRADIENT:   5.0865E+02  0.0000E+00  9.4687E+00  1.2813E+03  8.0498E+00  4.7696E+01  3.1060E-01  3.5737E+00  3.0618E+01  1.3613E+00
             7.5683E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1719.45231460841        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     2958             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0041E+00  1.0000E-02  1.9656E+00  1.6132E+00  9.5144E-01  9.7903E-01  6.6438E+00  1.7343E+00  5.7399E-01  1.0368E+00
             9.5833E-01
 PARAMETER:  1.0409E-01 -4.9948E+00  7.7578E-01  5.7825E-01  5.0226E-02  7.8805E-02  1.9937E+00  6.5063E-01 -4.5514E-01  1.3610E-01
             5.7432E-02
 GRADIENT:   5.0583E+02  0.0000E+00  9.1990E+00  1.2725E+03  9.4032E+00  4.7475E+01  2.8365E-01  3.3817E+00  3.0674E+01  1.0517E+00
             7.8397E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1719.45525649810        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     3151
 NPARAMETR:  1.0041E+00  1.0000E-02  1.9645E+00  1.6131E+00  9.5019E-01  9.7903E-01  7.0561E+00  1.7313E+00  5.7397E-01  1.0368E+00
             9.5829E-01
 PARAMETER:  1.0409E-01 -4.9948E+00  7.7523E-01  5.7818E-01  4.8902E-02  7.8807E-02  2.0539E+00  6.4890E-01 -4.5518E-01  1.3612E-01
             5.7392E-02
 GRADIENT:   1.9986E+00  0.0000E+00 -1.0778E-01 -2.4013E+01  1.0067E+00  1.0195E-01 -6.5221E-03  1.1311E-01  2.2301E-01 -1.6857E-01
            -3.3981E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1719.45777764087        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3348             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0041E+00  1.0000E-02  1.9616E+00  1.6130E+00  9.4871E-01  9.7903E-01  7.8322E+00  1.7283E+00  5.7389E-01  1.0375E+00
             9.5833E-01
 PARAMETER:  1.0408E-01 -4.9948E+00  7.7377E-01  5.7808E-01  4.7353E-02  7.8812E-02  2.1582E+00  6.4712E-01 -4.5532E-01  1.3681E-01
             5.7442E-02
 GRADIENT:   5.0572E+02  0.0000E+00  1.0194E+01  1.2722E+03  7.0418E+00  4.7500E+01  4.1711E-01  3.2669E+00  3.0728E+01  1.4924E+00
             8.1326E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1719.46029642698        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     3542             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0041E+00  1.0000E-02  1.9583E+00  1.6129E+00  9.4887E-01  9.7905E-01  9.4985E+00  1.7264E+00  5.7329E-01  1.0360E+00
             9.5825E-01
 PARAMETER:  1.0407E-01 -4.9948E+00  7.7209E-01  5.7802E-01  4.7513E-02  7.8823E-02  2.3511E+00  6.4604E-01 -4.5636E-01  1.3539E-01
             5.7354E-02
 GRADIENT:   5.0586E+02  0.0000E+00  9.7663E+00  1.2721E+03  8.1828E+00  4.7490E+01  6.9648E-01  3.2410E+00  3.0611E+01  1.2358E+00
             7.4692E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1719.46064031047        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     3733
 NPARAMETR:  1.0041E+00  1.0000E-02  1.9550E+00  1.6128E+00  9.4835E-01  9.7904E-01  9.5517E+00  1.7245E+00  5.7347E-01  1.0350E+00
             9.5830E-01
 PARAMETER:  1.0407E-01 -4.9948E+00  7.7038E-01  5.7794E-01  4.6964E-02  7.8818E-02  2.3567E+00  6.4491E-01 -4.5605E-01  1.3442E-01
             5.7407E-02
 GRADIENT:   1.9183E+00  0.0000E+00 -2.6776E-01 -2.3981E+01  1.2373E+00  1.0654E-01  1.5617E-02  1.2759E-01  1.8171E-01 -2.1266E-01
            -3.1756E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1719.46183399513        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     3927             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0041E+00  1.0000E-02  1.9526E+00  1.6126E+00  9.4677E-01  9.7904E-01  9.5603E+00  1.7208E+00  5.7354E-01  1.0360E+00
             9.5833E-01
 PARAMETER:  1.0406E-01 -4.9948E+00  7.6918E-01  5.7784E-01  4.5298E-02  7.8821E-02  2.3576E+00  6.4278E-01 -4.5592E-01  1.3541E-01
             5.7438E-02
 GRADIENT:   5.0557E+02  0.0000E+00  1.0202E+01  1.2713E+03  7.0257E+00  4.7494E+01  7.1053E-01  3.1967E+00  3.0688E+01  1.4724E+00
             8.0605E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1719.46224550624        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     4118
 NPARAMETR:  1.0041E+00  1.0000E-02  1.9501E+00  1.6125E+00  9.4607E-01  9.7905E-01  9.5639E+00  1.7186E+00  5.7358E-01  1.0356E+00
             9.5834E-01
 PARAMETER:  1.0405E-01 -4.9948E+00  7.6788E-01  5.7777E-01  4.4564E-02  7.8823E-02  2.3580E+00  6.4148E-01 -4.5586E-01  1.3503E-01
             5.7448E-02
 GRADIENT:   1.8488E+00  0.0000E+00  3.6136E-01 -2.3845E+01 -3.9898E-01  1.2323E-01  1.6193E-02  7.1158E-02  2.1051E-01  1.3001E-01
             9.8623E-03

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1719.46278603402        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     4312             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0041E+00  1.0000E-02  1.9447E+00  1.6123E+00  9.4591E-01  9.7905E-01  9.5586E+00  1.7161E+00  5.7360E-01  1.0339E+00
             9.5832E-01
 PARAMETER:  1.0405E-01 -4.9948E+00  7.6511E-01  5.7766E-01  4.4391E-02  7.8824E-02  2.3574E+00  6.4005E-01 -4.5582E-01  1.3335E-01
             5.7426E-02
 GRADIENT:   5.0528E+02  0.0000E+00  9.6994E+00  1.2704E+03  8.1217E+00  4.7449E+01  7.1026E-01  3.2137E+00  3.0661E+01  1.2190E+00
             7.8580E-01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1719.46292036070        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     4503
 NPARAMETR:  1.0041E+00  1.0000E-02  1.9422E+00  1.6122E+00  9.4555E-01  9.7905E-01  9.5597E+00  1.7144E+00  5.7362E-01  1.0333E+00
             9.5831E-01
 PARAMETER:  1.0404E-01 -4.9948E+00  7.6382E-01  5.7760E-01  4.4013E-02  7.8825E-02  2.3576E+00  6.3904E-01 -4.5578E-01  1.3273E-01
             5.7421E-02
 GRADIENT:   1.9965E+00  0.0000E+00 -3.1319E-01 -2.4042E+01  1.1280E+00  1.2859E-01  1.5679E-02  1.2518E-01  1.8272E-01 -1.8004E-01
            -2.1304E-02

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1719.46345061307        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     4697             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0040E+00  1.0000E-02  1.9416E+00  1.6121E+00  9.4404E-01  9.7906E-01  9.5737E+00  1.7114E+00  5.7369E-01  1.0345E+00
             9.5836E-01
 PARAMETER:  1.0404E-01 -4.9948E+00  7.6352E-01  5.7753E-01  4.2413E-02  7.8833E-02  2.3590E+00  6.3733E-01 -4.5567E-01  1.3390E-01
             5.7467E-02
 GRADIENT:   5.0527E+02  0.0000E+00  1.0360E+01  1.2698E+03  6.5838E+00  4.7447E+01  7.1351E-01  3.1190E+00  3.0697E+01  1.5148E+00
             8.1892E-01

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1719.46387631100        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     4888
 NPARAMETR:  1.0040E+00  1.0000E-02  1.9395E+00  1.6120E+00  9.4378E-01  9.7905E-01  9.5714E+00  1.7103E+00  5.7370E-01  1.0340E+00
             9.5835E-01
 PARAMETER:  1.0403E-01 -4.9948E+00  7.6241E-01  5.7748E-01  4.2142E-02  7.8831E-02  2.3588E+00  6.3666E-01 -4.5565E-01  1.3345E-01
             5.7456E-02
 GRADIENT:   1.9706E+00  0.0000E+00  3.0298E-01 -2.3920E+01 -4.2107E-01  1.2843E-01  1.6271E-02  7.0715E-02  2.1298E-01  1.2935E-01
             1.3205E-02

0ITERATION NO.:  140    OBJECTIVE VALUE:  -1719.46399481513        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     5082             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0040E+00  1.0000E-02  1.9347E+00  1.6119E+00  9.4387E-01  9.7905E-01  9.5653E+00  1.7082E+00  5.7371E-01  1.0321E+00
             9.5831E-01
 PARAMETER:  1.0403E-01 -4.9948E+00  7.5995E-01  5.7739E-01  4.2237E-02  7.8832E-02  2.3581E+00  6.3546E-01 -4.5563E-01  1.3157E-01
             5.7420E-02
 GRADIENT:   5.0506E+02  0.0000E+00  9.5964E+00  1.2692E+03  8.3164E+00  4.7435E+01  7.1176E-01  3.1592E+00  3.0646E+01  1.1362E+00
             7.7604E-01

0ITERATION NO.:  144    OBJECTIVE VALUE:  -1719.46435557991        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     5224
 NPARAMETR:  1.0040E+00  1.0000E-02  1.9331E+00  1.6118E+00  9.4344E-01  9.7905E-01  9.5673E+00  1.7067E+00  5.7374E-01  1.0321E+00
             9.5831E-01
 PARAMETER:  1.0403E-01 -4.9948E+00  7.6055E-01  5.7738E-01  4.1352E-02  7.8833E-02  2.3590E+00  6.3501E-01 -4.5558E-01  1.3288E-01
             5.7458E-02
 GRADIENT:   3.4234E-03  0.0000E+00  3.2042E-01  1.9238E-01 -3.7474E-01  3.3354E-04  1.4125E-04  4.1091E-02  2.5771E-03  6.9955E-02
             1.3719E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     5224
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.5748E-04 -4.6172E-05 -3.6799E-02 -9.5491E-03 -4.4489E-02
 SE:             2.9868E-02  1.9737E-03  1.8793E-02  2.8834E-02  2.0402E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9579E-01  9.8134E-01  5.0213E-02  7.4051E-01  2.9214E-02

 ETASHRINKSD(%)  1.0000E-10  9.3388E+01  3.7042E+01  3.4034E+00  3.1650E+01
 ETASHRINKVR(%)  1.0000E-10  9.9563E+01  6.0363E+01  6.6909E+00  5.3282E+01
 EBVSHRINKSD(%)  3.8582E-01  9.3562E+01  4.0713E+01  3.8215E+00  2.7570E+01
 EBVSHRINKVR(%)  7.7015E-01  9.9586E+01  6.4850E+01  7.4969E+00  4.7539E+01
 RELATIVEINF(%)  9.5661E+01  8.8720E-03  8.9385E+00  2.1815E+00  9.1105E+00
 EPSSHRINKSD(%)  4.5190E+01
 EPSSHRINKVR(%)  6.9959E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1719.4643555799114     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -984.31352901617322     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    73.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.27
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1719.464       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  1.94E+00  1.61E+00  9.43E-01  9.79E-01  9.57E+00  1.71E+00  5.74E-01  1.03E+00  9.58E-01
 


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
+        1.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.24E+00  0.00E+00  3.87E+01
 
 TH 4
+       -1.69E+01  0.00E+00 -2.62E+01  1.20E+03
 
 TH 5
+       -1.94E+00  0.00E+00 -1.14E+02 -9.77E+01  6.43E+02
 
 TH 6
+       -6.36E-02  0.00E+00  8.19E-02 -3.23E+00  2.44E-01  2.04E+02
 
 TH 7
+        9.28E-04  0.00E+00  2.23E-03  1.80E-03 -2.60E-03  6.15E-04  1.61E-03
 
 TH 8
+        1.07E-01  0.00E+00 -1.45E+01 -4.63E+00 -9.94E+00  2.27E-02  9.36E-05  1.98E+01
 
 TH 9
+        2.75E+00  0.00E+00  4.92E+00 -1.64E+00  4.65E-01 -2.74E+00  1.16E-01  7.86E-01  5.23E+02
 
 TH10
+        8.43E-01  0.00E+00 -8.22E-03 -3.51E+00 -7.65E+01  4.42E-01  4.06E-03  8.77E+00  2.87E+00  6.45E+01
 
 TH11
+       -8.26E+00  0.00E+00 -3.21E+00 -1.65E+01 -4.13E+00  3.25E+00  5.17E-03  4.56E+00  2.01E+01  1.02E+01  2.26E+02
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       80.038
Stop Time:
Wed Sep 29 19:14:38 CDT 2021
