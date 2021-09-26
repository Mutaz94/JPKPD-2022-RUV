Fri Sep 24 21:04:40 CDT 2021
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
$DATA ../../../../data/int/A2/dat1.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2638.21386521882        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.6065E+01  8.8097E+01  1.7666E+02 -1.1370E+01  7.4621E+01  1.5990E+01 -1.4041E+02 -1.7768E+02 -1.5017E+01 -1.0032E+02
            -1.9992E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3177.29972074635        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9329E-01  9.9810E-01  8.8666E-01  1.0288E+00  9.6249E-01  9.3195E-01  1.0683E+00  1.0202E+00  9.6226E-01  9.5871E-01
             1.8101E+00
 PARAMETER:  9.3269E-02  9.8094E-02 -2.0296E-02  1.2844E-01  6.1765E-02  2.9524E-02  1.6606E-01  1.2002E-01  6.1527E-02  5.7836E-02
             6.9336E-01
 GRADIENT:  -5.5293E+00  1.2757E+01 -2.8900E+00  9.5105E+00  2.6647E+00 -8.7433E+00 -2.2061E+00  3.0060E+00 -3.4224E+00 -1.8092E+00
            -9.6507E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3178.33178494518        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  9.9935E-01  8.4884E-01  7.5904E-01  1.1083E+00  8.1049E-01  9.4859E-01  1.1983E+00  7.9579E-01  9.4454E-01  8.5684E-01
             1.8111E+00
 PARAMETER:  9.9350E-02 -6.3883E-02 -1.7570E-01  2.0280E-01 -1.1011E-01  4.7223E-02  2.8094E-01 -1.2841E-01  4.2946E-02 -5.4506E-02
             6.9393E-01
 GRADIENT:   1.0323E+01  1.4376E+01 -1.0811E+01  3.3212E+01  9.5416E+00 -1.8799E+00  1.8075E-01  3.1153E+00 -3.6012E+00 -1.6512E+00
             2.8875E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3179.14412642102        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      238
 NPARAMETR:  9.9223E-01  7.0127E-01  6.4969E-01  1.1687E+00  6.6720E-01  9.5258E-01  1.3106E+00  5.5136E-01  9.5033E-01  7.9179E-01
             1.7990E+00
 PARAMETER:  9.2203E-02 -2.5487E-01 -3.3126E-01  2.5588E-01 -3.0466E-01  5.1417E-02  3.7045E-01 -4.9537E-01  4.9051E-02 -1.3346E-01
             6.8721E-01
 GRADIENT:  -5.3764E+00  1.1951E+01  2.2901E+01  4.0738E+01 -8.8555E+00 -5.4205E-01 -2.3938E-01 -1.3052E+00  1.6613E+00  6.5277E-01
             3.9982E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3179.14701312182        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      317
 NPARAMETR:  9.9133E-01  6.7946E-01  6.3210E-01  1.1784E+00  6.4606E-01  9.5296E-01  1.3277E+00  5.2618E-01  9.5020E-01  7.8158E-01
             1.7966E+00
 PARAMETER:  9.1292E-02 -2.8646E-01 -3.5871E-01  2.6418E-01 -3.3686E-01  5.1817E-02  3.8347E-01 -5.4211E-01  4.8913E-02 -1.4644E-01
             6.8589E-01
 GRADIENT:  -7.3774E+00  1.1897E+01  2.9343E+01  4.5746E+01 -1.2890E+01 -4.8541E-01 -1.7002E-01 -1.6311E+00  1.9002E+00  7.3337E-01
             8.9141E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3179.14783847962        NO. OF FUNC. EVALS.:  82
 CUMULATIVE NO. OF FUNC. EVALS.:      399
 NPARAMETR:  9.9089E-01  6.6924E-01  6.2363E-01  1.1828E+00  6.3618E-01  9.5317E-01  1.3356E+00  5.1418E-01  9.5021E-01  7.7692E-01
             1.7955E+00
 PARAMETER:  9.0852E-02 -3.0161E-01 -3.7220E-01  2.6785E-01 -3.5227E-01  5.2033E-02  3.8939E-01 -5.6518E-01  4.8931E-02 -1.5241E-01
             6.8526E-01
 GRADIENT:  -8.3405E+00  1.1773E+01  3.2378E+01  4.7915E+01 -1.4834E+01 -4.5127E-01 -1.3353E-01 -1.7816E+00  2.0094E+00  7.7103E-01
            -7.0245E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3179.14814120884        NO. OF FUNC. EVALS.:  84
 CUMULATIVE NO. OF FUNC. EVALS.:      483
 NPARAMETR:  9.9062E-01  6.6292E-01  6.1831E-01  1.1854E+00  6.3007E-01  9.5330E-01  1.3405E+00  5.0665E-01  9.5024E-01  7.7408E-01
             1.7948E+00
 PARAMETER:  9.0577E-02 -3.1111E-01 -3.8076E-01  2.7005E-01 -3.6193E-01  5.2175E-02  3.9301E-01 -5.7994E-01  4.8963E-02 -1.5608E-01
             6.8487E-01
 GRADIENT:  -8.9425E+00  1.1669E+01  3.4261E+01  4.9209E+01 -1.6056E+01 -4.2830E-01 -1.0917E-01 -1.8740E+00  2.0751E+00  7.9333E-01
            -1.7321E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3179.14840918803        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      568
 NPARAMETR:  9.9042E-01  6.5829E-01  6.1439E-01  1.1872E+00  6.2560E-01  9.5340E-01  1.3440E+00  5.0106E-01  9.5028E-01  7.7202E-01
             1.7943E+00
 PARAMETER:  9.0376E-02 -3.1811E-01 -3.8713E-01  2.7163E-01 -3.6905E-01  5.2281E-02  3.9562E-01 -5.9102E-01  4.8997E-02 -1.5875E-01
             6.8459E-01
 GRADIENT:  -9.3845E+00  1.1579E+01  3.5638E+01  5.0132E+01 -1.6955E+01 -4.1091E-01 -9.0691E-02 -1.9411E+00  2.1222E+00  8.0918E-01
            -2.5023E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3179.14856870877        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      653
 NPARAMETR:  9.9027E-01  6.5471E-01  6.1133E-01  1.1887E+00  6.2214E-01  9.5348E-01  1.3467E+00  4.9670E-01  9.5031E-01  7.7044E-01
             1.7939E+00
 PARAMETER:  9.0219E-02 -3.2356E-01 -3.9211E-01  2.7283E-01 -3.7458E-01  5.2365E-02  3.9763E-01 -5.9976E-01  4.9028E-02 -1.6080E-01
             6.8438E-01
 GRADIENT:  -9.7274E+00  1.1502E+01  3.6702E+01  5.0829E+01 -1.7653E+01 -3.9695E-01 -7.6048E-02 -1.9927E+00  2.1582E+00  8.2117E-01
            -3.1074E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3179.14867381326        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      738
 NPARAMETR:  9.9014E-01  6.5180E-01  6.0883E-01  1.1898E+00  6.1933E-01  9.5355E-01  1.3489E+00  4.9311E-01  9.5033E-01  7.6916E-01
             1.7936E+00
 PARAMETER:  9.0091E-02 -3.2802E-01 -3.9621E-01  2.7379E-01 -3.7911E-01  5.2435E-02  3.9925E-01 -6.0702E-01  4.9058E-02 -1.6246E-01
             6.8420E-01
 GRADIENT:  -1.0007E+01  1.1435E+01  3.7568E+01  5.1388E+01 -1.8224E+01 -3.8528E-01 -6.3909E-02 -2.0346E+00  2.1871E+00  8.3080E-01
            -3.6060E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3179.14874863531        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      823
 NPARAMETR:  9.9003E-01  6.4933E-01  6.0671E-01  1.1908E+00  6.1695E-01  9.5361E-01  1.3507E+00  4.9004E-01  9.5036E-01  7.6808E-01
             1.7933E+00
 PARAMETER:  8.9983E-02 -3.3181E-01 -3.9970E-01  2.7460E-01 -3.8296E-01  5.2495E-02  4.0062E-01 -6.1326E-01  4.9085E-02 -1.6386E-01
             6.8405E-01
 GRADIENT:  -1.0244E+01  1.1374E+01  3.8300E+01  5.1854E+01 -1.8708E+01 -3.7521E-01 -5.3493E-02 -2.0699E+00  2.2113E+00  8.3884E-01
            -4.0320E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3179.14880385647        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      908
 NPARAMETR:  9.8994E-01  6.4720E-01  6.0487E-01  1.1916E+00  6.1490E-01  9.5366E-01  1.3523E+00  4.8737E-01  9.5038E-01  7.6716E-01
             1.7931E+00
 PARAMETER:  8.9889E-02 -3.3510E-01 -4.0274E-01  2.7529E-01 -3.8630E-01  5.2547E-02  4.0179E-01 -6.1873E-01  4.9110E-02 -1.6506E-01
             6.8393E-01
 GRADIENT:  -1.0450E+01  1.1320E+01  3.8932E+01  5.2252E+01 -1.9127E+01 -3.6636E-01 -4.4407E-02 -2.1004E+00  2.2320E+00  8.4569E-01
            -4.4025E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3179.14884715954        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      993
 NPARAMETR:  9.8986E-01  6.4530E-01  6.0322E-01  1.1923E+00  6.1307E-01  9.5370E-01  1.3537E+00  4.8497E-01  9.5041E-01  7.6634E-01
             1.7929E+00
 PARAMETER:  8.9806E-02 -3.3803E-01 -4.0547E-01  2.7590E-01 -3.8928E-01  5.2595E-02  4.0283E-01 -6.2367E-01  4.9134E-02 -1.6613E-01
             6.8381E-01
 GRADIENT:  -1.0632E+01  1.1269E+01  3.9494E+01  5.2602E+01 -1.9501E+01 -3.5835E-01 -3.6222E-02 -2.1274E+00  2.2504E+00  8.5174E-01
            -4.7348E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3179.14888141728        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:     1078
 NPARAMETR:  9.8978E-01  6.4361E-01  6.0174E-01  1.1930E+00  6.1143E-01  9.5374E-01  1.3549E+00  4.8281E-01  9.5043E-01  7.6561E-01
             1.7927E+00
 PARAMETER:  8.9731E-02 -3.4067E-01 -4.0792E-01  2.7644E-01 -3.9196E-01  5.2637E-02  4.0376E-01 -6.2814E-01  4.9156E-02 -1.6709E-01
             6.8371E-01
 GRADIENT:  -1.0796E+01  1.1222E+01  3.9996E+01  5.2912E+01 -1.9835E+01 -3.5110E-01 -2.8855E-02 -2.1515E+00  2.2667E+00  8.5710E-01
            -5.0331E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3179.14890725277        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  9.8972E-01  6.4211E-01  6.0044E-01  1.1935E+00  6.0999E-01  9.5378E-01  1.3561E+00  4.8088E-01  9.5045E-01  7.6497E-01
             1.7925E+00
 PARAMETER:  8.9665E-02 -3.4299E-01 -4.1010E-01  2.7691E-01 -3.9432E-01  5.2675E-02  4.0458E-01 -6.3213E-01  4.9177E-02 -1.6792E-01
             6.8363E-01
 GRADIENT:  -1.0940E+01  1.1180E+01  4.0438E+01  5.3183E+01 -2.0130E+01 -3.4465E-01 -2.2325E-02 -2.1726E+00  2.2809E+00  8.6178E-01
            -5.2970E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -3179.14892882983        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:     1248
 NPARAMETR:  9.8966E-01  6.4070E-01  5.9921E-01  1.1941E+00  6.0863E-01  9.5381E-01  1.3571E+00  4.7906E-01  9.5047E-01  7.6437E-01
             1.7924E+00
 PARAMETER:  8.9603E-02 -3.4519E-01 -4.1215E-01  2.7735E-01 -3.9655E-01  5.2711E-02  4.0534E-01 -6.3592E-01  4.9197E-02 -1.6871E-01
             6.8354E-01
 GRADIENT:  -1.1075E+01  1.1140E+01  4.0853E+01  5.3435E+01 -2.0407E+01 -3.3858E-01 -1.6152E-02 -2.1925E+00  2.2943E+00  8.6615E-01
            -5.5460E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -3179.14893967027        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:     1333
 NPARAMETR:  9.8960E-01  6.3933E-01  5.9800E-01  1.1946E+00  6.0730E-01  9.5384E-01  1.3581E+00  4.7728E-01  9.5048E-01  7.6378E-01
             1.7922E+00
 PARAMETER:  8.9543E-02 -3.4734E-01 -4.1416E-01  2.7778E-01 -3.9873E-01  5.2745E-02  4.0609E-01 -6.3966E-01  4.9217E-02 -1.6947E-01
             6.8346E-01
 GRADIENT:  -1.1207E+01  1.1099E+01  4.1257E+01  5.3680E+01 -2.0678E+01 -3.3256E-01 -1.0077E-02 -2.2118E+00  2.3072E+00  8.7037E-01
            -5.7898E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -3179.28448553712        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1515
 NPARAMETR:  9.9461E-01  6.4315E-01  6.1187E-01  1.1856E+00  6.1922E-01  9.5777E-01  1.3690E+00  5.1851E-01  9.4692E-01  7.6015E-01
             1.7918E+00
 PARAMETER:  9.4595E-02 -3.4138E-01 -3.9123E-01  2.7022E-01 -3.7930E-01  5.6853E-02  4.1408E-01 -5.5680E-01  4.5456E-02 -1.7424E-01
             6.8322E-01
 GRADIENT:  -1.1712E+01 -4.0523E+00  4.2277E+01  3.0231E+01 -2.4445E+01  4.4187E-01  6.3184E-01 -1.8848E+00  1.2438E+00  1.1179E-01
            -1.2686E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -3187.69358232291        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1698
 NPARAMETR:  9.7084E-01  3.4284E-01  2.8714E-01  1.2512E+00  3.0287E-01  9.8557E-01  1.6009E+00  1.2205E+00  1.0504E+00  5.6083E-01
             1.6232E+00
 PARAMETER:  7.0410E-02 -9.7049E-01 -1.1478E+00  3.2408E-01 -1.0945E+00  8.5460E-02  5.7055E-01  2.9928E-01  1.4915E-01 -4.7834E-01
             5.8437E-01
 GRADIENT:  -6.6997E+01  3.1788E+01  1.1106E+02  1.6216E+02 -1.5075E+02  1.8325E+00  7.4164E+00  1.7792E+01 -1.3321E+01  2.4376E+00
             4.6530E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -3219.26440876730        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1877
 NPARAMETR:  1.0020E+00  2.3928E-01  1.6535E-01  1.0405E+00  2.2217E-01  9.7710E-01  1.5311E+00  1.1014E+00  1.1858E+00  5.5145E-01
             1.5250E+00
 PARAMETER:  1.0199E-01 -1.3301E+00 -1.6997E+00  1.3967E-01 -1.4043E+00  7.6830E-02  5.2597E-01  1.9654E-01  2.7043E-01 -4.9520E-01
             5.2202E-01
 GRADIENT:   6.7390E+00 -6.8993E+00  7.2979E+00  2.1584E+00  7.9408E+00  4.5408E-02 -5.5976E+00  1.6870E+00 -3.0851E+00 -5.9583E-01
            -2.1050E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -3220.01101565688        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2052
 NPARAMETR:  9.9914E-01  2.2381E-01  1.4651E-01  1.0038E+00  2.0568E-01  9.7622E-01  1.5454E+00  1.0935E+00  1.2375E+00  5.9031E-01
             1.5358E+00
 PARAMETER:  9.9141E-02 -1.3970E+00 -1.8206E+00  1.0382E-01 -1.4814E+00  7.5938E-02  5.3529E-01  1.8934E-01  3.1310E-01 -4.2711E-01
             5.2904E-01
 GRADIENT:   1.2695E-01  1.2797E-01 -2.6139E-01  3.7674E-01  3.9346E-01 -4.2768E-02  9.4557E-02  1.5014E-01 -5.5028E-02  1.1668E-01
             1.6692E-01

0ITERATION NO.:  103    OBJECTIVE VALUE:  -3220.01149237189        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     2144
 NPARAMETR:  9.9909E-01  2.2347E-01  1.4622E-01  1.0029E+00  2.0538E-01  9.7632E-01  1.5449E+00  1.0926E+00  1.2384E+00  5.9039E-01
             1.5357E+00
 PARAMETER:  9.9093E-02 -1.3985E+00 -1.8227E+00  1.0290E-01 -1.4829E+00  7.6038E-02  5.3498E-01  1.8851E-01  3.1386E-01 -4.2697E-01
             5.2897E-01
 GRADIENT:   2.4301E-02 -1.2564E-02 -3.5681E-02  7.3897E-02 -3.0092E-04  4.8777E-04  1.1294E-02  3.7969E-03 -3.2226E-02  2.5028E-03
             4.0606E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2144
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.9902E-04  1.0515E-02  1.1240E-02 -5.1568E-05  1.7301E-02
 SE:             2.9732E-02  2.6828E-02  2.2348E-02  2.8307E-02  2.0598E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8124E-01  6.9510E-01  6.1499E-01  9.9855E-01  4.0094E-01

 ETASHRINKSD(%)  3.9493E-01  1.0122E+01  2.5131E+01  5.1669E+00  3.0993E+01
 ETASHRINKVR(%)  7.8831E-01  1.9220E+01  4.3947E+01  1.0067E+01  5.2381E+01
 EBVSHRINKSD(%)  6.9901E-01  8.8226E+00  2.4147E+01  4.6028E+00  3.2065E+01
 EBVSHRINKVR(%)  1.3931E+00  1.6867E+01  4.2463E+01  8.9937E+00  5.3848E+01
 RELATIVEINF(%)  9.8595E+01  2.9376E+01  8.7422E+00  4.9420E+01  6.1041E+00
 EPSSHRINKSD(%)  2.3002E+01
 EPSSHRINKVR(%)  4.0713E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3220.0114923718861     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1565.9221326034753     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    45.40
 Elapsed covariance  time in seconds:    13.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3220.011       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  2.23E-01  1.46E-01  1.00E+00  2.05E-01  9.76E-01  1.54E+00  1.09E+00  1.24E+00  5.90E-01  1.54E+00
 


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
 
         2.94E-02  6.66E-02  7.35E-02  1.62E-01  6.29E-02  5.89E-02  1.52E-01  2.56E-01  2.23E-01  1.64E-01  1.02E-01
 


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
+        8.67E-04
 
 TH 2
+       -5.42E-04  4.44E-03
 
 TH 3
+       -6.09E-04  4.82E-03  5.40E-03
 
 TH 4
+       -1.49E-03  1.03E-02  1.15E-02  2.61E-02
 
 TH 5
+       -4.72E-04  4.13E-03  4.61E-03  9.83E-03  3.96E-03
 
 TH 6
+       -1.93E-04  1.28E-04  7.39E-05  6.49E-04  7.03E-05  3.47E-03
 
 TH 7
+       -1.58E-04  6.49E-03  7.01E-03  1.52E-02  6.11E-03  1.03E-04  2.32E-02
 
 TH 8
+        2.13E-03 -1.07E-02 -1.24E-02 -2.30E-02 -1.02E-02  2.65E-03 -1.30E-02  6.57E-02
 
 TH 9
+        1.71E-03 -1.31E-02 -1.44E-02 -3.09E-02 -1.23E-02  6.30E-04 -2.04E-02  4.11E-02  4.97E-02
 
 TH10
+        8.64E-04 -8.45E-03 -9.20E-03 -1.97E-02 -7.85E-03 -9.84E-04 -1.44E-02  2.21E-02  2.49E-02  2.69E-02
 
 TH11
+        3.79E-04  3.55E-03  3.95E-03  7.96E-03  3.42E-03 -5.28E-04  7.86E-03 -7.70E-03 -8.92E-03 -7.45E-03  1.03E-02
 
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
+        2.94E-02
 
 TH 2
+       -2.77E-01  6.66E-02
 
 TH 3
+       -2.81E-01  9.85E-01  7.35E-02
 
 TH 4
+       -3.13E-01  9.58E-01  9.73E-01  1.62E-01
 
 TH 5
+       -2.55E-01  9.85E-01  9.97E-01  9.67E-01  6.29E-02
 
 TH 6
+       -1.11E-01  3.27E-02  1.71E-02  6.82E-02  1.90E-02  5.89E-02
 
 TH 7
+       -3.53E-02  6.40E-01  6.27E-01  6.20E-01  6.39E-01  1.15E-02  1.52E-01
 
 TH 8
+        2.82E-01 -6.28E-01 -6.59E-01 -5.54E-01 -6.30E-01  1.75E-01 -3.33E-01  2.56E-01
 
 TH 9
+        2.60E-01 -8.82E-01 -8.80E-01 -8.58E-01 -8.74E-01  4.80E-02 -6.00E-01  7.19E-01  2.23E-01
 
 TH10
+        1.79E-01 -7.74E-01 -7.64E-01 -7.44E-01 -7.61E-01 -1.02E-01 -5.76E-01  5.25E-01  6.82E-01  1.64E-01
 
 TH11
+        1.27E-01  5.24E-01  5.29E-01  4.85E-01  5.35E-01 -8.81E-02  5.09E-01 -2.95E-01 -3.94E-01 -4.47E-01  1.02E-01
 
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
+        1.65E+03
 
 TH 2
+        4.71E+02  9.70E+03
 
 TH 3
+        1.29E+03 -8.78E+03  8.71E+04
 
 TH 4
+        2.35E+02  7.72E+02 -7.73E+03  1.49E+03
 
 TH 5
+       -2.12E+03 -7.43E+02 -7.11E+04  4.46E+03  7.10E+04
 
 TH 6
+        5.72E+01 -1.05E+02  4.91E+01 -4.97E+01  1.01E+02  3.31E+02
 
 TH 7
+       -1.62E+01 -4.93E+01  3.90E+02 -2.57E+01 -3.65E+02  7.42E-01  8.87E+01
 
 TH 8
+       -4.49E+01 -2.41E+02  1.86E+03 -2.40E+02 -1.30E+03 -1.92E+01 -3.57E+00  8.28E+01
 
 TH 9
+        5.36E+01  3.69E+02 -8.87E+02  1.68E+02  5.48E+02 -1.86E+01  2.02E+01 -5.75E+01  1.44E+02
 
 TH10
+        3.04E+01  1.97E+02 -7.32E+01  3.09E+01 -2.64E+01  2.79E+01  1.56E+01 -1.34E+01  9.83E+00  1.00E+02
 
 TH11
+       -1.45E+02 -2.57E+01 -5.86E+02  5.99E+01  3.48E+02  2.60E+01 -3.23E+01 -6.53E+00 -2.78E+01  4.89E+00  1.74E+02
 
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
 #CPUT: Total CPU Time in Seconds,       58.630
Stop Time:
Fri Sep 24 21:05:42 CDT 2021
