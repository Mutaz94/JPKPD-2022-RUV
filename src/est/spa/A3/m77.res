Wed Sep 29 13:47:40 CDT 2021
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
$DATA ../../../../data/spa/A3/dat77.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -596.659607398177        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7792E+02  5.8785E+01  1.4449E+02 -6.5455E+01  5.4216E+00  4.0993E+01  1.0967E+00 -5.5447E+01 -5.7126E+01 -5.1208E+01
            -1.9304E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1259.39791472558        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.9292E-01  9.6184E-01  7.8739E-01  1.0955E+00  8.5107E-01  9.1543E-01  8.6762E-01  1.0277E+00  1.0594E+00  9.7055E-01
             1.8747E+00
 PARAMETER:  9.2892E-02  6.1096E-02 -1.3903E-01  1.9123E-01 -6.1255E-02  1.1641E-02 -4.1997E-02  1.2730E-01  1.5771E-01  7.0103E-02
             7.2844E-01
 GRADIENT:   4.3524E+01  5.0812E+01  4.1041E+01  4.5989E+01 -3.0542E+01 -2.1700E+01  5.9861E+00 -3.8361E+00 -4.9962E+00  3.4289E+00
            -3.9415E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1285.31118158630        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      179
 NPARAMETR:  9.8633E-01  6.6638E-01  3.3939E-01  1.2429E+00  4.6199E-01  1.0191E+00  5.2072E-01  1.5728E-01  9.4958E-01  7.5753E-01
             1.9674E+00
 PARAMETER:  8.6233E-02 -3.0590E-01 -9.8060E-01  3.1741E-01 -6.7221E-01  1.1896E-01 -5.5254E-01 -1.7497E+00  4.8261E-02 -1.7769E-01
             7.7672E-01
 GRADIENT:  -7.9058E+01  1.1663E+01 -8.5311E+01  2.0504E+02  1.2712E+02  5.3886E+00 -2.7448E+00  6.0453E-02 -2.2552E+01 -5.9360E+00
            -3.1533E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1357.11398028483        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      357
 NPARAMETR:  1.0361E+00  8.6154E-01  4.8158E-01  1.0896E+00  6.3357E-01  9.8773E-01  4.1517E-01  4.4158E-01  1.0149E+00  4.8520E-01
             3.0866E+00
 PARAMETER:  1.3545E-01 -4.9028E-02 -6.3068E-01  1.8584E-01 -3.5638E-01  8.7658E-02 -7.7907E-01 -7.1739E-01  1.1476E-01 -6.2319E-01
             1.2271E+00
 GRADIENT:   1.0087E+01 -1.6821E+01 -4.6339E+00 -2.0597E+01  1.0172E+01  1.2008E+01  1.7655E-01  7.2821E-01  1.4391E+01  2.7271E+00
            -1.2034E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1361.18360738626        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  1.0372E+00  6.5513E-01  3.8807E-01  1.1828E+00  4.7979E-01  9.7886E-01  8.0699E-01  3.9929E-01  7.7283E-01  2.8989E-01
             3.1685E+00
 PARAMETER:  1.3654E-01 -3.2293E-01 -8.4657E-01  2.6791E-01 -6.3441E-01  7.8629E-02 -1.1444E-01 -8.1808E-01 -1.5770E-01 -1.1382E+00
             1.2533E+00
 GRADIENT:   4.5867E+00  3.6520E+00 -6.4522E+00  1.2788E+01  1.1646E+01  7.3429E+00 -4.3180E-01  8.6251E-01 -5.3334E+00  1.4408E+00
             1.5102E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1363.08298895852        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  1.0337E+00  4.6508E-01  3.6043E-01  1.2538E+00  4.0232E-01  9.5916E-01  1.1943E+00  3.8927E-01  7.5644E-01  1.7259E-01
             3.1126E+00
 PARAMETER:  1.3311E-01 -6.6555E-01 -9.2045E-01  3.2617E-01 -8.1050E-01  5.8301E-02  2.7756E-01 -8.4347E-01 -1.7913E-01 -1.6568E+00
             1.2354E+00
 GRADIENT:  -2.9795E+00  1.7990E+00 -3.8618E+00  1.6824E+01  7.2582E+00 -2.2666E+00 -5.1710E-01 -5.9816E-01 -1.7546E+00 -2.3447E-01
            -3.6978E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1364.16364656016        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  1.0342E+00  3.4628E-01  2.9959E-01  1.2455E+00  3.3072E-01  9.8382E-01  1.7835E+00  5.3421E-01  7.7197E-01  8.0089E-02
             2.9694E+00
 PARAMETER:  1.3364E-01 -9.6050E-01 -1.1053E+00  3.1952E-01 -1.0065E+00  8.3685E-02  6.7856E-01 -5.2697E-01 -1.5880E-01 -2.4246E+00
             1.1884E+00
 GRADIENT:  -2.4273E-01 -2.6348E+00 -5.9769E+00  8.6764E+00  8.8665E+00  2.6263E+00 -1.1635E-01 -2.7005E-01  2.1284E-01 -1.1024E-01
             3.6922E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1364.50545113329        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1062
 NPARAMETR:  1.0325E+00  4.3497E-01  2.6979E-01  1.1862E+00  3.2630E-01  9.7647E-01  1.4499E+00  5.2476E-01  7.9501E-01  1.2024E-01
             2.8370E+00
 PARAMETER:  1.3195E-01 -7.3247E-01 -1.2101E+00  2.7077E-01 -1.0199E+00  7.6191E-02  4.7146E-01 -5.4482E-01 -1.2940E-01 -2.0182E+00
             1.1428E+00
 GRADIENT:  -2.3898E+00  5.1212E+00  6.6623E+00  1.5811E+01 -1.0527E+01 -2.2088E+00 -6.8496E-01 -3.5127E-01 -3.0208E+00 -1.5173E-01
            -1.5677E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1365.17008502916        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1238
 NPARAMETR:  1.0317E+00  4.5793E-01  2.2722E-01  1.1231E+00  3.0302E-01  9.8618E-01  1.2734E+00  5.0435E-01  8.3608E-01  1.1219E-01
             2.8820E+00
 PARAMETER:  1.3125E-01 -6.8104E-01 -1.3818E+00  2.1607E-01 -1.0940E+00  8.6082E-02  3.4172E-01 -5.8449E-01 -7.9029E-02 -2.0876E+00
             1.1585E+00
 GRADIENT:   5.4593E-01  6.6981E-01  2.5096E+00 -1.1040E+00 -3.8820E+00  9.9515E-02  8.4878E-03 -1.4526E-01 -2.1398E-02 -1.1341E-01
             4.5856E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1365.19386966991        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1414            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0306E+00  4.7399E-01  2.2354E-01  1.1140E+00  3.0482E-01  9.8602E-01  1.2251E+00  5.0552E-01  8.4055E-01  1.6494E-01
             2.8781E+00
 PARAMETER:  1.3017E-01 -6.4658E-01 -1.3981E+00  2.0791E-01 -1.0880E+00  8.5922E-02  3.0303E-01 -5.8216E-01 -7.3702E-02 -1.7022E+00
             1.1571E+00
 GRADIENT:   4.9286E+01  7.1653E+00  1.5887E+01  1.7727E+01  5.4934E+01  3.6665E+00  1.6553E+00  6.4149E-01  8.8489E-01  2.8183E-01
             1.1182E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1365.23467750433        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1590
 NPARAMETR:  1.0308E+00  4.7239E-01  2.2426E-01  1.1155E+00  3.0465E-01  9.8577E-01  1.1893E+00  4.5753E-01  8.4147E-01  2.1455E-01
             2.8764E+00
 PARAMETER:  1.3037E-01 -6.4996E-01 -1.3950E+00  2.0928E-01 -1.0886E+00  8.5668E-02  2.7339E-01 -6.8191E-01 -7.2606E-02 -1.4392E+00
             1.1565E+00
 GRADIENT:   2.0865E-01  1.9497E-01 -9.9950E-01 -7.5789E-01  3.4378E+00  3.7927E-02  1.5133E-01  1.5511E-01  2.6094E-05  2.9153E-02
            -5.6297E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1365.33314678993        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1767
 NPARAMETR:  1.0360E+00  4.4697E-01  2.3048E-01  1.1341E+00  3.0171E-01  9.8515E-01  1.0918E+00  2.5929E-01  8.3710E-01  3.6348E-01
             2.9311E+00
 PARAMETER:  1.3532E-01 -7.0525E-01 -1.3676E+00  2.2584E-01 -1.0983E+00  8.5035E-02  1.8783E-01 -1.2498E+00 -7.7815E-02 -9.1203E-01
             1.1754E+00
 GRADIENT:   6.4405E+00  1.5385E+00 -1.9053E+00  5.9986E-01  5.7225E+00  3.2693E-01  3.8049E-01  1.2428E-01 -1.1360E-01  4.3268E-01
             4.2414E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1365.72147735995        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1942
 NPARAMETR:  1.0389E+00  3.6079E-01  2.5632E-01  1.1960E+00  2.9712E-01  9.8429E-01  1.0240E+00  6.6733E-02  8.2538E-01  4.9077E-01
             2.9320E+00
 PARAMETER:  1.3820E-01 -9.1946E-01 -1.2613E+00  2.7894E-01 -1.1136E+00  8.4166E-02  1.2373E-01 -2.6071E+00 -9.1914E-02 -6.1178E-01
             1.1757E+00
 GRADIENT:   5.5082E+00  6.7042E+00  1.0286E+01  4.5299E+00 -1.1474E+01  7.3066E-01  3.0138E-01  2.5009E-02  8.7280E-01  2.4937E-01
            -6.6329E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1366.14740468876        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2117
 NPARAMETR:  1.0374E+00  2.5756E-01  2.8381E-01  1.2601E+00  2.9885E-01  9.7692E-01  9.1289E-01  1.0000E-02  7.8979E-01  5.2276E-01
             2.9980E+00
 PARAMETER:  1.3674E-01 -1.2565E+00 -1.1595E+00  3.3121E-01 -1.1078E+00  7.6648E-02  8.8551E-03 -4.5227E+00 -1.3599E-01 -5.4863E-01
             1.1979E+00
 GRADIENT:   5.2997E+00  3.0772E+00  8.3325E+00  3.6885E+00 -8.4266E+00  1.3297E-01  2.4312E-02  3.9213E-04 -4.3846E-01 -1.4671E-02
             4.4196E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1366.38093015281        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2294
 NPARAMETR:  1.0288E+00  1.7205E-01  3.0904E-01  1.3141E+00  3.0775E-01  9.7230E-01  7.9520E-01  1.0000E-02  7.6799E-01  5.4439E-01
             2.9775E+00
 PARAMETER:  1.2839E-01 -1.6600E+00 -1.0743E+00  3.7313E-01 -1.0785E+00  7.1911E-02 -1.2916E-01 -6.5860E+00 -1.6398E-01 -5.0809E-01
             1.1911E+00
 GRADIENT:   4.5217E-01  1.2691E+00  5.5763E+00  3.2740E+00 -2.9788E+00  7.0061E-02  1.6283E-02  0.0000E+00  7.2157E-02 -1.1428E-01
            -1.0938E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1366.41263340362        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2469
 NPARAMETR:  1.0270E+00  1.3508E-01  3.1955E-01  1.3362E+00  3.1311E-01  9.7058E-01  7.2719E-01  1.0000E-02  7.5414E-01  5.3961E-01
             2.9998E+00
 PARAMETER:  1.2664E-01 -1.9019E+00 -1.0408E+00  3.8984E-01 -1.0612E+00  7.0137E-02 -2.1857E-01 -7.7278E+00 -1.8217E-01 -5.1691E-01
             1.1985E+00
 GRADIENT:   3.7558E+00  4.2005E-01  5.0801E-01  7.8446E-01  4.7181E+00  3.1532E-01  1.2673E-02  0.0000E+00  1.2549E-01 -4.5389E-02
             1.5365E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1367.05295113231        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2650
 NPARAMETR:  1.0236E+00  3.3641E-02  2.6961E-01  1.3003E+00  2.6992E-01  9.7907E-01  1.3888E+00  1.0000E-02  7.8455E-01  5.9836E-01
             2.8621E+00
 PARAMETER:  1.2335E-01 -3.2920E+00 -1.2108E+00  3.6258E-01 -1.2096E+00  7.8852E-02  4.2847E-01 -4.6404E+00 -1.4264E-01 -4.1357E-01
             1.1516E+00
 GRADIENT:   1.3976E+01 -2.2366E-01 -4.0098E+00 -1.0491E+01  1.4091E+01  2.0936E+00  3.4025E-03  0.0000E+00  2.4773E+00  3.1760E-01
            -7.0191E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1367.26593630432        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2829             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0176E+00  2.4150E-02  2.5711E-01  1.2882E+00  2.5951E-01  9.7383E-01  1.4054E+00  1.0000E-02  7.8318E-01  5.9137E-01
             2.8939E+00
 PARAMETER:  1.1749E-01 -3.6235E+00 -1.2583E+00  3.5324E-01 -1.2490E+00  7.3482E-02  4.4034E-01 -5.6784E+00 -1.4439E-01 -4.2531E-01
             1.1626E+00
 GRADIENT:   4.2595E+01  5.8247E-02  2.2137E+01  4.1899E+01  5.9227E+01  3.3022E+00  1.9329E-03  0.0000E+00  1.4850E+00  7.4332E-01
             9.6682E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1367.26873816099        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3013
 NPARAMETR:  1.0181E+00  3.3274E-02  2.5678E-01  1.2859E+00  2.5980E-01  9.7412E-01  6.2532E-01  1.0000E-02  7.8492E-01  5.8964E-01
             2.8900E+00
 PARAMETER:  1.1795E-01 -3.3030E+00 -1.2595E+00  3.5146E-01 -1.2479E+00  7.3776E-02 -3.6948E-01 -8.6126E+00 -1.4218E-01 -4.2824E-01
             1.1613E+00
 GRADIENT:  -4.1295E-01  1.6285E-02  8.2830E-01  8.5444E-01 -1.4453E+00 -1.0335E-02  2.2445E-04  0.0000E+00  4.1746E-02 -7.7671E-02
            -5.4000E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1367.27007834717        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     3189            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0185E+00  3.5760E-02  2.5656E-01  1.2844E+00  2.5982E-01  9.7416E-01  2.9503E-01  1.0000E-02  7.8491E-01  5.8922E-01
             2.8928E+00
 PARAMETER:  1.1838E-01 -3.2309E+00 -1.2604E+00  3.5032E-01 -1.2477E+00  7.3821E-02 -1.1207E+00 -9.8248E+00 -1.4219E-01 -4.2895E-01
             1.1622E+00
 GRADIENT:   4.2480E+01  1.1603E-01  2.2254E+01  4.1513E+01  5.9282E+01  3.2549E+00  3.6330E-04  0.0000E+00  1.2480E+00  6.3348E-01
             9.2927E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1367.27029859146        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3371
 NPARAMETR:  1.0186E+00  3.5839E-02  2.5643E-01  1.2843E+00  2.5983E-01  9.7421E-01  2.7023E-01  1.0000E-02  7.8515E-01  5.8908E-01
             2.8925E+00
 PARAMETER:  1.1842E-01 -3.2287E+00 -1.2609E+00  3.5024E-01 -1.2477E+00  7.3870E-02 -1.2085E+00 -9.8248E+00 -1.4188E-01 -4.2919E-01
             1.1621E+00
 GRADIENT:   9.2429E-02 -7.0243E-04  2.3414E-02 -1.3666E-01  4.6881E-02  4.1608E-03  5.0580E-05  0.0000E+00 -4.0808E-03  3.3730E-03
            -5.3583E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1367.27083147212        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3559
 NPARAMETR:  1.0188E+00  3.6824E-02  2.5613E-01  1.2832E+00  2.5967E-01  9.7428E-01  2.0890E-01  1.0000E-02  7.8550E-01  5.8913E-01
             2.8917E+00
 PARAMETER:  1.1866E-01 -3.2016E+00 -1.2621E+00  3.4940E-01 -1.2483E+00  7.3948E-02 -1.4659E+00 -9.8248E+00 -1.4143E-01 -4.2911E-01
             1.1618E+00
 GRADIENT:   4.6681E-01 -7.8992E-03  1.3722E-01 -8.7438E-01  1.9076E-01  1.3694E-02  3.2481E-05  0.0000E+00 -2.0776E-02 -8.5570E-03
            -1.4777E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1367.27104572321        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     3720
 NPARAMETR:  1.0186E+00  3.7011E-02  2.5606E-01  1.2832E+00  2.5960E-01  9.7421E-01  1.6452E-01  1.0000E-02  7.8567E-01  5.8898E-01
             2.8917E+00
 PARAMETER:  1.1838E-01 -3.1965E+00 -1.2624E+00  3.4933E-01 -1.2486E+00  7.3874E-02 -1.7047E+00 -9.8248E+00 -1.4121E-01 -4.2936E-01
             1.1618E+00
 GRADIENT:  -2.2238E-01  2.5312E-05  3.5850E-01 -6.2524E-01 -2.5333E-01 -1.9924E-02  2.0578E-05  0.0000E+00 -9.2516E-03 -2.8347E-02
            -1.2080E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1367.27129277130        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     3906             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0189E+00  3.7363E-02  2.5591E-01  1.2828E+00  2.5951E-01  9.7436E-01  1.3136E-01  1.0000E-02  7.8606E-01  5.8881E-01
             2.8916E+00
 PARAMETER:  1.1869E-01 -3.1871E+00 -1.2629E+00  3.4905E-01 -1.2490E+00  7.4030E-02 -1.9298E+00 -9.8248E+00 -1.4072E-01 -4.2966E-01
             1.1618E+00
 GRADIENT:   4.3014E+01  1.1270E-01  2.2061E+01  4.0722E+01  6.0180E+01  3.2869E+00  1.2736E-04  0.0000E+00  1.3137E+00  5.8040E-01
             9.1558E+00

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1367.27139155330        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     4086
 NPARAMETR:  1.0188E+00  3.7647E-02  2.5588E-01  1.2828E+00  2.5946E-01  9.7434E-01  1.2314E-01  1.0000E-02  7.8596E-01  5.8896E-01
             2.8917E+00
 PARAMETER:  1.1865E-01 -3.1795E+00 -1.2631E+00  3.4902E-01 -1.2492E+00  7.4005E-02 -1.9944E+00 -9.8248E+00 -1.4085E-01 -4.2939E-01
             1.1618E+00
 GRADIENT:   1.7098E-01  5.0263E-03  7.2917E-01 -5.0894E-01 -8.5299E-01  9.3349E-03  1.2332E-05  0.0000E+00 -6.7512E-03 -3.8629E-02
            -9.2992E-02

0ITERATION NO.:  121    OBJECTIVE VALUE:  -1367.27139155330        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     4108
 NPARAMETR:  1.0188E+00  3.7647E-02  2.5588E-01  1.2828E+00  2.5946E-01  9.7434E-01  1.2314E-01  1.0000E-02  7.8596E-01  5.8896E-01
             2.8917E+00
 PARAMETER:  1.1865E-01 -3.1795E+00 -1.2631E+00  3.4902E-01 -1.2492E+00  7.4005E-02 -1.9944E+00 -9.8248E+00 -1.4085E-01 -4.2939E-01
             1.1618E+00
 GRADIENT:   1.7098E-01  5.0263E-03  7.2917E-01 -5.0894E-01 -8.5299E-01  9.3349E-03  1.2332E-05  0.0000E+00 -6.7512E-03 -3.8629E-02
            -9.2992E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     4108
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.1479E-04 -5.8141E-05  1.1442E-04 -1.2453E-02  7.0130E-04
 SE:             2.9012E-02  5.8834E-05  2.4534E-04  2.5385E-02  2.0399E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8584E-01  3.2304E-01  6.4093E-01  6.2374E-01  9.7257E-01

 ETASHRINKSD(%)  2.8066E+00  9.9803E+01  9.9178E+01  1.4956E+01  3.1661E+01
 ETASHRINKVR(%)  5.5344E+00  1.0000E+02  9.9993E+01  2.7674E+01  5.3298E+01
 EBVSHRINKSD(%)  2.8096E+00  9.9807E+01  9.9161E+01  1.4282E+01  3.1634E+01
 EBVSHRINKVR(%)  5.5402E+00  1.0000E+02  9.9993E+01  2.6524E+01  5.3261E+01
 RELATIVEINF(%)  6.3659E+01  4.1174E-05  2.0789E-04  8.0605E+00  1.3161E+00
 EPSSHRINKSD(%)  2.9797E+01
 EPSSHRINKVR(%)  5.0716E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1367.2713915533050     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -632.12056498956679     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    51.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1367.271       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  3.76E-02  2.56E-01  1.28E+00  2.59E-01  9.74E-01  1.23E-01  1.00E-02  7.86E-01  5.89E-01  2.89E+00
 


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
+        1.07E+03
 
 TH 2
+       -1.38E+02  1.89E+02
 
 TH 3
+       -8.61E+01  6.03E+02  1.18E+04
 
 TH 4
+       -5.42E+01  2.54E+02 -8.74E+02  8.06E+02
 
 TH 5
+        2.96E+02 -1.44E+03 -1.44E+04 -5.44E+02  2.18E+04
 
 TH 6
+       -1.46E+00 -1.03E+01  2.39E+01 -1.26E+01  1.70E+01  1.87E+02
 
 TH 7
+       -6.25E-03 -1.06E-02 -9.58E-04  1.25E-03  2.29E-03 -2.25E-02  8.21E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.78E+00 -3.94E+01  1.54E+02 -3.19E+01  5.81E+01  8.59E-02 -1.42E-02  0.00E+00  1.71E+02
 
 TH10
+       -1.02E+01  2.56E+00 -1.83E+02  8.30E+00  2.52E+02  2.08E+00  1.16E-03  0.00E+00  4.30E+00  1.10E+02
 
 TH11
+       -1.30E+01 -1.03E+00 -4.53E+00 -8.93E+00 -8.31E+00  3.01E+00  1.52E-03  0.00E+00  1.45E+01  2.83E+01  3.58E+01
 
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
 #CPUT: Total CPU Time in Seconds,       57.312
Stop Time:
Wed Sep 29 13:48:39 CDT 2021
