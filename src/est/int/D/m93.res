Wed Sep 29 10:04:55 CDT 2021
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
$DATA ../../../../data/int/D/dat93.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m93.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   54352.9703050682        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0454E+03  8.9697E+02  7.7237E+01  8.6362E+02  9.8080E+01 -3.6139E+03 -1.9282E+03 -1.0779E+02 -2.7321E+03 -9.6046E+02
            -1.0563E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -490.592224906127        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0289E+00  2.0182E+00  7.8653E-01  1.7803E+00  1.0268E+00  4.6739E+00  4.4459E+00  9.8091E-01  1.6296E+00  1.2407E+00
             1.3072E+01
 PARAMETER:  1.2854E-01  8.0221E-01 -1.4012E-01  6.7676E-01  1.2641E-01  1.6420E+00  1.5920E+00  8.0723E-02  5.8835E-01  3.1565E-01
             2.6705E+00
 GRADIENT:  -2.0976E+01  4.3613E+01 -5.2753E+01  2.0354E+02  2.1260E+01  1.6919E+02  2.1757E+01  4.3487E+00 -1.6544E+01  2.1860E+01
             1.9511E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -591.633998858094        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  8.7802E-01  4.5698E+00  1.4977E+01  2.1656E+00  2.5088E+00  2.7000E+00  1.3610E+01  6.4174E-01  1.4452E+00  1.3863E+00
             1.3487E+01
 PARAMETER: -3.0091E-02  1.6195E+00  2.8065E+00  8.7269E-01  1.0198E+00  1.0933E+00  2.7108E+00 -3.4356E-01  4.6822E-01  4.2666E-01
             2.7018E+00
 GRADIENT:  -5.5484E+01  4.1012E+01 -1.0649E+01  8.6419E+01  1.4028E+01  9.3309E+01  7.6044E+01  8.1554E-02  6.8345E+00  2.1100E+01
             1.1128E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -670.614079796813        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.1049E+00  1.5504E+00  8.6295E+00  1.1725E+00  2.3180E+00  1.7310E+00  4.3328E+00  2.4409E+00  1.0445E+00  4.8809E-01
             1.3940E+01
 PARAMETER:  1.9978E-01  5.3851E-01  2.2552E+00  2.5914E-01  9.4072E-01  6.4871E-01  1.5662E+00  9.9238E-01  1.4358E-01 -6.1725E-01
             2.7347E+00
 GRADIENT:  -1.5804E+00 -4.7862E+00 -4.6200E+00 -2.6489E+00  8.3674E+00 -8.7957E+00 -7.9536E+00  1.3870E+00  1.0231E+01  3.1660E+00
             1.8120E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -694.832680816335        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0673E+00  1.0725E+00  2.3696E+01  1.4142E+00  2.5039E+00  1.9659E+00  6.3883E+00  9.6492E-02  8.8547E-01  3.9871E-01
             1.2445E+01
 PARAMETER:  1.6513E-01  1.7003E-01  3.2653E+00  4.4657E-01  1.0178E+00  7.7594E-01  1.9545E+00 -2.2383E+00 -2.1632E-02 -8.1951E-01
             2.6214E+00
 GRADIENT:   2.6952E+00  7.9579E+00 -1.0901E+00  1.2589E+01  7.2105E-01  1.4056E+01 -9.5502E-01  4.7698E-04 -1.4981E+00  1.6842E+00
            -4.1045E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -699.589500647236        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  1.0687E+00  5.8746E-01  5.5025E+01  1.6393E+00  2.6430E+00  1.8630E+00  7.3907E+00  1.2499E-01  1.0803E+00  1.0896E-01
             1.2676E+01
 PARAMETER:  1.6648E-01 -4.3194E-01  4.1078E+00  5.9429E-01  1.0719E+00  7.2218E-01  2.1002E+00 -1.9795E+00  1.7721E-01 -2.1168E+00
             2.6397E+00
 GRADIENT:  -2.0073E-01 -1.3654E-01 -9.0792E-01  3.6344E+00  7.6947E+00  8.8737E-01  2.3397E+00  1.9506E-04 -1.6269E+00  1.1993E-01
            -8.8660E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -700.569653380248        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      694
 NPARAMETR:  1.0635E+00  5.6774E-01  1.8039E+02  1.6745E+00  2.5904E+00  1.8566E+00  7.9880E+00  1.1331E-01  1.1225E+00  2.8982E-02
             1.2706E+01
 PARAMETER:  1.6160E-01 -4.6608E-01  5.2951E+00  6.1553E-01  1.0518E+00  7.1873E-01  2.1779E+00 -2.0777E+00  2.1553E-01 -3.4411E+00
             2.6421E+00
 GRADIENT:  -2.4895E+00  3.1234E+00 -1.1288E-01  1.1736E+00 -5.4275E+00  3.9348E-01  1.6061E+01  1.2885E-05 -1.1782E+00  8.0242E-03
            -4.5138E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -701.011818653204        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      873
 NPARAMETR:  1.0706E+00  4.9802E-01  1.2528E+03  1.7130E+00  2.6712E+00  1.8561E+00  7.9206E+00  1.1373E-01  1.1780E+00  1.0000E-02
             1.2751E+01
 PARAMETER:  1.6819E-01 -5.9711E-01  7.2331E+00  6.3823E-01  1.0825E+00  7.1847E-01  2.1695E+00 -2.0740E+00  2.6385E-01 -5.3893E+00
             2.6456E+00
 GRADIENT:  -3.7844E-01  5.9164E-01 -2.7602E-02 -1.1123E+00  1.2383E+00  4.7746E-01  9.4820E+00  1.1211E-06  4.9362E-01  0.0000E+00
             1.7088E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -701.232974989897        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1035
 NPARAMETR:  1.0721E+00  4.7853E-01  4.4763E+04  1.7220E+00  2.6623E+00  1.8528E+00  8.3076E+00  1.1855E-01  1.1806E+00  1.0000E-02
             1.2740E+01
 PARAMETER:  1.6957E-01 -6.3704E-01  1.0809E+01  6.4349E-01  1.0792E+00  7.1668E-01  2.2172E+00 -2.0324E+00  2.6598E-01 -5.5508E+00
             2.6447E+00
 GRADIENT:  -9.2748E-01  1.4054E+00 -7.4687E-04 -5.3251E+00  1.9302E-01  8.2246E-01  1.8539E+01  1.5824E-07  8.2333E-01  0.0000E+00
             2.3847E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -701.428853486399        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1226
 NPARAMETR:  1.0705E+00  4.2885E-01  1.1248E+05  1.7459E+00  2.6630E+00  1.8518E+00  8.5206E+00  1.1364E-01  1.1901E+00  1.0000E-02
             1.2734E+01
 PARAMETER:  1.6815E-01 -7.4664E-01  1.1731E+01  6.5726E-01  1.0794E+00  7.1616E-01  2.2425E+00 -2.0748E+00  2.7405E-01 -5.5508E+00
             2.6443E+00
 GRADIENT:  -2.0380E+01  8.1324E-01 -3.2578E-04 -2.9255E-01  1.3559E-01  8.2999E+00  1.9375E+01 -1.8756E-06 -1.8020E+00  0.0000E+00
             7.8207E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -701.602978049963        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:     1427             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0701E+00  3.7688E-01  1.5536E+05  1.7843E+00  2.6610E+00  1.8524E+00  8.6751E+00  1.1721E-01  1.2154E+00  1.0000E-02
             1.2729E+01
 PARAMETER:  1.6776E-01 -8.7583E-01  1.2054E+01  6.7905E-01  1.0787E+00  7.1646E-01  2.2605E+00 -2.0438E+00  2.9507E-01 -5.5508E+00
             2.6439E+00
 GRADIENT:  -2.4947E+01  3.6099E+00 -2.4233E-04  2.3539E+01  1.9369E+00  2.0830E+01  2.0053E+02 -1.4719E-06 -5.2559E-01  0.0000E+00
             6.2049E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -701.647194258112        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1612
 NPARAMETR:  1.0706E+00  3.5677E-01  1.6137E+05  1.8067E+00  2.6610E+00  1.8517E+00  8.8043E+00  1.1479E-01  1.2196E+00  1.0000E-02
             1.2732E+01
 PARAMETER:  1.6820E-01 -9.3067E-01  1.2091E+01  6.9148E-01  1.0787E+00  7.1609E-01  2.2752E+00 -2.0647E+00  2.9852E-01 -5.5508E+00
             2.6441E+00
 GRADIENT:  -4.7457E+02 -5.6515E+00 -2.1006E-04  3.3050E+01  1.2386E+01  1.5463E+02  1.2933E+01  7.8227E-05 -8.2536E-01  0.0000E+00
             1.8917E+02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -701.700478360587        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1811             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0694E+00  3.4542E-01  2.0596E+05  1.8176E+00  2.6624E+00  1.8511E+00  8.9608E+00  1.1839E-01  1.2329E+00  1.0000E-02
             1.2741E+01
 PARAMETER:  1.6714E-01 -9.6299E-01  1.2335E+01  6.9751E-01  1.0792E+00  7.1580E-01  2.2929E+00 -2.0338E+00  3.0934E-01 -5.5508E+00
             2.6449E+00
 GRADIENT:  -1.2651E+03  2.0047E+00 -2.1013E-04  3.9631E+02 -2.4531E+00  6.3846E+02  1.9309E+02 -2.4327E-04 -1.5397E+02  0.0000E+00
             4.2629E+02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -701.715599386100        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1969
 NPARAMETR:  1.0728E+00  3.3206E-01  2.0496E+05  1.8126E+00  2.6704E+00  1.8493E+00  9.0133E+00  1.1845E-01  1.2524E+00  1.0000E-02
             1.2764E+01
 PARAMETER:  1.7030E-01 -1.0024E+00  1.2331E+01  6.9476E-01  1.0822E+00  7.1478E-01  2.2987E+00 -2.0332E+00  3.2509E-01 -5.5508E+00
             2.6466E+00
 GRADIENT:  -6.5967E+01 -1.6163E+00 -2.3315E-04  6.6765E+00  1.1665E+00  2.5099E+01  2.1864E+01 -5.9126E-06 -2.6062E+00  0.0000E+00
             3.0888E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -701.744665913474        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2165             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0716E+00  3.2414E-01  2.6850E+05  1.8307E+00  2.6638E+00  1.8503E+00  9.1041E+00  1.1569E-01  1.2473E+00  1.0000E-02
             1.2746E+01
 PARAMETER:  1.6916E-01 -1.0266E+00  1.2601E+01  7.0470E-01  1.0798E+00  7.1535E-01  2.3087E+00 -2.0568E+00  3.2100E-01 -5.5508E+00
             2.6452E+00
 GRADIENT:  -5.5412E+01  5.3162E-02 -1.9122E-04  1.9527E+01  3.1731E+00  2.4432E+01  2.2089E+02  1.2251E-06  2.5008E+00  0.0000E+00
             7.6343E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -701.759000509716        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2351
 NPARAMETR:  1.0713E+00  3.1764E-01  2.6006E+05  1.8330E+00  2.6655E+00  1.8503E+00  9.1625E+00  1.0314E-01  1.2562E+00  1.0000E-02
             1.2751E+01
 PARAMETER:  1.6885E-01 -1.0468E+00  1.2569E+01  7.0597E-01  1.0804E+00  7.1533E-01  2.3151E+00 -2.1717E+00  3.2811E-01 -5.5508E+00
             2.6456E+00
 GRADIENT:  -3.6148E+02 -6.3663E+00 -1.7609E-04  1.3738E+01  1.0279E+01  1.1823E+02  2.0148E+01  3.7143E-06  4.1161E+00  0.0000E+00
             1.5261E+02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -701.767258770872        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     2546             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0714E+00  3.0610E-01  3.1779E+05  1.8370E+00  2.6685E+00  1.8501E+00  9.2265E+00  1.5964E-01  1.2649E+00  1.0000E-02
             1.2760E+01
 PARAMETER:  1.6900E-01 -1.0838E+00  1.2769E+01  7.0815E-01  1.0815E+00  7.1522E-01  2.3221E+00 -1.7348E+00  3.3496E-01 -5.5508E+00
             2.6463E+00
 GRADIENT:  -5.1631E+01  3.3264E+00 -1.4680E-04  3.7063E+01  2.6198E+00  1.3721E+01  2.2700E+02 -1.4835E-06  3.7636E-01  0.0000E+00
             7.1471E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -701.771735338311        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2734
 NPARAMETR:  1.0715E+00  3.0406E-01  3.8512E+05  1.8406E+00  2.6674E+00  1.8499E+00  9.2585E+00  1.6482E-01  1.2645E+00  1.0000E-02
             1.2757E+01
 PARAMETER:  1.6906E-01 -1.0905E+00  1.2961E+01  7.1011E-01  1.0811E+00  7.1514E-01  2.3255E+00 -1.7029E+00  3.3467E-01 -5.5508E+00
             2.6461E+00
 GRADIENT:  -3.3914E+01 -3.1627E+00 -1.3614E-04 -8.7031E-01  2.1795E+00 -1.8072E+00  2.5214E+01 -1.0993E-05  3.7471E+00  0.0000E+00
             1.9137E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -701.775160505523        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2928             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0717E+00  2.9934E-01  3.4691E+05  1.8468E+00  2.6665E+00  1.8498E+00  9.3182E+00  1.5834E-01  1.2659E+00  1.0000E-02
             1.2756E+01
 PARAMETER:  1.6926E-01 -1.1062E+00  1.2857E+01  7.1343E-01  1.0808E+00  7.1506E-01  2.3320E+00 -1.7430E+00  3.3574E-01 -5.5508E+00
             2.6460E+00
 GRADIENT:  -4.8506E+01  5.7452E+00 -1.3834E-04  3.7977E+01  2.2876E+00  7.4570E+00  2.3046E+02  1.1147E-05 -2.4864E-01  0.0000E+00
             6.3413E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -701.776331617707        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3113
 NPARAMETR:  1.0717E+00  2.9714E-01  4.3298E+05  1.8482E+00  2.6668E+00  1.8496E+00  9.3365E+00  1.7790E-01  1.2677E+00  1.0000E-02
             1.2757E+01
 PARAMETER:  1.6929E-01 -1.1136E+00  1.3078E+01  7.1419E-01  1.0809E+00  7.1498E-01  2.3339E+00 -1.6265E+00  3.3724E-01 -5.5508E+00
             2.6461E+00
 GRADIENT:  -3.2689E+00 -4.8350E+00 -1.0673E-04 -2.2066E+01  1.4514E+00  3.9844E+01  2.1912E+01  2.1512E-06 -4.7005E+00  0.0000E+00
             1.8680E+01

0ITERATION NO.:   98    OBJECTIVE VALUE:  -701.776576539362        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:     3217
 NPARAMETR:  1.0718E+00  2.9745E-01  3.9859E+05  1.8519E+00  2.6652E+00  1.8492E+00  9.3749E+00  1.7467E-01  1.2686E+00  1.0000E-02
             1.2753E+01
 PARAMETER:  1.6922E-01 -1.1237E+00  1.3019E+01  7.1389E-01  1.0816E+00  7.1501E-01  2.3342E+00 -1.6428E+00  3.3841E-01 -5.5508E+00
             2.6465E+00
 GRADIENT:  -3.9242E-02 -5.8332E-02  2.2670E-04 -6.8446E-01  1.4816E-01  2.3812E-02 -2.9199E-01  4.9870E-04  2.9435E-02  0.0000E+00
             4.8394E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3217
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7221E-02  7.0125E-02  6.9822E-09 -9.2015E-02  1.7006E-06
 SE:             2.7449E-02  2.0646E-02  7.8011E-09  1.6269E-02  9.0763E-05
 N:                     100         100         100         100         100

 P VAL.:         5.3040E-01  6.8271E-04  3.7077E-01  1.5540E-08  9.8505E-01

 ETASHRINKSD(%)  8.0430E+00  3.0832E+01  1.0000E+02  4.5498E+01  9.9696E+01
 ETASHRINKVR(%)  1.5439E+01  5.2157E+01  1.0000E+02  7.0296E+01  9.9999E+01
 EBVSHRINKSD(%)  9.8169E+00  3.0664E+01  1.0000E+02  3.7406E+01  9.9622E+01
 EBVSHRINKVR(%)  1.8670E+01  5.1926E+01  1.0000E+02  6.0820E+01  9.9999E+01
 RELATIVEINF(%)  8.1012E+01  3.0047E+01  0.0000E+00  2.3546E+01  2.7803E-04
 EPSSHRINKSD(%)  3.8225E+00
 EPSSHRINKVR(%)  7.4989E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -701.77657653936205     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       952.31278322904870     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   125.04
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    18.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -701.777       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.07E+00  2.94E-01  4.08E+05  1.85E+00  2.67E+00  1.85E+00  9.34E+00  1.75E-01  1.27E+00  1.00E-02  1.28E+01
 


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
+        4.95E+02
 
 TH 2
+       -1.57E+02  1.33E+02
 
 TH 3
+       -2.19E-07  2.20E-07  4.64E-16
 
 TH 4
+       -2.21E+02  9.36E+01  5.61E-08  2.03E+02
 
 TH 5
+        2.96E+01 -1.07E+01 -4.02E-09 -2.89E+01  4.35E+00
 
 TH 6
+       -2.37E+01 -3.30E+00  6.83E-09  1.27E+01 -3.85E+00  2.36E+01
 
 TH 7
+       -1.38E-01  6.41E+00  1.60E-08 -3.49E+00  7.59E-01 -1.43E+00  7.89E-01
 
 TH 8
+        7.32E+00 -5.42E+00 -8.62E-09 -3.94E+00  4.18E-01  5.00E-01 -2.58E-01  2.34E-01
 
 TH 9
+        1.69E+01 -6.90E+00  1.98E-08 -3.26E+01  4.15E+00  6.89E+00  1.10E+00  4.39E-01  1.06E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.03E-01 -2.81E+00  1.13E-09 -7.15E+00  9.38E-01  9.99E-01  1.19E-01  1.09E-01  2.06E+00  0.00E+00  7.37E-01
 
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
+        2.27E+02
 
 TH 2
+       -7.46E+00  8.65E+01
 
 TH 3
+       -2.58E-08  1.64E-07  3.00E-14
 
 TH 4
+       -8.72E+00  2.54E+01 -4.12E-08  1.14E+02
 
 TH 5
+       -8.78E-01 -3.96E+00  1.94E-08 -1.17E+01  2.09E+01
 
 TH 6
+       -3.14E+00 -2.27E+00 -4.29E-09 -5.27E-01 -3.67E-01  3.82E+01
 
 TH 7
+        4.38E-01  6.43E+00  6.45E-10 -3.93E+00  2.06E-01  6.55E-02  1.15E+00
 
 TH 8
+        1.47E+00 -3.12E+00  1.43E-07 -7.58E-01  2.15E-01  5.26E-01 -6.35E-03  4.26E+00
 
 TH 9
+       -2.81E+00  1.60E+00  8.77E-08 -2.72E+01  3.64E+00 -2.00E+00  3.33E-01  4.10E-01  3.54E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.75E+00 -2.43E+00 -4.30E-10 -8.61E+00  7.12E-01  1.57E+00  1.55E-01  1.22E-02  2.96E+00  0.00E+00  5.43E+00
 
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
+        2.31E+02
 
 TH 2
+        6.65E+01  8.78E+01
 
 TH 3
+       -8.85E-10  5.99E-10  5.66E-19
 
 TH 4
+        9.93E+01  2.29E+01 -1.02E-09  1.17E+02
 
 TH 5
+       -1.33E+01  2.10E+00 -2.45E-11 -2.05E+01  1.48E+01
 
 TH 6
+        1.78E+01  1.31E+01 -5.36E-11 -1.35E+01  5.72E-01  4.06E+01
 
 TH 7
+       -1.37E+00  5.98E+00  6.31E-11 -3.79E+00  1.49E+00  6.93E-01  9.31E-01
 
 TH 8
+       -3.44E-02  5.47E-03  2.19E-12 -3.11E-03 -3.23E-04 -3.11E-03  2.37E-04  1.87E-04
 
 TH 9
+       -2.24E+01 -6.09E-01 -5.11E-10 -2.51E+01  5.23E+00  1.58E+01  4.41E-01  3.12E-03  2.72E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.15E+01 -7.95E+00 -9.06E-10 -3.53E+01  5.61E+00  6.46E+00  1.43E+00 -7.53E-03  1.21E+01  0.00E+00  1.10E+02
 
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
 #CPUT: Total CPU Time in Seconds,      143.571
Stop Time:
Wed Sep 29 10:07:20 CDT 2021
