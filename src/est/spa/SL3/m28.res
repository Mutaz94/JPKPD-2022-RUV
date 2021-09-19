Sat Sep 18 12:45:14 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat28.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1597.51230122901        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.5890E+01 -1.2103E+02 -5.4728E+01 -1.1369E+02  1.2174E+02  1.0373E+01 -2.8860E+00 -1.9316E+00 -1.8755E+01  7.3305E+00
            -6.9161E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1613.66220235975        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.7894E-01  9.9230E-01  1.0551E+00  1.0825E+00  8.9903E-01  9.4937E-01  9.0915E-01  1.0515E+00  1.0735E+00  7.6853E-01
             1.3097E+00
 PARAMETER:  7.8715E-02  9.2271E-02  1.5359E-01  1.7923E-01 -6.4431E-03  4.8046E-02  4.7527E-03  1.5024E-01  1.7092E-01 -1.6328E-01
             3.6980E-01
 GRADIENT:   7.7941E+00 -2.3056E+00  1.2123E+01 -1.5478E+01 -1.2325E+01 -8.8013E+00  3.6941E+00 -3.5201E+00  9.1627E+00 -3.1737E+00
             4.5538E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1617.55349217157        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.7623E-01  9.0813E-01  1.3403E+00  1.1571E+00  9.7063E-01  9.6780E-01  7.5065E-01  1.4473E+00  1.0571E+00  8.9614E-01
             1.1868E+00
 PARAMETER:  7.5942E-02  3.6302E-03  3.9291E-01  2.4596E-01  7.0189E-02  6.7274E-02 -1.8681E-01  4.6967E-01  1.5553E-01 -9.6538E-03
             2.7126E-01
 GRADIENT:   1.2265E+01  4.0255E+00  3.8330E+00  8.8044E+00 -1.3329E+01 -7.6148E-01  3.1229E+00  3.7210E+00  6.5633E+00  2.9531E+00
             1.2102E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1617.98359394365        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.7149E-01  9.8565E-01  1.3034E+00  1.0993E+00  9.9638E-01  9.7066E-01  5.4155E-01  1.4184E+00  1.1147E+00  9.2169E-01
             1.1546E+00
 PARAMETER:  7.1079E-02  8.5543E-02  3.6501E-01  1.9468E-01  9.6371E-02  7.0217E-02 -5.1332E-01  4.4953E-01  2.0858E-01  1.8449E-02
             2.4378E-01
 GRADIENT:   1.7911E+00  2.0562E+00  8.4978E-01  2.1533E+00 -1.4361E+00  1.5683E-01  1.0541E+00 -1.8351E-01  8.2865E-01  1.3547E+00
             6.9094E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1618.53969968095        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.7235E-01  8.8051E-01  1.4568E+00  1.1689E+00  1.0082E+00  9.6901E-01  1.3575E-01  1.5167E+00  1.0875E+00  9.4841E-01
             1.1572E+00
 PARAMETER:  7.1957E-02 -2.7255E-02  4.7626E-01  2.5606E-01  1.0817E-01  6.8518E-02 -1.8969E+00  5.1656E-01  1.8388E-01  4.7035E-02
             2.4600E-01
 GRADIENT:   5.8642E+00  5.3705E+00  1.1332E+00  1.0604E+01 -1.1443E+00 -2.4224E-02  5.5381E-02  3.2451E-01 -2.1747E+00  2.1064E+00
             8.6270E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1619.09264714554        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      479
 NPARAMETR:  9.7899E-01  6.8236E-01  1.4183E+00  1.3048E+00  9.3101E-01  9.7239E-01  3.4838E-02  1.3527E+00  9.8469E-01  8.8879E-01
             1.1591E+00
 PARAMETER:  7.8770E-02 -2.8220E-01  4.4944E-01  3.6602E-01  2.8512E-02  7.2004E-02 -3.2571E+00  4.0208E-01  8.4572E-02 -1.7895E-02
             2.4760E-01
 GRADIENT:  -2.9773E+00  6.3375E+00  3.2190E+00  1.1667E+01 -5.7432E+00 -6.4117E-01  1.4910E-03 -1.1491E+00  5.1973E-01  2.2280E-01
             1.0424E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1619.35374702698        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      655
 NPARAMETR:  9.7894E-01  5.0052E-01  1.3678E+00  1.4134E+00  8.6380E-01  9.7256E-01  1.0000E-02  1.2722E+00  8.9552E-01  8.3914E-01
             1.1603E+00
 PARAMETER:  7.8720E-02 -5.9210E-01  4.1317E-01  4.4599E-01 -4.6419E-02  7.2179E-02 -5.0996E+00  3.4073E-01 -1.0355E-02 -7.5381E-02
             2.4866E-01
 GRADIENT:   1.0082E+00  1.8155E+00  7.7320E-01  4.5084E+00 -1.9053E+00 -2.2437E-01  0.0000E+00  1.4424E-01 -1.9762E-01 -3.4891E-01
             1.4404E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1619.47136640864        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  9.7722E-01  3.1464E-01  1.2695E+00  1.5232E+00  7.8204E-01  9.7266E-01  1.0000E-02  1.1728E+00  8.2070E-01  7.8785E-01
             1.1588E+00
 PARAMETER:  7.6959E-02 -1.0563E+00  3.3859E-01  5.2079E-01 -1.4585E-01  7.2277E-02 -8.2361E+00  2.5938E-01 -9.7602E-02 -1.3845E-01
             2.4738E-01
 GRADIENT:   2.1620E+00  1.6017E+00 -1.6267E-02  8.9689E+00 -1.3004E+00  1.3264E-01  0.0000E+00 -4.1835E-02 -9.6156E-01 -1.7679E-02
             7.2074E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1619.57637494611        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1007
 NPARAMETR:  9.7491E-01  1.9030E-01  1.2296E+00  1.5935E+00  7.3737E-01  9.7222E-01  1.0000E-02  1.1497E+00  7.7774E-01  7.5891E-01
             1.1556E+00
 PARAMETER:  7.4593E-02 -1.5591E+00  3.0671E-01  5.6595E-01 -2.0466E-01  7.1831E-02 -1.2038E+01  2.3946E-01 -1.5136E-01 -1.7587E-01
             2.4459E-01
 GRADIENT:   1.4553E+00  1.1132E+00  1.0460E+00  9.3110E+00 -2.7834E+00  2.3977E-01  0.0000E+00 -1.0794E-01 -1.1988E+00  2.1907E-01
            -5.3603E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1619.72243892659        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1183
 NPARAMETR:  9.7153E-01  8.1381E-02  1.1637E+00  1.6482E+00  6.8979E-01  9.6929E-01  1.0000E-02  1.1158E+00  7.4855E-01  7.1756E-01
             1.1580E+00
 PARAMETER:  7.1116E-02 -2.4086E+00  2.5161E-01  5.9971E-01 -2.7136E-01  6.8807E-02 -1.8881E+01  2.0955E-01 -1.8961E-01 -2.3191E-01
             2.4667E-01
 GRADIENT:  -1.8139E+00  2.8974E-01 -2.2990E-01  5.1789E+00 -4.4557E-01 -6.4236E-01  0.0000E+00  1.7738E-01  4.0657E-01 -7.9248E-01
             1.5416E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1619.80594745445        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1358
 NPARAMETR:  9.7069E-01  2.4756E-02  1.1242E+00  1.6750E+00  6.6286E-01  9.7027E-01  1.0000E-02  1.0799E+00  7.3068E-01  7.1078E-01
             1.1581E+00
 PARAMETER:  7.0248E-02 -3.5987E+00  2.1704E-01  6.1580E-01 -3.1119E-01  6.9820E-02 -2.8847E+01  1.7689E-01 -2.1378E-01 -2.4139E-01
             2.4677E-01
 GRADIENT:  -1.0136E+00  1.0371E-01  1.4877E+00  3.7022E+00 -3.3524E+00 -7.2343E-02  0.0000E+00 -2.2339E-01 -6.8455E-01  1.7742E-01
             2.2453E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1619.83663147906        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1533
 NPARAMETR:  9.7083E-01  1.0000E-02  1.1568E+00  1.6867E+00  6.7265E-01  9.7053E-01  1.0000E-02  1.1143E+00  7.2638E-01  7.1395E-01
             1.1568E+00
 PARAMETER:  7.0400E-02 -4.5613E+00  2.4563E-01  6.2278E-01 -2.9654E-01  7.0082E-02 -3.7112E+01  2.0827E-01 -2.1968E-01 -2.3695E-01
             2.4569E-01
 GRADIENT:   2.4533E-01  0.0000E+00  2.6106E-01  2.1462E+00 -4.8113E-01  8.2871E-02  0.0000E+00 -3.4574E-03 -3.2741E-01  8.4596E-02
            -1.4386E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1619.83852964143        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1708
 NPARAMETR:  9.7073E-01  1.0000E-02  1.1451E+00  1.6844E+00  6.6852E-01  9.7031E-01  1.0000E-02  1.1044E+00  7.2737E-01  7.1113E-01
             1.1572E+00
 PARAMETER:  7.0294E-02 -4.5144E+00  2.3546E-01  6.2143E-01 -3.0270E-01  6.9861E-02 -3.6700E+01  1.9929E-01 -2.1832E-01 -2.4090E-01
             2.4603E-01
 GRADIENT:  -1.5549E-02  0.0000E+00 -2.5407E-04  7.1756E-03 -5.5624E-04  1.1437E-03  0.0000E+00 -3.3803E-04  4.3002E-04  6.2175E-04
             1.4776E-04

0ITERATION NO.:   61    OBJECTIVE VALUE:  -1619.83852964143        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1730
 NPARAMETR:  9.7073E-01  1.0000E-02  1.1451E+00  1.6844E+00  6.6852E-01  9.7031E-01  1.0000E-02  1.1044E+00  7.2737E-01  7.1113E-01
             1.1572E+00
 PARAMETER:  7.0294E-02 -4.5144E+00  2.3546E-01  6.2143E-01 -3.0270E-01  6.9861E-02 -3.6700E+01  1.9929E-01 -2.1832E-01 -2.4090E-01
             2.4603E-01
 GRADIENT:  -1.5549E-02  0.0000E+00 -2.5407E-04  7.1756E-03 -5.5624E-04  1.1437E-03  0.0000E+00 -3.3803E-04  4.3002E-04  6.2175E-04
             1.4776E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1730
 NO. OF SIG. DIGITS IN FINAL EST.:  4.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0819E-04 -2.6700E-06 -1.9697E-02 -5.0086E-03 -2.6137E-02
 SE:             2.9800E-02  1.7925E-06  1.9425E-02  2.9053E-02  1.8995E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9443E-01  1.3635E-01  3.1058E-01  8.6313E-01  1.6882E-01

 ETASHRINKSD(%)  1.6726E-01  9.9994E+01  3.4925E+01  2.6682E+00  3.6365E+01
 ETASHRINKVR(%)  3.3425E-01  1.0000E+02  5.7652E+01  5.2651E+00  5.9506E+01
 EBVSHRINKSD(%)  5.5602E-01  9.9994E+01  3.5940E+01  3.1100E+00  3.5574E+01
 EBVSHRINKVR(%)  1.1089E+00  1.0000E+02  5.8963E+01  6.1232E+00  5.8493E+01
 RELATIVEINF(%)  9.5582E+01  1.7162E-08  4.7838E+00  7.5864E+00  2.8773E+00
 EPSSHRINKSD(%)  4.3802E+01
 EPSSHRINKVR(%)  6.8418E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1619.8385296414344     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -884.68770307769626     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.32
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1619.839       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.71E-01  1.00E-02  1.15E+00  1.68E+00  6.69E-01  9.70E-01  1.00E-02  1.10E+00  7.27E-01  7.11E-01  1.16E+00
 


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
+        1.23E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.55E-02  0.00E+00  2.91E+02
 
 TH 4
+       -1.26E+01  0.00E+00 -3.06E+01  6.95E+02
 
 TH 5
+        7.45E+00  0.00E+00 -6.91E+02 -1.36E+02  2.06E+03
 
 TH 6
+        1.37E+01  0.00E+00  1.99E+00 -3.95E+00 -3.79E+00  2.06E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        4.10E+00  0.00E+00 -3.59E+01 -2.80E+00 -1.71E+01 -2.29E+00  0.00E+00  4.32E+01
 
 TH 9
+        2.41E+00  0.00E+00  8.68E+00 -1.92E+00  2.37E+00 -2.12E+00  0.00E+00  1.77E-01  3.31E+02
 
 TH10
+        3.46E+00  0.00E+00  2.49E+00 -1.03E+00 -1.14E+02 -5.32E+00  0.00E+00  2.93E+01  5.51E+00  8.62E+01
 
 TH11
+       -5.77E+00  0.00E+00 -4.81E+00 -7.67E+00 -1.39E+01  1.09E+00  0.00E+00  1.05E+01  1.23E+01  1.55E+01  1.58E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.986
Stop Time:
Sat Sep 18 12:45:39 CDT 2021
