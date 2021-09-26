Sat Sep 25 00:00:58 CDT 2021
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
$DATA ../../../../data/int/SL1/dat16.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m16.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3339.51235828344        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.0271E+01 -6.1988E+01 -2.7076E+01  3.3601E+01  1.3628E+02  4.5213E+01 -7.3363E+01 -5.7769E+01 -1.8357E+01 -1.5482E+01
            -1.0572E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3525.21592407610        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0244E+00  1.1374E+00  1.0909E+00  9.5926E-01  1.0245E+00  7.9193E-01  1.2806E+00  1.0365E+00  9.3605E-01  1.0257E+00
             1.4061E+00
 PARAMETER:  1.2414E-01  2.2876E-01  1.8698E-01  5.8407E-02  1.2418E-01 -1.3329E-01  3.4737E-01  1.3589E-01  3.3911E-02  1.2533E-01
             4.4082E-01
 GRADIENT:  -1.6609E+01  3.2120E+01 -1.1360E+01  2.0584E+01  1.7014E+01 -4.3879E+01  1.2072E+01 -5.3990E+00 -6.3717E+00  2.1911E-01
            -5.8521E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3526.54213130000        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      249
 NPARAMETR:  1.0420E+00  1.2692E+00  1.2588E+00  8.9701E-01  1.1594E+00  8.0827E-01  1.1220E+00  1.3993E+00  9.8232E-01  1.1438E+00
             1.3993E+00
 PARAMETER:  1.4110E-01  3.3838E-01  3.3013E-01 -8.6914E-03  2.4787E-01 -1.1286E-01  2.1514E-01  4.3598E-01  8.2158E-02  2.3437E-01
             4.3597E-01
 GRADIENT:   1.3324E+01  3.2120E+01  1.4640E+00  1.4235E+01  1.6950E+01 -3.6298E+01  4.8587E+00 -8.5820E+00  7.0349E+00 -8.0229E+00
            -6.5878E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3531.28063494590        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      425
 NPARAMETR:  1.0379E+00  1.2375E+00  1.3484E+00  9.0453E-01  1.1694E+00  8.7921E-01  1.0743E+00  1.8295E+00  8.8572E-01  1.2043E+00
             1.4326E+00
 PARAMETER:  1.3719E-01  3.1306E-01  3.9890E-01 -3.3529E-04  2.5646E-01 -2.8729E-02  1.7168E-01  7.0406E-01 -2.1355E-02  2.8591E-01
             4.5949E-01
 GRADIENT:  -4.2289E-02  8.5273E-01 -4.6759E+00 -2.7858E+00  4.7389E+00 -2.8443E-01 -3.4776E+00 -1.0813E+00  1.2058E+00 -2.2655E+00
             4.4632E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3533.20443801642        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      607
 NPARAMETR:  1.0370E+00  1.2244E+00  1.5813E+00  9.1721E-01  1.2168E+00  8.8860E-01  1.0981E+00  2.4498E+00  7.9986E-01  1.2232E+00
             1.4247E+00
 PARAMETER:  1.3636E-01  3.0248E-01  5.5824E-01  1.3585E-02  2.9621E-01 -1.8111E-02  1.9360E-01  9.9600E-01 -1.2331E-01  3.0150E-01
             4.5398E-01
 GRADIENT:  -1.9752E+00 -1.9770E+00 -9.6147E+00 -1.6220E+01  4.4166E+00  3.7403E+00 -5.9031E+00  1.6385E+00  3.8827E+00 -2.1731E+00
             1.7020E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3533.53078501058        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  1.0367E+00  1.2218E+00  1.5864E+00  9.3482E-01  1.2054E+00  8.7247E-01  1.1766E+00  2.3941E+00  7.2156E-01  1.2509E+00
             1.4020E+00
 PARAMETER:  1.3601E-01  3.0032E-01  5.6146E-01  3.2595E-02  2.8683E-01 -3.6427E-02  2.6267E-01  9.7300E-01 -2.2633E-01  3.2390E-01
             4.3791E-01
 GRADIENT:   2.8717E+01  3.5749E+01 -9.8784E+00  1.1692E+01  5.2275E+00 -1.6904E+00  2.4275E+00 -5.1769E+00 -2.3859E+00  3.4577E+00
            -1.5897E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3533.84635003887        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  1.0334E+00  1.1890E+00  1.5867E+00  9.3604E-01  1.2055E+00  8.7873E-01  1.1682E+00  2.3940E+00  7.4921E-01  1.2237E+00
             1.4023E+00
 PARAMETER:  1.3281E-01  2.7314E-01  5.6164E-01  3.3906E-02  2.8686E-01 -2.9281E-02  2.5545E-01  9.7297E-01 -1.8873E-01  3.0187E-01
             4.3808E-01
 GRADIENT:   1.9258E+01  6.2363E+00 -1.0975E+01 -1.2378E+01  1.6231E+01  1.1480E+00  5.7254E-01 -3.9002E+00  1.7980E-01  1.6587E+00
            -1.2549E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3533.95419149007        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  1.0350E+00  1.1822E+00  1.5867E+00  9.4730E-01  1.2055E+00  8.8156E-01  1.1701E+00  2.3939E+00  7.5297E-01  1.2071E+00
             1.4023E+00
 PARAMETER:  1.3436E-01  2.6735E-01  5.6167E-01  4.5864E-02  2.8687E-01 -2.6059E-02  2.5706E-01  9.7293E-01 -1.8373E-01  2.8818E-01
             4.3810E-01
 GRADIENT:  -7.1662E+00 -2.2746E-01 -1.4585E+01 -1.0274E+00  1.7462E+01  5.5930E-01 -1.4449E+00 -5.3301E+00  3.6503E-01 -4.6722E-01
            -1.2148E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3533.98082130586        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1203
 NPARAMETR:  1.0363E+00  1.1782E+00  1.5873E+00  9.4935E-01  1.2053E+00  8.7842E-01  1.1819E+00  2.3945E+00  7.4939E-01  1.2084E+00
             1.4025E+00
 PARAMETER:  1.3564E-01  2.6398E-01  5.6202E-01  4.8027E-02  2.8675E-01 -2.9636E-02  2.6709E-01  9.7315E-01 -1.8849E-01  2.8931E-01
             4.3827E-01
 GRADIENT:  -3.5579E+00 -1.0329E+00 -1.5169E+01 -7.0759E-01  1.8825E+01 -8.2060E-01 -9.0752E-02 -5.6384E+00  4.6814E-01  1.1716E-01
            -1.1218E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3534.40566675805        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     1375
 NPARAMETR:  1.0336E+00  1.1675E+00  1.6056E+00  9.5446E-01  1.1928E+00  8.7445E-01  1.1888E+00  2.4503E+00  7.4909E-01  1.1891E+00
             1.4095E+00
 PARAMETER:  1.3305E-01  2.5484E-01  5.7349E-01  5.3392E-02  2.7629E-01 -3.4164E-02  2.7297E-01  9.9623E-01 -1.8890E-01  2.7319E-01
             4.4324E-01
 GRADIENT:   1.9390E+01  9.7647E+00 -1.2502E+01 -2.8326E+00  1.6775E+01 -8.0158E-01  1.1394E+00 -2.1081E+00  1.4832E+00  1.2672E-01
             3.3747E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3534.45683094665        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     1504             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0331E+00  1.1647E+00  1.6060E+00  9.5586E-01  1.1901E+00  8.8480E-01  1.1967E+00  2.4571E+00  7.3494E-01  1.1945E+00
             1.4093E+00
 PARAMETER:  1.3261E-01  2.5244E-01  5.7377E-01  5.4858E-02  2.7406E-01 -2.2394E-02  2.7958E-01  9.9899E-01 -2.0796E-01  2.7771E-01
             4.4306E-01
 GRADIENT:   1.8375E+01  9.5704E+00 -1.2802E+01 -4.0144E+00  1.4172E+01  3.8018E+00  1.4684E+00 -2.4482E+00  1.6995E-01  9.1510E-01
             3.5121E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3534.46265142802        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1638
 NPARAMETR:  1.0331E+00  1.1647E+00  1.6061E+00  9.5586E-01  1.1901E+00  8.8008E-01  1.1957E+00  2.4571E+00  7.3798E-01  1.1929E+00
             1.4093E+00
 PARAMETER:  1.3261E-01  2.5244E-01  5.7378E-01  5.4857E-02  2.7406E-01 -2.7740E-02  2.7872E-01  9.9899E-01 -2.0383E-01  2.7637E-01
             4.4306E-01
 GRADIENT:  -1.2137E+01 -1.1942E+00 -1.4157E+01 -8.5329E+00  8.2125E+00 -9.8771E-02 -6.9700E-02 -3.3196E+00  2.2830E-02 -2.0546E-02
             2.7405E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3534.49609954483        NO. OF FUNC. EVALS.: 206
 CUMULATIVE NO. OF FUNC. EVALS.:     1844             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0337E+00  1.1648E+00  1.6070E+00  9.5641E-01  1.1890E+00  8.8011E-01  1.1957E+00  2.4624E+00  7.3794E-01  1.1929E+00
             1.4090E+00
 PARAMETER:  1.3315E-01  2.5255E-01  5.7435E-01  5.5431E-02  2.7312E-01 -2.7714E-02  2.7877E-01  1.0011E+00 -2.0390E-01  2.7637E-01
             4.4290E-01
 GRADIENT:   1.9902E+01  1.0521E+01 -1.2687E+01 -3.3504E+00  1.3413E+01  1.7412E+00  1.4135E+00 -2.0994E+00  5.1591E-01  8.8140E-01
             3.4414E+00

0ITERATION NO.:   62    OBJECTIVE VALUE:  -3534.49609954483        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1910
 NPARAMETR:  1.0337E+00  1.1648E+00  1.6069E+00  9.5640E-01  1.1890E+00  8.8015E-01  1.1959E+00  2.4622E+00  7.3779E-01  1.1928E+00
             1.4091E+00
 PARAMETER:  1.3315E-01  2.5255E-01  5.7435E-01  5.5431E-02  2.7312E-01 -2.7714E-02  2.7877E-01  1.0011E+00 -2.0390E-01  2.7637E-01
             4.4290E-01
 GRADIENT:  -4.1107E+03  4.3225E+03  1.6626E+03  1.0913E+04 -3.9892E+03 -6.7252E-02 -8.7577E-02  1.0898E+03  5.9795E-02  6.9639E-02
            -2.4656E+03
 NUMSIGDIG:         3.3         3.3         3.8         3.3         3.3         2.4         2.5         3.3         2.1         2.5
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1910
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.8134E-03 -1.6009E-02 -1.5366E-02  1.6249E-02 -3.6385E-02
 SE:             2.9802E-02  2.4379E-02  2.2002E-02  2.2242E-02  2.3429E-02
 N:                     100         100         100         100         100

 P VAL.:         8.7169E-01  5.1141E-01  4.8495E-01  4.6505E-01  1.2043E-01

 ETASHRINKSD(%)  1.5825E-01  1.8326E+01  2.6291E+01  2.5486E+01  2.1510E+01
 ETASHRINKVR(%)  3.1626E-01  3.3294E+01  4.5670E+01  4.4477E+01  3.8393E+01
 EBVSHRINKSD(%)  5.8890E-01  1.8619E+01  2.9351E+01  2.9149E+01  1.7946E+01
 EBVSHRINKVR(%)  1.1743E+00  3.3771E+01  5.0087E+01  4.9801E+01  3.2672E+01
 RELATIVEINF(%)  9.8816E+01  2.7302E+01  3.8408E+01  2.0384E+01  4.1003E+01
 EPSSHRINKSD(%)  2.0639E+01
 EPSSHRINKVR(%)  3.7019E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3534.4960995448296     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1880.4067397764188     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    57.04
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3534.496       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.16E+00  1.61E+00  9.56E-01  1.19E+00  8.80E-01  1.20E+00  2.46E+00  7.38E-01  1.19E+00  1.41E+00
 


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
+        1.44E+07
 
 TH 2
+       -6.74E+06  3.15E+06
 
 TH 3
+        1.04E+02  1.70E+02  9.22E+05
 
 TH 4
+       -2.07E+07  9.70E+06 -2.07E+02  2.99E+07
 
 TH 5
+        6.56E+01  9.07E+02 -2.48E+03 -3.55E+03  2.59E+06
 
 TH 6
+        8.72E+01 -3.66E+01  2.49E+01 -9.15E+01  3.00E+01  2.51E+02
 
 TH 7
+       -4.75E+02  2.48E+02 -2.27E+02  6.65E+02 -1.99E+02 -6.09E-01  6.47E+01
 
 TH 8
+       -9.35E+00 -1.32E+02 -3.45E+02  4.77E+02  2.31E+02 -4.02E+00  2.61E+01  4.51E+04
 
 TH 9
+        4.97E+03 -2.33E+03  4.88E+02 -7.13E+03  2.13E+03  4.29E+00  2.24E+01 -2.69E+02  8.34E+01
 
 TH10
+       -6.02E+06  7.30E+01  9.59E+02  8.66E+06 -1.09E+02  9.23E-01  5.84E+00  6.85E+00  5.50E+06  2.51E+06
 
 TH11
+        2.66E+01  4.97E+02 -1.85E+03 -1.91E+03 -9.55E+02  1.85E+01 -9.78E+01 -4.09E+02  1.11E+03 -2.47E+01  7.03E+05
 
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
 #CPUT: Total CPU Time in Seconds,       71.847
Stop Time:
Sat Sep 25 00:02:12 CDT 2021
