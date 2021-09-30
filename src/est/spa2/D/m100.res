Thu Sep 30 10:13:54 CDT 2021
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
$DATA ../../../../data/spa2/D/dat100.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   45366.9421142132        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3570E+03  9.5566E+02  2.0900E+01  9.9092E+02 -4.3076E+01 -3.1731E+03 -1.7412E+03 -4.6027E+01 -2.3447E+03 -6.6126E+02
            -8.6450E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -234.315093460706        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0982E+00  1.2251E+00  1.0018E+00  1.2674E+00  1.1732E+00  2.0209E+00  1.3299E+00  9.6848E-01  1.0946E+00  8.8752E-01
             1.4588E+01
 PARAMETER:  1.9370E-01  3.0299E-01  1.0176E-01  3.3695E-01  2.5975E-01  8.0355E-01  3.8513E-01  6.7975E-02  1.9043E-01 -1.9321E-02
             2.7802E+00
 GRADIENT:  -1.7440E+01  2.0712E+01 -6.1186E+00  6.0290E+01  3.7234E-01  5.0609E+01 -3.0336E+01  2.7829E+00 -2.1176E+01  6.7427E+00
            -1.1967E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -328.427268533834        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.2050E+00  1.5530E+00  8.6014E+00  1.4587E+00  6.8153E+00  1.4764E+00  5.0246E+00  3.6829E-01  8.5850E-01  1.3262E-02
             1.6117E+01
 PARAMETER:  2.8650E-01  5.4016E-01  2.2519E+00  4.7752E-01  2.0192E+00  4.8963E-01  1.7143E+00 -8.9887E-01 -5.2568E-02 -4.2228E+00
             2.8799E+00
 GRADIENT:   6.8747E+01  6.8708E+00 -8.6867E-01  1.8702E+01  2.2202E+00 -1.2215E+01 -9.4667E+00  1.9008E-04  3.3201E+00  3.8103E-06
             1.4412E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -339.515773987596        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0400E+00  1.1289E+00  4.1380E+00  1.2742E+00  2.1263E+00  1.6513E+00  5.3765E+00  4.4339E-01  5.4950E-01  1.4950E-01
             1.5029E+01
 PARAMETER:  1.3919E-01  2.2127E-01  1.5202E+00  3.4231E-01  8.5439E-01  6.0158E-01  1.7820E+00 -7.1330E-01 -4.9875E-01 -1.8005E+00
             2.8100E+00
 GRADIENT:  -8.3862E+00 -7.7851E+00 -2.7566E+00 -5.6056E+00  3.4247E+00  1.9790E+01  1.0230E+01  6.6858E-02  1.2544E+00  1.1616E-01
             3.1756E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -342.240877734994        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      347
 NPARAMETR:  1.0504E+00  1.1144E+00  5.7875E+00  1.3446E+00  2.4706E+00  1.5500E+00  5.8642E+00  4.9827E-01  5.3577E-01  1.6734E-01
             1.5056E+01
 PARAMETER:  1.4917E-01  2.0835E-01  1.8557E+00  3.9609E-01  1.0044E+00  5.3828E-01  1.8689E+00 -5.9661E-01 -5.2405E-01 -1.6878E+00
             2.8118E+00
 GRADIENT:  -3.9580E+00 -2.7979E+00 -2.4154E+00 -3.4902E+00  5.6186E+00  4.3018E+00 -8.5317E+00  3.8133E-02 -6.2212E-02  1.0648E-01
            -1.7298E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -344.982024876608        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      522
 NPARAMETR:  1.0604E+00  7.2662E-01  1.5193E+01  1.6380E+00  2.4968E+00  1.4481E+00  7.1016E+00  9.0539E-01  9.0497E-01  2.3226E-02
             1.5529E+01
 PARAMETER:  1.5865E-01 -2.1935E-01  2.8208E+00  5.9349E-01  1.0150E+00  4.7027E-01  2.0603E+00  6.1191E-04  1.4940E-04 -3.6625E+00
             2.8427E+00
 GRADIENT:  -1.7412E+00 -6.8467E-01 -6.4143E-01 -8.7076E-01  7.3597E-01 -1.2268E+00  6.1357E-01  2.1416E-02  7.2025E-01  2.1594E-03
            -1.6001E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -345.406701047455        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  1.0664E+00  6.8930E-01  1.0305E+02  1.6720E+00  2.6208E+00  1.4493E+00  7.2428E+00  9.9186E-01  9.2564E-01  5.3478E-02
             1.5609E+01
 PARAMETER:  1.6431E-01 -2.7208E-01  4.7353E+00  6.1401E-01  1.0635E+00  4.7109E-01  2.0800E+00  9.1832E-02  2.2726E-02 -2.8285E+00
             2.8478E+00
 GRADIENT:  -4.2150E-01 -9.2031E-01 -5.1639E-02 -2.1936E+00  4.8931E-01 -5.7670E-01  1.2404E+00  4.7069E-04  5.6719E-01  1.1414E-02
             2.5899E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -345.478948541851        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  1.0683E+00  7.0452E-01  4.9336E+02  1.6767E+00  2.6384E+00  1.4557E+00  7.2987E+00  9.9180E-01  9.1648E-01  2.2513E-02
             1.5625E+01
 PARAMETER:  1.6611E-01 -2.5024E-01  6.3012E+00  6.1685E-01  1.0702E+00  4.7546E-01  2.0877E+00  9.1767E-02  1.2788E-02 -3.6936E+00
             2.8489E+00
 GRADIENT:  -1.7110E-01 -3.5264E-02 -9.4785E-03 -7.0077E-03  9.4294E-02  2.1331E-01  3.0318E+00  1.8902E-05  4.0046E-02  2.0163E-03
             1.5339E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -345.492420218763        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  1.0688E+00  7.0368E-01  6.1242E+03  1.6778E+00  2.6392E+00  1.4556E+00  7.3344E+00  9.9081E-01  9.1602E-01  1.3306E-02
             1.5622E+01
 PARAMETER:  1.6653E-01 -2.5144E-01  8.8200E+00  6.1746E-01  1.0705E+00  4.7542E-01  2.0926E+00  9.0764E-02  1.2280E-02 -4.2196E+00
             2.8487E+00
 GRADIENT:   1.6795E-01  1.2964E-01 -7.3314E-04 -2.2910E-01 -1.3349E-03  1.7987E-01  3.8955E+00  2.1022E-06  2.6884E-02  7.0614E-04
             1.1989E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -345.498282404179        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1225
 NPARAMETR:  1.0672E+00  6.9383E-01  2.4473E+04  1.6801E+00  2.6381E+00  1.4546E+00  7.3731E+00  9.8845E-01  9.1713E-01  1.0000E-02
             1.5591E+01
 PARAMETER:  1.6500E-01 -2.6552E-01  1.0205E+01  6.1885E-01  1.0701E+00  4.7472E-01  2.0978E+00  8.8384E-02  1.3491E-02 -4.5580E+00
             2.8467E+00
 GRADIENT:   4.2495E+00  8.2830E-01 -1.7283E-04  9.6265E+00  4.7433E-01  2.8496E+00  5.6787E+01  1.4511E-05 -8.4105E-04  0.0000E+00
             3.3842E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -345.499111284290        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1412
 NPARAMETR:  1.0671E+00  6.9073E-01  2.6640E+04  1.6824E+00  2.6361E+00  1.4544E+00  7.3796E+00  9.8752E-01  9.1986E-01  1.0000E-02
             1.5586E+01
 PARAMETER:  1.6494E-01 -2.7000E-01  1.0290E+01  6.2022E-01  1.0693E+00  4.7458E-01  2.0987E+00  8.7442E-02  1.6471E-02 -4.5580E+00
             2.8464E+00
 GRADIENT:  -4.4748E-01  1.1024E-01 -1.6791E-04  2.4132E-01 -1.0104E-01  1.1073E-01  4.2160E+00  3.5880E-07 -1.5959E-01  0.0000E+00
            -6.8160E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -345.499994617686        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1591
 NPARAMETR:  1.0668E+00  6.8500E-01  2.8187E+04  1.6851E+00  2.6362E+00  1.4548E+00  7.3881E+00  1.0898E+00  9.2258E-01  1.0000E-02
             1.5590E+01
 PARAMETER:  1.6470E-01 -2.7833E-01  1.0347E+01  6.2181E-01  1.0693E+00  4.7485E-01  2.0999E+00  1.8597E-01  1.9417E-02 -4.5580E+00
             2.8466E+00
 GRADIENT:  -8.0863E-01  1.2245E-02 -1.5867E-04  2.3491E-01 -8.3376E-02  1.9490E-01  4.0595E+00 -1.6415E-06 -1.7528E-01  0.0000E+00
            -2.4422E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -345.500791968550        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1782             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0671E+00  6.8249E-01  3.7215E+04  1.6861E+00  2.6386E+00  1.4543E+00  7.4021E+00  1.0934E+00  9.2519E-01  1.0000E-02
             1.5591E+01
 PARAMETER:  1.6497E-01 -2.8201E-01  1.0624E+01  6.2239E-01  1.0703E+00  4.7452E-01  2.1018E+00  1.8927E-01  2.2245E-02 -4.5580E+00
             2.8467E+00
 GRADIENT:   4.2117E+00  7.2238E-01 -1.1619E-04  9.3849E+00  5.3251E-01  2.8703E+00  5.7350E+01  1.1868E-05  1.1406E-01  0.0000E+00
             3.4060E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -345.501265105402        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1964             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0672E+00  6.8103E-01  3.2564E+04  1.6882E+00  2.6375E+00  1.4542E+00  7.4061E+00  1.0923E+00  9.2757E-01  1.0171E-02
             1.5591E+01
 PARAMETER:  1.6505E-01 -2.8414E-01  1.0491E+01  6.2369E-01  1.0698E+00  4.7448E-01  2.1023E+00  1.8825E-01  2.4814E-02 -4.4882E+00
             2.8467E+00
 GRADIENT:   4.2721E+00  7.9862E-01 -1.3154E-04  9.7980E+00  4.7305E-01  2.8746E+00  5.7310E+01  1.5339E-05  5.0729E-02  4.4204E-04
             3.3748E+01

0ITERATION NO.:   69    OBJECTIVE VALUE:  -345.501378935935        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:     2090
 NPARAMETR:  1.0673E+00  6.8042E-01  3.9561E+04  1.6883E+00  2.6387E+00  1.4541E+00  7.4127E+00  1.0916E+00  9.2799E-01  1.0000E-02
             1.5592E+01
 PARAMETER:  1.6509E-01 -2.8504E-01  1.0686E+01  6.2373E-01  1.0703E+00  4.7438E-01  2.1032E+00  1.8766E-01  2.5270E-02 -4.5293E+00
             2.8467E+00
 GRADIENT:  -5.3356E-01  3.1701E-03 -1.1775E-04 -2.3315E-01  2.0068E-02  1.2319E-01  4.4273E+00  4.4374E-06 -2.0098E-02  9.6090E-05
            -5.2157E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2090
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8715E-02  4.7368E-02  2.4138E-07 -7.1048E-02  9.6962E-06
 SE:             2.6361E-02  2.2134E-02  7.2603E-08  1.2222E-02  3.0262E-05
 N:                     100         100         100         100         100

 P VAL.:         2.7603E-01  3.2349E-02  8.8537E-04  6.1539E-09  7.4866E-01

 ETASHRINKSD(%)  1.1686E+01  2.5848E+01  1.0000E+02  5.9054E+01  9.9899E+01
 ETASHRINKVR(%)  2.2007E+01  4.5015E+01  1.0000E+02  8.3234E+01  1.0000E+02
 EBVSHRINKSD(%)  1.6460E+01  2.1841E+01  9.9999E+01  5.9535E+01  9.9826E+01
 EBVSHRINKVR(%)  3.0210E+01  3.8911E+01  1.0000E+02  8.3626E+01  1.0000E+02
 RELATIVEINF(%)  6.5362E+01  3.4004E+01  4.1702E-10  7.0996E+00  4.4506E-05
 EPSSHRINKSD(%)  3.1147E+00
 EPSSHRINKVR(%)  6.1325E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -345.50137893593455     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       757.22486090967254     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    54.64
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -345.501       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.07E+00  6.80E-01  3.96E+04  1.69E+00  2.64E+00  1.45E+00  7.41E+00  1.09E+00  9.28E-01  1.00E-02  1.56E+01
 


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
+        2.85E+02
 
 TH 2
+       -7.79E+00  2.88E+01
 
 TH 3
+       -1.43E-06 -3.46E-07  3.50E-13
 
 TH 4
+       -2.27E+01  2.47E+01 -1.00E-07  1.03E+02
 
 TH 5
+        5.41E-01 -2.24E+00  7.73E-08 -6.98E+00  4.78E+00
 
 TH 6
+       -1.70E+01 -3.89E+00 -1.41E-07  3.48E-01 -4.31E-01  5.04E+01
 
 TH 7
+        1.46E+00  3.25E+00 -2.10E-08 -4.60E+00  1.33E-01 -2.59E-01  1.76E+00
 
 TH 8
+        3.00E-01 -1.07E+00 -3.06E-07  1.53E-01  1.87E-01 -1.08E-01 -4.47E-02  2.98E-01
 
 TH 9
+        1.01E+01 -2.24E+00  4.99E-07 -3.64E+01  2.79E+00 -2.94E+00  1.32E+00 -5.92E+00  3.52E+01
 
 TH10
+       -1.61E+00  1.72E+01 -6.70E-06  1.70E+00  1.25E-01  5.72E+00  5.48E-04 -1.60E+01  2.37E+01  1.74E+01
 
 TH11
+       -8.70E+00 -2.00E+00  1.36E-09 -6.45E+00  3.86E-01  1.53E+00  1.03E-01  1.68E-02  2.32E+00 -3.06E-02  2.46E+00
 
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
 #CPUT: Total CPU Time in Seconds,       67.635
Stop Time:
Thu Sep 30 10:15:03 CDT 2021
