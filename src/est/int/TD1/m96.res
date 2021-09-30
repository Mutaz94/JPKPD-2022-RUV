Wed Sep 29 06:50:46 CDT 2021
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
$DATA ../../../../data/int/TD1/dat96.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3449.17796218511        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0404E+02  5.8970E+01  8.1367E+01  9.5182E+01  9.4937E+01  2.7117E+01  2.3415E+01 -4.5430E+02 -1.0978E+02  1.7076E+01
            -2.1143E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3570.39164143800        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0565E+00  9.3432E-01  1.0354E+00  9.0085E-01  1.0028E+00  9.5321E-01  9.5825E-01  1.6732E+00  9.3770E-01  9.7223E-01
             1.6650E+00
 PARAMETER:  1.5494E-01  3.2064E-02  1.3477E-01 -4.4130E-03  1.0284E-01  5.2079E-02  5.7355E-02  6.1477E-01  3.5677E-02  7.1841E-02
             6.0983E-01
 GRADIENT:   2.0211E+02 -1.9220E+02 -4.0329E+01 -2.0855E+02  2.2582E+01 -4.7647E+01  3.0472E-02  4.3856E+01 -1.9014E+01  1.6069E+01
             6.4327E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3607.37140832846        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      334
 NPARAMETR:  9.4207E-01  1.3050E+00  2.4803E+00  7.9115E-01  1.4339E+00  9.7880E-01  9.5666E-01  3.1661E+00  4.2358E-01  1.0750E+00
             1.5453E+00
 PARAMETER:  4.0320E-02  3.6620E-01  1.0084E+00 -1.3427E-01  4.6042E-01  7.8573E-02  5.5697E-02  1.2525E+00 -7.5901E-01  1.7236E-01
             5.3519E-01
 GRADIENT:  -7.1226E+01 -6.0259E+01  3.2010E+01 -9.1064E+01 -6.1398E+01 -2.4067E+01 -5.5825E+01  2.5804E+01 -2.0248E+01 -3.2957E+01
             5.1390E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3617.74323505030        NO. OF FUNC. EVALS.: 152
 CUMULATIVE NO. OF FUNC. EVALS.:      486
 NPARAMETR:  9.4132E-01  1.3022E+00  2.4674E+00  8.3072E-01  1.4328E+00  1.0211E+00  1.1604E+00  3.1685E+00  4.2244E-01  1.2534E+00
             1.5453E+00
 PARAMETER:  3.9529E-02  3.6406E-01  1.0032E+00 -8.5458E-02  4.5963E-01  1.2091E-01  2.4873E-01  1.2533E+00 -7.6170E-01  3.2582E-01
             5.3524E-01
 GRADIENT:   1.2026E+02  1.6989E+02  3.1981E+01  3.2163E+01  2.9930E+01  2.2029E+01  8.7644E+00  4.4465E+01 -4.0895E+00  1.1941E+01
             5.3307E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3639.95929510375        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      669
 NPARAMETR:  9.4205E-01  1.2992E+00  2.4284E+00  8.7478E-01  1.4427E+00  1.0940E+00  8.9646E-01  3.0515E+00  1.2363E+00  1.1172E+00
             1.4813E+00
 PARAMETER:  4.0305E-02  3.6172E-01  9.8724E-01 -3.3785E-02  4.6650E-01  1.8987E-01 -9.2964E-03  1.2156E+00  3.1215E-01  2.1084E-01
             4.9291E-01
 GRADIENT:  -5.5996E+01  4.4434E+01  1.8075E+01  4.6126E+01  1.4143E+01  1.9012E+01  2.0885E+01  3.9642E+01  4.6625E+01  1.1243E+00
             4.8192E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3663.00890252516        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      846
 NPARAMETR:  9.4381E-01  1.2913E+00  2.3432E+00  7.8657E-01  1.4604E+00  1.0600E+00  7.5522E-02  2.8129E+00  1.7379E+00  1.4084E+00
             1.3671E+00
 PARAMETER:  4.2173E-02  3.5562E-01  9.5151E-01 -1.4007E-01  4.7872E-01  1.5822E-01 -2.4833E+00  1.1342E+00  6.5267E-01  4.4243E-01
             4.1268E-01
 GRADIENT:  -5.3594E+01 -2.8120E+01  1.9228E+01 -1.4254E+01  7.8238E+00  7.6706E+00 -1.1323E-01  1.6947E+01  6.1406E+01  3.9262E+01
             3.5399E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3708.21137609537        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1023
 NPARAMETR:  9.4786E-01  1.2734E+00  2.1586E+00  8.1240E-01  1.4992E+00  1.0378E+00  1.0000E-02  2.3351E+00  1.4136E+00  1.2243E+00
             1.1394E+00
 PARAMETER:  4.6453E-02  3.4167E-01  8.6947E-01 -1.0776E-01  5.0495E-01  1.3710E-01 -8.2669E+00  9.4807E-01  4.4613E-01  3.0235E-01
             2.3046E-01
 GRADIENT:  -4.4453E+01  8.6216E+01  1.6141E+01  3.9765E-01  7.7513E+01 -7.0789E-01  0.0000E+00 -8.4535E+00 -1.1410E+00 -1.0115E+00
             2.9342E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3714.90479326322        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1121
 NPARAMETR:  9.3353E-01  1.1811E+00  1.8150E+00  8.1765E-01  1.2960E+00  1.0223E+00  1.0000E-02  2.2550E+00  1.3829E+00  1.2020E+00
             1.1370E+00
 PARAMETER:  3.1222E-02  2.6644E-01  6.9610E-01 -1.0132E-01  3.5924E-01  1.2201E-01 -8.2669E+00  9.1315E-01  4.2418E-01  2.8397E-01
             2.2842E-01
 GRADIENT:   2.6917E+02  1.6740E+02  3.0628E+01 -1.5456E+01  1.1713E+02  4.6258E+01  0.0000E+00  1.1112E+01  5.4876E+01  2.2605E+01
             2.9307E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3740.40897995645        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  9.7676E-01  1.1324E+00  1.3477E+00  8.9909E-01  1.1722E+00  1.0123E+00  1.0000E-02  1.5394E+00  1.2251E+00  1.1372E+00
             1.0960E+00
 PARAMETER:  7.6481E-02  2.2438E-01  3.9842E-01 -6.3720E-03  2.5891E-01  1.1218E-01 -8.2669E+00  5.3142E-01  3.0305E-01  2.2861E-01
             1.9163E-01
 GRADIENT:   1.9184E+01 -6.0829E+00  3.6025E+01  8.8016E+00  3.1665E+01 -9.3969E+00  0.0000E+00 -9.2169E+00 -5.7779E+00 -2.6903E+01
            -1.9285E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3746.38831522481        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1426
 NPARAMETR:  9.7559E-01  1.1330E+00  1.1066E+00  9.0253E-01  1.1170E+00  1.0166E+00  1.0000E-02  1.3247E+00  1.2111E+00  1.2166E+00
             1.1008E+00
 PARAMETER:  7.5288E-02  2.2486E-01  2.0131E-01 -2.5533E-03  2.1061E-01  1.1648E-01 -8.2669E+00  3.8121E-01  2.9156E-01  2.9602E-01
             1.9603E-01
 GRADIENT:   1.6416E+01  9.3296E+00 -1.2229E+00  2.6715E+01  2.2291E+01 -7.5269E+00  0.0000E+00  5.5798E+00 -1.1446E+01 -1.0642E+01
            -1.0613E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3747.14989636089        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1565
 NPARAMETR:  9.7608E-01  1.1343E+00  1.1072E+00  8.9791E-01  1.1158E+00  1.0158E+00  1.0000E-02  1.0750E+00  1.2384E+00  1.2544E+00
             1.1103E+00
 PARAMETER:  7.5792E-02  2.2600E-01  2.0185E-01 -7.6816E-03  2.0955E-01  1.1569E-01 -8.2669E+00  1.7229E-01  3.1381E-01  3.2663E-01
             2.0462E-01
 GRADIENT:   3.7798E+02  1.7837E+02  2.8323E+01  7.8949E+01  1.0338E+02  3.8960E+01  0.0000E+00 -2.8855E+00  2.2644E+01  1.6796E+01
            -1.7497E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3748.19996533141        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1743            RESET HESSIAN, TYPE II
 NPARAMETR:  9.6874E-01  1.1380E+00  1.0717E+00  8.8469E-01  1.0980E+00  1.0361E+00  1.0000E-02  1.1596E+00  1.2483E+00  1.2601E+00
             1.1100E+00
 PARAMETER:  6.8236E-02  2.2923E-01  1.6928E-01 -2.2518E-02  1.9349E-01  1.3549E-01 -8.2669E+00  2.4804E-01  3.2179E-01  3.3122E-01
             2.0440E-01
 GRADIENT:   3.6189E+02  1.8868E+02  1.2984E+01  6.2640E+01  8.2673E+01  5.6268E+01  0.0000E+00  4.1866E+00  2.4987E+01  1.8191E+01
             5.5688E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3748.51948849008        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1906
 NPARAMETR:  9.6787E-01  1.1389E+00  1.0362E+00  8.7885E-01  1.0842E+00  1.0350E+00  1.0000E-02  1.0551E+00  1.2547E+00  1.2644E+00
             1.1111E+00
 PARAMETER:  6.7337E-02  2.3007E-01  1.3556E-01 -2.9143E-02  1.8087E-01  1.3439E-01 -8.2669E+00  1.5364E-01  3.2691E-01  3.3459E-01
             2.0535E-01
 GRADIENT:  -3.4308E-01  2.0449E+00  1.3311E-01 -5.4093E-01  4.0771E-02 -7.7079E-02  0.0000E+00  2.4561E-02  1.2177E-01  4.8843E-02
             3.4684E-01

0ITERATION NO.:   61    OBJECTIVE VALUE:  -3748.51948849008        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1928
 NPARAMETR:  9.6787E-01  1.1389E+00  1.0362E+00  8.7885E-01  1.0842E+00  1.0350E+00  1.0000E-02  1.0551E+00  1.2547E+00  1.2644E+00
             1.1111E+00
 PARAMETER:  6.7337E-02  2.3007E-01  1.3556E-01 -2.9143E-02  1.8087E-01  1.3439E-01 -8.2669E+00  1.5364E-01  3.2691E-01  3.3459E-01
             2.0535E-01
 GRADIENT:  -3.4308E-01  2.0449E+00  1.3311E-01 -5.4093E-01  4.0771E-02 -7.7079E-02  0.0000E+00  2.4561E-02  1.2177E-01  4.8843E-02
             3.4684E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1928
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4570E-03 -1.1886E-03 -1.9780E-02  6.9291E-04 -5.7262E-03
 SE:             2.9928E-02  5.2590E-04  1.5845E-02  2.9590E-02  2.8870E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6117E-01  2.3808E-02  2.1191E-01  9.8132E-01  8.4278E-01

 ETASHRINKSD(%)  1.0000E-10  9.8238E+01  4.6916E+01  8.6902E-01  3.2806E+00
 ETASHRINKVR(%)  1.0000E-10  9.9969E+01  7.1821E+01  1.7305E+00  6.4536E+00
 EBVSHRINKSD(%)  2.8733E-01  9.8826E+01  4.8705E+01  1.2067E+00  3.2973E+00
 EBVSHRINKVR(%)  5.7384E-01  9.9986E+01  7.3688E+01  2.3988E+00  6.4859E+00
 RELATIVEINF(%)  9.9422E+01  5.2240E-03  2.1513E+01  5.4827E+01  4.5107E+01
 EPSSHRINKSD(%)  1.9375E+01
 EPSSHRINKVR(%)  3.4996E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3748.5194884900843     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2094.4301287216736     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    54.57
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3748.519       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  1.14E+00  1.04E+00  8.79E-01  1.08E+00  1.03E+00  1.00E-02  1.06E+00  1.25E+00  1.26E+00  1.11E+00
 


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
+        1.10E+03
 
 TH 2
+       -6.41E+00  1.24E+03
 
 TH 3
+       -1.26E-01  9.74E+00  2.30E+02
 
 TH 4
+       -4.28E+00  5.72E+02 -2.41E+01  8.92E+02
 
 TH 5
+       -1.51E-01 -4.34E+02 -1.32E+02 -7.53E+00  5.52E+02
 
 TH 6
+        2.85E+00 -2.64E+00  1.64E-01 -1.37E+00 -5.49E-02  1.84E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -5.29E-02  9.86E-01 -2.11E+01 -2.24E+00 -3.39E+00 -5.23E-02  0.00E+00  2.16E+01
 
 TH 9
+        1.12E+00 -1.26E+02  3.05E+00  4.20E+00  7.70E-01  9.99E-02  0.00E+00 -1.78E-01  1.19E+02
 
 TH10
+        1.14E-01 -2.76E+01  5.65E+00  1.75E+00 -2.58E+00  7.42E-02  0.00E+00 -1.08E+00  4.96E-01  1.07E+02
 
 TH11
+       -6.87E+00 -4.40E+01 -2.71E+01 -9.73E+00  5.58E-01  1.41E+00  0.00E+00  2.91E+01  5.36E+00  1.02E+01  8.96E+02
 
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
 #CPUT: Total CPU Time in Seconds,       66.918
Stop Time:
Wed Sep 29 06:51:55 CDT 2021
