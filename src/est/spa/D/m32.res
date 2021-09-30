Wed Sep 29 19:57:27 CDT 2021
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
$DATA ../../../../data/spa/D/dat32.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m32.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13731.2989880271        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0044E+02  1.4634E+02 -7.0766E+01 -1.0229E+01  2.0063E+02 -1.2313E+03 -5.8326E+02 -2.7017E+01 -1.1758E+03 -2.8521E+02
            -2.7174E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -665.582587131154        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.5214E+00  1.3470E+00  1.0182E+00  1.8633E+00  1.0464E+00  1.4286E+00  1.1361E+00  9.7777E-01  1.2769E+00  1.0080E+00
             1.4816E+01
 PARAMETER:  5.1964E-01  3.9785E-01  1.1801E-01  7.2233E-01  1.4537E-01  4.5671E-01  2.2757E-01  7.7521E-02  3.4445E-01  1.0800E-01
             2.7957E+00
 GRADIENT:   4.7951E+00  5.5975E+01  2.9741E+00  1.0091E+02 -2.2815E+01 -6.9799E+00 -4.2013E+00  3.6296E+00  7.4187E-01  3.5315E+00
             9.8980E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -690.017624646314        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3698E+00  8.6665E-01  1.5993E+00  1.7365E+00  1.2404E+00  1.3955E+00  2.2628E+00  3.1893E-01  1.0430E+00  9.7378E-01
             1.4065E+01
 PARAMETER:  4.1466E-01 -4.3125E-02  5.6959E-01  6.5185E-01  3.1542E-01  4.3323E-01  9.1660E-01 -1.0428E+00  1.4209E-01  7.3428E-02
             2.7437E+00
 GRADIENT:  -1.7320E+01  1.6392E+01  7.8115E+00  2.6975E+01 -1.2219E+01  7.9063E+00  7.2804E+00  1.9436E-01  1.0813E+01  2.3369E+00
             1.5079E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -712.112568415701        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1739E+00  6.7042E-01  1.5386E+00  1.4698E+00  2.9469E+00  1.4427E+00  1.6104E+00  5.5833E-02  8.1372E-01  4.9257E+00
             1.1293E+01
 PARAMETER:  2.6037E-01 -2.9985E-01  5.3088E-01  4.8511E-01  1.1808E+00  4.6650E-01  5.7647E-01 -2.7854E+00 -1.0614E-01  1.6945E+00
             2.5241E+00
 GRADIENT:  -1.8198E+01 -2.1627E+00  9.2081E+00 -2.9191E+01 -9.9172E+00  3.5217E+01  1.5050E+00  5.3080E-04  4.5191E+00 -8.4677E-01
             5.0267E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -726.061211021106        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.1030E+00  4.7457E-01  6.9406E-01  1.3880E+00  5.8230E+00  1.1136E+00  8.3061E-01  1.0000E-02  3.0466E-01  6.6324E+00
             1.0333E+01
 PARAMETER:  1.9806E-01 -6.4535E-01 -2.6519E-01  4.2786E-01  1.8618E+00  2.0759E-01 -8.5601E-02 -6.2014E+00 -1.0886E+00  1.9920E+00
             2.4353E+00
 GRADIENT:  -5.1657E+00  1.8673E+01 -3.1078E-01  1.4611E+01 -5.6301E+00 -3.0971E+01  3.8648E-01  0.0000E+00  1.2443E+00 -3.9000E+00
            -6.5874E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -744.920973097490        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.5816E-01  8.6932E-02  2.4594E-01  1.1873E+00  8.2064E+00  1.2723E+00  2.5505E-01  1.0000E-02  2.0767E-02  6.1345E+00
             1.0802E+01
 PARAMETER:  5.7255E-02 -2.3426E+00 -1.3026E+00  2.7168E-01  2.2049E+00  3.4083E-01 -1.2663E+00 -1.3188E+01 -3.7744E+00  1.9139E+00
             2.4797E+00
 GRADIENT:   1.0282E+01  3.2625E+00 -4.5596E+01  1.2589E+02 -1.4419E+01 -1.2915E+01  6.2436E-03  0.0000E+00  3.3957E-04  5.1791E+01
            -2.3274E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -757.779627826674        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  6.7301E-01  1.0000E-02  9.3932E-02  7.8270E-01  5.7459E+00  1.3619E+00  8.0594E-02  1.0000E-02  1.0000E-02  3.5969E+00
             1.0777E+01
 PARAMETER: -2.9599E-01 -4.6094E+00 -2.2652E+00 -1.4501E-01  1.8485E+00  4.0892E-01 -2.4183E+00 -1.9322E+01 -6.7251E+00  1.3801E+00
             2.4774E+00
 GRADIENT:  -1.0137E+02  0.0000E+00 -8.6128E+01  2.2891E+02 -4.9490E+01 -1.2563E+01  3.3242E-05  0.0000E+00  0.0000E+00  6.7908E+01
             1.8156E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -764.281598529848        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      520
 NPARAMETR:  5.0522E-01  1.0000E-02  4.1921E-02  4.8697E-01  5.8956E+00  1.4260E+00  2.3784E-02  1.0000E-02  1.0000E-02  2.2424E+00
             1.0308E+01
 PARAMETER: -5.8277E-01 -6.5644E+00 -3.0720E+00 -6.1955E-01  1.8742E+00  4.5491E-01 -3.6387E+00 -2.5042E+01 -9.5158E+00  9.0757E-01
             2.4329E+00
 GRADIENT:  -9.4837E+01  0.0000E+00 -1.5085E+02  3.2245E+02 -6.7956E+01  4.3980E+00  2.8899E-05  0.0000E+00  0.0000E+00  6.7657E+01
            -8.5333E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -764.313518296457        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      603
 NPARAMETR:  4.8592E-01  1.0000E-02  3.7858E-02  4.5593E-01  5.9223E+00  1.4315E+00  2.0021E-02  1.0000E-02  1.0000E-02  2.1072E+00
             1.0265E+01
 PARAMETER: -6.2170E-01 -6.8355E+00 -3.1739E+00 -6.8542E-01  1.8787E+00  4.5876E-01 -3.8110E+00 -2.5820E+01 -9.9068E+00  8.4534E-01
             2.4287E+00
 GRADIENT:  -9.3641E+01  0.0000E+00 -1.5491E+02  3.2907E+02 -6.7793E+01  6.4139E+00  2.7211E-05  0.0000E+00  0.0000E+00  6.5820E+01
            -1.0669E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -815.530678903942        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      784
 NPARAMETR:  5.4897E-01  1.0000E-02  3.5334E-02  3.8400E-01  7.1265E+00  1.3892E+00  1.4033E-02  1.0000E-02  1.0000E-02  1.3106E+00
             1.0081E+01
 PARAMETER: -4.9972E-01 -7.9895E+00 -3.2429E+00 -8.5712E-01  2.0638E+00  4.2873E-01 -4.1664E+00 -2.7800E+01 -1.1383E+01  3.7048E-01
             2.4106E+00
 GRADIENT:   3.0933E+01  0.0000E+00 -3.2507E+01  3.4147E+01  5.2442E-01  1.6807E+01 -8.8327E-08  0.0000E+00  0.0000E+00 -1.2434E+00
            -2.8039E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -817.315946203025        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      959
 NPARAMETR:  5.1607E-01  1.0000E-02  3.2051E-02  3.5434E-01  6.9605E+00  1.3220E+00  1.0878E-02  1.0000E-02  1.0000E-02  1.1527E+00
             1.0296E+01
 PARAMETER: -5.6151E-01 -8.4557E+00 -3.3404E+00 -9.3749E-01  2.0403E+00  3.7917E-01 -4.4210E+00 -2.8862E+01 -1.2043E+01  2.4213E-01
             2.4318E+00
 GRADIENT:  -7.7623E-01  0.0000E+00 -3.5956E+00  9.5138E+00  3.0564E-01  4.7126E+00 -8.6210E-07  0.0000E+00  0.0000E+00 -3.8789E+00
             1.0635E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -819.790746559007        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1136
 NPARAMETR:  4.5258E-01  1.0000E-02  2.3139E-02  2.7295E-01  8.4258E+00  1.2787E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7530E+00
             1.0096E+01
 PARAMETER: -6.9279E-01 -8.2217E+00 -3.6662E+00 -1.1985E+00  2.2313E+00  3.4585E-01 -5.3389E+00 -3.0629E+01 -1.2521E+01  6.6135E-01
             2.4122E+00
 GRADIENT:  -8.3796E-01  0.0000E+00  8.4492E+00 -1.1646E+01  3.1322E+00 -3.1932E-01  0.0000E+00  0.0000E+00  0.0000E+00 -3.8818E+00
            -9.4493E+00

0ITERATION NO.:   59    OBJECTIVE VALUE:  -820.678288669136        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1278
 NPARAMETR:  4.7747E-01  1.0000E-02  2.6190E-02  3.0057E-01  1.0569E+01  1.2913E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.0734E+00
             1.0280E+01
 PARAMETER: -6.3887E-01 -6.4415E+00 -3.5440E+00 -1.1024E+00  2.4560E+00  3.5648E-01 -5.4120E+00 -2.8663E+01 -1.0900E+01  1.4896E+00
             2.4311E+00
 GRADIENT:   2.8094E+00  0.0000E+00 -4.0921E+01 -6.1914E+00 -6.7472E+01  8.3609E-01  0.0000E+00  0.0000E+00  0.0000E+00 -2.7520E-01
             4.7034E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1278
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.7347E-04  3.2810E-06  7.8275E-05 -1.7818E-04 -4.0110E-03
 SE:             2.8498E-02  1.4475E-06  2.9864E-04  3.7455E-04  3.2651E-03
 N:                     100         100         100         100         100

 P VAL.:         9.8954E-01  2.3415E-02  7.9324E-01  6.3427E-01  2.1928E-01

 ETASHRINKSD(%)  4.5282E+00  9.9995E+01  9.9000E+01  9.8745E+01  8.9062E+01
 ETASHRINKVR(%)  8.8514E+00  1.0000E+02  9.9990E+01  9.9984E+01  9.8804E+01
 EBVSHRINKSD(%)  4.5267E+00  9.9995E+01  9.8949E+01  9.8681E+01  8.9889E+01
 EBVSHRINKVR(%)  8.8484E+00  1.0000E+02  9.9989E+01  9.9983E+01  9.8978E+01
 RELATIVEINF(%)  3.3051E+00  2.0840E-08  6.2904E-05  1.0058E-04  3.9427E-01
 EPSSHRINKSD(%)  7.6641E+00
 EPSSHRINKVR(%)  1.4741E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -820.67828866913601     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -85.527462105397831     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -820.678       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.78E-01  1.00E-02  2.61E-02  3.00E-01  1.05E+01  1.29E+00  1.00E-02  1.00E-02  1.00E-02  4.01E+00  1.03E+01
 


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
+        2.61E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.14E+04  0.00E+00  1.90E+06
 
 TH 4
+       -6.97E+02  0.00E+00 -1.20E+05  3.41E+04
 
 TH 5
+        2.12E+00  0.00E+00  3.47E+03  1.30E+01  2.46E+01
 
 TH 6
+       -1.02E+01  0.00E+00  3.00E+02 -7.64E+01 -4.77E-01  9.57E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        1.31E-02  0.00E+00 -1.45E+02  1.09E+03 -6.27E-01 -5.82E-02  0.00E+00  0.00E+00  0.00E+00  2.18E-01
 
 TH11
+       -2.73E+01  0.00E+00  3.75E+02 -2.84E+02 -3.68E-02  9.52E-01  0.00E+00  0.00E+00  0.00E+00 -1.80E-02  4.14E+00
 
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
 #CPUT: Total CPU Time in Seconds,       23.076
Stop Time:
Wed Sep 29 19:57:51 CDT 2021
