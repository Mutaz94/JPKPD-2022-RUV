Sat Sep 18 11:24:15 CDT 2021
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
$DATA ../../../../data/spa/S1/dat96.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1639.44483334136        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2254E+02 -7.7800E+00 -3.4213E+01  3.5968E+01  2.8954E+01 -9.7554E+00  8.8363E+00  1.0099E+01  2.4789E+01  1.3788E+01
            -4.1147E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1645.35739368741        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9095E-01  1.0481E+00  1.1108E+00  9.7688E-01  1.0483E+00  1.0264E+00  9.3104E-01  9.2210E-01  8.5636E-01  9.0111E-01
             1.1358E+00
 PARAMETER:  9.0908E-02  1.4700E-01  2.0507E-01  7.6609E-02  1.4720E-01  1.2602E-01  2.8548E-02  1.8902E-02 -5.5069E-02 -4.1287E-03
             2.2737E-01
 GRADIENT:   8.8662E+01  2.0071E+01 -1.8668E+00  2.9203E+01  1.9751E+01  2.9824E+00 -3.0960E+00  6.3628E-01 -9.3474E+00 -9.6900E+00
             6.4713E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1646.67105220913        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.8383E-01  1.0241E+00  1.0626E+00  9.8297E-01  1.0301E+00  1.0319E+00  9.6018E-01  5.8419E-01  8.8056E-01  9.7116E-01
             1.1293E+00
 PARAMETER:  8.3702E-02  1.2384E-01  1.6073E-01  8.2828E-02  1.2962E-01  1.3141E-01  5.9363E-02 -4.3753E-01 -2.7199E-02  7.0738E-02
             2.2163E-01
 GRADIENT:   7.3608E+01  8.9956E+00 -8.0576E-01  1.7817E+01  1.8958E+01  6.0563E+00  3.3525E-01 -7.4409E-01 -2.6691E+00 -1.7609E+00
             5.8354E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1647.76643122027        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.6077E-01  1.0469E+00  8.6404E-01  9.5065E-01  9.3873E-01  1.0158E+00  9.7472E-01  4.1316E-01  8.9801E-01  8.9030E-01
             1.1084E+00
 PARAMETER:  5.9979E-02  1.4579E-01 -4.6139E-02  4.9393E-02  3.6775E-02  1.1572E-01  7.4395E-02 -7.8391E-01 -7.5732E-03 -1.6196E-02
             2.0292E-01
 GRADIENT:   2.2574E+01 -2.1019E-01 -7.7012E+00  7.4993E+00  1.0782E+01 -1.7031E-01  7.0720E-01  6.9246E-01  1.8671E+00  1.5300E+00
             3.2068E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1647.76836567459        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.5795E-01  1.0343E+00  8.6522E-01  9.5745E-01  9.3200E-01  1.0158E+00  9.8446E-01  3.9929E-01  8.9055E-01  8.8650E-01
             1.1065E+00
 PARAMETER:  5.7044E-02  1.3377E-01 -4.4773E-02  5.6521E-02  2.9577E-02  1.1568E-01  8.4337E-02 -8.1806E-01 -1.5921E-02 -2.0477E-02
             2.0124E-01
 GRADIENT:   1.6552E+01 -3.2766E-01 -5.9912E+00  5.4776E+00  8.1351E+00 -2.0540E-01  5.9085E-01  6.0463E-01  1.4505E+00  1.3337E+00
             2.4082E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1647.76867078419        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.5695E-01  1.0296E+00  8.6412E-01  9.5994E-01  9.2867E-01  1.0159E+00  9.8849E-01  3.8793E-01  8.8759E-01  8.8436E-01
             1.1059E+00
 PARAMETER:  5.5994E-02  1.2919E-01 -4.6048E-02  5.9111E-02  2.6001E-02  1.1575E-01  8.8420E-02 -8.4694E-01 -1.9247E-02 -2.2894E-02
             2.0064E-01
 GRADIENT:   1.4372E+01 -3.3513E-01 -5.2954E+00  4.7484E+00  7.1219E+00 -1.9689E-01  5.3218E-01  5.5138E-01  1.2811E+00  1.2149E+00
             2.1087E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1647.76890410338        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  9.5630E-01  1.0266E+00  8.6278E-01  9.6149E-01  9.2621E-01  1.0159E+00  9.9122E-01  3.7803E-01  8.8562E-01  8.8271E-01
             1.1055E+00
 PARAMETER:  5.5320E-02  1.2622E-01 -4.7592E-02  6.0731E-02  2.3343E-02  1.1582E-01  9.1186E-02 -8.7277E-01 -2.1470E-02 -2.4760E-02
             2.0026E-01
 GRADIENT:   1.2961E+01 -3.2752E-01 -4.8216E+00  4.2784E+00  6.4510E+00 -1.8605E-01  4.8918E-01  5.1011E-01  1.1657E+00  1.1234E+00
             1.9108E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1647.76902533476        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      523
 NPARAMETR:  9.5587E-01  1.0245E+00  8.6158E-01  9.6251E-01  9.2440E-01  1.0160E+00  9.9312E-01  3.7011E-01  8.8427E-01  8.8148E-01
             1.1052E+00
 PARAMETER:  5.4870E-02  1.2424E-01 -4.8982E-02  6.1785E-02  2.1394E-02  1.1588E-01  9.3096E-02 -8.9395E-01 -2.2989E-02 -2.6151E-02
             2.0000E-01
 GRADIENT:   1.2012E+01 -3.1714E-01 -4.4931E+00  3.9631E+00  5.9935E+00 -1.7677E-01  4.5821E-01  4.7956E-01  1.0858E+00  1.0558E+00
             1.7759E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1647.76909370660        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      598
 NPARAMETR:  9.5557E-01  1.0231E+00  8.6058E-01  9.6320E-01  9.2306E-01  1.0160E+00  9.9448E-01  3.6383E-01  8.8332E-01  8.8055E-01
             1.1050E+00
 PARAMETER:  5.4553E-02  1.2285E-01 -5.0153E-02  6.2506E-02  1.9934E-02  1.1592E-01  9.4468E-02 -9.1107E-01 -2.4073E-02 -2.7205E-02
             1.9981E-01
 GRADIENT:   1.1343E+01 -3.0722E-01 -4.2566E+00  3.7409E+00  5.6677E+00 -1.6943E-01  4.3542E-01  4.5666E-01  1.0283E+00  1.0053E+00
             1.6797E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1647.76913510744        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      673
 NPARAMETR:  9.5535E-01  1.0221E+00  8.5974E-01  9.6370E-01  9.2201E-01  1.0161E+00  9.9551E-01  3.5876E-01  8.8260E-01  8.7983E-01
             1.1048E+00
 PARAMETER:  5.4319E-02  1.2183E-01 -5.1128E-02  6.3029E-02  1.8804E-02  1.1595E-01  9.5498E-02 -9.2509E-01 -2.4882E-02 -2.8026E-02
             1.9968E-01
 GRADIENT:   1.0846E+01 -2.9851E-01 -4.0786E+00  3.5763E+00  5.4244E+00 -1.6326E-01  4.1798E-01  4.3898E-01  9.8510E-01  9.6627E-01
             1.6079E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1647.76916389896        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      748
 NPARAMETR:  9.5517E-01  1.0212E+00  8.5901E-01  9.6410E-01  9.2115E-01  1.0161E+00  9.9634E-01  3.5445E-01  8.8203E-01  8.7923E-01
             1.1047E+00
 PARAMETER:  5.4133E-02  1.2102E-01 -5.1978E-02  6.3437E-02  1.7870E-02  1.1598E-01  9.6328E-02 -9.3720E-01 -2.5532E-02 -2.8707E-02
             1.9957E-01
 GRADIENT:   1.0449E+01 -2.9074E-01 -3.9348E+00  3.4447E+00  5.2290E+00 -1.5821E-01  4.0376E-01  4.2442E-01  9.5020E-01  9.3419E-01
             1.5502E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1647.76918308968        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      823
 NPARAMETR:  9.5502E-01  1.0205E+00  8.5833E-01  9.6443E-01  9.2039E-01  1.0161E+00  9.9706E-01  3.5052E-01  8.8153E-01  8.7870E-01
             1.1046E+00
 PARAMETER:  5.3972E-02  1.2033E-01 -5.2762E-02  6.3783E-02  1.7042E-02  1.1601E-01  9.7051E-02 -9.4833E-01 -2.6095E-02 -2.9313E-02
             1.9948E-01
 GRADIENT:   1.0106E+01 -2.8345E-01 -3.8094E+00  3.3310E+00  5.0594E+00 -1.5331E-01  3.9123E-01  4.1157E-01  9.1979E-01  9.0585E-01
             1.5001E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1648.01130802750        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1004
 NPARAMETR:  9.6508E-01  9.2324E-01  8.6773E-01  1.0265E+00  8.7815E-01  1.0267E+00  1.0850E+00  2.4904E-01  8.3749E-01  8.6183E-01
             1.1038E+00
 PARAMETER:  6.4453E-02  2.0130E-02 -4.1873E-02  1.2615E-01 -2.9943E-02  1.2633E-01  1.8154E-01 -1.2901E+00 -7.7346E-02 -4.8695E-02
             1.9875E-01
 GRADIENT:  -1.1028E+00  8.3873E-01 -3.2647E-01  1.5953E+00 -7.6747E-01 -2.8740E-01  4.0034E-02  1.0976E-01 -5.8361E-03  2.2914E-01
             5.1362E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1648.04221142214        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1180
 NPARAMETR:  9.6453E-01  9.4509E-01  8.5150E-01  1.0110E+00  8.7911E-01  1.0269E+00  1.0667E+00  1.4150E-01  8.4592E-01  8.6007E-01
             1.1038E+00
 PARAMETER:  6.3882E-02  4.3524E-02 -6.0760E-02  1.1090E-01 -2.8840E-02  1.2656E-01  1.6455E-01 -1.8554E+00 -6.7327E-02 -5.0742E-02
             1.9879E-01
 GRADIENT:  -2.7289E+00  4.6630E-01  1.4829E+00 -8.8323E-01 -2.3466E+00 -2.7194E-01 -6.2903E-02  1.9925E-02 -4.8989E-02 -2.3011E-01
             2.3890E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1648.06301859561        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1355
 NPARAMETR:  9.6620E-01  1.0027E+00  8.3014E-01  9.7476E-01  8.9562E-01  1.0280E+00  1.0202E+00  5.0573E-02  8.6938E-01  8.6532E-01
             1.1028E+00
 PARAMETER:  6.5618E-02  1.0274E-01 -8.6155E-02  7.4435E-02 -1.0235E-02  1.2764E-01  1.1999E-01 -2.8843E+00 -3.9971E-02 -4.4653E-02
             1.9785E-01
 GRADIENT:   1.2688E-01 -2.8549E-01 -2.3910E-01 -1.9123E-01  4.1228E-01  1.8928E-02  2.0274E-02  3.9494E-03  3.7502E-02  1.5537E-02
            -1.4916E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1648.06521058381        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1530
 NPARAMETR:  9.6613E-01  9.9316E-01  8.3180E-01  9.8084E-01  8.9183E-01  1.0280E+00  1.0281E+00  1.0566E-02  8.6502E-01  8.6402E-01
             1.1029E+00
 PARAMETER:  6.5542E-02  9.3139E-02 -8.4157E-02  8.0652E-02 -1.4475E-02  1.2758E-01  1.2770E-01 -4.4501E+00 -4.5006E-02 -4.6156E-02
             1.9790E-01
 GRADIENT:   4.5786E-02  4.9879E-02 -6.7151E-03  8.6360E-02  1.9038E-02  7.1763E-03  2.3113E-03  1.5902E-04 -4.1060E-03 -4.8771E-03
            -1.7712E-02

0ITERATION NO.:   78    OBJECTIVE VALUE:  -1648.06522166583        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1622
 NPARAMETR:  9.6611E-01  9.9306E-01  8.3179E-01  9.8085E-01  8.9177E-01  1.0279E+00  1.0281E+00  1.0000E-02  8.6500E-01  8.6400E-01
             1.1029E+00
 PARAMETER:  6.5520E-02  9.3033E-02 -8.4180E-02  8.0668E-02 -1.4545E-02  1.2756E-01  1.2772E-01 -4.7334E+00 -4.5023E-02 -4.6186E-02
             1.9794E-01
 GRADIENT:  -4.5789E-04  1.3449E-03  3.1227E-03 -1.1277E-03 -4.5664E-03 -6.2504E-04 -1.3386E-03  0.0000E+00 -1.1647E-03  2.6660E-04
             1.9486E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1622
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.9912E-04 -6.0727E-03 -3.2728E-04 -1.5599E-03 -1.8443E-02
 SE:             2.9808E-02  2.0032E-02  1.6562E-04  2.4074E-02  2.3158E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9467E-01  7.6177E-01  4.8144E-02  9.4834E-01  4.2581E-01

 ETASHRINKSD(%)  1.3985E-01  3.2891E+01  9.9445E+01  1.9348E+01  2.2417E+01
 ETASHRINKVR(%)  2.7951E-01  5.4963E+01  9.9997E+01  3.4953E+01  3.9809E+01
 EBVSHRINKSD(%)  4.8034E-01  3.2764E+01  9.9478E+01  1.9625E+01  2.1231E+01
 EBVSHRINKVR(%)  9.5837E-01  5.4793E+01  9.9997E+01  3.5399E+01  3.7955E+01
 RELATIVEINF(%)  9.8608E+01  1.5488E+00  2.5496E-04  2.6955E+00  5.1089E+00
 EPSSHRINKSD(%)  4.2086E+01
 EPSSHRINKVR(%)  6.6460E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1648.0652216658257     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -912.91439510208750     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.74
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1648.065       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.66E-01  9.93E-01  8.32E-01  9.81E-01  8.92E-01  1.03E+00  1.03E+00  1.00E-02  8.65E-01  8.64E-01  1.10E+00
 


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
+        1.11E+03
 
 TH 2
+       -1.23E+01  4.99E+02
 
 TH 3
+        1.68E+01  2.06E+02  4.65E+02
 
 TH 4
+       -1.85E+01  4.64E+02 -2.21E+02  1.00E+03
 
 TH 5
+       -4.70E+00 -3.88E+02 -6.42E+02  2.61E+02  1.16E+03
 
 TH 6
+       -1.45E+00 -2.01E+00  3.08E+00 -4.87E+00 -4.30E+00  1.84E+02
 
 TH 7
+       -6.51E-01  2.45E+01  6.13E+00 -9.87E+00 -1.22E+01  1.01E+00  4.45E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.24E+00 -2.56E+01 -2.20E+01  2.58E+01  3.93E+00  3.59E+00  2.34E+01  0.00E+00  1.21E+02
 
 TH10
+       -1.69E+00 -3.74E+00 -4.15E+01 -1.53E+01 -5.81E+01 -9.54E-01  1.38E+01  0.00E+00  4.66E+00  1.02E+02
 
 TH11
+       -6.43E+00 -1.55E+01 -3.51E+01 -4.62E+00  8.45E+00  2.49E+00  6.14E+00  0.00E+00  1.18E+01  2.55E+01  1.83E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.559
Stop Time:
Sat Sep 18 11:24:38 CDT 2021
