Thu Sep 30 09:13:30 CDT 2021
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
$DATA ../../../../data/spa2/D/dat51.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m51.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17750.0069718171        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.7552E+02  3.8539E+02 -3.8994E+01  2.9682E+02  2.5901E+02 -1.6724E+03 -7.7038E+02 -6.9843E+01 -1.2430E+03 -6.6830E+02
            -3.5493E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -633.957650487316        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3004E+00  1.1216E+00  1.0153E+00  1.3697E+00  9.9685E-01  2.0267E+00  1.5084E+00  9.9157E-01  1.5138E+00  1.1462E+00
             1.4076E+01
 PARAMETER:  3.6267E-01  2.1474E-01  1.1519E-01  4.1463E-01  9.6849E-02  8.0641E-01  5.1105E-01  9.1535E-02  5.1462E-01  2.3643E-01
             2.7445E+00
 GRADIENT:   3.3194E+00 -3.6766E+01 -2.4785E+01 -2.6270E+01  3.0434E+01  5.4690E+01 -6.1068E+00  4.5774E+00 -1.5118E-01  1.8958E+01
             3.3060E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -704.142574170908        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.3122E+00  4.1984E-01  1.4958E+01  3.2650E+00  2.2410E+00  2.3094E+00  8.6944E+00  5.3702E-01  3.1620E+00  1.9834E+00
             1.1452E+01
 PARAMETER:  3.7171E-01 -7.6787E-01  2.8052E+00  1.2833E+00  9.0690E-01  9.3697E-01  2.2627E+00 -5.2172E-01  1.2512E+00  7.8483E-01
             2.5382E+00
 GRADIENT:   3.2235E+01  8.0240E+00 -4.5374E+00  9.6608E+01 -1.5849E+01  6.1896E+01  8.7032E+00  5.7216E-02  4.6671E+01  2.7861E+01
             2.3594E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -713.403438030388        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.2709E+00  4.9860E-01  1.2462E+01  2.6583E+00  2.0873E+00  2.3181E+00  7.5716E+00  7.3915E-01  2.9644E+00  1.8715E+00
             1.1397E+01
 PARAMETER:  3.3974E-01 -5.9595E-01  2.6227E+00  1.0777E+00  8.3588E-01  9.4075E-01  2.1244E+00 -2.0226E-01  1.1867E+00  7.2674E-01
             2.5334E+00
 GRADIENT:   1.4647E+01  1.3109E+01 -6.5255E+00  5.4246E+01 -1.6484E+01  5.4725E+01  2.7882E+01  3.9121E-02  2.9237E+01  2.7911E+01
             2.5533E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -713.926573742517        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  1.2691E+00  4.9848E-01  1.2515E+01  2.6482E+00  2.0870E+00  2.3166E+00  7.5526E+00  7.6110E-01  2.9612E+00  1.8667E+00
             1.1385E+01
 PARAMETER:  3.3830E-01 -5.9620E-01  2.6269E+00  1.0739E+00  8.3572E-01  9.4011E-01  2.1219E+00 -1.7299E-01  1.1856E+00  7.2415E-01
             2.5323E+00
 GRADIENT:  -2.2264E+00  1.1475E+01 -6.6858E+00  2.2752E+01 -1.7193E+01  2.3751E+01  2.3313E+01  4.0026E-02  1.1207E+01  2.7464E+01
             2.3064E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -761.930750924542        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      497
 NPARAMETR:  1.2306E+00  4.7498E-01  1.4460E+01  1.3183E+00  2.2177E+00  2.2408E+00  7.2106E+00  7.6826E-01  3.1443E+00  1.6397E+00
             8.7347E+00
 PARAMETER:  3.0749E-01 -6.4448E-01  2.7714E+00  3.7631E-01  8.9646E-01  9.0683E-01  2.0755E+00 -1.6363E-01  1.2456E+00  5.9449E-01
             2.2673E+00
 GRADIENT:   5.0681E+01 -6.2818E+00 -9.5133E-01 -5.1944E+01  9.5811E+00  4.4784E+01  1.2427E+02  1.5173E-02  8.0894E+01  2.6652E+01
             8.1817E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -800.481955439583        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      570
 NPARAMETR:  1.0758E+00  8.6072E-01  2.6389E+01  1.4560E+00  1.9893E+00  1.9444E+00  3.9842E+00  7.7741E-01  2.0388E+00  4.2850E-01
             8.1241E+00
 PARAMETER:  1.7311E-01 -4.9991E-02  3.3730E+00  4.7567E-01  7.8781E-01  7.6497E-01  1.4823E+00 -1.5178E-01  8.1238E-01 -7.4745E-01
             2.1948E+00
 GRADIENT:  -1.0171E+01 -1.3289E+01  6.1230E-01  2.4705E+01 -1.4778E+00 -2.3949E+00 -1.8710E+01  6.0816E-03  3.1666E+01  1.8272E+00
            -1.7521E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -813.007141262142        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      643
 NPARAMETR:  1.1048E+00  1.3457E+00  2.3898E+01  9.8665E-01  2.0459E+00  1.9568E+00  3.7595E+00  7.9432E-01  1.1511E+00  1.4696E-01
             8.1592E+00
 PARAMETER:  1.9967E-01  3.9690E-01  3.2738E+00  8.6560E-02  8.1582E-01  7.7133E-01  1.4243E+00 -1.3027E-01  2.4068E-01 -1.8176E+00
             2.1991E+00
 GRADIENT:   7.5205E+00 -4.6331E+00  2.5073E-01 -1.9433E+01  5.2653E+00 -9.8709E-01  1.0494E+01  1.6511E-03  1.2025E+01  2.1944E-01
            -3.2533E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -819.170031784412        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      740
 NPARAMETR:  1.0820E+00  1.3085E+00  2.2769E+01  1.0114E+00  2.0483E+00  1.8437E+00  3.9116E+00  8.1320E-01  6.1752E-01  5.8547E-02
             8.3402E+00
 PARAMETER:  1.7879E-01  3.6886E-01  3.2254E+00  1.1137E-01  8.1699E-01  7.1176E-01  1.4639E+00 -1.0678E-01 -3.8204E-01 -2.7379E+00
             2.2211E+00
 GRADIENT:  -1.3586E+01 -6.8952E+00  2.8324E-01 -9.0662E+00 -2.9111E-01 -4.0884E+01 -3.8032E+01  2.8844E-03  2.5015E+00  3.3827E-02
            -2.1339E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -829.889399627572        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      918
 NPARAMETR:  1.1074E+00  7.3641E-01  1.8394E+01  1.3805E+00  2.0457E+00  2.0614E+00  5.8014E+00  7.5418E-01  8.1598E-01  1.0557E-01
             8.5152E+00
 PARAMETER:  2.0201E-01 -2.0597E-01  3.0120E+00  4.2244E-01  8.1576E-01  8.2337E-01  1.8581E+00 -1.8213E-01 -1.0337E-01 -2.1484E+00
             2.2418E+00
 GRADIENT:  -2.3488E+00  1.0911E+00  3.9840E-01  6.2766E+00 -1.4025E+00  1.4385E+00 -1.6729E+00  9.6544E-03 -2.7492E+00  9.8517E-02
             5.8924E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -830.293193646238        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1095
 NPARAMETR:  1.1118E+00  6.1197E-01  1.6694E+01  1.4640E+00  2.0333E+00  2.0409E+00  6.2181E+00  7.0009E-01  9.4726E-01  1.3754E-01
             8.5026E+00
 PARAMETER:  2.0602E-01 -3.9107E-01  2.9151E+00  4.8115E-01  8.0968E-01  8.1340E-01  1.9275E+00 -2.5654E-01  4.5821E-02 -1.8838E+00
             2.2404E+00
 GRADIENT:  -2.6649E-01  5.5735E-01  3.6204E-01  3.9498E+00 -1.7497E-01 -1.7592E+00  1.1524E+00  1.1880E-02 -1.3470E-01  1.6686E-01
            -1.1046E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -830.405027989255        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1274             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1122E+00  6.1551E-01  9.8613E+00  1.4501E+00  2.0325E+00  2.0659E+00  6.2659E+00  5.1968E-01  9.3940E-01  8.2070E-02
             8.5030E+00
 PARAMETER:  2.0636E-01 -3.8531E-01  2.3886E+00  4.7161E-01  8.0927E-01  8.2557E-01  1.9351E+00 -5.5455E-01  3.7485E-02 -2.4002E+00
             2.2404E+00
 GRADIENT:   1.5947E+01  3.2424E+00 -5.6990E-01  1.1338E+01  9.4641E+00  4.1836E+01  9.6163E+01  2.1433E-02  1.3603E+00  6.1624E-02
             2.8660E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -830.706827033440        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1450
 NPARAMETR:  1.1128E+00  5.8675E-01  7.2531E+00  1.4575E+00  1.8695E+00  2.0523E+00  6.2951E+00  5.8383E-02  9.5045E-01  2.0176E-02
             8.5111E+00
 PARAMETER:  2.0684E-01 -4.3315E-01  2.0814E+00  4.7670E-01  7.2565E-01  8.1897E-01  1.9398E+00 -2.7407E+00  4.9181E-02 -3.8033E+00
             2.2414E+00
 GRADIENT:   7.3196E-01 -8.4235E-03 -1.7788E-01  3.8149E-01  1.0819E+00  1.2930E-01  2.9669E+00  5.6739E-04  4.4868E-01  3.7972E-03
             1.4657E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -830.735577050242        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1628
 NPARAMETR:  1.1108E+00  5.9105E-01  7.2456E+00  1.4491E+00  1.8616E+00  2.0693E+00  6.3534E+00  3.1751E-02  9.3875E-01  1.0000E-02
             8.4982E+00
 PARAMETER:  2.0512E-01 -4.2585E-01  2.0804E+00  4.7096E-01  7.2145E-01  8.2722E-01  1.9490E+00 -3.3498E+00  3.6790E-02 -5.1359E+00
             2.2399E+00
 GRADIENT:   3.0723E-01  3.8479E-01 -5.2997E-02 -1.9691E+00  3.9681E-01  2.6974E+00  5.5479E+00  1.6585E-04  4.0452E-01  0.0000E+00
             1.1888E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -830.744741836783        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1806
 NPARAMETR:  1.1102E+00  5.7839E-01  7.2708E+00  1.4538E+00  1.8615E+00  2.0692E+00  6.3830E+00  1.0000E-02  9.4379E-01  1.0000E-02
             8.4976E+00
 PARAMETER:  2.0452E-01 -4.4750E-01  2.0839E+00  4.7418E-01  7.2136E-01  8.2717E-01  1.9536E+00 -4.7297E+00  4.2151E-02 -5.1318E+00
             2.2398E+00
 GRADIENT:   8.7137E-02  2.5375E-02 -7.1544E-02 -2.5630E+00  5.7380E-01  2.7244E+00  5.3379E+00  0.0000E+00  4.9399E-01  0.0000E+00
             1.4425E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -830.753914112180        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1995             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1102E+00  5.7003E-01  7.3016E+00  1.4594E+00  1.8599E+00  2.0696E+00  6.4220E+00  1.0000E-02  9.4629E-01  1.0000E-02
             8.4934E+00
 PARAMETER:  2.0451E-01 -4.6207E-01  2.0881E+00  4.7806E-01  7.2051E-01  8.2735E-01  1.9597E+00 -4.7226E+00  4.4789E-02 -5.1318E+00
             2.2393E+00
 GRADIENT:   1.6059E+01  3.7328E+00  4.6143E-02  1.3705E+01  1.5190E+00  4.3009E+01  1.0135E+02  0.0000E+00  5.1632E-01  0.0000E+00
             2.7367E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -830.757444345185        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2174             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1101E+00  5.6518E-01  7.3205E+00  1.4624E+00  1.8596E+00  2.0696E+00  6.4374E+00  1.0000E-02  9.4732E-01  1.0000E-02
             8.4927E+00
 PARAMETER:  2.0448E-01 -4.7061E-01  2.0907E+00  4.8010E-01  7.2038E-01  8.2735E-01  1.9621E+00 -4.7226E+00  4.5882E-02 -5.1318E+00
             2.2392E+00
 GRADIENT:   1.6057E+01  3.7375E+00  5.4313E-02  1.4240E+01  1.4399E+00  4.3034E+01  1.0173E+02  0.0000E+00  3.9729E-01  0.0000E+00
             2.7110E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -830.759742199894        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2354
 NPARAMETR:  1.1103E+00  5.6149E-01  7.3627E+00  1.4661E+00  1.8586E+00  2.0690E+00  6.4457E+00  1.0000E-02  9.4658E-01  1.0000E-02
             8.4926E+00
 PARAMETER:  2.0464E-01 -4.7717E-01  2.0964E+00  4.8260E-01  7.1980E-01  8.2708E-01  1.9634E+00 -4.7226E+00  4.5096E-02 -5.1318E+00
             2.2392E+00
 GRADIENT:   1.9387E-01  5.6321E-02  1.6281E-02 -5.1392E-02 -1.6540E-01  2.7846E+00  5.0599E+00  0.0000E+00 -1.6368E-01  0.0000E+00
            -3.3518E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -830.761546935177        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2539
 NPARAMETR:  1.1103E+00  5.5738E-01  7.3751E+00  1.4678E+00  1.8593E+00  2.0694E+00  6.4626E+00  1.0000E-02  9.4920E-01  1.0000E-02
             8.4931E+00
 PARAMETER:  2.0460E-01 -4.8451E-01  2.0981E+00  4.8377E-01  7.2020E-01  8.2728E-01  1.9660E+00 -4.7226E+00  4.7864E-02 -5.1318E+00
             2.2393E+00
 GRADIENT:   1.8342E-01 -2.3358E-02 -3.7803E-03 -5.0720E-01 -7.9336E-03  2.8529E+00  5.2318E+00  0.0000E+00 -5.0723E-02  0.0000E+00
            -6.4041E-02

0ITERATION NO.:   93    OBJECTIVE VALUE:  -830.762137526466        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:     2636
 NPARAMETR:  1.1101E+00  5.5759E-01  7.3789E+00  1.4686E+00  1.8593E+00  2.0694E+00  6.4745E+00  1.0000E-02  9.5032E-01  1.0000E-02
             8.4923E+00
 PARAMETER:  2.0448E-01 -4.8414E-01  2.0986E+00  4.8429E-01  7.2021E-01  8.2725E-01  1.9679E+00 -4.7226E+00  4.9045E-02 -5.1318E+00
             2.2392E+00
 GRADIENT:   1.3617E-01  1.0418E-01 -6.0264E-03 -4.7933E-01 -2.5070E-02  2.8251E+00  5.5882E+00  0.0000E+00 -1.3945E-02  0.0000E+00
            -1.9701E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2636
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2964E-02  5.0285E-02  3.6577E-06 -7.6351E-02  2.7014E-05
 SE:             2.8684E-02  2.2557E-02  9.9304E-06  1.3965E-02  7.2171E-05
 N:                     100         100         100         100         100

 P VAL.:         6.5129E-01  2.5795E-02  7.1262E-01  4.5771E-08  7.0818E-01

 ETASHRINKSD(%)  3.9057E+00  2.4433E+01  9.9967E+01  5.3217E+01  9.9758E+01
 ETASHRINKVR(%)  7.6589E+00  4.2896E+01  1.0000E+02  7.8113E+01  9.9999E+01
 EBVSHRINKSD(%)  4.5261E+00  2.0610E+01  9.9954E+01  5.3127E+01  9.9680E+01
 EBVSHRINKVR(%)  8.8474E+00  3.6972E+01  1.0000E+02  7.8029E+01  9.9999E+01
 RELATIVEINF(%)  9.0392E+01  3.1844E+01  2.5664E-06  1.0512E+01  1.1944E-04
 EPSSHRINKSD(%)  9.4041E+00
 EPSSHRINKVR(%)  1.7924E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -830.76213752646561     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       271.96410231914149     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    60.79
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    11.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -830.762       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.11E+00  5.58E-01  7.38E+00  1.47E+00  1.86E+00  2.07E+00  6.47E+00  1.00E-02  9.50E-01  1.00E-02  8.49E+00
 


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
+        4.94E+02
 
 TH 2
+       -1.49E+02  9.35E+01
 
 TH 3
+        2.99E-01 -3.71E-01  5.21E-03
 
 TH 4
+       -3.08E+02  1.27E+02  2.79E-01  3.60E+02
 
 TH 5
+        1.33E+01 -7.75E-01 -1.37E-01 -3.15E+01  5.30E+00
 
 TH 6
+       -7.17E+01  9.89E+00  2.55E-01  6.40E+01 -1.02E+01  4.90E+01
 
 TH 7
+        4.53E-01  5.45E+00 -8.91E-02 -8.91E+00  2.58E+00 -3.67E+00  1.73E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.45E+01 -2.10E+01 -2.26E-01 -1.02E+02  1.22E+01 -2.33E+01  5.11E+00  0.00E+00  3.32E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        6.03E+00 -4.44E+00 -2.08E-02 -1.51E+01  1.57E+00 -7.68E-01  5.32E-01  0.00E+00  4.48E+00  0.00E+00  1.27E+00
 
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
+        1.87E+02
 
 TH 2
+       -3.95E+00  5.59E+01
 
 TH 3
+        5.33E-02  1.57E-01  5.88E-02
 
 TH 4
+       -9.48E+00  4.41E+01  1.76E-01  1.70E+02
 
 TH 5
+       -1.83E+00 -8.57E+00 -1.11E+00 -1.67E+01  2.97E+01
 
 TH 6
+       -5.82E-01 -1.09E+00  2.97E-02  2.22E+00 -8.94E-01  3.72E+01
 
 TH 7
+        4.78E-01  5.85E+00 -2.62E-02 -7.47E+00  5.71E-01 -2.73E-01  2.47E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.10E-02 -3.92E+00 -1.86E-01 -4.65E+01  6.67E+00 -1.27E+00  1.78E+00  0.00E+00  3.39E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.07E+00 -3.33E+00 -4.15E-04 -1.20E+01  1.02E+00  1.63E+00  4.53E-01  0.00E+00  3.45E+00  0.00E+00  9.55E+00
 
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
+        1.93E+02
 
 TH 2
+        7.39E+01  5.77E+01
 
 TH 3
+        3.43E-01  1.21E-01  2.68E-02
 
 TH 4
+        9.79E+01  4.61E+01  2.65E-01  1.74E+02
 
 TH 5
+       -1.70E+01 -2.15E+00 -6.36E-01 -1.84E+01  1.81E+01
 
 TH 6
+        2.99E+01  1.28E+01 -4.80E-02 -1.61E+01 -5.87E-01  4.33E+01
 
 TH 7
+        4.37E-01  6.20E+00 -2.41E-02 -7.57E+00  2.63E+00  2.40E+00  2.85E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.15E+01 -2.96E+00 -1.92E-02 -3.96E+01  3.16E+00  6.75E+00  1.62E+00  0.00E+00  2.14E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.41E+01 -2.16E+01 -2.36E-01 -2.00E+01  7.62E+00  7.81E+00 -4.03E-01  0.00E+00  2.09E+00  0.00E+00  1.92E+02
 
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
 #CPUT: Total CPU Time in Seconds,       72.381
Stop Time:
Thu Sep 30 09:14:44 CDT 2021
