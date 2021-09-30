Wed Sep 29 09:22:57 CDT 2021
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
$DATA ../../../../data/int/D/dat62.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   29173.8076609797        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3767E+02  4.7131E+02 -5.6586E+01  2.3948E+02  3.8972E+02 -2.0397E+03 -1.0518E+03 -9.5309E+01 -2.0444E+03 -6.5667E+02
            -5.9720E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1009.16910242983        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.4666E+00  1.9031E+00  9.7521E-01  2.4306E+00  6.6064E-01  4.8796E+00  3.6227E+00  1.0010E+00  3.3497E+00  1.7302E+00
             1.2057E+01
 PARAMETER:  4.8294E-01  7.4349E-01  7.4901E-02  9.8815E-01 -3.1455E-01  1.6851E+00  1.3872E+00  1.0097E-01  1.3089E+00  6.4826E-01
             2.5896E+00
 GRADIENT:   1.4571E+01  6.7872E+01 -1.0364E+01  1.0195E+02 -8.4294E+01  1.2926E+02 -2.2436E+01  3.8979E+00  3.2743E+01  8.5281E+00
             4.1681E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1131.34594795465        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.4535E+00  3.3933E+00  2.9398E+01  3.5051E+00  2.6276E+00  3.4046E+00  1.0154E+01  6.7405E-01  3.0754E+00  2.2668E+00
             1.2074E+01
 PARAMETER:  4.7394E-01  1.3218E+00  3.4809E+00  1.3542E+00  1.0661E+00  1.3251E+00  2.4179E+00 -2.9446E-01  1.2234E+00  9.1838E-01
             2.5911E+00
 GRADIENT:   2.0464E+01  4.5586E+01 -7.3884E+00  1.4058E+02  3.0066E+01  1.0657E+02  5.7938E+01  3.0129E-02  7.5747E+00  6.0333E+01
             4.3262E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1299.25048853114        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.1336E+00  8.8276E-01  2.0332E+01  1.0371E+00  1.9019E+00  2.6353E+00  5.9523E+00  1.4953E+00  1.9718E+00  8.7374E-01
             9.7565E+00
 PARAMETER:  2.2536E-01 -2.4704E-02  3.1122E+00  1.3644E-01  7.4283E-01  1.0690E+00  1.8838E+00  5.0233E-01  7.7893E-01 -3.4973E-02
             2.3779E+00
 GRADIENT:  -2.6764E+01 -2.6082E+01  2.0064E+00 -5.1100E+01 -2.7086E+01  3.7364E+01  7.7198E+01  1.5918E-02  2.6981E+01  1.2651E+01
             3.3482E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1312.31921658416        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.1091E+00  1.5165E+00  1.2197E+01  8.0901E-01  2.1215E+00  2.3827E+00  4.4421E+00  3.2198E+00  2.0919E+00  7.9625E-01
             8.1588E+00
 PARAMETER:  2.0358E-01  5.1638E-01  2.6012E+00 -1.1194E-01  8.5210E-01  9.6823E-01  1.5911E+00  1.2693E+00  8.3809E-01 -1.2785E-01
             2.1991E+00
 GRADIENT:  -2.6150E+01 -1.4292E+01 -6.1168E-01 -1.5109E+01  2.4463E+01 -6.6401E+00 -1.7302E+00 -1.3644E-01  1.6629E+01  9.8307E+00
             3.9095E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1325.31412932334        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.1938E+00  1.5823E+00  5.7748E+00  8.6441E-01  1.8964E+00  2.4399E+00  4.6523E+00  1.9997E+00  1.3821E+00  3.5920E-01
             8.0821E+00
 PARAMETER:  2.7713E-01  5.5886E-01  1.8535E+00 -4.5712E-02  7.3996E-01  9.9194E-01  1.6374E+00  7.9297E-01  4.2358E-01 -9.2387E-01
             2.1897E+00
 GRADIENT:   4.6127E+00  3.1960E+00 -3.7171E-01  7.2340E-01  9.3688E-01  3.7769E+00  2.9553E+00 -3.4624E-01  4.2950E+00  1.6624E+00
            -8.4960E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1346.07031682379        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  1.2967E+00  1.2345E+00  1.0708E+01  1.1324E+00  1.9998E+00  2.7450E+00  7.0134E+00  2.0221E+00  1.1304E+00  2.6077E-01
             8.4231E+00
 PARAMETER:  3.5979E-01  3.1065E-01  2.4709E+00  2.2431E-01  7.9303E-01  1.1098E+00  2.0478E+00  8.0416E-01  2.2257E-01 -1.2441E+00
             2.2310E+00
 GRADIENT:   3.4123E+00  4.7975E+00 -1.5638E+00 -8.2435E-01  6.0884E+00  3.3883E+00  8.2701E-01 -3.1171E-01 -1.5457E+00  6.2979E-01
             3.2547E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1347.77419755860        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  1.2672E+00  8.1945E-01  2.2083E+01  1.2782E+00  2.0769E+00  2.7139E+00  7.6842E+00  3.6670E+00  1.3215E+00  2.4904E-01
             8.3694E+00
 PARAMETER:  3.3679E-01 -9.9124E-02  3.1948E+00  3.4544E-01  8.3087E-01  1.0984E+00  2.1392E+00  1.3994E+00  3.7879E-01 -1.2901E+00
             2.2246E+00
 GRADIENT:  -2.2861E+00 -1.3897E+00 -6.1991E-01  6.3845E-01 -7.4942E-01 -5.6469E-01 -1.8775E+00  1.9388E+00 -5.7911E-01  3.8978E-01
            -5.9321E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1347.85402905823        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  1.2797E+00  8.9641E-01  1.8934E+01  1.2548E+00  2.0577E+00  2.7180E+00  7.6087E+00  3.4360E+00  1.2903E+00  2.2176E-01
             8.4025E+00
 PARAMETER:  3.4663E-01 -9.3545E-03  3.0410E+00  3.2700E-01  8.2161E-01  1.0999E+00  2.1293E+00  1.3343E+00  3.5485E-01 -1.4062E+00
             2.2285E+00
 GRADIENT:   1.3553E-01  2.5250E-01 -2.3306E-01  2.9304E-01  2.3945E-01  1.8897E-01  3.5008E-01  2.3023E-01 -4.3930E-01  3.7251E-01
             9.8205E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1348.01467527770        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  1.2788E+00  8.7447E-01  1.7168E+01  1.2585E+00  2.0404E+00  2.7166E+00  7.6406E+00  3.3736E+00  1.3072E+00  5.5839E-02
             8.3985E+00
 PARAMETER:  3.4592E-01 -3.4139E-02  2.9431E+00  3.2995E-01  8.1313E-01  1.0994E+00  2.1335E+00  1.3160E+00  3.6792E-01 -2.7853E+00
             2.2281E+00
 GRADIENT:   4.6177E-02 -2.7595E-01 -3.4182E-01 -2.4111E-01  4.1548E-01  1.0464E-01  5.1668E-01  1.7155E-01  3.0688E-01  2.3245E-02
             4.1569E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1348.04311556457        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  1.2760E+00  8.7738E-01  1.9997E+01  1.2648E+00  2.0683E+00  2.7149E+00  7.6390E+00  3.4483E+00  1.3164E+00  2.9912E-02
             8.3975E+00
 PARAMETER:  3.4375E-01 -3.0819E-02  3.0956E+00  3.3492E-01  8.2675E-01  1.0988E+00  2.1333E+00  1.3379E+00  3.7490E-01 -3.4095E+00
             2.2279E+00
 GRADIENT:  -5.1018E-01 -2.7629E-02 -4.6694E-02  2.4749E-01  4.4595E-02 -1.0967E-01  1.4142E-01 -6.3178E-03  9.3165E-02  5.9450E-03
            -6.4257E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1348.09439741376        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1424             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2804E+00  8.5051E-01  2.0295E+01  1.2738E+00  2.0687E+00  2.7334E+00  7.8246E+00  3.4565E+00  1.3203E+00  1.0000E-02
             8.4031E+00
 PARAMETER:  3.4719E-01 -6.1920E-02  3.1104E+00  3.4201E-01  8.2691E-01  1.1055E+00  2.1573E+00  1.3403E+00  3.7785E-01 -4.6142E+00
             2.2286E+00
 GRADIENT:   2.9686E+01  1.4390E+00  1.1232E-01  1.2522E+01  4.4020E+00  4.9979E+01  1.5699E+02 -2.0836E-02  2.0050E+00  0.0000E+00
             4.1019E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1348.09893613724        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1581
 NPARAMETR:  1.2795E+00  8.4172E-01  2.0538E+01  1.2788E+00  2.0689E+00  2.7375E+00  7.8345E+00  3.4645E+00  1.3171E+00  1.0000E-02
             8.3994E+00
 PARAMETER:  3.4650E-01 -7.2312E-02  3.1223E+00  3.4592E-01  8.2704E-01  1.1070E+00  2.1585E+00  1.3426E+00  3.7543E-01 -4.6142E+00
             2.2282E+00
 GRADIENT:   2.9464E+01  1.3598E+00  1.4971E-01  1.4020E+01  3.9080E+00  5.0613E+01  1.5730E+02 -2.5828E-02  1.4772E+00  0.0000E+00
             3.9691E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1348.10027148017        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1768
 NPARAMETR:  1.2801E+00  8.3963E-01  2.0655E+01  1.2805E+00  2.0698E+00  2.7339E+00  7.8442E+00  3.4708E+00  1.3196E+00  1.0000E-02
             8.4004E+00
 PARAMETER:  3.4697E-01 -7.4792E-02  3.1280E+00  3.4727E-01  8.2744E-01  1.1057E+00  2.1598E+00  1.3444E+00  3.7735E-01 -4.6142E+00
             2.2283E+00
 GRADIENT:   3.7366E-01 -1.8263E-02 -4.3687E-02 -2.1568E-01 -1.5813E-01  2.2352E+00  3.3494E+00 -6.9823E-02 -7.0141E-02  0.0000E+00
            -3.1465E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1348.10127640038        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1957             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2802E+00  8.3877E-01  2.0841E+01  1.2822E+00  2.0711E+00  2.7338E+00  7.8585E+00  3.4795E+00  1.3227E+00  1.0000E-02
             8.4008E+00
 PARAMETER:  3.4699E-01 -7.5814E-02  3.1369E+00  3.4858E-01  8.2806E-01  1.1057E+00  2.1616E+00  1.3469E+00  3.7968E-01 -4.6142E+00
             2.2283E+00
 GRADIENT:   2.9631E+01  1.4286E+00  1.5158E-01  1.4192E+01  3.9444E+00  5.0055E+01  1.5810E+02 -4.4685E-02  1.6270E+00  0.0000E+00
             3.9896E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1348.10173106288        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2139
 NPARAMETR:  1.2802E+00  8.3774E-01  2.0946E+01  1.2832E+00  2.0718E+00  2.7338E+00  7.8624E+00  3.4850E+00  1.3236E+00  1.0000E-02
             8.4007E+00
 PARAMETER:  3.4699E-01 -7.7049E-02  3.1419E+00  3.4932E-01  8.2840E-01  1.1057E+00  2.1621E+00  1.3485E+00  3.8035E-01 -4.6142E+00
             2.2283E+00
 GRADIENT:   3.7862E-01  5.3784E-02 -2.9790E-02 -1.4854E-01 -1.6629E-01  2.2286E+00  3.6338E+00 -9.2311E-02 -1.6573E-02  0.0000E+00
            -2.7175E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1348.10213521966        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2328             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2802E+00  8.3365E-01  2.1065E+01  1.2838E+00  2.0730E+00  2.7338E+00  7.8637E+00  3.4929E+00  1.3253E+00  1.0000E-02
             8.4012E+00
 PARAMETER:  3.4699E-01 -8.1940E-02  3.1476E+00  3.4980E-01  8.2900E-01  1.1057E+00  2.1623E+00  1.3507E+00  3.8162E-01 -4.6142E+00
             2.2284E+00
 GRADIENT:   2.9624E+01  1.3065E+00  1.4055E-01  1.4173E+01  4.0971E+00  5.0041E+01  1.5857E+02 -3.4084E-02  1.6896E+00  0.0000E+00
             4.0059E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1348.10250105243        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2508            RESET HESSIAN, TYPE II
 NPARAMETR:  1.2802E+00  8.3433E-01  2.1162E+01  1.2844E+00  2.0737E+00  2.7337E+00  7.8669E+00  3.4979E+00  1.3261E+00  1.0000E-02
             8.4016E+00
 PARAMETER:  3.4701E-01 -8.1132E-02  3.1522E+00  3.5026E-01  8.2935E-01  1.1057E+00  2.1627E+00  1.3522E+00  3.8227E-01 -4.6142E+00
             2.2284E+00
 GRADIENT:   2.9624E+01  1.3551E+00  1.4728E-01  1.4262E+01  4.0872E+00  5.0027E+01  1.5856E+02 -4.5089E-02  1.6995E+00  0.0000E+00
             4.0120E+01

0ITERATION NO.:   87    OBJECTIVE VALUE:  -1348.10250105243        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:     2566
 NPARAMETR:  1.2802E+00  8.3433E-01  2.1162E+01  1.2844E+00  2.0737E+00  2.7337E+00  7.8669E+00  3.4979E+00  1.3261E+00  1.0000E-02
             8.4016E+00
 PARAMETER:  3.4701E-01 -8.1132E-02  3.1522E+00  3.5026E-01  8.2935E-01  1.1057E+00  2.1627E+00  1.3522E+00  3.8227E-01 -4.6142E+00
             2.2284E+00
 GRADIENT:   3.7659E-01 -2.2868E-02 -3.4996E-02 -2.7813E-01 -4.2258E-04  2.2229E+00  3.5889E+00 -9.1734E-02  5.2640E-02  0.0000E+00
             1.5468E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2566
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.9011E-03  4.0912E-02 -8.2340E-03 -7.9887E-02  7.8273E-06
 SE:             2.9361E-02  2.3710E-02  4.8852E-03  1.5000E-02  1.5267E-04
 N:                     100         100         100         100         100

 P VAL.:         7.6177E-01  8.4431E-02  9.1890E-02  1.0070E-07  9.5911E-01

 ETASHRINKSD(%)  1.6374E+00  2.0570E+01  8.3634E+01  4.9748E+01  9.9489E+01
 ETASHRINKVR(%)  3.2481E+00  3.6909E+01  9.7322E+01  7.4748E+01  9.9997E+01
 EBVSHRINKSD(%)  2.3597E+00  1.4236E+01  8.6928E+01  5.2937E+01  9.9457E+01
 EBVSHRINKVR(%)  4.6637E+00  2.6445E+01  9.8291E+01  7.7851E+01  9.9997E+01
 RELATIVEINF(%)  9.5205E+01  4.2725E+01  5.3419E-01  1.3040E+01  9.2150E-04
 EPSSHRINKSD(%)  7.3918E+00
 EPSSHRINKVR(%)  1.4237E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1348.1025010524299     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       305.98685871598082     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    95.37
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    19.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1348.103       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.28E+00  8.34E-01  2.12E+01  1.28E+00  2.07E+00  2.73E+00  7.87E+00  3.50E+00  1.33E+00  1.00E-02  8.40E+00
 


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
+        5.14E+01
 
 TH 2
+        2.84E+00  3.43E+00
 
 TH 3
+        4.15E-02  1.08E-02  2.75E-03
 
 TH 4
+       -1.26E+01  1.86E+01 -8.15E-02  1.22E+02
 
 TH 5
+       -3.14E+00 -5.85E+00 -4.25E-01 -1.44E+01  7.17E+01
 
 TH 6
+        6.35E+00 -4.53E-01  2.70E-03 -6.32E+00  1.07E+00  9.90E-01
 
 TH 7
+        1.37E+00 -8.19E-01  2.02E-03 -5.77E+00  9.84E-01  3.92E-01  2.86E-01
 
 TH 8
+       -9.65E-02  3.33E-01  1.09E-02  1.57E+00 -2.12E+00 -9.85E-02 -8.35E-02  7.39E-02
 
 TH 9
+        5.17E+00 -4.96E+00  1.86E-02 -3.35E+01  4.59E+00  1.95E+00  1.61E+00 -4.52E-01  9.28E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.83E-01 -8.89E-01  6.14E-03 -5.52E+00  2.61E-01  1.81E-01  2.41E-01 -6.75E-02  1.53E+00  0.00E+00  6.35E-01
 
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
+        8.43E+01
 
 TH 2
+       -1.76E-01  2.10E+01
 
 TH 3
+        1.14E-02  3.68E-02  1.13E-02
 
 TH 4
+       -1.81E+00  1.84E+01 -1.35E-01  1.23E+02
 
 TH 5
+       -1.48E+00 -4.53E+00 -4.29E-01 -9.45E+00  7.06E+01
 
 TH 6
+       -7.34E-01  5.41E-02  3.52E-03  7.30E-01  1.69E-01  2.17E+01
 
 TH 7
+        1.21E-01  2.29E+00 -6.04E-03 -5.91E+00  7.58E-01  1.28E-03  1.75E+00
 
 TH 8
+        9.03E-02 -4.50E-02 -1.05E-01  1.56E+00 -2.05E+00 -8.72E-02 -2.19E-02  2.10E+00
 
 TH 9
+        3.81E-01 -2.14E+00 -2.96E-02 -2.96E+01  4.60E+00 -4.14E-01  1.33E+00 -3.32E-01  1.86E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.42E+00 -1.56E+00  1.34E-02 -7.90E+00  5.78E-01  1.40E+00  1.03E-01 -4.11E-01  3.65E+00  0.00E+00  1.39E+01
 
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
+        8.68E+01
 
 TH 2
+        3.11E+01  1.98E+01
 
 TH 3
+        2.33E-02  2.72E-02  4.21E-03
 
 TH 4
+        4.50E+01  2.11E+01  1.20E-01  1.24E+02
 
 TH 5
+       -2.03E+00 -5.11E-01 -4.74E-01 -1.10E+01  6.87E+01
 
 TH 6
+        2.04E+01  7.34E+00 -3.35E-02 -1.70E+01  4.58E+00  2.88E+01
 
 TH 7
+        1.28E-01  4.30E+00 -1.47E-02 -6.44E+00  2.44E+00  2.64E+00  5.03E+00
 
 TH 8
+       -1.42E-02 -2.83E-02 -2.89E-03 -1.31E-01  2.19E-01  1.21E-02  1.17E-02  7.63E-03
 
 TH 9
+       -9.41E+00 -2.58E+00  1.54E-02 -4.75E+01  4.48E-02  9.79E+00  1.87E+00 -1.04E-02  3.91E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.64E+01  7.03E+00 -4.19E-02 -1.09E+02  2.17E+01  1.57E+01  2.03E+01 -9.89E-02  6.96E+01  0.00E+00  5.14E+02
 
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
 #CPUT: Total CPU Time in Seconds,      114.551
Stop Time:
Wed Sep 29 09:24:55 CDT 2021
