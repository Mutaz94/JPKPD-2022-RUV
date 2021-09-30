Thu Sep 30 09:45:12 CDT 2021
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
$DATA ../../../../data/spa2/D/dat76.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m76.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   24260.3871063173        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0003E+02  3.7865E+02 -2.7477E+01  2.0693E+02  2.4033E+02 -2.9434E+03 -1.1649E+03 -1.3480E+02 -1.8949E+03 -7.0528E+02
            -4.6075E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -568.505711438352        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.4005E+00  1.2649E+00  8.1899E-01  2.0256E+00  1.0714E+00  3.4949E+00  1.6903E+00  9.7982E-01  1.8013E+00  1.0477E+00
             1.2737E+01
 PARAMETER:  4.3684E-01  3.3499E-01 -9.9679E-02  8.0587E-01  1.6895E-01  1.3513E+00  6.2491E-01  7.9616E-02  6.8849E-01  1.4660E-01
             2.6445E+00
 GRADIENT:   1.5853E+01  2.7494E+01 -2.9737E+01  1.0922E+02  3.6610E+01  1.0773E+02 -1.5450E+01  4.2462E+00 -5.7311E+01  1.0667E+01
             6.1709E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -665.702380132141        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.2550E+00  1.5949E+00  2.5997E+00  1.4951E+00  1.8241E+00  3.5334E+00  4.6012E+00  3.7052E-01  2.9454E+00  8.7860E-01
             1.1702E+01
 PARAMETER:  3.2714E-01  5.6683E-01  1.0554E+00  5.0219E-01  7.0108E-01  1.3623E+00  1.6263E+00 -8.9285E-01  1.1802E+00 -2.9421E-02
             2.5597E+00
 GRADIENT:  -2.2359E+00 -1.8842E-01 -1.7316E+01  3.0633E+01  8.2441E+00  1.1137E+02  2.8313E+01  9.6151E-02  5.1160E+01  6.1601E+00
             1.8457E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -680.336713411185        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.1970E+00  1.1046E+00  2.7963E+00  1.4114E+00  1.3028E+00  3.0904E+00  5.0376E+00  4.6747E-01  1.9829E+00  1.0296E+00
             1.1605E+01
 PARAMETER:  2.7983E-01  1.9946E-01  1.1283E+00  4.4460E-01  3.6454E-01  1.2283E+00  1.7169E+00 -6.6042E-01  7.8456E-01  1.2921E-01
             2.5514E+00
 GRADIENT:  -1.5126E+01 -3.9030E+00  5.4344E+00 -3.5750E+00 -3.6257E+01  6.8487E+01  3.0346E+01  1.8039E-01  3.6528E+01  1.1846E+01
             1.8784E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -699.306788280492        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      417
 NPARAMETR:  1.2349E+00  9.7918E-01  4.8602E+00  1.3999E+00  1.8519E+00  2.9434E+00  5.9210E+00  5.6074E-01  1.7345E+00  9.0232E-01
             1.1270E+01
 PARAMETER:  3.1100E-01  7.8956E-02  1.6811E+00  4.3641E-01  7.1621E-01  1.1796E+00  1.8785E+00 -4.7849E-01  6.5074E-01 -2.7846E-03
             2.5222E+00
 GRADIENT:  -1.6027E+01 -7.1724E+00 -1.1364E+00 -3.4449E+01 -2.2965E+00  2.5084E+00 -7.1528E-03  8.2176E-02  2.9344E+01  6.8835E+00
             1.4521E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -713.092564658729        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  1.2369E+00  1.0488E+00  5.4482E+00  1.3980E+00  1.8470E+00  2.9473E+00  5.8191E+00  1.3355E-01  1.7353E+00  2.0956E-01
             9.8622E+00
 PARAMETER:  3.1264E-01  1.4762E-01  1.7953E+00  4.3501E-01  7.1354E-01  1.1809E+00  1.8611E+00 -1.9133E+00  6.5119E-01 -1.4627E+00
             2.3887E+00
 GRADIENT:  -3.8117E+00  4.4288E-01  1.3395E-01 -9.1025E+00 -6.5022E+00  9.1080E-01 -5.1095E+00  3.3445E-03  2.4338E+01  3.9019E-01
             8.5329E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -714.635463419991        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:      779             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2407E+00  9.4615E-01  5.7358E+00  1.4097E+00  1.8741E+00  2.9032E+00  6.2236E+00  1.0000E-02  1.6515E+00  1.0000E-02
             9.8129E+00
 PARAMETER:  3.1569E-01  4.4648E-02  1.8467E+00  4.4340E-01  7.2812E-01  1.1658E+00  1.9284E+00 -6.1311E+00  6.0167E-01 -6.9100E+00
             2.3837E+00
 GRADIENT:   1.6941E+01  5.7108E-01  2.9365E-02 -6.3184E+00 -2.3481E+00  5.5118E+01  9.3779E+01  0.0000E+00  2.6072E+01  0.0000E+00
             2.6077E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -715.103273193905        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      958
 NPARAMETR:  1.2433E+00  1.0054E+00  6.0529E+00  1.4171E+00  1.8776E+00  3.0041E+00  6.1576E+00  1.0000E-02  1.6260E+00  1.0000E-02
             9.9500E+00
 PARAMETER:  3.1779E-01  1.0539E-01  1.9005E+00  4.4863E-01  7.3001E-01  1.2000E+00  1.9177E+00 -6.1311E+00  5.8612E-01 -6.9100E+00
             2.3976E+00
 GRADIENT:  -3.0543E+00  1.6322E+00  5.5404E-01 -1.2732E+01 -6.1860E+00  7.7481E+00  1.5906E+00  0.0000E+00  2.2685E+01  0.0000E+00
             1.1309E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -715.714202077656        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1135
 NPARAMETR:  1.2486E+00  9.5829E-01  6.2099E+00  1.4276E+00  1.8821E+00  3.0182E+00  6.3743E+00  1.0000E-02  1.5849E+00  1.0000E-02
             9.9505E+00
 PARAMETER:  3.2200E-01  5.7398E-02  1.9261E+00  4.5602E-01  7.3240E-01  1.2046E+00  1.9523E+00 -6.1311E+00  5.6052E-01 -6.9100E+00
             2.3976E+00
 GRADIENT:  -1.9372E+00  1.5298E+00  5.8793E-01 -1.5621E+01 -5.7478E+00  9.3459E+00  4.8992E+00  0.0000E+00  2.2343E+01  0.0000E+00
             1.1976E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -717.830866908789        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1314
 NPARAMETR:  1.2823E+00  6.5089E-01  5.7746E+00  1.4887E+00  1.8977E+00  3.0743E+00  7.4223E+00  1.0000E-02  1.3688E+00  1.0000E-02
             9.8472E+00
 PARAMETER:  3.4862E-01 -3.2942E-01  1.8535E+00  4.9790E-01  7.4065E-01  1.2231E+00  2.1045E+00 -6.1311E+00  4.1394E-01 -6.9100E+00
             2.3872E+00
 GRADIENT:   5.4461E+00 -2.8175E+00 -7.1414E-01 -2.5113E+01  1.1694E+00  1.5275E+01  9.2108E+00  0.0000E+00  1.8139E+01  0.0000E+00
             1.1426E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -718.615673936791        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     1517
 NPARAMETR:  1.2731E+00  6.1135E-01  5.9245E+00  1.5034E+00  1.8944E+00  3.0546E+00  7.4300E+00  1.0000E-02  1.3149E+00  1.0000E-02
             9.8412E+00
 PARAMETER:  3.4148E-01 -3.9208E-01  1.8791E+00  5.0773E-01  7.3890E-01  1.2167E+00  2.1055E+00 -6.1311E+00  3.7375E-01 -6.9100E+00
             2.3866E+00
 GRADIENT:   3.9965E+00 -3.9223E+00 -3.2743E-01 -2.0980E+01 -1.3679E-01  1.3443E+01  4.6254E+00  0.0000E+00  1.5252E+01  0.0000E+00
            -2.7826E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -719.016688114309        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     1647
 NPARAMETR:  1.2754E+00  7.0764E-01  5.6209E+00  1.4992E+00  1.8945E+00  3.0606E+00  7.0263E+00  1.0000E-02  1.3175E+00  1.0000E-02
             9.7109E+00
 PARAMETER:  3.4330E-01 -2.4582E-01  1.8265E+00  5.0496E-01  7.3895E-01  1.2186E+00  2.0497E+00 -6.1311E+00  3.7575E-01 -6.9100E+00
             2.3732E+00
 GRADIENT:   5.0487E+00 -1.6746E+00 -4.8697E-01 -7.0054E+00 -1.2830E+00  1.3605E+01  6.5481E-01  0.0000E+00  1.3403E+01  0.0000E+00
            -2.4691E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -719.288148510609        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     1850             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2637E+00  7.4326E-01  5.9024E+00  1.4999E+00  1.9144E+00  3.0623E+00  7.1017E+00  1.0000E-02  1.3092E+00  1.0000E-02
             9.7743E+00
 PARAMETER:  3.3402E-01 -1.9671E-01  1.8754E+00  5.0539E-01  7.4942E-01  1.2192E+00  2.0603E+00 -6.1311E+00  3.6940E-01 -6.9100E+00
             2.3798E+00
 GRADIENT:   2.3038E+01  1.7031E+00 -1.2154E-01  5.7711E+00 -7.8180E-01  7.6023E+01  1.1996E+02  0.0000E+00  1.4401E+01  0.0000E+00
             1.2202E+01

0ITERATION NO.:   64    OBJECTIVE VALUE:  -719.294683676713        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:     1947
 NPARAMETR:  1.2535E+00  7.4098E-01  6.0160E+00  1.4986E+00  1.9170E+00  3.0676E+00  7.0644E+00  1.0000E-02  1.3100E+00  1.0000E-02
             9.7472E+00
 PARAMETER:  3.2559E-01 -2.0180E-01  1.8863E+00  5.0509E-01  7.5158E-01  1.2197E+00  2.0575E+00 -6.1311E+00  3.6961E-01 -6.9100E+00
             2.3787E+00
 GRADIENT:  -6.2090E+02 -1.4315E-01 -1.6452E-01  3.9577E+02  2.6720E+02 -1.4912E+02  5.3541E+01  0.0000E+00 -5.3668E+02  0.0000E+00
             6.1555E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1947
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.6624E-03  4.6153E-02  2.7511E-06 -7.8383E-02  1.4169E-05
 SE:             2.8125E-02  2.2082E-02  9.1891E-06  1.2988E-02  5.7624E-05
 N:                     100         100         100         100         100

 P VAL.:         7.5809E-01  3.6611E-02  7.6465E-01  1.5949E-09  8.0577E-01

 ETASHRINKSD(%)  5.7761E+00  2.6023E+01  9.9969E+01  5.6489E+01  9.9807E+01
 ETASHRINKVR(%)  1.1219E+01  4.5274E+01  1.0000E+02  8.1068E+01  1.0000E+02
 EBVSHRINKSD(%)  2.9367E+00  1.9293E+01  9.9951E+01  5.1041E+01  9.9725E+01
 EBVSHRINKVR(%)  5.7872E+00  3.4864E+01  1.0000E+02  7.6030E+01  9.9999E+01
 RELATIVEINF(%)  9.3860E+01  3.3833E+01  3.3688E-06  1.2218E+01  1.0667E-04
 EPSSHRINKSD(%)  8.2649E+00
 EPSSHRINKVR(%)  1.5847E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -719.29468367671291     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       383.43155616889419     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    46.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.49
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -719.295       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.25E+00  7.39E-01  5.97E+00  1.50E+00  1.92E+00  3.06E+00  7.08E+00  1.00E-02  1.31E+00  1.00E-02  9.76E+00
 


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
+        3.03E+04
 
 TH 2
+        2.07E+01  2.87E+01
 
 TH 3
+       -6.88E-01  1.76E-01  1.25E-01
 
 TH 4
+       -1.87E+00  1.01E+01  2.72E-01  8.95E+03
 
 TH 5
+        1.05E+02 -1.10E+01 -1.08E+00 -7.38E+00  2.44E+03
 
 TH 6
+       -3.28E+01  2.17E+00 -3.91E-02  3.70E-01  1.07E+01  3.45E+02
 
 TH 7
+        8.54E+00  2.86E+00 -1.22E-02 -4.83E+00 -2.51E+00  9.89E+01  2.53E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.49E-02  1.62E+01 -9.35E-01 -2.24E+01  3.96E+00 -1.64E-01  1.40E+00  0.00E+00  2.17E+04
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        2.09E+00 -2.28E+00  2.57E-02 -7.50E+00 -1.17E+00  6.68E+01 -6.23E-02  0.00E+00  2.16E+00  0.00E+00  1.46E+01
 
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
 #CPUT: Total CPU Time in Seconds,       58.423
Stop Time:
Thu Sep 30 09:46:21 CDT 2021
