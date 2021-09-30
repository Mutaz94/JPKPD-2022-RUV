Wed Sep 29 15:19:20 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat75.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m75.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1650.28490165972        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2105E+02 -4.7701E+01  9.6429E+00 -7.4655E+01  4.6364E-01  4.1895E+01 -2.4567E+00 -6.2609E-01 -9.9754E+00  6.7642E+00
            -1.4018E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1657.30902102020        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9686E-01  1.0422E+00  1.0054E+00  1.0693E+00  1.0033E+00  1.0226E+00  1.0099E+00  1.0047E+00  1.0156E+00  9.8224E-01
             1.0432E+00
 PARAMETER:  9.6850E-02  1.4130E-01  1.0540E-01  1.6703E-01  1.0331E-01  1.2232E-01  1.0982E-01  1.0470E-01  1.1548E-01  8.2082E-02
             1.4233E-01
 GRADIENT:   1.3507E-01  3.3490E+00  5.2518E+00 -3.1469E+00 -1.1086E+01  1.3856E+00 -1.0591E+00 -1.6925E+00 -2.0135E-01  3.0853E+00
             2.6203E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1657.42653317574        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  9.9736E-01  1.0397E+00  1.0183E+00  1.0717E+00  1.0068E+00  1.0194E+00  1.0612E+00  1.0616E+00  9.9671E-01  9.4827E-01
             1.0412E+00
 PARAMETER:  9.7354E-02  1.3893E-01  1.1812E-01  1.6929E-01  1.0678E-01  1.1922E-01  1.5942E-01  1.5981E-01  9.6706E-02  4.6888E-02
             1.4035E-01
 GRADIENT:   1.4472E+00  3.7934E+00  5.4969E+00 -3.3869E+00 -8.9530E+00  2.3641E-01 -1.9491E-01 -2.8702E-01 -6.8951E-01 -2.3516E-01
             1.3179E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1657.60273237331        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.9218E-01  7.9285E-01  1.3288E+00  1.2481E+00  1.0590E+00  1.0125E+00  9.8289E-01  1.3255E+00  9.7727E-01  1.0618E+00
             1.0336E+00
 PARAMETER:  9.2147E-02 -1.3212E-01  3.8431E-01  3.2160E-01  1.5731E-01  1.1241E-01  8.2743E-02  3.8176E-01  7.7006E-02  1.5993E-01
             1.3302E-01
 GRADIENT:  -2.3194E+00  5.3718E+00 -2.4134E+00  1.3192E+01  1.6368E+00 -1.3337E+00  3.3221E+00  2.1867E+00  5.1381E+00  2.3195E+00
            -7.9600E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1658.77295929893        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  9.9144E-01  5.6128E-01  1.3947E+00  1.3712E+00  1.0065E+00  1.0163E+00  5.3658E-01  1.2704E+00  9.1884E-01  1.0443E+00
             1.0392E+00
 PARAMETER:  9.1402E-02 -4.7754E-01  4.3271E-01  4.1565E-01  1.0650E-01  1.1615E-01 -5.2254E-01  3.3935E-01  1.5357E-02  1.4334E-01
             1.3842E-01
 GRADIENT:   2.0005E+00 -5.9446E+00 -6.8773E-01 -1.8840E+01  4.2328E+00  8.5167E-01  9.4002E-01 -1.3697E+00  1.4509E+00 -9.8392E-01
             1.6174E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1658.83334006426        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.9017E-01  5.1086E-01  1.4396E+00  1.4127E+00  1.0030E+00  1.0148E+00  4.4932E-01  1.3231E+00  8.9545E-01  1.0476E+00
             1.0371E+00
 PARAMETER:  9.0121E-02 -5.7166E-01  4.6435E-01  4.4553E-01  1.0296E-01  1.1472E-01 -7.0003E-01  3.8001E-01 -1.0430E-02  1.4646E-01
             1.3643E-01
 GRADIENT:   7.6111E-01 -1.0326E+00  6.4057E-01 -3.9434E+00 -1.2386E+00  4.3177E-01  5.5118E-01 -1.4790E-02 -4.8938E-01 -2.1866E-02
             6.2601E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1658.83434758469        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  9.8979E-01  4.9550E-01  1.4500E+00  1.4232E+00  1.0028E+00  1.0143E+00  4.0177E-01  1.3324E+00  8.9140E-01  1.0494E+00
             1.0365E+00
 PARAMETER:  8.9736E-02 -6.0219E-01  4.7154E-01  4.5288E-01  1.0282E-01  1.1421E-01 -8.1186E-01  3.8698E-01 -1.4965E-02  1.4825E-01
             1.3587E-01
 GRADIENT:   4.4419E-01 -1.1338E+00  2.4924E-01 -3.1423E+00 -6.1557E-01  2.7116E-01  4.3369E-01  1.2225E-01 -9.7309E-02  7.4057E-02
             4.0128E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1658.83675435958        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  9.8945E-01  4.8179E-01  1.4591E+00  1.4324E+00  1.0028E+00  1.0138E+00  3.4286E-01  1.3406E+00  8.8793E-01  1.0511E+00
             1.0360E+00
 PARAMETER:  8.9392E-02 -6.3024E-01  4.7783E-01  4.5932E-01  1.0281E-01  1.1375E-01 -9.7042E-01  3.9313E-01 -1.8862E-02  1.4987E-01
             1.3536E-01
 GRADIENT:   1.6619E-01 -1.2984E+00 -1.6902E-01 -2.6019E+00  1.0180E-01  1.2299E-01  3.1308E-01  2.2882E-01  2.6696E-01  1.3685E-01
             1.9931E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1658.84898489745        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1411
 NPARAMETR:  9.8904E-01  4.6701E-01  1.4701E+00  1.4423E+00  1.0036E+00  1.0132E+00  2.3219E-01  1.3511E+00  8.8465E-01  1.0537E+00
             1.0353E+00
 PARAMETER:  8.8980E-02 -6.6141E-01  4.8535E-01  4.6621E-01  1.0359E-01  1.1315E-01 -1.3602E+00  4.0095E-01 -2.2568E-02  1.5227E-01
             1.3469E-01
 GRADIENT:  -2.0786E-01 -1.5804E+00 -7.7498E-01 -2.0909E+00  1.1696E+00 -7.6462E-02  1.4590E-01  3.6532E-01  7.4925E-01  2.0039E-01
            -6.6319E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1659.57102467703        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1593
 NPARAMETR:  9.9141E-01  6.5846E-01  1.4219E+00  1.3308E+00  1.0431E+00  1.0146E+00  1.0000E-02  1.3367E+00  9.6313E-01  1.0770E+00
             1.0350E+00
 PARAMETER:  9.1374E-02 -3.1786E-01  4.5199E-01  3.8574E-01  1.4215E-01  1.1448E-01 -7.7702E+00  3.9023E-01  6.2438E-02  1.7415E-01
             1.3438E-01
 GRADIENT:  -7.9281E-01  7.7334E+00  3.6763E+00  2.0976E+01 -6.4947E+00 -1.7717E-01  0.0000E+00 -1.1952E+00 -2.2901E+00  4.7570E-01
            -3.8894E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1660.21776527827        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1769
 NPARAMETR:  9.9440E-01  9.0457E-01  1.3555E+00  1.1604E+00  1.1140E+00  1.0171E+00  1.0000E-02  1.4325E+00  1.1028E+00  1.1060E+00
             1.0358E+00
 PARAMETER:  9.4388E-02 -2.9609E-04  4.0420E-01  2.4878E-01  2.0792E-01  1.1697E-01 -1.1662E+01  4.5944E-01  1.9783E-01  2.0074E-01
             1.3517E-01
 GRADIENT:   5.2109E-01  5.2189E+00  1.9593E+00  6.5338E+00 -2.1546E+00  2.1732E-01  0.0000E+00 -1.0752E+00 -1.3705E+00 -1.6525E-01
             1.1048E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1660.25645446114        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1951
 NPARAMETR:  9.9533E-01  9.3324E-01  1.3439E+00  1.1360E+00  1.1236E+00  1.0172E+00  1.0000E-02  1.4633E+00  1.1267E+00  1.1095E+00
             1.0348E+00
 PARAMETER:  9.5316E-02  3.0903E-02  3.9555E-01  2.2754E-01  2.1654E-01  1.1709E-01 -1.1928E+01  4.8068E-01  2.1927E-01  2.0392E-01
             1.3418E-01
 GRADIENT:   2.1251E+00 -7.9386E-01 -8.0954E-02 -1.2117E+00  3.6944E-01  2.2205E-01  0.0000E+00  2.7517E-02  1.8885E-01 -2.6398E-02
             7.3917E-03

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1660.25654197462        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     2008
 NPARAMETR:  9.9508E-01  9.3417E-01  1.3440E+00  1.1366E+00  1.1235E+00  1.0172E+00  1.0000E-02  1.4631E+00  1.1265E+00  1.1096E+00
             1.0348E+00
 PARAMETER:  9.5066E-02  3.1903E-02  3.9568E-01  2.2801E-01  2.1640E-01  1.1709E-01 -1.1928E+01  4.8058E-01  2.1913E-01  2.0402E-01
             1.3418E-01
 GRADIENT:   1.5523E+00  5.4537E-01  2.1090E-01  5.0670E-01 -2.3062E-01  2.1243E-01  0.0000E+00 -5.6895E-02 -5.7332E-02  3.6705E-03
            -3.5510E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2008
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0135E-03 -6.4819E-04 -3.8565E-02 -3.3356E-03 -4.1071E-02
 SE:             2.9871E-02  1.8600E-04  1.7326E-02  2.9293E-02  2.1250E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7294E-01  4.9256E-04  2.6026E-02  9.0934E-01  5.3261E-02

 ETASHRINKSD(%)  1.0000E-10  9.9377E+01  4.1955E+01  1.8656E+00  2.8811E+01
 ETASHRINKVR(%)  1.0000E-10  9.9996E+01  6.6308E+01  3.6964E+00  4.9321E+01
 EBVSHRINKSD(%)  4.3613E-01  9.9433E+01  4.4926E+01  2.5251E+00  2.6214E+01
 EBVSHRINKVR(%)  8.7036E-01  9.9997E+01  6.9669E+01  4.9864E+00  4.5557E+01
 RELATIVEINF(%)  9.8812E+01  3.1440E-04  9.5225E+00  1.1376E+01  1.7226E+01
 EPSSHRINKSD(%)  4.5092E+01
 EPSSHRINKVR(%)  6.9851E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1660.2565419746181     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -925.10571541087995     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.61
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1660.257       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  9.34E-01  1.34E+00  1.14E+00  1.12E+00  1.02E+00  1.00E-02  1.46E+00  1.13E+00  1.11E+00  1.03E+00
 


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
+        1.08E+03
 
 TH 2
+       -1.72E+01  4.52E+02
 
 TH 3
+        2.82E+00  7.81E+01  9.63E+01
 
 TH 4
+       -9.65E+00  4.37E+02 -6.74E+00  6.44E+02
 
 TH 5
+        2.50E+00 -1.71E+02 -1.35E+02 -2.77E+01  4.02E+02
 
 TH 6
+       -1.36E-01 -2.41E+00  7.86E-01 -2.35E+00 -3.26E-01  1.90E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        7.22E-02 -2.34E+01 -2.45E+01 -2.05E+00 -1.11E+00  1.71E-01  0.00E+00  2.18E+01
 
 TH 9
+        2.81E+00 -9.40E+01  3.30E+00  2.42E+00  3.07E+00 -5.02E-02  0.00E+00 -2.36E-01  1.42E+02
 
 TH10
+        9.44E-01  4.42E-01 -4.19E+00  7.25E-03 -5.92E+01  6.95E-02  0.00E+00  1.05E+01  1.35E+00  6.07E+01
 
 TH11
+       -7.26E+00 -1.68E+01 -8.31E+00 -8.77E+00 -6.72E+00  1.63E+00  0.00E+00  4.50E+00  8.74E+00  9.04E+00  1.97E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.702
Stop Time:
Wed Sep 29 15:19:54 CDT 2021
