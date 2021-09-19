Sat Sep 18 13:23:51 CDT 2021
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
$DATA ../../../../data/spa/S2/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1697.81647087872        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.6187E+01 -4.1483E+01 -5.6016E+01  6.6882E+00  9.2997E+01  1.6703E+01 -2.1828E+00  1.0730E+01 -2.9558E+00 -1.2572E+01
             3.5317E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1706.28562396412        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      176
 NPARAMETR:  1.0600E+00  1.0188E+00  1.0927E+00  1.0020E+00  9.7680E-01  9.2415E-01  9.9944E-01  9.3074E-01  1.0286E+00  1.0298E+00
             1.0042E+00
 PARAMETER:  1.5827E-01  1.1864E-01  1.8861E-01  1.0204E-01  7.6527E-02  2.1124E-02  9.9444E-02  2.8223E-02  1.2815E-01  1.2940E-01
             1.0423E-01
 GRADIENT:   8.3647E+00  1.7281E+00 -1.9406E+00 -3.1026E+00 -4.5092E+00 -1.2700E+01  1.9276E+00  5.6667E+00  2.2058E-01 -6.0002E+00
             4.8701E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1707.56805756571        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  1.0566E+00  8.9711E-01  1.2193E+00  1.0910E+00  9.8067E-01  9.3380E-01  8.5817E-01  7.8525E-01  1.0116E+00  1.1202E+00
             1.0036E+00
 PARAMETER:  1.5503E-01 -8.5782E-03  2.9827E-01  1.8714E-01  8.0477E-02  3.1503E-02 -5.2949E-02 -1.4175E-01  1.1153E-01  2.1355E-01
             1.0361E-01
 GRADIENT:   3.3336E+00  1.1711E+01  4.5939E+00  1.3096E+01 -9.8417E+00 -7.6223E+00 -1.7957E+00 -7.0014E-01 -1.2539E+00 -3.7447E+00
             6.7222E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1708.02337171369        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  1.0530E+00  7.4826E-01  1.3437E+00  1.1906E+00  9.7809E-01  9.5075E-01  1.0106E+00  8.6650E-01  9.4160E-01  1.1543E+00
             9.9329E-01
 PARAMETER:  1.5164E-01 -1.9000E-01  3.9541E-01  2.7447E-01  7.7849E-02  4.9494E-02  1.1058E-01 -4.3290E-02  3.9827E-02  2.4350E-01
             9.3263E-02
 GRADIENT:  -6.7527E-01  9.1461E+00  2.9358E+00  1.6040E+01 -4.8132E+00  1.9472E-01  1.5006E-01 -2.0402E-01 -8.3854E-02 -6.7441E-01
            -2.2635E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1708.44783971535        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  1.0499E+00  5.1623E-01  1.4133E+00  1.3264E+00  9.3044E-01  9.5185E-01  1.1268E+00  9.0154E-01  8.6217E-01  1.1529E+00
             9.9632E-01
 PARAMETER:  1.4873E-01 -5.6121E-01  4.4593E-01  3.8249E-01  2.7900E-02  5.0648E-02  2.1937E-01 -3.6534E-03 -4.8303E-02  2.4232E-01
             9.6308E-02
 GRADIENT:  -1.1973E+00  1.1491E+00  1.7969E-01  4.5971E-01 -2.3802E+00  1.6768E+00 -2.1852E-01  1.8176E-01 -9.9320E-01  1.2091E+00
            -8.1114E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1708.61292518041        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  1.0476E+00  3.3083E-01  1.4998E+00  1.4472E+00  9.0391E-01  9.4739E-01  1.3321E+00  9.7502E-01  8.0720E-01  1.1428E+00
             9.9650E-01
 PARAMETER:  1.4649E-01 -1.0061E+00  5.0534E-01  4.6962E-01 -1.0226E-03  4.5957E-02  3.8676E-01  7.4703E-02 -1.1418E-01  2.3352E-01
             9.6498E-02
 GRADIENT:  -2.6797E-01  1.8733E+00  1.1980E+00  1.0049E+01 -2.8260E+00  7.4818E-01 -1.8123E-01 -3.2648E-01 -6.4246E-01  1.9021E-01
            -2.8125E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1708.70362482505        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  1.0452E+00  1.9469E-01  1.5523E+00  1.5312E+00  8.8296E-01  9.4464E-01  2.0445E+00  1.0376E+00  7.6237E-01  1.1274E+00
             9.9675E-01
 PARAMETER:  1.4424E-01 -1.5363E+00  5.3971E-01  5.2604E-01 -2.4472E-02  4.3054E-02  8.1514E-01  1.3692E-01 -1.7133E-01  2.1995E-01
             9.6746E-02
 GRADIENT:  -5.5901E-01  9.7106E-01  1.0165E+00  8.5354E+00 -2.0438E+00  2.6772E-01 -6.9237E-02 -2.3518E-01 -1.6478E+00 -1.6604E-01
            -2.8254E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1708.74859785997        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  1.0444E+00  1.3933E-01  1.5501E+00  1.5599E+00  8.6887E-01  9.4259E-01  2.8868E+00  1.0426E+00  7.4899E-01  1.1173E+00
             9.9699E-01
 PARAMETER:  1.4345E-01 -1.8709E+00  5.3830E-01  5.4464E-01 -4.0563E-02  4.0874E-02  1.1602E+00  1.4169E-01 -1.8903E-01  2.1092E-01
             9.6984E-02
 GRADIENT:  -1.5453E-01  1.6069E-01  5.8128E-01 -9.5527E-01 -6.1084E-01 -3.2082E-01  1.3487E-01 -1.8863E-01 -7.3404E-02 -1.1927E-02
             4.2591E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1708.81415252460        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1411
 NPARAMETR:  1.0429E+00  5.8713E-02  1.6310E+00  1.6132E+00  8.7222E-01  9.4153E-01  4.1951E+00  1.1172E+00  7.3158E-01  1.1260E+00
             9.9736E-01
 PARAMETER:  1.4205E-01 -2.7351E+00  5.8922E-01  5.7820E-01 -3.6717E-02  3.9751E-02  1.5339E+00  2.1083E-01 -2.1255E-01  2.1865E-01
             9.7353E-02
 GRADIENT:  -3.2749E-01  9.0570E-02  5.0190E-01  3.7398E-01 -4.1907E-01 -2.9961E-01  9.6756E-03 -1.6466E-01 -7.6288E-01 -7.4685E-02
             6.2885E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1708.85018360604        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1588
 NPARAMETR:  1.0423E+00  1.6781E-02  1.6569E+00  1.6403E+00  8.6916E-01  9.4139E-01  5.1104E+00  1.1450E+00  7.2625E-01  1.1280E+00
             9.9722E-01
 PARAMETER:  1.4142E-01 -3.9875E+00  6.0493E-01  5.9488E-01 -4.0226E-02  3.9603E-02  1.7313E+00  2.3536E-01 -2.1987E-01  2.2047E-01
             9.7213E-02
 GRADIENT:  -1.6847E-01  3.2481E-02  1.0015E-01  2.3837E+00 -7.1809E-01 -1.3449E-01  4.2231E-04  4.5540E-03  2.8490E-01  2.0006E-01
             6.9885E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1708.85742641040        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1765
 NPARAMETR:  1.0422E+00  1.0000E-02  1.6671E+00  1.6440E+00  8.7071E-01  9.4163E-01  5.2706E+00  1.1541E+00  7.2404E-01  1.1281E+00
             9.9710E-01
 PARAMETER:  1.4137E-01 -4.5278E+00  6.1109E-01  5.9712E-01 -3.8451E-02  3.9853E-02  1.7621E+00  2.4334E-01 -2.2291E-01  2.2051E-01
             9.7093E-02
 GRADIENT:   6.7013E-02  0.0000E+00 -5.9557E-05 -8.3857E-02  2.7257E-02  1.4623E-02  4.2499E-04  1.7500E-03 -2.6353E-03 -5.3534E-03
             5.2370E-04

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1708.85744068467        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1949
 NPARAMETR:  1.0422E+00  1.0000E-02  1.6667E+00  1.6439E+00  8.7057E-01  9.4162E-01  5.2539E+00  1.1538E+00  7.2406E-01  1.1280E+00
             9.9708E-01
 PARAMETER:  1.4138E-01 -4.5279E+00  6.1083E-01  5.9705E-01 -3.8605E-02  3.9848E-02  1.7590E+00  2.4303E-01 -2.2289E-01  2.2047E-01
             9.7080E-02
 GRADIENT:   1.3758E-02  0.0000E+00  1.2989E-02 -2.5262E-01  5.5342E-03  4.8138E-03  4.3461E-04  2.3196E-03  1.7121E-03  2.5266E-03
            -2.5087E-03

0ITERATION NO.:   56    OBJECTIVE VALUE:  -1708.85744068467        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:     2000
 NPARAMETR:  1.0422E+00  1.0000E-02  1.6666E+00  1.6439E+00  8.7056E-01  9.4162E-01  5.2447E+00  1.1535E+00  7.2405E-01  1.1280E+00
             9.9708E-01
 PARAMETER:  1.4138E-01 -4.5279E+00  6.1083E-01  5.9705E-01 -3.8605E-02  3.9848E-02  1.7590E+00  2.4303E-01 -2.2289E-01  2.2047E-01
             9.7080E-02
 GRADIENT:  -1.9499E-03  0.0000E+00  1.2273E-02  5.2160E-03  1.6681E-02 -2.1743E-03  9.6763E-04  4.5919E-03  7.7349E-04  1.8861E-03
             5.4414E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2000
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.5736E-04 -1.4570E-03 -2.6047E-02 -6.7706E-03 -3.3883E-02
 SE:             2.9856E-02  8.9105E-04  1.5067E-02  2.9205E-02  2.2540E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9579E-01  1.0202E-01  8.3860E-02  8.1667E-01  1.3278E-01

 ETASHRINKSD(%)  1.0000E-10  9.7015E+01  4.9523E+01  2.1597E+00  2.4487E+01
 ETASHRINKVR(%)  1.0000E-10  9.9911E+01  7.4521E+01  4.2728E+00  4.2978E+01
 EBVSHRINKSD(%)  4.4296E-01  9.7145E+01  5.3143E+01  2.4395E+00  2.0423E+01
 EBVSHRINKVR(%)  8.8395E-01  9.9918E+01  7.8044E+01  4.8195E+00  3.6676E+01
 RELATIVEINF(%)  9.5836E+01  3.0120E-03  4.6206E+00  4.1206E+00  8.5402E+00
 EPSSHRINKSD(%)  4.4348E+01
 EPSSHRINKVR(%)  6.9029E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1708.8574406846685     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -973.70661412093034     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.47
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1708.857       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  1.67E+00  1.64E+00  8.71E-01  9.42E-01  5.25E+00  1.15E+00  7.24E-01  1.13E+00  9.97E-01
 


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
+        1.13E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        4.45E-01  0.00E+00  7.10E+01
 
 TH 4
+       -1.44E+01  0.00E+00 -2.90E+01  7.47E+02
 
 TH 5
+        3.87E+01  0.00E+00 -1.82E+02 -5.07E+01  7.32E+02
 
 TH 6
+       -2.73E+01  0.00E+00 -1.87E+00  2.74E-01 -2.82E+01  2.12E+02
 
 TH 7
+        7.66E-02  0.00E+00  8.97E-03 -3.56E-02 -1.19E-01 -1.23E-01 -9.57E-04
 
 TH 8
+        6.79E+00  0.00E+00 -1.95E+01 -3.35E+00 -4.60E+00 -1.31E+00  9.07E-02  1.96E+01
 
 TH 9
+        1.08E+01  0.00E+00  7.17E+00 -2.23E+00 -5.87E+00  1.44E+01  3.05E-02  1.00E+00  3.43E+02
 
 TH10
+        1.03E+01  0.00E+00 -6.04E+00 -5.72E-01 -6.36E+01 -4.48E+00  9.90E-02  1.39E+01 -4.30E+00  6.41E+01
 
 TH11
+       -1.61E+01  0.00E+00 -9.38E+00 -9.81E+00 -1.16E+01  2.40E+01 -2.43E-01  7.10E+00 -6.78E+00  1.46E+01  1.84E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.307
Stop Time:
Sat Sep 18 13:24:38 CDT 2021
