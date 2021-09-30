Thu Sep 30 09:54:42 CDT 2021
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
$DATA ../../../../data/spa2/D/dat85.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   29689.9986508782        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.7189E+02  6.2660E+02  3.3793E+00  5.1347E+02  1.3487E+02 -2.7580E+03 -1.3969E+03 -5.6583E+01 -2.1628E+03 -6.8829E+02
            -5.6445E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -445.615904241197        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2050E+00  1.2890E+00  1.0240E+00  1.4050E+00  9.6237E-01  1.5993E+00  1.2393E+00  9.8342E-01  1.1793E+00  1.0212E+00
             1.4862E+01
 PARAMETER:  2.8648E-01  3.5385E-01  1.2370E-01  4.4005E-01  6.1647E-02  5.6954E-01  3.1455E-01  8.3282E-02  2.6494E-01  1.2099E-01
             2.7988E+00
 GRADIENT:  -3.9867E+01  3.6450E+01  4.5423E+00  5.9595E+01 -2.2718E+01 -1.0027E+00 -1.9990E+01  2.7118E+00 -5.7916E+00  9.3918E+00
             6.0969E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -517.178837369861        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.1402E+00  7.9051E-01  4.4543E+00  1.8279E+00  1.3292E+00  2.1815E+00  5.5122E+00  3.4712E-01  1.1015E+00  1.5089E-01
             1.4469E+01
 PARAMETER:  2.3123E-01 -1.3507E-01  1.5939E+00  7.0314E-01  3.8461E-01  8.8001E-01  1.8070E+00 -9.5808E-01  1.9666E-01 -1.7912E+00
             2.7720E+00
 GRADIENT:  -3.0370E+01  6.4559E+00  1.1814E+01  7.1757E+01 -6.4181E+01  5.6577E+01  2.0207E-01  8.6308E-02 -7.1566E+00  1.2255E-01
             1.6616E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -546.399512942502        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.1527E+00  1.0728E+00  1.1212E+01  1.5247E+00  2.5821E+00  1.7341E+00  4.9059E+00  4.8941E-01  1.0888E+00  9.0170E-01
             1.3215E+01
 PARAMETER:  2.4211E-01  1.7026E-01  2.5170E+00  5.2177E-01  1.0486E+00  6.5046E-01  1.6904E+00 -6.1456E-01  1.8503E-01 -3.4716E-03
             2.6814E+00
 GRADIENT:   1.7780E+00 -3.0377E+00 -8.7794E-01  2.0450E+00  4.6599E+00  8.5624E+00  7.7829E+00  1.1370E-02  9.2220E+00  3.7901E+00
             1.2259E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -554.686427176923        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.0993E+00  9.6258E-01  1.3765E+01  1.4424E+00  2.6838E+00  1.7665E+00  5.5752E+00  1.1756E+00  8.6991E-01  7.1631E-01
             1.1972E+01
 PARAMETER:  1.9464E-01  6.1862E-02  2.7222E+00  4.6630E-01  1.0872E+00  6.6898E-01  1.8183E+00  2.6175E-01 -3.9361E-02 -2.3364E-01
             2.5826E+00
 GRADIENT:  -7.8250E+00 -1.3621E+00 -8.4421E-01 -1.5058E+01  7.9272E+00 -4.9456E+00 -2.9060E+00  3.8473E-02  3.5863E+00  2.3529E+00
             1.0167E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -557.622713594639        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.1138E+00  6.6653E-01  1.6921E+01  1.6116E+00  2.3886E+00  1.7790E+00  6.4793E+00  1.1535E+00  9.3735E-01  1.9481E-01
             1.1967E+01
 PARAMETER:  2.0780E-01 -3.0567E-01  2.9286E+00  5.7721E-01  9.7072E-01  6.7605E-01  1.9686E+00  2.4282E-01  3.5303E-02 -1.5357E+00
             2.5821E+00
 GRADIENT:  -1.9620E-01 -1.1663E+00 -4.5743E-01  7.0045E-01  6.9256E-01 -2.9854E-01  8.5718E-01  3.9931E-02 -1.1124E-01  2.0547E-01
            -1.6761E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -557.915694613572        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      688
 NPARAMETR:  1.1141E+00  7.1874E-01  4.5465E+01  1.5970E+00  2.4356E+00  1.7817E+00  6.4079E+00  7.7924E-01  9.3728E-01  5.4755E-02
             1.1968E+01
 PARAMETER:  2.0808E-01 -2.3026E-01  3.9170E+00  5.6815E-01  9.9017E-01  6.7755E-01  1.9575E+00 -1.4944E-01  3.5222E-02 -2.8049E+00
             2.5822E+00
 GRADIENT:  -6.5476E-01  3.2039E-01 -4.7256E-02  4.6234E-01 -8.2496E-01 -3.7586E-01  2.4459E+00  2.1816E-03  2.4621E-01  1.6778E-02
            -5.9907E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -557.978163030246        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      854             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1161E+00  7.1063E-01  2.3752E+02  1.5999E+00  2.4870E+00  1.7885E+00  6.4812E+00  7.7713E-01  9.2110E-01  1.0000E-02
             1.1962E+01
 PARAMETER:  2.0986E-01 -2.4160E-01  5.5702E+00  5.6994E-01  1.0111E+00  6.8139E-01  1.9689E+00 -1.5215E-01  1.7810E-02 -4.5765E+00
             2.5818E+00
 GRADIENT:   9.3041E+00  1.3166E+00 -3.2360E-03  1.2060E+01  1.7998E-01  1.2630E+01  6.1398E+01  7.3832E-05 -2.6361E-01  0.0000E+00
             2.9972E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -557.983208824729        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1029
 NPARAMETR:  1.1173E+00  7.0493E-01  4.1527E+02  1.6071E+00  2.5016E+00  1.7881E+00  6.4690E+00  7.7696E-01  9.3423E-01  1.0000E-02
             1.1986E+01
 PARAMETER:  2.1091E-01 -2.4965E-01  6.1289E+00  5.7442E-01  1.0169E+00  6.8117E-01  1.9670E+00 -1.5237E-01  3.1963E-02 -4.6175E+00
             2.5837E+00
 GRADIENT:   3.1047E-01  3.9612E-02 -2.5414E-03 -1.9731E-01 -1.1886E-01  7.3328E-01  2.9504E+00  2.4121E-05 -2.3509E-02  0.0000E+00
             1.0835E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -557.994468404962        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1212             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1167E+00  6.9813E-01  2.3007E+03  1.6111E+00  2.5088E+00  1.7860E+00  6.5289E+00  7.7833E-01  9.3606E-01  1.0000E-02
             1.1971E+01
 PARAMETER:  2.1037E-01 -2.5935E-01  7.8410E+00  5.7689E-01  1.0198E+00  6.8000E-01  1.9762E+00 -1.5061E-01  3.3926E-02 -4.6263E+00
             2.5825E+00
 GRADIENT:   9.3496E+00  1.3348E+00 -4.2166E-04  1.1839E+01  6.8142E-01  1.2280E+01  6.2477E+01  7.0864E-07  1.1218E-01  0.0000E+00
             3.0779E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -557.997719429741        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1390
 NPARAMETR:  1.1162E+00  6.8789E-01  4.7270E+03  1.6143E+00  2.5101E+00  1.7863E+00  6.5439E+00  7.7898E-01  9.3748E-01  1.0000E-02
             1.1969E+01
 PARAMETER:  2.0995E-01 -2.7413E-01  8.5611E+00  5.7888E-01  1.0203E+00  6.8014E-01  1.9785E+00 -1.4977E-01  3.5438E-02 -4.6263E+00
             2.5823E+00
 GRADIENT:   1.1653E-01 -4.1368E-02 -2.1540E-04 -5.7481E-01  3.6914E-02  3.9517E-01  3.8449E+00 -1.7703E-06 -8.6521E-02  0.0000E+00
            -1.9512E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -557.999496714385        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1583             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1159E+00  6.8634E-01  5.7665E+05  1.6184E+00  2.5088E+00  1.7863E+00  6.5678E+00  7.7800E-01  9.4081E-01  1.0000E-02
             1.1966E+01
 PARAMETER:  2.0969E-01 -2.7638E-01  1.3365E+01  5.8146E-01  1.0198E+00  6.8017E-01  1.9822E+00 -1.5103E-01  3.8986E-02 -4.6263E+00
             2.5820E+00
 GRADIENT:  -9.5834E+01  7.9485E+00  8.3646E-07  5.5429E+01 -1.4566E+00  2.3906E+01  5.8880E+01 -6.5819E-07 -9.9390E+00  0.0000E+00
             5.0128E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -557.999881967397        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     1780
 NPARAMETR:  1.1160E+00  6.8432E-01  2.1029E+06  1.6182E+00  2.5096E+00  1.7863E+00  6.5689E+00  7.7799E-01  9.4094E-01  1.0000E-02
             1.1967E+01
 PARAMETER:  2.0971E-01 -2.7932E-01  1.4659E+01  5.8134E-01  1.0201E+00  6.8017E-01  1.9823E+00 -1.5104E-01  3.9127E-02 -4.6263E+00
             2.5822E+00
 GRADIENT:  -2.5870E+01  1.5909E+00 -7.1993E-07  4.1391E+00 -2.7704E-01 -1.2029E+00  4.4452E+00 -1.2595E-04 -2.1186E-01  0.0000E+00
             6.7064E+00

0ITERATION NO.:   61    OBJECTIVE VALUE:  -557.999881967397        NO. OF FUNC. EVALS.:  32
 CUMULATIVE NO. OF FUNC. EVALS.:     1812
 NPARAMETR:  1.1160E+00  6.8242E-01  2.1022E+06  1.6182E+00  2.5096E+00  1.7863E+00  6.5689E+00  7.7800E-01  9.4095E-01  1.0000E-02
             1.1968E+01
 PARAMETER:  2.0971E-01 -2.7932E-01  1.4659E+01  5.8134E-01  1.0201E+00  6.8017E-01  1.9823E+00 -1.5104E-01  3.9127E-02 -4.6263E+00
             2.5822E+00
 GRADIENT:   8.4574E-02  1.7068E-01  1.1806E-03 -8.0655E-02 -4.1784E-02  6.2859E-03 -1.4265E-01 -1.7787E-01 -6.2661E-02  0.0000E+00
            -2.3574E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1812
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1410E-02  4.4211E-02  2.9849E-09 -7.0931E-02  6.5777E-06
 SE:             2.7625E-02  2.2234E-02  1.0335E-09  1.3021E-02  3.7860E-05
 N:                     100         100         100         100         100

 P VAL.:         4.3834E-01  4.6762E-02  3.8761E-03  5.1210E-08  8.6207E-01

 ETASHRINKSD(%)  7.4519E+00  2.5513E+01  1.0000E+02  5.6378E+01  9.9873E+01
 ETASHRINKVR(%)  1.4348E+01  4.4517E+01  1.0000E+02  8.0971E+01  1.0000E+02
 EBVSHRINKSD(%)  9.2139E+00  2.1604E+01  1.0000E+02  5.7189E+01  9.9790E+01
 EBVSHRINKVR(%)  1.7579E+01  3.8540E+01  1.0000E+02  8.1672E+01  1.0000E+02
 RELATIVEINF(%)  8.0614E+01  3.0986E+01  0.0000E+00  8.0067E+00  6.1308E-05
 EPSSHRINKSD(%)  5.9365E+00
 EPSSHRINKVR(%)  1.1521E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -557.99988196739662     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       544.72635787821048     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    42.91
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -558.000       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.12E+00  6.84E-01  2.10E+06  1.62E+00  2.51E+00  1.79E+00  6.57E+00  7.78E-01  9.41E-01  1.00E-02  1.20E+01
 


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
+        2.18E+02
 
 TH 2
+        3.98E+00  4.11E+01
 
 TH 3
+       -7.11E-08  7.80E-08  3.92E-16
 
 TH 4
+       -9.34E+00  2.85E+01  5.58E-10  1.19E+02
 
 TH 5
+       -2.14E+00 -3.33E+00 -2.40E-09 -9.38E+00  7.79E+00
 
 TH 6
+       -4.36E+00  4.03E-01 -1.07E-08  3.45E+00 -7.22E-01  4.11E+01
 
 TH 7
+        7.17E-01  4.45E+00 -2.20E-09 -5.69E+00  1.99E-01 -5.04E-01  2.23E+00
 
 TH 8
+       -1.03E+01  1.13E+01  4.07E-07  5.20E+00 -2.37E+00  3.29E+00 -1.80E-01 -1.85E+01
 
 TH 9
+       -1.19E+01 -4.87E+00  1.85E-07 -3.68E+01  7.63E+00  2.27E+00  1.16E+00  8.53E+01 -4.79E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.46E+00 -2.49E+00 -1.40E-09 -9.17E+00  5.85E-01  1.91E+00  2.91E-01  1.82E-01  3.16E+00  0.00E+00  4.46E+00
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       55.841
Stop Time:
Thu Sep 30 09:55:39 CDT 2021
