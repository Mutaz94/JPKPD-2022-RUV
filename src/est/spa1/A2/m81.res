Wed Sep 29 23:42:43 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat81.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -396.957893411364        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0445E+02  7.5484E+01  8.7641E+01  7.7956E+01  2.4292E+02  2.3433E+01 -4.9633E+01 -9.3471E+01  1.2677E+00 -1.4811E+02
            -3.1180E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1607.42537616745        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0492E+00  9.2074E-01  9.6860E-01  1.0827E+00  8.9735E-01  1.1454E+00  8.5961E-01  8.9269E-01  7.9572E-01  7.9518E-01
             3.2109E+00
 PARAMETER:  1.4806E-01  1.7418E-02  6.8097E-02  1.7943E-01 -8.3096E-03  2.3576E-01 -5.1274E-02 -1.3514E-02 -1.2851E-01 -1.2919E-01
             1.2666E+00
 GRADIENT:   6.8587E+01 -2.3363E+01 -2.2251E+01 -1.4171E+01  4.0614E+01  4.2660E+01 -4.3057E+00  4.7634E+00 -1.7084E+01  9.5858E+00
             3.5327E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1616.21471765352        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  1.0455E+00  7.4469E-01  5.5452E-01  1.1626E+00  5.7254E-01  1.0759E+00  7.4058E-01  3.9399E-01  9.3413E-01  3.1432E-01
             3.2331E+00
 PARAMETER:  1.4450E-01 -1.9478E-01 -4.8965E-01  2.5065E-01 -4.5768E-01  1.7317E-01 -2.0032E-01 -8.3142E-01  3.1865E-02 -1.0574E+00
             1.2735E+00
 GRADIENT:   4.9898E+01  2.4450E+01 -7.1176E+00  7.2911E+01  1.5910E+01  1.5647E+01 -5.3582E+00 -1.3772E-01  1.4241E+00 -1.8831E+00
             4.0512E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1616.33247829540        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  1.0444E+00  7.3653E-01  4.6757E-01  1.1556E+00  5.1621E-01  1.0850E+00  8.1424E-01  3.8271E-01  9.1855E-01  2.9340E-01
             3.1980E+00
 PARAMETER:  1.4342E-01 -2.0581E-01 -6.6020E-01  2.4461E-01 -5.6123E-01  1.8155E-01 -1.0550E-01 -8.6049E-01  1.5041E-02 -1.1262E+00
             1.2625E+00
 GRADIENT:  -1.3918E+01  2.3330E+01 -2.0490E+01  7.3805E+01  1.8579E+01  7.6837E+00 -8.1947E+00 -6.4573E-01 -4.4302E+00 -2.6774E+00
             2.7686E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1627.02691952617        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.0449E+00  5.6447E-01  3.6510E-01  1.1954E+00  3.8964E-01  1.0122E+00  1.2312E+00  1.5447E-01  9.5869E-01  7.1522E-01
             2.7092E+00
 PARAMETER:  1.4392E-01 -4.7186E-01 -9.0759E-01  2.7846E-01 -8.4254E-01  1.1213E-01  3.0801E-01 -1.7678E+00  5.7811E-02 -2.3516E-01
             1.0967E+00
 GRADIENT:  -4.2463E+00  4.6816E+01 -2.3281E+01  1.1674E+02  3.6851E+00 -2.3876E+01  3.3447E+00  2.3607E-01 -1.4446E+01  9.5415E+00
            -1.9590E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1638.66202384857        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      621
 NPARAMETR:  1.0433E+00  3.6342E-01  3.3220E-01  1.1776E+00  3.2609E-01  1.0795E+00  9.8692E-01  3.1090E-02  1.0233E+00  6.5969E-01
             2.7013E+00
 PARAMETER:  1.4235E-01 -9.1218E-01 -1.0020E+00  2.6351E-01 -1.0206E+00  1.7650E-01  8.6834E-02 -3.3709E+00  1.2298E-01 -3.1599E-01
             1.0937E+00
 GRADIENT:   2.4379E+00  5.5066E+00  6.3799E-01  1.1197E+01 -1.2849E+00  3.2217E+00 -6.2494E-01  3.2396E-03  4.2702E+00 -2.4223E+00
            -1.1891E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1640.63239656789        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      796
 NPARAMETR:  1.0300E+00  1.4838E-01  4.2217E-01  1.3220E+00  3.4624E-01  1.0486E+00  1.4003E+00  1.0000E-02  9.3295E-01  7.0971E-01
             2.7471E+00
 PARAMETER:  1.2958E-01 -1.8080E+00 -7.6235E-01  3.7911E-01 -9.6063E-01  1.4749E-01  4.3665E-01 -8.0499E+00  3.0601E-02 -2.4290E-01
             1.1105E+00
 GRADIENT:  -4.1950E-01  3.6455E+00  2.0053E+01  1.6258E+01 -2.6973E+01 -2.6686E+00 -1.1914E-01  0.0000E+00  3.5911E+00  2.6272E-01
            -3.0900E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1641.81168169461        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      971
 NPARAMETR:  1.0214E+00  3.4419E-02  4.1953E-01  1.3636E+00  3.3516E-01  1.0506E+00  2.3517E+00  1.0000E-02  8.9746E-01  6.9872E-01
             2.7591E+00
 PARAMETER:  1.2119E-01 -3.2692E+00 -7.6863E-01  4.1009E-01 -9.9315E-01  1.4933E-01  9.5513E-01 -1.7194E+01 -8.1838E-03 -2.5850E-01
             1.1149E+00
 GRADIENT:  -4.0523E+00  5.1279E-01  5.4406E+00  1.6324E+01 -1.2853E+01 -6.0486E-01 -1.2501E-02  0.0000E+00 -5.7964E-01 -3.7047E-01
             8.5670E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1642.03414272042        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1154             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0218E+00  1.0000E-02  4.2021E-01  1.3637E+00  3.3462E-01  1.0510E+00  4.1156E+00  1.0000E-02  8.9450E-01  7.0161E-01
             2.7563E+00
 PARAMETER:  1.2152E-01 -4.7975E+00 -7.6700E-01  4.1021E-01 -9.9477E-01  1.4972E-01  1.5148E+00 -2.7272E+01 -1.1492E-02 -2.5437E-01
             1.1139E+00
 GRADIENT:   6.1683E+01  0.0000E+00  1.2411E+01  9.5495E+01  5.2438E+01  8.5386E+00  1.0806E-02  0.0000E+00  2.3235E+00  7.1072E-01
             1.1174E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1642.03478651999        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1347             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0217E+00  1.0000E-02  4.2012E-01  1.3637E+00  3.3460E-01  1.0509E+00  4.9191E+00  1.0000E-02  8.9363E-01  7.0187E-01
             2.7562E+00
 PARAMETER:  1.2148E-01 -4.7975E+00 -7.6720E-01  4.1022E-01 -9.9481E-01  1.4969E-01  1.6931E+00 -2.7272E+01 -1.2468E-02 -2.5401E-01
             1.1138E+00
 GRADIENT:   6.1568E+01  0.0000E+00  1.2170E+01  9.5623E+01  5.2644E+01  8.5236E+00  1.6809E-02  0.0000E+00  2.0556E+00  7.6093E-01
             1.1123E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1642.04879448748        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1504
 NPARAMETR:  1.0213E+00  1.0000E-02  4.1979E-01  1.3644E+00  3.3445E-01  1.0505E+00  1.2188E+01  1.0000E-02  8.9257E-01  7.0103E-01
             2.7561E+00
 PARAMETER:  1.2105E-01 -4.7975E+00 -7.6801E-01  4.1074E-01 -9.9526E-01  1.4926E-01  2.6004E+00 -2.7272E+01 -1.3652E-02 -2.5521E-01
             1.1138E+00
 GRADIENT:  -7.7656E-01  0.0000E+00  7.7609E+00  7.1493E+00 -1.2837E+01  1.0799E-02  9.5270E-01  0.0000E+00  3.7673E+00  3.1527E+00
             1.4789E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1642.18671903517        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1604
 NPARAMETR:  1.0217E+00  1.0000E-02  4.1709E-01  1.3602E+00  3.3483E-01  1.0507E+00  1.2395E+01  1.0000E-02  8.7958E-01  6.7817E-01
             2.7569E+00
 PARAMETER:  1.2144E-01 -4.7975E+00 -7.7446E-01  4.0766E-01 -9.9412E-01  1.4949E-01  2.6173E+00 -2.7272E+01 -2.8309E-02 -2.8835E-01
             1.1141E+00
 GRADIENT:   6.1503E+01  0.0000E+00  1.1834E+01  9.7246E+01  5.3056E+01  8.4821E+00  1.5387E+01  0.0000E+00  1.4681E+00  4.6634E-01
             9.9600E+00

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1642.18671903517        NO. OF FUNC. EVALS.:  67
 CUMULATIVE NO. OF FUNC. EVALS.:     1671
 NPARAMETR:  1.0217E+00  1.0000E-02  4.1761E-01  1.3611E+00  3.3433E-01  1.0507E+00  1.2342E+01  1.0000E-02  8.8046E-01  6.7922E-01
             2.7618E+00
 PARAMETER:  1.2144E-01 -4.7975E+00 -7.7446E-01  4.0766E-01 -9.9412E-01  1.4949E-01  2.6173E+00 -2.7272E+01 -2.8309E-02 -2.8835E-01
             1.1141E+00
 GRADIENT:   1.3252E-01  0.0000E+00 -1.1611E+02 -2.1933E+02  8.8937E+01  7.3286E-03  3.3776E+01  0.0000E+00 -4.4929E-01 -2.6382E-01
            -8.3086E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1671
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.9164E-04  5.5018E-03  8.5656E-05 -1.0256E-02 -4.8512E-03
 SE:             2.9243E-02  4.2444E-03  2.1672E-04  2.7239E-02  2.0816E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8931E-01  1.9489E-01  6.9267E-01  7.0653E-01  8.1572E-01

 ETASHRINKSD(%)  2.0314E+00  8.5781E+01  9.9274E+01  8.7454E+00  3.0265E+01
 ETASHRINKVR(%)  4.0215E+00  9.7978E+01  9.9995E+01  1.6726E+01  5.1371E+01
 EBVSHRINKSD(%)  2.0058E+00  8.8344E+01  9.9255E+01  7.8034E+00  2.9810E+01
 EBVSHRINKVR(%)  3.9715E+00  9.8641E+01  9.9994E+01  1.4998E+01  5.0734E+01
 RELATIVEINF(%)  9.5416E+01  1.0116E+00  3.6364E-04  2.8797E+01  3.3347E+00
 EPSSHRINKSD(%)  2.4932E+01
 EPSSHRINKVR(%)  4.3648E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1642.1867190351695     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -723.24818583049682     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.70
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1642.187       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.00E-02  4.17E-01  1.36E+00  3.35E-01  1.05E+00  1.24E+01  1.00E-02  8.80E-01  6.78E-01  2.76E+00
 


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
+        9.27E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.31E+00  0.00E+00  2.46E+04
 
 TH 4
+       -1.09E+01  0.00E+00  1.21E+04  7.94E+03
 
 TH 5
+        6.45E+01  0.00E+00 -4.65E+03 -1.21E+04  2.96E+04
 
 TH 6
+        3.57E-01  0.00E+00 -4.94E+00 -1.72E+01  1.62E+01  1.65E+02
 
 TH 7
+       -2.02E-01  0.00E+00  8.19E+00  2.60E+00 -1.75E+01  1.47E-01  2.17E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.37E+00  0.00E+00  7.90E+04 -1.30E+02  2.17E+02  4.90E+00  1.92E+00  0.00E+00  1.98E+02
 
 TH10
+       -6.27E+00  0.00E+00  3.54E+04 -9.37E+01  1.60E+02  1.06E+00  1.67E+00  0.00E+00  1.10E+01  1.26E+02
 
 TH11
+       -1.12E+01  0.00E+00  2.25E+03  1.32E+03  2.74E+01  1.64E-01  7.47E-02  0.00E+00  8.45E+03  3.81E+03  3.05E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.124
Stop Time:
Wed Sep 29 23:43:20 CDT 2021
