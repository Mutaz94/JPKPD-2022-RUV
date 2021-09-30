Thu Sep 30 09:26:01 CDT 2021
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
$DATA ../../../../data/spa2/D/dat63.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m63.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   22994.7867578763        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1381E+02  2.7464E+02 -3.0737E+01  8.8902E+01  2.5895E+02 -1.9877E+03 -9.2874E+02 -4.3874E+01 -1.7398E+03 -6.1066E+02
            -4.5268E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -638.565794323769        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2811E+00  1.2197E+00  9.5324E-01  1.6047E+00  1.0716E+00  2.2529E+00  1.5217E+00  9.6607E-01  1.5469E+00  1.0706E+00
             1.3793E+01
 PARAMETER:  3.4771E-01  2.9860E-01  5.2116E-02  5.7295E-01  1.6911E-01  9.1223E-01  5.1982E-01  6.5481E-02  5.3626E-01  1.6823E-01
             2.7242E+00
 GRADIENT:  -4.6782E+01 -4.4206E+00 -2.2728E+01  2.3133E+01  2.6714E+01  3.3030E+01 -1.4923E+01  5.6739E+00 -3.1093E+01  1.4680E+01
             2.0273E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -710.261355204843        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.2194E+00  2.0448E+00  8.2988E+00  1.6242E+00  4.8025E+00  3.8825E+00  4.8735E+00  3.7640E-01  6.9824E+00  7.3139E-01
             1.1422E+01
 PARAMETER:  2.9836E-01  8.1528E-01  2.2161E+00  5.8500E-01  1.6691E+00  1.4565E+00  1.6838E+00 -8.7709E-01  2.0434E+00 -2.1280E-01
             2.5355E+00
 GRADIENT:  -7.7218E+00  7.0198E+00 -1.7449E+01  1.8905E+01  8.5912E+00  1.3957E+02  5.9234E+01  1.4895E-01  1.0267E+02  9.5305E-01
             1.9490E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -775.631713220098        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.2772E+00  7.3348E-01  1.8858E+01  1.9794E+00  2.3660E+00  2.4601E+00  3.3857E+00  1.2601E+00  3.1237E+00  1.0786E+00
             1.0736E+01
 PARAMETER:  3.4469E-01 -2.0996E-01  3.0369E+00  7.8278E-01  9.6119E-01  1.0002E+00  1.3196E+00  3.3116E-01  1.2390E+00  1.7566E-01
             2.4736E+00
 GRADIENT:  -2.6082E+01 -4.4741E+00 -2.3581E+00  4.3419E-01  7.2403E+00 -1.2805E+01  1.1675E+01  5.9988E-01  2.1250E+01  7.6953E+00
             1.7856E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -793.132922436325        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.2800E+00  5.7144E-01  2.8143E+01  1.7969E+00  2.0365E+00  2.5358E+00  2.4864E+00  8.4823E-01  2.8814E+00  3.2949E-01
             8.8394E+00
 PARAMETER:  3.4686E-01 -4.5959E-01  3.4373E+00  6.8607E-01  8.1123E-01  1.0305E+00  1.0108E+00 -6.4607E-02  1.1583E+00 -1.0102E+00
             2.2792E+00
 GRADIENT:   6.7930E+00 -4.8203E+00 -1.6003E+00 -7.8445E+00  4.0184E+00  5.7809E+00  6.7753E+00  2.0008E-01 -2.7517E+00  1.1698E+00
            -1.6694E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -810.276782026329        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.2452E+00  1.2461E+00  3.0561E+01  1.2924E+00  2.1321E+00  2.3249E+00  7.7130E-01  3.6817E-01  3.9631E+00  6.2854E-02
             8.9778E+00
 PARAMETER:  3.1929E-01  3.1998E-01  3.5197E+00  3.5649E-01  8.5711E-01  9.4367E-01 -1.5968E-01 -8.9922E-01  1.4770E+00 -2.6669E+00
             2.2948E+00
 GRADIENT:  -2.7799E-01 -5.1833E+00 -2.2421E+00 -1.9346E+00  6.1851E+00 -1.3965E+01  3.5637E+00  5.0329E-02  6.8490E+00  4.1263E-02
            -2.1146E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -813.874833202542        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      464
 NPARAMETR:  1.2544E+00  1.2004E+00  1.1352E+02  1.3417E+00  2.1162E+00  2.4329E+00  3.1521E-01  1.7353E-01  3.8249E+00  1.0000E-02
             9.2021E+00
 PARAMETER:  3.2664E-01  2.8267E-01  4.8320E+00  3.9394E-01  8.4960E-01  9.8909E-01 -1.0545E+00 -1.6514E+00  1.4415E+00 -4.5174E+00
             2.3194E+00
 GRADIENT:  -2.1926E+01 -2.4385E+00  2.6364E-01 -6.9345E+00 -1.9710E+00 -2.6822E+01  2.1642E-01  1.1726E-03 -3.2067E+01  4.0106E-04
            -2.2022E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -820.749162462500        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      641
 NPARAMETR:  1.3684E+00  1.5366E+00  6.6617E+01  1.1291E+00  2.1703E+00  2.6768E+00  3.9464E-01  4.3137E-01  5.1342E+00  1.2675E-02
             9.4917E+00
 PARAMETER:  4.1362E-01  5.2957E-01  4.2990E+00  2.2138E-01  8.7486E-01  1.0846E+00 -8.2977E-01 -7.4078E-01  1.7359E+00 -4.2682E+00
             2.3504E+00
 GRADIENT:   2.3505E+00  5.4508E+00  8.3574E-02  3.1063E+00 -7.9529E-01  1.4409E+00 -9.8209E-01  2.0632E-02 -1.2472E-01  1.4754E-03
            -3.6325E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -820.989163381661        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      816
 NPARAMETR:  1.3445E+00  1.5918E+00  5.1665E+01  1.0173E+00  2.1553E+00  2.6500E+00  5.1622E-01  4.0744E-01  5.3306E+00  2.0318E-02
             9.4805E+00
 PARAMETER:  3.9602E-01  5.6486E-01  4.0448E+00  1.1717E-01  8.6794E-01  1.0746E+00 -5.6123E-01 -7.9786E-01  1.7735E+00 -3.7962E+00
             2.3492E+00
 GRADIENT:  -1.8982E+00 -2.0268E+00 -2.7339E-01 -4.0145E-01  7.6968E-01 -1.0254E+00 -4.0266E-01  2.9561E-02  1.5287E+00  3.8866E-03
             1.2647E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -821.084105125239        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      996
 NPARAMETR:  1.3561E+00  1.6256E+00  5.9920E+01  9.9824E-01  2.1572E+00  2.6713E+00  5.8756E-01  1.4606E-01  5.3713E+00  1.0000E-02
             9.4722E+00
 PARAMETER:  4.0464E-01  5.8585E-01  4.1930E+00  9.8237E-02  8.6883E-01  1.0826E+00 -4.3177E-01 -1.8238E+00  1.7811E+00 -4.8938E+00
             2.3484E+00
 GRADIENT:   3.2558E-01 -1.9610E-01  2.9731E-03  8.4868E-02  3.8664E-02  1.1642E+00  6.1471E-03  2.9278E-03  4.7121E-01  0.0000E+00
            -1.0902E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -821.084254353616        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1171
 NPARAMETR:  1.3561E+00  1.6276E+00  5.9786E+01  9.9580E-01  2.1571E+00  2.6658E+00  5.8836E-01  1.1461E-01  5.3765E+00  1.0000E-02
             9.4756E+00
 PARAMETER:  4.0462E-01  5.8710E-01  4.1908E+00  9.5793E-02  8.6876E-01  1.0805E+00 -4.3041E-01 -2.0662E+00  1.7820E+00 -4.8928E+00
             2.3487E+00
 GRADIENT:   3.3661E-01 -1.5367E-01  1.7521E-03  3.9436E-02  2.0482E-02  6.0397E-01 -5.0033E-03  1.8104E-03  4.5997E-01  0.0000E+00
             1.9141E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -821.085682952707        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1350
 NPARAMETR:  1.3563E+00  1.6284E+00  5.9698E+01  9.9467E-01  2.1568E+00  2.6714E+00  5.8945E-01  8.3140E-02  5.3817E+00  1.0000E-02
             9.4728E+00
 PARAMETER:  4.0474E-01  5.8758E-01  4.1893E+00  9.4656E-02  8.6863E-01  1.0826E+00 -4.2857E-01 -2.3872E+00  1.7830E+00 -4.8948E+00
             2.3484E+00
 GRADIENT:   3.5619E-01 -1.6133E-01  1.2044E-03  4.3339E-02  1.7181E-02  1.1788E+00 -2.8252E-03  9.5624E-04  5.2747E-01  0.0000E+00
            -6.6726E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -821.086352709161        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1521
 NPARAMETR:  1.3561E+00  1.6298E+00  5.9563E+01  9.9311E-01  2.1561E+00  2.6717E+00  5.9063E-01  1.8853E-02  5.3872E+00  1.0000E-02
             9.4719E+00
 PARAMETER:  4.0460E-01  5.8844E-01  4.1870E+00  9.3081E-02  8.6830E-01  1.0827E+00 -4.2656E-01 -3.8711E+00  1.7840E+00 -4.8948E+00
             2.3483E+00
 GRADIENT:   3.2895E-01 -1.0856E-01  4.9717E-04  4.8747E-02 -3.3016E-02  1.2052E+00 -7.8894E-03  4.9453E-05  5.7027E-01  0.0000E+00
            -1.9590E-01

0ITERATION NO.:   63    OBJECTIVE VALUE:  -821.086386109958        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1617
 NPARAMETR:  1.3563E+00  1.6299E+00  5.9537E+01  9.9211E-01  2.1565E+00  2.6717E+00  5.9112E-01  1.8680E-02  5.3880E+00  1.0000E-02
             9.4718E+00
 PARAMETER:  4.0477E-01  5.8850E-01  4.1866E+00  9.2081E-02  8.6847E-01  1.0827E+00 -4.2574E-01 -3.8803E+00  1.7842E+00 -4.8948E+00
             2.3483E+00
 GRADIENT:   3.8092E-01 -2.4195E-01 -8.2509E-04 -1.0916E-02  2.1840E-02  1.2042E+00  4.1216E-03  4.8581E-05  5.8199E-01  0.0000E+00
            -1.4391E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1617
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2402E-02 -4.4341E-02  9.4590E-06  1.7069E-02  2.7498E-05
 SE:             2.8986E-02  8.7071E-03  3.8604E-06  2.5090E-02  5.3360E-05
 N:                     100         100         100         100         100

 P VAL.:         6.6877E-01  3.5401E-07  1.4275E-02  4.9630E-01  6.0632E-01

 ETASHRINKSD(%)  2.8919E+00  7.0830E+01  9.9987E+01  1.5945E+01  9.9821E+01
 ETASHRINKVR(%)  5.7002E+00  9.1491E+01  1.0000E+02  2.9348E+01  1.0000E+02
 EBVSHRINKSD(%)  2.8731E+00  7.7774E+01  9.9950E+01  9.5170E+00  9.9725E+01
 EBVSHRINKVR(%)  5.6636E+00  9.5060E+01  1.0000E+02  1.8128E+01  9.9999E+01
 RELATIVEINF(%)  9.3871E+01  2.5450E+00  2.2878E-05  4.3547E+01  6.2933E-04
 EPSSHRINKSD(%)  8.4479E+00
 EPSSHRINKVR(%)  1.6182E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -821.08638610995763     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       281.63985373564947     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.99
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.07
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -821.086       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.36E+00  1.63E+00  5.95E+01  9.92E-01  2.16E+00  2.67E+00  5.91E-01  1.87E-02  5.39E+00  1.00E-02  9.47E+00
 


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
+        7.51E+01
 
 TH 2
+       -2.66E+00  9.32E+01
 
 TH 3
+       -7.14E-04  5.23E-05  2.08E-04
 
 TH 4
+       -8.70E-01  2.64E+01  1.45E-03  3.20E+01
 
 TH 5
+       -7.94E-01 -1.28E+01 -1.33E-02 -4.99E+00  1.81E+01
 
 TH 6
+       -3.00E+00  1.85E+00 -1.30E-04 -2.68E-01  8.42E-03  1.96E+01
 
 TH 7
+        3.23E-01 -2.46E+01  2.52E-04 -2.37E+00  3.15E+00 -2.37E-01  1.60E+01
 
 TH 8
+        5.46E-02 -1.87E-02 -4.08E-05 -4.42E-01 -1.13E-02  1.18E-02 -6.56E-02  7.81E-03
 
 TH 9
+        6.61E-01 -1.08E+01  1.41E-03  2.30E+00  1.12E+00 -5.30E-01  2.19E+00 -2.77E-03  3.87E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.15E+00 -8.24E+00  5.94E-04 -3.13E+00  6.34E-01  8.04E-01  2.86E+00  1.57E-04  4.96E-01  0.00E+00  7.22E+00
 
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
 #CPUT: Total CPU Time in Seconds,       47.141
Stop Time:
Thu Sep 30 09:26:50 CDT 2021
