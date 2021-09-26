Sat Sep 25 08:02:04 CDT 2021
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
$DATA ../../../../data/spa/A1/dat39.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m39.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1323.47639664066        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.3796E+01  3.2433E+01  1.8368E+01 -3.6527E+00  2.1043E+01  2.3859E+01 -1.8792E+01  1.4814E+01 -3.4115E+01 -4.1130E+01
            -6.1616E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1492.37218648212        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0132E+00  8.6364E-01  1.1292E+00  1.0753E+00  9.9621E-01  8.9891E-01  1.0077E+00  7.9533E-01  1.1249E+00  9.9966E-01
             1.9715E+00
 PARAMETER:  1.1307E-01 -4.6597E-02  2.2153E-01  1.7263E-01  9.6200E-02 -6.5764E-03  1.0766E-01 -1.2900E-01  2.1773E-01  9.9660E-02
             7.7877E-01
 GRADIENT:   6.8634E+01 -8.8844E+00  4.6655E-01 -2.5328E+01  9.5734E+00 -1.4944E+01  3.4858E+00  7.6226E+00  1.3318E+01 -5.7760E+00
             1.4078E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1498.41712507511        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0094E+00  7.0982E-01  8.2389E-01  1.1868E+00  7.9677E-01  9.3592E-01  1.0988E+00  1.6533E-01  8.8884E-01  9.1888E-01
             1.9176E+00
 PARAMETER:  1.0940E-01 -2.4274E-01 -9.3714E-02  2.7122E-01 -1.2719E-01  3.3779E-02  1.9420E-01 -1.6998E+00 -1.7833E-02  1.5401E-02
             7.5109E-01
 GRADIENT:   4.6842E+01  6.1179E+00 -4.2526E+01  3.0069E+01  5.4888E+01 -1.2306E+00 -5.7005E+00  6.7495E-01 -2.8153E+01  9.1128E+00
             7.2315E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1508.30719287844        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.8122E-01  4.7078E-01  4.7739E-01  1.2450E+00  4.6040E-01  9.4457E-01  1.4655E+00  4.9347E-02  9.4388E-01  5.3269E-01
             1.8308E+00
 PARAMETER:  8.1041E-02 -6.5337E-01 -6.3942E-01  3.1913E-01 -6.7567E-01  4.2976E-02  4.8220E-01 -2.9089E+00  4.2249E-02 -5.2981E-01
             7.0473E-01
 GRADIENT:  -2.9409E+01  3.0311E+01  3.2265E+01  3.8649E+01 -4.3618E+01 -8.6532E-02 -1.1229E+00  6.3195E-02  8.4938E+00 -5.3768E+00
            -6.2872E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1511.12840784952        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.9370E-01  2.9527E-01  3.7239E-01  1.2461E+00  3.5816E-01  9.4839E-01  1.9058E+00  1.0000E-02  8.9315E-01  5.3262E-01
             1.7974E+00
 PARAMETER:  9.3678E-02 -1.1199E+00 -8.8781E-01  3.2000E-01 -9.2676E-01  4.7006E-02  7.4491E-01 -4.7401E+00 -1.2997E-02 -5.2994E-01
             6.8637E-01
 GRADIENT:   2.0324E+00  6.5313E+00  9.4068E+00  4.3725E+00 -1.7012E+01  1.1289E+00  1.0557E+00  0.0000E+00 -2.7622E+00 -4.2772E-01
            -3.0673E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1516.88919552394        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.8105E-01  7.0545E-02  4.6913E-01  1.3840E+00  3.9342E-01  9.6100E-01  2.5179E+00  1.0000E-02  8.6952E-01  6.0916E-01
             1.8277E+00
 PARAMETER:  8.0866E-02 -2.5515E+00 -6.5689E-01  4.2500E-01 -8.3289E-01  6.0224E-02  1.0234E+00 -8.4088E+00 -3.9811E-02 -3.9568E-01
             7.0308E-01
 GRADIENT:   3.4597E+00  1.7687E-01 -1.2873E+01 -1.2011E+01  2.5625E+01  9.5511E+00 -5.4692E-01  0.0000E+00 -2.1145E-01 -3.3774E+00
             1.4141E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1517.43760713560        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  9.7958E-01  4.8711E-02  4.6630E-01  1.4103E+00  3.8484E-01  9.4846E-01  2.7920E+00  1.0000E-02  8.6499E-01  6.1713E-01
             1.8272E+00
 PARAMETER:  7.9365E-02 -2.9218E+00 -6.6292E-01  4.4377E-01 -8.5494E-01  4.7082E-02  1.1268E+00 -9.4838E+00 -4.5038E-02 -3.8267E-01
             7.0278E-01
 GRADIENT:   5.9895E-01  1.2753E+00  9.3971E+00  2.7028E+01 -1.5431E+01  4.6794E+00 -3.2307E-01  0.0000E+00 -6.5706E-01 -1.3475E+00
             2.2066E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1517.44047084387        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  9.7930E-01  4.4150E-02  4.6045E-01  1.4058E+00  3.8104E-01  9.4546E-01  2.8964E+00  1.0000E-02  8.6545E-01  6.1739E-01
             1.8232E+00
 PARAMETER:  7.9083E-02 -3.0202E+00 -6.7555E-01  4.4058E-01 -8.6484E-01  4.3911E-02  1.1635E+00 -9.8333E+00 -4.4511E-02 -3.8225E-01
             7.0059E-01
 GRADIENT:   4.9925E-01  1.0346E+00  7.5927E+00  2.1247E+01 -1.2416E+01  3.5248E+00 -2.8366E-01  0.0000E+00 -5.8913E-01 -1.0431E+00
             1.6609E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1517.44092756488        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      594
 NPARAMETR:  9.7911E-01  4.1137E-02  4.5753E-01  1.4038E+00  3.7907E-01  9.4412E-01  2.9734E+00  1.0000E-02  8.6563E-01  6.1749E-01
             1.8213E+00
 PARAMETER:  7.8893E-02 -3.0908E+00 -6.8192E-01  4.3921E-01 -8.7004E-01  4.2499E-02  1.1897E+00 -1.0076E+01 -4.4298E-02 -3.8209E-01
             6.9957E-01
 GRADIENT:   4.2889E-01  9.0693E-01  6.6594E+00  1.8456E+01 -1.0881E+01  3.0169E+00 -2.5817E-01  0.0000E+00 -5.3397E-01 -8.9983E-01
             1.4185E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1518.80763248739        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      751
 NPARAMETR:  9.7927E-01  1.1157E-02  5.4859E-01  1.4601E+00  4.3089E-01  9.3515E-01  4.6438E+00  1.0000E-02  8.4484E-01  6.6918E-01
             1.8216E+00
 PARAMETER:  7.9057E-02 -4.3957E+00 -5.0041E-01  4.7849E-01 -7.4190E-01  3.2947E-02  1.6355E+00 -1.3881E+01 -6.8613E-02 -3.0170E-01
             6.9974E-01
 GRADIENT:  -3.5916E+00  1.7414E-01  4.6164E+00 -6.8442E+00 -5.7609E+00 -2.9752E-01 -3.8877E-02  0.0000E+00 -1.0942E+00  2.2870E+00
             1.4463E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1518.87224369853        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      929
 NPARAMETR:  9.7994E-01  1.0000E-02  5.4513E-01  1.4632E+00  4.2883E-01  9.3480E-01  6.8237E+00  1.0000E-02  8.4614E-01  6.6000E-01
             1.8219E+00
 PARAMETER:  7.9734E-02 -5.3893E+00 -5.0673E-01  4.8066E-01 -7.4670E-01  3.2581E-02  2.0204E+00 -1.7060E+01 -6.7075E-02 -3.1552E-01
             6.9989E-01
 GRADIENT:  -1.9683E+00  0.0000E+00  3.2712E+00  1.0061E+00 -5.1821E+00 -4.8698E-01 -6.0157E-02  0.0000E+00 -6.0592E-01  6.5032E-01
             6.0278E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1518.91041827725        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1109
 NPARAMETR:  9.8063E-01  1.0000E-02  5.4538E-01  1.4631E+00  4.2969E-01  9.3597E-01  1.0670E+01  1.0000E-02  8.4717E-01  6.5726E-01
             1.8206E+00
 PARAMETER:  8.0442E-02 -6.4796E+00 -5.0627E-01  4.8055E-01 -7.4470E-01  3.3831E-02  2.4674E+00 -2.0569E+01 -6.5849E-02 -3.1967E-01
             6.9919E-01
 GRADIENT:  -2.6358E-02  0.0000E+00 -1.0218E-01  9.9077E-03  1.5191E-01  9.6492E-03  2.4190E-02  0.0000E+00  1.9424E-02 -3.0480E-02
            -1.2881E-02

0ITERATION NO.:   59    OBJECTIVE VALUE:  -1518.91048778611        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     1238
 NPARAMETR:  9.8065E-01  1.0000E-02  5.4539E-01  1.4630E+00  4.2966E-01  9.3595E-01  1.0614E+01  1.0000E-02  8.4711E-01  6.5736E-01
             1.8206E+00
 PARAMETER:  8.0456E-02 -6.4585E+00 -5.0626E-01  4.8052E-01 -7.4475E-01  3.3802E-02  2.4622E+00 -2.0501E+01 -6.5928E-02 -3.1953E-01
             6.9919E-01
 GRADIENT:   3.4611E-03  0.0000E+00  4.2285E-02 -4.7587E-02 -4.3288E-02 -5.6253E-04  1.3384E-03  0.0000E+00 -3.2806E-03 -1.7937E-03
             7.2642E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1238
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.6857E-05  8.3250E-04 -1.8209E-05 -6.6573E-03 -9.4486E-03
 SE:             2.9471E-02  1.9784E-03  2.4654E-04  2.8236E-02  2.2871E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9927E-01  6.7390E-01  9.4112E-01  8.1361E-01  6.7952E-01

 ETASHRINKSD(%)  1.2695E+00  9.3372E+01  9.9174E+01  5.4067E+00  2.3379E+01
 ETASHRINKVR(%)  2.5229E+00  9.9561E+01  9.9993E+01  1.0521E+01  4.1292E+01
 EBVSHRINKSD(%)  1.4332E+00  9.4226E+01  9.9113E+01  5.0828E+00  2.2774E+01
 EBVSHRINKVR(%)  2.8458E+00  9.9667E+01  9.9992E+01  9.9072E+00  4.0362E+01
 RELATIVEINF(%)  8.6456E+01  2.2370E-02  3.2645E-04  1.1369E+01  1.9561E+00
 EPSSHRINKSD(%)  3.6691E+01
 EPSSHRINKVR(%)  5.9920E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1518.9104877861073     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -783.75966122236912     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1518.910       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.00E-02  5.45E-01  1.46E+00  4.30E-01  9.36E-01  1.06E+01  1.00E-02  8.47E-01  6.57E-01  1.82E+00
 


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
+        1.28E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.05E+01  0.00E+00  2.95E+03
 
 TH 4
+       -2.70E+01  0.00E+00 -2.06E+02  6.53E+02
 
 TH 5
+        6.69E+01  0.00E+00 -4.74E+03 -1.43E+02  8.29E+03
 
 TH 6
+        5.14E+00  0.00E+00  6.52E+00 -7.15E+00 -8.73E-01  2.14E+02
 
 TH 7
+       -4.55E-04  0.00E+00 -9.29E-02 -6.31E-02  2.06E-01  2.34E-02  1.62E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -8.50E+00  0.00E+00  3.80E+01 -7.83E+00 -8.52E-01  6.51E+00 -1.05E-02  0.00E+00  2.30E+02
 
 TH10
+       -5.74E+00  0.00E+00 -3.90E+01  8.66E-01 -6.61E+01 -3.46E+00 -8.14E-02  0.00E+00  1.14E+00  1.68E+02
 
 TH11
+       -1.25E+01  0.00E+00 -2.17E+01 -6.77E+00  7.24E+00  1.32E+00 -1.30E-02  0.00E+00  6.65E+00  3.59E+01  7.89E+01
 
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
 #CPUT: Total CPU Time in Seconds,       18.504
Stop Time:
Sat Sep 25 08:02:24 CDT 2021
