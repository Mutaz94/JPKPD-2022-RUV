Sat Sep 18 03:33:43 CDT 2021
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
$DATA ../../../../data/int/SL1/dat64.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3191.58495806543        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.0672E+01 -3.3067E+01 -9.3178E-01  7.9888E+01  1.6080E+02 -2.7587E+01 -6.9194E+01 -1.1385E+02 -6.8632E+01  4.1773E+00
            -1.1827E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3425.67553008276        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0167E+00  1.1759E+00  1.1412E+00  8.9964E-01  1.0248E+00  1.1645E+00  1.2248E+00  1.1446E+00  1.2322E+00  8.3759E-01
             1.4465E+00
 PARAMETER:  1.1660E-01  2.6203E-01  2.3204E-01 -5.7572E-03  1.2447E-01  2.5225E-01  3.0274E-01  2.3502E-01  3.0876E-01 -7.7224E-02
             4.6915E-01
 GRADIENT:   3.6416E+01  5.9576E+01  8.6515E-01 -1.7406E+00 -2.0241E+01  3.5500E+01  3.6602E+00 -1.2845E+01  2.0741E+00 -7.7109E+00
            -6.1916E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3430.88612750135        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0211E+00  1.0959E+00  1.2211E+00  9.5592E-01  9.9848E-01  1.1065E+00  1.4008E+00  1.8950E+00  1.2041E+00  7.0779E-01
             1.4040E+00
 PARAMETER:  1.2091E-01  1.9162E-01  2.9975E-01  5.4914E-02  9.8476E-02  2.0120E-01  4.3703E-01  7.3920E-01  2.8572E-01 -2.4561E-01
             4.3932E-01
 GRADIENT:   4.8344E+01  5.3305E+01 -7.8603E+00  8.8597E+00  6.5934E+00  1.5230E+01  1.5053E+01 -7.2007E+00  8.6236E+00 -9.8748E+00
            -5.3114E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3432.25684777355        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      278
 NPARAMETR:  1.0212E+00  1.0960E+00  1.2949E+00  9.6260E-01  1.0161E+00  1.1214E+00  1.4051E+00  2.0860E+00  1.1614E+00  7.2371E-01
             1.4013E+00
 PARAMETER:  1.2095E-01  1.9168E-01  3.5844E-01  6.1887E-02  1.1599E-01  2.1455E-01  4.4011E-01  8.3525E-01  2.4965E-01 -2.2336E-01
             4.3741E-01
 GRADIENT:   1.8535E+01  4.4148E+01 -8.7005E+00  2.1282E+00  7.2912E+00  1.3615E+01  1.5397E+01 -8.1014E+00  7.1525E+00 -1.0731E+01
            -5.1964E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3434.95917497914        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  1.0199E+00  1.0856E+00  1.4550E+00  9.8206E-01  1.0476E+00  1.0646E+00  1.3654E+00  2.5577E+00  1.0075E+00  8.8865E-01
             1.4001E+00
 PARAMETER:  1.1971E-01  1.8216E-01  4.7503E-01  8.1892E-02  1.4650E-01  1.6261E-01  4.1146E-01  1.0391E+00  1.0747E-01 -1.8049E-02
             4.3652E-01
 GRADIENT:   4.7526E+01  4.4550E+01 -1.3876E+00  3.2407E+00 -2.7449E+00 -9.4252E-01  2.1181E+01 -8.1559E-01  9.2982E-01 -7.3278E-01
            -3.8451E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3434.97966405254        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  1.0199E+00  1.0856E+00  1.4550E+00  9.8206E-01  1.0476E+00  1.0659E+00  1.3654E+00  2.5690E+00  1.0032E+00  8.9295E-01
             1.4001E+00
 PARAMETER:  1.1971E-01  1.8216E-01  4.7503E-01  8.1892E-02  1.4650E-01  1.6382E-01  4.1146E-01  1.0435E+00  1.0324E-01 -1.3223E-02
             4.3652E-01
 GRADIENT:   4.7493E+01  4.4487E+01 -1.6672E+00  2.6201E+00 -3.6750E+00 -4.1368E-01  2.1296E+01 -3.7475E-01  4.2960E-01 -3.4509E-01
            -3.8058E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3435.03580183473        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  1.0199E+00  1.0855E+00  1.4547E+00  9.8206E-01  1.0476E+00  1.0816E+00  1.3651E+00  2.5679E+00  1.0024E+00  8.9300E-01
             1.4004E+00
 PARAMETER:  1.1971E-01  1.8206E-01  4.7479E-01  8.1892E-02  1.4650E-01  1.7841E-01  4.1125E-01  1.0431E+00  1.0235E-01 -1.3167E-02
             4.3674E-01
 GRADIENT:   4.7018E+01  4.4349E+01 -1.7030E+00  2.4909E+00 -3.6590E+00  5.8470E+00  2.1274E+01 -4.4284E-01  2.6533E-01 -3.4474E-01
            -3.7604E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3438.12337015010        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      805
 NPARAMETR:  1.0042E+00  1.0080E+00  1.4927E+00  9.8606E-01  1.0559E+00  1.0810E+00  1.1247E+00  2.6253E+00  1.0011E+00  9.2522E-01
             1.4417E+00
 PARAMETER:  1.0422E-01  1.0799E-01  5.0060E-01  8.5965E-02  1.5440E-01  1.7788E-01  2.1750E-01  1.0652E+00  1.0107E-01  2.2272E-02
             4.6584E-01
 GRADIENT:   1.1699E+01 -4.0609E+01 -2.1528E+01 -3.2756E+01  5.8463E+00  6.4744E+00 -1.1919E+01 -4.1036E+00  7.0929E+00  1.5554E+01
             1.0778E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3442.88806267594        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  1.0006E+00  9.9845E-01  1.6532E+00  1.0117E+00  1.0583E+00  1.0713E+00  1.1817E+00  2.8397E+00  9.5660E-01  8.2255E-01
             1.4344E+00
 PARAMETER:  1.0063E-01  9.8445E-02  6.0273E-01  1.1161E-01  1.5665E-01  1.6885E-01  2.6692E-01  1.1437E+00  5.5631E-02 -9.5347E-02
             4.6074E-01
 GRADIENT:   3.8132E+00 -9.7235E+00 -1.0283E+01 -2.2191E+01  8.4397E+00  2.3925E+00 -6.4897E+00 -3.4274E+00  3.2348E+00  3.7756E+00
             9.4354E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3443.47378549751        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:      984
 NPARAMETR:  1.0002E+00  9.9253E-01  1.6850E+00  1.0187E+00  1.0564E+00  1.0700E+00  1.1959E+00  2.8858E+00  9.4850E-01  8.0753E-01
             1.4319E+00
 PARAMETER:  1.0018E-01  9.2498E-02  6.2177E-01  1.1852E-01  1.5485E-01  1.6764E-01  2.7891E-01  1.1598E+00  4.7128E-02 -1.1378E-01
             4.5903E-01
 GRADIENT:  -2.0645E+01 -9.2952E+00 -1.0213E+01 -2.5641E+01  4.8124E+00 -3.7162E+00 -6.6927E+00 -4.6143E+00  2.3428E+00  2.5798E+00
             7.3797E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3443.52524403931        NO. OF FUNC. EVALS.: 202
 CUMULATIVE NO. OF FUNC. EVALS.:     1186
 NPARAMETR:  1.0003E+00  9.9195E-01  1.6866E+00  1.0192E+00  1.0561E+00  1.0782E+00  1.1969E+00  2.8883E+00  9.4805E-01  8.0652E-01
             1.4318E+00
 PARAMETER:  1.0026E-01  9.1918E-02  6.2270E-01  1.1903E-01  1.5461E-01  1.7534E-01  2.7975E-01  1.1607E+00  4.6651E-02 -1.1503E-01
             4.5893E-01
 GRADIENT:  -2.0157E+01 -9.1472E+00 -1.0123E+01 -2.5494E+01  4.7941E+00 -6.7566E-01 -6.6425E+00 -4.5905E+00  2.3276E+00  2.5444E+00
             7.3716E+00

0ITERATION NO.:   54    OBJECTIVE VALUE:  -3443.52647340228        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     1330
 NPARAMETR:  1.0003E+00  9.9195E-01  1.6866E+00  1.0192E+00  1.0561E+00  1.0801E+00  1.1969E+00  2.8883E+00  9.4798E-01  8.0636E-01
             1.4318E+00
 PARAMETER:  1.0026E-01  9.1918E-02  6.2270E-01  1.1903E-01  1.5461E-01  1.7705E-01  2.7975E-01  1.1607E+00  4.6582E-02 -1.1522E-01
             4.5893E-01
 GRADIENT:  -1.2812E+05  4.8601E+05 -2.0631E+04  2.1577E+05 -3.1435E+05 -9.2599E-03  9.1805E+04  2.2087E+04 -2.5686E+05 -2.2292E+05
            -5.5967E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1330
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0488E-02 -1.4511E-02 -1.3729E-02  2.2032E-02 -4.7650E-02
 SE:             2.9815E-02  2.3317E-02  2.4795E-02  2.5124E-02  1.9584E-02
 N:                     100         100         100         100         100

 P VAL.:         7.2502E-01  5.3373E-01  5.7979E-01  3.8052E-01  1.4971E-02

 ETASHRINKSD(%)  1.1543E-01  2.1887E+01  1.6934E+01  1.5833E+01  3.4391E+01
 ETASHRINKVR(%)  2.3072E-01  3.8983E+01  3.1000E+01  2.9159E+01  5.6955E+01
 EBVSHRINKSD(%)  4.6020E-01  2.4523E+01  1.8754E+01  1.6219E+01  3.2098E+01
 EBVSHRINKVR(%)  9.1827E-01  4.3032E+01  3.3990E+01  2.9807E+01  5.3894E+01
 RELATIVEINF(%)  9.9076E+01  2.9509E+01  5.2225E+01  4.1737E+01  2.7419E+01
 EPSSHRINKSD(%)  2.1170E+01
 EPSSHRINKVR(%)  3.7859E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3443.5264734022799     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1789.4371136338691     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.27
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3443.526       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  9.92E-01  1.69E+00  1.02E+00  1.06E+00  1.08E+00  1.20E+00  2.89E+00  9.48E-01  8.06E-01  1.43E+00
 


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
+        1.78E+09
 
 TH 2
+       -1.15E+09  1.82E+09
 
 TH 3
+       -1.43E+02  2.88E+03  1.10E+07
 
 TH 4
+       -9.42E+08  9.52E+08  2.19E+03  1.21E+09
 
 TH 5
+        7.00E+08 -7.07E+08  1.13E+03 -5.78E+08  6.71E+08
 
 TH 6
+       -1.74E+03  3.47E+03 -1.68E+02  1.45E+03 -2.12E+03  1.68E+02
 
 TH 7
+       -3.41E+08  3.45E+08  1.95E+03  2.82E+08 -2.10E+08  5.28E+02  1.59E+08
 
 TH 8
+        1.71E+07 -1.73E+07  6.98E+02 -7.06E+02  1.05E+07  5.27E+01 -5.11E+06  5.69E+05
 
 TH 9
+        6.03E+08 -6.09E+08  6.46E+02 -5.51E+03  3.70E+08 -1.86E+03 -1.80E+08 -1.98E+02  7.15E+08
 
 TH10
+        6.15E+08 -1.24E+09  1.65E+02 -5.08E+08  7.55E+08 -1.90E+03 -3.68E+08 -5.07E+01  7.27E+03  7.44E+08
 
 TH11
+       -1.74E+08  1.76E+08  1.69E+03  1.44E+08 -1.07E+08 -2.68E+02  5.21E+07 -2.61E+06 -9.20E+07 -9.38E+07  4.14E+07
 
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
 #CPUT: Total CPU Time in Seconds,       54.680
Stop Time:
Sat Sep 18 03:34:40 CDT 2021
