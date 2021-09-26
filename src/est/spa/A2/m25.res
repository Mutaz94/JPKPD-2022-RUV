Sat Sep 25 08:33:03 CDT 2021
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
$DATA ../../../../data/spa/A2/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1283.59240857231        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.5669E+01 -3.6072E+01  1.5876E+01 -9.1215E+01  1.5997E+01  2.7195E+01 -7.4844E+00  4.6055E+00 -3.1663E+01 -1.0936E+01
            -7.2175E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1489.24271211806        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0381E+00  1.0097E+00  1.1113E+00  1.0649E+00  1.0538E+00  8.8949E-01  9.5897E-01  8.8813E-01  1.0468E+00  8.4852E-01
             1.8925E+00
 PARAMETER:  1.3743E-01  1.0965E-01  2.0557E-01  1.6284E-01  1.5237E-01 -1.7109E-02  5.8106E-02 -1.8636E-02  1.4577E-01 -6.4259E-02
             7.3789E-01
 GRADIENT:   1.1542E+02 -5.2410E+00 -2.3636E-01 -5.8922E+00  2.1602E+01 -1.4467E+01  2.4537E+00  2.5317E+00 -2.7279E-01 -6.0402E+00
            -6.0053E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1493.61619799720        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0233E+00  8.7267E-01  8.9397E-01  1.1510E+00  8.8197E-01  9.1505E-01  8.5247E-01  2.6955E-01  9.9381E-01  9.0644E-01
             1.9176E+00
 PARAMETER:  1.2304E-01 -3.6203E-02 -1.2081E-02  2.4063E-01 -2.5597E-02  1.1218E-02 -5.9614E-02 -1.2110E+00  9.3790E-02  1.7642E-03
             7.5110E-01
 GRADIENT:   6.5975E+01  1.1081E+00 -1.7818E+01  1.8290E+01  1.9541E+01 -1.6586E+00 -3.7087E+00  5.3833E-01 -7.6221E+00  6.7004E+00
            -3.7685E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1497.31828624717        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0008E+00  7.8993E-01  8.8071E-01  1.1888E+00  8.2355E-01  9.1070E-01  1.1431E+00  2.4177E-01  9.3918E-01  7.1748E-01
             2.0788E+00
 PARAMETER:  1.0081E-01 -1.3581E-01 -2.7026E-02  2.7297E-01 -9.4135E-02  6.4542E-03  2.3370E-01 -1.3198E+00  3.7256E-02 -2.3200E-01
             8.3181E-01
 GRADIENT:   3.6839E+00  1.9487E+00 -4.6966E-01  2.9501E+00  3.6795E-01 -2.2870E-01  3.9458E-01  3.2508E-01 -6.5562E-02  7.3410E-01
            -8.4340E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1497.73110113021        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.9817E-01  6.4284E-01  8.5878E-01  1.2679E+00  7.5558E-01  9.1104E-01  1.3155E+00  7.7016E-02  8.8992E-01  7.1271E-01
             2.0642E+00
 PARAMETER:  9.8164E-02 -3.4185E-01 -5.2239E-02  3.3733E-01 -1.8026E-01  6.8367E-03  3.7418E-01 -2.4637E+00 -1.6626E-02 -2.3869E-01
             8.2473E-01
 GRADIENT:   3.1479E-01  7.5867E-01  1.4265E+00 -1.0496E+00 -2.3200E+00  1.5615E-01 -1.3240E-03  3.1132E-02 -1.8510E-01  7.6720E-02
            -1.8854E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1497.89129987945        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.9722E-01  5.7680E-01  8.7031E-01  1.3064E+00  7.4079E-01  9.1004E-01  1.3903E+00  3.7300E-02  8.7385E-01  7.2645E-01
             2.0680E+00
 PARAMETER:  9.7212E-02 -4.5026E-01 -3.8906E-02  3.6728E-01 -2.0003E-01  5.7321E-03  4.2952E-01 -3.1888E+00 -3.4841E-02 -2.1959E-01
             8.2659E-01
 GRADIENT:   1.2941E-01  3.3356E-01  7.4027E-01 -1.2346E+00 -1.2148E+00  1.0611E-01  7.9745E-02  7.0917E-03  1.1099E-01  1.8630E-01
            -2.5379E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1498.59364090704        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  9.9783E-01  3.1721E-01  9.4620E-01  1.4753E+00  7.0865E-01  9.0804E-01  1.8382E+00  1.0000E-02  8.1189E-01  7.8636E-01
             2.0741E+00
 PARAMETER:  9.7830E-02 -1.0482E+00  4.4698E-02  4.8886E-01 -2.4440E-01  3.5324E-03  7.0877E-01 -8.1150E+00 -1.0839E-01 -1.4034E-01
             8.2954E-01
 GRADIENT:   5.2096E+00  3.3002E+00  4.7739E+00  1.3824E+01 -6.1173E+00  2.2845E-01 -9.1143E-02  0.0000E+00 -8.6442E-01  5.6123E-01
             2.9504E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1499.41214075656        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      712
 NPARAMETR:  9.9381E-01  1.1284E-01  8.1868E-01  1.5612E+00  5.9681E-01  9.0559E-01  3.2780E+00  1.0000E-02  7.8487E-01  7.7680E-01
             2.0384E+00
 PARAMETER:  9.3796E-02 -2.0818E+00 -1.0006E-01  5.4546E-01 -4.1615E-01  8.3089E-04  1.2872E+00 -1.9760E+01 -1.4223E-01 -1.5257E-01
             8.1217E-01
 GRADIENT:   6.4115E+00  9.0853E-01  3.6581E-01  1.0770E+01 -3.6850E+00 -6.3183E-01 -1.2418E-01  0.0000E+00  9.3862E-01  1.3137E+00
             2.2475E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1499.73026724082        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  9.8868E-01  4.0999E-02  8.2564E-01  1.5958E+00  5.8689E-01  9.0674E-01  4.9087E+00  1.0000E-02  7.7159E-01  7.7848E-01
             2.0311E+00
 PARAMETER:  8.8615E-02 -3.0942E+00 -9.1595E-02  5.6737E-01 -4.3292E-01  2.0974E-03  1.6910E+00 -3.1767E+01 -1.5930E-01 -1.5042E-01
             8.0859E-01
 GRADIENT:  -7.7206E-01  2.1281E-01  1.9493E+00  4.7746E+00 -3.2482E+00  3.8427E-01 -1.3145E-01  0.0000E+00 -3.4256E-02  2.9583E-02
             4.0107E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1499.86689511088        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  9.8760E-01  1.0000E-02  8.1026E-01  1.6042E+00  5.7472E-01  9.0580E-01  9.5355E+00  1.0000E-02  7.6839E-01  7.7946E-01
             2.0227E+00
 PARAMETER:  8.7526E-02 -4.6475E+00 -1.1040E-01  5.7265E-01 -4.5387E-01  1.0680E-03  2.3550E+00 -5.0884E+01 -1.6346E-01 -1.4915E-01
             8.0441E-01
 GRADIENT:  -6.1073E-01  0.0000E+00 -1.1640E+00 -3.3578E+00  2.1503E+00  1.2393E-01 -3.3866E-02  0.0000E+00  2.7223E-01 -5.4504E-02
            -3.1364E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1499.87435384805        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1241
 NPARAMETR:  9.8790E-01  1.0000E-02  8.0844E-01  1.6053E+00  5.7315E-01  9.0555E-01  1.1221E+01  1.0000E-02  7.6761E-01  7.7888E-01
             2.0228E+00
 PARAMETER:  8.7831E-02 -5.0180E+00 -1.1265E-01  5.7333E-01 -4.5661E-01  7.8653E-04  2.5178E+00 -5.5456E+01 -1.6448E-01 -1.4989E-01
             8.0446E-01
 GRADIENT:   3.3848E-04  0.0000E+00  3.9176E-02 -9.7591E-02 -7.9167E-02 -2.8151E-03 -4.2599E-03  0.0000E+00  1.8668E-02  2.1112E-02
             2.1768E-02

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1499.87437743170        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1341
 NPARAMETR:  9.8793E-01  1.0000E-02  8.0913E-01  1.6055E+00  5.7352E-01  9.0546E-01  1.1254E+01  1.0000E-02  7.6749E-01  7.7883E-01
             2.0229E+00
 PARAMETER:  8.7827E-02 -5.0300E+00 -1.1181E-01  5.7345E-01 -4.5595E-01  7.8302E-04  2.5229E+00 -5.5603E+01 -1.6459E-01 -1.4987E-01
             8.0454E-01
 GRADIENT:  -1.2106E-02  0.0000E+00 -2.7037E-03  1.4777E-02  4.1237E-03  5.6513E-03  3.1116E-04  0.0000E+00  2.2878E-03  1.9086E-03
             2.3759E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1341
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -7.5629E-05  8.3780E-04  6.6166E-06 -9.7126E-03 -1.8308E-02
 SE:             2.9302E-02  1.9064E-03  1.9721E-04  2.7682E-02  2.1544E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9794E-01  6.6032E-01  9.7324E-01  7.2570E-01  3.9544E-01

 ETASHRINKSD(%)  1.8334E+00  9.3613E+01  9.9339E+01  7.2602E+00  2.7825E+01
 ETASHRINKVR(%)  3.6331E+00  9.9592E+01  9.9996E+01  1.3993E+01  4.7908E+01
 EBVSHRINKSD(%)  1.8598E+00  9.4134E+01  9.9310E+01  7.0838E+00  2.7561E+01
 EBVSHRINKVR(%)  3.6850E+00  9.9656E+01  9.9995E+01  1.3666E+01  4.7526E+01
 RELATIVEINF(%)  8.3616E+01  1.0307E-02  2.3192E-04  4.2495E+00  1.6465E+00
 EPSSHRINKSD(%)  3.4445E+01
 EPSSHRINKVR(%)  5.7026E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1499.8743774316999     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -764.72355086796176     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.06
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1499.874       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  1.00E-02  8.09E-01  1.61E+00  5.74E-01  9.06E-01  1.13E+01  1.00E-02  7.68E-01  7.79E-01  2.02E+00
 


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
+        1.35E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        2.99E+00  0.00E+00  8.28E+02
 
 TH 4
+       -3.69E+01  0.00E+00 -8.04E+01  6.32E+02
 
 TH 5
+        3.41E+01  0.00E+00 -1.46E+03 -1.24E+02  2.87E+03
 
 TH 6
+        1.60E-01  0.00E+00 -2.60E+00 -9.64E+00  5.15E-01  2.44E+02
 
 TH 7
+        1.23E-02  0.00E+00  1.08E-02 -1.65E-02 -7.43E-03 -2.58E-02  3.75E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.13E+00  0.00E+00  3.88E+00 -9.54E+00 -2.93E+00 -5.91E+00  3.44E-02  0.00E+00  2.60E+02
 
 TH10
+       -2.01E+01  0.00E+00  1.17E+01 -1.66E+00 -5.78E+01  2.60E+00 -2.01E-02  0.00E+00 -1.64E+00  9.05E+01
 
 TH11
+       -1.47E+01  0.00E+00 -1.69E+01 -8.35E+00 -5.44E+00  5.02E+00  5.57E-03  0.00E+00  1.22E+01  3.08E+01  6.49E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.202
Stop Time:
Sat Sep 25 08:33:26 CDT 2021
