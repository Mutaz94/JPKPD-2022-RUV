Wed Sep 29 12:31:44 CDT 2021
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
$DATA ../../../../data/spa/A2/dat3.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1019.18965975403        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8563E+02  3.1000E+00  2.2721E+01  9.1865E+00  1.2519E+02  5.1199E+01 -1.8009E+01 -1.8907E+01 -1.5801E+01 -5.8481E+01
            -1.1796E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1427.82074681209        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0826E+00  1.0373E+00  9.8466E-01  1.0693E+00  9.3629E-01  9.0479E-01  9.2195E-01  9.9821E-01  9.1552E-01  8.5637E-01
             2.6628E+00
 PARAMETER:  1.7933E-01  1.3657E-01  8.4546E-02  1.6699E-01  3.4171E-02 -5.0573E-05  1.8732E-02  9.8210E-02  1.1734E-02 -5.5058E-02
             1.0794E+00
 GRADIENT:   2.1753E+02  2.1808E+01  2.7196E-01  3.4812E+01 -5.6811E+00 -2.0083E+01  4.7595E+00  1.3053E+00  8.8590E-01  1.2444E+01
             8.9681E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1431.98617469793        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0692E+00  1.1355E+00  7.3243E-01  1.0038E+00  8.5483E-01  9.6283E-01  8.4883E-01  1.1056E+00  1.0906E+00  3.7594E-01
             2.7110E+00
 PARAMETER:  1.6691E-01  2.2704E-01 -2.1139E-01  1.0376E-01 -5.6855E-02  6.2119E-02 -6.3900E-02  2.0039E-01  1.8671E-01 -8.7832E-01
             1.0973E+00
 GRADIENT:   1.5792E+02  2.0270E+01 -1.1343E+01  4.6842E+01  1.1777E+01  3.3733E+00  6.6783E+00  5.1524E+00  2.0491E+01  2.2204E+00
             1.4698E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1435.14619910106        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      277
 NPARAMETR:  1.0450E+00  1.0629E+00  7.1844E-01  1.0104E+00  8.1481E-01  9.6693E-01  8.8125E-01  9.0184E-01  9.3373E-01  3.5491E-01
             2.7364E+00
 PARAMETER:  1.4405E-01  1.6097E-01 -2.3067E-01  1.1033E-01 -1.0480E-01  6.6373E-02 -2.6409E-02 -3.3148E-03  3.1427E-02 -9.3590E-01
             1.1066E+00
 GRADIENT:   1.3854E+01 -9.6423E+00 -6.6431E+00 -2.1111E+00  9.1405E+00  2.5336E+00  2.4700E+00  1.4114E+00  2.7961E+00  1.7297E+00
             7.9982E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1436.32576966675        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      454
 NPARAMETR:  1.0385E+00  1.3426E+00  7.1465E-01  8.4781E-01  9.3714E-01  9.6001E-01  7.1945E-01  1.2342E+00  1.0457E+00  8.7401E-02
             2.7570E+00
 PARAMETER:  1.3780E-01  3.9463E-01 -2.3597E-01 -6.5094E-02  3.5079E-02  5.9187E-02 -2.2928E-01  3.1042E-01  1.4469E-01 -2.3372E+00
             1.1141E+00
 GRADIENT:  -1.1119E+00  1.0826E+01  5.2711E+00  6.4973E+00 -1.1890E+01  4.7764E-01 -1.7429E+00  1.2064E-01  1.5641E+00  7.4856E-02
            -7.4647E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1437.54459055114        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      629
 NPARAMETR:  1.0361E+00  1.6669E+00  4.6750E-01  6.3108E-01  9.6819E-01  9.5113E-01  6.7703E-01  1.3799E+00  1.1035E+00  1.0000E-02
             2.7525E+00
 PARAMETER:  1.3551E-01  6.1098E-01 -6.6035E-01 -3.6031E-01  6.7669E-02  4.9891E-02 -2.9004E-01  4.2198E-01  1.9845E-01 -4.6253E+00
             1.1125E+00
 GRADIENT:  -8.0846E+00  3.1171E+01  9.6505E+00  4.9325E+00 -2.0829E+01 -2.0011E+00 -1.9265E+00 -2.1985E+00 -3.1873E+00  0.0000E+00
            -1.7124E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1448.26735262043        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      814
 NPARAMETR:  1.0015E+00  1.9998E+00  1.1412E-01  3.5021E-01  9.0785E-01  9.2563E-01  6.1289E-01  1.5757E+00  1.2580E+00  1.0000E-02
             2.6349E+00
 PARAMETER:  1.0155E-01  7.9305E-01 -2.0705E+00 -9.4923E-01  3.3214E-03  2.2716E-02 -3.8957E-01  5.5471E-01  3.2951E-01 -9.9913E+00
             1.0689E+00
 GRADIENT:  -6.4982E+01  5.5127E+01  1.3209E+01 -1.0627E+01 -3.1587E+01 -1.0955E+01 -6.3724E+00 -5.6852E+00 -1.1116E+01  0.0000E+00
             6.4271E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1453.70577190775        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      989
 NPARAMETR:  1.0286E+00  2.2204E+00  4.4342E-02  2.2199E-01  1.0415E+00  9.5728E-01  6.0418E-01  9.8882E-01  2.1916E+00  1.0000E-02
             2.6050E+00
 PARAMETER:  1.2816E-01  8.9770E-01 -3.0158E+00 -1.4051E+00  1.4067E-01  5.6341E-02 -4.0388E-01  8.8756E-02  8.8461E-01 -1.3833E+01
             1.0574E+00
 GRADIENT:   8.0827E+00  1.7487E+01 -7.2465E+00  1.7703E+01  1.1393E+01  2.5924E+00  2.6639E+00  2.0131E+00 -2.0309E+00  0.0000E+00
             2.5704E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1455.65859659019        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1168
 NPARAMETR:  1.0238E+00  2.2811E+00  4.1270E-02  1.8675E-01  1.0935E+00  9.4781E-01  5.9307E-01  3.0581E-01  2.7012E+00  1.0000E-02
             2.5769E+00
 PARAMETER:  1.2352E-01  9.2467E-01 -3.0876E+00 -1.5780E+00  1.8936E-01  4.6403E-02 -4.2244E-01 -1.0848E+00  1.0937E+00 -1.4287E+01
             1.0466E+00
 GRADIENT:  -5.9539E+00  2.4768E+01 -3.6402E-01  9.0655E+00  3.9563E+00 -1.3611E+00  4.5976E-01  3.3939E-01  4.2445E+00  0.0000E+00
            -6.9660E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1456.22098322496        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1349             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0261E+00  2.2762E+00  4.0004E-02  1.7364E-01  1.0973E+00  9.4915E-01  5.9280E-01  5.5503E-02  2.6407E+00  1.0000E-02
             2.5830E+00
 PARAMETER:  1.2572E-01  9.2249E-01 -3.1188E+00 -1.6508E+00  1.9282E-01  4.7807E-02 -4.2289E-01 -2.7913E+00  1.0711E+00 -1.4409E+01
             1.0490E+00
 GRADIENT:   7.9737E+01  2.1504E+02  4.4265E+00  1.0240E+01  3.8852E+00  6.2247E+00  4.1038E+00  1.1434E-02  3.5646E+00  0.0000E+00
             8.5412E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1456.22832730139        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1528
 NPARAMETR:  1.0259E+00  2.2756E+00  4.0319E-02  1.7531E-01  1.0947E+00  9.4932E-01  5.9294E-01  1.1263E-02  2.6286E+00  1.0000E-02
             2.5816E+00
 PARAMETER:  1.2558E-01  9.2225E-01 -3.1109E+00 -1.6412E+00  1.9049E-01  4.7991E-02 -4.2266E-01 -4.3863E+00  1.0664E+00 -1.4409E+01
             1.0484E+00
 GRADIENT:   1.7558E-01  7.3397E-01 -1.7195E-02  3.3194E-01  4.5004E-01  5.1273E-03  2.2088E-02  4.4264E-04  4.9816E-02  0.0000E+00
            -4.3124E-02

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1456.23036223089        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1593
 NPARAMETR:  1.0261E+00  2.2708E+00  4.0061E-02  1.7550E-01  1.0945E+00  9.4915E-01  5.9329E-01  1.0000E-02  2.5972E+00  1.0000E-02
             2.6088E+00
 PARAMETER:  1.2543E-01  9.2180E-01 -3.1106E+00 -1.6441E+00  1.8972E-01  4.7969E-02 -4.2275E-01 -4.5957E+00  1.0651E+00 -1.4409E+01
             1.0484E+00
 GRADIENT:  -3.6529E-01  3.6980E+00  2.6733E-01 -2.6610E-01 -3.2955E-01  2.8965E-02 -6.4638E-02  0.0000E+00  1.1047E+03  0.0000E+00
            -1.1513E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1593
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0389E-03 -2.0198E-02  5.1773E-05  2.1734E-02 -2.3824E-04
 SE:             2.9223E-02  2.5082E-02  3.1847E-05  1.7592E-02  1.9781E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7164E-01  4.2067E-01  1.0401E-01  2.1666E-01  2.2843E-01

 ETASHRINKSD(%)  2.1001E+00  1.5971E+01  9.9893E+01  4.1065E+01  9.9337E+01
 ETASHRINKVR(%)  4.1562E+00  2.9392E+01  1.0000E+02  6.5266E+01  9.9996E+01
 EBVSHRINKSD(%)  2.3553E+00  1.6923E+01  9.9834E+01  3.7897E+01  9.9260E+01
 EBVSHRINKVR(%)  4.6551E+00  3.0981E+01  1.0000E+02  6.1432E+01  9.9995E+01
 RELATIVEINF(%)  9.2331E+01  1.8773E+01  1.9975E-04  1.1633E+01  1.3958E-03
 EPSSHRINKSD(%)  2.8077E+01
 EPSSHRINKVR(%)  4.8270E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1456.2303622308903     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -721.07953566715207     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.15
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.01
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1456.230       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.27E+00  4.03E-02  1.75E-01  1.09E+00  9.49E-01  5.93E-01  1.00E-02  2.62E+00  1.00E-02  2.58E+00
 


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
+        1.12E+03
 
 TH 2
+       -5.43E+01  4.24E+02
 
 TH 3
+       -2.78E+02  3.57E+02  2.43E+04
 
 TH 4
+       -3.00E+01  2.27E+02 -2.92E+03  2.17E+03
 
 TH 5
+       -2.72E+01 -2.62E+02 -6.92E+02  4.27E+02  4.93E+02
 
 TH 6
+        1.13E-02 -8.91E+00 -4.96E+01 -2.46E+01 -2.75E+00  2.01E+02
 
 TH 7
+        4.27E+00 -1.06E+01 -1.86E+01  2.93E+01  2.05E+01  6.01E+00  2.73E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  2.82E-01
 
 TH 9
+        4.05E+01 -2.57E+01  4.75E+03 -5.67E+02 -9.42E+00  3.92E+01  4.68E+01  0.00E+00  3.68E+03
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        8.66E+04  5.31E+03  8.65E+04 -3.85E+04 -3.74E+00 -3.80E+01 -2.34E+01  0.00E+00 -3.89E+03  0.00E+00  4.18E+03
 
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
 #CPUT: Total CPU Time in Seconds,       28.221
Stop Time:
Wed Sep 29 12:32:26 CDT 2021
