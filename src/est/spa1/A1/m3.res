Wed Sep 29 21:46:29 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat3.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1361.34033608525        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9509E+02  2.1814E+01  1.0211E+01  5.7524E+01  1.6519E+02  5.9727E+01 -3.0690E+01 -1.6937E+01 -1.1121E+01 -6.5173E+01
            -1.4028E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1775.01201687150        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0548E+00  9.7291E-01  1.0322E+00  1.0711E+00  9.4987E-01  1.0205E+00  9.5371E-01  8.2786E-01  8.7375E-01  9.0972E-01
             2.1158E+00
 PARAMETER:  1.5333E-01  7.2539E-02  1.3172E-01  1.6868E-01  4.8572E-02  1.2030E-01  5.2605E-02 -8.8912E-02 -3.4958E-02  5.3847E-03
             8.4944E-01
 GRADIENT:   1.9477E+02  2.7245E+00 -8.4934E+00  2.6687E+01  4.0089E+01  2.8359E+01 -4.0033E+00  2.0569E+00 -1.2481E+01 -7.3337E+00
            -6.8799E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1788.33908045754        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0385E+00  8.0166E-01  6.2815E-01  1.1562E+00  6.5239E-01  9.3737E-01  9.4765E-01  1.8605E-01  9.8670E-01  7.3798E-01
             2.0903E+00
 PARAMETER:  1.3782E-01 -1.2107E-01 -3.6497E-01  2.4515E-01 -3.2711E-01  3.5328E-02  4.6232E-02 -1.5817E+00  8.6614E-02 -2.0384E-01
             8.3730E-01
 GRADIENT:   1.4443E+02  2.1751E+01 -3.0884E+01  1.2550E+02  8.9227E+01 -6.0149E+00 -9.6331E+00  4.6119E-01  2.8243E-01 -2.4300E+01
            -5.2727E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1793.69145285579        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      336
 NPARAMETR:  1.0246E+00  7.0038E-01  6.6138E-01  1.1848E+00  6.3609E-01  9.7124E-01  1.2089E+00  1.0704E-01  9.0368E-01  8.1854E-01
             2.2149E+00
 PARAMETER:  1.2430E-01 -2.5614E-01 -3.1342E-01  2.6955E-01 -3.5241E-01  7.0814E-02  2.8969E-01 -2.1346E+00 -1.2818E-03 -1.0024E-01
             8.9522E-01
 GRADIENT:  -2.1359E+01  2.5886E+00 -2.7169E+01  1.4251E+01  5.1369E+01 -2.1863E+00 -2.7976E+00  1.8901E-01 -7.7238E+00 -5.7620E+00
             4.9215E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1803.23181478656        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  1.0396E+00  2.9791E-01  6.3947E-01  1.3828E+00  4.9968E-01  9.6964E-01  1.4815E+00  2.8528E-02  8.5035E-01  7.7900E-01
             2.0814E+00
 PARAMETER:  1.3884E-01 -1.1110E+00 -3.4711E-01  4.2413E-01 -5.9379E-01  6.9165E-02  4.9304E-01 -3.4569E+00 -6.2104E-02 -1.4975E-01
             8.3304E-01
 GRADIENT:   3.3281E+01  7.3890E+00  1.0334E+01  2.2755E+01 -1.5689E+01 -1.1957E+00 -2.8071E+00  1.0016E-02 -9.6962E+00 -8.1699E+00
            -5.3957E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1804.91280888745        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      691
 NPARAMETR:  1.0219E+00  2.4153E-01  6.4450E-01  1.4042E+00  4.9336E-01  9.7171E-01  2.2624E+00  1.1601E-02  8.5668E-01  8.3788E-01
             2.0557E+00
 PARAMETER:  1.2171E-01 -1.3208E+00 -3.3928E-01  4.3949E-01 -6.0651E-01  7.1302E-02  9.1642E-01 -4.3567E+00 -5.4688E-02 -7.6880E-02
             8.2064E-01
 GRADIENT:  -1.6174E-01  5.1822E+00  6.5480E+00  9.0034E+00 -1.5498E+01  1.2411E+00  6.5799E-01  2.0319E-03  3.1675E-03  3.4449E+00
            -8.5005E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1806.12736458045        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      870
 NPARAMETR:  1.0188E+00  1.0845E-01  6.7878E-01  1.4765E+00  4.9218E-01  9.6562E-01  3.4866E+00  1.0000E-02  8.2229E-01  8.3832E-01
             2.0788E+00
 PARAMETER:  1.1858E-01 -2.1214E+00 -2.8746E-01  4.8966E-01 -6.0892E-01  6.5011E-02  1.3489E+00 -5.3981E+00 -9.5657E-02 -7.6353E-02
             8.3180E-01
 GRADIENT:   3.5045E+00  1.3833E+00  2.8806E+00  5.9841E+00 -6.0247E+00  2.2068E-01  2.3627E-01  0.0000E+00 -3.1442E+00  3.3803E-01
            -1.1901E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1806.45363529440        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1046
 NPARAMETR:  1.0135E+00  3.3233E-02  7.0216E-01  1.5229E+00  4.9262E-01  9.6040E-01  5.2217E+00  1.0000E-02  8.1852E-01  8.4958E-01
             2.0805E+00
 PARAMETER:  1.1340E-01 -3.3042E+00 -2.5359E-01  5.2060E-01 -6.0802E-01  5.9592E-02  1.7528E+00 -6.9353E+00 -1.0026E-01 -6.3008E-02
             8.3262E-01
 GRADIENT:  -2.0422E+00  3.5003E-01  2.4158E+00  1.2712E+01 -5.2756E+00 -1.1252E+00 -1.1553E-01  0.0000E+00 -4.5307E-02  2.5650E-02
            -1.3507E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1806.59085121602        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1221
 NPARAMETR:  1.0133E+00  1.0000E-02  7.0224E-01  1.5276E+00  4.9055E-01  9.6215E-01  7.5795E+00  1.0000E-02  8.1472E-01  8.4970E-01
             2.0833E+00
 PARAMETER:  1.1321E-01 -4.6683E+00 -2.5348E-01  5.2372E-01 -6.1223E-01  6.1412E-02  2.1254E+00 -8.6423E+00 -1.0491E-01 -6.2870E-02
             8.3395E-01
 GRADIENT:  -1.1031E-01  0.0000E+00 -5.4883E-01 -2.8434E-01  8.6248E-01 -1.5526E-01 -3.5293E-02  0.0000E+00  5.5521E-02 -8.9682E-02
            -9.7700E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1806.59779946580        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1407
 NPARAMETR:  1.0136E+00  1.0000E-02  7.0233E-01  1.5245E+00  4.9016E-01  9.6266E-01  8.7691E+00  1.0000E-02  8.1455E-01  8.4994E-01
             2.0833E+00
 PARAMETER:  1.1352E-01 -4.6683E+00 -2.5336E-01  5.2165E-01 -6.1303E-01  6.1943E-02  2.2712E+00 -8.6423E+00 -1.0512E-01 -6.2587E-02
             8.3394E-01
 GRADIENT:   7.3054E-01  0.0000E+00  1.2853E+00 -6.5080E+00 -4.5676E-01  8.4099E-02 -3.5134E-02  0.0000E+00  9.2814E-02 -2.0138E-02
             5.6743E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1806.60513544300        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1549
 NPARAMETR:  1.0135E+00  1.0000E-02  7.0184E-01  1.5261E+00  4.9017E-01  9.6257E-01  1.1018E+01  1.0000E-02  8.1441E-01  8.4995E-01
             2.0832E+00
 PARAMETER:  1.1341E-01 -4.6683E+00 -2.5405E-01  5.2273E-01 -6.1299E-01  6.1852E-02  2.4995E+00 -8.6423E+00 -1.0530E-01 -6.2581E-02
             8.3392E-01
 GRADIENT:   4.1592E-01  0.0000E+00  1.7658E-01 -3.0185E+00  2.8128E-01  4.3870E-02  6.6693E-03  0.0000E+00  1.2609E-01  3.8266E-02
             4.7827E-03

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1806.60513544300        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1571
 NPARAMETR:  1.0135E+00  1.0000E-02  7.0184E-01  1.5261E+00  4.9017E-01  9.6257E-01  1.1018E+01  1.0000E-02  8.1441E-01  8.4995E-01
             2.0832E+00
 PARAMETER:  1.1341E-01 -4.6683E+00 -2.5405E-01  5.2273E-01 -6.1299E-01  6.1852E-02  2.4995E+00 -8.6423E+00 -1.0530E-01 -6.2581E-02
             8.3392E-01
 GRADIENT:   4.1592E-01  0.0000E+00  1.7658E-01 -3.0185E+00  2.8128E-01  4.3870E-02  6.6693E-03  0.0000E+00  1.2609E-01  3.8266E-02
             4.7827E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1571
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.9687E-04  5.5398E-04 -2.1522E-05 -6.8878E-03 -1.3012E-02
 SE:             2.9496E-02  1.7287E-03  1.9944E-04  2.8229E-02  2.3143E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8386E-01  7.4862E-01  9.1407E-01  8.0723E-01  5.7396E-01

 ETASHRINKSD(%)  1.1862E+00  9.4209E+01  9.9332E+01  5.4294E+00  2.2469E+01
 ETASHRINKVR(%)  2.3584E+00  9.9665E+01  9.9996E+01  1.0564E+01  3.9889E+01
 EBVSHRINKSD(%)  1.3810E+00  9.4490E+01  9.9315E+01  4.8968E+00  2.1705E+01
 EBVSHRINKVR(%)  2.7429E+00  9.9696E+01  9.9995E+01  9.5538E+00  3.8698E+01
 RELATIVEINF(%)  8.6147E+01  1.4712E-02  3.3622E-04  6.2429E+00  4.0525E+00
 EPSSHRINKSD(%)  2.7892E+01
 EPSSHRINKVR(%)  4.8004E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1806.6051354430042     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -887.66660223833151     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.19
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1806.605       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.00E-02  7.02E-01  1.53E+00  4.90E-01  9.63E-01  1.10E+01  1.00E-02  8.14E-01  8.50E-01  2.08E+00
 


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
+        1.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.56E+01  0.00E+00  1.07E+03
 
 TH 4
+       -1.91E+01  0.00E+00 -1.57E+02  6.52E+02
 
 TH 5
+        3.59E+01  0.00E+00 -1.87E+03 -1.14E+02  3.96E+03
 
 TH 6
+        5.78E-01  0.00E+00  2.56E+00 -5.36E+00  9.58E-01  2.03E+02
 
 TH 7
+        1.39E-03  0.00E+00  6.50E-03 -7.03E-03 -6.69E-03  2.69E-03  2.27E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.84E+00  0.00E+00  3.04E+01 -4.98E+00 -5.34E+00  2.85E+00  2.96E-02  0.00E+00  2.41E+02
 
 TH10
+        9.43E-02  0.00E+00 -2.78E+01  5.86E+00 -6.84E+01  1.32E+00  1.94E-03  0.00E+00 -2.64E+00  1.17E+02
 
 TH11
+       -1.35E+01  0.00E+00 -9.92E+00 -1.25E+01  1.12E+01  9.28E-01  3.03E-03  0.00E+00  9.22E+00  1.95E+01  1.06E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.057
Stop Time:
Wed Sep 29 21:47:02 CDT 2021
