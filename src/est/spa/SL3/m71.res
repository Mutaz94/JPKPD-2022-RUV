Sat Sep 25 11:52:04 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat71.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1556.99358302184        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3702E+02 -7.8589E+01  1.0612E+01 -1.2089E+02  7.7152E+00  1.5440E+00 -2.0433E+01 -7.9889E+00 -1.7613E+01 -4.9136E+00
            -9.6249E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1574.51066879873        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.6792E-01  1.1156E+00  9.9049E-01  1.0021E+00  1.0363E+00  9.6658E-01  1.2190E+00  1.0712E+00  1.0045E+00  9.5323E-01
             1.2835E+00
 PARAMETER:  6.7393E-02  2.0939E-01  9.0444E-02  1.0205E-01  1.3570E-01  6.6010E-02  2.9803E-01  1.6874E-01  1.0449E-01  5.2102E-02
             3.4957E-01
 GRADIENT:   4.5914E+01 -5.5791E+00  4.3941E+00 -1.3317E+01  2.6935E+00 -8.0244E+00  3.4619E+00 -2.6549E+00  4.3851E+00 -4.6390E-01
             2.2607E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1574.87232427330        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  9.6863E-01  1.0924E+00  9.4137E-01  1.0165E+00  9.9223E-01  9.9028E-01  1.2694E+00  1.1747E+00  9.7843E-01  8.8397E-01
             1.2718E+00
 PARAMETER:  6.8123E-02  1.8838E-01  3.9582E-02  1.1635E-01  9.2203E-02  9.0237E-02  3.3852E-01  2.6099E-01  7.8196E-02 -2.3331E-02
             3.4042E-01
 GRADIENT:   4.6048E+01 -4.7484E+00 -1.9294E+00 -7.5850E+00  1.3007E+00  1.4951E+00  6.1765E+00  3.5473E+00  4.9628E+00  5.5359E-01
             2.0422E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1575.77497046532        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.4971E-01  1.1391E+00  7.5905E-01  9.8268E-01  9.0426E-01  9.8420E-01  1.2074E+00  9.3432E-01  9.6719E-01  7.8833E-01
             1.2106E+00
 PARAMETER:  4.8399E-02  2.3023E-01 -1.7568E-01  8.2529E-02 -6.4001E-04  8.4078E-02  2.8850E-01  3.2066E-02  6.6635E-02 -1.3784E-01
             2.9115E-01
 GRADIENT:   9.5093E-01  6.1962E+00  2.2783E+00  3.5036E+00 -3.5464E+00 -1.4060E+00  6.2590E-01  3.0345E-01  1.4841E+00 -3.1438E-01
             9.6053E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1575.97691390411        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  9.6111E-01  1.1667E+00  6.9670E-01  9.6037E-01  8.8333E-01  9.9704E-01  1.1993E+00  8.5155E-01  9.5610E-01  7.6793E-01
             1.2054E+00
 PARAMETER:  6.0333E-02  2.5415E-01 -2.6141E-01  5.9563E-02 -2.4052E-02  9.7031E-02  2.8170E-01 -6.0701E-02  5.5104E-02 -1.6406E-01
             2.8679E-01
 GRADIENT:   1.9072E+00 -4.6529E-01 -4.0718E-01  1.2752E+00  4.1818E-01  8.0971E-01  4.4721E-01 -3.6787E-02 -5.6747E-01  2.2843E-01
            -9.7633E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1576.13385212720        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  9.5941E-01  1.3284E+00  5.1815E-01  8.4263E-01  8.4609E-01  9.9367E-01  1.0794E+00  6.8537E-01  1.0141E+00  6.7980E-01
             1.2108E+00
 PARAMETER:  5.8560E-02  3.8395E-01 -5.5749E-01 -7.1221E-02 -6.7131E-02  9.3655E-02  1.7639E-01 -2.7779E-01  1.1396E-01 -2.8595E-01
             2.9130E-01
 GRADIENT:  -3.7946E+00  4.2265E+00  2.6950E+00 -2.4485E-01 -5.1283E+00 -7.1258E-01  3.1311E-01 -1.5113E-01  1.7054E-01 -6.3225E-02
             7.1184E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1576.28744576651        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  9.5883E-01  1.4530E+00  3.4500E-01  7.3362E-01  7.8969E-01  9.9070E-01  9.7483E-01  5.0972E-01  1.0693E+00  5.5843E-01
             1.2094E+00
 PARAMETER:  5.7962E-02  4.7362E-01 -9.6422E-01 -2.0976E-01 -1.3611E-01  9.0653E-02  7.4503E-02 -5.7389E-01  1.6702E-01 -4.8262E-01
             2.9010E-01
 GRADIENT:  -3.9004E-01 -1.4763E+00 -8.1181E-01  1.9251E+00  1.8813E+00 -1.2009E+00 -8.3782E-01  3.5123E-01  8.3497E-01  3.2847E-02
             1.5420E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1576.37095641128        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  9.5733E-01  1.5275E+00  2.8287E-01  6.7290E-01  7.8404E-01  9.9172E-01  9.3190E-01  3.8259E-01  1.1158E+00  5.2768E-01
             1.2046E+00
 PARAMETER:  5.6391E-02  5.2363E-01 -1.1628E+00 -2.9616E-01 -1.4330E-01  9.1687E-02  2.9467E-02 -8.6080E-01  2.0961E-01 -5.3926E-01
             2.8611E-01
 GRADIENT:   9.6662E-01 -2.0060E+00 -8.7349E-01  2.6116E-01  1.4662E+00 -1.7276E-01  4.9142E-02  2.5893E-01  6.2974E-01 -6.9250E-02
             7.6553E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1576.41974154180        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1081
 NPARAMETR:  9.5731E-01  1.5575E+00  2.7503E-01  6.5743E-01  7.9360E-01  9.9160E-01  9.2161E-01  2.5463E-01  1.1310E+00  5.4365E-01
             1.2023E+00
 PARAMETER:  5.6367E-02  5.4309E-01 -1.1909E+00 -3.1942E-01 -1.3117E-01  9.1563E-02  1.8366E-02 -1.2679E+00  2.2308E-01 -5.0946E-01
             2.8422E-01
 GRADIENT:   1.1299E+00  5.1447E+00  3.3437E-01  1.5767E+00 -1.8656E+00 -1.9954E-01  1.1594E-01  8.0209E-02 -7.3401E-03 -1.8041E-01
            -1.8923E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1576.46200364855        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1256
 NPARAMETR:  9.5662E-01  1.5428E+00  2.7327E-01  6.6219E-01  7.8535E-01  9.9203E-01  9.2461E-01  6.5778E-02  1.1250E+00  5.4019E-01
             1.2018E+00
 PARAMETER:  5.5647E-02  5.3359E-01 -1.1973E+00 -3.1220E-01 -1.4163E-01  9.2000E-02  2.1612E-02 -2.6215E+00  2.1782E-01 -5.1584E-01
             2.8381E-01
 GRADIENT:  -1.5430E-02  3.8614E-02 -1.6741E-02 -1.7147E-02 -7.1401E-02 -2.0684E-03 -2.8300E-02  3.9214E-03  8.3630E-03 -1.1449E-04
            -1.6475E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1576.46390711388        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1431
 NPARAMETR:  9.5665E-01  1.5419E+00  2.7383E-01  6.6289E-01  7.8532E-01  9.9206E-01  9.2529E-01  1.0000E-02  1.1245E+00  5.4088E-01
             1.2018E+00
 PARAMETER:  5.5686E-02  5.3300E-01 -1.1953E+00 -3.1115E-01 -1.4166E-01  9.2029E-02  2.2347E-02 -4.9210E+00  2.1731E-01 -5.1455E-01
             2.8383E-01
 GRADIENT:   5.9076E-03  1.1503E-02 -6.3275E-03  1.6975E-02  1.4631E-02 -8.6588E-04  6.3475E-03  0.0000E+00  2.3497E-03 -1.5147E-03
            -3.3275E-03

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1576.46390927043        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1488
 NPARAMETR:  9.5665E-01  1.5415E+00  2.7389E-01  6.6307E-01  7.8517E-01  9.9206E-01  9.2537E-01  1.0000E-02  1.1242E+00  5.4077E-01
             1.2018E+00
 PARAMETER:  5.5681E-02  5.3278E-01 -1.1950E+00 -3.1088E-01 -1.4186E-01  9.2033E-02  2.2441E-02 -4.9090E+00  2.1711E-01 -5.1477E-01
             2.8385E-01
 GRADIENT:  -4.1528E-03 -7.6497E-03  4.1307E-03 -1.0991E-02 -8.0853E-03  2.3731E-04 -4.9640E-03  0.0000E+00 -1.7428E-03  3.9709E-04
             3.3644E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1488
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3108E-04 -1.2000E-02 -2.0866E-04  1.1402E-02 -2.0964E-02
 SE:             2.9800E-02  2.6650E-02  1.1749E-04  2.4893E-02  1.5745E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9114E-01  6.5250E-01  7.5740E-02  6.4692E-01  1.8302E-01

 ETASHRINKSD(%)  1.6579E-01  1.0718E+01  9.9606E+01  1.6605E+01  4.7253E+01
 ETASHRINKVR(%)  3.3131E-01  2.0288E+01  9.9998E+01  3.0454E+01  7.2177E+01
 EBVSHRINKSD(%)  5.8581E-01  1.0890E+01  9.9629E+01  1.6344E+01  4.8215E+01
 EBVSHRINKVR(%)  1.1682E+00  2.0594E+01  9.9999E+01  3.0017E+01  7.3183E+01
 RELATIVEINF(%)  9.7192E+01  6.7785E+00  1.1854E-04  4.7497E+00  2.5490E+00
 EPSSHRINKSD(%)  4.3569E+01
 EPSSHRINKVR(%)  6.8155E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1576.4639092704283     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -841.31308270669012     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1576.464       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.57E-01  1.54E+00  2.74E-01  6.63E-01  7.85E-01  9.92E-01  9.25E-01  1.00E-02  1.12E+00  5.41E-01  1.20E+00
 


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
+        1.22E+03
 
 TH 2
+       -9.53E+00  4.38E+02
 
 TH 3
+       -5.09E+01  4.07E+02  2.10E+03
 
 TH 4
+       -1.18E+01  1.77E+02 -1.24E+03  1.41E+03
 
 TH 5
+       -1.78E+01 -5.09E+02 -1.62E+03  8.45E+02  1.66E+03
 
 TH 6
+       -4.35E-01 -1.50E+00 -3.93E+00 -5.88E+00 -2.58E+00  2.00E+02
 
 TH 7
+       -1.63E-01  1.13E+01 -1.09E+02  1.08E+00  1.92E+01  1.38E+00  1.48E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.25E+00 -1.20E+01 -1.92E+01  5.01E+01 -3.05E+01  4.94E-02  1.24E+01  0.00E+00  7.73E+01
 
 TH10
+       -2.94E-01 -1.33E+01 -3.05E+01 -3.05E+01 -9.44E+01  2.00E-01  2.18E+01  0.00E+00  1.72E+01  7.02E+01
 
 TH11
+       -9.07E+00 -1.11E+01 -2.44E+01 -6.49E+00 -2.80E+01  2.59E+00  8.43E+00  0.00E+00  1.25E+01  2.19E+01  1.46E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.850
Stop Time:
Sat Sep 25 11:52:28 CDT 2021
