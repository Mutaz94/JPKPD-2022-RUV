Thu Sep 30 00:19:33 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat50.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   549.006959316516        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5822E+02  7.5699E+01  1.6257E+02 -3.2456E+01  1.7896E+02  3.6548E+01 -6.9362E+01 -7.5349E+01 -5.2938E+01 -1.3798E+02
            -4.8936E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1424.80490488691        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0007E+00  1.0234E+00  9.5165E-01  1.1694E+00  9.6900E-01  8.5127E-01  9.7474E-01  9.6161E-01  9.5859E-01  9.9041E-01
             5.3494E+00
 PARAMETER:  1.0069E-01  1.2313E-01  5.0437E-02  2.5653E-01  6.8511E-02 -6.1032E-02  7.4417E-02  6.0859E-02  5.7704E-02  9.0363E-02
             1.7770E+00
 GRADIENT:  -6.3209E+01  6.1544E-01 -1.4546E+01  1.7098E+01 -1.3613E+01 -1.1264E+01  1.2656E+01  7.7533E+00  2.6175E+01  2.2919E+01
             3.3695E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1474.30607055701        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.8240E-01  7.1856E-01  2.3077E-01  1.2449E+00  3.7592E-01  9.6106E-01  4.9021E-01  1.3583E-02  1.2454E+00  3.1244E-01
             4.2573E+00
 PARAMETER:  8.2242E-02 -2.3050E-01 -1.3663E+00  3.1902E-01 -8.7837E-01  6.0285E-02 -6.1292E-01 -4.1989E+00  3.1946E-01 -1.0633E+00
             1.5486E+00
 GRADIENT:  -4.8719E+01 -2.0541E+01 -9.8547E+01  1.3376E+02  1.3199E+02  2.6354E+00  3.9558E+00  2.5347E-03  1.5490E+01  5.0958E+00
             1.8732E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1509.98371420242        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      257
 NPARAMETR:  9.7696E-01  6.4057E-01  2.6220E-01  1.1972E+00  3.4870E-01  9.3559E-01  5.2662E-01  1.0000E-02  1.1733E+00  3.5755E-01
             3.3638E+00
 PARAMETER:  7.6695E-02 -3.4539E-01 -1.2387E+00  2.7997E-01 -9.5354E-01  3.3423E-02 -5.4128E-01 -5.3009E+00  2.5985E-01 -9.2849E-01
             1.3131E+00
 GRADIENT:  -5.5037E+01  8.2264E+01  1.9014E+01  7.8482E+01 -7.4878E+01 -8.7136E+00 -2.8879E+00  0.0000E+00 -2.7357E+00 -5.0874E+00
            -2.1920E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1516.19085729928        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      433
 NPARAMETR:  9.9866E-01  5.2464E-01  3.0055E-01  1.2046E+00  3.5338E-01  9.4570E-01  6.8304E-01  1.0000E-02  1.0359E+00  2.1498E-01
             3.4714E+00
 PARAMETER:  9.8656E-02 -5.4504E-01 -1.1021E+00  2.8614E-01 -9.4021E-01  4.4174E-02 -2.8120E-01 -4.9390E+00  1.3532E-01 -1.4372E+00
             1.3446E+00
 GRADIENT:   2.4927E+00  9.9477E+00  5.0979E-01  9.0271E+00 -1.1376E+00  1.1781E+00 -1.4518E+00  0.0000E+00 -5.7618E+00 -2.2483E+00
            -8.2533E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1518.00597448871        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      610
 NPARAMETR:  9.9540E-01  3.9923E-01  2.3263E-01  1.1869E+00  2.7702E-01  9.3993E-01  9.4332E-01  1.0000E-02  1.1311E+00  2.3661E-01
             3.3269E+00
 PARAMETER:  9.5386E-02 -8.1822E-01 -1.3583E+00  2.7132E-01 -1.1837E+00  3.8045E-02  4.1653E-02 -6.0818E+00  2.2320E-01 -1.3414E+00
             1.3020E+00
 GRADIENT:   1.4242E+00  6.2063E+00 -3.6628E+00  2.4897E+01 -1.4429E+01 -4.1446E+00  3.2227E-01  0.0000E+00 -8.6882E+00 -6.2384E+00
            -2.3462E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1524.06945169082        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      789
 NPARAMETR:  9.9317E-01  3.7755E-01  1.9937E-01  1.1157E+00  2.4975E-01  9.5816E-01  4.9265E-01  1.0000E-02  1.2607E+00  5.2763E-01
             3.2244E+00
 PARAMETER:  9.3151E-02 -8.7406E-01 -1.5126E+00  2.0944E-01 -1.2873E+00  5.7256E-02 -6.0795E-01 -6.9174E+00  3.3163E-01 -5.3936E-01
             1.2707E+00
 GRADIENT:  -3.2314E+00 -6.4758E+00 -8.4303E+00 -1.2174E+01  1.9738E+01  1.1345E+00  1.5799E+00  0.0000E+00 -2.5808E-01 -2.8536E-01
             9.0664E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1525.03791942934        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      964
 NPARAMETR:  9.9307E-01  3.7979E-01  1.8864E-01  1.1098E+00  2.4139E-01  9.6236E-01  1.2147E-01  1.0000E-02  1.2925E+00  5.7082E-01
             3.1650E+00
 PARAMETER:  9.3050E-02 -8.6814E-01 -1.5679E+00  2.0414E-01 -1.3213E+00  6.1637E-02 -2.0081E+00 -6.0173E+00  3.5660E-01 -4.6068E-01
             1.2522E+00
 GRADIENT:  -3.1712E+00  1.6244E+00  1.6803E+00  5.4179E-02 -4.4356E+00  1.5861E+00  1.1122E-01  0.0000E+00 -9.2085E-01  2.4390E-01
            -7.4412E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1525.10158727840        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1139
 NPARAMETR:  9.9445E-01  3.8014E-01  1.8947E-01  1.1113E+00  2.4236E-01  9.5784E-01  1.7512E-02  1.0676E-02  1.2954E+00  5.6847E-01
             3.1692E+00
 PARAMETER:  9.4439E-02 -8.6723E-01 -1.5635E+00  2.0554E-01 -1.3173E+00  5.6929E-02 -3.9449E+00 -4.4397E+00  3.5880E-01 -4.6481E-01
             1.2535E+00
 GRADIENT:   6.3096E-02 -2.2481E-02 -2.9197E-01  2.0049E-01  4.5796E-01  1.4778E-02  2.2577E-03 -1.1067E-03 -5.6462E-02 -2.6470E-02
            -9.4843E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1526.76018474200        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1324
 NPARAMETR:  9.9369E-01  3.7748E-01  1.7466E-01  1.0963E+00  2.3224E-01  9.6704E-01  1.0000E-02  7.8682E-01  1.4004E+00  4.4846E-01
             3.0891E+00
 PARAMETER:  9.3675E-02 -8.7423E-01 -1.6449E+00  1.9191E-01 -1.3600E+00  6.6485E-02 -9.3847E+00 -1.3976E-01  4.3674E-01 -7.0195E-01
             1.2279E+00
 GRADIENT:   8.8694E-01  9.8770E-01  4.8448E+00  3.9700E+00 -5.0102E+00  1.8826E+00  0.0000E+00 -4.4614E-01  7.8858E+00  1.5968E+00
             1.1531E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1527.00342550977        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1486
 NPARAMETR:  9.9354E-01  3.7695E-01  1.7217E-01  1.0883E+00  2.3132E-01  9.6221E-01  1.0000E-02  8.8213E-01  1.3696E+00  4.0287E-01
             3.0808E+00
 PARAMETER:  9.3519E-02 -8.7563E-01 -1.6593E+00  1.8461E-01 -1.3639E+00  6.1479E-02 -9.5063E+00 -2.5412E-02  4.1454E-01 -8.0915E-01
             1.2252E+00
 GRADIENT:   5.1961E-02 -4.4890E-02 -9.8474E-02 -4.2309E-02  1.2831E-01  1.4092E-02  0.0000E+00 -3.2681E-03  1.1072E-02  1.1136E-02
            -5.4661E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1486
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9357E-03 -2.2135E-04  1.4484E-02 -8.0704E-03  1.3157E-02
 SE:             2.8814E-02  1.3680E-04  1.7172E-02  2.6634E-02  1.5360E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4644E-01  1.0566E-01  3.9897E-01  7.6188E-01  3.9168E-01

 ETASHRINKSD(%)  3.4685E+00  9.9542E+01  4.2472E+01  1.0772E+01  4.8542E+01
 ETASHRINKVR(%)  6.8167E+00  9.9998E+01  6.6905E+01  2.0384E+01  7.3521E+01
 EBVSHRINKSD(%)  3.2011E+00  9.9516E+01  4.2167E+01  8.8257E+00  4.9260E+01
 EBVSHRINKVR(%)  6.2998E+00  9.9998E+01  6.6553E+01  1.6873E+01  7.4254E+01
 RELATIVEINF(%)  9.3510E+01  2.9389E-04  2.7409E+00  4.6968E+01  9.0567E-01
 EPSSHRINKSD(%)  2.5145E+01
 EPSSHRINKVR(%)  4.3967E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1527.0034255097685     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -608.06489230509578     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.46
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1527.003       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  3.77E-01  1.72E-01  1.09E+00  2.31E-01  9.62E-01  1.00E-02  8.82E-01  1.37E+00  4.03E-01  3.08E+00
 


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
+       -5.52E+01  1.96E+03
 
 TH 3
+       -1.25E+02  3.50E+03  1.55E+04
 
 TH 4
+       -1.52E+01  1.50E+02 -2.33E+02  4.17E+02
 
 TH 5
+        2.63E+02 -6.92E+03 -2.03E+04 -5.87E+02  3.46E+04
 
 TH 6
+        5.10E+00 -1.67E+01  4.00E+01 -1.10E+01  2.97E+01  1.87E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -3.51E+00 -2.88E+01  7.24E+01  8.14E-01  6.79E+01  1.38E+00  0.00E+00  2.49E+01
 
 TH 9
+        1.22E+01 -4.78E+01  1.08E+02 -5.80E+00  1.36E+02 -3.60E-01  0.00E+00 -9.82E-01  6.73E+01
 
 TH10
+       -1.76E+00 -9.50E+01 -5.50E+01 -4.16E-01  4.12E+02  3.44E+00  0.00E+00  3.66E+01  9.81E+00  9.88E+01
 
 TH11
+       -1.97E+01 -1.67E+01 -5.69E+01 -4.70E+00  5.31E+01  2.85E+00  0.00E+00  1.11E+01  5.64E+00  1.30E+01  4.96E+01
 
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
 #CPUT: Total CPU Time in Seconds,       31.588
Stop Time:
Thu Sep 30 00:20:16 CDT 2021
