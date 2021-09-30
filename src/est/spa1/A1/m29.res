Wed Sep 29 21:59:15 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1414.10904393281        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6264E+02  1.7354E+01 -1.7036E+00  4.3603E+01  1.2696E+02  2.4551E+01 -3.7729E+01  2.5175E-01 -1.6269E+01 -6.0020E+01
            -1.2343E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1736.63279835512        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1445E+00  1.0031E+00  1.1805E+00  1.0798E+00  9.2218E-01  1.1685E+00  1.1861E+00  8.0521E-01  1.1519E+00  9.8234E-01
             2.3895E+00
 PARAMETER:  2.3497E-01  1.0312E-01  2.6596E-01  1.7675E-01  1.8982E-02  2.5570E-01  2.7064E-01 -1.1666E-01  2.4140E-01  8.2186E-02
             9.7108E-01
 GRADIENT:   3.9962E+02  4.3422E+01  2.5748E+01  3.9611E+01 -7.0858E+01  3.2238E+01  1.0479E+01  5.1196E+00  2.5437E+01  3.8865E+00
             1.2639E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1759.38229685189        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0438E+00  9.3343E-01  9.4015E-01  1.0747E+00  8.8470E-01  1.0314E+00  6.1001E-01  4.1181E-02  1.0468E+00  1.0635E+00
             2.2627E+00
 PARAMETER:  1.4282E-01  3.1113E-02  3.8284E-02  1.7202E-01 -2.2508E-02  1.3091E-01 -3.9428E-01 -3.0898E+00  1.4570E-01  1.6155E-01
             9.1654E-01
 GRADIENT:   2.1289E+02  1.4724E+01 -2.9166E+01  2.3414E+01  2.4472E+01  1.0505E+01 -8.6395E+00  2.4822E-02 -2.5217E+01  1.1639E+01
             9.2772E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1770.38543040111        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.4681E-01  9.0490E-01  1.0702E+00  1.0737E+00  9.1753E-01  9.8829E-01  1.1444E+00  1.3371E-01  1.0446E+00  1.1118E+00
             1.9777E+00
 PARAMETER:  4.5348E-02  6.8894E-05  1.6781E-01  1.7107E-01  1.3927E-02  8.8220E-02  2.3485E-01 -1.9120E+00  1.4361E-01  2.0596E-01
             7.8193E-01
 GRADIENT:   7.6906E+00  4.8965E+00  2.7974E+00  5.5531E-02 -7.2068E+00 -5.1649E+00 -7.5763E-03  2.1434E-01  3.3597E+00  9.5495E+00
            -1.9566E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1770.68052602038        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.3903E-01  8.0553E-01  7.7504E-01  1.1048E+00  7.3643E-01  1.0125E+00  1.3243E+00  9.5505E-02  9.6066E-01  8.9102E-01
             1.9510E+00
 PARAMETER:  3.7096E-02 -1.1625E-01 -1.5485E-01  1.9968E-01 -2.0594E-01  1.1246E-01  3.8091E-01 -2.2486E+00  5.9870E-02 -1.5388E-02
             7.6836E-01
 GRADIENT:  -7.0467E+00  1.1703E+00 -2.3569E+00  7.8632E+00  1.2611E+01  3.7987E+00  1.5715E-01  1.7735E-01  6.9344E-01 -6.6458E-01
            -4.1502E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1777.80665266196        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      462
 NPARAMETR:  1.0042E+00  6.1907E-01  9.1247E-01  1.2853E+00  7.3213E-01  1.0281E+00  1.5271E+00  1.2279E-01  9.0849E-01  9.4542E-01
             1.9874E+00
 PARAMETER:  1.0419E-01 -3.7954E-01  8.3953E-03  3.5096E-01 -2.1180E-01  1.2770E-01  5.2339E-01 -1.9973E+00  4.0235E-03  4.3875E-02
             7.8681E-01
 GRADIENT:   2.9057E+01  2.4053E+01 -1.8587E+00  4.1776E+01 -7.7272E+00  1.9839E+00 -2.2947E+00  2.7039E-01 -6.5696E+00  1.4250E+00
             5.7301E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1784.99013759387        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      639
 NPARAMETR:  9.7064E-01  2.2499E-01  1.0010E+00  1.5147E+00  6.6496E-01  1.0093E+00  2.6352E+00  1.3781E-02  8.6497E-01  9.4536E-01
             1.9914E+00
 PARAMETER:  7.0199E-02 -1.3917E+00  1.0104E-01  5.1522E-01 -3.0802E-01  1.0924E-01  1.0690E+00 -4.1844E+00 -4.5057E-02  4.3815E-02
             7.8885E-01
 GRADIENT:  -2.5785E+01  1.0040E+01  1.3273E+01  4.8943E+01 -2.3987E+01 -2.6458E+00  5.8659E-01  3.3125E-03  1.4643E+00 -2.1310E+00
             1.0138E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1788.37473448955        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      814
 NPARAMETR:  9.8125E-01  6.6273E-02  9.4288E-01  1.5586E+00  6.1656E-01  1.0012E+00  3.9258E+00  1.0000E-02  8.4849E-01  9.5613E-01
             1.9310E+00
 PARAMETER:  8.1076E-02 -2.6140E+00  4.1178E-02  5.4376E-01 -3.8360E-01  1.0122E-01  1.4676E+00 -7.9323E+00 -6.4294E-02  5.5136E-02
             7.5806E-01
 GRADIENT:   8.6463E+00  8.7401E-01  1.2039E+00 -3.1899E+00 -1.6462E+00 -4.2540E+00 -4.8567E-01  0.0000E+00  3.2756E+00  2.3646E-01
            -1.0585E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1789.32340146907        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      989
 NPARAMETR:  9.7573E-01  1.0000E-02  9.5903E-01  1.5928E+00  6.1388E-01  1.0150E+00  8.1499E+00  1.0000E-02  8.3011E-01  9.5454E-01
             1.9539E+00
 PARAMETER:  7.5426E-02 -4.7119E+00  5.8167E-02  5.6547E-01 -3.8796E-01  1.1493E-01  2.1980E+00 -1.4717E+01 -8.6201E-02  5.3470E-02
             7.6984E-01
 GRADIENT:  -8.0676E-01  0.0000E+00 -6.3010E-01 -1.3747E+00  1.3225E+00  1.3816E+00 -6.2238E-02  0.0000E+00  2.2662E-01 -4.0650E-01
             7.1702E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1789.33722605800        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1167
 NPARAMETR:  9.7687E-01  1.0000E-02  9.5897E-01  1.5948E+00  6.1328E-01  1.0062E+00  1.1362E+01  1.0000E-02  8.2862E-01  9.5837E-01
             1.9494E+00
 PARAMETER:  7.6596E-02 -5.6098E+00  5.8104E-02  5.6674E-01 -3.8894E-01  1.0614E-01  2.5303E+00 -1.7603E+01 -8.7990E-02  5.7474E-02
             7.6751E-01
 GRADIENT:   1.6633E+00  0.0000E+00 -8.5340E-03  2.7670E+00 -5.2637E-01 -1.9370E+00 -3.6729E-02  0.0000E+00 -2.4254E-01  2.2540E-01
            -1.3469E+00

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1789.35121790203        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1260
 NPARAMETR:  9.7611E-01  1.0000E-02  9.5738E-01  1.5928E+00  6.1275E-01  1.0114E+00  1.1427E+01  1.0000E-02  8.2943E-01  9.5598E-01
             1.9531E+00
 PARAMETER:  7.5817E-02 -5.6320E+00  5.6442E-02  5.6551E-01 -3.8981E-01  1.1129E-01  2.5360E+00 -1.7691E+01 -8.7021E-02  5.4982E-02
             7.6940E-01
 GRADIENT:   5.2956E-03  0.0000E+00 -9.1188E-02 -5.6863E-01  1.5113E-01  3.1103E-02 -2.5043E-02  0.0000E+00  1.1156E-01  1.1251E-02
             4.6828E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1260
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.7160E-04 -6.5586E-05 -4.0483E-05 -8.2750E-03 -1.8583E-02
 SE:             2.9530E-02  1.8802E-03  1.4839E-04  2.8636E-02  2.3237E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7375E-01  9.7217E-01  7.8499E-01  7.7260E-01  4.2389E-01

 ETASHRINKSD(%)  1.0707E+00  9.3701E+01  9.9503E+01  4.0669E+00  2.2151E+01
 ETASHRINKVR(%)  2.1299E+00  9.9603E+01  9.9998E+01  7.9684E+00  3.9396E+01
 EBVSHRINKSD(%)  1.1832E+00  9.4362E+01  9.9417E+01  3.9800E+00  2.1292E+01
 EBVSHRINKVR(%)  2.3524E+00  9.9682E+01  9.9997E+01  7.8015E+00  3.8050E+01
 RELATIVEINF(%)  9.1598E+01  1.4026E-02  2.6633E-04  5.2718E+00  4.9384E+00
 EPSSHRINKSD(%)  2.8487E+01
 EPSSHRINKVR(%)  4.8859E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1789.3512179020297     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -870.41268469735701     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.54
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1789.351       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  1.00E-02  9.57E-01  1.59E+00  6.13E-01  1.01E+00  1.14E+01  1.00E-02  8.29E-01  9.56E-01  1.95E+00
 


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
+        1.11E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.07E+01  0.00E+00  4.09E+02
 
 TH 4
+       -1.41E+01  0.00E+00 -8.15E+01  5.88E+02
 
 TH 5
+        2.19E+01  0.00E+00 -8.45E+02 -6.73E+01  2.01E+03
 
 TH 6
+        2.94E+00  0.00E+00  3.40E-01 -4.97E+00 -2.01E+00  1.83E+02
 
 TH 7
+        7.79E-03  0.00E+00 -8.43E-02 -4.24E+01  2.28E-01 -3.10E-03  1.35E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.61E+00  0.00E+00  1.45E+01 -6.71E+00 -1.78E+00 -4.81E-01 -1.04E-01  0.00E+00  2.48E+02
 
 TH10
+       -6.00E-01  0.00E+00 -6.65E+00  3.94E+00 -6.27E+01  1.27E+00 -1.28E-01  0.00E+00  2.86E-01  9.28E+01
 
 TH11
+       -1.27E+01  0.00E+00 -9.91E+00 -1.14E+01  4.88E+00  2.69E+00 -1.90E-02  0.00E+00  6.93E+00  1.78E+01  1.19E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       25.036
Stop Time:
Wed Sep 29 21:59:41 CDT 2021
