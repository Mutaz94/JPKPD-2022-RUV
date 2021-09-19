Sat Sep 18 14:30:07 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat15.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1649.40331969551        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.6154E+01 -1.0013E+02 -1.9293E+01 -1.1518E+02  7.4839E+01  2.4136E+01 -1.4710E+01 -9.2569E+00 -2.8035E+01 -6.5921E+00
            -4.2090E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1663.99861961902        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  1.0443E+00  1.0812E+00  9.6191E-01  1.0309E+00  9.4719E-01  9.1141E-01  1.0783E+00  1.1209E+00  1.0987E+00  9.2797E-01
             1.0803E+00
 PARAMETER:  1.4332E-01  1.7811E-01  6.1171E-02  1.3043E-01  4.5741E-02  7.2424E-03  1.7538E-01  2.1413E-01  1.9411E-01  2.5245E-02
             1.7728E-01
 GRADIENT:   4.3892E+01 -2.4404E+00 -5.2846E-01  3.6312E+00  1.6084E+00 -7.4685E+00  5.6420E-01 -3.2451E+00  6.0122E+00  1.3048E+00
            -5.3367E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1664.55995666599        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0418E+00  1.1141E+00  1.0341E+00  1.0155E+00  9.8397E-01  9.3354E-01  1.0385E+00  1.3624E+00  1.0935E+00  8.7670E-01
             1.1093E+00
 PARAMETER:  1.4096E-01  2.0809E-01  1.3353E-01  1.1540E-01  8.3838E-02  3.1233E-02  1.3777E-01  4.0924E-01  1.8936E-01 -3.1591E-02
             2.0372E-01
 GRADIENT:   3.5252E+01  1.6135E+00  8.1764E-01  1.4592E+00  2.7399E+00  2.3763E+00 -2.3108E-02 -1.1596E+00  4.5414E-01 -4.9570E+00
             4.0914E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1665.13421536852        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0505E+00  1.2431E+00  1.1199E+00  9.4277E-01  1.0772E+00  9.3370E-01  8.9050E-01  1.6201E+00  1.1900E+00  1.0152E+00
             1.0912E+00
 PARAMETER:  1.4924E-01  3.1757E-01  2.1326E-01  4.1068E-02  1.7434E-01  3.1404E-02 -1.5967E-02  5.8246E-01  2.7399E-01  1.1507E-01
             1.8725E-01
 GRADIENT:   7.5661E+00  2.5941E-01  4.7429E-01  1.2267E+00 -1.9443E+00 -1.8341E-01 -9.3716E-02 -9.8371E-02 -8.8964E-01  1.7316E+00
            -2.0499E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1665.51702416071        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      474
 NPARAMETR:  1.0477E+00  1.5186E+00  8.9752E-01  7.5968E-01  1.1224E+00  9.3486E-01  8.3946E-01  1.6565E+00  1.3646E+00  1.0159E+00
             1.0987E+00
 PARAMETER:  1.4656E-01  5.1781E-01 -8.1238E-03 -1.7486E-01  2.1543E-01  3.2647E-02 -7.4992E-02  6.0473E-01  4.1087E-01  1.1573E-01
             1.9415E-01
 GRADIENT:  -2.6212E+00  3.8677E+00  1.4235E+00  1.9885E+00 -1.8586E+00 -1.9870E-01  7.6903E-01 -2.9921E-01  5.3566E-03  7.4525E-01
             2.0573E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1665.70205994193        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      650
 NPARAMETR:  1.0491E+00  1.7109E+00  6.5724E-01  6.2448E-01  1.1264E+00  9.3660E-01  8.1093E-01  1.5165E+00  1.5132E+00  9.8918E-01
             1.1000E+00
 PARAMETER:  1.4792E-01  6.3700E-01 -3.1970E-01 -3.7083E-01  2.1904E-01  3.4497E-02 -1.0957E-01  5.1639E-01  5.1423E-01  8.9124E-02
             1.9531E-01
 GRADIENT:  -1.3420E+00  7.9260E-01 -5.1421E-01  1.4795E+00  7.7225E-01  3.5620E-02 -6.5348E-02  7.5581E-02 -1.3719E-01 -6.2433E-02
             3.9122E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1665.71867877270        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      825
 NPARAMETR:  1.0493E+00  1.7760E+00  6.0576E-01  5.7881E-01  1.1418E+00  9.3664E-01  7.9711E-01  1.5092E+00  1.5849E+00  9.9827E-01
             1.1001E+00
 PARAMETER:  1.4813E-01  6.7435E-01 -4.0127E-01 -4.4679E-01  2.3260E-01  3.4545E-02 -1.2677E-01  5.1156E-01  5.6053E-01  9.8267E-02
             1.9538E-01
 GRADIENT:  -9.4320E-01  1.8447E-01 -5.0882E-02  2.1334E-01  2.2147E-01 -5.1083E-03 -6.1402E-02 -5.7870E-02 -1.5100E-02 -1.5444E-02
             1.6424E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1665.71905527777        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1005
 NPARAMETR:  1.0498E+00  1.7770E+00  6.0669E-01  5.7785E-01  1.1430E+00  9.3665E-01  7.9679E-01  1.5167E+00  1.5870E+00  9.9918E-01
             1.1000E+00
 PARAMETER:  1.4864E-01  6.7495E-01 -3.9973E-01 -4.4843E-01  2.3364E-01  3.4550E-02 -1.2716E-01  5.1654E-01  5.6183E-01  9.9181E-02
             1.9529E-01
 GRADIENT:   3.4839E-01 -2.0302E-01 -4.6307E-03 -5.4991E-02  1.5134E-01  3.2867E-03  1.6610E-02 -3.6566E-02 -3.1616E-03 -1.3555E-02
             1.1622E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1665.71913884731        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1138
 NPARAMETR:  1.0497E+00  1.7771E+00  6.0678E-01  5.7789E-01  1.1429E+00  9.3664E-01  7.9661E-01  1.5193E+00  1.5868E+00  9.9932E-01
             1.0997E+00
 PARAMETER:  1.4851E-01  6.7500E-01 -3.9973E-01 -4.4840E-01  2.3353E-01  3.4546E-02 -1.2727E-01  5.1841E-01  5.6156E-01  9.9290E-02
             1.9502E-01
 GRADIENT:   3.2583E+05 -6.1195E-03 -1.2104E+05 -2.5148E-02 -2.0721E+05  5.8100E-05  1.3841E-02  9.3272E+04 -4.2079E-02 -4.8389E+05
            -2.4814E+05
 NUMSIGDIG:         3.3         5.4         3.3         4.0         3.3         5.8         2.8         3.3         3.4         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1138
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.3553E-04 -3.5602E-02 -2.9399E-02  2.9040E-02 -4.8086E-02
 SE:             2.9817E-02  2.2525E-02  1.1736E-02  2.2915E-02  2.0672E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8299E-01  1.1398E-01  1.2244E-02  2.0505E-01  2.0011E-02

 ETASHRINKSD(%)  1.1051E-01  2.4538E+01  6.0683E+01  2.3233E+01  3.0746E+01
 ETASHRINKVR(%)  2.2089E-01  4.3055E+01  8.4542E+01  4.1069E+01  5.2040E+01
 EBVSHRINKSD(%)  6.0746E-01  2.4195E+01  6.2990E+01  2.4216E+01  2.8589E+01
 EBVSHRINKVR(%)  1.2112E+00  4.2536E+01  8.6303E+01  4.2568E+01  4.9005E+01
 RELATIVEINF(%)  9.8741E+01  3.2940E+00  1.1929E+00  3.4466E+00  1.3929E+01
 EPSSHRINKSD(%)  4.5264E+01
 EPSSHRINKVR(%)  7.0040E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1665.7191388473093     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -930.56831228357112     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.40
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1665.719       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.78E+00  6.07E-01  5.78E-01  1.14E+00  9.37E-01  7.97E-01  1.52E+00  1.59E+00  9.99E-01  1.10E+00
 


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
+        4.98E+08
 
 TH 2
+       -7.21E+00  3.41E+02
 
 TH 3
+       -2.18E+03  1.05E+02  2.06E+08
 
 TH 4
+       -1.58E+01  2.94E+02 -1.69E+02  7.63E+02
 
 TH 5
+       -1.99E+03 -1.51E+02 -5.65E+04  1.59E+02  1.70E+08
 
 TH 6
+        4.49E+03 -1.06E+00 -2.89E+03 -4.11E+00 -2.63E+03  2.18E+02
 
 TH 7
+        5.56E+01 -1.61E-01 -3.55E+01 -9.66E+00 -4.20E+01 -4.59E+00  1.13E+02
 
 TH 8
+        6.72E+02 -9.56E+00  1.90E+04  8.36E+00  9.93E+02  8.89E+02  1.66E+01  1.95E+07
 
 TH 9
+        1.33E+00 -1.50E+01 -1.38E+01 -5.24E+07 -5.97E+00 -4.43E-01  1.42E+01 -1.72E+07  3.49E+01
 
 TH10
+       -5.29E+03 -7.58E+00 -1.50E+05 -5.00E+00  4.54E+08 -7.01E+03 -7.03E+01  4.57E+02  2.12E+00  1.21E+09
 
 TH11
+       -2.48E+03 -1.25E+01  1.79E+04 -1.42E+00  1.63E+04 -3.26E+03 -3.18E+01 -5.52E+03  3.22E+00  4.36E+04  2.63E+08
 
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
 #CPUT: Total CPU Time in Seconds,       20.807
Stop Time:
Sat Sep 18 14:30:29 CDT 2021
