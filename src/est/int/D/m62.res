Sat Sep 18 07:15:33 CDT 2021
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
$DATA ../../../../data/int/D/dat62.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   29174.6601461978        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6234E+01  4.0857E+02 -6.4052E+01  1.5209E+02  3.3372E+02 -2.2545E+03 -1.0731E+03 -9.6607E+01 -2.1083E+03 -7.1890E+02
            -5.9776E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1013.45944772694        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  2.0452E+00  1.7077E+00  9.3229E-01  2.2833E+00  7.6956E-01  6.4879E+00  3.9412E+00  1.0225E+00  3.6251E+00  1.7967E+00
             1.1588E+01
 PARAMETER:  8.1547E-01  6.3515E-01  2.9885E-02  9.2560E-01 -1.6194E-01  1.9699E+00  1.4715E+00  1.2228E-01  1.3879E+00  6.8595E-01
             2.5500E+00
 GRADIENT:   1.7086E+01  2.3557E+01 -2.8385E+01  7.8644E+01 -5.0144E+01  1.2682E+02 -8.9916E+00  4.1581E+00  4.1137E+01  1.8703E+01
             3.7121E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1032.44549366570        NO. OF FUNC. EVALS.:  84
 CUMULATIVE NO. OF FUNC. EVALS.:      167
 NPARAMETR:  2.1363E+00  3.6235E+00  1.0139E+00  1.5833E+00  1.4945E+00  6.3946E+00  3.5505E+00  1.5403E+00  4.1521E+00  2.5475E+00
             1.1554E+01
 PARAMETER:  8.5906E-01  1.3874E+00  1.1378E-01  5.5949E-01  5.0177E-01  1.9555E+00  1.3671E+00  5.3200E-01  1.5236E+00  1.0351E+00
             2.5470E+00
 GRADIENT:   1.9581E+01  5.5612E+01 -2.2576E+01  8.4843E+01 -1.7380E+01  1.2202E+02 -1.7648E+01  2.4198E+00  2.2382E+00  5.3116E+01
             3.6509E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1240.55923900641        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      246
 NPARAMETR:  1.0090E+00  1.0070E+00  5.9317E+02  1.6819E+00  3.3527E+00  2.1179E+00  5.6513E+00  5.0703E+01  2.8399E+00  1.1940E+00
             9.8733E+00
 PARAMETER:  1.0897E-01  1.0701E-01  6.4855E+00  6.1991E-01  1.3098E+00  8.5041E-01  1.8319E+00  4.0260E+00  1.1438E+00  2.7727E-01
             2.3898E+00
 GRADIENT:  -8.8121E+01 -1.4199E+01 -4.5323E-01  1.0066E+01  7.7772E+01 -4.2841E+01 -7.1846E+01  1.0811E+01  2.8530E+01  1.6812E+01
             2.2865E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1342.53358783162        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      319
 NPARAMETR:  1.2928E+00  4.7036E-01  5.2941E+02  1.5226E+00  2.2138E+00  2.7421E+00  8.5508E+00  4.5278E+00  1.6835E+00  7.2268E-01
             8.5091E+00
 PARAMETER:  3.5684E-01 -6.5426E-01  6.3718E+00  5.2041E-01  8.9472E-01  1.1087E+00  2.2460E+00  1.6102E+00  6.2087E-01 -2.2479E-01
             2.2411E+00
 GRADIENT:   3.8479E+00 -5.9264E+00  5.5549E-03  4.0553E+00  3.0353E+00  5.1358E+00  6.8234E+00  1.7116E-03  5.4140E+00  6.1364E+00
             4.3218E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1344.52380347299        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      393
 NPARAMETR:  1.2443E+00  6.7477E-01  2.6459E+02  1.3320E+00  2.1927E+00  2.6946E+00  7.4011E+00  7.3682E+00  1.5515E+00  3.6806E-01
             8.2565E+00
 PARAMETER:  3.1860E-01 -2.9338E-01  5.6782E+00  3.8666E-01  8.8515E-01  1.0913E+00  2.1016E+00  2.0972E+00  5.3923E-01 -8.9950E-01
             2.2110E+00
 GRADIENT:  -3.9746E+00 -7.0050E+00  9.5869E-02 -5.0767E+00 -3.0610E+00 -2.0690E+00 -5.2134E+00 -1.4980E-02  4.6088E+00  8.9214E-01
            -1.4802E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1344.88488505838        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      497
 NPARAMETR:  1.2439E+00  7.6941E-01  2.0986E+02  1.2967E+00  2.1993E+00  2.6872E+00  7.2178E+00  8.4318E+00  1.5190E+00  3.1771E-01
             8.2569E+00
 PARAMETER:  3.1823E-01 -1.6214E-01  5.4465E+00  3.5981E-01  8.8814E-01  1.0885E+00  2.0765E+00  2.2320E+00  5.1807E-01 -1.0466E+00
             2.2110E+00
 GRADIENT:  -6.8339E+00 -5.5161E+00  1.4286E-01 -6.2056E+00 -2.1598E+00 -7.7822E+00 -2.0166E+01 -5.6692E-02  4.0127E+00  6.3273E-01
            -2.0858E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1346.27389267415        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      666
 NPARAMETR:  1.2437E+00  8.1653E-01  2.1047E+02  1.2969E+00  2.2004E+00  2.6857E+00  7.8801E+00  8.4219E+00  1.3799E+00  2.3256E-01
             8.2705E+00
 PARAMETER:  3.1806E-01 -1.0269E-01  5.4493E+00  3.6000E-01  8.8862E-01  1.0879E+00  2.1643E+00  2.2308E+00  4.2202E-01 -1.3586E+00
             2.2127E+00
 GRADIENT:  -6.7478E+00  3.7347E-02  1.8415E-01 -4.3448E+00 -4.8221E+00 -7.1501E+00 -1.6979E-01 -7.8520E-02  7.2306E-02  2.7215E-01
            -2.7401E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1346.29452205366        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      843
 NPARAMETR:  1.2437E+00  8.1604E-01  2.1047E+02  1.2969E+00  2.2004E+00  2.6857E+00  7.8780E+00  8.4218E+00  1.3774E+00  2.1836E-01
             8.2719E+00
 PARAMETER:  3.1806E-01 -1.0329E-01  5.4493E+00  3.6000E-01  8.8863E-01  1.0879E+00  2.1641E+00  2.2308E+00  4.2017E-01 -1.4216E+00
             2.2129E+00
 GRADIENT:   2.4627E+03 -7.6048E+03 -1.4339E+02 -1.0955E+03  8.8301E+02  7.1500E+02  3.6346E+02  3.5072E+02 -9.3467E+02 -5.5217E+02
            -2.0313E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      843
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7465E-02  4.3585E-02 -1.0647E-03 -8.1779E-02  8.9433E-05
 SE:             2.9972E-02  2.3565E-02  8.6936E-04  1.5086E-02  3.4331E-03
 N:                     100         100         100         100         100

 P VAL.:         5.6009E-01  6.4374E-02  2.2067E-01  5.9408E-08  9.7922E-01

 ETASHRINKSD(%)  1.0000E-10  2.1055E+01  9.7088E+01  4.9461E+01  8.8499E+01
 ETASHRINKVR(%)  1.0000E-10  3.7677E+01  9.9915E+01  7.4458E+01  9.8677E+01
 EBVSHRINKSD(%)  2.3869E+00  1.4775E+01  9.7086E+01  5.2220E+01  8.8063E+01
 EBVSHRINKVR(%)  4.7168E+00  2.7367E+01  9.9915E+01  7.7171E+01  9.8575E+01
 RELATIVEINF(%)  9.5148E+01  4.1301E+01  2.1735E-02  1.2930E+01  3.6033E-01
 EPSSHRINKSD(%)  6.6490E+00
 EPSSHRINKVR(%)  1.2856E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1346.2945220536590     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       307.79483771475179     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.86
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1346.295       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.24E+00  8.16E-01  2.10E+02  1.30E+00  2.20E+00  2.69E+00  7.88E+00  8.42E+00  1.38E+00  2.18E-01  8.27E+00
 


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
+        1.26E+06
 
 TH 2
+       -1.15E+02  2.76E+07
 
 TH 3
+       -4.31E+02  2.84E-02  1.48E-01
 
 TH 4
+       -9.61E+02  1.23E+02  3.65E+02  9.02E+05
 
 TH 5
+        2.38E+03 -2.64E+01 -8.30E-01 -2.03E+03  5.19E+04
 
 TH 6
+        6.84E+01 -1.13E+01 -2.40E-02 -5.85E+01  1.48E+01  2.30E+04
 
 TH 7
+        1.36E+02 -6.33E+02 -4.68E-02 -1.21E+02  2.77E+01  1.84E+01  6.80E+02
 
 TH 8
+        2.63E+04 -2.88E+00  9.18E+00 -2.23E+04  5.07E+01  1.50E+00  2.88E+00  5.52E+02
 
 TH 9
+        1.48E+02 -7.26E+02 -5.36E-02 -1.56E+02  3.40E+01  2.03E+01  4.87E+00  3.07E+00  5.86E+05
 
 TH10
+        7.16E+02 -3.36E+03 -2.45E-01 -6.14E+02  1.48E+02  9.65E+01  1.67E+01  1.50E+01  1.09E+06  2.04E+06
 
 TH11
+       -2.70E+04 -8.93E+02 -6.60E-02  2.29E+04 -5.49E+03  2.77E+01  4.44E+00  4.08E+00  2.63E-01 -1.47E+01  5.96E+02
 
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
 #CPUT: Total CPU Time in Seconds,       43.958
Stop Time:
Sat Sep 18 07:16:19 CDT 2021
