Sat Sep 18 09:33:57 CDT 2021
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
$DATA ../../../../data/spa/A1/dat98.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m98.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1366.65691151039        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2042E+02 -1.0182E+01 -1.0512E+01 -2.6934E+00  7.6629E+01  2.1225E+01 -2.5654E+01  1.6566E+00 -2.8443E+01 -2.8750E+01
            -4.6916E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1486.84769664123        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.6176E-01  9.9899E-01  1.1501E+00  1.0429E+00  9.8738E-01  8.9583E-01  1.0833E+00  8.8616E-01  1.0826E+00  9.1748E-01
             1.8306E+00
 PARAMETER:  6.1012E-02  9.8993E-02  2.3987E-01  1.4204E-01  8.7298E-02 -1.0004E-02  1.7999E-01 -2.0855E-02  1.7937E-01  1.3871E-02
             7.0466E-01
 GRADIENT:  -9.9219E+00  2.7230E+01  8.0301E+00  2.3681E+01 -1.7794E+01 -1.4649E+01  1.3692E+00  3.8169E+00  4.3772E+00 -5.2487E+00
            -1.3813E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1489.91214178090        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.5397E-01  8.0514E-01  1.3034E+00  1.1592E+00  9.9457E-01  9.1874E-01  1.1797E+00  3.6614E-01  1.0051E+00  1.1111E+00
             1.8173E+00
 PARAMETER:  5.2880E-02 -1.1674E-01  3.6495E-01  2.4776E-01  9.4555E-02  1.5253E-02  2.6528E-01 -9.0473E-01  1.0510E-01  2.0537E-01
             6.9737E-01
 GRADIENT:  -2.4054E+01  1.3358E+01  6.0189E+00  1.4857E+01 -8.7633E+00 -3.9081E+00  1.2191E+00  5.2758E-01  1.5472E+00  7.9940E+00
             3.5966E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1493.86237725746        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.6072E-01  5.0819E-01  8.8000E-01  1.3048E+00  6.9771E-01  9.2454E-01  1.8051E+00  8.0166E-02  8.7065E-01  8.6774E-01
             1.7582E+00
 PARAMETER:  5.9928E-02 -5.7690E-01 -2.7839E-02  3.6608E-01 -2.5995E-01  2.1539E-02  6.9062E-01 -2.4237E+00 -3.8513E-02 -4.1868E-02
             6.6426E-01
 GRADIENT:  -5.0312E+00  1.4521E+01  8.2620E+00  1.8886E+01 -1.9546E+01 -2.5349E+00  2.1455E+00  7.6368E-02 -1.8599E+00  2.2396E+00
             7.6307E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1496.96980563292        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.5940E-01  2.4236E-01  7.8664E-01  1.4184E+00  5.8526E-01  9.2956E-01  2.7161E+00  1.0000E-02  8.1810E-01  8.3320E-01
             1.7326E+00
 PARAMETER:  5.8557E-02 -1.3173E+00 -1.3998E-01  4.4951E-01 -4.3569E-01  2.6960E-02  1.0992E+00 -5.6071E+00 -1.0077E-01 -8.2486E-02
             6.4960E-01
 GRADIENT:   2.5214E+00  6.3472E+00  7.0721E+00 -4.8265E+00 -1.2127E+01  5.4239E-01  4.4338E+00  0.0000E+00 -5.4983E+00 -6.6023E-01
            -1.7719E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1499.80915032232        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.5249E-01  6.2622E-02  7.6139E-01  1.5093E+00  5.3752E-01  9.2527E-01  4.8140E+00  1.0000E-02  8.1364E-01  8.1665E-01
             1.7206E+00
 PARAMETER:  5.1327E-02 -2.6706E+00 -1.7261E-01  5.1165E-01 -5.2080E-01  2.2327E-02  1.6715E+00 -1.2207E+01 -1.0623E-01 -1.0255E-01
             6.4268E-01
 GRADIENT:  -4.0542E+00 -1.0248E+00  1.7104E+01  2.3232E+01 -2.2792E+01 -3.0624E-01 -6.5419E+00  0.0000E+00  6.2520E+00 -1.7712E+00
            -3.0507E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1501.37072721775        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      463
 NPARAMETR:  9.5362E-01  2.1824E-02  7.1322E-01  1.5152E+00  5.1010E-01  9.2695E-01  8.1985E+00  1.0000E-02  7.8073E-01  8.3126E-01
             1.7048E+00
 PARAMETER:  5.2513E-02 -3.7247E+00 -2.3797E-01  5.1555E-01 -5.7314E-01  2.4145E-02  2.2040E+00 -1.7833E+01 -1.4752E-01 -8.4809E-02
             6.3345E-01
 GRADIENT:  -1.1059E+01 -8.3014E-01  4.0364E+00 -7.3308E+00 -1.1900E+01 -7.5706E-01 -3.8775E+00  0.0000E+00 -5.3954E+00  2.7588E+00
            -1.1190E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1502.42726508895        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      644
 NPARAMETR:  9.5796E-01  1.0000E-02  8.0432E-01  1.5478E+00  5.5473E-01  9.2760E-01  1.2505E+01  1.0000E-02  7.9034E-01  8.5529E-01
             1.7233E+00
 PARAMETER:  5.7055E-02 -4.5164E+00 -1.1776E-01  5.3683E-01 -4.8927E-01  2.4841E-02  2.6261E+00 -2.1769E+01 -1.3529E-01 -5.6311E-02
             6.4424E-01
 GRADIENT:   2.3127E+00  0.0000E+00  4.1240E+00 -2.6443E-01 -8.1204E+00  3.5306E-02  5.7848E-01  0.0000E+00  3.8746E-01  9.2567E-01
             4.8822E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1502.48008980200        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      822
 NPARAMETR:  9.5693E-01  1.0000E-02  8.4270E-01  1.5561E+00  5.7581E-01  9.2710E-01  1.2525E+01  1.0000E-02  7.8744E-01  8.6327E-01
             1.7270E+00
 PARAMETER:  5.5973E-02 -4.5096E+00 -7.1149E-02  5.4217E-01 -4.5199E-01  2.4307E-02  2.6277E+00 -2.1662E+01 -1.3897E-01 -4.7023E-02
             6.4638E-01
 GRADIENT:  -4.8631E-02  6.3892E-03 -1.5458E-01 -2.5677E-01  2.4896E-01 -8.1409E-03  4.2026E-01  0.0000E+00 -1.5119E-01 -1.8264E-01
            -1.1061E-01

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1502.48008980200        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      874
 NPARAMETR:  9.5693E-01  1.0000E-02  8.4273E-01  1.5562E+00  5.7579E-01  9.2711E-01  1.2521E+01  1.0000E-02  7.8750E-01  8.6336E-01
             1.7271E+00
 PARAMETER:  5.5973E-02 -4.5096E+00 -7.1149E-02  5.4217E-01 -4.5199E-01  2.4307E-02  2.6277E+00 -2.1662E+01 -1.3897E-01 -4.7023E-02
             6.4638E-01
 GRADIENT:  -5.0534E-02  3.5114E-03 -1.2214E-01 -1.7148E+03  5.3427E-01 -7.0956E-03  3.5315E+02  0.0000E+00 -1.2154E-01 -1.5826E-01
            -1.3061E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      874
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.0809E-06  5.9662E-03 -6.0393E-05 -9.4025E-03 -1.7362E-02
 SE:             2.9537E-02  4.7361E-03  1.8122E-04  2.8097E-02  2.3490E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9994E-01  2.0777E-01  7.3894E-01  7.3789E-01  4.5985E-01

 ETASHRINKSD(%)  1.0460E+00  8.4134E+01  9.9393E+01  5.8729E+00  2.1304E+01
 ETASHRINKVR(%)  2.0810E+00  9.7483E+01  9.9996E+01  1.1401E+01  3.8069E+01
 EBVSHRINKSD(%)  1.2956E+00  8.7831E+01  9.9338E+01  5.5874E+00  2.0075E+01
 EBVSHRINKVR(%)  2.5744E+00  9.8519E+01  9.9996E+01  1.0863E+01  3.6120E+01
 RELATIVEINF(%)  9.7174E+01  1.0148E+00  2.4514E-04  4.5271E+01  3.6032E+00
 EPSSHRINKSD(%)  3.7744E+01
 EPSSHRINKVR(%)  6.1242E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1502.4800898020044     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -767.32926323826621     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.41
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1502.480       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.57E-01  1.00E-02  8.43E-01  1.56E+00  5.76E-01  9.27E-01  1.25E+01  1.00E-02  7.87E-01  8.63E-01  1.73E+00
 


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
+        1.40E+03
 
 TH 2
+        8.97E+01  1.15E+08
 
 TH 3
+       -4.18E+01 -5.65E+02  1.01E+03
 
 TH 4
+        1.05E+02  3.99E+02 -1.04E+03  3.27E+05
 
 TH 5
+        4.63E+02  1.99E+07 -1.06E+07 -1.06E+06  3.43E+06
 
 TH 6
+       -3.45E+00 -8.76E+00 -2.49E+00 -1.24E+01 -2.44E+01  2.16E+02
 
 TH 7
+       -3.07E+00 -6.94E+00  2.25E+01  3.35E+01 -6.01E+01  1.92E-01  2.14E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.30E+01 -6.92E+02  3.82E+02 -3.77E+02 -8.17E+06  1.28E+00  7.17E+00  0.00E+00  6.96E+02
 
 TH10
+       -7.01E+01 -7.36E+02  3.76E+02 -7.37E+02 -1.04E+07 -1.21E+00  1.66E+01  0.00E+00  4.41E+02  5.61E+02
 
 TH11
+       -3.50E+01 -2.51E+02  1.14E+02 -2.35E+02 -8.01E+05  2.15E+00  5.11E+00  0.00E+00  1.68E+02  1.86E+02  1.87E+05
 
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
 #CPUT: Total CPU Time in Seconds,       15.848
Stop Time:
Sat Sep 18 09:34:15 CDT 2021
