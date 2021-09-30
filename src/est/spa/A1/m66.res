Wed Sep 29 12:16:50 CDT 2021
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
$DATA ../../../../data/spa/A1/dat66.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1295.63412895042        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8953E+02 -1.3617E+01  1.9098E+01 -2.6945E+01  5.9841E+01  4.2910E+01 -2.7997E+01  2.7414E+00 -2.7758E+01 -6.1480E+01
            -5.4145E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1435.97426675712        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0021E+00  1.0318E+00  1.0325E+00  1.0864E+00  1.0226E+00  9.6999E-01  1.3237E+00  7.7604E-01  1.0259E+00  1.5268E+00
             2.1983E+00
 PARAMETER:  1.0208E-01  1.3134E-01  1.3195E-01  1.8284E-01  1.2238E-01  6.9535E-02  3.8047E-01 -1.5356E-01  1.2558E-01  5.2316E-01
             8.8769E-01
 GRADIENT:   1.7921E+02  2.2218E+01 -1.8982E+01  4.7412E+01 -1.1818E+01  1.3571E+01  1.2398E+01  5.3591E+00  1.1513E+01  3.9247E+01
             9.3401E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1445.20328391861        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.6930E-01  8.7795E-01  1.1309E+00  1.1681E+00  9.8991E-01  9.4260E-01  1.0524E+00  1.2833E-01  9.8994E-01  1.3026E+00
             2.1640E+00
 PARAMETER:  6.8814E-02 -3.0163E-02  2.2301E-01  2.5535E-01  8.9858E-02  4.0885E-02  1.5110E-01 -1.9532E+00  8.9893E-02  3.6436E-01
             8.7195E-01
 GRADIENT:   1.1136E+02  1.5952E+01 -3.0567E-01  5.3818E+01 -4.6079E+00  8.5820E+00 -6.0000E+00  1.0605E-01  3.6072E-01  7.2880E+00
             7.8901E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1450.27365582355        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.2421E-01  7.9112E-01  8.4024E-01  1.1743E+00  7.8706E-01  9.2276E-01  1.3995E+00  1.0512E-01  9.5656E-01  1.0581E+00
             1.8771E+00
 PARAMETER:  2.1182E-02 -1.3430E-01 -7.4072E-02  2.6071E-01 -1.3945E-01  1.9616E-02  4.3611E-01 -2.1526E+00  5.5591E-02  1.5646E-01
             7.2971E-01
 GRADIENT:   2.0720E+01  1.3035E+01  2.6721E+00  4.5238E+01 -9.8272E+00 -3.8237E-01  1.2521E+00  1.7895E-01  5.5567E+00  1.2942E+00
             2.3075E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1454.30493521816        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      385
 NPARAMETR:  9.5651E-01  6.8789E-01  1.1208E+00  1.2675E+00  9.2469E-01  9.4662E-01  1.5706E+00  2.1027E-01  9.1886E-01  1.2864E+00
             1.8309E+00
 PARAMETER:  5.5540E-02 -2.7413E-01  2.1407E-01  3.3707E-01  2.1707E-02  4.5147E-02  5.5146E-01 -1.4594E+00  1.5383E-02  3.5184E-01
             7.0481E-01
 GRADIENT:   1.4044E+01  4.2083E+00  2.8118E-01  1.3509E+00 -1.6793E+00  1.6749E+00 -4.7035E-01  2.6724E-01  3.3915E-01 -6.8840E-02
             3.3956E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1455.62066296613        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      562
 NPARAMETR:  9.4738E-01  4.0946E-01  1.1141E+00  1.4420E+00  8.2024E-01  9.4191E-01  2.2122E+00  4.1275E-02  8.5404E-01  1.2422E+00
             1.8339E+00
 PARAMETER:  4.5942E-02 -7.9292E-01  2.0807E-01  4.6604E-01 -9.8156E-02  4.0157E-02  8.9399E-01 -3.0875E+00 -5.7778E-02  3.1689E-01
             7.0642E-01
 GRADIENT:  -1.4863E+00  6.9795E+00  3.1726E+00  1.8941E+01 -8.7847E+00  6.4346E-01  1.0575E+00  1.0457E-02 -2.8898E-03  2.2373E+00
             3.2240E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1457.12654454001        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      740
 NPARAMETR:  9.3909E-01  1.2990E-01  1.1610E+00  1.6009E+00  7.7322E-01  9.3398E-01  3.9773E+00  1.0000E-02  8.1057E-01  1.2317E+00
             1.8186E+00
 PARAMETER:  3.7152E-02 -1.9410E+00  2.4929E-01  5.7054E-01 -1.5719E-01  3.1702E-02  1.4806E+00 -9.3844E+00 -1.1002E-01  3.0841E-01
             6.9809E-01
 GRADIENT:  -9.6626E+00  1.7014E+00  1.1331E+00  1.5661E+01 -2.7820E+00 -1.3388E+00 -6.0023E-01  0.0000E+00  7.9027E-01  4.0733E-01
            -3.7186E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1458.31952492601        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  9.4096E-01  3.5162E-02  1.1761E+00  1.6485E+00  7.5704E-01  9.3645E-01  6.9900E+00  1.0000E-02  7.9316E-01  1.2183E+00
             1.8474E+00
 PARAMETER:  3.9148E-02 -3.2478E+00  2.6220E-01  5.9985E-01 -1.7833E-01  3.4346E-02  2.0445E+00 -1.8006E+01 -1.3173E-01  2.9745E-01
             7.1381E-01
 GRADIENT:  -3.8762E-01 -1.0623E+00  2.3800E+00  9.0783E+00 -3.6067E+00  3.5216E-01 -3.0150E+00  0.0000E+00  2.5637E+00  1.8389E+00
             3.1786E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1458.59409058553        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1092
 NPARAMETR:  9.4081E-01  3.2046E-02  1.1317E+00  1.6410E+00  7.3749E-01  9.3560E-01  7.4939E+00  1.0000E-02  7.8977E-01  1.1915E+00
             1.8365E+00
 PARAMETER:  3.8990E-02 -3.3406E+00  2.2371E-01  5.9530E-01 -2.0450E-01  3.3428E-02  2.1141E+00 -1.8712E+01 -1.3602E-01  2.7520E-01
             7.0787E-01
 GRADIENT:  -1.6972E-01 -6.7935E-02  1.6064E-01 -1.2626E-01 -4.2988E-01 -4.0679E-02 -1.6543E-01  0.0000E+00 -1.6621E-02  2.4943E-01
            -2.7783E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1458.59409058553        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1114
 NPARAMETR:  9.4081E-01  3.2046E-02  1.1317E+00  1.6410E+00  7.3749E-01  9.3560E-01  7.4939E+00  1.0000E-02  7.8977E-01  1.1915E+00
             1.8365E+00
 PARAMETER:  3.8990E-02 -3.3406E+00  2.2371E-01  5.9530E-01 -2.0450E-01  3.3428E-02  2.1141E+00 -1.8712E+01 -1.3602E-01  2.7520E-01
             7.0787E-01
 GRADIENT:  -1.6972E-01 -6.7935E-02  1.6064E-01 -1.2626E-01 -4.2988E-01 -4.0679E-02 -1.6543E-01  0.0000E+00 -1.6621E-02  2.4943E-01
            -2.7783E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1114
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9575E-04  7.7001E-03 -6.3688E-05 -1.4879E-02 -2.9101E-02
 SE:             2.9433E-02  6.6023E-03  1.3174E-04  2.7510E-02  2.2898E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9198E-01  2.4350E-01  6.2878E-01  5.8861E-01  2.0377E-01

 ETASHRINKSD(%)  1.3960E+00  7.7882E+01  9.9559E+01  7.8384E+00  2.3288E+01
 ETASHRINKVR(%)  2.7726E+00  9.5108E+01  9.9998E+01  1.5062E+01  4.1152E+01
 EBVSHRINKSD(%)  1.4648E+00  8.4498E+01  9.9526E+01  6.8667E+00  2.0155E+01
 EBVSHRINKVR(%)  2.9081E+00  9.7597E+01  9.9998E+01  1.3262E+01  3.6247E+01
 RELATIVEINF(%)  9.6735E+01  1.0176E+00  2.0368E-04  3.4976E+01  5.7376E+00
 EPSSHRINKSD(%)  3.6614E+01
 EPSSHRINKVR(%)  5.9823E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1458.5940905855275     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -723.44326402178933     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1458.594       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.41E-01  3.20E-02  1.13E+00  1.64E+00  7.37E-01  9.36E-01  7.49E+00  1.00E-02  7.90E-01  1.19E+00  1.84E+00
 


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
+        1.39E+03
 
 TH 2
+       -8.35E+01  2.82E+05
 
 TH 3
+       -1.19E+00  1.89E+02  2.12E+02
 
 TH 4
+       -2.10E+01  8.89E+02 -5.97E+01  2.28E+03
 
 TH 5
+        1.68E+01 -7.68E+02 -4.24E+02  2.19E+01  1.04E+03
 
 TH 6
+        7.33E-01 -3.21E+00  2.53E+00 -7.42E+00 -2.76E+00  2.14E+02
 
 TH 7
+       -2.07E-01  9.30E+02  7.13E-01 -1.03E+02 -3.36E+00  1.60E-02  1.33E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.44E-01  2.62E+02  4.69E+01 -5.30E+01 -1.31E+02  1.57E+00  1.74E+00  0.00E+00  3.76E+02
 
 TH10
+       -6.90E+00  2.19E+02  2.12E+01 -3.85E+01 -3.44E+04  6.30E-01  1.10E+00  0.00E+00  1.04E+02  1.58E+04
 
 TH11
+       -1.44E+01  1.93E+01 -7.82E+00 -1.56E+01 -2.11E+00  2.59E+00  1.20E-01  0.00E+00  2.28E+01  2.19E+01  7.97E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.578
Stop Time:
Wed Sep 29 12:17:13 CDT 2021
