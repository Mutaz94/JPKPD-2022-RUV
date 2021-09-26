Sat Sep 25 12:39:34 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat5.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m5.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1657.45775574289        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5602E+02 -3.9145E+01  1.4499E-01 -5.4197E+01 -3.5993E+01  1.5182E+01 -1.2111E+01  8.0573E+00  9.0150E+00  9.7308E+00
             3.7947E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1666.35502238922        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.6360E-01  1.1668E+00  1.0880E+00  9.3654E-01  1.1776E+00  9.3189E-01  1.3111E+00  9.2311E-01  8.8679E-01  1.0117E+00
             9.0351E-01
 PARAMETER:  6.2925E-02  2.5428E-01  1.8437E-01  3.4438E-02  2.6345E-01  2.9462E-02  3.7084E-01  1.9988E-02 -2.0145E-02  1.1162E-01
            -1.4705E-03
 GRADIENT:   9.0048E+01  2.0741E+01  8.4304E+00 -3.8955E-01  2.6527E+01 -6.0769E+00  1.6834E+01 -4.1284E+00 -1.9730E+00 -1.3663E+01
            -6.5302E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1667.21682608878        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.6130E-01  1.0396E+00  1.0975E+00  1.0264E+00  1.1084E+00  9.5761E-01  1.4174E+00  8.0280E-01  8.3427E-01  1.0060E+00
             8.9239E-01
 PARAMETER:  6.0527E-02  1.3882E-01  1.9300E-01  1.2601E-01  2.0288E-01  5.6686E-02  4.4879E-01 -1.1965E-01 -8.1201E-02  1.0594E-01
            -1.3851E-02
 GRADIENT:   8.4512E+01  2.2454E+01  1.1191E+01  1.9410E+01  1.7188E+01  5.1800E+00  1.3270E+01 -5.2650E+00 -2.6726E+00 -8.6489E+00
            -1.1173E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1669.37657645335        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.3885E-01  1.0367E+00  7.9320E-01  9.9725E-01  9.3200E-01  9.4808E-01  1.3321E+00  4.8693E-01  8.4429E-01  8.5052E-01
             9.0423E-01
 PARAMETER:  3.6897E-02  1.3602E-01 -1.3168E-01  9.7249E-02  2.9583E-02  4.6686E-02  3.8675E-01 -6.1964E-01 -6.9264E-02 -6.1904E-02
            -6.7591E-04
 GRADIENT:   1.8900E+01  2.8121E+00 -8.3470E+00  1.2081E+01  1.3268E+01  8.3551E-01  4.3450E+00  1.2573E+00  4.8350E-01 -8.0794E-01
            -1.1883E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1669.37791394167        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.3634E-01  1.0365E+00  7.6427E-01  9.9268E-01  9.1149E-01  9.4786E-01  1.3201E+00  4.3474E-01  8.4424E-01  8.3153E-01
             9.0531E-01
 PARAMETER:  3.4227E-02  1.3585E-01 -1.6883E-01  9.2657E-02  7.3256E-03  4.6457E-02  3.7767E-01 -7.3300E-01 -6.9316E-02 -8.4482E-02
             5.2411E-04
 GRADIENT:   1.1454E+01  1.3071E+00 -6.7553E+00  8.3108E+00  9.2506E+00  5.1077E-01  2.8957E+00  1.1867E+00  5.2209E-01 -1.7930E-01
            -4.8791E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1669.83013059463        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      480
 NPARAMETR:  9.4750E-01  1.0060E+00  7.3813E-01  1.0090E+00  8.7559E-01  9.5468E-01  1.3594E+00  2.7497E-01  8.3203E-01  8.0564E-01
             9.0650E-01
 PARAMETER:  4.6068E-02  1.0596E-01 -2.0363E-01  1.0898E-01 -3.2862E-02  5.3618E-02  4.0705E-01 -1.1911E+00 -8.3883E-02 -1.1612E-01
             1.8353E-03
 GRADIENT:  -5.0333E+00  5.3663E-01  1.3806E+00  2.6353E-01  5.3991E-01 -1.4127E-02  7.6575E-01  1.1866E-01  4.5977E-02 -1.0087E+00
            -6.0091E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1669.96426524385        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      655
 NPARAMETR:  9.4962E-01  1.0159E+00  6.8568E-01  9.9666E-01  8.4543E-01  9.5502E-01  1.3426E+00  7.4127E-02  8.3182E-01  7.7727E-01
             9.0671E-01
 PARAMETER:  4.8312E-02  1.1579E-01 -2.7734E-01  9.6653E-02 -6.7904E-02  5.3973E-02  3.9459E-01 -2.5020E+00 -8.4137E-02 -1.5197E-01
             2.0712E-03
 GRADIENT:  -9.4307E-01 -5.3267E-01 -1.2191E+00 -2.2071E-01 -1.6842E-01 -4.2850E-02  4.2467E-01  2.7720E-02  2.6516E-01  6.3900E-01
             1.0360E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1669.98137113388        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  9.4995E-01  1.0445E+00  6.8672E-01  9.8072E-01  8.6091E-01  9.5503E-01  1.3100E+00  2.1476E-02  8.4359E-01  7.8770E-01
             9.0705E-01
 PARAMETER:  4.8653E-02  1.4352E-01 -2.7583E-01  8.0531E-02 -4.9768E-02  5.3987E-02  3.7003E-01 -3.7408E+00 -7.0092E-02 -1.3864E-01
             2.4466E-03
 GRADIENT:  -2.6349E-01 -2.1103E-01 -3.9291E-01  1.6730E-01  4.0489E-01 -5.9391E-02 -4.4233E-03  1.7042E-03  1.7211E-01  1.4835E-01
             3.8943E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1669.98256135707        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      993
 NPARAMETR:  9.5004E-01  1.0385E+00  6.8662E-01  9.8414E-01  8.5753E-01  9.5516E-01  1.3168E+00  1.0000E-02  8.4038E-01  7.8444E-01
             9.0706E-01
 PARAMETER:  4.8748E-02  1.3776E-01 -2.7597E-01  8.4015E-02 -5.3694E-02  5.4125E-02  3.7520E-01 -4.5675E+00 -7.3906E-02 -1.4279E-01
             2.4579E-03
 GRADIENT:   2.5821E-04  7.5080E-03  2.4800E-02 -1.5070E-02 -2.2108E-02  1.3024E-03 -6.5521E-05  0.0000E+00 -3.8225E-03 -3.5613E-03
            -1.3859E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      993
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.5168E-05  6.6166E-04 -5.1462E-04 -3.6963E-03 -1.0604E-02
 SE:             2.9847E-02  2.3811E-02  2.0294E-04  2.3684E-02  2.1784E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9879E-01  9.7783E-01  1.1218E-02  8.7598E-01  6.2643E-01

 ETASHRINKSD(%)  1.0252E-02  2.0232E+01  9.9320E+01  2.0656E+01  2.7020E+01
 ETASHRINKVR(%)  2.0503E-02  3.6370E+01  9.9995E+01  3.7045E+01  4.6739E+01
 EBVSHRINKSD(%)  3.9006E-01  1.9446E+01  9.9401E+01  2.1296E+01  2.6475E+01
 EBVSHRINKVR(%)  7.7860E-01  3.5110E+01  9.9996E+01  3.8057E+01  4.5940E+01
 RELATIVEINF(%)  9.9019E+01  5.8818E+00  4.3548E-04  5.7317E+00  5.0555E+00
 EPSSHRINKSD(%)  4.4275E+01
 EPSSHRINKVR(%)  6.8947E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1669.9825613570658     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -934.83173479332766     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1669.983       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.50E-01  1.04E+00  6.87E-01  9.84E-01  8.58E-01  9.55E-01  1.32E+00  1.00E-02  8.40E-01  7.84E-01  9.07E-01
 


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
+        1.34E+03
 
 TH 2
+       -5.13E+00  3.98E+02
 
 TH 3
+        1.96E+01  2.43E+02  9.61E+02
 
 TH 4
+       -1.09E+01  2.95E+02 -4.79E+02  1.02E+03
 
 TH 5
+       -5.70E+00 -3.52E+02 -9.60E+02  5.13E+02  1.32E+03
 
 TH 6
+       -3.14E-01 -9.45E-01  3.44E+00 -3.00E+00 -3.00E+00  2.16E+02
 
 TH 7
+        1.11E+00  2.97E+01 -2.54E+01 -1.50E+01  7.89E+00 -9.23E-02  5.09E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.51E+00 -2.09E+01 -4.44E+01  3.61E+01  5.84E+00 -4.93E-01  1.66E+01  0.00E+00  1.25E+02
 
 TH10
+       -1.77E+00 -1.41E+01 -7.88E+01 -2.40E+01 -6.68E+01  3.95E-01  1.19E+01  0.00E+00  1.69E+01  1.09E+02
 
 TH11
+       -8.56E+00 -1.28E+01 -4.31E+01 -5.72E-01  6.88E+00  2.43E+00  5.70E+00  0.00E+00  1.16E+01  2.25E+01  2.61E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.085
Stop Time:
Sat Sep 25 12:39:51 CDT 2021
