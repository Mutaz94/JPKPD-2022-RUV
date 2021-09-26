Sat Sep 25 02:32:46 CDT 2021
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
$DATA ../../../../data/int/SL3/dat65.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      978
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      878
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
 RAW OUTPUT FILE (FILE): m65.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -962.070932168160        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.6663E+01 -2.6818E+01  1.8868E+01  1.2433E+02  1.6102E+02  4.4867E+01 -1.3388E+02 -1.6267E+02 -1.0491E+02 -4.9100E+01
            -5.3074E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2783.16228564734        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0458E+00  1.1687E+00  1.0841E+00  9.1337E-01  1.0446E+00  8.0546E-01  1.1430E+00  9.6165E-01  1.0244E+00  1.1144E+00
             2.3256E+00
 PARAMETER:  1.4474E-01  2.5593E-01  1.8074E-01  9.3809E-03  1.4367E-01 -1.1634E-01  2.3366E-01  6.0897E-02  1.2414E-01  2.0835E-01
             9.4398E-01
 GRADIENT:   1.8923E+01 -3.0701E+00 -1.3985E+01  1.2517E+01 -5.1900E-01 -2.8390E+01  4.4906E+00  3.8867E+00 -8.5899E+00 -1.2498E+01
            -1.8637E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2791.62654355431        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0457E+00  1.5235E+00  1.6281E+00  7.2675E-01  1.3725E+00  8.2565E-01  8.1176E-01  7.0203E-01  1.1851E+00  1.5461E+00
             2.3267E+00
 PARAMETER:  1.4470E-01  5.2100E-01  5.8742E-01 -2.1917E-01  4.1661E-01 -9.1581E-02 -1.0855E-01 -2.5378E-01  2.6985E-01  5.3574E-01
             9.4446E-01
 GRADIENT:   1.7772E+01  4.3167E+01  1.4203E+01  2.3467E+01 -2.7203E+01 -1.8120E+01 -1.0468E+01 -1.3479E+00 -3.8909E+00 -4.7638E+00
            -1.8722E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2802.06245184844        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0447E+00  1.6299E+00  1.2468E+00  6.3794E-01  1.4257E+00  8.5541E-01  8.1520E-01  3.5860E-01  1.2792E+00  1.5458E+00
             2.5014E+00
 PARAMETER:  1.4372E-01  5.8851E-01  3.2056E-01 -3.4952E-01  4.5466E-01 -5.6171E-02 -1.0432E-01 -9.2555E-01  3.4621E-01  5.3551E-01
             1.0169E+00
 GRADIENT:   5.7336E+00  8.9569E+00  2.3272E+00  3.6144E+00  3.2003E+00 -3.5871E+00 -7.9248E-02 -2.6568E-02 -1.3789E+00 -2.8980E+00
            -2.3163E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2803.23855041221        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      409
 NPARAMETR:  1.0468E+00  1.9004E+00  9.2489E-01  4.6366E-01  1.5523E+00  8.6490E-01  7.3049E-01  2.4422E-01  1.6054E+00  1.6629E+00
             2.4993E+00
 PARAMETER:  1.4578E-01  7.4204E-01  2.1915E-02 -6.6860E-01  5.3972E-01 -4.5136E-02 -2.1404E-01 -1.3097E+00  5.7340E-01  6.0858E-01
             1.0160E+00
 GRADIENT:  -6.9495E-01  5.0929E+00 -7.1056E-01  3.0025E+00 -6.3792E-01 -4.9381E-01 -1.2719E+00  2.8282E-02  1.2621E-01 -5.9786E-02
            -8.4161E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2803.26983535927        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.0470E+00  1.9267E+00  9.1267E-01  4.4294E-01  1.5740E+00  8.6591E-01  7.3125E-01  2.3513E-01  1.6431E+00  1.6783E+00
             2.4994E+00
 PARAMETER:  1.4594E-01  7.5583E-01  8.6194E-03 -7.1432E-01  5.5363E-01 -4.3979E-02 -2.1300E-01 -1.3476E+00  5.9659E-01  6.1781E-01
             1.0161E+00
 GRADIENT:  -7.7191E-02 -2.1259E-01 -4.0986E-03 -1.3831E-01 -2.4555E-01 -4.3124E-02 -7.7098E-03  2.4385E-02  9.5185E-02 -4.8116E-02
             7.0171E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2803.28191285585        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      763
 NPARAMETR:  1.0470E+00  1.9244E+00  9.1479E-01  4.4458E-01  1.5734E+00  8.6601E-01  7.3181E-01  3.3582E-02  1.6385E+00  1.6783E+00
             2.4996E+00
 PARAMETER:  1.4596E-01  7.5462E-01  1.0936E-02 -7.1062E-01  5.5322E-01 -4.3854E-02 -2.1224E-01 -3.2938E+00  5.9380E-01  6.1776E-01
             1.0161E+00
 GRADIENT:  -2.3302E-02 -9.1451E-02 -1.3523E-02 -2.9272E-02  1.6105E-02  1.0864E-03  8.2574E-03  4.8739E-04 -2.0373E-03  1.1379E-02
            -1.5939E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2803.28215355933        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      944             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0470E+00  1.9236E+00  9.1592E-01  4.4513E-01  1.5730E+00  8.6601E-01  7.3192E-01  1.0000E-02  1.6374E+00  1.6779E+00
             2.4996E+00
 PARAMETER:  1.4597E-01  7.5422E-01  1.2171E-02 -7.0938E-01  5.5297E-01 -4.3863E-02 -2.1208E-01 -4.7959E+00  5.9312E-01  6.1757E-01
             1.0161E+00
 GRADIENT:   1.1541E+01  2.5961E+01  8.6099E-03  2.4323E+00  2.8657E+00  5.8305E-01  3.7964E-01  0.0000E+00  5.0712E-01  8.5578E-01
             1.7854E+00

0ITERATION NO.:   38    OBJECTIVE VALUE:  -2803.28215370687        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1037
 NPARAMETR:  1.0470E+00  1.9236E+00  9.1601E-01  4.4514E-01  1.5730E+00  8.6601E-01  7.3192E-01  1.0000E-02  1.6374E+00  1.6779E+00
             2.4996E+00
 PARAMETER:  1.4597E-01  7.5422E-01  1.2271E-02 -7.0937E-01  5.5296E-01 -4.3861E-02 -2.1209E-01 -4.7959E+00  5.9311E-01  6.1756E-01
             1.0161E+00
 GRADIENT:   5.2697E-03 -6.7885E-03  1.2449E-03 -2.7136E-03 -1.8770E-03 -7.5596E-04 -3.9224E-05  0.0000E+00  6.3681E-04  2.0663E-04
            -2.2077E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1037
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1430E-03 -2.8943E-02 -6.8280E-05  2.4028E-02 -1.9757E-02
 SE:             2.9330E-02  2.2499E-02  4.9406E-05  2.0562E-02  2.6468E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6891E-01  1.9830E-01  1.6697E-01  2.4259E-01  4.5540E-01

 ETASHRINKSD(%)  1.7414E+00  2.4625E+01  9.9834E+01  3.1113E+01  1.1330E+01
 ETASHRINKVR(%)  3.4525E+00  4.3186E+01  1.0000E+02  5.2546E+01  2.1376E+01
 EBVSHRINKSD(%)  1.8293E+00  2.3780E+01  9.9853E+01  3.4998E+01  9.2671E+00
 EBVSHRINKVR(%)  3.6251E+00  4.1905E+01  1.0000E+02  5.7748E+01  1.7675E+01
 RELATIVEINF(%)  9.6291E+01  6.9777E+00  1.0743E-04  5.0092E+00  3.4781E+01
 EPSSHRINKSD(%)  1.6727E+01
 EPSSHRINKVR(%)  3.0656E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          878
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1613.6560643074051     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2803.2821537068676     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1189.6260893994624     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2803.282       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.92E+00  9.16E-01  4.45E-01  1.57E+00  8.66E-01  7.32E-01  1.00E-02  1.64E+00  1.68E+00  2.50E+00
 


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
+        1.30E+03
 
 TH 2
+       -1.56E+01  3.37E+02
 
 TH 3
+        1.15E+00  1.76E+01  2.82E+01
 
 TH 4
+       -2.42E+01  4.13E+02 -4.95E+01  8.88E+02
 
 TH 5
+       -3.18E+00 -5.30E+01 -3.06E+01  8.03E+01  1.30E+02
 
 TH 6
+        5.48E+00 -4.10E+00  1.71E-01 -7.76E+00 -1.59E+00  2.49E+02
 
 TH 7
+        4.21E+00 -6.04E+00  3.86E-01 -1.49E+01 -3.07E+00 -1.15E-01  1.49E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.92E-01 -4.23E+00 -4.93E-01  3.67E+01  2.31E+00  2.80E-01  1.61E+01  0.00E+00  1.89E+01
 
 TH10
+        9.90E-01 -1.14E+01 -6.16E+00  2.16E+01 -8.63E+00  2.94E-01  2.23E+00  0.00E+00  3.93E+00  4.57E+01
 
 TH11
+       -1.91E+01 -1.44E+01 -2.24E+00 -1.52E+01  3.11E-01  3.38E+00  6.81E+00  0.00E+00  3.27E+00  3.21E+00  1.83E+02
 
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
 #CPUT: Total CPU Time in Seconds,       37.432
Stop Time:
Sat Sep 25 02:33:25 CDT 2021
