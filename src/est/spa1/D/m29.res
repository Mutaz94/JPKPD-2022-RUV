Thu Sep 30 02:53:37 CDT 2021
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
$DATA ../../../../data/spa1/D/dat29.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12369.7774988889        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1766E+02  7.6234E+01 -1.0824E+02  8.9749E+00  2.5899E+02 -1.1875E+03 -6.1947E+02 -1.9163E+01 -1.0570E+03 -3.6103E+02
            -2.5451E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -725.399609240840        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.5355E+00  1.1634E+00  1.1493E+00  2.2278E+00  1.1024E+00  3.3039E+00  1.6166E+00  9.3581E-01  3.4383E+00  1.0734E+00
             1.2169E+01
 PARAMETER:  5.2885E-01  2.5138E-01  2.3915E-01  9.0100E-01  1.9751E-01  1.2951E+00  5.8035E-01  3.3655E-02  1.3350E+00  1.7080E-01
             2.5989E+00
 GRADIENT:   4.0373E+01 -4.1605E-01 -3.5815E+01  5.8068E+01  1.2419E+00  1.3663E+02  7.8338E+00  5.4644E+00  7.0947E+01  6.1592E+00
             2.6575E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -788.783408346096        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.4007E+00  1.2670E+00  6.9360E+00  2.1384E+00  8.0593E+00  2.1012E+00  2.7093E+00  3.4767E-01  2.8509E+00  1.1615E+01
             1.1272E+01
 PARAMETER:  4.3701E-01  3.3662E-01  2.0367E+00  8.6005E-01  2.1868E+00  8.4249E-01  1.0967E+00 -9.5651E-01  1.1476E+00  2.5523E+00
             2.5223E+00
 GRADIENT:   3.4381E+01  1.5255E+01  2.2863E+00  5.8711E+01 -4.7039E+00  1.6083E+01  1.3609E+01 -1.1831E-02  3.9858E+01  2.0980E+01
             2.3388E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -848.767566195941        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0981E+00  1.5453E+00  1.0817E+01  1.0426E+00  4.5743E+00  1.8906E+00  1.0178E+00  4.5573E+00  2.9872E+00  5.5459E+00
             8.6431E+00
 PARAMETER:  1.9361E-01  5.3520E-01  2.4811E+00  1.4174E-01  1.6205E+00  7.3690E-01  1.1766E-01  1.6167E+00  1.1943E+00  1.8131E+00
             2.2568E+00
 GRADIENT:  -5.5786E+01  4.5704E+00  4.5412E+00 -1.5119E+01 -5.2931E+00 -1.7461E+01  2.3629E+00  3.9994E+00 -1.3774E+01  3.3387E+00
             1.3004E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -873.984907685224        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.1897E+00  8.5481E-01  1.2446E+00  1.5743E+00  9.6194E+00  1.8306E+00  1.1302E+00  1.4448E+00  1.9909E+00  4.0057E+00
             7.6480E+00
 PARAMETER:  2.7373E-01 -5.6874E-02  3.1881E-01  5.5383E-01  2.3638E+00  7.0462E-01  2.2237E-01  4.6800E-01  7.8860E-01  1.4877E+00
             2.1344E+00
 GRADIENT:   4.4540E+01  4.7464E+01 -3.8562E+00  5.6202E+01 -2.5270E+00 -1.1989E+01  2.7782E-01  2.8787E+00 -2.6479E+01 -4.4145E-01
            -2.6236E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -931.476166207111        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  5.7859E-01  4.9920E-02  1.0998E-01  8.0839E-01  1.3662E+03  1.8248E+00  1.0000E-02  1.2056E-01  1.2338E+00  5.6525E+00
             6.9595E+00
 PARAMETER: -4.4715E-01 -2.8973E+00 -2.1075E+00 -1.1271E-01  7.3198E+00  7.0145E-01 -6.7471E+00 -2.0156E+00  3.1006E-01  1.8321E+00
             2.0401E+00
 GRADIENT:  -8.5003E+01  3.8023E+00  6.3147E+01  3.7979E+01 -3.8278E-03  3.6982E+01  0.0000E+00 -1.1466E-01  1.8034E+01  2.4030E-06
            -4.4939E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -960.348791549664        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  6.1214E-01  2.7762E-02  5.9416E-02  5.4164E-01  4.7427E+03  1.6055E+00  1.0000E-02  6.1393E-02  1.1340E+00  6.6137E+00
             7.6056E+00
 PARAMETER: -3.9080E-01 -3.4841E+00 -2.7232E+00 -5.1316E-01  8.5644E+00  5.7343E-01 -8.7219E+00 -2.6905E+00  2.2572E-01  1.9891E+00
             2.1289E+00
 GRADIENT:   1.1223E+01  5.8961E-01  1.7913E+00 -5.3094E+00  2.1435E-04  2.4117E+01  0.0000E+00 -2.8558E-02  1.5074E+01  1.0988E-07
             6.4444E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -963.475144005048        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      710
 NPARAMETR:  5.3635E-01  1.7058E-02  4.2479E-02  4.2973E-01  1.3779E+04  1.4404E+00  1.0000E-02  4.1610E-02  1.0513E+00  7.7338E+00
             7.5777E+00
 PARAMETER: -5.2297E-01 -3.9711E+00 -3.0588E+00 -7.4459E-01  9.6309E+00  4.6493E-01 -1.0135E+01 -3.0794E+00  1.5005E-01  2.1456E+00
             2.1252E+00
 GRADIENT:  -1.5235E+00  4.9763E-02  1.9536E+00 -1.7593E+00  3.7165E-05 -1.5951E+00  0.0000E+00 -2.3070E-02  5.0383E-01 -3.9210E-10
             1.6960E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -963.495668730485        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      867
 NPARAMETR:  5.3331E-01  1.6446E-02  4.1940E-02  4.2816E-01  1.3778E+04  1.4428E+00  1.0000E-02  4.3378E-02  1.0495E+00  1.2540E-02
             7.5536E+00
 PARAMETER: -5.2866E-01 -4.0077E+00 -3.0715E+00 -7.4825E-01  9.6308E+00  4.6660E-01 -1.0135E+01 -3.0378E+00  1.4834E-01 -4.2788E+00
             2.1220E+00
 GRADIENT:  -1.9473E+00  4.3130E-02 -4.1242E+00  6.6348E+00  5.1269E-05 -8.2098E-01  0.0000E+00 -2.5118E-02 -3.4118E-01  0.0000E+00
            -1.8366E+00

0ITERATION NO.:   42    OBJECTIVE VALUE:  -963.495674022529        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      924
 NPARAMETR:  5.3613E-01  1.6572E-02  4.2282E-02  4.2929E-01  1.3777E+04  1.4457E+00  1.0000E-02  4.3880E-02  1.0493E+00  1.2540E-02
             7.5643E+00
 PARAMETER: -5.2337E-01 -4.0001E+00 -3.0634E+00 -7.4563E-01  9.6308E+00  4.6861E-01 -1.0135E+01 -3.0263E+00  1.4813E-01 -4.2788E+00
             2.1234E+00
 GRADIENT:   2.6041E-01  1.6334E-02 -1.3271E+00  1.7389E+00  5.2979E-05 -2.3812E-01  0.0000E+00 -2.5758E-02 -1.1819E-01  0.0000E+00
            -5.1326E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      924
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4994E-03 -8.6216E-08  1.0247E-04 -1.8980E-02 -8.6740E-10
 SE:             2.9046E-02  1.5683E-06  8.9950E-04  2.6043E-02  1.1786E-09
 N:                     100         100         100         100         100

 P VAL.:         9.3143E-01  9.5616E-01  9.0930E-01  4.6612E-01  4.6176E-01

 ETASHRINKSD(%)  2.6920E+00  9.9995E+01  9.6987E+01  1.2753E+01  1.0000E+02
 ETASHRINKVR(%)  5.3115E+00  1.0000E+02  9.9909E+01  2.3879E+01  1.0000E+02
 EBVSHRINKSD(%)  2.5790E+00  9.9993E+01  9.7198E+01  1.2558E+01  1.0000E+02
 EBVSHRINKVR(%)  5.0915E+00  1.0000E+02  9.9922E+01  2.3538E+01  1.0000E+02
 RELATIVEINF(%)  4.1706E+00  5.5104E-08  1.0092E-03  9.7762E-01  0.0000E+00
 EPSSHRINKSD(%)  1.3656E+01
 EPSSHRINKVR(%)  2.5448E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -963.49567402252910     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -44.557140817856407     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -963.496       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.36E-01  1.66E-02  4.23E-02  4.29E-01  1.38E+04  1.45E+00  1.00E-02  4.39E-02  1.05E+00  1.25E-02  7.56E+00
 


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
+        1.74E+03
 
 TH 2
+       -4.90E+02  2.08E+03
 
 TH 3
+       -7.36E+03 -3.88E+01  2.47E+05
 
 TH 4
+       -4.43E+01  1.96E+02 -2.90E+04  4.04E+03
 
 TH 5
+        2.75E-07 -2.78E-07 -2.42E-06  2.02E-07 -1.61E-13
 
 TH 6
+        6.56E+00  1.77E+01 -2.11E+02 -1.29E+01  8.44E-09  8.48E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        4.32E-01 -7.62E-01 -1.16E+01  1.14E-01 -1.38E-08  5.52E-01  0.00E+00 -1.35E+01
 
 TH 9
+        6.88E+00 -2.15E+01  4.58E+02 -7.30E+01 -2.36E-08 -6.19E-01  0.00E+00  2.70E+00  1.06E+02
 
 TH10
+       -9.32E-03  7.14E-03  1.75E-03  3.78E-04 -1.55E-08 -5.45E-03  0.00E+00  2.34E-02  3.33E-03  7.03E-02
 
 TH11
+       -2.24E+01  5.38E+00  1.91E+02 -1.50E+01 -4.38E-09  1.79E+00  0.00E+00  1.05E-01  3.04E+00  1.82E-04  9.04E+00
 
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
 #CPUT: Total CPU Time in Seconds,       24.153
Stop Time:
Thu Sep 30 02:54:03 CDT 2021
