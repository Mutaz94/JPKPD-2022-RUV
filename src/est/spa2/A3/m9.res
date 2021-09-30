Thu Sep 30 06:11:05 CDT 2021
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
$DATA ../../../../data/spa2/A3/dat9.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m9.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -780.865566363762        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3272E+02  7.6473E+01  1.8079E+02 -4.2318E+00  2.1589E+02  2.0161E+01 -5.3904E+01 -1.6891E+02 -2.9554E+01 -1.2612E+02
            -2.8230E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1770.94694935464        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.3418E-01  9.4106E-01  8.5272E-01  1.1368E+00  8.2312E-01  8.8087E-01  9.6856E-01  9.9607E-01  9.6135E-01  1.1659E+00
             3.4673E+00
 PARAMETER:  3.1917E-02  3.9248E-02 -5.9320E-02  2.2821E-01 -9.4654E-02 -2.6843E-02  6.8056E-02  9.6059E-02  6.0584E-02  2.5352E-01
             1.3434E+00
 GRADIENT:  -7.1748E+00  1.9777E+01 -2.8608E+00  3.5056E+01 -1.2136E+01 -1.2740E+01  1.2398E+01  1.0195E+01  1.8117E+01  2.6605E+01
             2.7777E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1781.03411541714        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.2253E-01  5.7990E-01  4.6055E-01  1.3758E+00  4.3473E-01  9.8961E-01  6.5407E-01  6.1354E-01  9.6080E-01  1.0915E+00
             3.1841E+00
 PARAMETER:  1.9364E-02 -4.4490E-01 -6.7532E-01  4.1906E-01 -7.3303E-01  8.9557E-02 -3.2454E-01 -3.8850E-01  6.0006E-02  1.8759E-01
             1.2582E+00
 GRADIENT:  -2.9271E+01  1.0646E+02  2.4877E+01  2.4835E+02 -6.3772E+01  1.6670E+01  1.4335E+00  8.1302E+00 -1.1479E+01  4.1386E+01
             2.2794E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1819.81336415843        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      271
 NPARAMETR:  9.1813E-01  5.5543E-01  4.8692E-01  1.2974E+00  4.7115E-01  9.6001E-01  3.8030E-01  7.3679E-02  8.5660E-01  1.0028E+00
             2.8789E+00
 PARAMETER:  1.4588E-02 -4.8802E-01 -6.1965E-01  3.6037E-01 -6.5258E-01  5.9190E-02 -8.6681E-01 -2.5080E+00 -5.4788E-02  1.0279E-01
             1.1574E+00
 GRADIENT:  -5.8814E+01  3.5370E+01 -3.9651E+00  1.1676E+02 -1.4118E+01  4.5352E+00 -8.0403E-01  9.1875E-02 -3.0883E+01  1.3805E+01
             1.2435E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1835.17609130860        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  9.3346E-01  5.0219E-01  4.2023E-01  1.2211E+00  4.1943E-01  9.4601E-01  1.6539E-01  1.7515E-02  9.8117E-01  9.3807E-01
             2.5377E+00
 PARAMETER:  3.1146E-02 -5.8878E-01 -7.6695E-01  2.9972E-01 -7.6886E-01  4.4501E-02 -1.6995E+00 -3.9447E+00  8.0993E-02  3.6074E-02
             1.0313E+00
 GRADIENT:  -2.2405E+00  1.6983E+00  3.9170E+00 -6.8707E+00 -2.3274E+00 -1.9716E-02 -3.8606E-01  3.9612E-03 -3.5205E+00 -2.8787E+00
             1.6811E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1838.03398305924        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      626
 NPARAMETR:  9.3581E-01  4.8289E-01  3.9247E-01  1.2499E+00  3.9150E-01  9.3259E-01  9.3960E-01  1.5545E-02  1.0091E+00  8.7591E-01
             2.5721E+00
 PARAMETER:  3.3655E-02 -6.2796E-01 -8.3529E-01  3.2307E-01 -8.3777E-01  3.0214E-02  3.7697E-02 -4.0640E+00  1.0902E-01 -3.2486E-02
             1.0447E+00
 GRADIENT:   1.1100E+00  1.9335E+01  1.5857E+01  3.5232E+01 -2.9483E+01 -5.9827E+00 -1.1764E+00  3.7848E-03  6.6100E+00  1.3512E+01
             3.9833E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1846.79237518163        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      806
 NPARAMETR:  9.3350E-01  3.3556E-01  2.3449E-01  1.1660E+00  2.6652E-01  9.5479E-01  1.5414E+00  1.0000E-02  1.0890E+00  6.0099E-01
             2.3892E+00
 PARAMETER:  3.1184E-02 -9.9197E-01 -1.3504E+00  2.5358E-01 -1.2223E+00  5.3740E-02  5.3270E-01 -6.7301E+00  1.8525E-01 -4.0918E-01
             9.7095E-01
 GRADIENT:   1.6909E+00  9.0999E+00 -1.1275E+01  1.2203E+01  6.5993E+00 -8.5733E-01  5.7608E+00  0.0000E+00 -1.8257E+00 -2.0336E+00
            -2.8390E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1847.15336823361        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      982
 NPARAMETR:  9.3294E-01  3.2385E-01  2.3073E-01  1.1542E+00  2.6113E-01  9.5653E-01  1.4853E+00  1.0000E-02  1.0966E+00  6.1691E-01
             2.3876E+00
 PARAMETER:  3.0590E-02 -1.0275E+00 -1.3665E+00  2.4340E-01 -1.2428E+00  5.5553E-02  4.9561E-01 -7.2771E+00  1.9220E-01 -3.8303E-01
             9.7028E-01
 GRADIENT:  -4.4192E-02 -3.7840E-02 -1.2631E-01  3.4966E-01 -2.3379E-01  3.2707E-02 -1.1171E-01  0.0000E+00 -1.9705E-01 -2.5909E-01
             1.5215E-01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1847.15336823361        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1004
 NPARAMETR:  9.3294E-01  3.2385E-01  2.3073E-01  1.1542E+00  2.6113E-01  9.5653E-01  1.4853E+00  1.0000E-02  1.0966E+00  6.1691E-01
             2.3876E+00
 PARAMETER:  3.0590E-02 -1.0275E+00 -1.3665E+00  2.4340E-01 -1.2428E+00  5.5553E-02  4.9561E-01 -7.2771E+00  1.9220E-01 -3.8303E-01
             9.7028E-01
 GRADIENT:  -4.4192E-02 -3.7840E-02 -1.2631E-01  3.4966E-01 -2.3379E-01  3.2707E-02 -1.1171E-01  0.0000E+00 -1.9705E-01 -2.5909E-01
             1.5215E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1004
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.9688E-04  1.7701E-02 -1.3394E-04 -9.7575E-03  1.1171E-02
 SE:             2.9258E-02  2.1202E-02  2.5848E-04  2.7449E-02  2.2094E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7555E-01  4.0378E-01  6.0434E-01  7.2223E-01  6.1314E-01

 ETASHRINKSD(%)  1.9818E+00  2.8971E+01  9.9134E+01  8.0418E+00  2.5982E+01
 ETASHRINKVR(%)  3.9244E+00  4.9548E+01  9.9993E+01  1.5437E+01  4.5213E+01
 EBVSHRINKSD(%)  1.9431E+00  2.7658E+01  9.9196E+01  7.2277E+00  2.6309E+01
 EBVSHRINKVR(%)  3.8484E+00  4.7667E+01  9.9994E+01  1.3933E+01  4.5697E+01
 RELATIVEINF(%)  9.6071E+01  1.0555E+01  4.0141E-04  5.0545E+01  2.3409E+00
 EPSSHRINKSD(%)  2.6445E+01
 EPSSHRINKVR(%)  4.5896E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1847.1533682336103     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -744.42712838800321     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.83
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1847.153       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.33E-01  3.24E-01  2.31E-01  1.15E+00  2.61E-01  9.57E-01  1.49E+00  1.00E-02  1.10E+00  6.17E-01  2.39E+00
 


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
+        9.52E+00  2.56E+03
 
 TH 3
+       -3.51E+01  1.46E+03  1.24E+04
 
 TH 4
+       -1.57E+01  4.79E+01 -5.33E+02  5.99E+02
 
 TH 5
+        7.47E+01 -4.63E+03 -1.51E+04 -1.09E+02  2.37E+04
 
 TH 6
+        4.92E+00 -5.60E+00  3.54E+01 -9.07E+00 -1.63E+01  1.99E+02
 
 TH 7
+        6.91E-02  7.08E+01 -2.70E+01  1.28E-01 -1.24E+01  1.46E-01  2.77E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.04E+01 -1.90E+01  1.64E+02 -3.86E+00 -1.19E+01 -2.03E-01 -1.92E+00  0.00E+00  1.14E+02
 
 TH10
+       -6.69E+00  5.81E+00 -1.96E+02  4.37E+00  1.76E+02  8.33E-01  1.59E+01  0.00E+00  3.57E+00  1.69E+02
 
 TH11
+       -1.81E+01 -1.45E+01 -1.03E+02 -9.39E+00  9.49E+01  3.13E+00  7.22E+00  0.00E+00  1.10E+01  1.79E+01  9.68E+01
 
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
 #CPUT: Total CPU Time in Seconds,       26.116
Stop Time:
Thu Sep 30 06:11:32 CDT 2021
