Wed Sep 29 21:09:57 CDT 2021
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
$DATA ../../../../data/spa1/B/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1491.03187275261        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9679E+02 -2.3845E+01  4.4039E+01 -1.0302E+00  6.1430E+01  2.0097E+01 -3.3608E+01 -2.2512E+02 -6.2910E+01 -1.1222E+00
            -8.1857E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1942.12287053627        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      140
 NPARAMETR:  8.1286E-01  1.0146E+00  9.7444E-01  9.8225E-01  9.8052E-01  1.0178E+00  1.0699E+00  1.1699E+00  1.0489E+00  9.9501E-01
             1.6475E+00
 PARAMETER: -1.0720E-01  1.1454E-01  7.4105E-02  8.2092E-02  8.0330E-02  1.1760E-01  1.6756E-01  2.5693E-01  1.4776E-01  9.5002E-02
             5.9928E-01
 GRADIENT:  -3.4521E+02 -5.6363E+01 -3.4358E+01 -3.5188E+01  8.5620E+00 -5.4978E+01 -2.8882E+01  9.9013E+00 -6.4096E+00  1.8482E+01
             2.7124E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1971.19358703608        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      318
 NPARAMETR:  8.6596E-01  7.7253E-01  1.0094E+00  1.1947E+00  8.9998E-01  9.5838E-01  2.5497E+00  8.8390E-01  5.9812E-01  9.6975E-01
             1.5271E+00
 PARAMETER: -4.3914E-02 -1.5809E-01  1.0935E-01  2.7792E-01 -5.3778E-03  5.7491E-02  1.0360E+00 -2.3407E-02 -4.1397E-01  6.9288E-02
             5.2339E-01
 GRADIENT:  -2.2945E+02  3.3013E+01 -3.9684E+01  9.1068E+01  2.0929E+01 -4.5720E+01  2.4731E+01  2.4575E+00 -3.6648E+01  2.0992E+01
             2.3515E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2043.68258792167        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      507             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1313E-01  7.2269E-01  9.8565E-01  1.2141E+00  8.5682E-01  9.5182E-01  2.1338E+00  6.9051E-01  9.0202E-01  8.1364E-01
             1.1268E+00
 PARAMETER:  9.1248E-03 -2.2477E-01  8.5545E-02  2.9401E-01 -5.4531E-02  5.0621E-02  8.5790E-01 -2.7033E-01 -3.1153E-03 -1.0624E-01
             2.1940E-01
 GRADIENT:   2.3157E+02  5.1184E+01 -1.2152E+01  3.6036E+02  3.1215E+01  2.1176E+01  6.2478E+01 -7.2066E+00  1.8000E+01  2.9082E-02
             8.0374E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2047.57806412581        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      664
 NPARAMETR:  9.5115E-01  6.7275E-01  9.9790E-01  1.2142E+00  8.3998E-01  1.0172E+00  2.0417E+00  6.9054E-01  9.0204E-01  8.1892E-01
             1.1268E+00
 PARAMETER:  4.9916E-02 -2.9639E-01  9.7894E-02  2.9406E-01 -7.4371E-02  1.1705E-01  8.1378E-01 -2.7028E-01 -3.0990E-03 -9.9775E-02
             2.1937E-01
 GRADIENT:   3.7283E+00 -9.2711E-01  2.3576E+00  2.5373E+01  8.9724E-02 -1.4561E-01 -9.0362E-01 -8.5289E+00  7.3506E+00  1.6831E-01
             7.8124E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2047.61426596854        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      839
 NPARAMETR:  9.4962E-01  6.7204E-01  9.7863E-01  1.2142E+00  8.2993E-01  1.0177E+00  2.0632E+00  6.9054E-01  9.0204E-01  8.0147E-01
             1.1268E+00
 PARAMETER:  4.8305E-02 -2.9743E-01  7.8397E-02  2.9405E-01 -8.6410E-02  1.1758E-01  8.2425E-01 -2.7028E-01 -3.1000E-03 -1.2131E-01
             2.1936E-01
 GRADIENT:   1.5553E-01  1.5359E-02 -3.6662E-02  2.8486E+01 -4.0036E-02  2.0271E-02  3.8933E-03 -7.5078E+00  7.8472E+00  7.9848E-03
             7.8246E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2049.81439885969        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1012
 NPARAMETR:  9.4953E-01  6.5911E-01  9.7782E-01  1.1831E+00  8.2958E-01  9.9902E-01  1.9763E+00  6.9223E-01  8.8454E-01  8.0132E-01
             1.0863E+00
 PARAMETER:  4.8210E-02 -3.1687E-01  7.7571E-02  2.6812E-01 -8.6837E-02  9.9015E-02  7.8123E-01 -2.6784E-01 -2.2693E-02 -1.2149E-01
             1.8274E-01
 GRADIENT:   5.5153E-01 -1.9572E+01  9.2229E+00 -2.5860E+01 -6.1821E+00 -7.5175E+00 -8.3830E+00 -8.0747E+00  2.4732E+00 -8.2162E-01
             5.3764E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2050.81868118271        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1182
 NPARAMETR:  9.4969E-01  7.1179E-01  9.5823E-01  1.1823E+00  8.3444E-01  1.0271E+00  2.0197E+00  6.9176E-01  8.8431E-01  7.9123E-01
             1.0868E+00
 PARAMETER:  4.8383E-02 -2.3997E-01  5.7333E-02  2.6744E-01 -8.0996E-02  1.2670E-01  8.0296E-01 -2.6851E-01 -2.2944E-02 -1.3416E-01
             1.8320E-01
 GRADIENT:  -1.1560E-01  2.3911E-01  2.3342E+00  1.8270E+01 -1.2925E+00  3.3206E+00  1.3741E+00 -6.9671E+00  2.4262E+00 -6.2877E-01
             5.3872E+01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2050.81868118271        NO. OF FUNC. EVALS.:  26
 CUMULATIVE NO. OF FUNC. EVALS.:     1208
 NPARAMETR:  9.4969E-01  7.1179E-01  9.5823E-01  1.1823E+00  8.3444E-01  1.0271E+00  2.0197E+00  6.9176E-01  8.8431E-01  7.9123E-01
             1.0868E+00
 PARAMETER:  4.8383E-02 -2.3997E-01  5.7333E-02  2.6744E-01 -8.0996E-02  1.2670E-01  8.0296E-01 -2.6851E-01 -2.2944E-02 -1.3416E-01
             1.8320E-01
 GRADIENT:   2.6992E+00  1.2211E+05 -2.9297E+05  5.4784E+04 -8.5267E-01 -2.3127E+05  1.8249E+04  1.0896E+05  1.4649E+05 -3.2776E+00
            -1.6026E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1208
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6272E-03  1.8469E-02 -2.9084E-02 -2.7456E-02 -9.0866E-03
 SE:             2.9645E-02  2.1984E-02  1.3430E-02  2.3477E-02  2.0331E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5623E-01  4.0084E-01  3.0335E-02  2.4222E-01  6.5493E-01

 ETASHRINKSD(%)  6.8469E-01  2.6351E+01  5.5009E+01  2.1348E+01  3.1888E+01
 ETASHRINKVR(%)  1.3647E+00  4.5758E+01  7.9758E+01  3.8138E+01  5.3607E+01
 EBVSHRINKSD(%)  4.0003E-01  2.6805E+01  6.1253E+01  1.9725E+01  3.0441E+01
 EBVSHRINKVR(%)  7.9846E-01  4.6425E+01  8.4987E+01  3.5559E+01  5.1615E+01
 RELATIVEINF(%)  9.8780E+01  8.3971E+00  2.4477E+00  1.0752E+01  8.0609E+00
 EPSSHRINKSD(%)  3.7554E+01
 EPSSHRINKVR(%)  6.1005E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2050.8186811827086     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1131.8801479780359     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.88
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2050.819       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.50E-01  7.12E-01  9.58E-01  1.18E+00  8.34E-01  1.03E+00  2.02E+00  6.92E-01  8.84E-01  7.91E-01  1.09E+00
 


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
+        8.12E+07
 
 TH 2
+        2.28E+03  3.08E+02
 
 TH 3
+       -3.40E+04  4.48E+07  7.97E+07
 
 TH 4
+        3.53E+02 -4.43E+02  1.00E+04  7.32E+06
 
 TH 5
+        9.24E+07 -2.65E+02 -7.17E+03  2.27E+03  1.05E+08
 
 TH 6
+        5.93E+07  6.59E+07 -1.63E+03  4.94E+02 -6.75E+07  4.33E+07
 
 TH 7
+       -4.76E+06  5.29E+06 -5.59E+02  1.60E+02  5.41E+06  9.62E+01  2.79E+05
 
 TH 8
+        5.93E+02 -1.18E+03  1.73E+04 -8.49E+03  3.41E+03  8.35E+02  2.89E+02  2.12E+07
 
 TH 9
+        1.26E+03 -2.48E+03  3.65E+04 -1.07E+04  7.13E+03  1.77E+03 -5.11E+06 -1.80E+04  9.36E+07
 
 TH10
+        7.27E+07  4.04E+07  7.20E+07 -2.95E+02  8.27E+07  7.07E-01 -4.26E+06 -4.28E+02 -9.39E+02  6.50E+07
 
 TH11
+       -5.66E+02  1.09E+03 -1.62E+04  7.94E+03 -3.17E+03 -7.83E+02 -2.67E+02  5.92E+04  1.70E+04  4.38E+02  1.86E+07
 
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
 #CPUT: Total CPU Time in Seconds,       27.448
Stop Time:
Wed Sep 29 21:10:26 CDT 2021
