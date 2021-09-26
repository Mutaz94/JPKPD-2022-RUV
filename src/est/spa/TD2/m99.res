Sat Sep 25 13:56:56 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1682.72927170826        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.0270E+01 -6.3804E+01 -2.7476E+01 -4.4392E+01  5.2158E+01  3.4619E+01 -1.2854E+01  2.7427E+00  2.1131E+00  4.9919E+00
             7.8454E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1687.96889890505        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9981E-01  1.0979E+00  9.6691E-01  9.8200E-01  9.9230E-01  8.9942E-01  1.2078E+00  1.0034E+00  9.5829E-01  8.8846E-01
             1.0016E+00
 PARAMETER:  9.9807E-02  1.9338E-01  6.6349E-02  8.1834E-02  9.2274E-02 -6.0093E-03  2.8878E-01  1.0338E-01  5.7394E-02 -1.8261E-02
             1.0159E-01
 GRADIENT:   2.4279E+01  8.7490E+00 -4.3177E+00  1.5967E+01  2.6963E+01 -4.9704E+00  6.0606E+00 -5.5876E-01  1.9448E+00 -1.8185E+00
             5.3669E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1688.49516194493        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0016E+00  1.0497E+00  8.1170E-01  1.0021E+00  8.6973E-01  9.0854E-01  1.2862E+00  9.7868E-01  9.1508E-01  6.8482E-01
             9.5985E-01
 PARAMETER:  1.0163E-01  1.4847E-01 -1.0863E-01  1.0214E-01 -3.9574E-02  4.0886E-03  3.5169E-01  7.8452E-02  1.1260E-02 -2.7860E-01
             5.9027E-02
 GRADIENT:   3.0923E+01  1.4115E+01 -2.1128E+00  2.3120E+01  1.1433E+01 -1.4008E+00  7.1270E+00  2.0487E-01  1.9016E+00 -5.0705E+00
            -1.0345E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1689.02310879319        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.9502E-01  1.1017E+00  6.9068E-01  9.5048E-01  8.3553E-01  9.1047E-01  1.2052E+00  6.9625E-01  9.2192E-01  7.1488E-01
             9.7617E-01
 PARAMETER:  9.5007E-02  1.9687E-01 -2.7007E-01  4.9216E-02 -7.9693E-02  6.2008E-03  2.8668E-01 -2.6204E-01  1.8705E-02 -2.3563E-01
             7.5885E-02
 GRADIENT:   7.7234E+00  5.5428E+00 -3.6471E+00  1.0598E+01  6.2245E+00 -1.4174E+00  2.3417E+00  5.0825E-01  8.1677E-02 -8.8759E-01
            -1.9628E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1689.36892815544        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      347
 NPARAMETR:  1.0040E+00  1.1139E+00  7.0994E-01  9.4525E-01  8.5326E-01  9.1951E-01  1.1859E+00  6.8403E-01  9.3229E-01  7.5544E-01
             9.8497E-01
 PARAMETER:  1.0397E-01  2.0784E-01 -2.4258E-01  4.3696E-02 -5.8692E-02  1.6089E-02  2.7050E-01 -2.7976E-01  2.9889E-02 -1.8046E-01
             8.4857E-02
 GRADIENT:  -1.1374E+01 -2.8925E+00 -2.4743E+00  9.8228E-01  1.9449E+00 -1.6172E+00  3.4776E-01  4.7874E-01 -3.7082E-01  6.6364E-01
             1.6712E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1689.64604383144        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      522
 NPARAMETR:  1.0085E+00  1.3021E+00  5.9317E-01  8.2280E-01  8.7978E-01  9.2391E-01  1.0539E+00  4.5978E-01  1.0223E+00  7.6154E-01
             9.8247E-01
 PARAMETER:  1.0851E-01  3.6395E-01 -4.2227E-01 -9.5036E-02 -2.8082E-02  2.0864E-02  1.5252E-01 -6.7700E-01  1.2203E-01 -1.7241E-01
             8.2313E-02
 GRADIENT:  -1.6331E+00  9.2127E-01 -3.7510E-01  1.7311E+00 -4.9207E-01 -2.7869E-01  1.9574E-01  1.0156E-01  1.2242E-01  3.8898E-01
             4.7181E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1689.68746912342        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      697
 NPARAMETR:  1.0093E+00  1.3982E+00  5.3906E-01  7.5804E-01  9.0068E-01  9.2486E-01  9.9773E-01  3.3063E-01  1.0775E+00  7.6437E-01
             9.8330E-01
 PARAMETER:  1.0921E-01  4.3515E-01 -5.1793E-01 -1.7702E-01 -4.6089E-03  2.1882E-02  9.7731E-02 -1.0067E+00  1.7467E-01 -1.6870E-01
             8.3157E-02
 GRADIENT:  -1.7554E-01 -4.7297E-01 -2.7040E-01 -2.4238E-01  3.0718E-01  1.6217E-03  6.3878E-03  7.2652E-02  6.3096E-02  1.2463E-01
             5.2302E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1689.70910185745        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      875
 NPARAMETR:  1.0093E+00  1.3986E+00  5.3010E-01  7.5721E-01  8.9586E-01  9.2491E-01  9.9926E-01  1.3834E-01  1.0772E+00  7.6607E-01
             9.8305E-01
 PARAMETER:  1.0925E-01  4.3550E-01 -5.3469E-01 -1.7812E-01 -9.9746E-03  2.1946E-02  9.9259E-02 -1.8780E+00  1.7435E-01 -1.6648E-01
             8.2906E-02
 GRADIENT:  -1.6559E-01  5.1569E-01 -9.3359E-02  6.7754E-01  2.6712E-01 -1.8362E-02  2.1263E-03  6.6583E-03 -5.9981E-02 -4.2622E-02
            -1.4766E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1689.71339422962        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1050
 NPARAMETR:  1.0093E+00  1.4092E+00  5.2296E-01  7.4954E-01  8.9758E-01  9.2499E-01  9.9308E-01  2.4645E-02  1.0845E+00  7.6624E-01
             9.8338E-01
 PARAMETER:  1.0930E-01  4.4303E-01 -5.4824E-01 -1.8830E-01 -8.0573E-03  2.2022E-02  9.3052E-02 -3.6032E+00  1.8110E-01 -1.6626E-01
             8.3243E-02
 GRADIENT:  -3.3549E-02 -1.0891E-01 -6.7181E-02  7.6506E-02  7.7855E-02  3.2826E-03 -6.4935E-02  2.6584E-04  5.1136E-02  9.3525E-03
            -1.4620E-02

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1689.71352832005        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1177
 NPARAMETR:  1.0094E+00  1.4093E+00  5.2306E-01  7.4952E-01  8.9763E-01  9.2498E-01  9.9345E-01  1.0000E-02  1.0842E+00  7.6629E-01
             9.8343E-01
 PARAMETER:  1.0931E-01  4.4307E-01 -5.4807E-01 -1.8832E-01 -7.9972E-03  2.2011E-02  9.3425E-02 -4.5852E+00  1.8086E-01 -1.6620E-01
             8.3290E-02
 GRADIENT:  -3.8690E-03 -2.8173E-03 -6.5782E-04  5.2220E-03  6.0397E-03 -1.4730E-04 -4.9432E-03  0.0000E+00  1.2477E-03 -1.9582E-03
            -2.6766E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1177
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0238E-04 -1.4352E-02 -3.5840E-04  1.1974E-02 -2.2824E-02
 SE:             2.9833E-02  2.4719E-02  1.4771E-04  2.3810E-02  2.0888E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9726E-01  5.6151E-01  1.5249E-02  6.1505E-01  2.7452E-01

 ETASHRINKSD(%)  5.5650E-02  1.7187E+01  9.9505E+01  2.0233E+01  3.0023E+01
 ETASHRINKVR(%)  1.1127E-01  3.1420E+01  9.9998E+01  3.6373E+01  5.1032E+01
 EBVSHRINKSD(%)  4.8113E-01  1.7020E+01  9.9570E+01  2.0773E+01  2.9879E+01
 EBVSHRINKVR(%)  9.5994E-01  3.1143E+01  9.9998E+01  3.7231E+01  5.0830E+01
 RELATIVEINF(%)  9.9012E+01  5.0260E+00  1.6656E-04  4.3499E+00  5.6591E+00
 EPSSHRINKSD(%)  4.4698E+01
 EPSSHRINKVR(%)  6.9417E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1689.7135283200539     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -954.56270175631573     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.91
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1689.714       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.41E+00  5.23E-01  7.50E-01  8.98E-01  9.25E-01  9.93E-01  1.00E-02  1.08E+00  7.66E-01  9.83E-01
 


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
+        1.26E+03
 
 TH 2
+       -5.89E+00  4.02E+02
 
 TH 3
+        9.69E+00  2.47E+02  8.02E+02
 
 TH 4
+       -1.57E+01  2.66E+02 -5.31E+02  1.07E+03
 
 TH 5
+       -6.20E+00 -3.40E+02 -8.08E+02  5.15E+02  1.16E+03
 
 TH 6
+       -1.08E+00 -8.50E-01  1.66E+00 -4.00E+00 -2.50E+00  2.26E+02
 
 TH 7
+        5.00E-01  2.08E+01 -3.48E+01 -9.86E+00  7.11E+00  1.38E+00  9.69E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.12E+00 -2.00E+01 -4.82E+01  5.69E+01 -1.08E+01 -3.40E-01  1.46E+01  0.00E+00  7.61E+01
 
 TH10
+       -1.88E+00 -1.48E+01 -5.09E+01 -2.30E+01 -7.85E+01  3.84E-01  1.95E+01  0.00E+00  1.59E+01  9.00E+01
 
 TH11
+       -9.81E+00 -1.35E+01 -2.95E+01  8.53E-01 -6.10E+00  2.66E+00  8.92E+00  0.00E+00  8.69E+00  2.22E+01  2.11E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.546
Stop Time:
Sat Sep 25 13:57:18 CDT 2021
