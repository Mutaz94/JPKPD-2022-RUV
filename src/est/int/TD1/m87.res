Wed Sep 29 06:45:43 CDT 2021
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
$DATA ../../../../data/int/TD1/dat87.csv ignore=@
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3332.69093674638        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1591E+02  2.4188E+01  4.8668E+00  4.9189E+01  1.7872E+02  3.4484E+01 -3.5924E+01 -2.1103E+02 -4.3976E+01 -1.8346E+01
            -8.4205E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3845.95331497124        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      150
 NPARAMETR:  9.2803E-01  1.0387E+00  1.1060E+00  9.9040E-01  9.9059E-01  9.6110E-01  1.0171E+00  1.0011E+00  1.0856E+00  9.5218E-01
             1.2142E+00
 PARAMETER:  2.5311E-02  1.3796E-01  2.0075E-01  9.0353E-02  9.0547E-02  6.0321E-02  1.1691E-01  1.0107E-01  1.8217E-01  5.1004E-02
             2.9406E-01
 GRADIENT:  -2.1694E+02 -4.0790E+01 -6.5982E+00 -2.4618E+01  2.3639E-01 -5.1475E+01  2.9532E+00  1.0242E+01  1.8044E+00  2.8754E+00
             2.7591E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3848.38668066443        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      328
 NPARAMETR:  9.5736E-01  1.0004E+00  9.6321E-01  1.0011E+00  9.2651E-01  8.9219E-01  9.3974E-01  5.8542E-01  1.0539E+00  1.0062E+00
             1.1903E+00
 PARAMETER:  5.6428E-02  1.0040E-01  6.2520E-02  1.0113E-01  2.3674E-02 -1.4077E-02  3.7846E-02 -4.3543E-01  1.5253E-01  1.0615E-01
             2.7420E-01
 GRADIENT:  -1.6524E+02 -3.8993E+01 -4.5762E+01 -3.8565E+01  1.4664E+01 -7.7931E+01 -5.5282E+00 -4.1506E+00 -1.1167E+01  1.9458E+01
             2.1873E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3877.46452191072        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  9.9472E-01  1.0402E+00  1.0749E+00  9.9962E-01  9.7819E-01  9.8936E-01  1.0782E+00  8.6009E-01  1.1050E+00  9.2962E-01
             1.0677E+00
 PARAMETER:  9.4702E-02  1.3939E-01  1.7225E-01  9.9616E-02  7.7951E-02  8.9300E-02  1.7532E-01 -5.0717E-02  1.9986E-01  2.7025E-02
             1.6552E-01
 GRADIENT:  -4.5629E+01 -1.3319E+01 -5.0648E+00 -5.0396E+00  4.3128E+00 -2.0206E+01  6.3101E+00 -1.8115E+00 -2.2330E-01  9.0629E-01
             4.1594E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3878.17708792072        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      691
 NPARAMETR:  9.9848E-01  1.0503E+00  1.0865E+00  9.9675E-01  9.8667E-01  9.9771E-01  1.0637E+00  8.8369E-01  1.1048E+00  9.3846E-01
             1.0629E+00
 PARAMETER:  9.8479E-02  1.4909E-01  1.8294E-01  9.6747E-02  8.6585E-02  9.7706E-02  1.6179E-01 -2.3650E-02  1.9963E-01  3.6485E-02
             1.6096E-01
 GRADIENT:  -3.6400E+01 -1.0093E+01 -3.8355E+00 -3.5152E+00  3.2627E+00 -1.6278E+01  4.8785E+00 -1.3931E+00 -3.0984E-01  6.5448E-01
             3.2135E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3878.44591617558        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:      863
 NPARAMETR:  9.9895E-01  1.0503E+00  1.0875E+00  9.9622E-01  9.8720E-01  1.0565E+00  1.0628E+00  8.8715E-01  1.1036E+00  9.3843E-01
             1.0638E+00
 PARAMETER:  9.8948E-02  1.4904E-01  1.8391E-01  9.6215E-02  8.7117E-02  1.5500E-01  1.6093E-01 -1.9736E-02  1.9857E-01  3.6454E-02
             1.6181E-01
 GRADIENT:   3.8667E+02  1.0234E+02  7.6876E+00  8.5901E+01  4.9322E+01  1.0443E+02  1.6304E+01 -7.5045E-01  1.9452E+01  4.4855E+00
             3.6933E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3878.50051120920        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1004
 NPARAMETR:  9.9907E-01  1.0501E+00  1.0873E+00  9.9752E-01  9.8635E-01  1.0488E+00  1.0626E+00  8.8704E-01  1.1039E+00  9.3855E-01
             1.0635E+00
 PARAMETER:  9.9074E-02  1.4885E-01  1.8368E-01  9.7513E-02  8.6261E-02  1.4762E-01  1.6073E-01 -1.9861E-02  1.9881E-01  3.6579E-02
             1.6161E-01
 GRADIENT:  -3.1505E+01 -1.0056E+01 -3.4241E+00 -2.8058E+00  2.4787E+00  4.4177E+00  4.7125E+00 -1.1682E+00 -4.1062E-01  6.7469E-01
             3.3634E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3879.09829972683        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1172
 NPARAMETR:  1.0146E+00  1.0542E+00  1.0903E+00  9.9805E-01  9.8619E-01  1.0358E+00  1.0307E+00  9.0245E-01  1.1034E+00  9.4104E-01
             1.0565E+00
 PARAMETER:  1.1448E-01  1.5278E-01  1.8643E-01  9.8051E-02  8.6091E-02  1.3513E-01  1.3022E-01 -2.6430E-03  1.9837E-01  3.9232E-02
             1.5495E-01
 GRADIENT:   4.8523E+02  1.1101E+02  9.6718E+00  9.0481E+01  4.3770E+01  8.2218E+01  9.2039E+00 -4.6607E-02  1.9728E+01  3.6292E+00
             2.1446E+01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -3879.09829972683        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  1.0151E+00  1.0539E+00  1.0906E+00  9.9819E-01  9.8619E-01  1.0364E+00  1.0306E+00  9.0259E-01  1.1030E+00  9.4198E-01
             1.0567E+00
 PARAMETER:  1.1448E-01  1.5278E-01  1.8643E-01  9.8051E-02  8.6091E-02  1.3513E-01  1.3022E-01 -2.6430E-03  1.9837E-01  3.9232E-02
             1.5495E-01
 GRADIENT:  -1.6437E+00  1.7105E+05 -7.0065E+04 -6.9620E-01  3.5214E+01 -3.9639E-01  2.1586E-02 -3.1863E+04  1.3170E+05 -3.2312E-01
            -1.6897E+05
 NUMSIGDIG:         1.9         2.3         2.3         2.3         6.5         1.8         3.0         2.3         2.3         1.5
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1235
 NO. OF SIG. DIGITS IN FINAL EST.:  1.5
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3811E-03 -2.2594E-02 -1.8779E-02  1.3631E-02 -2.6006E-02
 SE:             2.9931E-02  2.3141E-02  1.5850E-02  2.7657E-02  2.4270E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6320E-01  3.2888E-01  2.3610E-01  6.2212E-01  2.8393E-01

 ETASHRINKSD(%)  1.0000E-10  2.2475E+01  4.6900E+01  7.3449E+00  1.8694E+01
 ETASHRINKVR(%)  1.0000E-10  3.9898E+01  7.1804E+01  1.4150E+01  3.3893E+01
 EBVSHRINKSD(%)  2.7447E-01  2.2624E+01  4.7980E+01  8.0828E+00  1.8190E+01
 EBVSHRINKVR(%)  5.4820E-01  4.0130E+01  7.2939E+01  1.5512E+01  3.3072E+01
 RELATIVEINF(%)  9.9450E+01  3.2558E+01  1.7483E+01  6.1985E+01  2.8479E+01
 EPSSHRINKSD(%)  2.1512E+01
 EPSSHRINKVR(%)  3.8396E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3879.0982997268325     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2225.0089399584217     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.31
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3879.098       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.05E+00  1.09E+00  9.98E-01  9.86E-01  1.04E+00  1.03E+00  9.02E-01  1.10E+00  9.41E-01  1.06E+00
 


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
+        4.84E+07
 
 TH 2
+        6.19E+02  2.52E+07
 
 TH 3
+        2.77E+07  1.08E+04  3.16E+07
 
 TH 4
+       -5.63E+07 -4.06E+07 -6.44E+07  1.31E+08
 
 TH 5
+       -8.32E+03  4.11E+07  6.52E+07 -1.33E+08  1.34E+08
 
 TH 6
+        9.80E+06  1.02E+03  2.30E+07 -9.03E-01 -6.07E-01  1.84E+02
 
 TH 7
+       -4.19E+07 -4.85E+03 -2.39E+07  4.88E+07 -8.67E+07 -1.54E-01  3.63E+07
 
 TH 8
+       -6.23E+07 -3.69E+04 -3.56E+07  7.25E+07  1.07E+04  1.26E+07  5.39E+07  8.01E+07
 
 TH 9
+       -2.57E+07  1.85E+07 -2.94E+07  5.98E+07 -6.05E+07  7.57E+02  2.22E+07  3.31E+07  1.36E+07
 
 TH10
+       -5.98E+07 -6.02E+03 -3.41E+07  6.95E+07 -1.41E+08  9.74E-02  5.17E+07  7.69E+07  3.17E+07  7.38E+07
 
 TH11
+       -6.17E+02 -3.53E+03  1.96E+07 -3.99E+07  4.04E+07 -1.01E+03  4.82E+03  3.64E+04 -1.82E+07  5.90E+03  4.88E+07
 
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
 #CPUT: Total CPU Time in Seconds,       47.446
Stop Time:
Wed Sep 29 06:46:32 CDT 2021
