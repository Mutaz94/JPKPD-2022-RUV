Thu Sep 30 01:48:42 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat1.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2103.57737860323        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3108E+02  1.0305E+01  1.7473E+01  4.6304E+01 -3.4637E+01  5.0940E+01 -1.5235E+00  5.1090E-01  3.0831E+01 -1.7049E+01
            -3.9870E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2110.49134155025        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0047E+00  1.1903E+00  1.0314E+00  9.1193E-01  1.2058E+00  9.8659E-01  1.0628E+00  9.7426E-01  8.0194E-01  1.2627E+00
             1.0039E+00
 PARAMETER:  1.0473E-01  2.7420E-01  1.3096E-01  7.8104E-03  2.8717E-01  8.6497E-02  1.6092E-01  7.3919E-02 -1.2072E-01  3.3329E-01
             1.0390E-01
 GRADIENT:   1.1774E+00 -3.9499E+00 -6.5239E+00  5.9116E+00  2.1443E+01 -2.4548E+00 -9.2288E-01 -2.2508E+00 -1.6058E+01 -8.4793E+00
            -3.4337E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2111.23686084107        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0074E+00  1.2925E+00  1.0513E+00  8.6604E-01  1.2800E+00  9.9100E-01  1.0043E+00  1.0097E+00  8.9864E-01  1.4698E+00
             1.0053E+00
 PARAMETER:  1.0737E-01  3.5654E-01  1.4998E-01 -4.3823E-02  3.4689E-01  9.0955E-02  1.0429E-01  1.0968E-01 -6.8717E-03  4.8513E-01
             1.0530E-01
 GRADIENT:   7.2216E+00  2.0057E+01 -4.3687E+00  2.7990E+01  9.3643E+00 -6.6080E-01  3.3989E+00 -2.2918E+00 -6.5073E+00  1.0446E+01
            -3.0676E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2113.82486467874        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0073E+00  1.5157E+00  8.1522E-01  7.1572E-01  1.2430E+00  9.9543E-01  8.3119E-01  9.1512E-01  1.0843E+00  1.3177E+00
             1.0366E+00
 PARAMETER:  1.0728E-01  5.1589E-01 -1.0430E-01 -2.3446E-01  3.1750E-01  9.5418E-02 -8.4901E-02  1.1303E-02  1.8094E-01  3.7588E-01
             1.3594E-01
 GRADIENT:   9.7950E-01  2.6487E+01  5.7216E+00  1.6084E+01 -1.1773E+01 -1.0696E-02 -2.9781E+00 -5.7594E-02 -2.5474E+00  9.7592E-01
            -4.4054E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2115.40075901977        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0087E+00  1.8685E+00  4.4999E-01  4.7537E-01  1.2891E+00  1.0008E+00  7.4327E-01  4.7910E-01  1.3795E+00  1.2945E+00
             1.0288E+00
 PARAMETER:  1.0863E-01  7.2513E-01 -6.9852E-01 -6.4366E-01  3.5394E-01  1.0083E-01 -1.9670E-01 -6.3584E-01  4.2172E-01  3.5815E-01
             1.2840E-01
 GRADIENT:  -8.8595E-01  1.7435E+01 -5.4041E-02  8.2996E+00 -5.5248E+00  8.2971E-01 -8.8374E-01  8.0431E-01 -3.1617E-01  6.8176E-01
            -4.9579E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2115.63285175763        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      908
 NPARAMETR:  1.0106E+00  1.8801E+00  4.3081E-01  4.4073E-01  1.3112E+00  9.9915E-01  7.4038E-01  1.6719E-01  1.4252E+00  1.3011E+00
             1.0375E+00
 PARAMETER:  1.1053E-01  7.3130E-01 -7.4209E-01 -7.1931E-01  3.7094E-01  9.9152E-02 -2.0059E-01 -1.6886E+00  4.5432E-01  3.6325E-01
             1.3677E-01
 GRADIENT:   3.5974E+00 -3.5896E+01  2.1287E+00 -1.5543E+01 -4.5732E+00  3.1830E-01  3.6645E-01  9.6487E-02 -1.9320E+00 -1.0997E+00
             1.8788E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2115.89505563965        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1088             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0100E+00  1.8920E+00  4.2840E-01  4.4674E-01  1.3177E+00  9.9890E-01  7.3791E-01  5.3208E-02  1.4373E+00  1.3159E+00
             1.0353E+00
 PARAMETER:  1.0994E-01  7.3762E-01 -7.4769E-01 -7.0578E-01  3.7586E-01  9.8900E-02 -2.0393E-01 -2.8335E+00  4.6275E-01  3.7454E-01
             1.3474E-01
 GRADIENT:   4.6266E+02  9.0221E+02  2.5113E+00  1.1437E+02  1.9708E+01  4.6835E+01  1.3120E+01  1.3853E-02  1.3374E+01  4.4162E+00
             1.6782E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2115.90706471643        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1274             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0096E+00  1.8917E+00  4.3445E-01  4.5035E-01  1.3169E+00  9.9875E-01  7.3817E-01  1.0000E-02  1.4325E+00  1.3150E+00
             1.0350E+00
 PARAMETER:  1.0958E-01  7.3748E-01 -7.3367E-01 -6.9774E-01  3.7525E-01  9.8753E-02 -2.0358E-01 -4.8579E+00  4.5939E-01  3.7381E-01
             1.3443E-01
 GRADIENT:   4.6064E+02  9.0750E+02  3.5889E+00  1.1605E+02  1.8144E+01  4.6865E+01  1.2841E+01  0.0000E+00  1.2973E+01  3.9350E+00
             9.6342E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -2115.91408196410        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:     1406
 NPARAMETR:  1.0099E+00  1.8865E+00  4.3386E-01  4.4976E-01  1.3184E+00  9.9883E-01  7.3885E-01  1.0000E-02  1.4326E+00  1.3199E+00
             1.0355E+00
 PARAMETER:  1.0985E-01  7.3475E-01 -7.3504E-01 -6.9905E-01  3.7641E-01  9.8829E-02 -2.0266E-01 -4.8579E+00  4.5950E-01  3.7755E-01
             1.3490E-01
 GRADIENT:   2.0291E+00 -1.3189E+01 -8.0206E-01 -1.4734E+00  5.5784E-01  1.6155E-01  4.8803E-02  0.0000E+00  4.0196E-01  6.1982E-01
             3.7813E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1406
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.6739E-04 -3.2444E-02 -2.6861E-04  3.7017E-02 -4.2139E-02
 SE:             2.9886E-02  2.4683E-02  8.6590E-05  2.1849E-02  2.2803E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8485E-01  1.8871E-01  1.9220E-03  9.0224E-02  6.4607E-02

 ETASHRINKSD(%)  1.0000E-10  1.7307E+01  9.9710E+01  2.6802E+01  2.3607E+01
 ETASHRINKVR(%)  1.0000E-10  3.1619E+01  9.9999E+01  4.6421E+01  4.1641E+01
 EBVSHRINKSD(%)  3.5597E-01  1.6758E+01  9.9742E+01  3.0655E+01  1.9199E+01
 EBVSHRINKVR(%)  7.1067E-01  3.0707E+01  9.9999E+01  5.1912E+01  3.4711E+01
 RELATIVEINF(%)  9.9220E+01  6.4671E+00  1.1798E-04  3.9084E+00  2.0548E+01
 EPSSHRINKSD(%)  3.3866E+01
 EPSSHRINKVR(%)  5.6263E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2115.9140819641020     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1196.9755487594293     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.87
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2115.914       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.89E+00  4.34E-01  4.50E-01  1.32E+00  9.99E-01  7.39E-01  1.00E-02  1.43E+00  1.32E+00  1.04E+00
 


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
+        1.08E+03
 
 TH 2
+       -4.57E+00  4.18E+02
 
 TH 3
+        8.20E+00  1.25E+02  3.75E+02
 
 TH 4
+       -9.06E+00  4.11E+02 -3.73E+02  1.30E+03
 
 TH 5
+        1.89E+00 -8.77E+01 -2.01E+02  2.22E+02  2.44E+02
 
 TH 6
+        4.03E-01 -9.09E-01  2.88E+00 -3.24E+00  2.88E-01  1.97E+02
 
 TH 7
+        7.11E-01  7.53E+00 -1.14E+01 -1.25E+01 -6.94E+00 -4.51E-01  1.90E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.81E-01 -1.66E+01 -3.77E+01  7.12E+01  2.63E+00 -3.44E-01  2.00E+01  0.00E+00  3.40E+01
 
 TH10
+        1.54E+00 -1.17E+01 -2.86E+01  1.53E+01 -3.82E+01  6.45E-01  3.32E+00  0.00E+00  6.60E+00  5.16E+01
 
 TH11
+       -7.98E+00 -2.17E+01 -3.98E+01  1.09E+01 -1.43E+00  1.67E+00  1.27E+01  0.00E+00  4.28E+00  9.57E+00  3.78E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.482
Stop Time:
Thu Sep 30 01:49:14 CDT 2021
