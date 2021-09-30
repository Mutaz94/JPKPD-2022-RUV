Wed Sep 29 23:33:26 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat63.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m63.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1484.81424152712        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3916E+02 -5.5808E+01  3.5959E+01 -6.2744E+01  1.0523E+02  1.7041E+01 -2.2330E+01 -4.5049E+01 -2.6360E+01 -1.7227E+01
            -1.1592E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1838.79533352114        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0364E+00  1.0733E+00  1.2353E+00  1.0440E+00  1.0726E+00  1.1002E+00  1.0330E+00  8.1665E-01  1.1512E+00  6.4243E-01
             1.9755E+00
 PARAMETER:  1.3580E-01  1.7077E-01  3.1132E-01  1.4308E-01  1.7010E-01  1.9552E-01  1.3246E-01 -1.0254E-01  2.4083E-01 -3.4249E-01
             7.8081E-01
 GRADIENT:   1.3359E+02 -7.3086E-01  1.1973E+01 -7.7041E+00  9.3248E+00  4.0757E+01  3.6452E+00  8.0032E-01  1.3140E+01 -1.3811E+00
            -2.9606E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1843.89097457382        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0467E+00  9.0474E-01  5.5017E-01  1.1338E+00  6.5634E-01  1.0109E+00  1.2018E+00  4.4911E-01  1.0333E+00  2.7996E-01
             1.9169E+00
 PARAMETER:  1.4565E-01 -1.0808E-04 -4.9754E-01  2.2557E-01 -3.2108E-01  1.1084E-01  2.8384E-01 -7.0050E-01  1.3273E-01 -1.1731E+00
             7.5069E-01
 GRADIENT:   1.5242E+02 -2.1373E+01 -4.7928E+01  1.1674E+02  1.0141E+02 -3.1519E+00 -1.0417E+01  3.2598E+00  1.9111E+01  2.8407E-01
            -2.0475E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1849.49185649961        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      338
 NPARAMETR:  1.0489E+00  7.1971E-01  6.4535E-01  1.2436E+00  6.1781E-01  1.0507E+00  1.5813E+00  1.5253E-01  8.9707E-01  4.2092E-01
             1.9782E+00
 PARAMETER:  1.4770E-01 -2.2891E-01 -3.3797E-01  3.1799E-01 -3.8157E-01  1.4947E-01  5.5825E-01 -1.7804E+00 -8.6208E-03 -7.6532E-01
             7.8219E-01
 GRADIENT:   2.5000E+00  2.7980E+01  2.8384E+01  2.7666E+01 -4.2838E+01  3.1072E+00 -1.1885E+00  2.3899E-01  1.0324E-01 -6.8406E-01
             1.6335E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1851.89136069830        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      513
 NPARAMETR:  1.0454E+00  4.9001E-01  5.6786E-01  1.3092E+00  5.1703E-01  1.0384E+00  1.9628E+00  3.7211E-02  8.5993E-01  4.7464E-01
             1.9383E+00
 PARAMETER:  1.4440E-01 -6.1332E-01 -4.6588E-01  3.6942E-01 -5.5965E-01  1.3768E-01  7.7436E-01 -3.1911E+00 -5.0899E-02 -6.4520E-01
             7.6180E-01
 GRADIENT:   2.3291E+00 -5.8849E-01  2.5165E+00 -7.5026E+00 -2.2274E+00 -6.2044E-01 -1.9134E+00  3.0424E-02  9.1557E-01 -1.2112E-01
             9.2100E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1851.98549127581        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      689
 NPARAMETR:  1.0438E+00  4.4508E-01  5.4696E-01  1.3265E+00  4.9506E-01  1.0398E+00  2.1088E+00  1.5900E-02  8.5213E-01  4.8581E-01
             1.9296E+00
 PARAMETER:  1.4288E-01 -7.0951E-01 -5.0339E-01  3.8253E-01 -6.0307E-01  1.3907E-01  8.4610E-01 -4.0414E+00 -6.0022E-02 -6.2193E-01
             7.5734E-01
 GRADIENT:   7.3807E-01 -1.0310E+00 -7.2768E+00 -5.1875E-01  7.7589E+00  4.9934E-03 -8.9696E-02  6.4174E-03  5.4960E-01  4.2217E-01
             7.6097E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1852.02709090919        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      865
 NPARAMETR:  1.0385E+00  4.2123E-01  5.5950E-01  1.3384E+00  4.9594E-01  1.0369E+00  2.1813E+00  1.0000E-02  8.4835E-01  4.9261E-01
             1.9306E+00
 PARAMETER:  1.3781E-01 -7.6458E-01 -4.8072E-01  3.9151E-01 -6.0130E-01  1.3628E-01  8.7991E-01 -7.4099E+00 -6.4462E-02 -6.0803E-01
             7.5784E-01
 GRADIENT:  -7.8791E+00 -3.3840E-01  5.8906E-01 -7.1163E+00 -1.7847E+00 -6.8720E-01 -3.8171E-01  0.0000E+00 -2.9036E-01 -1.3607E-01
            -2.3027E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1852.04465126302        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1045
 NPARAMETR:  1.0324E+00  3.5417E-01  5.7947E-01  1.3740E+00  4.9268E-01  1.0327E+00  2.3835E+00  1.0000E-02  8.4950E-01  5.2432E-01
             1.9291E+00
 PARAMETER:  1.3187E-01 -9.3797E-01 -4.4565E-01  4.1774E-01 -6.0790E-01  1.3215E-01  9.6857E-01 -1.5459E+01 -6.3105E-02 -5.4565E-01
             7.5705E-01
 GRADIENT:  -1.5501E+01 -2.1145E-01  3.8225E+00 -1.1197E+01 -6.1946E+00 -1.2475E+00 -1.0195E+00  0.0000E+00 -4.4612E-01  4.1462E-01
            -1.6611E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1852.14869057009        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1223
 NPARAMETR:  1.0345E+00  3.0158E-01  5.8302E-01  1.4001E+00  4.8721E-01  1.0324E+00  2.6910E+00  1.0000E-02  8.4353E-01  5.2387E-01
             1.9292E+00
 PARAMETER:  1.3391E-01 -1.0987E+00 -4.3954E-01  4.3652E-01 -6.1906E-01  1.3187E-01  1.0899E+00 -2.1827E+01 -7.0155E-02 -5.4651E-01
             7.5709E-01
 GRADIENT:  -7.7901E+00 -2.0099E+00 -1.1534E+00 -5.7581E+00  1.5157E+00 -7.7294E-01 -1.0409E+00  0.0000E+00  5.0582E-01  1.9059E-01
            -9.9690E-01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1852.18051562024        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1351
 NPARAMETR:  1.0388E+00  3.0608E-01  5.8275E-01  1.4009E+00  4.8665E-01  1.0342E+00  2.6777E+00  1.0000E-02  8.4408E-01  5.2207E-01
             1.9319E+00
 PARAMETER:  1.3806E-01 -1.0839E+00 -4.3999E-01  4.3709E-01 -6.2021E-01  1.3364E-01  1.0849E+00 -2.1383E+01 -6.9509E-02 -5.4996E-01
             7.5850E-01
 GRADIENT:   3.9470E-01  1.7184E-01  1.9282E+00 -2.0105E+00 -4.3314E+00  1.5272E-03  2.9789E-01  0.0000E+00 -3.6765E-03 -1.2865E-01
            -5.0201E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1351
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2815E-03  4.5402E-02 -3.7699E-04 -2.5877E-02  1.4534E-02
 SE:             2.9645E-02  1.9196E-02  2.1955E-04  2.5994E-02  1.6967E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3865E-01  1.8021E-02  8.5961E-02  3.1948E-01  3.9167E-01

 ETASHRINKSD(%)  6.8412E-01  3.5691E+01  9.9264E+01  1.2918E+01  4.3157E+01
 ETASHRINKVR(%)  1.3636E+00  5.8643E+01  9.9995E+01  2.4168E+01  6.7689E+01
 EBVSHRINKSD(%)  1.1273E+00  4.4410E+01  9.9136E+01  9.8816E+00  3.9629E+01
 EBVSHRINKVR(%)  2.2419E+00  6.9097E+01  9.9993E+01  1.8787E+01  6.3553E+01
 RELATIVEINF(%)  9.7103E+01  9.1451E+00  4.1745E-04  3.4000E+01  1.9400E+00
 EPSSHRINKSD(%)  2.8321E+01
 EPSSHRINKVR(%)  4.8621E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1852.1805156202370     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -933.24198241556428     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.20
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1852.181       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  3.06E-01  5.83E-01  1.40E+00  4.87E-01  1.03E+00  2.68E+00  1.00E-02  8.44E-01  5.22E-01  1.93E+00
 


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
+        9.43E+02
 
 TH 2
+       -1.28E+01  1.00E+03
 
 TH 3
+       -1.46E+01  7.81E+01  2.91E+03
 
 TH 4
+       -1.35E+01  1.39E+02 -2.26E+02  6.77E+02
 
 TH 5
+        4.15E+01 -5.51E+02 -4.37E+03  5.14E+01  7.00E+03
 
 TH 6
+        2.65E+00  1.50E+01 -1.64E+01 -8.85E+00  2.34E+01  1.78E+02
 
 TH 7
+        4.13E+00  1.06E+02 -6.17E+01 -1.70E+01  6.53E+01  3.19E+00  1.68E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -3.13E+00 -1.46E+02  1.16E+02  2.28E+01 -1.13E+02 -8.23E+00 -1.61E+01  0.00E+00  2.43E+02
 
 TH10
+       -2.87E+00 -2.61E+01 -7.13E+01 -1.29E+01  9.58E+01 -1.56E+00  1.17E-01  0.00E+00  4.06E+00  1.16E+02
 
 TH11
+       -1.14E+01 -2.58E+01 -2.71E+01 -1.16E+01  1.42E+01  2.08E+00 -9.93E-01  0.00E+00  1.06E+01  3.08E+01  1.23E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       29.508
Stop Time:
Wed Sep 29 23:33:58 CDT 2021
