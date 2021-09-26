Fri Sep 24 19:33:31 CDT 2021
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
$DATA ../../../../data/int/B/dat49.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3205.14342768162        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4753E+02 -4.9362E+01  9.3900E+01 -5.6056E+01  3.3608E+01 -7.6426E+00 -4.5677E+01 -2.5291E+02 -6.2296E+01 -1.7863E-02
            -8.8755E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3571.52696551940        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  8.8131E-01  1.0599E+00  1.1299E+00  9.7902E-01  1.0110E+00  1.0390E+00  1.3821E+00  8.7649E-01  9.7666E-01  9.2419E-01
             1.7217E+00
 PARAMETER: -2.6350E-02  1.5819E-01  2.2211E-01  7.8797E-02  1.1096E-01  1.3830E-01  4.2358E-01 -3.1833E-02  7.6385E-02  2.1163E-02
             6.4330E-01
 GRADIENT:  -1.5629E+02  9.6254E+00  6.2676E+01 -2.4856E+01 -7.8282E+01 -4.1917E+00  1.5262E+00  5.8792E-01  9.6796E+00  8.5087E+00
             7.1292E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3642.75923136511        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  8.9835E-01  9.2405E-01  6.5882E-01  1.0651E+00  7.6737E-01  1.0170E+00  1.4443E+00  2.4485E-01  9.1841E-01  4.1448E-01
             1.4938E+00
 PARAMETER: -7.1937E-03  2.1010E-02 -3.1731E-01  1.6306E-01 -1.6479E-01  1.1684E-01  4.6763E-01 -1.3071E+00  1.4892E-02 -7.8073E-01
             5.0133E-01
 GRADIENT:  -1.1406E+02  7.2646E+01 -1.6905E+02  5.1059E+01  3.2706E+01 -8.5402E+00 -4.0954E+01 -1.7719E+00 -2.9508E+01 -1.9181E+01
             5.0838E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3776.45930405966        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.2023E-01  8.5955E-01  7.6438E-01  1.0906E+00  7.7253E-01  1.0231E+00  1.5682E+00  4.3182E-01  1.0811E+00  8.2037E-01
             1.0880E+00
 PARAMETER:  1.6871E-02 -5.1344E-02 -1.6869E-01  1.8675E-01 -1.5808E-01  1.2281E-01  5.4994E-01 -7.3975E-01  1.7798E-01 -9.8002E-02
             1.8430E-01
 GRADIENT:  -3.9297E+01  2.9293E+01 -2.2186E+01  3.6781E+01 -4.2173E+01  3.1037E+00  3.9240E+01 -4.1927E+00  2.4422E+01  8.7146E+00
             9.3436E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3783.38851648989        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:      361
 NPARAMETR:  9.3762E-01  8.4819E-01  7.9154E-01  1.0659E+00  7.9645E-01  1.0080E+00  1.3439E+00  4.4865E-01  1.0008E+00  8.1305E-01
             1.0818E+00
 PARAMETER:  3.5592E-02 -6.4646E-02 -1.3378E-01  1.6383E-01 -1.2759E-01  1.0801E-01  3.9557E-01 -7.0152E-01  1.0076E-01 -1.0696E-01
             1.7859E-01
 GRADIENT:  -5.3973E-01 -4.6890E+00 -4.4992E+00 -4.3844E+00 -2.7180E+01 -1.7170E+00  1.9616E+00 -5.8113E+00 -1.5181E+00 -1.5346E+00
             8.2995E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3785.26456143159        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      436
 NPARAMETR:  9.3735E-01  9.0700E-01  8.3992E-01  1.0438E+00  8.5820E-01  1.0128E+00  1.3100E+00  4.4865E-01  1.0055E+00  8.5817E-01
             1.0818E+00
 PARAMETER:  3.5301E-02  2.3843E-03 -7.4443E-02  1.4291E-01 -5.2921E-02  1.1267E-01  3.7004E-01 -7.0152E-01  1.0551E-01 -5.2950E-02
             1.7859E-01
 GRADIENT:  -9.0652E-01 -2.1362E-01  1.3629E-01 -4.1972E-03 -3.0268E-03  4.9605E-01 -1.6644E-01 -8.8131E+00  3.6535E-02 -1.6092E-01
             7.6710E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3785.75798891049        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      596
 NPARAMETR:  9.5407E-01  9.2112E-01  8.4761E-01  1.0453E+00  8.6899E-01  1.0287E+00  1.3308E+00  4.4865E-01  1.0062E+00  8.6207E-01
             1.0818E+00
 PARAMETER:  5.2982E-02  1.7831E-02 -6.5335E-02  1.4435E-01 -4.0423E-02  1.2827E-01  3.8578E-01 -7.0152E-01  1.0614E-01 -4.8416E-02
             1.7859E-01
 GRADIENT:   3.2296E-01  4.7713E-01 -7.0627E-01 -4.1321E-01 -9.5401E-02  5.5894E-01  1.0490E-01 -9.5043E+00 -2.6642E-01  1.6760E-01
             7.3982E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3787.68777012623        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  9.5339E-01  9.2036E-01  8.4803E-01  1.0454E+00  8.6888E-01  1.0265E+00  1.3295E+00  5.3297E-01  1.0068E+00  8.6082E-01
             1.0750E+00
 PARAMETER:  5.2273E-02  1.7013E-02 -6.4833E-02  1.4441E-01 -4.0550E-02  1.2620E-01  3.8477E-01 -5.2930E-01  1.0679E-01 -4.9873E-02
             1.7236E-01
 GRADIENT:   3.6265E+01  5.8127E+00 -9.3474E+00  1.2446E+01  3.8351E+00  7.0932E+00  3.7019E+00 -7.1674E+00  1.3403E+00  2.5603E+00
             7.0933E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3787.88251342235        NO. OF FUNC. EVALS.: 112
 CUMULATIVE NO. OF FUNC. EVALS.:      866
 NPARAMETR:  9.5324E-01  9.2036E-01  8.4855E-01  1.0454E+00  8.6901E-01  1.0264E+00  1.3295E+00  5.4061E-01  1.0067E+00  8.5971E-01
             1.0738E+00
 PARAMETER:  5.2109E-02  1.7010E-02 -6.4227E-02  1.4436E-01 -4.0397E-02  1.2608E-01  3.8481E-01 -5.1505E-01  1.0670E-01 -5.1165E-02
             1.7117E-01
 GRADIENT:  -1.3768E+00 -5.7329E-01 -1.0247E+01 -9.7707E-01 -1.4359E+00 -3.9935E-01 -2.5428E-01 -7.0853E+00  1.0910E-02  2.1835E+00
             6.8925E+01

0ITERATION NO.:   41    OBJECTIVE VALUE:  -3787.88251342235        NO. OF FUNC. EVALS.:  42
 CUMULATIVE NO. OF FUNC. EVALS.:      908
 NPARAMETR:  9.5324E-01  9.2036E-01  8.4856E-01  1.0454E+00  8.6901E-01  1.0264E+00  1.3295E+00  5.4089E-01  1.0067E+00  8.5969E-01
             1.0737E+00
 PARAMETER:  5.2109E-02  1.7010E-02 -6.4227E-02  1.4436E-01 -4.0397E-02  1.2608E-01  3.8481E-01 -5.1505E-01  1.0670E-01 -5.1165E-02
             1.7117E-01
 GRADIENT:  -1.3391E+00 -5.5010E-01 -1.0212E+01 -9.0655E-01 -1.3902E+00 -3.9154E-01 -2.3583E-01 -7.0822E+00  9.7778E-03  2.1796E+00
             6.8839E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      908
 NO. OF SIG. DIGITS UNREPORTABLE

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0485E-03 -1.8704E-02 -1.0100E-02  1.1590E-02 -2.4026E-02
 SE:             2.9906E-02  2.5285E-02  1.4315E-02  2.7970E-02  2.3408E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7203E-01  4.5948E-01  4.8044E-01  6.7859E-01  3.0472E-01

 ETASHRINKSD(%)  1.0000E-10  1.5292E+01  5.2044E+01  6.2972E+00  2.1579E+01
 ETASHRINKVR(%)  1.0000E-10  2.8246E+01  7.7003E+01  1.2198E+01  3.8502E+01
 EBVSHRINKSD(%)  2.9651E-01  1.5035E+01  5.6858E+01  7.2703E+00  2.2162E+01
 EBVSHRINKVR(%)  5.9213E-01  2.7810E+01  8.1388E+01  1.4012E+01  3.9412E+01
 RELATIVEINF(%)  9.9406E+01  3.5960E+01  8.5046E+00  6.0895E+01  1.8975E+01
 EPSSHRINKSD(%)  2.3395E+01
 EPSSHRINKVR(%)  4.1317E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3787.8825134223498     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2133.7931536539390     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.45
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.35
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3787.883       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.53E-01  9.20E-01  8.49E-01  1.05E+00  8.69E-01  1.03E+00  1.33E+00  5.41E-01  1.01E+00  8.60E-01  1.07E+00
 


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
+        1.15E+03
 
 TH 2
+       -1.58E+00  5.32E+02
 
 TH 3
+        9.51E-01 -3.41E+01  9.61E+02
 
 TH 4
+        6.01E+09  1.51E+02 -3.20E+01  8.55E+02
 
 TH 5
+       -1.60E+00 -3.65E+02 -6.32E+02  2.60E+02  1.18E+03
 
 TH 6
+        2.34E+00  3.75E-01 -1.76E+00 -1.27E+00 -1.15E+00  1.87E+02
 
 TH 7
+        9.95E-01  9.42E+00 -1.61E-01  1.01E+01 -3.56E+00 -1.81E-01  6.03E+01
 
 TH 8
+       -3.26E+09  3.37E+09 -3.66E+09  2.06E+09 -3.57E+09 -5.82E-01  6.07E+08  1.11E+09
 
 TH 9
+       -3.47E-01 -1.76E+00  2.28E+01  1.53E+01  6.34E+00  2.54E+00  1.11E+01  2.89E+09  1.46E+02
 
 TH10
+       -1.98E+00 -1.34E+01 -2.65E+01 -1.06E+01 -2.37E+01  1.56E-01  2.92E+01  3.61E+09 -4.44E+00  9.81E+01
 
 TH11
+        4.93E+09 -5.11E+09 -1.61E+01 -3.12E+09 -3.52E+01  2.17E+00 -9.19E+08 -1.69E+09 -4.38E+09  1.71E+01  2.56E+09
 
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
 #CPUT: Total CPU Time in Seconds,       32.908
Stop Time:
Fri Sep 24 19:34:07 CDT 2021
