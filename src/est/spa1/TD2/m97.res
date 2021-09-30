Thu Sep 30 02:30:37 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat97.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2098.04424518874        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.5430E+02 -3.0532E+01 -4.6317E+01  3.2154E+01  7.1891E+01  3.3893E+01  1.9062E+00  1.6853E+01  2.7394E+01  1.4169E+01
             3.3669E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2104.97449719805        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.9577E-01  1.0600E+00  1.0710E+00  1.0300E+00  9.8048E-01  1.0906E+00  1.0034E+00  9.3723E-01  9.2921E-01  9.3907E-01
             9.3863E-01
 PARAMETER:  9.5765E-02  1.5828E-01  1.6857E-01  1.2951E-01  8.0291E-02  1.8672E-01  1.0336E-01  3.5171E-02  2.6581E-02  3.7139E-02
             3.6671E-02
 GRADIENT:   1.2478E+02  3.7966E+01  9.2824E+00  3.7890E+01 -2.0089E+01  1.0982E+01 -4.6239E+00  9.2026E+00  5.3006E+00  3.3557E+00
            -2.4106E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2111.06941035373        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  9.8190E-01  8.0563E-01  8.8650E-01  1.2194E+00  7.9355E-01  1.0582E+00  1.7013E+00  3.1569E-01  7.8543E-01  7.7848E-01
             9.4402E-01
 PARAMETER:  8.1734E-02 -1.1614E-01 -2.0470E-02  2.9840E-01 -1.3124E-01  1.5655E-01  6.3137E-01 -1.0530E+00 -1.4152E-01 -1.5041E-01
             4.2392E-02
 GRADIENT:   1.0773E+02  6.5751E+01 -1.1200E+01  1.3443E+02  1.0288E+01  2.5290E+00  1.4919E+01  9.7070E-01  2.4774E+00 -1.4085E+00
            -9.3447E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2120.86898384183        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.1335E-01  7.0790E-01  8.2834E-01  1.1837E+00  7.3720E-01  1.0389E+00  1.5657E+00  2.3408E-01  7.4111E-01  7.3797E-01
             9.5262E-01
 PARAMETER:  9.3654E-03 -2.4545E-01 -8.8335E-02  2.6863E-01 -2.0490E-01  1.3815E-01  5.4831E-01 -1.3521E+00 -1.9960E-01 -2.0385E-01
             5.1462E-02
 GRADIENT:  -3.6419E+01  8.5258E+00  1.5758E+01 -1.4141E+01 -2.4344E+01  7.6272E-01 -2.1592E+00  4.9077E-01 -5.9749E+00 -7.7044E-01
             1.3809E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2122.03000731385        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      704             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2913E-01  5.7469E-01  8.7471E-01  1.2633E+00  7.2763E-01  1.0355E+00  1.8252E+00  6.3789E-02  7.3397E-01  7.8312E-01
             9.4983E-01
 PARAMETER:  2.6494E-02 -4.5392E-01 -3.3863E-02  3.3371E-01 -2.1796E-01  1.3488E-01  7.0171E-01 -2.6522E+00 -2.0929E-01 -1.4447E-01
             4.8527E-02
 GRADIENT:   4.2572E+02  6.4794E+01  6.6846E+00  4.4594E+02  2.5779E+01  9.6161E+01  3.0964E+01  5.1488E-02  1.7034E+01  3.3900E-01
             6.3537E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2122.03392934822        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      860
 NPARAMETR:  9.2798E-01  5.7463E-01  8.7407E-01  1.2640E+00  7.2771E-01  1.0337E+00  1.8221E+00  5.5017E-02  7.3393E-01  7.8638E-01
             9.4996E-01
 PARAMETER:  2.5260E-02 -4.5403E-01 -3.4589E-02  3.3424E-01 -2.1785E-01  1.3310E-01  6.9999E-01 -2.8001E+00 -2.0934E-01 -1.4031E-01
             4.8664E-02
 GRADIENT:   4.2238E-01  2.2915E-02  3.9306E-01 -1.5081E+00  3.2426E-02  1.8936E-01  6.6232E-02  1.0425E-02  8.0530E-03 -1.6144E-01
            -4.4407E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2122.03918578299        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1019
 NPARAMETR:  9.2851E-01  5.7419E-01  8.7391E-01  1.2631E+00  7.2750E-01  1.0350E+00  1.8256E+00  1.8785E-02  7.3408E-01  7.8747E-01
             9.4998E-01
 PARAMETER:  2.5828E-02 -4.5479E-01 -3.4781E-02  3.3359E-01 -2.1814E-01  1.3436E-01  7.0191E-01 -3.8747E+00 -2.0913E-01 -1.3894E-01
             4.8681E-02
 GRADIENT:   1.6149E+00 -3.2251E-01  1.0829E+00 -4.0682E+00 -7.3535E-01  7.0267E-01  3.0560E-01  1.2713E-03  1.0211E-01 -4.2296E-03
            -3.0416E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2122.04152691600        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1204
 NPARAMETR:  9.2848E-01  5.7587E-01  8.7292E-01  1.2625E+00  7.2746E-01  1.0346E+00  1.8227E+00  1.0000E-02  7.3398E-01  7.8654E-01
             9.4991E-01
 PARAMETER:  2.5792E-02 -4.5188E-01 -3.5909E-02  3.3309E-01 -2.1819E-01  1.3398E-01  7.0032E-01 -5.8223E+00 -2.0928E-01 -1.4011E-01
             4.8614E-02
 GRADIENT:   1.4801E+00 -5.8830E-02  9.5292E-01 -3.0918E+00 -6.2597E-01  5.3589E-01  3.0656E-01  0.0000E+00 -1.9825E-02 -7.9383E-02
            -8.3559E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -2122.04192958948        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     1309
 NPARAMETR:  9.2852E-01  5.7747E-01  8.7264E-01  1.2616E+00  7.2701E-01  1.0346E+00  1.8208E+00  1.0000E-02  7.3415E-01  7.8582E-01
             9.4984E-01
 PARAMETER:  2.5824E-02 -4.5139E-01 -3.6909E-02  3.3293E-01 -2.1796E-01  1.3402E-01  7.0007E-01 -5.8223E+00 -2.0889E-01 -1.3963E-01
             4.8715E-02
 GRADIENT:  -2.3266E-02 -4.4912E-01 -5.4346E-01  1.0821E+00  1.2590E+00 -4.6934E-03  6.3694E-02  0.0000E+00  2.1002E-02  7.4594E-02
             9.3531E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1309
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1457E-04  2.1159E-02 -4.6215E-04 -1.7947E-02  1.4001E-03
 SE:             2.9887E-02  2.0364E-02  2.2885E-04  2.5066E-02  2.3686E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8092E-01  2.9878E-01  4.3438E-02  4.7399E-01  9.5286E-01

 ETASHRINKSD(%)  1.0000E-10  3.1778E+01  9.9233E+01  1.6027E+01  2.0649E+01
 ETASHRINKVR(%)  1.0000E-10  5.3457E+01  9.9994E+01  2.9486E+01  3.7035E+01
 EBVSHRINKSD(%)  3.0378E-01  3.3702E+01  9.9266E+01  1.4990E+01  1.8731E+01
 EBVSHRINKVR(%)  6.0664E-01  5.6046E+01  9.9995E+01  2.7733E+01  3.3953E+01
 RELATIVEINF(%)  9.8929E+01  6.0744E+00  6.3125E-04  1.1422E+01  8.5635E+00
 EPSSHRINKSD(%)  3.3210E+01
 EPSSHRINKVR(%)  5.5391E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2122.0419295894762     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1203.1033963848035     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.03
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2122.042       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.29E-01  5.76E-01  8.72E-01  1.26E+00  7.28E-01  1.03E+00  1.82E+00  1.00E-02  7.34E-01  7.87E-01  9.50E-01
 


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
+        1.20E+03
 
 TH 2
+       -1.22E+01  4.69E+02
 
 TH 3
+        8.97E+00  2.39E+02  8.41E+02
 
 TH 4
+        1.09E+00  3.78E+02 -3.12E+02  9.48E+02
 
 TH 5
+       -1.88E+00 -5.04E+02 -1.23E+03  3.51E+02  2.22E+03
 
 TH 6
+        6.47E-01 -2.61E+00  1.65E+00 -8.73E-01 -2.43E+00  1.84E+02
 
 TH 7
+        1.26E+00  3.61E+01  2.60E-01 -7.29E+00 -4.52E+00  9.92E-02  1.88E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.16E+00 -1.98E+01 -1.88E+01  1.66E+00  3.05E+01 -4.41E-01  1.22E+01  0.00E+00  2.10E+02
 
 TH10
+       -1.96E-03  4.18E+00 -9.21E+01 -3.22E+01 -2.67E+01  3.16E-01  8.83E+00  0.00E+00  7.00E+00  1.37E+02
 
 TH11
+       -7.51E+00 -1.19E+01 -3.74E+01 -1.03E+01  1.98E+01  1.30E+00  9.06E-01  0.00E+00  1.04E+01  3.08E+01  4.56E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.480
Stop Time:
Thu Sep 30 02:31:07 CDT 2021
