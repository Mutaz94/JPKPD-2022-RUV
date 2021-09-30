Wed Sep 29 21:47:43 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat6.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1650.98661401023        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4147E+02 -9.2530E+00 -3.8837E+01  7.6564E+01  2.5039E+02  5.2788E+01 -2.5412E+01 -5.6693E+00 -4.7965E+00 -7.1581E+01
            -8.0529E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1860.11479439199        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9838E-01  9.0252E-01  1.0153E+00  1.0570E+00  8.2873E-01  9.2219E-01  9.7795E-01  8.4268E-01  1.0834E+00  8.8267E-01
             1.6146E+00
 PARAMETER:  9.8382E-02 -2.5666E-03  1.1522E-01  1.5547E-01 -8.7864E-02  1.8996E-02  7.7699E-02 -7.1169E-02  1.8006E-01 -2.4807E-02
             5.7908E-01
 GRADIENT:   1.5566E+02  2.4913E+00  8.2076E+00  1.7902E+01  2.5557E+01 -4.4054E+00  9.0841E-01  6.4341E+00  1.8027E+01 -1.4645E+01
            -7.6314E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1870.06131342298        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      202
 NPARAMETR:  9.9862E-01  6.3390E-01  8.0123E-01  1.2123E+00  6.4254E-01  9.4082E-01  1.1910E+00  3.9194E-01  8.6952E-01  8.8781E-01
             1.5718E+00
 PARAMETER:  9.8617E-02 -3.5586E-01 -1.2161E-01  2.9249E-01 -3.4232E-01  3.8999E-02  2.7479E-01 -8.3666E-01 -3.9811E-02 -1.8992E-02
             5.5223E-01
 GRADIENT:  -9.5214E+00  8.8147E+00  6.0906E+00 -9.7521E+00 -2.3484E+00 -1.4166E+01 -9.0868E+00  2.4212E+00 -2.1199E+01  4.1541E+00
            -8.3139E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1877.89021712572        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  1.0048E+00  5.3012E-01  7.3388E-01  1.2781E+00  5.8129E-01  9.8250E-01  1.5494E+00  1.1680E-01  8.9929E-01  8.4385E-01
             1.7256E+00
 PARAMETER:  1.0475E-01 -5.3464E-01 -2.0941E-01  3.4535E-01 -4.4250E-01  8.2348E-02  5.3786E-01 -2.0473E+00 -6.1461E-03 -6.9777E-02
             6.4555E-01
 GRADIENT:   3.9367E+00  1.1945E+01 -3.5797E+00  2.5655E+01 -3.7861E+00  3.8093E+00 -1.0190E+00  2.4143E-01  2.9496E+00  7.3701E+00
             9.3182E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1880.18698295671        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      556
 NPARAMETR:  9.9943E-01  2.9468E-01  7.3532E-01  1.3782E+00  5.2602E-01  9.7052E-01  2.2310E+00  2.7278E-02  8.4072E-01  7.9732E-01
             1.7129E+00
 PARAMETER:  9.9430E-02 -1.1219E+00 -2.0745E-01  4.2078E-01 -5.4242E-01  7.0078E-02  9.0243E-01 -3.5017E+00 -7.3499E-02 -1.2649E-01
             6.3816E-01
 GRADIENT:   4.0857E+00  3.2390E+00  9.0121E+00 -1.1854E+01 -1.9461E+01  8.5094E-01  1.2172E+00  6.1953E-03 -2.2329E+00 -4.7223E+00
            -6.3235E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1881.87954486091        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  9.9366E-01  1.1040E-01  8.7249E-01  1.5137E+00  5.5934E-01  9.6827E-01  3.1327E+00  1.0000E-02  8.1696E-01  8.9204E-01
             1.7008E+00
 PARAMETER:  9.3640E-02 -2.1037E+00 -3.6403E-02  5.1457E-01 -4.8100E-01  6.7759E-02  1.2419E+00 -5.1291E+00 -1.0217E-01 -1.4240E-02
             6.3113E-01
 GRADIENT:   3.4525E+00  9.0577E-01 -6.0939E-01  9.2588E+00  1.8129E+00  9.3928E-01 -6.9662E-01  0.0000E+00  5.2173E-01  5.6889E-02
            -6.7999E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1882.32735690779        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      910
 NPARAMETR:  9.8920E-01  3.6780E-02  8.7588E-01  1.5504E+00  5.4598E-01  9.6286E-01  5.9485E+00  1.0000E-02  8.0630E-01  8.9211E-01
             1.7194E+00
 PARAMETER:  8.9143E-02 -3.2028E+00 -3.2527E-02  5.3852E-01 -5.0516E-01  6.2149E-02  1.8831E+00 -8.4006E+00 -1.1530E-01 -1.4165E-02
             6.4198E-01
 GRADIENT:  -2.6196E+00  5.9386E-01  2.0047E+00  4.6825E+00 -3.7865E+00 -6.1498E-01  4.0653E-01  0.0000E+00  1.1441E+00  5.1132E-01
             2.6250E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1882.44185624586        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1088
 NPARAMETR:  9.8849E-01  1.0050E-02  8.8155E-01  1.5647E+00  5.4371E-01  9.6240E-01  1.1324E+01  1.0000E-02  7.9943E-01  8.9291E-01
             1.7193E+00
 PARAMETER:  8.8427E-02 -4.5002E+00 -2.6073E-02  5.4770E-01 -5.0934E-01  6.1677E-02  2.5269E+00 -1.2138E+01 -1.2385E-01 -1.3267E-02
             6.4195E-01
 GRADIENT:  -2.5250E+00  4.4278E-01  2.0295E+00  3.5785E+00 -3.7908E+00 -6.2537E-01  1.1277E-01  0.0000E+00  6.7980E-02  2.4663E-01
             2.2930E+00

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1882.46646169467        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1222
 NPARAMETR:  9.8984E-01  1.0000E-02  8.7968E-01  1.5612E+00  5.4429E-01  9.6446E-01  1.1181E+01  1.0000E-02  7.9936E-01  8.9207E-01
             1.7147E+00
 PARAMETER:  8.9787E-02 -4.5233E+00 -2.8196E-02  5.4547E-01 -5.0827E-01  6.3810E-02  2.5143E+00 -1.2138E+01 -1.2394E-01 -1.4208E-02
             6.3924E-01
 GRADIENT:   9.6176E-01  1.4144E-01 -4.4132E-01 -3.0350E+00  1.1902E+00  1.9617E-01  5.4738E-03  0.0000E+00  3.2037E-03 -1.1723E-01
             6.9864E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1222
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.1286E-04  8.6387E-04 -5.7470E-05 -6.8288E-03 -1.4339E-02
 SE:             2.9589E-02  1.9772E-03  1.9226E-04  2.8793E-02  2.4462E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8887E-01  6.6217E-01  7.6501E-01  8.1253E-01  5.5776E-01

 ETASHRINKSD(%)  8.7225E-01  9.3376E+01  9.9356E+01  3.5399E+00  1.8049E+01
 ETASHRINKVR(%)  1.7369E+00  9.9561E+01  9.9996E+01  6.9545E+00  3.2841E+01
 EBVSHRINKSD(%)  1.0067E+00  9.4038E+01  9.9347E+01  3.4890E+00  1.7167E+01
 EBVSHRINKVR(%)  2.0032E+00  9.9645E+01  9.9996E+01  6.8564E+00  3.1387E+01
 RELATIVEINF(%)  9.2256E+01  1.7736E-02  3.8859E-04  6.1423E+00  6.9797E+00
 EPSSHRINKSD(%)  2.9683E+01
 EPSSHRINKVR(%)  5.0555E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1882.4664616946723     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -963.52792848999957     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.20
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1882.466       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.90E-01  1.00E-02  8.80E-01  1.56E+00  5.44E-01  9.64E-01  1.12E+01  1.00E-02  7.99E-01  8.92E-01  1.71E+00
 


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
+       -1.16E+01  2.09E+05
 
 TH 3
+       -1.13E+01  3.10E+01  6.07E+02
 
 TH 4
+       -1.34E+01  5.85E+01 -1.15E+02  6.66E+02
 
 TH 5
+        2.05E+01 -1.16E+02 -1.28E+03 -8.24E+01  3.26E+03
 
 TH 6
+        1.58E+00 -1.59E+00  1.37E+00 -4.71E+00 -2.57E+00  2.04E+02
 
 TH 7
+       -2.82E-01 -1.84E+01  5.08E+00  6.02E+00 -1.33E+01  3.30E-02  4.64E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.35E+00 -4.90E+01  2.41E+01 -6.85E+00 -3.90E+00 -1.12E-01  2.43E+00  0.00E+00  2.72E+02
 
 TH10
+        5.40E-01 -3.10E+01 -2.00E+01  4.46E+00 -7.37E+01  5.92E-01  4.11E+00  0.00E+00  3.38E-01  1.21E+02
 
 TH11
+       -1.33E+01 -1.00E+01 -1.26E+01 -1.09E+01  1.32E+01  2.58E+00  9.00E-01  0.00E+00  7.72E+00  2.35E+01  1.49E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.582
Stop Time:
Wed Sep 29 21:48:10 CDT 2021
