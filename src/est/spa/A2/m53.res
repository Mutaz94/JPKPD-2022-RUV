Wed Sep 29 12:53:05 CDT 2021
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
$DATA ../../../../data/spa/A2/dat53.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m53.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1065.85233309200        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5328E+02 -8.5641E+00  2.3696E+01 -8.8039E+00  1.0382E+02  3.2852E+01  1.0858E+00 -1.4892E+01  1.7376E+00 -2.5443E+01
            -1.0691E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1391.79509473938        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1976E+00  1.0528E+00  9.8541E-01  1.0769E+00  9.5140E-01  1.1663E+00  8.0878E-01  9.4699E-01  7.4830E-01  7.0861E-01
             2.7443E+00
 PARAMETER:  2.8032E-01  1.5145E-01  8.5298E-02  1.7408E-01  5.0177E-02  2.5383E-01 -1.1223E-01  4.5536E-02 -1.8995E-01 -2.4445E-01
             1.1095E+00
 GRADIENT:   4.4658E+02  2.8170E+01  7.0787E+00  3.4633E+01 -1.2977E+01  4.5375E+00  5.1815E-01  2.0614E+00 -8.1134E+00  1.0641E+01
             4.2456E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1407.41860866304        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.1631E+00  9.0073E-01  4.1518E-01  1.0967E+00  5.3974E-01  1.1948E+00  9.1954E-01  7.9421E-01  9.3488E-01  2.2779E-01
             2.3907E+00
 PARAMETER:  2.5107E-01 -4.5520E-03 -7.7903E-01  1.9227E-01 -5.1667E-01  2.7799E-01  1.6116E-02 -1.3040E-01  3.2666E-02 -1.3793E+00
             9.7159E-01
 GRADIENT:   4.0024E+02  4.0518E+01  3.9579E+00  9.7884E+01  7.7280E+00  2.8690E+01  2.1141E+00  6.5648E+00  2.5954E+01  1.9294E+00
             1.7289E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1411.11147430051        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.1455E+00  1.0210E+00  5.0170E-01  1.0248E+00  6.4487E-01  1.1716E+00  7.9194E-01  5.9666E-01  9.9794E-01  2.8425E-01
             2.4323E+00
 PARAMETER:  2.3587E-01  1.2082E-01 -5.8974E-01  1.2446E-01 -3.3871E-01  2.5840E-01 -1.3327E-01 -4.1641E-01  9.7934E-02 -1.1579E+00
             9.8882E-01
 GRADIENT:   2.2953E+02  2.9909E+01  2.2454E+01  1.8599E+01 -3.9663E+01  9.0405E+00 -5.0245E+00 -6.9557E-01  2.3248E+01  1.8085E-01
            -1.8343E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1440.60554166951        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      473
 NPARAMETR:  9.7265E-01  7.8652E-01  2.4960E-01  1.0197E+00  4.0025E-01  1.0210E+00  9.9954E-01  3.7549E-01  7.8922E-01  1.7456E-01
             2.2265E+00
 PARAMETER:  7.2270E-02 -1.4013E-01 -1.2879E+00  1.1951E-01 -8.1568E-01  1.2079E-01  9.9536E-02 -8.7953E-01 -1.3672E-01 -1.6455E+00
             9.0042E-01
 GRADIENT:  -2.2186E+01  1.0378E+01 -1.6066E+01  5.2130E+01  2.2026E+01 -1.0305E+01  3.4445E+00 -5.3345E-01 -4.2239E+00  7.4235E-01
             6.7753E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1442.70722675081        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  9.8052E-01  6.9500E-01  2.2518E-01  1.0128E+00  3.5388E-01  1.0533E+00  1.0338E+00  5.8390E-01  7.9846E-01  1.1842E-01
             2.1294E+00
 PARAMETER:  8.0323E-02 -2.6384E-01 -1.3908E+00  1.1270E-01 -9.3880E-01  1.5193E-01  1.3322E-01 -4.3803E-01 -1.2508E-01 -2.0335E+00
             8.5584E-01
 GRADIENT:  -7.0983E-01 -2.8402E-01 -3.3202E-01 -3.3947E-01  2.0353E-01  2.0036E-01  1.3704E-02  7.7640E-02 -1.9782E-01  1.8122E-01
             8.3769E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1442.76292061042        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      825
 NPARAMETR:  9.8082E-01  7.0264E-01  2.2679E-01  1.0101E+00  3.5758E-01  1.0526E+00  1.0432E+00  5.9975E-01  7.9934E-01  4.4852E-02
             2.1319E+00
 PARAMETER:  8.0638E-02 -2.5291E-01 -1.3838E+00  1.1001E-01 -9.2839E-01  1.5125E-01  1.4234E-01 -4.1124E-01 -1.2397E-01 -3.0044E+00
             8.5699E-01
 GRADIENT:   9.4078E-02 -1.1740E+00  2.0360E-01 -2.2689E+00  1.0230E+00  1.4041E-01  1.2203E+00  1.9954E-03  6.6530E-01  2.2075E-02
             6.7952E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1442.78348649390        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1002
 NPARAMETR:  9.8017E-01  7.1928E-01  2.2497E-01  1.0027E+00  3.6151E-01  1.0517E+00  1.0162E+00  6.1338E-01  7.9903E-01  1.1663E-02
             2.1338E+00
 PARAMETER:  7.9972E-02 -2.2951E-01 -1.3918E+00  1.0267E-01 -9.1746E-01  1.5045E-01  1.1609E-01 -3.8877E-01 -1.2436E-01 -4.3513E+00
             8.5788E-01
 GRADIENT:   1.1426E-03 -1.8202E-01 -1.1058E-01  2.5445E-02  2.8466E-01  2.0598E-02 -1.4744E-01 -1.5776E-02 -2.1270E-01  1.9414E-03
             1.5707E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1442.78450272033        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  9.8019E-01  7.1829E-01  2.2496E-01  1.0030E+00  3.6115E-01  1.0517E+00  1.0182E+00  6.1454E-01  7.9959E-01  1.0000E-02
             2.1326E+00
 PARAMETER:  7.9992E-02 -2.3089E-01 -1.3918E+00  1.0302E-01 -9.1845E-01  1.5044E-01  1.1805E-01 -3.8688E-01 -1.2365E-01 -4.9267E+00
             8.5734E-01
 GRADIENT:   1.3350E-02 -1.5700E-02 -9.6849E-03 -4.9339E-02  4.3694E-02  2.3977E-03  2.1065E-02  1.1590E-02 -2.4363E-02  0.0000E+00
             6.0843E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1059
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.4407E-04  6.1081E-03 -1.1799E-02 -1.0996E-02  1.8718E-04
 SE:             2.9514E-02  2.3621E-02  1.1103E-02  2.4893E-02  3.7015E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8529E-01  7.9596E-01  2.8792E-01  6.5868E-01  6.1308E-01

 ETASHRINKSD(%)  1.1253E+00  2.0865E+01  6.2804E+01  1.6606E+01  9.8760E+01
 ETASHRINKVR(%)  2.2379E+00  3.7377E+01  8.6165E+01  3.0454E+01  9.9985E+01
 EBVSHRINKSD(%)  1.3909E+00  2.0304E+01  6.3302E+01  1.6549E+01  9.8767E+01
 EBVSHRINKVR(%)  2.7625E+00  3.6486E+01  8.6532E+01  3.0359E+01  9.9985E+01
 RELATIVEINF(%)  9.4642E+01  3.7729E+00  1.4789E+00  1.6417E+01  6.1114E-04
 EPSSHRINKSD(%)  3.5911E+01
 EPSSHRINKVR(%)  5.8926E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1442.7845027203314     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -707.63367615659320     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1442.785       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  7.18E-01  2.25E-01  1.00E+00  3.61E-01  1.05E+00  1.02E+00  6.15E-01  8.00E-01  1.00E-02  2.13E+00
 


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
+        1.02E+03
 
 TH 2
+       -2.86E+01  1.33E+03
 
 TH 3
+       -1.59E+02  2.06E+03  8.23E+03
 
 TH 4
+       -4.03E+01  2.34E+02 -1.60E+03  1.20E+03
 
 TH 5
+        1.27E+02 -3.71E+03 -8.96E+03  7.39E+02  1.31E+04
 
 TH 6
+        4.61E-01 -5.01E+00  6.72E+00 -1.16E+01  2.88E+01  1.70E+02
 
 TH 7
+       -4.46E-01  4.14E+01 -1.39E+02 -7.11E+00  8.18E+01  8.67E-01  7.61E+01
 
 TH 8
+        1.08E+00 -1.90E+01 -4.83E+01  4.32E+00  7.40E+01  9.84E-01  8.03E+00  1.60E+01
 
 TH 9
+        4.23E+00 -2.15E+01  6.21E+01 -6.53E+00  4.92E+01  1.38E+00  1.51E+01  8.40E+00  1.59E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.06E+01 -1.20E+01 -4.98E+01 -1.40E+01  1.55E+00  1.36E+00  1.19E+01  8.86E+00  1.06E+01  0.00E+00  5.89E+01
 
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
 #CPUT: Total CPU Time in Seconds,       18.150
Stop Time:
Wed Sep 29 12:53:25 CDT 2021
