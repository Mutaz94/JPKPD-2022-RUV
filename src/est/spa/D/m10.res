Sat Sep 25 14:02:45 CDT 2021
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
$DATA ../../../../data/spa/D/dat10.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m10.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1161.87373039605        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.3040E+01 -1.1739E+02 -4.4561E+01 -1.3088E+02  8.0496E+01 -3.5886E+02 -2.7945E+02 -1.1422E+01 -3.0189E+02 -3.4180E+01
            -5.9435E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1437.77346970572        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0421E+00  1.3730E+00  1.1710E+00  8.4202E-01  1.2927E+00  1.4408E+00  2.8860E+00  1.0543E+00  1.3835E+00  8.9890E-01
             1.0300E+00
 PARAMETER:  1.4126E-01  4.1703E-01  2.5788E-01 -7.1947E-02  3.5673E-01  4.6519E-01  1.1599E+00  1.5285E-01  4.2464E-01 -6.5841E-03
             1.2957E-01
 GRADIENT:   5.4408E+01  3.2806E+01 -5.1463E+00 -2.3940E+00  1.1340E+02 -2.8171E+01  8.8108E+01 -9.1689E+00 -3.7072E+00 -3.1435E+01
            -2.4828E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1450.65091442105        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      244
 NPARAMETR:  1.0533E+00  1.0035E+00  1.7094E+00  1.0718E+00  1.2967E+00  1.6163E+00  3.5246E+00  2.0546E+00  1.1788E+00  1.0572E+00
             1.0909E+00
 PARAMETER:  1.5191E-01  1.0353E-01  6.3612E-01  1.6937E-01  3.5980E-01  5.8014E-01  1.3598E+00  8.2010E-01  2.6447E-01  1.5565E-01
             1.8703E-01
 GRADIENT:   8.9063E+00  5.7865E+00 -1.7948E+01 -1.0042E+01  5.7633E+01 -2.2771E+01  3.7287E+01  6.5355E+00  6.3554E+00 -6.3293E+00
             1.1879E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1454.06349261120        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      426
 NPARAMETR:  1.0543E+00  9.6396E-01  1.9337E+00  1.1046E+00  1.3173E+00  1.6729E+00  3.3404E+00  1.9942E+00  1.2216E+00  1.1531E+00
             1.0372E+00
 PARAMETER:  1.5284E-01  6.3294E-02  7.5946E-01  1.9948E-01  3.7556E-01  6.1454E-01  1.3061E+00  7.9026E-01  3.0016E-01  2.4250E-01
             1.3656E-01
 GRADIENT:   9.6618E+00  4.9799E+00 -3.3477E+00 -7.8588E+00  4.1202E+01 -8.0144E+00  2.1507E+01 -3.1822E+00  3.0611E+00  2.0127E+00
            -6.0202E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1455.20081227847        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:      575
 NPARAMETR:  1.0541E+00  9.6270E-01  1.9313E+00  1.1283E+00  1.2105E+00  1.6551E+00  3.3432E+00  1.9924E+00  1.2213E+00  1.0910E+00
             1.0630E+00
 PARAMETER:  1.5268E-01  6.1989E-02  7.5820E-01  2.2071E-01  2.9104E-01  6.0386E-01  1.3069E+00  7.8933E-01  2.9988E-01  1.8709E-01
             1.6107E-01
 GRADIENT:   9.8688E+00  1.3729E+01  1.0029E+01 -1.0189E+01 -2.5631E+01 -1.2270E+01  1.6609E+01  7.5404E-01  1.0849E+00  8.8716E+00
             4.4844E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1455.91218703841        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      750
 NPARAMETR:  1.0541E+00  9.6270E-01  1.9313E+00  1.1395E+00  1.2299E+00  1.7057E+00  3.3432E+00  1.9924E+00  1.2213E+00  1.0322E+00
             1.0569E+00
 PARAMETER:  1.5268E-01  6.1989E-02  7.5821E-01  2.3061E-01  3.0693E-01  6.3398E-01  1.3069E+00  7.8932E-01  2.9988E-01  1.3168E-01
             1.5537E-01
 GRADIENT:   9.2328E+00  1.4542E+01  5.6036E+00 -1.0064E-01  9.4542E-03  1.2437E-02  1.6178E+01  3.6508E-01  4.5937E-01  6.5807E-04
             1.3548E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1456.02388069711        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:      920
 NPARAMETR:  1.0535E+00  9.5963E-01  1.9292E+00  1.1396E+00  1.2298E+00  1.7048E+00  3.3317E+00  1.9921E+00  1.2210E+00  1.0322E+00
             1.0566E+00
 PARAMETER:  1.5210E-01  5.8793E-02  7.5710E-01  2.3066E-01  3.0688E-01  6.3348E-01  1.3035E+00  7.8919E-01  2.9969E-01  1.3170E-01
             1.5506E-01
 GRADIENT:   8.8033E+00  1.3952E+01  5.4037E+00 -7.2955E-02  3.1587E-01 -1.7512E-01  1.5409E+01  4.8319E-01  3.7915E-01  1.7263E-02
             1.3667E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1456.02388069711        NO. OF FUNC. EVALS.:  33
 CUMULATIVE NO. OF FUNC. EVALS.:      953
 NPARAMETR:  1.0535E+00  9.5965E-01  1.9289E+00  1.1396E+00  1.2296E+00  1.7059E+00  3.3326E+00  1.9924E+00  1.2210E+00  1.0321E+00
             1.0566E+00
 PARAMETER:  1.5210E-01  5.8793E-02  7.5710E-01  2.3066E-01  3.0688E-01  6.3348E-01  1.3035E+00  7.8919E-01  2.9969E-01  1.3170E-01
             1.5506E-01
 GRADIENT:   1.0184E+05 -3.0978E+05  2.0471E+04 -6.2377E-02  2.3797E-01 -2.9292E-01 -2.3752E+04 -3.9305E+04  5.1684E+04  1.4926E-02
             1.2828E-02
 NUMSIGDIG:         3.3         3.3         3.3         3.4         2.8         2.6         3.3         3.3         3.3         2.7
                    3.7

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      953
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.1449E-03  2.4412E-03 -5.9736E-02 -3.5315E-02 -3.9327E-02
 SE:             2.9965E-02  2.2667E-02  1.6441E-02  1.9316E-02  1.8829E-02
 N:                     100         100         100         100         100

 P VAL.:         8.3752E-01  9.1424E-01  2.7995E-04  6.7505E-02  3.6737E-02

 ETASHRINKSD(%)  1.0000E-10  2.4063E+01  4.4919E+01  3.5289E+01  3.6921E+01
 ETASHRINKVR(%)  1.0000E-10  4.2336E+01  6.9661E+01  5.8125E+01  6.0210E+01
 EBVSHRINKSD(%)  1.7702E-01  1.7382E+01  5.2629E+01  3.6373E+01  3.1726E+01
 EBVSHRINKVR(%)  3.5372E-01  3.1743E+01  7.7560E+01  5.9516E+01  5.3387E+01
 RELATIVEINF(%)  9.9542E+01  1.5667E+01  6.9404E+00  8.1718E+00  1.9917E+01
 EPSSHRINKSD(%)  4.5453E+01
 EPSSHRINKVR(%)  7.0246E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1456.0238806971117     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -720.87305413337356     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1456.024       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  9.60E-01  1.93E+00  1.14E+00  1.23E+00  1.70E+00  3.33E+00  1.99E+00  1.22E+00  1.03E+00  1.06E+00
 


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
+        3.02E+08
 
 TH 2
+       -5.04E+08  8.41E+08
 
 TH 3
+        3.31E+07 -5.53E+07  3.63E+06
 
 TH 4
+        1.84E+08 -3.07E+08  2.02E+07  1.12E+08
 
 TH 5
+        4.53E+03 -2.14E+08  4.50E+02  7.81E+07  3.42E+02
 
 TH 6
+        9.05E+02 -1.52E+03  1.00E+02 -2.72E-01 -2.91E-01  6.83E+01
 
 TH 7
+       -1.78E+01  3.48E+01 -7.00E+02  2.64E+02 -1.62E+02 -3.37E+01  4.11E+05
 
 TH 8
+       -5.00E+01  7.77E+01 -1.95E+03  7.83E+02 -4.70E+02 -9.27E+01  2.95E+03  3.14E+06
 
 TH 9
+        2.15E+02 -3.54E+02  8.29E+03 -3.33E+03  1.99E+03  3.98E+02  2.42E+02 -3.50E+04  5.78E+07
 
 TH10
+        2.40E+03 -4.00E+03  2.63E+02 -2.17E+08 -6.26E+01 -2.53E-01 -8.92E+01 -2.43E+02  1.05E+03  5.46E+01
 
 TH11
+       -1.99E+04 -4.93E+08 -2.19E+03  1.80E+08 -1.02E+01  1.17E+00  7.35E+02  2.03E+03 -8.71E+03  1.67E+01  1.88E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.091
Stop Time:
Sat Sep 25 14:03:13 CDT 2021
