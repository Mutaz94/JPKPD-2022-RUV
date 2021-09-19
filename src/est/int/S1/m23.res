Sat Sep 18 02:26:10 CDT 2021
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
$DATA ../../../../data/int/S1/dat23.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3654.04526815202        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1574E+01  6.2287E+01 -1.4367E+01  8.5272E+01  9.1465E+00  1.8061E+01 -3.1283E+01 -2.5787E+01 -3.8793E+00 -2.3759E+00
            -3.7325E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3686.97957057311        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  1.0178E+00  9.8878E-01  1.0584E+00  9.7257E-01  1.0174E+00  9.2202E-01  1.1609E+00  1.1196E+00  9.5542E-01  9.9334E-01
             1.1486E+00
 PARAMETER:  1.1762E-01  8.8719E-02  1.5677E-01  7.2186E-02  1.1727E-01  1.8810E-02  2.4917E-01  2.1301E-01  5.4393E-02  9.3317E-02
             2.3855E-01
 GRADIENT:   4.3659E+01  1.6111E+01 -1.0460E+00  2.0959E+01  5.2523E-01 -1.4735E+01  1.6546E+00 -1.1471E+01 -1.6721E+00  1.6982E+00
            -2.5459E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3689.57121768221        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      167
 NPARAMETR:  1.0185E+00  9.7493E-01  1.1842E+00  9.8303E-01  1.0561E+00  9.8249E-01  1.1919E+00  1.5588E+00  9.1112E-01  9.8721E-01
             1.1349E+00
 PARAMETER:  1.1837E-01  7.4608E-02  2.6904E-01  8.2880E-02  1.5462E-01  8.2335E-02  2.7557E-01  5.4389E-01  6.9201E-03  8.7125E-02
             2.2658E-01
 GRADIENT:   4.6920E+01 -7.3653E-01 -2.7158E+00  2.3524E+01  7.0525E+00  1.1095E+01  2.6535E+00 -3.4820E+00  6.0048E-01  1.7818E+00
            -2.6226E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3689.76830909263        NO. OF FUNC. EVALS.:  82
 CUMULATIVE NO. OF FUNC. EVALS.:      249
 NPARAMETR:  1.0170E+00  9.8358E-01  1.1989E+00  9.7780E-01  1.0662E+00  9.8061E-01  1.1916E+00  1.5850E+00  9.2492E-01  1.0007E+00
             1.1359E+00
 PARAMETER:  1.1688E-01  8.3445E-02  2.8141E-01  7.7553E-02  1.6406E-01  8.0422E-02  2.7527E-01  5.6060E-01  2.1957E-02  1.0074E-01
             2.2744E-01
 GRADIENT:   4.2877E+01  3.2806E-01 -1.1958E+00  2.0125E+01  7.7710E+00  1.0376E+01  3.6419E+00 -2.8259E+00  4.6094E+00  2.9400E+00
            -2.3394E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3689.81429633952        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:      444             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0170E+00  9.8394E-01  1.1995E+00  9.7740E-01  1.0661E+00  9.7923E-01  1.1907E+00  1.5850E+00  9.2305E-01  9.9950E-01
             1.1368E+00
 PARAMETER:  1.1684E-01  8.3808E-02  2.8188E-01  7.7140E-02  1.6397E-01  7.9011E-02  2.7453E-01  5.6059E-01  1.9927E-02  9.9501E-02
             2.2821E-01
 GRADIENT:   4.2693E+01  5.4039E-01 -9.9759E-01  1.9242E+01  7.1917E+00  9.8352E+00  3.5273E+00 -2.8643E+00  4.1877E+00  2.6607E+00
            -2.2010E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3689.84889056012        NO. OF FUNC. EVALS.: 152
 CUMULATIVE NO. OF FUNC. EVALS.:      596
 NPARAMETR:  1.0170E+00  9.8389E-01  1.1996E+00  9.7745E-01  1.0662E+00  9.6096E-01  1.1908E+00  1.5846E+00  9.2302E-01  9.9949E-01
             1.1367E+00
 PARAMETER:  1.1690E-01  8.3758E-02  2.8202E-01  7.7190E-02  1.6406E-01  6.0181E-02  2.7466E-01  5.6030E-01  1.9896E-02  9.9492E-02
             2.2809E-01
 GRADIENT:  -8.5719E-01 -5.4449E+00 -2.4906E+00  1.0939E+01  1.1535E+00 -1.6348E+00  1.6196E+00 -3.5757E+00  3.6035E+00  2.2945E+00
            -2.2748E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3689.86321430402        NO. OF FUNC. EVALS.: 205
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  1.0170E+00  9.8400E-01  1.1998E+00  9.7733E-01  1.0661E+00  9.6490E-01  1.1906E+00  1.5857E+00  9.2267E-01  9.9924E-01
             1.1369E+00
 PARAMETER:  1.1690E-01  8.3875E-02  2.8212E-01  7.7072E-02  1.6404E-01  6.4264E-02  2.7445E-01  5.6103E-01  1.9519E-02  9.9237E-02
             2.2827E-01
 GRADIENT:  -8.2758E-01 -5.3773E+00 -2.5141E+00  1.0676E+01  9.6634E-01 -1.8458E-02  1.5908E+00 -3.5224E+00  3.5478E+00  2.2386E+00
            -2.2356E+01

0ITERATION NO.:   33    OBJECTIVE VALUE:  -3689.86356386988        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      921
 NPARAMETR:  1.0170E+00  9.8400E-01  1.1998E+00  9.7733E-01  1.0661E+00  9.6490E-01  1.1906E+00  1.5857E+00  9.2258E-01  9.9924E-01
             1.1369E+00
 PARAMETER:  1.1690E-01  8.3875E-02  2.8212E-01  7.7072E-02  1.6404E-01  6.4264E-02  2.7445E-01  5.6103E-01  1.9419E-02  9.9237E-02
             2.2827E-01
 GRADIENT:  -3.0228E+05  7.0676E+05  3.7300E+01 -7.0677E+05 -1.3218E+01 -1.9491E-02 -2.5752E+05  1.2586E+05 -7.0677E+05 -7.0677E+05
            -3.7827E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      921
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.8703E-04 -2.0079E-02 -3.0726E-02  9.9945E-03 -3.8476E-02
 SE:             2.9890E-02  2.3543E-02  2.1323E-02  2.6174E-02  2.2545E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7899E-01  3.9373E-01  1.4959E-01  7.0258E-01  8.7882E-02

 ETASHRINKSD(%)  1.0000E-10  2.1128E+01  2.8565E+01  1.2312E+01  2.4473E+01
 ETASHRINKVR(%)  1.0000E-10  3.7792E+01  4.8971E+01  2.3109E+01  4.2956E+01
 EBVSHRINKSD(%)  3.2268E-01  2.0466E+01  3.2250E+01  1.3027E+01  2.1949E+01
 EBVSHRINKVR(%)  6.4433E-01  3.6743E+01  5.4100E+01  2.4357E+01  3.9081E+01
 RELATIVEINF(%)  9.9352E+01  3.6358E+01  3.4618E+01  5.0915E+01  3.1139E+01
 EPSSHRINKSD(%)  2.0804E+01
 EPSSHRINKVR(%)  3.7280E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3689.8635638698802     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2035.7742041014694     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.62
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3689.864       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  9.84E-01  1.20E+00  9.77E-01  1.07E+00  9.65E-01  1.19E+00  1.59E+00  9.23E-01  9.99E-01  1.14E+00
 


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
+        1.25E+09
 
 TH 2
+        2.03E+03  1.82E+09
 
 TH 3
+        1.40E+00  2.49E+01  3.08E+08
 
 TH 4
+        1.52E+09 -7.50E+04  5.34E+08  3.70E+09
 
 TH 5
+        8.50E+08 -1.03E+09  2.99E+08  1.03E+09  1.16E+09
 
 TH 6
+       -2.51E+03  3.03E+03 -8.28E-01 -3.06E+03 -6.24E-01  2.18E+02
 
 TH 7
+        4.55E+08 -1.48E+04  1.60E+08  5.53E+08  3.09E+08 -9.13E+02  3.31E+08
 
 TH 8
+        2.25E+02 -3.56E+03  5.87E+07 -8.34E+03 -3.75E+03  3.36E+02 -1.64E+03  2.23E+07
 
 TH 9
+        3.95E+04 -4.78E+04  1.13E+01 -1.96E+09  1.09E+01 -3.24E+03 -5.86E+08 -5.26E+03  2.08E+09
 
 TH10
+        2.51E+04 -3.03E+04  1.12E+00  1.81E+09 -3.93E+01 -2.98E+03  9.13E+03 -3.35E+03  4.71E+04  1.77E+09
 
 TH11
+       -8.79E+00 -3.66E+01 -2.01E+08 -5.49E+00 -3.89E+08  5.57E-01  2.17E+01 -7.66E+07  2.02E+01  8.92E+00  5.25E+08
 
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
 #CPUT: Total CPU Time in Seconds,       43.055
Stop Time:
Sat Sep 18 02:26:55 CDT 2021
