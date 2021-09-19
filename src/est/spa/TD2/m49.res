Sat Sep 18 14:42:47 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat49.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1587.36841288052        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4995E+02 -4.8395E+01 -3.4504E+01 -7.2144E+00  4.3856E+01 -6.8546E+00 -4.1915E+01 -1.9408E+00 -2.3540E+01  9.2240E-01
            -2.4231E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1602.39674239152        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.3967E-01  1.1382E+00  1.0765E+00  9.1989E-01  1.0556E+00  1.0054E+00  1.3842E+00  1.0430E+00  1.0852E+00  9.4471E-01
             1.0496E+00
 PARAMETER:  3.7772E-02  2.2942E-01  1.7372E-01  1.6499E-02  1.5409E-01  1.0542E-01  4.2512E-01  1.4206E-01  1.8177E-01  4.3122E-02
             1.4837E-01
 GRADIENT:   1.0093E+01  5.6266E+00  1.0532E+01 -5.9372E+00 -6.3712E+00  7.6204E-01  4.8236E+00 -5.5215E+00  8.6631E+00  8.4461E-01
            -1.9272E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1602.83412640496        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  9.4000E-01  1.1380E+00  1.0600E+00  9.4130E-01  1.0352E+00  1.0215E+00  1.4340E+00  1.1971E+00  1.0376E+00  8.7960E-01
             1.0724E+00
 PARAMETER:  3.8121E-02  2.2923E-01  1.5831E-01  3.9511E-02  1.3456E-01  1.2129E-01  4.6044E-01  2.7992E-01  1.3692E-01 -2.8287E-02
             1.6989E-01
 GRADIENT:  -2.7682E+01  1.2729E+01  2.3226E+00  9.7248E+00 -8.2712E+00  5.0588E-01  4.9625E+00  4.2862E-01  4.7749E+00 -2.0890E+00
             6.8388E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1603.32640743184        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      432
 NPARAMETR:  9.4758E-01  1.0943E+00  1.1440E+00  9.6526E-01  1.0560E+00  1.0192E+00  1.4401E+00  1.2707E+00  1.0062E+00  9.2391E-01
             1.0604E+00
 PARAMETER:  4.6152E-02  1.9013E-01  2.3449E-01  6.4645E-02  1.5451E-01  1.1904E-01  4.6475E-01  3.3957E-01  1.0618E-01  2.0861E-02
             1.5862E-01
 GRADIENT:  -9.5785E+00  4.9525E+00  1.1436E+00  3.9583E+00 -3.7564E+00  1.3099E-01  1.8510E+00  1.3544E-01  1.7268E+00 -4.5162E-01
             2.2886E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1603.38573214837        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      617             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5016E-01  1.0927E+00  1.1546E+00  9.6410E-01  1.0622E+00  1.0188E+00  1.4311E+00  1.2818E+00  9.9462E-01  9.3566E-01
             1.0570E+00
 PARAMETER:  4.8877E-02  1.8864E-01  2.4371E-01  6.3437E-02  1.6038E-01  1.1858E-01  4.5842E-01  3.4826E-01  9.4601E-02  3.3496E-02
             1.5544E-01
 GRADIENT:   3.4443E+01  9.4251E+00  1.1492E+00  8.0504E+00 -1.2150E+00  6.9817E+00  3.7800E+00  3.1858E-02  6.8874E-01  5.0115E-02
             9.9992E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1603.38684126704        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  9.5026E-01  1.0927E+00  1.1544E+00  9.6402E-01  1.0624E+00  1.0187E+00  1.4308E+00  1.2818E+00  9.9536E-01  9.3557E-01
             1.0569E+00
 PARAMETER:  4.8985E-02  1.8863E-01  2.4354E-01  6.3362E-02  1.6048E-01  1.1855E-01  4.5820E-01  3.4830E-01  9.5346E-02  3.3396E-02
             1.5531E-01
 GRADIENT:  -3.4625E+00  2.2209E+00  7.5597E-01  1.4030E+00 -1.9191E+00 -5.6784E-04  5.0253E-01 -2.6192E-02 -6.5087E-03 -2.3679E-02
             8.2671E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1603.39116343727        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      984
 NPARAMETR:  9.5181E-01  1.0927E+00  1.1544E+00  9.6405E-01  1.0623E+00  1.0187E+00  1.4253E+00  1.2817E+00  9.9696E-01  9.3638E-01
             1.0548E+00
 PARAMETER:  5.0608E-02  1.8868E-01  2.4360E-01  6.3387E-02  1.6044E-01  1.1855E-01  4.5441E-01  3.4821E-01  9.6958E-02  3.4264E-02
             1.5338E-01
 GRADIENT:  -1.4245E-03  2.0031E+00  7.9268E-01  1.5763E+00 -2.0750E+00 -2.5883E-04 -3.0692E-03 -4.9794E-02  1.7723E-03 -3.3549E-03
            -1.5135E-03

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1603.39116343727        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:     1039
 NPARAMETR:  9.5181E-01  1.0927E+00  1.1542E+00  9.6395E-01  1.0624E+00  1.0187E+00  1.4253E+00  1.2818E+00  9.9696E-01  9.3638E-01
             1.0548E+00
 PARAMETER:  5.0608E-02  1.8868E-01  2.4360E-01  6.3387E-02  1.6044E-01  1.1855E-01  4.5441E-01  3.4821E-01  9.6958E-02  3.4264E-02
             1.5338E-01
 GRADIENT:  -1.7997E-02 -1.2267E+06  7.7446E-01  1.5492E+00 -2.0229E+00 -4.6220E-03 -5.2658E-03 -4.5284E-02  1.4618E-03 -3.4165E-03
            -1.1265E-03
 NUMSIGDIG:         4.1         3.3         1.7         1.6         1.7         4.0         4.1         2.3         4.0         3.6
                    4.8

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1039
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.0857E-04 -3.1866E-03 -3.7304E-02 -4.2468E-03 -3.4191E-02
 SE:             2.9835E-02  2.2641E-02  1.4321E-02  2.1694E-02  2.0213E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8907E-01  8.8807E-01  9.1901E-03  8.4480E-01  9.0740E-02

 ETASHRINKSD(%)  5.0039E-02  2.4150E+01  5.2024E+01  2.7324E+01  3.2283E+01
 ETASHRINKVR(%)  1.0005E-01  4.2467E+01  7.6983E+01  4.7182E+01  5.4144E+01
 EBVSHRINKSD(%)  4.6994E-01  2.4543E+01  5.5863E+01  2.8200E+01  2.9471E+01
 EBVSHRINKVR(%)  9.3767E-01  4.3063E+01  8.0519E+01  4.8447E+01  5.0257E+01
 RELATIVEINF(%)  9.8594E+01  3.3306E+00  2.5350E+00  3.0140E+00  1.1533E+01
 EPSSHRINKSD(%)  4.4643E+01
 EPSSHRINKVR(%)  6.9356E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1603.3911634372748     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -868.24033687353665     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.79
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1603.391       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.52E-01  1.09E+00  1.15E+00  9.64E-01  1.06E+00  1.02E+00  1.43E+00  1.28E+00  9.97E-01  9.36E-01  1.05E+00
 


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
+        1.17E+03
 
 TH 2
+        2.95E+09  2.72E+09
 
 TH 3
+        5.23E+00  9.98E+08  7.32E+08
 
 TH 4
+       -6.98E+00  2.91E+09  2.13E+09  6.26E+02
 
 TH 5
+       -2.22E+00 -1.65E+09 -1.21E+09 -3.52E+09  1.99E+09
 
 TH 6
+        7.99E-01 -5.50E+03  1.88E+00 -1.21E+00 -1.97E+00  1.89E+02
 
 TH 7
+        7.66E-01  4.33E+08  3.42E-01 -1.32E+01 -3.32E-01 -1.22E-01  3.81E+01
 
 TH 8
+       -7.78E-01 -6.29E+08 -4.61E+08  9.07E+00  3.87E-01 -5.64E-02  2.15E+00  1.58E+01
 
 TH 9
+        1.32E-02  2.05E+04 -1.11E+01  9.87E+00  5.51E+00 -9.25E-01  1.68E+01  4.17E+00  6.78E+01
 
 TH10
+        2.09E+00  7.39E+03 -7.71E+00 -1.32E+01 -7.30E+01  3.48E+00  4.81E+00  7.25E+00  3.14E+00  6.77E+01
 
 TH11
+       -7.81E+00  1.74E+09 -7.68E+00 -2.28E+00 -1.01E+01  2.56E+00  1.83E+00  3.41E+00  6.18E+00  1.31E+01  1.93E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.361
Stop Time:
Sat Sep 18 14:43:10 CDT 2021
