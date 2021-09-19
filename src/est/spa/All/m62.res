Sat Sep 18 15:56:25 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/All/dat62.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK

ET1 = EXP(ETA(1)*THETA(4))
ET2 = EXP(ETA(2)*THETA(5))
ET3 = EXP(ETA(3)*THETA(6))


CL = 5.0 * THETA(1) * ET1
V = 85  * THETA(2) * ET2
KA = 0.7 * THETA(3) * ET3

SC = V
$ERROR
CVERR 	= 0.05
W  	= THETA(7)*F*CVERR
Y  	= F + W * ERR(1)
$THETA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvKA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvK
(0,1) ; RUV
$OMEGA
0.9 FIX ;     IIV CL
0.9 FIX  ;     IIV V
0.9 FIX ;      IIV KA
$SIGMA  1  FIX;        [P]
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
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
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
0INITIAL ESTIMATE OF OMEGA:
 0.9000E+00
 0.0000E+00   0.9000E+00
 0.0000E+00   0.0000E+00   0.9000E+00
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

 ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
   ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3

 TRANSLATOR WILL CONVERT PARAMETERS
 CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            4           *           *           *           *
    3            *           -           -           -           -
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   9668.65034163109        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   1.2512E+02 -1.8252E+01  1.0855E+01 -2.8297E+02  5.6204E+01 -1.0901E+02 -2.1156E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -533.474052885074        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:       60
 NPARAMETR:  1.6976E+00  1.7435E+00  2.3857E+00  5.6924E-01  6.4024E-01  5.0439E-01  1.5420E+01
 PARAMETER:  6.2924E-01  6.5589E-01  9.6949E-01 -4.6345E-01 -3.4592E-01 -5.8441E-01  2.8357E+00
 GRADIENT:   6.8293E+01  1.4439E+01 -4.5927E+00 -1.2082E+01  3.5221E+01  1.3866E+00  1.3226E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -581.154171883666        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:      115
 NPARAMETR:  1.1547E+00  1.2465E+00  3.5422E+00  3.4985E-01  1.9911E-01  4.5763E-01  1.3155E+01
 PARAMETER:  2.4388E-01  3.2033E-01  1.3648E+00 -9.5026E-01 -1.5139E+00 -6.8170E-01  2.6768E+00
 GRADIENT:  -8.9059E+01  1.3520E+02 -2.6063E+00 -3.8171E+01  1.1978E+00  6.5382E-02  1.1066E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -593.587292443287        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      166
 NPARAMETR:  1.1271E+00  1.0764E+00  1.1289E+01  3.4502E-01  9.6785E-02  9.1959E-01  1.2681E+01
 PARAMETER:  2.1963E-01  1.7361E-01  2.5238E+00 -9.6416E-01 -2.2353E+00  1.6169E-02  2.6401E+00
 GRADIENT:   3.1502E+00  2.6609E+01 -5.5516E-01 -7.8747E+00 -3.5427E-02 -7.2093E-02 -1.6642E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -594.953001340909        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      218
 NPARAMETR:  1.0962E+00  1.0323E+00  2.0346E+02  3.4189E-01  2.8038E-02  8.0575E+00  1.2813E+01
 PARAMETER:  1.9185E-01  1.3174E-01  5.4154E+00 -9.7328E-01 -3.4742E+00  2.1866E+00  2.6505E+00
 GRADIENT:  -8.6122E-01  8.0420E-01  9.6631E-02 -7.2633E-01  6.8131E-03 -1.1549E-01  1.3375E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -594.972088903979        NO. OF FUNC. EVALS.:  90
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  1.1002E+00  1.0353E+00  2.1892E+02  3.4389E-01  2.7477E-02  8.9434E+00  1.2857E+01
 PARAMETER:  1.9545E-01  1.3472E-01  5.4887E+00 -9.6743E-01 -3.4944E+00  2.2909E+00  2.6539E+00
 GRADIENT:   4.4959E-01  7.4684E-02  6.5235E-02  1.1275E-01 -7.0480E-03 -6.5637E-02 -3.0962E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -595.051975035888        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      427
 NPARAMETR:  1.1050E+00  1.0408E+00  1.0269E+02  3.4530E-01  3.7832E-02  7.8006E+00  1.2881E+01
 PARAMETER:  1.9981E-01  1.3994E-01  4.7317E+00 -9.6333E-01 -3.1746E+00  2.1542E+00  2.6558E+00
 GRADIENT:   1.6439E+00  1.3968E+00  1.2404E-01  5.7849E-01 -4.3538E-03 -2.0395E-01  6.4764E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -595.065013320201        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      544
 NPARAMETR:  1.1014E+00  1.0373E+00  9.2438E+01  3.4424E-01  3.9141E-02  7.6123E+00  1.2847E+01
 PARAMETER:  1.9655E-01  1.3662E-01  4.6265E+00 -9.6642E-01 -3.1406E+00  2.1298E+00  2.6531E+00
 GRADIENT:  -1.0515E-01 -2.2127E-01  6.5548E-02 -2.7030E-02  1.6789E-02 -1.1670E-01  1.6039E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -595.065758367351        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      622
 NPARAMETR:  1.1008E+00  1.0368E+00  9.1827E+01  3.4356E-01  3.5088E-02  7.5905E+00  1.2851E+01
 PARAMETER:  1.9600E-01  1.3610E-01  4.6199E+00 -9.6840E-01 -3.2499E+00  2.1269E+00  2.6534E+00
 GRADIENT:   6.3777E-01  8.3901E-01  5.6641E-02  2.8118E-01  1.7819E-02  4.7606E-02  2.1672E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -595.066262167791        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.1007E+00  1.0362E+00  9.1835E+01  3.4389E-01  3.2599E-02  7.5900E+00  1.2850E+01
 PARAMETER:  1.9597E-01  1.3560E-01  4.6200E+00 -9.6743E-01 -3.3235E+00  2.1268E+00  2.6534E+00
 GRADIENT:   8.9336E-03 -4.6370E-02 -4.8302E-03  2.6239E-02  4.2980E-04  1.4392E-02 -1.7545E-01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -595.066288799270        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      805
 NPARAMETR:  1.1007E+00  1.0362E+00  9.1921E+01  3.4382E-01  3.1937E-02  7.5867E+00  1.2860E+01
 PARAMETER:  1.9597E-01  1.3558E-01  4.6199E+00 -9.6762E-01 -3.3406E+00  2.1269E+00  2.6536E+00
 GRADIENT:  -4.8253E-02 -7.0082E-02 -5.4711E+00  8.5809E-03  4.4224E-04  1.1965E+01 -9.3141E+00
 NUMSIGDIG:         3.8         3.6         3.3         4.5         2.6         3.3         3.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      805
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.5859E-02 -1.2360E-02 -8.8702E-03
 SE:             8.6082E-02  9.0971E-03  6.5331E-03
 N:                     100         100         100

 P VAL.:         6.7699E-01  1.7427E-01  1.7455E-01

 ETASHRINKSD(%)  8.8049E+00  9.0363E+01  9.3079E+01
 ETASHRINKVR(%)  1.6835E+01  9.9071E+01  9.9521E+01
 EBVSHRINKSD(%)  8.7277E+00  9.0729E+01  9.4943E+01
 EBVSHRINKVR(%)  1.6694E+01  9.9141E+01  9.9744E+01
 RELATIVEINF(%)  5.0960E+01  5.2865E-01  1.9642E-01
 EPSSHRINKSD(%)  4.1525E+00
 EPSSHRINKVR(%)  8.1325E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -595.06628879927018     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       140.08453776446800     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     7.05
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     2.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -595.066       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.10E+00  1.04E+00  9.18E+01  3.44E-01  3.20E-02  7.59E+00  1.29E+01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-01
 
 ETA2
+        0.00E+00  9.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.49E-01
 
 ETA2
+        0.00E+00  9.49E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.49E-01
 


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
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        6.72E+02
 
 TH 2
+       -3.30E+02  9.20E+02
 
 TH 3
+        1.07E-01  3.57E-01  2.89E-02
 
 TH 4
+        6.95E+01 -2.46E+02  1.41E-01  1.16E+03
 
 TH 5
+       -3.24E+01 -1.27E+02 -1.48E-01 -3.35E+01  5.55E+01
 
 TH 6
+       -2.72E+00 -9.45E+00 -7.69E-01 -3.81E+00  3.96E+00  2.04E+01
 
 TH 7
+       -1.01E+01 -1.25E+01  4.29E-01  1.09E+01  3.77E+00  1.34E-01  8.68E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,        9.274
Stop Time:
Sat Sep 18 15:56:36 CDT 2021
