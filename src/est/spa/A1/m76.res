Wed Sep 29 12:20:50 CDT 2021
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
$DATA ../../../../data/spa/A1/dat76.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m76.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -956.948498452628        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1985E+02  7.6481E+01  4.1634E+01  7.1061E+01  2.5671E+01 -5.8053E+00 -1.2127E+01 -9.0349E+00 -2.0210E+01 -5.1286E+01
            -1.2400E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1345.97033823036        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.2086E+00  9.9473E-01  1.0878E+00  1.0281E+00  9.6553E-01  1.3217E+00  9.5776E-01  9.3853E-01  1.0465E+00  8.5503E-01
             2.7027E+00
 PARAMETER:  2.8950E-01  9.4717E-02  1.8414E-01  1.2769E-01  6.4922E-02  3.7888E-01  5.6845E-02  3.6562E-02  1.4543E-01 -5.6617E-02
             1.0943E+00
 GRADIENT:   4.3467E+02  2.4967E+01  2.7916E+01  8.5302E+00 -5.1516E+01  3.1391E+01  5.5132E+00  5.0883E-01  1.0304E+01  6.9370E+00
             2.6655E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1378.44094618094        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0076E+00  7.9703E-01  6.1690E-01  1.1240E+00  7.0336E-01  1.0233E+00  1.2031E+00  3.8626E-01  9.1899E-01  5.6948E-01
             2.4092E+00
 PARAMETER:  1.0753E-01 -1.2687E-01 -3.8305E-01  2.1688E-01 -2.5189E-01  1.2306E-01  2.8489E-01 -8.5124E-01  1.5518E-02 -4.6303E-01
             9.7928E-01
 GRADIENT:   1.2707E+02 -2.6288E+00 -5.2868E+01  9.2319E+01  8.1165E+01 -1.3997E+01  3.4131E+00  1.4016E-01 -3.1971E+00 -1.3737E+00
            -3.2020E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1381.29326156624        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  9.6956E-01  7.9899E-01  6.2442E-01  1.0895E+00  6.8222E-01  1.0152E+00  1.1081E+00  5.2017E-01  9.5758E-01  5.4642E-01
             2.4028E+00
 PARAMETER:  6.9089E-02 -1.2440E-01 -3.7094E-01  1.8573E-01 -2.8240E-01  1.1505E-01  2.0264E-01 -5.5359E-01  5.6656E-02 -5.0436E-01
             9.7665E-01
 GRADIENT:  -2.4861E+01 -4.4490E-01 -9.7846E+00  4.6165E+00  7.1257E+00 -2.3682E+01 -1.4956E+00  5.3391E-01  2.8755E-01 -6.8854E-01
            -3.5621E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1384.93938324206        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      430
 NPARAMETR:  9.7801E-01  6.7267E-01  9.1225E-01  1.2005E+00  8.0562E-01  1.0678E+00  1.0866E+00  2.8279E-01  9.2615E-01  6.8946E-01
             2.6410E+00
 PARAMETER:  7.7769E-02 -2.9650E-01  8.1617E-03  2.8271E-01 -1.1615E-01  1.6562E-01  1.8306E-01 -1.1631E+00  2.3285E-02 -2.7184E-01
             1.0712E+00
 GRADIENT:  -2.0161E+00  2.7367E+00 -4.5613E-01  2.5152E+00 -8.7189E-01 -3.3937E-01 -2.5727E-01  3.1704E-01  9.0085E-01  1.4325E+00
             5.6156E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1385.91437312287        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      607
 NPARAMETR:  9.7493E-01  3.6645E-01  8.9986E-01  1.3724E+00  7.0314E-01  1.0669E+00  1.6557E+00  1.4066E-02  8.3700E-01  6.9239E-01
             2.6051E+00
 PARAMETER:  7.4609E-02 -9.0389E-01 -5.5159E-03  4.1654E-01 -2.5219E-01  1.6472E-01  6.0420E-01 -4.1640E+00 -7.7935E-02 -2.6760E-01
             1.0575E+00
 GRADIENT:   1.2743E+00  1.3506E+00  1.3296E+00  8.5009E-01 -2.3324E+00  5.8028E-02  1.7404E-01  9.9291E-04  2.8557E-01  2.1968E-01
             1.2869E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1386.07987765833        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      785
 NPARAMETR:  9.6992E-01  2.1279E-01  8.9089E-01  1.4587E+00  6.6003E-01  1.0663E+00  2.3601E+00  1.0000E-02  8.0319E-01  7.0421E-01
             2.5842E+00
 PARAMETER:  6.9463E-02 -1.4475E+00 -1.5530E-02  4.7757E-01 -3.1547E-01  1.6421E-01  9.5869E-01 -7.7091E+00 -1.1916E-01 -2.5068E-01
             1.0494E+00
 GRADIENT:  -2.1086E+00  1.5144E+00  2.5328E+00  7.7364E+00 -4.7462E+00  2.5983E-01  3.9468E-01  0.0000E+00 -1.5523E-01  1.6688E-02
            -2.9246E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1386.13624506920        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      968             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6906E-01  1.3966E-01  9.0979E-01  1.4959E+00  6.5576E-01  1.0642E+00  2.8085E+00  1.0000E-02  7.9094E-01  7.1513E-01
             2.5847E+00
 PARAMETER:  6.8572E-02 -1.8686E+00  5.4606E-03  5.0271E-01 -3.2196E-01  1.6226E-01  1.1326E+00 -1.0536E+01 -1.3453E-01 -2.3530E-01
             1.0496E+00
 GRADIENT:   6.2003E+01  3.1895E+00 -4.8145E-02  1.2437E+02  7.0251E+00  1.0188E+01  7.4074E-01  0.0000E+00  3.1965E+00  2.2284E-01
             8.6371E+00

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1386.13658202798        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  9.6896E-01  1.3987E-01  9.1070E-01  1.4965E+00  6.5539E-01  1.0641E+00  2.8146E+00  1.0000E-02  7.9006E-01  7.1481E-01
             2.5842E+00
 PARAMETER:  6.8466E-02 -1.8670E+00  6.4606E-03  5.0313E-01 -3.2252E-01  1.6216E-01  1.1348E+00 -1.0536E+01 -1.3565E-01 -2.3574E-01
             1.0494E+00
 GRADIENT:   1.0545E-01  1.5919E-01  5.5433E-01 -7.4640E-01 -3.7900E-01  2.1270E-02 -6.7174E-04  0.0000E+00  1.0794E-02 -1.8507E-02
             1.2968E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1061
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.0292E-05  4.4598E-04  8.1472E-05 -1.2813E-02 -1.7440E-02
 SE:             2.9152E-02  6.4062E-03  1.5906E-04  2.6138E-02  1.7615E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9890E-01  9.4450E-01  6.0850E-01  6.2400E-01  3.2213E-01

 ETASHRINKSD(%)  2.3359E+00  7.8539E+01  9.9467E+01  1.2434E+01  4.0987E+01
 ETASHRINKVR(%)  4.6173E+00  9.5394E+01  9.9997E+01  2.3321E+01  6.5175E+01
 EBVSHRINKSD(%)  2.1960E+00  7.9390E+01  9.9438E+01  1.1980E+01  4.1491E+01
 EBVSHRINKVR(%)  4.3437E+00  9.5752E+01  9.9997E+01  2.2525E+01  6.5767E+01
 RELATIVEINF(%)  8.6864E+01  8.7539E-02  1.7363E-04  2.4086E+00  1.2941E+00
 EPSSHRINKSD(%)  2.9787E+01
 EPSSHRINKVR(%)  5.0702E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1386.1365820279830     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -650.98575546424479     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.20
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1386.137       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.69E-01  1.40E-01  9.11E-01  1.50E+00  6.55E-01  1.06E+00  2.81E+00  1.00E-02  7.90E-01  7.15E-01  2.58E+00
 


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
+        1.00E+03
 
 TH 2
+       -5.40E+01  3.31E+02
 
 TH 3
+        3.37E+00  1.59E+02  4.14E+02
 
 TH 4
+       -4.23E+01  3.55E+02 -3.27E+01  6.16E+02
 
 TH 5
+        2.25E+01 -4.25E+02 -7.91E+02 -1.20E+02  1.65E+03
 
 TH 6
+        1.49E+00 -7.51E+00  5.72E+00 -1.22E+01 -2.76E+00  1.58E+02
 
 TH 7
+        3.26E-02  3.37E+00  4.89E-01 -3.97E-01  1.58E-01  9.08E-02  3.69E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.32E+00 -2.43E+01  2.10E+01 -1.31E+01 -6.79E+00  1.98E+00  1.49E+00  0.00E+00  1.97E+02
 
 TH10
+       -6.37E+00  3.92E+00 -3.93E+00 -8.62E+00 -1.73E+01 -5.36E-01  2.25E-01  0.00E+00  1.41E+00  4.27E+01
 
 TH11
+       -1.29E+01 -1.90E+00 -5.14E+00 -8.37E+00 -1.05E+01  2.87E+00  2.73E-01  0.00E+00  1.06E+01  2.41E+01  4.69E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.692
Stop Time:
Wed Sep 29 12:21:11 CDT 2021
