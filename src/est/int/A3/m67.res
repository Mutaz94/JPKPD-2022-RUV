Wed Sep 29 00:20:41 CDT 2021
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
$DATA ../../../../data/int/A3/dat67.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1466.90892959400        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2206E+02  1.6150E+02  3.0717E+02 -1.1309E+02  2.7784E+02  2.5768E+01 -1.9996E+02 -2.5622E+02 -5.1499E+01 -1.3433E+02
            -4.0934E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2877.97376262349        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0268E+00  1.0105E+00  8.6066E-01  1.1877E+00  8.8269E-01  9.6173E-01  9.9739E-01  9.2378E-01  9.8028E-01  8.8285E-01
             2.2889E+00
 PARAMETER:  1.2649E-01  1.1045E-01 -5.0050E-02  2.7202E-01 -2.4778E-02  6.0974E-02  9.7385E-02  2.0721E-02  8.0084E-02 -2.4597E-02
             9.2808E-01
 GRADIENT:   1.7995E+02  9.1309E+01  6.2789E+00  1.2984E+02  9.5387E-01 -5.8002E+00 -8.5700E+00  8.9818E+00 -6.9140E+00 -3.0199E+00
            -8.9173E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2893.25060311433        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      181
 NPARAMETR:  1.0196E+00  5.7432E-01  4.7285E-01  1.5099E+00  4.8588E-01  1.0003E+00  1.4099E+00  3.4365E-01  1.0887E+00  7.6941E-01
             2.1048E+00
 PARAMETER:  1.1946E-01 -4.5458E-01 -6.4897E-01  5.1203E-01 -6.2179E-01  1.0025E-01  4.4351E-01 -9.6813E-01  1.8502E-01 -1.6213E-01
             8.4422E-01
 GRADIENT:   7.9157E+01  8.2544E+01 -5.0819E+01  2.8442E+02  8.4597E+01 -1.3301E-01 -1.0827E-01 -5.3266E-01 -3.1538E+01  7.6799E+00
            -1.4766E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2920.21198406224        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      360
 NPARAMETR:  9.8085E-01  4.3938E-01  4.1087E-01  1.4886E+00  3.8750E-01  9.5389E-01  1.6059E+00  9.6593E-02  1.0901E+00  7.6021E-01
             2.2377E+00
 PARAMETER:  8.0659E-02 -7.2239E-01 -7.8948E-01  4.9782E-01 -8.4804E-01  5.2790E-02  5.7370E-01 -2.2373E+00  1.8623E-01 -1.7416E-01
             9.0544E-01
 GRADIENT:  -1.0139E+01  2.5575E+01  9.7863E+01  2.4165E+02 -7.0688E+01 -1.4687E+01  2.5739E+01 -2.9933E-02 -2.0507E+01  1.4252E+01
             2.6675E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2994.99266045415        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      537
 NPARAMETR:  9.8578E-01  2.0207E-01  1.3509E-01  1.0580E+00  1.7577E-01  1.0085E+00  1.5214E+00  4.2539E-02  1.4056E+00  7.5871E-01
             1.9895E+00
 PARAMETER:  8.5678E-02 -1.4992E+00 -1.9018E+00  1.5637E-01 -1.6386E+00  1.0851E-01  5.1961E-01 -3.0573E+00  4.4045E-01 -1.7613E-01
             7.8786E-01
 GRADIENT:   7.7639E+00  6.3922E+00  2.2591E+01  3.1031E+01 -4.5904E+01  2.2446E+00  1.3348E+00 -1.7569E-01 -6.0672E+00 -2.1556E+00
            -1.1972E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2996.84787810373        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      712
 NPARAMETR:  9.8139E-01  1.8353E-01  1.1529E-01  9.6770E-01  1.6228E-01  1.0019E+00  1.5296E+00  3.2446E-02  1.5443E+00  7.8561E-01
             1.9828E+00
 PARAMETER:  8.1219E-02 -1.5954E+00 -2.0603E+00  6.7168E-02 -1.7184E+00  1.0191E-01  5.2498E-01 -3.3282E+00  5.3456E-01 -1.4130E-01
             7.8452E-01
 GRADIENT:  -2.3460E-01  2.0207E-01 -2.5682E-01  1.2108E-01  2.8738E-01  8.0877E-02 -1.2835E-01 -1.0733E-01 -1.2340E-02 -5.5467E-02
            -9.1110E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3012.29767293415        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      894
 NPARAMETR:  9.7970E-01  1.7428E-01  1.1319E-01  9.6701E-01  1.6620E-01  1.0027E+00  1.5078E+00  9.6768E-01  1.5348E+00  7.8086E-01
             1.9726E+00
 PARAMETER:  7.9486E-02 -1.6471E+00 -2.0787E+00  6.6451E-02 -1.6946E+00  1.0271E-01  5.1064E-01  6.7142E-02  5.2838E-01 -1.4736E-01
             7.7935E-01
 GRADIENT:  -3.7744E-01 -3.9910E+01 -3.1506E+01 -2.3435E+01  9.0862E+01  2.9535E+00  1.5566E+00  5.2216E+00 -1.2305E+01  2.5829E+01
             7.4590E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3018.60894573435        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  9.7914E-01  1.8209E-01  1.1423E-01  9.9500E-01  1.6374E-01  9.9711E-01  1.4911E+00  9.3652E-01  1.5901E+00  7.1434E-01
             1.9150E+00
 PARAMETER:  7.8921E-02 -1.6032E+00 -2.0695E+00  9.4989E-02 -1.7095E+00  9.7103E-02  4.9954E-01  3.4413E-02  5.6383E-01 -2.3639E-01
             7.4972E-01
 GRADIENT:   2.6141E+00 -5.6958E+00 -2.2007E+01 -4.3455E+00  2.5058E+01 -2.7885E-01  7.2534E-01  3.5980E-01 -6.5422E+00  1.3718E+01
             1.5950E+01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -3018.66724881173        NO. OF FUNC. EVALS.:  68
 CUMULATIVE NO. OF FUNC. EVALS.:     1144
 NPARAMETR:  9.7909E-01  1.8234E-01  1.1418E-01  9.9534E-01  1.6376E-01  9.9722E-01  1.4910E+00  9.3609E-01  1.5906E+00  7.1265E-01
             1.9154E+00
 PARAMETER:  7.8866E-02 -1.6028E+00 -2.0688E+00  9.5385E-02 -1.7093E+00  9.7151E-02  4.9946E-01  3.4955E-02  5.6447E-01 -2.3876E-01
             7.4950E-01
 GRADIENT:   4.8074E+00 -8.0921E+02  1.2428E+03  2.5789E+04  4.4737E+01 -1.2903E+04 -1.2381E-01  3.3004E-01  4.5448E+03  2.3694E+01
            -3.4359E+03
 NUMSIGDIG:         6.3         2.3         2.3         2.3         4.1         2.3         7.2         1.1         2.3         5.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1144
 NO. OF SIG. DIGITS IN FINAL EST.:  1.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.0623E-04  2.1171E-02  3.1329E-02  2.1735E-03  1.2850E-02
 SE:             2.9475E-02  2.5819E-02  1.8486E-02  2.7552E-02  2.2501E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7818E-01  4.1225E-01  9.0124E-02  9.3712E-01  5.6795E-01

 ETASHRINKSD(%)  1.2559E+00  1.3502E+01  3.8070E+01  7.6979E+00  2.4618E+01
 ETASHRINKVR(%)  2.4960E+00  2.5180E+01  6.1647E+01  1.4803E+01  4.3176E+01
 EBVSHRINKSD(%)  1.3659E+00  1.1952E+01  3.7220E+01  7.3826E+00  2.1031E+01
 EBVSHRINKVR(%)  2.7131E+00  2.2475E+01  6.0587E+01  1.4220E+01  3.7639E+01
 RELATIVEINF(%)  9.7192E+01  3.3524E+01  7.6655E+00  5.2862E+01  1.1690E+01
 EPSSHRINKSD(%)  2.2317E+01
 EPSSHRINKVR(%)  3.9653E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3018.6672488117283     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1364.5778890433176     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.63
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3018.667       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.82E-01  1.14E-01  9.95E-01  1.64E-01  9.97E-01  1.49E+00  9.37E-01  1.59E+00  7.13E-01  1.91E+00
 


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
+        1.35E+07
 
 TH 2
+        3.91E+00  7.68E+05
 
 TH 3
+        2.83E+06  9.18E+05  1.22E+06
 
 TH 4
+        9.61E+00  2.07E+03 -3.38E+03  6.51E+06
 
 TH 5
+       -2.37E+06  7.87E+05 -2.92E+03 -6.37E+02  8.67E+05
 
 TH 6
+        1.89E+00 -2.02E+02  2.43E+02  5.64E+02 -1.46E+01  6.49E+06
 
 TH 7
+        8.85E+05 -2.97E+05  7.33E+05 -2.45E+00 -3.11E+05 -2.22E-01  1.16E+05
 
 TH 8
+       -1.76E+00  5.62E+03 -6.93E+03 -1.64E+04 -1.19E+02  1.64E+04  2.32E+00  7.37E+06
 
 TH 9
+        7.99E+01  2.47E+05 -6.07E+05 -7.25E+02 -2.56E+05  6.22E+01  9.69E+04 -1.81E+03  7.93E+04
 
 TH10
+       -7.75E+06 -1.30E+06 -1.63E+06  1.34E+01  1.36E+06  4.95E-01  3.26E+00  1.34E+01  4.11E+02  4.46E+06
 
 TH11
+       -6.25E+01 -1.54E+05  1.93E+05  4.61E+02 -1.54E+03 -3.56E+01  6.02E+04  1.16E+03 -2.35E+02 -2.50E+02  3.18E+04
 
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
 #CPUT: Total CPU Time in Seconds,       46.636
Stop Time:
Wed Sep 29 00:21:30 CDT 2021
