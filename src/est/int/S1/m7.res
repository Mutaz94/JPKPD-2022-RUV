Sat Sep 18 02:17:25 CDT 2021
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
$DATA ../../../../data/int/S1/dat7.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3325.61521025157        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.2864E+01 -2.0802E+01  1.1942E+02  9.5369E+01  8.6429E+01 -3.7641E+01 -6.5065E+01 -4.7314E+02 -1.3345E+02 -4.5087E+01
            -3.3742E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3501.21018077319        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:       92
 NPARAMETR:  8.8311E-01  1.0501E+00  8.8693E-01  8.6951E-01  9.2793E-01  1.0450E+00  1.0757E+00  1.9649E+00  1.1608E+00  1.0611E+00
             1.4686E+00
 PARAMETER: -2.4306E-02  1.4884E-01 -1.9993E-02 -3.9827E-02  2.5203E-02  1.4405E-01  1.7301E-01  7.7545E-01  2.4912E-01  1.5933E-01
             4.8430E-01
 GRADIENT:  -2.0052E+02 -4.9058E+01 -5.6968E+00 -9.6756E+01 -1.2682E+02 -3.9621E+01  1.7227E+01 -1.0469E+02  1.9782E+01 -2.8173E+00
             4.4683E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3503.63009015362        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  8.8312E-01  1.0501E+00  8.9394E-01  8.6949E-01  9.2797E-01  1.1472E+00  1.0757E+00  1.9648E+00  1.1608E+00  1.0612E+00
             1.4685E+00
 PARAMETER: -2.4292E-02  1.4886E-01 -1.2119E-02 -3.9846E-02  2.5248E-02  2.3732E-01  1.7298E-01  7.7537E-01  2.4911E-01  1.5938E-01
             4.8422E-01
 GRADIENT:  -1.6254E+02 -4.9060E+01 -3.6175E+00 -9.7733E+01 -1.2899E+02  6.0100E+00  1.7149E+01 -1.0535E+02  1.9818E+01 -2.8163E+00
             4.4692E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3503.93503786556        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      421
 NPARAMETR:  8.8313E-01  1.0501E+00  9.1334E-01  8.6949E-01  9.2798E-01  1.1422E+00  1.0757E+00  1.9652E+00  1.1608E+00  1.0612E+00
             1.4676E+00
 PARAMETER: -2.4285E-02  1.4887E-01  9.3561E-03 -3.9844E-02  2.5256E-02  2.3292E-01  1.7298E-01  7.7560E-01  2.4910E-01  1.5937E-01
             4.8360E-01
 GRADIENT:  -1.8675E+02 -5.5334E+01  1.5837E+00 -1.0370E+02 -1.3691E+02 -6.9410E+00  1.6275E+01 -1.0907E+02  1.8990E+01 -3.2985E+00
             4.4573E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3512.26027786210        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      602
 NPARAMETR:  8.8335E-01  1.0505E+00  9.5327E-01  8.6956E-01  9.2827E-01  1.0396E+00  1.0753E+00  1.9825E+00  1.1602E+00  1.0609E+00
             1.4337E+00
 PARAMETER: -2.4029E-02  1.4925E-01  5.2141E-02 -3.9766E-02  2.5572E-02  1.3881E-01  1.7264E-01  7.8437E-01  2.4857E-01  1.5913E-01
             4.6027E-01
 GRADIENT:  -2.2444E+02 -5.1741E+01  1.1473E+01 -1.0588E+02 -1.4895E+02 -4.8811E+01  1.5712E+01 -1.1123E+02  1.9108E+01 -4.2889E+00
             4.1701E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3524.27906397431        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      765
 NPARAMETR:  9.1699E-01  1.0505E+00  9.3746E-01  8.9445E-01  9.2828E-01  1.0737E+00  1.0269E+00  1.9826E+00  1.1190E+00  1.0734E+00
             1.4335E+00
 PARAMETER:  1.3338E-02  1.4924E-01  3.5418E-02 -1.1549E-02  2.5583E-02  1.7115E-01  1.2658E-01  7.8439E-01  2.1243E-01  1.7080E-01
             4.6009E-01
 GRADIENT:  -1.3928E+02 -4.5717E+01  4.8986E+00 -7.5168E+01 -1.3634E+02 -2.0540E+01  1.1920E+01 -1.1166E+02  1.8285E+01 -5.0345E+00
             4.1405E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3531.63396565203        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:      933
 NPARAMETR:  9.8395E-01  1.0504E+00  9.3457E-01  9.4359E-01  9.2826E-01  1.0971E+00  1.0270E+00  1.9821E+00  1.1191E+00  1.0944E+00
             1.4333E+00
 PARAMETER:  8.3815E-02  1.4920E-01  3.2332E-02  4.1934E-02  2.5554E-02  1.9265E-01  1.2661E-01  7.8416E-01  2.1249E-01  1.9025E-01
             4.5995E-01
 GRADIENT:  -4.8203E+00 -2.0452E+01  1.0960E-02 -1.5572E+00 -1.1949E+02 -1.7253E+00  1.3608E+01 -1.1124E+02  2.4195E+01 -2.2654E-01
             4.1181E+02

0ITERATION NO.:   33    OBJECTIVE VALUE:  -3531.64051602315        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1032
 NPARAMETR:  9.8479E-01  1.0504E+00  9.3466E-01  9.4389E-01  9.2826E-01  1.0989E+00  1.0270E+00  1.9821E+00  1.1191E+00  1.0949E+00
             1.4333E+00
 PARAMETER:  8.4677E-02  1.4920E-01  3.2429E-02  4.2261E-02  2.5554E-02  1.9409E-01  1.2661E-01  7.8416E-01  2.1249E-01  1.9064E-01
             4.5995E-01
 GRADIENT:   7.1714E+05 -2.4034E+05  7.1716E+05  7.1713E+05 -3.5868E+05 -1.1361E+00  2.8321E+05 -4.5878E+04  1.6877E+05 -1.8810E+05
            -1.5552E+05
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         1.5         3.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1032
 NO. OF SIG. DIGITS IN FINAL EST.:  1.5
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3230E-03 -2.6130E-02 -3.4975E-02  1.8240E-02  1.4725E-02
 SE:             2.9927E-02  2.0438E-02  3.0955E-02  2.4641E-02  2.2976E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3813E-01  2.0108E-01  2.5853E-01  4.5916E-01  5.2161E-01

 ETASHRINKSD(%)  1.0000E-10  3.1529E+01  1.0000E-10  1.7450E+01  2.3028E+01
 ETASHRINKVR(%)  1.0000E-10  5.3117E+01  1.0000E-10  3.1855E+01  4.0753E+01
 EBVSHRINKSD(%)  4.2849E-01  2.7461E+01  2.8377E+01  1.0832E+01  2.1591E+01
 EBVSHRINKVR(%)  8.5514E-01  4.7381E+01  4.8702E+01  2.0491E+01  3.8521E+01
 RELATIVEINF(%)  9.9140E+01  2.1671E+01  3.8383E+01  4.6917E+01  2.4085E+01
 EPSSHRINKSD(%)  3.7274E+01
 EPSSHRINKVR(%)  6.0655E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3531.6405160231502     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1877.5511562547395     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.24
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3531.641       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  1.05E+00  9.35E-01  9.44E-01  9.28E-01  1.10E+00  1.03E+00  1.98E+00  1.12E+00  1.09E+00  1.43E+00
 


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
+        1.85E+09
 
 TH 2
+       -1.16E+09  7.30E+08
 
 TH 3
+        1.95E+09 -1.22E+09  2.05E+09
 
 TH 4
+        1.93E+09  6.53E+04  2.03E+09  2.01E+09
 
 TH 5
+       -1.96E+09 -6.64E+04 -2.07E+09 -2.05E+09  2.08E+09
 
 TH 6
+        1.84E+03 -5.36E+08  1.94E+03  1.91E+03 -1.95E+03  1.63E+02
 
 TH 7
+        1.05E+03  4.72E+04  1.48E+09  1.46E+09 -1.49E+09  1.39E+03  1.06E+09
 
 TH 8
+       -8.77E+01 -3.95E+03 -5.37E+03  3.28E+03 -5.41E+03 -1.17E+02  8.02E+03  7.43E+06
 
 TH 9
+        7.66E+08  2.58E+04  8.07E+08  7.99E+08 -8.12E+08  7.62E+02  5.80E+08  1.78E+03  3.17E+08
 
 TH10
+       -8.72E+08  5.48E+08 -9.19E+08 -9.10E+08  9.25E+08 -8.68E+02 -5.24E+04  4.39E+03 -3.61E+08  4.12E+08
 
 TH11
+       -2.14E+02 -9.32E+03 -1.26E+04  7.72E+03 -1.27E+04 -2.73E+02  1.89E+04  3.42E+03  4.17E+03  1.04E+04  4.13E+07
 
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
 #CPUT: Total CPU Time in Seconds,       45.516
Stop Time:
Sat Sep 18 02:18:12 CDT 2021
