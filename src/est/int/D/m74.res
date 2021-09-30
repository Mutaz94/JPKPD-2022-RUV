Wed Sep 29 09:40:28 CDT 2021
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
$DATA ../../../../data/int/D/dat74.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   30491.6445758188        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4756E+02  5.3845E+02 -5.1623E+01  2.8153E+02  3.2423E+02 -3.0567E+03 -1.2382E+03 -8.9727E+01 -2.1712E+03 -8.9020E+02
            -6.0637E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -886.535665382960        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.4547E+00  1.5348E+00  8.2592E-01  2.6720E+00  9.4284E-01  5.2408E+00  3.6694E+00  1.0191E+00  3.4736E+00  1.6556E+00
             1.1622E+01
 PARAMETER:  4.7480E-01  5.2837E-01 -9.1255E-02  1.0828E+00  4.1140E-02  1.7565E+00  1.4000E+00  1.1888E-01  1.3452E+00  6.0417E-01
             2.5529E+00
 GRADIENT:   2.2494E+01 -2.8555E+00 -4.6168E+01  1.2381E+02 -1.8492E+01  2.2387E+02 -1.4080E+01  4.7784E+00  3.6253E+01  3.9539E+01
             3.8362E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -980.206232932822        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      179
 NPARAMETR:  9.4248E-01  3.8816E+00  3.1458E+01  5.3139E+00  4.4108E+00  3.7198E+00  2.1432E+00  1.1944E+00  2.9465E+01  9.4015E-01
             1.0981E+01
 PARAMETER:  4.0758E-02  1.4562E+00  3.5487E+00  1.7703E+00  1.5841E+00  1.4137E+00  8.6231E-01  2.7762E-01  3.4832E+00  3.8289E-02
             2.4962E+00
 GRADIENT:  -4.4888E+01 -2.1664E+01 -1.7367E+01  1.4261E+01  7.6393E+01  8.1221E+01  1.3509E+01 -9.4921E-01  4.0151E+01  7.1327E+00
             3.3875E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1045.21554021518        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      357
 NPARAMETR:  1.1103E+00  3.1136E+00  1.3783E+01  3.8281E+00  2.0815E+00  2.8039E+00  2.1929E+00  4.2725E+00  2.2180E+01  8.6072E-01
             1.0547E+01
 PARAMETER:  2.0461E-01  1.2358E+00  2.7234E+00  1.4424E+00  8.3310E-01  1.1310E+00  8.8520E-01  1.5522E+00  3.1992E+00 -4.9984E-02
             2.4559E+00
 GRADIENT:  -3.0363E+01 -3.8313E+01 -1.9277E+01  1.8490E+01 -5.3771E+01  1.9251E+01  3.7312E+01 -2.0732E+00  4.8110E+01  9.5849E+00
             3.3156E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1189.09482366985        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  1.2150E+00  1.6756E+00  5.2608E+01  1.0310E+00  2.3669E+00  3.1306E+00  5.9457E+00  3.7942E+00  1.1017E+00  4.0739E-01
             8.2999E+00
 PARAMETER:  2.9475E-01  6.1617E-01  4.0629E+00  1.3050E-01  9.6159E-01  1.2412E+00  1.8827E+00  1.4335E+00  1.9688E-01 -7.9798E-01
             2.2162E+00
 GRADIENT:  -1.8283E+00  1.2931E+01  6.1872E-01 -8.2008E+00 -2.6785E+01  1.0806E+01  2.0377E+01  3.8197E-02  2.4783E+00  2.1992E+00
            -4.0208E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1192.95863178830        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      714
 NPARAMETR:  1.2370E+00  1.4469E+00  6.1934E+01  1.0592E+00  2.5166E+00  3.0258E+00  5.5303E+00  3.8312E+00  1.0093E+00  2.4661E-01
             8.5509E+00
 PARAMETER:  3.1270E-01  4.6942E-01  4.2261E+00  1.5750E-01  1.0229E+00  1.2072E+00  1.8102E+00  1.4432E+00  1.0928E-01 -1.3000E+00
             2.2460E+00
 GRADIENT:  -4.8372E-01  3.2843E+00  5.4007E-02 -1.4601E+00  2.1906E+00 -1.1237E+00 -1.8899E+00  4.8937E-02  1.6654E-01  9.1347E-01
             8.8248E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1193.03545411321        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  1.2388E+00  1.4496E+00  6.0647E+01  1.0584E+00  2.5052E+00  3.0429E+00  5.5165E+00  4.1934E-01  1.0098E+00  2.4491E-01
             8.5349E+00
 PARAMETER:  3.1412E-01  4.7132E-01  4.2051E+00  1.5677E-01  1.0184E+00  1.2128E+00  1.8077E+00 -7.6907E-01  1.0978E-01 -1.3069E+00
             2.2442E+00
 GRADIENT:   4.4592E-02  3.4795E+00  1.4592E-01 -7.0988E-01 -6.2769E-02  9.9780E-01 -2.6385E+00  6.5092E-04  1.7154E-02  8.8966E-01
             5.1054E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1193.42709304312        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1077
 NPARAMETR:  1.2387E+00  1.3949E+00  5.2293E+01  1.0643E+00  2.4951E+00  3.0429E+00  5.6109E+00  8.4958E-02  1.0068E+00  1.8135E-01
             8.4738E+00
 PARAMETER:  3.1409E-01  4.3285E-01  4.0569E+00  1.6228E-01  1.0143E+00  1.2128E+00  1.8247E+00 -2.3656E+00  1.0681E-01 -1.6073E+00
             2.2370E+00
 GRADIENT:   6.3714E-01  2.3642E+00  1.6791E-01 -2.7078E+00 -1.8171E-01  7.2130E-01 -5.5537E-01  3.6172E-05 -2.7921E-01  4.7871E-01
            -7.8832E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1193.92602511813        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1268
 NPARAMETR:  1.2328E+00  1.2956E+00  2.7089E+01  1.0840E+00  2.4608E+00  3.0537E+00  5.6165E+00  6.9416E-02  1.0242E+00  1.9460E-02
             8.4746E+00
 PARAMETER:  3.0927E-01  3.5901E-01  3.3991E+00  1.8070E-01  1.0005E+00  1.2164E+00  1.8257E+00 -2.5676E+00  1.2392E-01 -3.8394E+00
             2.2371E+00
 GRADIENT:  -4.0053E-01 -6.9437E-01  2.1637E-02 -3.3521E+00  3.7010E+00  2.2514E+00 -4.6203E+00  1.0053E-04 -6.6867E-01  5.5470E-03
            -6.2102E+00

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1193.93587854261        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1368
 NPARAMETR:  1.2328E+00  1.2957E+00  2.7083E+01  1.0840E+00  2.4608E+00  3.0534E+00  5.6195E+00  6.3747E-02  1.0242E+00  1.4861E-02
             8.4845E+00
 PARAMETER:  3.0927E-01  3.5902E-01  3.3991E+00  1.8070E-01  1.0004E+00  1.2163E+00  1.8263E+00 -2.6266E+00  1.2392E-01 -4.0984E+00
             2.2383E+00
 GRADIENT:   6.4498E+00 -3.7069E+02  3.9089E+01  1.4664E+03 -2.6520E+02  1.8697E+01  1.3866E+02  6.1424E-05 -1.0724E+03  3.2285E-03
             1.0459E+02
 NUMSIGDIG:         4.7         2.3         2.3         2.3         2.3         3.7         2.3         0.0         2.3         0.6
                    2.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1368
 NO. OF SIG. DIGITS IN FINAL EST.:  0.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.8502E-03  2.6357E-02 -9.0398E-06 -6.3336E-02 -2.0293E-05
 SE:             2.9526E-02  2.6170E-02  2.8497E-05  1.2322E-02  1.9157E-04
 N:                     100         100         100         100         100

 P VAL.:         7.6437E-01  3.1386E-01  7.5108E-01  2.7511E-07  9.1564E-01

 ETASHRINKSD(%)  1.0850E+00  1.2328E+01  9.9905E+01  5.8720E+01  9.9358E+01
 ETASHRINKVR(%)  2.1582E+00  2.3137E+01  1.0000E+02  8.2959E+01  9.9996E+01
 EBVSHRINKSD(%)  1.8558E+00  8.8687E+00  9.9874E+01  6.6014E+01  9.9252E+01
 EBVSHRINKVR(%)  3.6773E+00  1.6951E+01  1.0000E+02  8.8449E+01  9.9994E+01
 RELATIVEINF(%)  9.6226E+01  4.2474E+01  3.5502E-05  5.9760E+00  1.2318E-03
 EPSSHRINKSD(%)  6.7116E+00
 EPSSHRINKVR(%)  1.2973E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1193.9358785426134     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       460.15348122579735     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    49.09
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.58
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1193.936       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.23E+00  1.30E+00  2.71E+01  1.08E+00  2.46E+00  3.05E+00  5.62E+00  6.54E-02  1.02E+00  1.50E-02  8.49E+00
 


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
+        9.14E+04
 
 TH 2
+       -3.14E+02  3.07E+04
 
 TH 3
+        1.59E+00 -1.55E+02  7.84E-01
 
 TH 4
+        7.43E+02 -7.46E-01  5.44E-02  1.73E+05
 
 TH 5
+       -7.02E+03 -1.67E+00 -1.63E-01  1.36E+04  1.17E+03
 
 TH 6
+        4.92E+03 -4.05E+03  2.05E+01  9.62E+03 -7.57E+02  4.48E+02
 
 TH 7
+        1.43E+01  1.71E+00 -3.00E-03 -1.17E+01 -5.99E+00  1.85E+02  6.63E+01
 
 TH 8
+        1.84E-02  9.53E-03 -1.68E-04 -2.78E-02 -1.78E-04 -2.51E-03  2.74E-04  5.27E-02
 
 TH 9
+       -1.15E+03 -1.83E+00 -5.68E+02 -3.08E+01  2.35E+00 -1.48E+04  2.16E+00  2.27E-02  4.12E+05
 
 TH10
+       -1.16E-01  9.76E-01 -4.83E-03 -2.13E+00  2.30E-01 -1.93E-03 -4.37E-02 -2.02E-01  3.06E+00  1.50E+01
 
 TH11
+        4.26E+00 -2.34E+00 -5.70E-04 -9.64E+00 -3.18E+00  1.12E+01  1.05E-01  4.60E-04  3.85E+00  2.96E-03  2.91E+01
 
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
 #CPUT: Total CPU Time in Seconds,       65.792
Stop Time:
Wed Sep 29 09:41:36 CDT 2021
