Sun Oct 24 00:20:18 CDT 2021
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
$DATA ../../../../data/SD3/SL3/dat86.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1998.85269263524        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6838E+02 -4.8083E+01 -3.9167E+01  2.5752E+01  1.3042E+02  5.5506E+01 -2.3014E+00 -3.7785E+00  1.1538E+01 -1.9068E+01
            -2.6775E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2038.76563744795        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  1.0284E+00  1.0274E+00  9.9436E-01  1.0099E+00  9.2374E-01  9.7743E-01  9.9596E-01  1.0084E+00  9.7810E-01  9.9735E-01
             1.2984E+00
 PARAMETER:  1.2798E-01  1.2698E-01  9.4348E-02  1.0983E-01  2.0678E-02  7.7167E-02  9.5949E-02  1.0832E-01  7.7858E-02  9.7342E-02
             3.6117E-01
 GRADIENT:   3.1542E+02  6.3259E-01  2.9577E+00  1.1608E+01  9.8172E+00  2.6787E+01  1.4985E+00 -1.0266E+00  6.9182E+00  4.3512E+00
             9.8865E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2039.55617113039        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      207
 NPARAMETR:  1.0263E+00  1.0544E+00  9.2815E-01  1.0183E+00  8.9192E-01  9.6841E-01  1.0134E+00  9.8990E-01  9.4397E-01  9.3316E-01
             1.2905E+00
 PARAMETER:  1.2597E-01  1.5298E-01  2.5433E-02  1.1811E-01 -1.4378E-02  6.7900E-02  1.1332E-01  8.9847E-02  4.2343E-02  3.0823E-02
             3.5501E-01
 GRADIENT:  -3.3080E+01  1.4148E+00  2.6815E+00  1.5395E-02 -5.4128E+00 -6.1123E+00 -1.9449E+00 -2.2255E-01 -2.6628E+00  5.4173E-01
            -3.5695E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2039.96964158710        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      386
 NPARAMETR:  1.0426E+00  1.0856E+00  8.4729E-01  9.9518E-01  8.7198E-01  9.7940E-01  1.0212E+00  8.8838E-01  9.6226E-01  9.1085E-01
             1.2926E+00
 PARAMETER:  1.4169E-01  1.8215E-01 -6.5712E-02  9.5166E-02 -3.6990E-02  7.9184E-02  1.2098E-01 -1.8358E-02  6.1533E-02  6.6227E-03
             3.5667E-01
 GRADIENT:   2.6332E+00 -4.7124E-01 -1.2937E+00  4.9605E+00  3.8620E-01 -1.4405E+00 -5.8435E-01  2.4373E-01  1.6523E-02  8.3540E-01
             1.0804E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2040.48432520096        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      562
 NPARAMETR:  1.0420E+00  1.3279E+00  6.3035E-01  8.3116E-01  8.7455E-01  9.9057E-01  9.1049E-01  5.9744E-01  1.0837E+00  8.7301E-01
             1.2899E+00
 PARAMETER:  1.4116E-01  3.8361E-01 -3.6148E-01 -8.4933E-02 -3.4043E-02  9.0524E-02  6.2277E-03 -4.1511E-01  1.8037E-01 -3.5813E-02
             3.5458E-01
 GRADIENT:  -3.7317E+00  1.8006E+00  7.0275E-01  1.6022E+00 -1.5131E+00  1.5631E+00  1.0920E-01  1.1712E-01  2.7204E-01 -4.0707E-01
            -7.3887E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2040.51348705030        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      740             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0453E+00  1.3798E+00  5.9971E-01  7.9607E-01  8.8609E-01  9.8725E-01  8.8733E-01  5.0753E-01  1.1154E+00  8.8074E-01
             1.2913E+00
 PARAMETER:  1.4432E-01  4.2194E-01 -4.1132E-01 -1.2807E-01 -2.0938E-02  8.7166E-02 -1.9540E-02 -5.7820E-01  2.0922E-01 -2.6988E-02
             3.5567E-01
 GRADIENT:   3.8938E+02  1.8115E+02  5.6401E+00  3.9312E+01  5.0371E+00  2.8085E+01  2.8314E+00  6.2575E-02  6.8303E+00  1.3548E-01
             3.3130E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2040.51561334045        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      920             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0450E+00  1.3795E+00  5.9885E-01  7.9627E-01  8.8590E-01  9.8714E-01  8.8697E-01  5.1301E-01  1.1150E+00  8.8172E-01
             1.2914E+00
 PARAMETER:  1.4402E-01  4.2173E-01 -4.1275E-01 -1.2781E-01 -2.1150E-02  8.7053E-02 -1.9942E-02 -5.6745E-01  2.0888E-01 -2.5876E-02
             3.5574E-01
 GRADIENT:   3.8779E+02  1.8040E+02  4.8726E+00  3.9690E+01  5.8516E+00  2.8055E+01  2.8396E+00  1.3007E-01  6.8513E+00  4.1187E-01
             3.5031E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2040.51575907893        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1087
 NPARAMETR:  1.0447E+00  1.3794E+00  5.9846E-01  7.9643E-01  8.8576E-01  9.8706E-01  8.8720E-01  5.1304E-01  1.1151E+00  8.8180E-01
             1.2914E+00
 PARAMETER:  1.4371E-01  4.2166E-01 -4.1340E-01 -1.2762E-01 -2.1309E-02  8.6980E-02 -1.9688E-02 -5.6739E-01  2.0890E-01 -2.5795E-02
             3.5576E-01
 GRADIENT:   1.3953E+00 -2.3946E-01  1.4469E-01 -3.8128E-02  2.6168E-01  5.3180E-02 -1.0103E-02  1.3132E-04  2.4779E-02  1.1001E-02
            -1.5349E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1087
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.0091E-04 -1.8237E-02 -1.2764E-02  1.1588E-02 -2.4299E-02
 SE:             2.9797E-02  2.1949E-02  6.7529E-03  2.4346E-02  2.1774E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7856E-01  4.0606E-01  5.8742E-02  6.3408E-01  2.6444E-01

 ETASHRINKSD(%)  1.7695E-01  2.6467E+01  7.7377E+01  1.8439E+01  2.7055E+01
 ETASHRINKVR(%)  3.5359E-01  4.5929E+01  9.4882E+01  3.3479E+01  4.6791E+01
 EBVSHRINKSD(%)  5.6155E-01  2.6412E+01  7.9417E+01  1.8880E+01  2.6559E+01
 EBVSHRINKVR(%)  1.1199E+00  4.5848E+01  9.5763E+01  3.4195E+01  4.6065E+01
 RELATIVEINF(%)  9.8742E+01  3.6004E+00  6.5508E-01  5.1026E+00  7.3442E+00
 EPSSHRINKSD(%)  3.2847E+01
 EPSSHRINKVR(%)  5.4904E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2040.5157590789313     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1121.5772258742586     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2040.516       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.38E+00  5.98E-01  7.96E-01  8.86E-01  9.87E-01  8.87E-01  5.13E-01  1.12E+00  8.82E-01  1.29E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       34.904
Stop Time:
Sun Oct 24 00:20:27 CDT 2021
