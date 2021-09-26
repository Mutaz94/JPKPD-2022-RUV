Sat Sep 25 07:58:57 CDT 2021
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
$DATA ../../../../data/spa/A1/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1136.81121757569        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.8796E+01 -2.3317E+01 -1.3904E+01 -3.1768E+01  1.0221E+02 -1.4957E+01 -3.6178E+01  4.2315E+00 -4.8684E+01 -5.4492E+01
            -8.8867E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1413.84869538272        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0508E+00  9.6120E-01  1.2945E+00  1.0661E+00  1.0604E+00  9.7312E-01  1.1262E+00  8.2538E-01  1.2079E+00  9.1031E-01
             2.4165E+00
 PARAMETER:  1.4951E-01  6.0426E-02  3.5816E-01  1.6399E-01  1.5867E-01  7.2747E-02  2.1887E-01 -9.1910E-02  2.8891E-01  6.0312E-03
             9.8232E-01
 GRADIENT:   1.3455E+02 -1.8382E+01 -1.1953E+01 -1.4020E+01  3.0680E+01 -2.1836E+01  6.6655E+00  3.7120E+00  1.7964E+01 -8.1540E+00
             2.7212E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1418.82742184954        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0400E+00  7.2185E-01  1.8651E+00  1.2710E+00  1.1215E+00  1.0113E+00  1.2514E+00  3.3958E-01  1.0535E+00  1.3764E+00
             2.3015E+00
 PARAMETER:  1.3921E-01 -2.2594E-01  7.2330E-01  3.3984E-01  2.1465E-01  1.1126E-01  3.2430E-01 -9.8004E-01  1.5209E-01  4.1949E-01
             9.3358E-01
 GRADIENT:   1.0914E+02  1.2176E+01  5.9535E+00  2.5970E+01 -1.3573E+01 -5.2740E+00  2.8982E+00  3.3119E-01  1.0156E+01  1.7541E+01
             2.4080E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1425.65572624973        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.7789E-01  5.8509E-01  1.2263E+00  1.2972E+00  8.6398E-01  1.0264E+00  1.3340E+00  1.3184E-01  9.9063E-01  9.8519E-01
             2.1808E+00
 PARAMETER:  7.7643E-02 -4.3598E-01  3.0400E-01  3.6022E-01 -4.6203E-02  1.2602E-01  3.8821E-01 -1.9261E+00  9.0590E-02  8.5080E-02
             8.7971E-01
 GRADIENT:  -1.5697E+01  6.0641E+00  7.2821E+00 -1.4351E+00 -1.0523E+01  3.3153E+00 -1.7016E-01  1.4654E-01  3.9486E-02 -2.7836E+00
            -3.9340E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1428.64065889154        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.8315E-01  2.9037E-01  1.0405E+00  1.4571E+00  7.1202E-01  1.0150E+00  2.1001E+00  1.6133E-02  8.9835E-01  9.2939E-01
             2.1559E+00
 PARAMETER:  8.3010E-02 -1.1366E+00  1.3973E-01  4.7642E-01 -2.3966E-01  1.1492E-01  8.4197E-01 -4.0269E+00 -7.1952E-03  2.6773E-02
             8.6819E-01
 GRADIENT:   1.2979E+00  3.0562E+00  5.1899E-01  6.9564E+00 -8.3725E-01 -3.8405E-01 -2.5854E-01  3.3454E-03 -1.4436E+00 -1.7704E+00
            -1.3090E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1430.12744361561        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  9.7854E-01  8.3931E-02  8.6334E-01  1.5326E+00  5.8631E-01  1.0133E+00  4.0086E+00  1.0000E-02  8.6485E-01  8.9065E-01
             2.1276E+00
 PARAMETER:  7.8303E-02 -2.3778E+00 -4.6945E-02  5.2698E-01 -4.3390E-01  1.1320E-01  1.4884E+00 -8.5333E+00 -4.5198E-02 -1.5803E-02
             8.5497E-01
 GRADIENT:  -3.9016E-01  6.2674E-01 -1.4120E+00 -6.1511E+00  1.2557E+00 -6.4113E-01  3.9737E-01  0.0000E+00 -8.2563E-01 -2.5802E-01
             2.5977E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1431.01806566109        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  9.7601E-01  1.6016E-02  8.9371E-01  1.5834E+00  5.8619E-01  1.0130E+00  8.2045E+00  1.0000E-02  8.6087E-01  9.0432E-01
             2.1174E+00
 PARAMETER:  7.5722E-02 -4.0342E+00 -1.2372E-02  5.5957E-01 -4.3411E-01  1.1289E-01  2.2047E+00 -1.4519E+01 -4.9812E-02 -5.7504E-04
             8.5020E-01
 GRADIENT:  -1.8486E+00  1.5306E-01  2.2340E+00  1.3404E+01 -4.4063E+00 -5.6983E-01 -7.2984E-02  0.0000E+00  2.4843E+00 -1.7161E-01
            -1.4820E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1431.27754714620        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      560
 NPARAMETR:  9.7941E-01  1.0000E-02  9.4598E-01  1.5986E+00  6.0947E-01  1.0164E+00  1.1145E+01  1.0000E-02  8.5038E-01  9.2325E-01
             2.1290E+00
 PARAMETER:  7.9198E-02 -4.6626E+00  4.4464E-02  5.6912E-01 -3.9517E-01  1.1626E-01  2.5110E+00 -1.6768E+01 -6.2077E-02  2.0150E-02
             8.5563E-01
 GRADIENT:  -3.1895E+00  0.0000E+00  1.3821E+00 -5.5235E+00 -3.8553E+00 -8.1786E-02 -3.8305E-02  0.0000E+00 -2.5227E-01  2.7060E-01
            -2.9571E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1431.36277876512        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      735
 NPARAMETR:  9.8092E-01  1.0000E-02  1.0255E+00  1.6151E+00  6.4590E-01  1.0163E+00  1.1621E+01  1.0000E-02  8.4762E-01  9.4753E-01
             2.1365E+00
 PARAMETER:  8.0733E-02 -4.7289E+00  1.2515E-01  5.7942E-01 -3.3711E-01  1.1614E-01  2.5528E+00 -1.6874E+01 -6.5323E-02  4.6108E-02
             8.5917E-01
 GRADIENT:  -1.3912E-03  0.0000E+00 -3.9150E-03 -4.9692E-02  1.6973E-02  3.4040E-03 -1.6211E-03  0.0000E+00 -1.5936E-02 -6.0464E-03
            -2.4952E-03

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1431.36277876512        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      757
 NPARAMETR:  9.8092E-01  1.0000E-02  1.0255E+00  1.6151E+00  6.4590E-01  1.0163E+00  1.1621E+01  1.0000E-02  8.4762E-01  9.4753E-01
             2.1365E+00
 PARAMETER:  8.0733E-02 -4.7289E+00  1.2515E-01  5.7942E-01 -3.3711E-01  1.1614E-01  2.5528E+00 -1.6874E+01 -6.5323E-02  4.6108E-02
             8.5917E-01
 GRADIENT:  -1.3912E-03  0.0000E+00 -3.9150E-03 -4.9692E-02  1.6973E-02  3.4040E-03 -1.6211E-03  0.0000E+00 -1.5936E-02 -6.0464E-03
            -2.4952E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      757
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.8869E-04 -1.3124E-04  9.5518E-06 -1.0608E-02 -2.3529E-02
 SE:             2.9340E-02  1.7813E-03  1.3490E-04  2.7833E-02  2.1907E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9487E-01  9.4126E-01  9.4355E-01  7.0311E-01  2.8281E-01

 ETASHRINKSD(%)  1.7074E+00  9.4033E+01  9.9548E+01  6.7557E+00  2.6608E+01
 ETASHRINKVR(%)  3.3857E+00  9.9644E+01  9.9998E+01  1.3055E+01  4.6137E+01
 EBVSHRINKSD(%)  1.6933E+00  9.4459E+01  9.9475E+01  6.4943E+00  2.5795E+01
 EBVSHRINKVR(%)  3.3579E+00  9.9693E+01  9.9997E+01  1.2567E+01  4.4936E+01
 RELATIVEINF(%)  9.0173E+01  9.7768E-03  1.7656E-04  3.9634E+00  2.6759E+00
 EPSSHRINKSD(%)  3.4926E+01
 EPSSHRINKVR(%)  5.7654E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1431.3627787651214     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -696.21195220138327     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1431.363       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.00E-02  1.03E+00  1.62E+00  6.46E-01  1.02E+00  1.16E+01  1.00E-02  8.48E-01  9.48E-01  2.14E+00
 


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
+        1.07E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.02E+00  0.00E+00  2.80E+02
 
 TH 4
+       -2.74E+01  0.00E+00 -4.09E+01  5.19E+02
 
 TH 5
+        1.58E+01  0.00E+00 -6.46E+02 -8.60E+01  1.61E+03
 
 TH 6
+        3.30E+00  0.00E+00  1.58E+01 -8.24E+00 -6.29E+00  1.63E+02
 
 TH 7
+        1.48E-02  0.00E+00  6.53E-03 -1.81E-02  5.93E-03  3.95E-02  4.32E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.65E+01  0.00E+00  1.14E+01 -8.23E+00 -6.90E+00 -6.05E+00 -7.73E-03  0.00E+00  2.08E+02
 
 TH10
+        1.07E+00  0.00E+00 -3.90E+00 -1.58E+00 -7.25E+01  4.24E+00 -3.60E-03  0.00E+00  1.10E+01  6.15E+01
 
 TH11
+       -1.24E+01  0.00E+00 -8.56E+00 -8.73E+00 -1.88E-01  2.81E+00  1.59E-03  0.00E+00  8.57E+00  1.93E+01  5.99E+01
 
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
 #CPUT: Total CPU Time in Seconds,       13.347
Stop Time:
Sat Sep 25 07:59:12 CDT 2021
