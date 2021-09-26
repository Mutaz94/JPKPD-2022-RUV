Sat Sep 25 01:36:53 CDT 2021
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
$DATA ../../../../data/int/SL2/dat80.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      996
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

 TOT. NO. OF OBS RECS:      896
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2125.83239426925        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5797E+00  6.8046E+01  1.1792E+02  1.1686E+01  4.7690E+01  4.2485E+01 -7.0362E+01 -1.2928E+02 -3.2930E+01 -2.4732E+00
            -3.2817E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3080.34142382231        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0263E+00  1.0638E+00  9.8820E-01  9.7406E-01  1.0385E+00  8.6472E-01  9.9090E-01  1.0035E+00  9.1715E-01  9.3548E-01
             2.0266E+00
 PARAMETER:  1.2600E-01  1.6189E-01  8.8131E-02  7.3713E-02  1.3774E-01 -4.5350E-02  9.0860E-02  1.0352E-01  1.3511E-02  3.3309E-02
             8.0637E-01
 GRADIENT:   6.0481E+00 -9.1394E+00 -5.7052E+00 -1.2034E+01  1.0856E-01 -9.0416E+00  2.1694E+00  1.0786E+00 -1.5261E+00 -3.5230E-02
            -6.1072E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3081.33065664749        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0244E+00  1.2265E+00  1.1236E+00  9.0882E-01  1.2062E+00  8.9534E-01  8.3111E-01  9.5093E-01  9.2130E-01  1.0435E+00
             2.0530E+00
 PARAMETER:  1.2410E-01  3.0419E-01  2.1652E-01  4.3940E-03  2.8747E-01 -1.0557E-02 -8.4995E-02  4.9686E-02  1.8030E-02  1.4262E-01
             8.1931E-01
 GRADIENT:  -8.2564E-01  2.2420E+01  4.0087E+00  2.7469E+01  1.1054E+01  4.3534E+00 -5.8730E+00 -5.3964E+00 -8.8360E+00 -5.4024E+00
            -4.9373E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3082.84657545596        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  1.0250E+00  1.2050E+00  1.1265E+00  9.0823E-01  1.1871E+00  8.8479E-01  8.6903E-01  1.1882E+00  9.5207E-01  1.0331E+00
             2.0836E+00
 PARAMETER:  1.2467E-01  2.8650E-01  2.1908E-01  3.7399E-03  2.7152E-01 -2.2405E-02 -4.0383E-02  2.7248E-01  5.0879E-02  1.3253E-01
             8.3412E-01
 GRADIENT:  -1.3719E+01 -4.4705E+00 -4.7442E-01  1.8097E+00 -1.6883E+00 -9.3442E-01 -4.4851E-01 -4.3085E-03 -2.2175E-01 -8.5048E-01
            -3.8329E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3083.81119944059        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  1.0294E+00  1.3978E+00  1.5240E+00  7.9964E-01  1.4488E+00  8.8492E-01  7.9235E-01  2.1112E+00  9.9541E-01  1.1361E+00
             2.0789E+00
 PARAMETER:  1.2899E-01  4.3491E-01  5.2135E-01 -1.2360E-01  4.7072E-01 -2.2261E-02 -1.3275E-01  8.4725E-01  9.5402E-02  2.2758E-01
             8.3186E-01
 GRADIENT:  -1.4176E+00 -2.8795E+00 -2.3870E+00 -3.7368E+00  3.6197E+00 -6.0621E-01 -3.8503E-02  8.5022E-01  4.9415E-01 -6.2378E-02
            -5.6722E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3085.14219253855        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  1.0305E+00  1.2453E+00  2.1103E+00  9.1238E-01  1.4372E+00  8.8548E-01  8.9540E-01  2.5760E+00  8.6384E-01  1.0435E+00
             2.0711E+00
 PARAMETER:  1.3004E-01  3.1935E-01  8.4681E-01  8.2984E-03  4.6270E-01 -2.1625E-02 -1.0484E-02  1.0462E+00 -4.6370E-02  1.4260E-01
             8.2807E-01
 GRADIENT:   2.1678E+00  5.3248E+00  1.4467E+00  2.7889E+00  3.9955E+00 -1.8028E-01  6.3697E-01 -4.9425E+00  3.5496E-01 -8.7365E-01
             6.6140E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3085.16085093612        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      795
 NPARAMETR:  1.0297E+00  1.2395E+00  2.1140E+00  9.1350E-01  1.4365E+00  8.8584E-01  8.9263E-01  2.5773E+00  8.6100E-01  1.0482E+00
             2.0712E+00
 PARAMETER:  1.2924E-01  3.1471E-01  8.4858E-01  9.5332E-03  4.6223E-01 -2.1221E-02 -1.3586E-02  1.0467E+00 -4.9657E-02  1.4712E-01
             8.2811E-01
 GRADIENT:   4.2190E-02  4.5635E-01  1.1772E+00 -7.5976E-01  4.8704E+00 -3.6373E-03  3.6194E-02 -4.8355E+00 -5.6068E-02 -1.6085E-02
             7.2322E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3085.16194900807        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      971
 NPARAMETR:  1.0297E+00  1.2372E+00  2.1140E+00  9.1520E-01  1.4365E+00  8.8585E-01  8.9344E-01  2.5773E+00  8.6026E-01  1.0474E+00
             2.0712E+00
 PARAMETER:  1.2923E-01  3.1282E-01  8.4858E-01  1.1389E-02  4.6223E-01 -2.1206E-02 -1.2682E-02  1.0467E+00 -5.0516E-02  1.4628E-01
             8.2811E-01
 GRADIENT:   2.2056E-02  6.6098E-03  7.8444E-01 -1.6000E-03  6.2180E+00 -1.6142E-03 -1.3367E-03 -4.8286E+00  3.4979E-03 -1.1729E-02
             7.3360E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3085.18898467007        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1146
 NPARAMETR:  1.0295E+00  1.2370E+00  2.1097E+00  9.1518E-01  1.4339E+00  8.8578E-01  8.9334E-01  2.5821E+00  8.6017E-01  1.0472E+00
             2.0690E+00
 PARAMETER:  1.2910E-01  3.1268E-01  8.4656E-01  1.1363E-02  4.6042E-01 -2.1283E-02 -1.2792E-02  1.0486E+00 -5.0622E-02  1.4613E-01
             8.2705E-01
 GRADIENT:  -2.6370E-01  4.8529E-01  8.6097E-01 -3.7007E-01  4.7054E+00 -4.2777E-02 -2.1980E-02 -4.6772E+00 -7.0796E-02  1.6732E-01
             5.3336E+00

0ITERATION NO.:   41    OBJECTIVE VALUE:  -3085.18898467007        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     1176
 NPARAMETR:  1.0295E+00  1.2370E+00  2.1096E+00  9.1518E-01  1.4340E+00  8.8580E-01  8.9337E-01  2.5819E+00  8.6026E-01  1.0472E+00
             2.0691E+00
 PARAMETER:  1.2910E-01  3.1268E-01  8.4656E-01  1.1363E-02  4.6042E-01 -2.1283E-02 -1.2792E-02  1.0486E+00 -5.0622E-02  1.4613E-01
             8.2705E-01
 GRADIENT:  -1.3999E+05 -1.1561E+05  2.1355E+04 -1.8074E+05 -7.8504E+04 -3.0663E-02 -1.9933E-02  3.4402E+04 -6.3228E-02 -2.4737E+05
            -4.3709E+04
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         2.8         2.6         3.3         2.2         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1176
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5862E-03 -2.1769E-02 -3.3102E-02  1.4480E-02 -3.7800E-02
 SE:             2.9616E-02  2.0950E-02  2.0571E-02  2.2968E-02  2.1906E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5729E-01  2.9875E-01  1.0759E-01  5.2841E-01  8.4426E-02

 ETASHRINKSD(%)  7.8178E-01  2.9815E+01  3.1084E+01  2.3053E+01  2.6613E+01
 ETASHRINKVR(%)  1.5575E+00  5.0741E+01  5.2506E+01  4.0792E+01  4.6144E+01
 EBVSHRINKSD(%)  1.1973E+00  3.0038E+01  3.6650E+01  2.5402E+01  2.3818E+01
 EBVSHRINKVR(%)  2.3802E+00  5.1054E+01  5.9867E+01  4.4352E+01  4.1962E+01
 RELATIVEINF(%)  9.7576E+01  8.8354E+00  2.1445E+01  1.0876E+01  2.4006E+01
 EPSSHRINKSD(%)  1.8592E+01
 EPSSHRINKVR(%)  3.3727E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          896
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1646.7378515027735     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3085.1889846700747     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1438.4511331673011     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.84
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3085.189       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.24E+00  2.11E+00  9.15E-01  1.43E+00  8.86E-01  8.93E-01  2.58E+00  8.60E-01  1.05E+00  2.07E+00
 


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
+        5.12E+08
 
 TH 2
+        1.76E+08  6.04E+07
 
 TH 3
+        1.77E+02 -2.20E+03  2.83E+06
 
 TH 4
+        4.82E+04  2.55E+08 -3.62E+03  1.08E+09
 
 TH 5
+       -4.79E+02  5.88E+03 -4.83E+02  9.80E+03  2.07E+07
 
 TH 6
+       -2.43E+03 -8.39E+02  1.81E+02 -3.54E+03 -4.91E+02  2.47E+02
 
 TH 7
+       -1.65E+04 -5.67E+03  1.23E+03 -2.40E+04 -3.33E+03  4.99E-01  8.23E+01
 
 TH 8
+        1.16E+02 -1.47E+03  8.65E+02 -2.37E+03 -2.93E+02  1.19E+02  8.13E+02  1.23E+06
 
 TH 9
+       -3.58E+03 -1.24E+03  2.68E+02 -5.17E+03 -7.13E+02  1.46E+01  2.97E+01  1.81E+02  9.81E+01
 
 TH10
+       -2.05E+03  2.59E+04 -3.48E+02  4.19E+04  9.05E+02 -2.11E+03 -1.44E+04 -2.25E+02 -3.11E+03  3.86E+08
 
 TH11
+       -1.98E+02  2.30E+03 -1.21E+03  3.74E+03  3.26E+03 -1.86E+02 -1.28E+03 -7.89E+02 -2.69E+02  3.70E+02  3.09E+06
 
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
 #CPUT: Total CPU Time in Seconds,       45.998
Stop Time:
Sat Sep 25 01:37:40 CDT 2021
