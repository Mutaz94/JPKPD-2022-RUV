Sat Sep 18 09:51:36 CDT 2021
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
$DATA ../../../../data/spa/A2/dat43.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1177.90169368204        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.9344E+01  1.4927E+01  5.1682E+01 -2.3288E+01  9.1870E+01  2.4010E+01 -3.6865E+00 -1.1655E+01  1.7211E+00 -7.4812E+01
            -8.3857E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1443.52026970061        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0073E+00  9.3457E-01  9.0662E-01  1.0429E+00  8.7923E-01  8.7578E-01  8.8737E-01  9.6333E-01  8.5044E-01  1.0160E+00
             2.3276E+00
 PARAMETER:  1.0729E-01  3.2335E-02  1.9726E-03  1.4197E-01 -2.8707E-02 -3.2642E-02 -1.9489E-02  6.2637E-02 -6.1998E-02  1.1587E-01
             9.4485E-01
 GRADIENT:   6.0765E+01 -1.6212E+01  5.7213E+00 -4.6553E+01 -1.8230E+00 -1.7492E+01  6.3209E+00  5.2786E+00  1.1610E-01  9.7759E+00
             2.0820E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1455.22570273424        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  1.0106E+00  6.7800E-01  3.3921E-01  1.1651E+00  4.2542E-01  9.5367E-01  7.1782E-01  2.4809E-01  8.2977E-01  4.5226E-01
             2.2337E+00
 PARAMETER:  1.1055E-01 -2.8861E-01 -9.8114E-01  2.5282E-01 -7.5467E-01  5.2559E-02 -2.3153E-01 -1.2940E+00 -8.6608E-02 -6.9350E-01
             9.0368E-01
 GRADIENT:   3.9455E+01  5.4937E+01 -1.2305E+01  1.7025E+02  7.6360E+00  6.6148E+00 -1.3978E+01 -2.2506E-01 -4.0711E+00 -1.0434E+01
             3.1249E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1467.51788428519        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  9.7251E-01  6.6139E-01  2.2104E-01  1.0468E+00  3.3936E-01  9.3200E-01  1.1195E+00  2.5927E-02  8.7446E-01  3.8436E-01
             2.0100E+00
 PARAMETER:  7.2130E-02 -3.1341E-01 -1.4094E+00  1.4569E-01 -9.8068E-01  2.9574E-02  2.1286E-01 -3.5525E+00 -3.4154E-02 -8.5618E-01
             7.9811E-01
 GRADIENT:  -2.8224E+01  6.4274E+01 -2.3515E+01  1.0040E+02 -1.4766E+00 -5.6020E+00  1.2987E+01 -4.4153E-04  6.9754E+00  4.0928E+00
             1.1348E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1476.88328280533        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.8076E-01  4.6032E-01  2.0683E-01  1.0388E+00  2.8131E-01  9.5853E-01  1.2509E+00  1.0000E-02  8.3822E-01  3.4841E-01
             1.8923E+00
 PARAMETER:  8.0574E-02 -6.7583E-01 -1.4758E+00  1.3807E-01 -1.1683E+00  5.7650E-02  3.2385E-01 -5.2719E+00 -7.6475E-02 -9.5438E-01
             7.3778E-01
 GRADIENT:  -1.9190E+01  4.3049E-01  3.7417E+00 -2.4595E+01 -3.1174E+00  1.9767E+00 -4.8646E+00  0.0000E+00 -1.2870E+00 -3.2428E+00
            -1.0344E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1478.83145483052        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  9.9163E-01  3.5696E-01  2.6244E-01  1.1474E+00  2.9171E-01  9.4634E-01  1.6596E+00  1.0000E-02  7.8588E-01  4.6843E-01
             1.9549E+00
 PARAMETER:  9.1595E-02 -9.3012E-01 -1.2377E+00  2.3751E-01 -1.1320E+00  4.4850E-02  6.0659E-01 -6.1988E+00 -1.4095E-01 -6.5838E-01
             7.7034E-01
 GRADIENT:  -6.4074E+00  7.9148E+00  1.5011E+01 -2.0147E+00 -2.7898E+01 -5.1507E-01  7.1269E-01  0.0000E+00  8.2998E-01 -5.3100E-01
             9.0564E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1480.41699074230        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.8647E-01  1.8885E-01  3.5773E-01  1.2923E+00  3.2902E-01  9.3248E-01  2.6948E+00  1.0000E-02  7.2726E-01  6.0373E-01
             2.0131E+00
 PARAMETER:  8.6375E-02 -1.5668E+00 -9.2797E-01  3.5642E-01 -1.0116E+00  3.0089E-02  1.0913E+00 -1.0401E+01 -2.1847E-01 -4.0462E-01
             7.9967E-01
 GRADIENT:  -2.8862E+00 -1.9346E+00  1.8207E+01 -1.9542E+00 -1.6953E+01 -2.1832E+00 -7.3164E+00  0.0000E+00  4.5252E+00  3.0152E+00
             5.4563E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1481.14525132971        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  9.8533E-01  1.4801E-01  3.4675E-01  1.2986E+00  3.1925E-01  9.3852E-01  3.3360E+00  1.0000E-02  7.1334E-01  5.9195E-01
             1.9930E+00
 PARAMETER:  8.5221E-02 -1.8105E+00 -9.5916E-01  3.6131E-01 -1.0418E+00  3.6544E-02  1.3048E+00 -1.2489E+01 -2.3780E-01 -4.2434E-01
             7.8966E-01
 GRADIENT:  -1.6075E-01  1.7057E-01  3.7981E-01  1.8867E-01 -6.3530E-01  1.8328E-02  4.3581E-02  0.0000E+00 -8.2695E-02  1.7719E-01
            -1.0474E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1481.14603876344        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      935
 NPARAMETR:  9.8530E-01  1.4640E-01  3.4641E-01  1.2988E+00  3.1896E-01  9.3842E-01  3.3587E+00  1.0000E-02  7.1338E-01  5.9076E-01
             1.9934E+00
 PARAMETER:  8.5187E-02 -1.8214E+00 -9.6013E-01  3.6147E-01 -1.0427E+00  3.6447E-02  1.3116E+00 -1.2581E+01 -2.3774E-01 -4.2634E-01
             7.8986E-01
 GRADIENT:  -2.4316E-02  5.1983E-02 -2.0366E-01  8.1707E-02  1.5064E-01 -1.8042E-02  8.1735E-02  0.0000E+00 -5.0427E-03 -4.2609E-02
             3.3547E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      935
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4298E-03  4.2800E-02 -2.3072E-04 -2.2801E-02  1.5909E-02
 SE:             2.9360E-02  1.5403E-02  2.8139E-04  2.5970E-02  2.0323E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6116E-01  5.4574E-03  4.1224E-01  3.7996E-01  4.3376E-01

 ETASHRINKSD(%)  1.6403E+00  4.8399E+01  9.9057E+01  1.2997E+01  3.1914E+01
 ETASHRINKVR(%)  3.2537E+00  7.3373E+01  9.9991E+01  2.4305E+01  5.3643E+01
 EBVSHRINKSD(%)  1.6604E+00  6.0701E+01  9.9052E+01  1.1743E+01  2.6588E+01
 EBVSHRINKVR(%)  3.2932E+00  8.4556E+01  9.9991E+01  2.2107E+01  4.6106E+01
 RELATIVEINF(%)  9.5194E+01  6.5554E+00  3.3384E-04  2.1350E+01  2.0559E+00
 EPSSHRINKSD(%)  3.6794E+01
 EPSSHRINKVR(%)  6.0050E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1481.1460387634377     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -745.99521219969949     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.69
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1481.146       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  1.46E-01  3.46E-01  1.30E+00  3.19E-01  9.38E-01  3.36E+00  1.00E-02  7.13E-01  5.91E-01  1.99E+00
 


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
+        1.26E+03
 
 TH 2
+       -4.73E+01  2.32E+03
 
 TH 3
+       -1.72E+01 -2.60E+02  9.56E+03
 
 TH 4
+       -4.37E+01  8.23E+01 -6.32E+02  1.07E+03
 
 TH 5
+        1.42E+02 -3.02E+02 -1.26E+04 -3.40E+02  1.90E+04
 
 TH 6
+       -5.26E+00 -4.29E+00  6.23E+00 -8.46E+00  7.64E+00  2.09E+02
 
 TH 7
+        2.83E+00  1.27E+02 -6.91E+01 -1.74E+01  8.59E+01  3.11E-02  8.75E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.49E+00 -5.75E+01  1.09E+02 -1.68E+01  3.79E+01  1.08E+00 -8.95E-01  0.00E+00  2.65E+02
 
 TH10
+       -1.21E+01 -2.20E+02 -1.18E+01  3.74E+01 -1.11E+02 -1.77E+00 -1.42E+01  0.00E+00 -1.13E+01  2.31E+02
 
 TH11
+       -1.26E+01 -3.19E+01 -2.08E+01 -9.54E+00  1.09E+01  2.37E+00 -1.52E+00  0.00E+00  1.36E+01  3.60E+01  6.55E+01
 
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
 #CPUT: Total CPU Time in Seconds,       17.342
Stop Time:
Sat Sep 18 09:51:55 CDT 2021
