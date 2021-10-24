Sat Oct 23 22:52:17 CDT 2021
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
$DATA ../../../../data/SD3/S1/dat3.csv ignore=@
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
Current Date:       23 OCT 2021
Days until program expires : 176
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2143.52240034256        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0314E+02 -2.9056E+01 -3.0867E+01 -4.6173E+00  3.8619E+01  6.9165E+01  5.3044E+00  1.7706E+01  5.3991E+00  1.7859E+01
            -2.6196E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2150.63661350165        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0270E+00  9.8642E-01  1.1482E+00  1.0853E+00  1.0321E+00  1.0208E+00  9.6171E-01  8.5030E-01  9.9191E-01  8.5778E-01
             1.0209E+00
 PARAMETER:  1.2661E-01  8.6324E-02  2.3820E-01  1.8182E-01  1.3164E-01  1.2056E-01  6.0959E-02 -6.2168E-02  9.1874E-02 -5.3408E-02
             1.2071E-01
 GRADIENT:   1.0615E+01  7.6675E+00  4.3298E+00  7.6901E+00  3.2080E+01  3.5756E+00 -1.2092E+00  1.6985E+00 -9.4295E-01 -1.6650E+01
            -2.2515E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2154.12904438484        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0376E+00  8.7355E-01  9.4442E-01  1.1466E+00  9.0153E-01  9.6136E-01  1.0436E+00  3.4886E-01  9.4308E-01  9.4365E-01
             1.0159E+00
 PARAMETER:  1.3691E-01 -3.5187E-02  4.2821E-02  2.3681E-01 -3.6658E-03  6.0594E-02  1.4267E-01 -9.5308E-01  4.1399E-02  4.1998E-02
             1.1575E-01
 GRADIENT:   3.7032E+01  5.9722E-01 -2.7809E+01  2.5593E+01  3.8903E+01 -2.1400E+01 -1.8855E+00  8.8424E-01 -3.5010E-01  5.5170E+00
            -1.1319E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2156.70431548140        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  1.0185E+00  7.0232E-01  9.4142E-01  1.2382E+00  8.0881E-01  1.0161E+00  1.2748E+00  3.0294E-01  8.7834E-01  8.7641E-01
             1.0303E+00
 PARAMETER:  1.1834E-01 -2.5337E-01  3.9636E-02  3.1362E-01 -1.1219E-01  1.1598E-01  3.4281E-01 -1.0942E+00 -2.9725E-02 -3.1923E-02
             1.2982E-01
 GRADIENT:  -3.4246E+00  7.4311E+00  4.1571E+00  5.9868E+00 -9.3715E+00  2.9569E+00 -8.3970E-02  1.7616E-02 -6.5829E-01  9.4904E-01
             1.7402E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2156.85030644489        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  1.0198E+00  6.4032E-01  9.4358E-01  1.2682E+00  7.9295E-01  1.0079E+00  1.3485E+00  2.9950E-01  8.6491E-01  8.7677E-01
             1.0272E+00
 PARAMETER:  1.1963E-01 -3.4578E-01  4.1930E-02  3.3763E-01 -1.3199E-01  1.0791E-01  3.9896E-01 -1.1057E+00 -4.5126E-02 -3.1506E-02
             1.2687E-01
 GRADIENT:   1.3207E+00  3.3003E-01 -2.4792E+00 -2.4321E+00 -6.1767E-02  1.9336E-01 -8.3501E-02 -2.7628E-02  2.6380E-01  7.2466E-01
             3.0603E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2156.86716708215        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      910
 NPARAMETR:  1.0162E+00  6.2290E-01  9.5852E-01  1.2780E+00  7.9410E-01  1.0053E+00  1.3779E+00  3.4244E-01  8.5841E-01  8.6960E-01
             1.0271E+00
 PARAMETER:  1.1607E-01 -3.7337E-01  5.7631E-02  3.4532E-01 -1.3054E-01  1.0526E-01  4.2058E-01 -9.7165E-01 -5.2671E-02 -3.9722E-02
             1.2672E-01
 GRADIENT:  -5.8463E+00 -9.5065E-02 -1.9509E-01 -6.2803E+00 -1.2405E+00 -7.4550E-01  1.3080E-02 -2.8945E-02 -2.3546E-01 -7.7882E-01
            -1.2966E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2156.88976505992        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1068
 NPARAMETR:  1.0205E+00  6.2129E-01  9.6056E-01  1.2805E+00  7.9482E-01  1.0081E+00  1.3772E+00  3.4249E-01  8.5897E-01  8.7530E-01
             1.0271E+00
 PARAMETER:  1.2028E-01 -3.7595E-01  5.9759E-02  3.4728E-01 -1.2964E-01  1.0807E-01  4.2002E-01 -9.7153E-01 -5.2026E-02 -3.3192E-02
             1.2670E-01
 GRADIENT:   3.6066E+00  4.1833E-01 -8.3100E-01 -3.7657E+00 -9.3365E-01  4.2122E-01  9.9967E-02 -5.3361E-03  5.2105E-02  8.3880E-03
            -1.3102E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -2156.88976505992        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1090
 NPARAMETR:  1.0205E+00  6.2129E-01  9.6056E-01  1.2805E+00  7.9482E-01  1.0081E+00  1.3772E+00  3.4249E-01  8.5897E-01  8.7530E-01
             1.0271E+00
 PARAMETER:  1.2028E-01 -3.7595E-01  5.9759E-02  3.4728E-01 -1.2964E-01  1.0807E-01  4.2002E-01 -9.7153E-01 -5.2026E-02 -3.3192E-02
             1.2670E-01
 GRADIENT:   3.6066E+00  4.1833E-01 -8.3100E-01 -3.7657E+00 -9.3365E-01  4.2122E-01  9.9967E-02 -5.3361E-03  5.2105E-02  8.3880E-03
            -1.3102E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1090
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.6765E-04  1.1940E-03 -1.2170E-02 -3.6506E-03 -1.4239E-02
 SE:             2.9872E-02  1.5495E-02  7.3746E-03  2.6999E-02  2.4024E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9552E-01  9.3858E-01  9.8893E-02  8.9245E-01  5.5340E-01

 ETASHRINKSD(%)  1.0000E-10  4.8090E+01  7.5294E+01  9.5488E+00  1.9515E+01
 ETASHRINKVR(%)  1.0000E-10  7.3054E+01  9.3896E+01  1.8186E+01  3.5222E+01
 EBVSHRINKSD(%)  3.4669E-01  4.9503E+01  7.6541E+01  9.4110E+00  1.7735E+01
 EBVSHRINKVR(%)  6.9218E-01  7.4501E+01  9.4497E+01  1.7936E+01  3.2325E+01
 RELATIVEINF(%)  9.8117E+01  1.6940E+00  8.4144E-01  7.5112E+00  8.2872E+00
 EPSSHRINKSD(%)  3.3040E+01
 EPSSHRINKVR(%)  5.5164E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2156.8897650599215     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1237.9512318552488     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2156.890       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  6.21E-01  9.61E-01  1.28E+00  7.95E-01  1.01E+00  1.38E+00  3.42E-01  8.59E-01  8.75E-01  1.03E+00
 


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
 #CPUT: Total CPU Time in Seconds,       87.792
Stop Time:
Sat Oct 23 22:52:32 CDT 2021
