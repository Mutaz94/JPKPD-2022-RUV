Sun Oct 24 03:08:38 CDT 2021
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
$DATA ../../../../data/SD4/SL1/dat94.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1660.02200252862        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3035E+02 -2.8903E+01 -3.2660E+01 -1.3050E+01  5.5970E+01  3.5277E+01 -5.0859E+00  1.0272E+01 -3.3591E+01  2.6405E+01
             1.2659E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1669.75904372550        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.9553E-01  9.9074E-01  1.0690E+00  1.0569E+00  9.7870E-01  1.0491E+00  1.0250E+00  9.5203E-01  1.2035E+00  8.5221E-01
             9.5467E-01
 PARAMETER:  9.5516E-02  9.0698E-02  1.6675E-01  1.5534E-01  7.8465E-02  1.4789E-01  1.2473E-01  5.0845E-02  2.8521E-01 -5.9928E-02
             5.3614E-02
 GRADIENT:  -3.4702E-01 -3.4622E+00  2.0082E+00 -4.7543E+00  1.5332E+01  7.5409E+00  4.7256E+00  1.7337E+00  9.5796E+00  5.6278E+00
            -8.7211E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1672.03643213107        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  9.9810E-01  8.7327E-01  8.1845E-01  1.1303E+00  7.8797E-01  9.9722E-01  1.2748E+00  7.6636E-01  1.0749E+00  5.6999E-01
             1.0115E+00
 PARAMETER:  9.8101E-02 -3.5505E-02 -1.0034E-01  2.2252E-01 -1.3830E-01  9.7213E-02  3.4282E-01 -1.6610E-01  1.7220E-01 -4.6213E-01
             1.1148E-01
 GRADIENT:   4.9090E-01  1.3652E+01 -1.2493E+00  1.8618E+01 -6.8729E+00 -1.3261E+01  4.3623E+00  4.0163E+00  7.2766E+00 -6.6099E-01
             1.4916E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1673.55217908254        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.9801E-01  8.2794E-01  8.1617E-01  1.1461E+00  7.7827E-01  1.0328E+00  1.3003E+00  5.7370E-01  1.0319E+00  6.4402E-01
             9.6962E-01
 PARAMETER:  9.8003E-02 -8.8809E-02 -1.0313E-01  2.3638E-01 -1.5068E-01  1.3230E-01  3.6260E-01 -4.5564E-01  1.3141E-01 -3.4002E-01
             6.9144E-02
 GRADIENT:   1.5471E+00  4.9414E+00  5.2406E-01  5.6371E+00 -1.0866E+00  8.0491E-01  8.8197E-01  4.8212E-01  1.0206E+00  1.4289E-01
            -1.2516E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1673.64806169670        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      700
 NPARAMETR:  9.9701E-01  7.5606E-01  8.0958E-01  1.1843E+00  7.4944E-01  1.0307E+00  1.3783E+00  5.0416E-01  9.9492E-01  6.4756E-01
             9.7238E-01
 PARAMETER:  9.7006E-02 -1.7964E-01 -1.1124E-01  2.6917E-01 -1.8843E-01  1.3026E-01  4.2082E-01 -5.8486E-01  9.4911E-02 -3.3455E-01
             7.1992E-02
 GRADIENT:   1.2763E-01  2.0372E+00  1.9484E+00  1.3165E+00 -2.8040E+00  8.2588E-02 -4.0528E-01 -4.0084E-01 -5.6293E-01 -3.9177E-01
            -1.4780E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1673.65140300838        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      875
 NPARAMETR:  9.9653E-01  7.0439E-01  8.3859E-01  1.2200E+00  7.4607E-01  1.0300E+00  1.4427E+00  5.3966E-01  9.8501E-01  6.5603E-01
             9.7214E-01
 PARAMETER:  9.6521E-02 -2.5042E-01 -7.6035E-02  2.9886E-01 -1.9294E-01  1.2956E-01  4.6653E-01 -5.1682E-01  8.4893E-02 -3.2154E-01
             7.1748E-02
 GRADIENT:   6.1265E-01  3.0436E+00  2.6933E+00  3.4167E+00 -3.7047E+00  3.3676E-02  4.1272E-01 -6.1559E-01  1.0035E+00 -4.4949E-01
            -2.2796E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1673.65215296419        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1052
 NPARAMETR:  9.9603E-01  6.6742E-01  8.6110E-01  1.2455E+00  7.4463E-01  1.0294E+00  1.4881E+00  5.7542E-01  9.7713E-01  6.6079E-01
             9.7201E-01
 PARAMETER:  9.6024E-02 -3.0434E-01 -4.9546E-02  3.1952E-01 -1.9487E-01  1.2895E-01  4.9750E-01 -4.5266E-01  7.6864E-02 -3.1432E-01
             7.1615E-02
 GRADIENT:   7.8454E-01  3.6227E+00  3.0762E+00  4.3944E+00 -4.5030E+00 -7.6621E-03  8.2103E-01 -5.9924E-01  1.6499E+00 -3.3706E-01
            -2.0342E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1673.68742471565        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  9.9639E-01  6.6356E-01  8.6099E-01  1.2426E+00  7.4474E-01  1.0299E+00  1.4785E+00  5.9352E-01  9.7166E-01  6.5926E-01
             9.7190E-01
 PARAMETER:  9.6387E-02 -3.1014E-01 -4.9667E-02  3.1717E-01 -1.9473E-01  1.2945E-01  4.9105E-01 -4.2169E-01  7.1251E-02 -3.1663E-01
             7.1497E-02
 GRADIENT:   1.6037E+00 -2.0454E-01 -6.3659E-01 -2.7922E+00  1.8459E-01  2.3274E-01  8.7141E-02 -2.4893E-03  5.7124E-02  2.3228E-01
             5.1541E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1673.68810421189        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1403
 NPARAMETR:  9.9638E-01  6.6368E-01  8.6132E-01  1.2426E+00  7.4482E-01  1.0299E+00  1.4787E+00  5.9724E-01  9.7161E-01  6.5750E-01
             9.7194E-01
 PARAMETER:  9.6375E-02 -3.0995E-01 -4.9290E-02  3.1721E-01 -1.9462E-01  1.2944E-01  4.9118E-01 -4.1544E-01  7.1194E-02 -3.1931E-01
             7.1536E-02
 GRADIENT:   1.5842E+00 -1.2813E-01 -6.3050E-01 -2.6452E+00  1.6148E-01  2.3118E-01  9.7332E-02  2.8407E-02  2.6889E-02  1.4889E-01
             5.1431E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1403
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.4740E-04  8.1820E-03 -2.2508E-02 -7.0578E-03 -1.4778E-02
 SE:             2.9846E-02  1.8252E-02  1.2879E-02  2.6548E-02  2.0237E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9071E-01  6.5396E-01  8.0521E-02  7.9036E-01  4.6523E-01

 ETASHRINKSD(%)  1.3585E-02  3.8852E+01  5.6855E+01  1.1061E+01  3.2203E+01
 ETASHRINKVR(%)  2.7168E-02  6.2609E+01  8.1385E+01  2.0899E+01  5.4036E+01
 EBVSHRINKSD(%)  4.0783E-01  4.0973E+01  5.7969E+01  1.0491E+01  3.0793E+01
 EBVSHRINKVR(%)  8.1400E-01  6.5159E+01  8.2334E+01  1.9882E+01  5.2104E+01
 RELATIVEINF(%)  9.8537E+01  2.3201E+00  2.0462E+00  8.8482E+00  3.6818E+00
 EPSSHRINKSD(%)  4.4722E+01
 EPSSHRINKVR(%)  6.9444E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1673.6881042118885     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -938.53727764815028     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1673.688       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  6.64E-01  8.61E-01  1.24E+00  7.45E-01  1.03E+00  1.48E+00  5.97E-01  9.72E-01  6.58E-01  9.72E-01
 


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
 #CPUT: Total CPU Time in Seconds,       41.443
Stop Time:
Sun Oct 24 03:08:47 CDT 2021
