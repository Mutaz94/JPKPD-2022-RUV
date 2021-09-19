Sat Sep 18 05:49:46 CDT 2021
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
$DATA ../../../../data/int/TD2/dat20.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3304.19006284087        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1096E+02 -1.6004E+01  7.2231E+01 -1.4454E+01  7.1790E+01  1.8275E+00  6.6160E+00 -2.4762E+02 -5.0008E+01  2.6233E+00
            -8.0430E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3732.36431964960        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  8.7912E-01  9.6966E-01  1.0314E+00  1.0393E+00  9.0213E-01  1.0448E+00  8.6420E-01  8.8714E-01  9.9125E-01  9.9232E-01
             1.4742E+00
 PARAMETER: -2.8835E-02  6.9186E-02  1.3087E-01  1.3857E-01 -2.9990E-03  1.4378E-01 -4.5950E-02 -1.9748E-02  9.1211E-02  9.2288E-02
             4.8808E-01
 GRADIENT:  -1.8656E+02  3.2001E+01  4.8356E+01  4.4936E+01 -7.0903E+01 -2.0469E+00  1.0237E+00  7.0097E+00 -4.9617E+00  7.3794E+00
             5.9756E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3792.48127502770        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.1088E-01  1.1857E+00  1.0105E+00  8.9423E-01  1.1205E+00  1.0111E+00  1.1923E+00  2.5697E-01  1.2300E+00  1.1234E+00
             1.2737E+00
 PARAMETER:  6.6583E-03  2.7035E-01  1.1040E-01 -1.1797E-02  2.1379E-01  1.1102E-01  2.7585E-01 -1.2588E+00  3.0704E-01  2.1634E-01
             3.4192E-01
 GRADIENT:  -1.1265E+02  5.6427E+01  3.2573E+00 -2.3056E+00  6.9217E+01 -3.4318E+00  4.5507E+01 -2.9067E+00  4.3591E+01  8.8056E+00
             3.6676E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3871.04993066433        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.4203E-01  9.7086E-01  9.5587E-01  9.9495E-01  9.2605E-01  9.4392E-01  1.0587E+00  6.9854E-01  9.7910E-01  1.0261E+00
             1.0402E+00
 PARAMETER:  4.0282E-02  7.0424E-02  5.4865E-02  9.4938E-02  2.3174E-02  4.2283E-02  1.5708E-01 -2.5876E-01  7.8878E-02  1.2573E-01
             1.3944E-01
 GRADIENT:  -3.2221E+01  5.3499E+00  1.7168E+01 -8.1672E+00 -2.2653E+00 -2.4349E+01  1.8747E+01 -1.0250E+01 -1.0874E+01  1.4680E+01
             9.0553E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3871.67436295052        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:      342
 NPARAMETR:  9.4217E-01  9.7071E-01  9.5680E-01  9.9520E-01  9.2604E-01  9.4512E-01  1.0548E+00  7.1883E-01  9.7959E-01  1.0234E+00
             1.0381E+00
 PARAMETER:  4.0426E-02  7.0275E-02  5.5841E-02  9.5190E-02  2.3164E-02  4.3552E-02  1.5335E-01 -2.3014E-01  7.9375E-02  1.2311E-01
             1.3735E-01
 GRADIENT:  -7.5502E+01 -2.4643E+00  1.5785E+01 -1.7876E+01 -7.7332E+00 -3.2235E+01  1.7152E+01 -9.9058E+00 -1.1886E+01  1.3839E+01
             8.8143E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3873.13379197039        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:      476
 NPARAMETR:  9.4216E-01  9.7071E-01  9.4762E-01  9.9520E-01  9.2521E-01  9.7422E-01  1.0548E+00  7.1882E-01  9.9883E-01  9.8439E-01
             1.0381E+00
 PARAMETER:  4.0421E-02  7.0270E-02  4.6200E-02  9.5186E-02  2.2265E-02  7.3877E-02  1.5334E-01 -2.3015E-01  9.8830E-02  8.4269E-02
             1.3736E-01
 GRADIENT:  -7.1016E+01  5.3398E-01  8.5534E+00 -1.7095E+01 -5.8805E-01 -1.8792E+01  1.6989E+01 -9.6311E+00 -6.6576E+00  6.4698E+00
             8.8606E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3873.39341568971        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      653
 NPARAMETR:  9.4216E-01  9.7071E-01  9.4593E-01  9.9520E-01  9.2405E-01  9.8130E-01  1.0548E+00  7.1882E-01  1.0025E+00  9.8069E-01
             1.0372E+00
 PARAMETER:  4.0422E-02  7.0270E-02  4.4415E-02  9.5186E-02  2.1006E-02  8.1119E-02  1.5334E-01 -2.3015E-01  1.0247E-01  8.0496E-02
             1.3650E-01
 GRADIENT:  -2.5876E+01  9.4593E+00  8.5818E+00 -7.4053E+00  3.9608E+00 -7.1591E+00  1.8006E+01 -9.5134E+00 -4.3293E+00  6.2953E+00
             8.7329E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3873.72341433781        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      842
 NPARAMETR:  9.4263E-01  9.7066E-01  9.4596E-01  9.9517E-01  9.2403E-01  1.0194E+00  1.0548E+00  7.1878E-01  1.0024E+00  9.8023E-01
             1.0372E+00
 PARAMETER:  4.0922E-02  7.0225E-02  4.4440E-02  9.5161E-02  2.0992E-02  1.1923E-01  1.5338E-01 -2.3021E-01  1.0244E-01  8.0036E-02
             1.3650E-01
 GRADIENT:  -1.9614E+01  9.3375E+00  8.6786E+00 -7.5819E+00  3.9014E+00  1.0156E+01  1.7993E+01 -9.5267E+00 -4.3208E+00  6.2251E+00
             8.7365E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3873.73972245147        NO. OF FUNC. EVALS.: 221
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  9.4268E-01  9.7066E-01  9.4595E-01  9.9518E-01  9.2403E-01  1.0194E+00  1.0546E+00  7.1894E-01  1.0025E+00  9.8019E-01
             1.0371E+00
 PARAMETER:  4.0972E-02  7.0223E-02  4.4436E-02  9.5165E-02  2.0993E-02  1.1923E-01  1.5315E-01 -2.2998E-01  1.0248E-01  7.9988E-02
             1.3643E-01
 GRADIENT:  -6.3596E+01  1.4908E+00  7.9536E+00 -1.7418E+01 -1.2231E+00 -6.2575E-03  1.6957E+01 -9.6058E+00 -5.6714E+00  5.7692E+00
             8.6996E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3873.96608309866        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1229
 NPARAMETR:  9.4270E-01  9.7066E-01  9.4595E-01  9.9518E-01  9.2399E-01  1.0190E+00  1.0545E+00  7.1911E-01  1.0034E+00  9.8017E-01
             1.0344E+00
 PARAMETER:  4.0972E-02  7.0222E-02  4.4436E-02  9.5165E-02  2.0941E-02  1.1884E-01  1.5315E-01 -2.2998E-01  1.0335E-01  7.9988E-02
             1.3380E-01
 GRADIENT:  -6.3572E+01  1.6442E+00  8.1635E+00 -1.7295E+01 -1.2088E+00 -1.5884E-01  1.6896E+01 -9.7851E+00 -5.4489E+00  5.7233E+00
            -2.3544E+06

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1229
 NO. OF SIG. DIGITS UNREPORTABLE

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9499E-02 -2.5334E-02 -1.9168E-02  2.0616E-02 -2.3294E-02
 SE:             2.9759E-02  2.1499E-02  1.6421E-02  2.8246E-02  2.4277E-02
 N:                     100         100         100         100         100

 P VAL.:         3.2155E-01  2.3865E-01  2.4309E-01  4.6545E-01  3.3730E-01

 ETASHRINKSD(%)  3.0523E-01  2.7976E+01  4.4989E+01  5.3733E+00  1.8669E+01
 ETASHRINKVR(%)  6.0953E-01  4.8125E+01  6.9738E+01  1.0458E+01  3.3853E+01
 EBVSHRINKSD(%)  2.7017E-01  2.1929E+01  5.0851E+01  7.6751E+00  1.6624E+01
 EBVSHRINKVR(%)  5.3961E-01  3.9049E+01  7.5844E+01  1.4761E+01  3.0485E+01
 RELATIVEINF(%)  9.9458E+01  3.2830E+01  1.5431E+01  6.2508E+01  2.9999E+01
 EPSSHRINKSD(%)  2.3957E+01
 EPSSHRINKVR(%)  4.2175E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3873.9660830986600     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2219.8767233302492     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.28
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3873.966       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.43E-01  9.71E-01  9.46E-01  9.95E-01  9.24E-01  1.02E+00  1.05E+00  7.19E-01  1.00E+00  9.80E-01  1.03E+00
 


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
+        1.20E+03
 
 TH 2
+        3.84E+00  6.52E+02
 
 TH 3
+       -7.12E+08 -2.19E+01  1.76E+10
 
 TH 4
+       -5.50E+00  8.15E+09 -1.67E+10  1.59E+10
 
 TH 5
+        2.62E+00 -3.58E+02  9.01E+09  1.90E+02  9.38E+02
 
 TH 6
+        6.97E+01  2.44E-01  2.81E+00  7.01E-01 -2.11E+00  1.90E+02
 
 TH 7
+        4.10E-01  1.16E+01  5.15E+09 -9.80E+09 -3.86E+00  8.31E-01  3.02E+09
 
 TH 8
+       -5.05E+09  4.91E+09 -1.01E+10  9.57E+09 -7.58E+00 -6.93E-01 -2.95E+09  2.88E+09
 
 TH 9
+       -3.10E-01 -5.91E+00  8.03E+09 -7.63E+09 -6.63E+08 -1.43E+00  4.70E+09  8.32E+00  7.32E+09
 
 TH10
+       -2.80E+00 -6.67E+08 -1.03E+01 -8.07E+09 -4.39E+01  3.67E+00  1.17E+01 -4.86E+09  3.18E+00  1.05E+02
 
 TH11
+       -6.04E+09  5.86E+09 -6.02E+09  5.72E+09 -6.16E+09  4.70E+09 -3.52E+09  7.54E+05 -5.49E+09 -5.81E+09  4.11E+09
 
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
 #CPUT: Total CPU Time in Seconds,       47.476
Stop Time:
Sat Sep 18 05:50:35 CDT 2021
