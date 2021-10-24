Sat Oct 23 22:25:08 CDT 2021
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
$DATA ../../../../data/SD3/A2/dat96.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1344.66817688840        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1029E+02  7.8841E+01  8.4366E-02  1.1848E+02  1.2572E+02  3.0289E+01 -3.2656E+00  1.8266E+01  1.3580E+01 -7.6048E+01
            -1.3870E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1724.75978963816        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  1.1438E+00  9.0349E-01  1.1417E+00  1.0531E+00  9.2510E-01  1.2096E+00  8.3681E-01  7.3331E-01  8.4616E-01  1.1382E+00
             2.3041E+00
 PARAMETER:  2.3438E-01 -1.4892E-03  2.3256E-01  1.5178E-01  2.2148E-02  2.9025E-01 -7.8161E-02 -2.1019E-01 -6.7048E-02  2.2948E-01
             9.3469E-01
 GRADIENT:   4.3210E+02  1.3438E+01  1.2125E+01 -3.4094E+00 -3.5302E+01  4.1380E+01 -5.3121E-01  6.2545E+00 -1.3499E+01  1.3676E+01
             8.4712E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1751.24760443678        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0492E+00  8.0724E-01  5.1119E-01  1.0714E+00  6.1256E-01  1.0550E+00  9.0557E-01  7.5138E-02  1.0250E+00  7.6846E-01
             2.1109E+00
 PARAMETER:  1.4803E-01 -1.1413E-01 -5.7101E-01  1.6898E-01 -3.9010E-01  1.5350E-01  8.0672E-04 -2.4884E+00  1.2468E-01 -1.6336E-01
             8.4711E-01
 GRADIENT:   2.7940E+02 -1.2274E+01 -5.5880E+01  8.0617E+01  1.2167E+02  1.8433E+01  6.9691E+00  1.7088E-01  2.7831E+01  2.6020E+00
             3.9235E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1752.85640689715        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  1.0226E+00  7.2441E-01  3.6347E-01  1.1080E+00  4.8497E-01  1.0136E+00  9.0125E-01  5.2340E-02  9.9178E-01  6.1286E-01
             2.0486E+00
 PARAMETER:  1.2231E-01 -2.2240E-01 -9.1206E-01  2.0253E-01 -6.2366E-01  1.1349E-01 -3.9771E-03 -2.8500E+00  9.1747E-02 -3.8961E-01
             8.1714E-01
 GRADIENT:   2.1219E+02 -1.7850E+01 -1.0170E+02  1.6778E+02  2.1746E+02  5.1193E+00  4.2153E+00  5.8550E-02  1.0706E+01 -7.0180E+00
             2.0867E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1767.33969836877        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      396
 NPARAMETR:  9.6551E-01  5.9524E-01  5.2856E-01  1.1781E+00  5.1821E-01  9.7986E-01  1.0410E+00  7.2063E-02  8.5868E-01  7.0575E-01
             1.9855E+00
 PARAMETER:  6.4896E-02 -4.1879E-01 -5.3760E-01  2.6390E-01 -5.5738E-01  7.9654E-02  1.4017E-01 -2.5302E+00 -5.2359E-02 -2.4850E-01
             7.8586E-01
 GRADIENT:  -1.2105E+01  3.0856E+01  4.7505E+00  5.2929E+01 -9.1289E+00 -1.0354E+01 -7.0081E-01  1.6859E-01 -5.0309E+00 -3.1142E+00
            -1.7380E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1772.13403042815        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      571
 NPARAMETR:  9.6627E-01  3.5854E-01  5.3917E-01  1.2708E+00  4.6062E-01  1.0029E+00  1.3320E+00  3.7552E-02  8.2805E-01  7.1346E-01
             2.0346E+00
 PARAMETER:  6.5688E-02 -9.2571E-01 -5.1772E-01  3.3963E-01 -6.7518E-01  1.0286E-01  3.8668E-01 -3.1820E+00 -8.8684E-02 -2.3763E-01
             8.1031E-01
 GRADIENT:   8.4050E-01  1.1508E+01  2.3344E+01  6.9567E+00 -4.0148E+01  1.0770E+00 -4.8012E-01  5.1247E-02  9.5826E-01  1.2237E+00
             8.0411E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1772.91461793632        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      747
 NPARAMETR:  9.6676E-01  2.6634E-01  5.1347E-01  1.2901E+00  4.3199E-01  1.0109E+00  1.8410E+00  1.7108E-02  8.2474E-01  7.0728E-01
             2.0138E+00
 PARAMETER:  6.6192E-02 -1.2230E+00 -5.6657E-01  3.5473E-01 -7.3935E-01  1.1082E-01  7.1030E-01 -3.9682E+00 -9.2688E-02 -2.4633E-01
             8.0003E-01
 GRADIENT:   9.1197E+00  3.6475E+00  9.2932E+00 -2.1577E+01 -2.1421E+01  5.0258E+00  2.3003E+00  1.0290E-02  4.9520E+00  2.5760E+00
             1.6133E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1775.77180946184        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      925
 NPARAMETR:  9.5246E-01  8.7288E-02  7.0894E-01  1.4572E+00  5.1209E-01  9.9078E-01  3.7280E+00  1.0000E-02  7.3861E-01  7.9179E-01
             2.0277E+00
 PARAMETER:  5.1297E-02 -2.3385E+00 -2.4399E-01  4.7650E-01 -5.6925E-01  9.0739E-02  1.4159E+00 -5.2066E+00 -2.0299E-01 -1.3346E-01
             8.0689E-01
 GRADIENT:  -7.0444E+00  2.4319E+00  5.4485E+00  1.9199E+01 -9.3235E+00 -1.0810E+00  9.8170E-01  0.0000E+00 -4.3606E+00 -1.2466E+00
             4.1897E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1776.51967093094        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1101
 NPARAMETR:  9.5200E-01  2.8224E-02  6.9699E-01  1.4732E+00  4.9889E-01  9.9253E-01  2.9519E+00  1.0000E-02  7.4453E-01  7.9653E-01
             2.0186E+00
 PARAMETER:  5.0806E-02 -3.4676E+00 -2.6098E-01  4.8744E-01 -5.9538E-01  9.2503E-02  1.1825E+00 -6.6842E+00 -1.9500E-01 -1.2749E-01
             8.0239E-01
 GRADIENT:  -2.3620E+00  1.9571E-01 -1.0549E+00 -2.6266E+00  1.2411E+00  1.9764E-01 -3.7350E-04  0.0000E+00 -2.1584E-01  8.4863E-02
             1.8977E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1776.63380834352        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1276
 NPARAMETR:  9.5260E-01  1.0000E-02  7.1219E-01  1.4873E+00  5.0325E-01  9.9163E-01  2.3568E+00  1.0000E-02  7.3989E-01  8.0430E-01
             2.0139E+00
 PARAMETER:  5.1440E-02 -4.5706E+00 -2.3942E-01  4.9696E-01 -5.8668E-01  9.1596E-02  9.5730E-01 -8.0639E+00 -2.0126E-01 -1.1778E-01
             8.0008E-01
 GRADIENT:   4.9874E-01  0.0000E+00  7.0931E-02 -5.6268E-01  2.5886E-02 -1.8391E-02 -7.4838E-05  0.0000E+00  2.8170E-03 -3.5671E-02
            -1.4033E-01

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1776.63380834352        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1298
 NPARAMETR:  9.5260E-01  1.0000E-02  7.1219E-01  1.4873E+00  5.0325E-01  9.9163E-01  2.3568E+00  1.0000E-02  7.3989E-01  8.0430E-01
             2.0139E+00
 PARAMETER:  5.1440E-02 -4.5706E+00 -2.3942E-01  4.9696E-01 -5.8668E-01  9.1596E-02  9.5730E-01 -8.0639E+00 -2.0126E-01 -1.1778E-01
             8.0008E-01
 GRADIENT:   4.9874E-01  0.0000E+00  7.0931E-02 -5.6268E-01  2.5886E-02 -1.8391E-02 -7.4838E-05  0.0000E+00  2.8170E-03 -3.5671E-02
            -1.4033E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1298
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.2810E-04 -6.0941E-04 -1.2601E-05 -7.6819E-03 -1.3046E-02
 SE:             2.9522E-02  4.1597E-04  1.8319E-04  2.8184E-02  2.2894E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8573E-01  1.4292E-01  9.4516E-01  7.8519E-01  5.6878E-01

 ETASHRINKSD(%)  1.0976E+00  9.8606E+01  9.9386E+01  5.5800E+00  2.3303E+01
 ETASHRINKVR(%)  2.1831E+00  9.9981E+01  9.9996E+01  1.0849E+01  4.1176E+01
 EBVSHRINKSD(%)  1.2445E+00  9.8670E+01  9.9297E+01  5.3360E+00  2.2797E+01
 EBVSHRINKVR(%)  2.4734E+00  9.9982E+01  9.9995E+01  1.0387E+01  4.0397E+01
 RELATIVEINF(%)  9.0098E+01  9.7856E-04  3.4863E-04  6.8703E+00  3.6910E+00
 EPSSHRINKSD(%)  2.7680E+01
 EPSSHRINKVR(%)  4.7698E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1776.6338083435153     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -857.69527513884259     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1776.634       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.53E-01  1.00E-02  7.12E-01  1.49E+00  5.03E-01  9.92E-01  2.36E+00  1.00E-02  7.40E-01  8.04E-01  2.01E+00
 


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
 #CPUT: Total CPU Time in Seconds,       99.723
Stop Time:
Sat Oct 23 22:25:25 CDT 2021
