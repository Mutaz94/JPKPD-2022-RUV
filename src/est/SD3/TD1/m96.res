Sun Oct 24 00:36:54 CDT 2021
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
$DATA ../../../../data/SD3/TD1/dat96.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2125.63221581077        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0044E+02  1.6235E+01 -6.6293E+01  1.1473E+02  5.8080E+01  3.3336E+01  1.4980E+01  1.9453E+01  2.7602E+01  3.6372E+01
             6.8723E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2135.04269555756        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  1.0026E+00  1.0770E+00  1.2661E+00  9.5450E-01  1.0639E+00  1.0707E+00  9.2076E-01  8.7856E-01  8.9754E-01  7.9237E-01
             9.8850E-01
 PARAMETER:  1.0263E-01  1.7417E-01  3.3596E-01  5.3435E-02  1.6195E-01  1.6828E-01  1.7440E-02 -2.9474E-02 -8.0981E-03 -1.3272E-01
             8.8436E-02
 GRADIENT:   6.6794E+01  2.8173E+01  2.4324E+01 -7.3070E+00 -2.2773E+01  1.6312E+01  2.8369E+00  3.3112E+00 -1.7137E+01 -1.8181E+01
            -1.8865E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2138.78607075749        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      347
 NPARAMETR:  1.0021E+00  8.6791E-01  1.5386E+00  1.0805E+00  1.0622E+00  1.0587E+00  6.1323E-01  4.9728E-01  9.4918E-01  9.8040E-01
             1.0231E+00
 PARAMETER:  1.0210E-01 -4.1672E-02  5.3088E-01  1.7745E-01  1.6033E-01  1.5702E-01 -3.8902E-01 -5.9860E-01  4.7838E-02  8.0204E-02
             1.2288E-01
 GRADIENT:   7.0833E+01  1.6824E+01  3.6938E+01 -1.0710E+01 -5.6372E+01  1.2909E+01  4.7987E-01 -5.9319E-01 -4.5841E+00 -2.0739E+00
             1.1682E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2143.11102868473        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      523
 NPARAMETR:  9.6797E-01  8.9090E-01  1.2041E+00  1.0656E+00  9.9966E-01  1.0187E+00  7.7428E-01  2.8953E-01  9.3819E-01  9.0740E-01
             1.0003E+00
 PARAMETER:  6.7450E-02 -1.5522E-02  2.8571E-01  1.6350E-01  9.9656E-02  1.1856E-01 -1.5582E-01 -1.1395E+00  3.6192E-02  2.8302E-03
             1.0027E-01
 GRADIENT:   1.5285E+00  4.1661E+00 -4.2376E+00  1.0211E+01  5.9398E+00  2.6802E-01 -7.7525E-01  5.1918E-01 -1.8295E+00  7.0771E-01
            -6.2320E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2143.30737081074        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      699
 NPARAMETR:  9.6632E-01  7.7702E-01  1.1820E+00  1.1297E+00  9.4641E-01  1.0164E+00  9.6554E-01  1.6446E-01  8.7840E-01  8.7178E-01
             1.0013E+00
 PARAMETER:  6.5744E-02 -1.5229E-01  2.6723E-01  2.2192E-01  4.4916E-02  1.1623E-01  6.4937E-02 -1.7051E+00 -2.9652E-02 -3.7216E-02
             1.0129E-01
 GRADIENT:   7.6082E-02 -2.4126E+00 -2.9343E+00 -9.6236E-01  4.1598E+00 -2.7759E-01  1.4644E-01  1.4609E-01  4.9513E-01  3.9159E-01
             3.6158E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2143.41924873024        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      876
 NPARAMETR:  9.6683E-01  8.3972E-01  1.1590E+00  1.0908E+00  9.5926E-01  1.0175E+00  9.1626E-01  5.1317E-02  9.0317E-01  8.7347E-01
             1.0009E+00
 PARAMETER:  6.6269E-02 -7.4683E-02  2.4758E-01  1.8691E-01  5.8404E-02  1.1733E-01  1.2541E-02 -2.8697E+00 -1.8492E-03 -3.5283E-02
             1.0091E-01
 GRADIENT:  -2.4136E-01  4.6909E-02 -1.7720E-01  1.1216E-01  9.3861E-02 -8.6174E-02  1.4043E-01  1.4706E-02  2.0634E-01  1.7663E-01
             1.5140E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2143.42819954532        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1054             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6768E-01  8.4263E-01  1.1608E+00  1.0888E+00  9.6106E-01  1.0181E+00  9.0668E-01  1.1045E-02  9.0517E-01  8.7512E-01
             1.0008E+00
 PARAMETER:  6.7150E-02 -7.1223E-02  2.4913E-01  1.8509E-01  6.0278E-02  1.1794E-01  2.0355E-03 -4.4058E+00  3.6241E-04 -3.3389E-02
             1.0084E-01
 GRADIENT:   4.3061E+02  2.9714E+01  3.8623E+00  1.7259E+02  7.1972E+00  5.5444E+01  2.1696E+00  1.4111E-03  8.4009E+00  6.1260E-01
             1.0080E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2143.42831348358        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1215             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6771E-01  8.4292E-01  1.1607E+00  1.0886E+00  9.6105E-01  1.0182E+00  9.0655E-01  1.0000E-02  9.0533E-01  8.7516E-01
             1.0008E+00
 PARAMETER:  6.7173E-02 -7.0885E-02  2.4900E-01  1.8494E-01  6.0269E-02  1.1804E-01  1.8868E-03 -5.4060E+00  5.4548E-04 -3.3346E-02
             1.0084E-01
 GRADIENT:   4.3065E+02  2.9776E+01  3.9123E+00  1.7239E+02  7.0694E+00  5.5525E+01  2.1761E+00  0.0000E+00  8.4172E+00  6.3248E-01
             1.0167E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2143.42838958468        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1401
 NPARAMETR:  9.6771E-01  8.4307E-01  1.1605E+00  1.0885E+00  9.6111E-01  1.0182E+00  9.0585E-01  1.0000E-02  9.0532E-01  8.7508E-01
             1.0008E+00
 PARAMETER:  6.7174E-02 -7.0705E-02  2.4881E-01  1.8484E-01  6.0335E-02  1.1804E-01  1.1145E-03 -5.4060E+00  5.3440E-04 -3.3443E-02
             1.0084E-01
 GRADIENT:   1.6366E+00 -9.1518E-02  5.7707E-02 -3.4776E-01 -7.2417E-02  1.9040E-01  6.4277E-03  0.0000E+00  5.6688E-03  7.6755E-03
             2.8834E-05

0ITERATION NO.:   41    OBJECTIVE VALUE:  -2143.42838958468        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1423
 NPARAMETR:  9.6771E-01  8.4307E-01  1.1605E+00  1.0885E+00  9.6111E-01  1.0182E+00  9.0585E-01  1.0000E-02  9.0532E-01  8.7508E-01
             1.0008E+00
 PARAMETER:  6.7174E-02 -7.0705E-02  2.4881E-01  1.8484E-01  6.0335E-02  1.1804E-01  1.1145E-03 -5.4060E+00  5.3440E-04 -3.3443E-02
             1.0084E-01
 GRADIENT:   1.6366E+00 -9.1518E-02  5.7707E-02 -3.4776E-01 -7.2417E-02  1.9040E-01  6.4277E-03  0.0000E+00  5.6688E-03  7.6755E-03
             2.8834E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1423
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1364E-04 -1.2305E-02 -2.4394E-04  2.6752E-03 -1.9725E-02
 SE:             2.9870E-02  1.4771E-02  1.6065E-04  2.6767E-02  2.4612E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9429E-01  4.0481E-01  1.2890E-01  9.2039E-01  4.2287E-01

 ETASHRINKSD(%)  1.0000E-10  5.0516E+01  9.9462E+01  1.0326E+01  1.7546E+01
 ETASHRINKVR(%)  1.0000E-10  7.5514E+01  9.9997E+01  1.9586E+01  3.2013E+01
 EBVSHRINKSD(%)  3.1493E-01  5.0577E+01  9.9464E+01  1.0469E+01  1.6150E+01
 EBVSHRINKVR(%)  6.2887E-01  7.5573E+01  9.9997E+01  1.9843E+01  2.9691E+01
 RELATIVEINF(%)  9.8596E+01  8.3477E-01  3.5658E-04  3.4798E+00  8.5092E+00
 EPSSHRINKSD(%)  3.2084E+01
 EPSSHRINKVR(%)  5.3875E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2143.4283895846811     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1224.4898563800084     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2143.428       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  8.43E-01  1.16E+00  1.09E+00  9.61E-01  1.02E+00  9.06E-01  1.00E-02  9.05E-01  8.75E-01  1.00E+00
 


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
 #CPUT: Total CPU Time in Seconds,       46.818
Stop Time:
Sun Oct 24 00:37:04 CDT 2021
