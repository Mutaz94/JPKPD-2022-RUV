Sun Oct 24 03:02:45 CDT 2021
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
$DATA ../../../../data/SD4/SL1/dat57.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1763.86552959112        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5832E+02 -6.6328E+01 -5.9214E+01 -7.4285E+00  6.5869E+01  6.6875E+01  1.0816E+01  1.4326E+01  3.2111E+01  2.2387E+01
             5.7449E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1777.96621703272        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0204E+00  1.1164E+00  1.1793E+00  1.0000E+00  1.0522E+00  8.7917E-01  9.2346E-01  9.2187E-01  8.5519E-01  8.8415E-01
             8.0365E-01
 PARAMETER:  1.2020E-01  2.1015E-01  2.6491E-01  1.0003E-01  1.5088E-01 -2.8772E-02  2.0369E-02  1.8651E-02 -5.6427E-02 -2.3127E-02
            -1.1860E-01
 GRADIENT:  -3.3287E+01  1.4481E+01  1.8002E+01 -6.0384E+00 -6.6563E+00 -2.1073E+01 -1.3611E+00 -3.6453E+00 -1.3082E+01 -1.8265E+01
            -4.5902E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1779.70955483823        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0266E+00  1.0380E+00  1.3464E+00  1.0505E+00  1.0713E+00  8.6839E-01  6.3949E-01  8.8597E-01  1.0074E+00  9.8221E-01
             8.2007E-01
 PARAMETER:  1.2627E-01  1.3733E-01  3.9746E-01  1.4928E-01  1.6887E-01 -4.1118E-02 -3.4708E-01 -2.1072E-02  1.0734E-01  8.2054E-02
            -9.8360E-02
 GRADIENT:  -1.3047E+01  1.4514E+01  2.5567E+01  2.0884E+00 -2.6837E+01 -2.5741E+01  2.2497E+00 -9.9979E+00  8.0921E+00 -1.3004E+01
            -3.6566E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1785.08371254298        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0331E+00  1.3122E+00  1.4661E+00  8.6988E-01  1.2458E+00  9.2776E-01  3.4316E-01  1.4472E+00  1.2523E+00  1.1379E+00
             9.0103E-01
 PARAMETER:  1.3260E-01  3.7171E-01  4.8260E-01 -3.9402E-02  3.1975E-01  2.5023E-02 -9.6955E-01  4.6966E-01  3.2500E-01  2.2923E-01
            -4.2212E-03
 GRADIENT:   8.0240E-01 -5.3385E+00 -5.0000E+00  7.8342E+00  6.1992E+00  1.3794E+00  2.1135E+00  1.5631E+00  6.5517E+00  1.2470E+00
             7.3562E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1785.68719585118        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      725             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0356E+00  1.4293E+00  1.4861E+00  7.7568E-01  1.2822E+00  9.2603E-01  1.4846E-01  1.5675E+00  1.4056E+00  1.1625E+00
             8.8802E-01
 PARAMETER:  1.3500E-01  4.5717E-01  4.9616E-01 -1.5401E-01  3.4858E-01  2.3156E-02 -1.8074E+00  5.4951E-01  4.4048E-01  2.5057E-01
            -1.8762E-02
 GRADIENT:   7.6958E+02  4.7475E+02  4.3799E+00  6.6374E+01  1.9311E+01  4.7685E+01  3.8111E+00  1.2999E-01  3.2873E+01  2.3845E+00
             3.8383E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1785.71303206987        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  1.0358E+00  1.4296E+00  1.4855E+00  7.7625E-01  1.2863E+00  9.2674E-01  1.1659E-01  1.5943E+00  1.4266E+00  1.1647E+00
             8.8842E-01
 PARAMETER:  1.3520E-01  4.5741E-01  4.9573E-01 -1.5328E-01  3.5178E-01  2.3916E-02 -2.0491E+00  5.6646E-01  4.5532E-01  2.5248E-01
            -1.8307E-02
 GRADIENT:   6.6342E+00  2.8296E+00 -4.0101E-01  1.2191E+00  3.3849E-01  4.5244E-01  6.9803E-02  9.7300E-03  4.7871E-01 -4.9312E-02
             1.0371E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1785.73731923408        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1083
 NPARAMETR:  1.0322E+00  1.4236E+00  1.4939E+00  7.7763E-01  1.2859E+00  9.2580E-01  6.6000E-02  1.5979E+00  1.4283E+00  1.1655E+00
             8.8821E-01
 PARAMETER:  1.3166E-01  4.5316E-01  5.0137E-01 -1.5150E-01  3.5143E-01  2.2905E-02 -2.6181E+00  5.6866E-01  4.5651E-01  2.5319E-01
            -1.8551E-02
 GRADIENT:  -2.5820E+00  1.0029E+00 -4.4564E-01 -1.0813E+00  5.6244E-01  3.3315E-02  3.0519E-02  4.9218E-02  8.3376E-02 -6.0472E-03
             1.3826E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1785.96162573704        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1261
 NPARAMETR:  1.0331E+00  1.3592E+00  1.5840E+00  8.2630E-01  1.2796E+00  9.2776E-01  1.0000E-02  1.6000E+00  1.3461E+00  1.1646E+00
             8.8718E-01
 PARAMETER:  1.3256E-01  4.0687E-01  5.5996E-01 -9.0800E-02  3.4655E-01  2.5016E-02 -8.9222E+00  5.6999E-01  3.9718E-01  2.5238E-01
            -1.9703E-02
 GRADIENT:   6.5900E-01  5.3780E+00  9.5344E-01  2.2081E+00 -8.1536E-01  1.0621E+00  0.0000E+00 -2.0801E-01 -1.1514E+00  3.5583E-02
            -2.8078E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1786.00959568827        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1438
 NPARAMETR:  1.0328E+00  1.3403E+00  1.5767E+00  8.3808E-01  1.2727E+00  9.2432E-01  1.0000E-02  1.5763E+00  1.3287E+00  1.1604E+00
             8.8709E-01
 PARAMETER:  1.3228E-01  3.9292E-01  5.5534E-01 -7.6640E-02  3.4117E-01  2.1305E-02 -1.0565E+01  5.5511E-01  3.8422E-01  2.4881E-01
            -1.9805E-02
 GRADIENT:   1.4297E-01  2.8489E-01  8.4750E-02  2.9225E-01  3.3665E-01 -3.3624E-01  0.0000E+00  7.2179E-02 -3.7165E-02  4.0977E-02
            -2.0598E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1786.01352227176        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     1586
 NPARAMETR:  1.0335E+00  1.3395E+00  1.5761E+00  8.3860E-01  1.2717E+00  9.2510E-01  1.0000E-02  1.5725E+00  1.3278E+00  1.1600E+00
             8.8696E-01
 PARAMETER:  1.3291E-01  3.9227E-01  5.5492E-01 -7.6016E-02  3.4037E-01  2.2145E-02 -1.0565E+01  5.5268E-01  3.8355E-01  2.4840E-01
            -1.9956E-02
 GRADIENT:   1.8339E+00  3.5399E-01  2.6945E-01  1.9040E-01 -1.4073E-01  2.4548E-04  0.0000E+00 -9.8849E-03 -1.2751E-02  4.1984E-02
            -8.0591E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1586
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6336E-03 -1.2374E-03 -3.0890E-02 -2.4568E-03 -3.3047E-02
 SE:             2.9899E-02  3.3688E-04  1.3040E-02  2.9238E-02  2.3423E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5643E-01  2.3957E-04  1.7845E-02  9.3304E-01  1.5827E-01

 ETASHRINKSD(%)  1.0000E-10  9.8871E+01  5.6314E+01  2.0479E+00  2.1531E+01
 ETASHRINKVR(%)  1.0000E-10  9.9987E+01  8.0915E+01  4.0539E+00  3.8426E+01
 EBVSHRINKSD(%)  3.8457E-01  9.9049E+01  6.0762E+01  2.5117E+00  1.8135E+01
 EBVSHRINKVR(%)  7.6766E-01  9.9991E+01  8.4604E+01  4.9602E+00  3.2982E+01
 RELATIVEINF(%)  9.9132E+01  1.2121E-03  6.6492E+00  1.4050E+01  2.9782E+01
 EPSSHRINKSD(%)  4.4007E+01
 EPSSHRINKVR(%)  6.8648E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1786.0135222717615     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1050.8626957080232     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1786.014       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.34E+00  1.58E+00  8.39E-01  1.27E+00  9.25E-01  1.00E-02  1.57E+00  1.33E+00  1.16E+00  8.87E-01
 


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
 #CPUT: Total CPU Time in Seconds,       44.971
Stop Time:
Sun Oct 24 03:02:55 CDT 2021
