Sun Oct 24 02:38:54 CDT 2021
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
$DATA ../../../../data/SD4/A3/dat100.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -204.567122797485        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3444E+02  8.3176E+01  8.5340E+01  2.4041E+01  2.1701E+02  3.9158E+01 -5.9904E+01 -6.0984E+01 -1.4908E+02 -1.2252E+02
            -2.4083E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1217.44947353242        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0209E+00  9.2721E-01  8.6673E-01  1.1055E+00  8.5675E-01  8.4718E-01  1.0346E+00  1.0585E+00  1.2058E+00  9.7581E-01
             4.3919E+00
 PARAMETER:  1.2073E-01  2.4428E-02 -4.3028E-02  2.0029E-01 -5.4613E-02 -6.5840E-02  1.3402E-01  1.5688E-01  2.8710E-01  7.5513E-02
             1.5798E+00
 GRADIENT:   1.3545E+02 -2.2618E+01 -1.3456E+01 -1.5146E+01  2.6653E+00 -1.0542E+01  1.0015E+01  7.0396E+00  2.3329E+01  2.4601E+01
             1.3073E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1232.75901608725        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      181
 NPARAMETR:  1.0166E+00  6.3156E-01  3.0320E-01  1.2315E+00  3.5812E-01  8.8632E-01  1.0122E+00  8.2941E-01  1.3243E+00  4.0575E-01
             3.7265E+00
 PARAMETER:  1.1647E-01 -3.5957E-01 -1.0934E+00  3.0823E-01 -9.2688E-01 -2.0676E-02  1.1214E-01 -8.7041E-02  3.8085E-01 -8.0201E-01
             1.4155E+00
 GRADIENT:   9.3032E+01  6.5052E+01  1.9858E+01  9.3031E+01 -6.8887E+01 -1.4503E+01  9.7783E-01  5.8609E+00  1.2007E+01  6.6359E+00
             8.0330E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1248.06250297882        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      358
 NPARAMETR:  9.7309E-01  5.7704E-01  2.4498E-01  1.1554E+00  3.1974E-01  9.1456E-01  1.0483E+00  3.4408E-01  1.3470E+00  3.0430E-01
             3.4614E+00
 PARAMETER:  7.2721E-02 -4.4984E-01 -1.3066E+00  2.4448E-01 -1.0403E+00  1.0689E-02  1.4720E-01 -9.6689E-01  3.9787E-01 -1.0897E+00
             1.3417E+00
 GRADIENT:   8.3307E+00  4.0650E+01  1.3833E+01  6.1191E+01 -4.9130E+01 -4.8446E+00 -1.3153E+00 -1.0875E+00  1.3078E+01 -1.6103E-01
             5.2548E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1255.10173858407        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  9.6165E-01  5.6654E-01  2.2004E-01  1.0610E+00  3.0936E-01  9.2504E-01  1.1826E+00  2.9619E-01  1.2766E+00  2.6381E-01
             3.0528E+00
 PARAMETER:  6.0891E-02 -4.6821E-01 -1.4139E+00  1.5920E-01 -1.0732E+00  2.2083E-02  2.6771E-01 -1.1168E+00  3.4422E-01 -1.2325E+00
             1.2161E+00
 GRADIENT:   7.2761E+00  7.9585E+00  3.9156E+00  7.6847E+00 -4.8742E+00 -2.7400E+00 -8.5099E-01 -1.5272E+00 -1.4272E+00 -9.4596E-01
            -3.6859E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1259.77536985413        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      712
 NPARAMETR:  9.4943E-01  5.0613E-01  1.7945E-01  1.0220E+00  2.6607E-01  9.3159E-01  1.0617E+00  1.2024E+00  1.3837E+00  2.4883E-01
             2.8243E+00
 PARAMETER:  4.8107E-02 -5.8095E-01 -1.6179E+00  1.2174E-01 -1.2240E+00  2.9134E-02  1.5985E-01  2.8428E-01  4.2474E-01 -1.2910E+00
             1.1383E+00
 GRADIENT:  -6.0192E+00  4.0298E-02 -2.8685E+00 -7.6583E+00 -1.3874E+01  1.2806E+00  7.4942E-01  8.4149E-01 -8.8621E+00 -7.5502E-01
             1.1342E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1260.79472214229        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      888
 NPARAMETR:  9.4484E-01  6.0410E-01  1.8928E-01  1.0120E+00  3.0023E-01  9.2184E-01  1.0104E+00  1.3019E+00  1.4371E+00  1.7876E-01
             2.7647E+00
 PARAMETER:  4.3265E-02 -4.0401E-01 -1.5645E+00  1.1193E-01 -1.1032E+00  1.8621E-02  1.1039E-01  3.6379E-01  4.6263E-01 -1.6217E+00
             1.1169E+00
 GRADIENT:  -3.9649E+00 -7.7073E-01 -3.8305E+00 -2.8380E-01  7.6762E+00  5.6676E-01  1.6201E+00  3.9055E-01  1.4259E-01  9.1592E-01
             2.8257E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1261.16790794422        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  9.4696E-01  6.0524E-01  1.9396E-01  1.0159E+00  3.0304E-01  9.2107E-01  1.0143E+00  1.3166E+00  1.4309E+00  4.6687E-02
             2.7588E+00
 PARAMETER:  4.5501E-02 -4.0213E-01 -1.5401E+00  1.1573E-01 -1.0939E+00  1.7777E-02  1.1424E-01  3.7507E-01  4.5829E-01 -2.9643E+00
             1.1148E+00
 GRADIENT:   1.3949E-01 -2.9707E-01 -2.5307E-01 -3.4750E-01  1.0660E+00  4.3998E-01 -3.5171E-01 -1.8556E-02 -2.9931E-01  4.4456E-02
            -3.2081E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1261.19073960495        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1239
 NPARAMETR:  9.4713E-01  6.0332E-01  1.9429E-01  1.0173E+00  3.0263E-01  9.1993E-01  1.0222E+00  1.3150E+00  1.4316E+00  1.0000E-02
             2.7602E+00
 PARAMETER:  4.5686E-02 -4.0530E-01 -1.5384E+00  1.1716E-01 -1.0953E+00  1.6542E-02  1.2191E-01  3.7380E-01  4.5882E-01 -4.5774E+00
             1.1153E+00
 GRADIENT:   1.3138E-01 -8.6897E-03 -2.1746E-03  4.4172E-02 -5.0510E-03 -1.6024E-02 -3.7080E-03 -3.9223E-02  4.3789E-02  0.0000E+00
             2.4403E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1261.19073960495        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1261
 NPARAMETR:  9.4713E-01  6.0332E-01  1.9429E-01  1.0173E+00  3.0263E-01  9.1993E-01  1.0222E+00  1.3150E+00  1.4316E+00  1.0000E-02
             2.7602E+00
 PARAMETER:  4.5686E-02 -4.0530E-01 -1.5384E+00  1.1716E-01 -1.0953E+00  1.6542E-02  1.2191E-01  3.7380E-01  4.5882E-01 -4.5774E+00
             1.1153E+00
 GRADIENT:   1.3138E-01 -8.6897E-03 -2.1746E-03  4.4172E-02 -5.0510E-03 -1.6024E-02 -3.7080E-03 -3.9223E-02  4.3789E-02  0.0000E+00
             2.4403E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1261
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2818E-03 -4.0052E-03 -2.8582E-03 -7.1862E-03  1.9353E-04
 SE:             2.8789E-02  1.9328E-02  1.7852E-02  2.6454E-02  3.2059E-04
 N:                     100         100         100         100         100

 P VAL.:         9.6449E-01  8.3584E-01  8.7280E-01  7.8589E-01  5.4605E-01

 ETASHRINKSD(%)  3.5538E+00  3.5248E+01  4.0192E+01  1.1376E+01  9.8926E+01
 ETASHRINKVR(%)  6.9814E+00  5.8071E+01  6.4230E+01  2.1457E+01  9.9988E+01
 EBVSHRINKSD(%)  3.4797E+00  3.4978E+01  4.0582E+01  1.0520E+01  9.8943E+01
 EBVSHRINKVR(%)  6.8383E+00  5.7721E+01  6.4695E+01  1.9933E+01  9.9989E+01
 RELATIVEINF(%)  8.8354E+01  2.0512E+00  8.5742E+00  3.9001E+01  4.5868E-04
 EPSSHRINKSD(%)  3.6256E+01
 EPSSHRINKVR(%)  5.9367E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1261.1907396049535     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -526.03991304121530     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.94
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1261.191       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.47E-01  6.03E-01  1.94E-01  1.02E+00  3.03E-01  9.20E-01  1.02E+00  1.31E+00  1.43E+00  1.00E-02  2.76E+00
 


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
 #CPUT: Total CPU Time in Seconds,       36.950
Stop Time:
Sun Oct 24 02:39:02 CDT 2021
