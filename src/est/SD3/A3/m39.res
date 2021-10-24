Sat Oct 23 22:36:25 CDT 2021
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
$DATA ../../../../data/SD3/A3/dat39.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m39.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   671.594351493371        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3346E+02  1.2844E+02  2.1394E+02  5.9565E+01  2.1331E+02  5.7213E+01 -8.7627E+01 -1.4969E+02 -4.9666E+01 -1.7109E+02
            -5.0539E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1438.14650375504        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0380E+00  9.9740E-01  9.2450E-01  1.1251E+00  9.8984E-01  8.1939E-01  9.6669E-01  9.7794E-01  8.9064E-01  1.0036E+00
             5.3524E+00
 PARAMETER:  1.3728E-01  9.7393E-02  2.1500E-02  2.1783E-01  8.9784E-02 -9.9200E-02  6.6119E-02  7.7693E-02 -1.5819E-02  1.0356E-01
             1.7775E+00
 GRADIENT:  -2.5816E+01 -1.1982E+01 -1.7744E+01 -5.0120E-01 -5.7283E+00 -1.4612E+01  1.3099E+01  7.4821E+00  2.6332E+01  2.1988E+01
             3.1546E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1491.49523605513        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0214E+00  6.7208E-01  3.2178E-01  1.2067E+00  3.9501E-01  8.5509E-01  1.0405E+00  1.0658E-02  9.8728E-01  1.9899E-01
             4.3332E+00
 PARAMETER:  1.2116E-01 -2.9738E-01 -1.0339E+00  2.8787E-01 -8.2885E-01 -5.6548E-02  1.3969E-01 -4.4415E+00  8.7198E-02 -1.5145E+00
             1.5663E+00
 GRADIENT:  -4.9665E+01  8.1867E+01  2.8498E+01  1.0446E+02 -7.4161E+01 -2.0811E+01  8.9724E+00  1.8138E-03  1.2889E+01  1.7030E+00
             1.7211E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1512.63286693594        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0110E+00  5.9669E-01  2.5764E-01  1.1276E+00  3.4045E-01  8.7267E-01  3.9630E-01  1.0000E-02  1.2404E+00  2.6167E-01
             3.5631E+00
 PARAMETER:  1.1090E-01 -4.1636E-01 -1.2562E+00  2.2006E-01 -9.7750E-01 -3.6194E-02 -8.2557E-01 -6.5486E+00  3.1547E-01 -1.2407E+00
             1.3706E+00
 GRADIENT:  -1.4183E+01  9.7731E+01  5.2352E+01  6.1697E+01 -6.1452E+01 -1.6695E+01 -4.0120E+00  0.0000E+00  1.2586E+01 -6.7018E+00
             5.4579E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1516.05392747612        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      345
 NPARAMETR:  1.0126E+00  5.3978E-01  2.6432E-01  1.1329E+00  3.2901E-01  9.0800E-01  3.8684E-01  1.0000E-02  1.1996E+00  2.6582E-01
             3.4198E+00
 PARAMETER:  1.1247E-01 -5.1659E-01 -1.2306E+00  2.2480E-01 -1.0117E+00  3.4899E-03 -8.4974E-01 -6.9825E+00  2.8198E-01 -1.2249E+00
             1.3296E+00
 GRADIENT:  -3.0018E+01  8.5715E+01  6.7159E+01  3.5765E+01 -1.1319E+02 -4.2186E+00 -4.1021E+00  0.0000E+00  4.7474E+00 -9.8776E+00
            -4.9867E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1529.75037643835        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      522
 NPARAMETR:  1.0252E+00  3.5675E-01  1.9861E-01  1.0973E+00  2.4697E-01  9.2656E-01  1.1569E-01  1.0000E-02  1.2386E+00  5.0748E-01
             3.2578E+00
 PARAMETER:  1.2491E-01 -9.3071E-01 -1.5164E+00  1.9286E-01 -1.2985E+00  2.3719E-02 -2.0569E+00 -9.2707E+00  3.1400E-01 -5.7830E-01
             1.2810E+00
 GRADIENT:   5.4996E+00  1.7815E+01  2.1507E+01  2.4847E+01 -5.3318E+01  3.4954E-01 -1.0234E-02  0.0000E+00 -1.6480E+00  3.9837E-01
            -2.6690E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1530.75776401007        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      699
 NPARAMETR:  1.0228E+00  3.4808E-01  2.0178E-01  1.0821E+00  2.5072E-01  9.2264E-01  8.3393E-02  1.0000E-02  1.2321E+00  4.7509E-01
             3.3535E+00
 PARAMETER:  1.2259E-01 -9.5533E-01 -1.5006E+00  1.7892E-01 -1.2834E+00  1.9480E-02 -2.3842E+00 -9.2412E+00  3.0873E-01 -6.4424E-01
             1.3100E+00
 GRADIENT:   4.1648E-01  1.2889E+00  1.4764E+00 -2.6030E+00 -2.6206E+00 -1.0205E-01 -4.4673E-03  0.0000E+00  7.8394E-01  5.7712E-02
             4.8261E-01

0ITERATION NO.:   32    OBJECTIVE VALUE:  -1530.76837964812        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      756
 NPARAMETR:  1.0228E+00  3.4132E-01  2.0375E-01  1.0867E+00  2.5110E-01  9.2309E-01  7.9041E-02  1.0000E-02  1.2191E+00  4.7010E-01
             3.3569E+00
 PARAMETER:  1.2255E-01 -9.7493E-01 -1.4909E+00  1.8315E-01 -1.2819E+00  1.9968E-02 -2.4378E+00 -9.2167E+00  2.9815E-01 -6.5480E-01
             1.3110E+00
 GRADIENT:   6.8581E-01 -6.2690E-01 -4.6866E-01 -2.1017E+00  8.9167E-01  2.5465E-01 -5.0037E-03  0.0000E+00 -5.6251E-01 -1.2603E-01
             2.4932E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      756
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8309E-04 -6.8448E-04  1.9038E-04 -1.0672E-02  3.1001E-03
 SE:             2.8710E-02  9.0943E-04  2.2943E-04  2.6451E-02  1.8370E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8380E-01  4.5166E-01  4.0666E-01  6.8660E-01  8.6598E-01

 ETASHRINKSD(%)  3.8189E+00  9.6953E+01  9.9231E+01  1.1386E+01  3.8459E+01
 ETASHRINKVR(%)  7.4920E+00  9.9907E+01  9.9994E+01  2.1475E+01  6.2127E+01
 EBVSHRINKSD(%)  3.5892E+00  9.7052E+01  9.9246E+01  9.7669E+00  3.8790E+01
 EBVSHRINKVR(%)  7.0495E+00  9.9913E+01  9.9994E+01  1.8580E+01  6.2533E+01
 RELATIVEINF(%)  9.2514E+01  1.1815E-02  3.0905E-04  3.6480E+01  1.1259E+00
 EPSSHRINKSD(%)  2.2005E+01
 EPSSHRINKVR(%)  3.9168E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1530.7683796481226     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -611.82984644344992     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1530.768       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  3.41E-01  2.04E-01  1.09E+00  2.51E-01  9.23E-01  7.90E-02  1.00E-02  1.22E+00  4.70E-01  3.36E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       62.223
Stop Time:
Sat Oct 23 22:36:36 CDT 2021
