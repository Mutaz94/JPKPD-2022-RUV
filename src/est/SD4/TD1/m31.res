Sun Oct 24 03:44:47 CDT 2021
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
$DATA ../../../../data/SD4/TD1/dat31.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1636.16363893270        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6223E+02 -1.9788E+01 -1.8072E+01  1.6058E+01  7.0393E+01  2.7543E+01  7.2583E+00 -8.4692E-01  9.2957E+00 -1.2443E+01
            -1.2085E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1639.89423489087        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.7746E-01  1.0168E+00  9.8933E-01  1.0344E+00  9.5978E-01  1.0301E+00  9.7402E-01  1.0078E+00  9.8180E-01  1.0275E+00
             1.0337E+00
 PARAMETER:  7.7199E-02  1.1666E-01  8.9268E-02  1.3380E-01  5.8951E-02  1.2968E-01  7.3673E-02  1.0775E-01  8.1635E-02  1.2711E-01
             1.3311E-01
 GRADIENT:  -1.1629E+00  3.1124E+00 -4.8724E+00  9.6621E+00  1.0335E+01  1.2603E+00  4.5582E+00  4.3825E-01  2.0649E+00 -1.3498E+00
             2.4895E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1640.17432147255        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.7883E-01  9.7841E-01  1.0455E+00  1.0612E+00  9.6550E-01  1.0185E+00  8.1032E-01  9.9098E-01  1.0182E+00  1.0656E+00
             1.0333E+00
 PARAMETER:  7.8600E-02  7.8172E-02  1.4450E-01  1.5936E-01  6.4891E-02  1.1836E-01 -1.1033E-01  9.0942E-02  1.1808E-01  1.6351E-01
             1.3275E-01
 GRADIENT:   2.8607E+00  5.9746E+00  8.2586E-01  1.3189E+01  3.4851E+00 -3.1277E+00  1.0215E+00 -2.5890E+00  4.2389E+00 -9.1815E-01
             1.4134E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1640.57212206522        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.7737E-01  1.1108E+00  1.1181E+00  9.7167E-01  1.0482E+00  1.0276E+00  6.3408E-01  1.3230E+00  1.1058E+00  1.1104E+00
             1.0294E+00
 PARAMETER:  7.7108E-02  2.0510E-01  2.1165E-01  7.1258E-02  1.4703E-01  1.2721E-01 -3.5559E-01  3.7987E-01  2.0057E-01  2.0470E-01
             1.2894E-01
 GRADIENT:  -4.4081E-01  2.8070E+00  1.0859E-01  2.8824E+00 -7.3107E-01  4.5923E-01  2.1482E-01  5.5276E-02 -3.6491E-01  1.9432E-01
            -5.5166E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1640.61417464903        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.7696E-01  1.0237E+00  1.2155E+00  1.0272E+00  1.0476E+00  1.0259E+00  5.5430E-01  1.3769E+00  1.0813E+00  1.1209E+00
             1.0305E+00
 PARAMETER:  7.6686E-02  1.2342E-01  2.9519E-01  1.2685E-01  1.4649E-01  1.2557E-01 -4.9006E-01  4.1980E-01  1.7816E-01  2.1411E-01
             1.3001E-01
 GRADIENT:   3.0160E-01  1.8786E+00  9.2816E-01  1.5022E+00 -1.3853E+00  5.6970E-02  1.2767E-01 -4.3873E-02  1.7813E-01 -5.8672E-02
             1.5562E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1640.63560643381        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  9.7580E-01  8.9503E-01  1.2942E+00  1.1102E+00  1.0269E+00  1.0249E+00  4.3766E-01  1.3664E+00  1.0326E+00  1.1251E+00
             1.0285E+00
 PARAMETER:  7.5499E-02 -1.0902E-02  3.5788E-01  2.0456E-01  1.2659E-01  1.2458E-01 -7.2631E-01  4.1218E-01  1.3209E-01  2.1790E-01
             1.2807E-01
 GRADIENT:  -1.6183E-02  1.8498E+00  5.6466E-01  2.6094E+00 -1.0157E+00 -6.5324E-02  6.0911E-02 -1.3286E-01  2.9884E-01  8.5112E-02
             2.0903E-03

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1640.64358730508        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  9.7594E-01  8.5388E-01  1.3148E+00  1.1340E+00  1.0207E+00  1.0250E+00  3.5374E-01  1.3659E+00  1.0175E+00  1.1253E+00
             1.0279E+00
 PARAMETER:  7.5643E-02 -5.7959E-02  3.7372E-01  2.2577E-01  1.2049E-01  1.2473E-01 -9.3919E-01  4.1184E-01  1.1731E-01  2.1806E-01
             1.2752E-01
 GRADIENT:   1.0570E+00  2.2664E-02 -3.2344E-01 -8.6749E-01  5.4008E-01  1.0383E-01  2.2102E-02 -1.4744E-02 -7.4280E-01 -1.9482E-01
            -1.2231E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1640.64750556232        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1239
 NPARAMETR:  9.7505E-01  8.4800E-01  1.3214E+00  1.1369E+00  1.0205E+00  1.0246E+00  2.3814E-01  1.3729E+00  1.0297E+00  1.1312E+00
             1.0277E+00
 PARAMETER:  7.4735E-02 -6.4870E-02  3.7865E-01  2.2829E-01  1.2032E-01  1.2427E-01 -1.3349E+00  4.1695E-01  1.2928E-01  2.2329E-01
             1.2729E-01
 GRADIENT:  -6.9154E-01  1.9946E-01 -2.0994E-01 -6.3183E-02 -6.1676E-01 -7.3091E-02  3.7874E-02 -4.6722E-02  8.3066E-01  9.6580E-02
             5.0208E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1640.66041047098        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1419
 NPARAMETR:  9.7484E-01  8.2584E-01  1.3533E+00  1.1499E+00  1.0245E+00  1.0244E+00  5.5002E-02  1.4008E+00  1.0279E+00  1.1397E+00
             1.0271E+00
 PARAMETER:  7.4521E-02 -9.1355E-02  4.0252E-01  2.3971E-01  1.2421E-01  1.2414E-01 -2.8004E+00  4.3704E-01  1.2748E-01  2.3076E-01
             1.2669E-01
 GRADIENT:  -5.2926E-01 -1.6076E+00 -3.7587E-01 -2.2029E+00  1.0130E-01 -4.9835E-02  9.2785E-03  1.6330E-01  1.9294E+00  6.0146E-01
             1.1543E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1640.67017956867        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1596
 NPARAMETR:  9.7568E-01  8.1403E-01  1.3713E+00  1.1595E+00  1.0260E+00  1.0249E+00  1.0000E-02  1.4137E+00  1.0157E+00  1.1369E+00
             1.0275E+00
 PARAMETER:  7.5378E-02 -1.0575E-01  4.1579E-01  2.4800E-01  1.2570E-01  1.2462E-01 -4.7783E+00  4.4621E-01  1.1560E-01  2.2828E-01
             1.2717E-01
 GRADIENT:   1.4939E+00  2.6221E-01  1.1091E-02 -5.2149E-03  3.4219E-01  1.8583E-01  0.0000E+00  4.9021E-02  1.7293E-01  4.1647E-02
             2.0294E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1640.67017956867        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1618
 NPARAMETR:  9.7568E-01  8.1403E-01  1.3713E+00  1.1595E+00  1.0260E+00  1.0249E+00  1.0000E-02  1.4137E+00  1.0157E+00  1.1369E+00
             1.0275E+00
 PARAMETER:  7.5378E-02 -1.0575E-01  4.1579E-01  2.4800E-01  1.2570E-01  1.2462E-01 -4.7783E+00  4.4621E-01  1.1560E-01  2.2828E-01
             1.2717E-01
 GRADIENT:   1.4939E+00  2.6221E-01  1.1091E-02 -5.2149E-03  3.4219E-01  1.8583E-01  0.0000E+00  4.9021E-02  1.7293E-01  4.1647E-02
             2.0294E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1618
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.3740E-04 -5.3841E-04 -3.7903E-02 -4.1206E-03 -3.7265E-02
 SE:             2.9833E-02  1.7823E-04  1.6426E-02  2.9180E-02  2.2288E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7761E-01  2.5212E-03  2.1029E-02  8.8770E-01  9.4529E-02

 ETASHRINKSD(%)  5.7051E-02  9.9403E+01  4.4970E+01  2.2435E+00  2.5332E+01
 ETASHRINKVR(%)  1.1407E-01  9.9996E+01  6.9717E+01  4.4366E+00  4.4246E+01
 EBVSHRINKSD(%)  4.4327E-01  9.9487E+01  4.8684E+01  2.7349E+00  2.2361E+01
 EBVSHRINKVR(%)  8.8458E-01  9.9997E+01  7.3667E+01  5.3951E+00  3.9722E+01
 RELATIVEINF(%)  9.8765E+01  2.2153E-04  7.7719E+00  9.4185E+00  1.6529E+01
 EPSSHRINKSD(%)  4.5228E+01
 EPSSHRINKVR(%)  7.0000E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1640.6701795686661     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -905.51935300492789     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1640.670       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  8.14E-01  1.37E+00  1.16E+00  1.03E+00  1.02E+00  1.00E-02  1.41E+00  1.02E+00  1.14E+00  1.03E+00
 


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
 #CPUT: Total CPU Time in Seconds,       48.763
Stop Time:
Sun Oct 24 03:44:58 CDT 2021
