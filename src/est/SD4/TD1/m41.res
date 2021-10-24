Sun Oct 24 03:46:43 CDT 2021
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
$DATA ../../../../data/SD4/TD1/dat41.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1705.77390111686        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5033E+02 -4.3136E+01  1.6658E+00 -3.3155E+01  2.8201E+01  4.8172E+01  2.3343E+00  1.6564E+00  2.7923E+01 -9.8824E+00
             2.4940E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1714.37747505650        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0354E+00  1.0686E+00  9.4674E-01  1.0376E+00  9.7912E-01  9.6330E-01  1.0028E+00  9.8964E-01  8.2222E-01  1.0852E+00
             9.1989E-01
 PARAMETER:  1.3477E-01  1.6639E-01  4.5273E-02  1.3690E-01  7.8901E-02  6.2612E-02  1.0283E-01  8.9589E-02 -9.5747E-02  1.8173E-01
             1.6504E-02
 GRADIENT:   2.7623E+00  1.2387E+01  2.2029E+00  1.3921E+01 -7.5208E+00 -4.6258E+00 -6.8539E+00  3.4458E+00 -6.6538E+00  3.3477E+00
            -9.6617E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1715.47580612760        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0358E+00  8.9727E-01  8.2491E-01  1.1426E+00  8.4334E-01  9.4211E-01  1.3161E+00  6.1800E-01  7.1112E-01  9.8139E-01
             9.3307E-01
 PARAMETER:  1.3516E-01 -8.3938E-03 -9.2487E-02  2.3332E-01 -7.0382E-02  4.0365E-02  3.7468E-01 -3.8127E-01 -2.4092E-01  8.1214E-02
             3.0721E-02
 GRADIENT:   1.9046E+00  2.2240E+01 -1.2751E+01  5.0212E+01  1.1015E+01 -1.4148E+01  9.0979E-02  2.3813E+00 -9.7155E+00  1.1052E+00
            -1.5202E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1717.31767060665        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0341E+00  8.0008E-01  7.6786E-01  1.1713E+00  7.6620E-01  9.8387E-01  1.3882E+00  3.2816E-01  7.3714E-01  9.5110E-01
             9.3076E-01
 PARAMETER:  1.3353E-01 -1.2304E-01 -1.6415E-01  2.5811E-01 -1.6632E-01  8.3744E-02  4.2803E-01 -1.0143E+00 -2.0498E-01  4.9869E-02
             2.8243E-02
 GRADIENT:  -1.6577E+00  4.7385E+00  4.2404E+00  2.6234E+00 -6.7160E+00  3.2783E+00  1.9450E-01  4.2179E-01  2.5895E-01  7.5998E-02
            -7.1319E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1717.48003370497        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  1.0347E+00  7.2983E-01  7.5113E-01  1.2043E+00  7.3415E-01  9.7516E-01  1.4892E+00  2.0076E-01  7.1895E-01  9.4799E-01
             9.3060E-01
 PARAMETER:  1.3415E-01 -2.1494E-01 -1.8617E-01  2.8593E-01 -2.0904E-01  7.4844E-02  4.9822E-01 -1.5056E+00 -2.2996E-01  4.6584E-02
             2.8077E-02
 GRADIENT:   3.6088E-01 -9.5088E-01 -2.4992E-01 -3.0062E+00 -1.5045E-01 -1.2113E-01  3.2931E-01  1.7323E-01  2.2084E-01  3.4891E-01
             1.6393E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1717.48018577888        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      899
 NPARAMETR:  1.0346E+00  7.2037E-01  7.5118E-01  1.2100E+00  7.3093E-01  9.7506E-01  1.5015E+00  1.8537E-01  7.1679E-01  9.4920E-01
             9.3052E-01
 PARAMETER:  1.3404E-01 -2.2800E-01 -1.8611E-01  2.9060E-01 -2.1344E-01  7.4743E-02  5.0645E-01 -1.5854E+00 -2.3297E-01  4.7859E-02
             2.7992E-02
 GRADIENT:   2.5445E-01 -9.0653E-01 -3.6496E-01 -2.0203E+00  1.3621E-01 -1.4171E-01  1.8186E-01  1.3754E-01  1.3489E-01  2.6803E-01
             1.6766E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1717.48192632409        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  1.0345E+00  7.0971E-01  7.4987E-01  1.2162E+00  7.2655E-01  9.7502E-01  1.5148E+00  1.5616E-01  7.1446E-01  9.5049E-01
             9.3042E-01
 PARAMETER:  1.3391E-01 -2.4290E-01 -1.8786E-01  2.9570E-01 -2.1945E-01  7.4698E-02  5.1527E-01 -1.7569E+00 -2.3623E-01  4.9217E-02
             2.7877E-02
 GRADIENT:   9.3517E-02 -8.9374E-01 -5.6541E-01 -7.1729E-01  5.7136E-01 -1.4503E-01 -2.6128E-02  8.8378E-02  3.1922E-02  1.7054E-01
             1.7699E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1717.53692793119        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1254
 NPARAMETR:  1.0346E+00  7.2891E-01  7.4534E-01  1.2058E+00  7.3010E-01  9.7599E-01  1.4824E+00  1.6442E-02  7.1914E-01  9.4841E-01
             9.3017E-01
 PARAMETER:  1.3401E-01 -2.1620E-01 -1.9392E-01  2.8717E-01 -2.1457E-01  7.5698E-02  4.9364E-01 -4.0079E+00 -2.2970E-01  4.7032E-02
             2.7614E-02
 GRADIENT:  -1.4974E-01  5.4264E-01  3.2271E+00  1.0700E+00 -1.3718E+00  1.5971E-01 -3.7017E-01  3.7890E-04 -1.4296E-01 -1.3823E+00
            -5.6864E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1717.60165644579        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1434
 NPARAMETR:  1.0367E+00  7.7532E-01  7.2740E-01  1.1768E+00  7.3602E-01  9.7655E-01  1.4201E+00  1.0000E-02  7.2942E-01  9.3945E-01
             9.3058E-01
 PARAMETER:  1.3601E-01 -1.5448E-01 -2.1828E-01  2.6281E-01 -2.0650E-01  7.6275E-02  4.5071E-01 -6.0413E+00 -2.1551E-01  3.7539E-02
             2.8054E-02
 GRADIENT:   3.1882E+00  3.0296E-01  1.9341E-01 -2.4463E-01 -4.9231E-01  2.0101E-01  8.5337E-02  0.0000E+00  8.9132E-02  1.8949E-02
            -1.8427E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1717.60165644579        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1456
 NPARAMETR:  1.0367E+00  7.7532E-01  7.2740E-01  1.1768E+00  7.3602E-01  9.7655E-01  1.4201E+00  1.0000E-02  7.2942E-01  9.3945E-01
             9.3058E-01
 PARAMETER:  1.3601E-01 -1.5448E-01 -2.1828E-01  2.6281E-01 -2.0650E-01  7.6275E-02  4.5071E-01 -6.0413E+00 -2.1551E-01  3.7539E-02
             2.8054E-02
 GRADIENT:   3.1882E+00  3.0296E-01  1.9341E-01 -2.4463E-01 -4.9231E-01  2.0101E-01  8.5337E-02  0.0000E+00  8.9132E-02  1.8949E-02
            -1.8427E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1456
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.9542E-05  7.4263E-03 -5.3524E-04 -1.0482E-02 -5.4977E-03
 SE:             2.9845E-02  2.0596E-02  2.2057E-04  2.4403E-02  2.4288E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9894E-01  7.1842E-01  1.5238E-02  6.6752E-01  8.2093E-01

 ETASHRINKSD(%)  1.5239E-02  3.1000E+01  9.9261E+01  1.8248E+01  1.8631E+01
 ETASHRINKVR(%)  3.0475E-02  5.2391E+01  9.9995E+01  3.3166E+01  3.3791E+01
 EBVSHRINKSD(%)  3.8495E-01  3.1533E+01  9.9342E+01  1.8178E+01  1.6734E+01
 EBVSHRINKVR(%)  7.6842E-01  5.3123E+01  9.9996E+01  3.3052E+01  3.0668E+01
 RELATIVEINF(%)  9.8746E+01  3.4894E+00  5.1179E-04  5.9055E+00  5.7031E+00
 EPSSHRINKSD(%)  4.4466E+01
 EPSSHRINKVR(%)  6.9160E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1717.6016564457925     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -982.45082988205434     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.76
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1717.602       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  7.75E-01  7.27E-01  1.18E+00  7.36E-01  9.77E-01  1.42E+00  1.00E-02  7.29E-01  9.39E-01  9.31E-01
 


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
 #CPUT: Total CPU Time in Seconds,       42.568
Stop Time:
Sun Oct 24 03:46:53 CDT 2021
