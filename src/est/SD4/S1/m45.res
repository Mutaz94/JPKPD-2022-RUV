Sun Oct 24 02:45:25 CDT 2021
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
$DATA ../../../../data/SD4/S1/dat45.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1615.27452180136        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1534E+02 -3.1782E+00 -3.7837E+01  6.6039E+01  7.8601E+01  4.5438E+00  6.1136E+00  6.7657E+00  1.1448E+01  3.7691E+00
            -3.7933E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1625.65651615821        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8648E-01  1.0241E+00  1.0309E+00  1.0107E+00  9.7554E-01  1.1647E+00  9.8544E-01  9.8162E-01  9.9050E-01  9.7096E-01
             1.0900E+00
 PARAMETER:  8.6391E-02  1.2377E-01  1.3048E-01  1.1062E-01  7.5233E-02  2.5247E-01  8.5335E-02  8.1452E-02  9.0459E-02  7.0527E-02
             1.8614E-01
 GRADIENT:   2.4962E+01  4.0694E-01 -9.2520E+00  1.1329E+01  1.4098E+01 -3.2139E-01  3.4568E+00  4.2557E+00  1.5128E+00  4.9778E+00
            -6.7095E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1626.72363432898        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.8487E-01  1.0014E+00  9.0318E-01  1.0144E+00  9.0826E-01  1.1684E+00  9.3294E-01  6.3788E-01  9.8699E-01  9.2603E-01
             1.0855E+00
 PARAMETER:  8.4754E-02  1.0137E-01 -1.8297E-03  1.1433E-01  3.7716E-03  2.5563E-01  3.0590E-02 -3.4961E-01  8.6910E-02  2.3146E-02
             1.8204E-01
 GRADIENT:   2.0795E+01 -4.8244E+00 -1.0626E+01  9.8853E+00  1.5431E+01  6.7895E-01 -2.2278E+00  1.5517E+00 -5.2981E-01  1.5160E+00
            -3.0396E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1627.48107669440        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.6858E-01  1.1343E+00  7.8726E-01  9.1887E-01  9.0358E-01  1.1674E+00  9.2327E-01  3.6471E-01  1.0451E+00  9.0127E-01
             1.0968E+00
 PARAMETER:  6.8080E-02  2.2606E-01 -1.3920E-01  1.5392E-02 -1.3880E-03  2.5476E-01  2.0166E-02 -9.0865E-01  1.4410E-01 -3.9450E-03
             1.9237E-01
 GRADIENT:  -8.0316E+00 -1.6603E+00 -5.1510E-01 -1.5289E+00 -2.3201E+00  9.2364E-02  5.7945E-01  4.8799E-01  1.0885E+00  1.7072E+00
             1.3337E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1627.71762088928        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  9.7388E-01  1.2884E+00  7.1297E-01  8.2034E-01  9.4393E-01  1.1668E+00  8.4641E-01  1.1161E-01  1.1317E+00  9.0725E-01
             1.0962E+00
 PARAMETER:  7.3533E-02  3.5343E-01 -2.3831E-01 -9.8039E-02  4.2296E-02  2.5429E-01 -6.6750E-02 -2.0928E+00  2.2376E-01  2.6626E-03
             1.9184E-01
 GRADIENT:  -8.2231E-02  6.0628E-02 -9.3428E-01  1.0703E+00  5.5988E-01 -3.7686E-01 -6.4316E-02  4.7790E-02 -3.7725E-03  1.9875E-01
             1.8329E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1627.74274432052        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      868
 NPARAMETR:  9.7457E-01  1.2939E+00  7.1295E-01  8.1604E-01  9.4670E-01  1.1742E+00  8.4434E-01  3.8408E-02  1.1369E+00  9.0856E-01
             1.0960E+00
 PARAMETER:  7.4238E-02  3.5764E-01 -2.3834E-01 -1.0329E-01  4.5230E-02  2.6061E-01 -6.9206E-02 -3.1595E+00  2.2829E-01  4.1021E-03
             1.9163E-01
 GRADIENT:   1.1101E+00 -2.6232E-01  1.5713E-01 -4.4010E-01 -3.5522E-01  2.2153E+00  3.9408E-03  5.3963E-03  9.6623E-04 -1.2064E-01
            -5.9077E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1627.74575217291        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1047             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7484E-01  1.2934E+00  7.1337E-01  8.1650E-01  9.4662E-01  1.1733E+00  8.4430E-01  1.0000E-02  1.1368E+00  9.0946E-01
             1.0960E+00
 PARAMETER:  7.4516E-02  3.5727E-01 -2.3776E-01 -1.0273E-01  4.5141E-02  2.5980E-01 -6.9245E-02 -5.1854E+00  2.2818E-01  5.1003E-03
             1.9163E-01
 GRADIENT:   3.7678E+02  2.0093E+02  3.3752E+00  4.4530E+01  6.5331E+00  1.8711E+02  3.4417E+00  0.0000E+00  1.1413E+01  5.1648E-01
             1.4865E+00

0ITERATION NO.:   32    OBJECTIVE VALUE:  -1627.74575217291        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:     1102
 NPARAMETR:  9.7484E-01  1.2934E+00  7.1337E-01  8.1650E-01  9.4662E-01  1.1733E+00  8.4430E-01  1.0000E-02  1.1368E+00  9.0946E-01
             1.0960E+00
 PARAMETER:  7.4516E-02  3.5727E-01 -2.3776E-01 -1.0273E-01  4.5141E-02  2.5980E-01 -6.9245E-02 -5.1854E+00  2.2818E-01  5.1003E-03
             1.9163E-01
 GRADIENT:   1.5570E+00 -6.5188E-02  1.8767E-01 -2.9461E-01 -5.0956E-01  1.8806E+00  4.5534E-03  0.0000E+00  2.4482E-02 -2.7794E-02
            -4.6611E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1102
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.6750E-05 -2.1393E-02 -3.1859E-04  1.0578E-02 -2.5396E-02
 SE:             2.9853E-02  2.0332E-02  1.3999E-04  2.4584E-02  2.3117E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9795E-01  2.9270E-01  2.2859E-02  6.6700E-01  2.7194E-01

 ETASHRINKSD(%)  1.0000E-10  3.1886E+01  9.9531E+01  1.7641E+01  2.2556E+01
 ETASHRINKVR(%)  1.0000E-10  5.3605E+01  9.9998E+01  3.2170E+01  4.0024E+01
 EBVSHRINKSD(%)  3.8132E-01  3.1279E+01  9.9568E+01  1.8097E+01  2.1732E+01
 EBVSHRINKVR(%)  7.6119E-01  5.2774E+01  9.9998E+01  3.2920E+01  3.8741E+01
 RELATIVEINF(%)  9.9118E+01  1.9827E+00  2.0131E-04  3.5807E+00  7.2071E+00
 EPSSHRINKSD(%)  4.3179E+01
 EPSSHRINKVR(%)  6.7714E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1627.7457521729075     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -892.59492560916931     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1627.746       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.29E+00  7.13E-01  8.17E-01  9.47E-01  1.17E+00  8.44E-01  1.00E-02  1.14E+00  9.09E-01  1.10E+00
 


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
 #CPUT: Total CPU Time in Seconds,       33.127
Stop Time:
Sun Oct 24 02:45:33 CDT 2021
