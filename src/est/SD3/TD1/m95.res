Sun Oct 24 00:36:45 CDT 2021
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
$DATA ../../../../data/SD3/TD1/dat95.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m95.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2106.07036695821        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9255E+02 -4.6056E+01 -4.2186E+01 -4.8666E+01  6.9134E+01  4.4317E+01 -8.3669E+00  1.9192E+01 -6.3900E+01  9.3042E+00
             3.7426E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2121.42067097978        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  1.0045E+00  9.5383E-01  1.1018E+00  1.0903E+00  9.7421E-01  9.7163E-01  1.0276E+00  9.0074E-01  1.3385E+00  9.4059E-01
             9.1672E-01
 PARAMETER:  1.0448E-01  5.2726E-02  1.9693E-01  1.8644E-01  7.3875E-02  7.1222E-02  1.2718E-01 -4.5422E-03  3.9152E-01  3.8750E-02
             1.3051E-02
 GRADIENT:   7.0776E+00 -1.0775E+01 -2.9943E+00 -1.2319E+01  1.0633E+01 -3.4963E+00  8.2868E+00  1.2100E+01  2.5063E+01 -4.8054E+00
            -3.3414E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2126.06197125100        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      347
 NPARAMETR:  1.0014E+00  8.5348E-01  9.7616E-01  1.1609E+00  8.8762E-01  1.0009E+00  9.4321E-01  3.5713E-01  1.3036E+00  1.0745E+00
             9.1627E-01
 PARAMETER:  1.0135E-01 -5.8429E-02  7.5868E-02  2.4918E-01 -1.9210E-02  1.0088E-01  4.1535E-02 -9.2967E-01  3.6515E-01  1.7183E-01
             1.2558E-02
 GRADIENT:   1.0457E+00 -1.0598E+01 -3.8379E+01  7.5184E+00  9.2834E+00  8.6375E+00  5.1965E+00  2.8199E+00  2.8949E+01  1.7577E+01
            -2.6514E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2129.94861839103        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  1.0007E+00  7.7261E-01  1.0605E+00  1.2216E+00  8.9100E-01  9.7783E-01  1.0007E+00  3.3246E-01  1.1266E+00  1.0391E+00
             9.5028E-01
 PARAMETER:  1.0073E-01 -1.5798E-01  1.5878E-01  3.0018E-01 -1.5406E-02  7.7577E-02  1.0066E-01 -1.0012E+00  2.1920E-01  1.3839E-01
             4.9000E-02
 GRADIENT:   1.8799E+00  7.7454E+00 -4.7586E-01  9.0644E+00 -2.8422E+00  4.7647E-01  1.6462E+00  1.4513E+00 -2.8992E+00  2.0126E+00
             2.1237E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2131.40727177812        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  9.9847E-01  6.7681E-01  1.0881E+00  1.2652E+00  8.8002E-01  9.7343E-01  3.1646E-01  1.3823E-01  1.1453E+00  1.1182E+00
             9.4774E-01
 PARAMETER:  9.8474E-02 -2.9037E-01  1.8440E-01  3.3521E-01 -2.7810E-02  7.3074E-02 -1.0506E+00 -1.8788E+00  2.3568E-01  2.1168E-01
             4.6326E-02
 GRADIENT:  -1.4801E-02 -1.1343E+00 -8.5718E-01  6.0652E-01  3.6961E-01 -4.1198E-01  2.4095E-01  1.9241E-01  1.0924E+00  1.1073E+00
            -2.8852E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2131.64456199538        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.9881E-01  7.1458E-01  1.0867E+00  1.2397E+00  8.9376E-01  9.7370E-01  3.8431E-02  3.8215E-02  1.1699E+00  1.1276E+00
             9.5118E-01
 PARAMETER:  9.8808E-02 -2.3606E-01  1.8311E-01  3.1488E-01 -1.2317E-02  7.3350E-02 -3.1589E+00 -3.1645E+00  2.5695E-01  2.2010E-01
             4.9945E-02
 GRADIENT:  -4.4234E-01  6.2494E-01 -3.1054E-01  3.2685E-01  9.9054E-01 -4.9125E-01  5.2241E-03  1.4857E-02 -9.0588E-01 -1.2445E-02
            -1.6674E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2131.64596376319        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  9.9886E-01  7.0513E-01  1.0865E+00  1.2453E+00  8.8983E-01  9.7461E-01  2.8422E-02  3.0868E-02  1.1671E+00  1.1262E+00
             9.5144E-01
 PARAMETER:  9.8862E-02 -2.4938E-01  1.8293E-01  3.1934E-01 -1.6721E-02  7.4286E-02 -3.4606E+00 -3.3780E+00  2.5452E-01  2.1885E-01
             5.0221E-02
 GRADIENT:  -3.7153E-02 -9.4307E-02 -7.7359E-02 -1.8109E-01  2.7561E-01 -6.9251E-02  3.3182E-03  9.7465E-03  1.2973E-01  8.5870E-02
             6.2262E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2131.64685067888        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1231
 NPARAMETR:  9.9885E-01  7.0159E-01  1.0862E+00  1.2475E+00  8.8825E-01  9.7481E-01  1.2387E-02  1.7801E-02  1.1649E+00  1.1253E+00
             9.5145E-01
 PARAMETER:  9.8850E-02 -2.5440E-01  1.8267E-01  3.2116E-01 -1.8502E-02  7.4491E-02 -4.2911E+00 -3.9285E+00  2.5267E-01  2.1807E-01
             5.0227E-02
 GRADIENT:   3.4860E-02 -9.4979E-02 -6.1844E-03 -8.8736E-02 -1.8691E-02  3.0691E-02  8.1845E-04  3.2674E-03  1.3492E-01  3.5900E-02
             3.2925E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -2131.64792560741        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1325
 NPARAMETR:  9.9884E-01  7.0179E-01  1.0861E+00  1.2475E+00  8.8826E-01  9.7474E-01  1.0000E-02  1.0000E-02  1.1644E+00  1.1251E+00
             9.5141E-01
 PARAMETER:  9.8837E-02 -2.5412E-01  1.8258E-01  3.2114E-01 -1.8487E-02  7.4412E-02 -5.7021E+00 -4.8510E+00  2.5222E-01  2.1786E-01
             5.0188E-02
 GRADIENT:  -6.7645E-03  5.0477E-02 -4.3880E-03  5.3202E-02 -2.0413E-02 -2.0752E-03  0.0000E+00  0.0000E+00 -8.1041E-02 -2.6103E-02
            -1.9203E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1325
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4010E-03 -3.7844E-04 -4.0090E-04 -2.7072E-03 -1.7185E-02
 SE:             2.9910E-02  1.3708E-04  2.1368E-04  2.9716E-02  2.6303E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6264E-01  5.7676E-03  6.0633E-02  9.2741E-01  5.1354E-01

 ETASHRINKSD(%)  1.0000E-10  9.9541E+01  9.9284E+01  4.4676E-01  1.1880E+01
 ETASHRINKVR(%)  1.0000E-10  9.9998E+01  9.9995E+01  8.9152E-01  2.2349E+01
 EBVSHRINKSD(%)  3.1221E-01  9.9594E+01  9.9318E+01  9.3502E-01  8.6433E+00
 EBVSHRINKVR(%)  6.2344E-01  9.9998E+01  9.9995E+01  1.8613E+00  1.6540E+01
 RELATIVEINF(%)  9.8922E+01  2.1940E-04  1.6895E-03  1.9031E+01  1.8121E+01
 EPSSHRINKSD(%)  3.2612E+01
 EPSSHRINKVR(%)  5.4589E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2131.6479256074113     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1212.7093924027386     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.33
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2131.648       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  7.02E-01  1.09E+00  1.25E+00  8.88E-01  9.75E-01  1.00E-02  1.00E-02  1.16E+00  1.13E+00  9.51E-01
 


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
 #CPUT: Total CPU Time in Seconds,       41.478
Stop Time:
Sun Oct 24 00:36:54 CDT 2021
