Sun Oct 24 03:07:59 CDT 2021
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
$DATA ../../../../data/SD4/SL1/dat90.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m90.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1624.35042472422        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6393E+02 -3.3697E+00 -2.5945E+01  4.4769E+01  6.0858E+01  2.9910E+01  4.9406E-01  4.3208E+00 -5.2739E-01 -1.9781E+00
            -3.4505E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1627.77626292866        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8031E-01  1.0180E+00  1.0195E+00  1.0148E+00  9.8090E-01  1.0139E+00  1.0039E+00  9.8843E-01  1.0243E+00  9.9250E-01
             1.0825E+00
 PARAMETER:  8.0110E-02  1.1780E-01  1.1927E-01  1.1469E-01  8.0711E-02  1.1384E-01  1.0386E-01  8.8362E-02  1.2401E-01  9.2472E-02
             1.7930E-01
 GRADIENT:  -3.3293E-01 -1.9788E+00 -8.3204E+00  7.0149E+00  1.1205E+01 -4.2593E-02  6.4100E-02  2.8392E+00 -1.8209E+00  1.1420E+00
             4.7507E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1627.90912739803        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      355
 NPARAMETR:  9.8223E-01  1.0578E+00  9.7085E-01  9.8663E-01  9.7612E-01  1.0000E+00  9.6889E-01  8.3521E-01  1.0437E+00  9.8293E-01
             1.0852E+00
 PARAMETER:  8.2070E-02  1.5622E-01  7.0417E-02  8.6542E-02  7.5828E-02  9.9998E-02  6.8395E-02 -8.0074E-02  1.4281E-01  8.2779E-02
             1.8172E-01
 GRADIENT:   2.9022E+00 -3.6071E-01 -2.0671E+00  5.6381E+00  9.8115E+00 -5.6592E+00 -1.5898E+00 -9.5154E-01 -2.7867E+00 -1.4056E+00
             3.5552E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1628.21371051311        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.8156E-01  1.1845E+00  8.4487E-01  8.9826E-01  9.6576E-01  1.0198E+00  9.3953E-01  7.5819E-01  1.1142E+00  9.5692E-01
             1.0811E+00
 PARAMETER:  8.1387E-02  2.6928E-01 -6.8575E-02 -7.2921E-03  6.5157E-02  1.1964E-01  3.7621E-02 -1.7682E-01  2.0814E-01  5.5966E-02
             1.7795E-01
 GRADIENT:  -8.7705E-01  1.7328E+00  1.7313E+00  8.7832E-01 -4.3694E+00  1.5967E+00  5.9327E-01  6.4219E-02  3.2239E-01  6.8600E-01
            -8.4802E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1628.43505398068        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  9.8372E-01  1.4190E+00  6.8080E-01  7.4857E-01  1.0041E+00  1.0120E+00  8.4537E-01  6.3236E-01  1.2536E+00  9.4138E-01
             1.0832E+00
 PARAMETER:  8.3590E-02  4.4993E-01 -2.8448E-01 -1.8960E-01  1.0405E-01  1.1195E-01 -6.7982E-02 -3.5830E-01  3.2600E-01  3.9589E-02
             1.7990E-01
 GRADIENT:   1.2655E+00  1.0579E+01  8.7429E-01  8.5370E+00 -2.5918E+00 -1.9195E+00 -6.2592E-01  2.6793E-01 -3.9665E-01 -1.1671E+00
            -3.9235E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1628.57393425620        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.8392E-01  1.6033E+00  5.4876E-01  6.2371E-01  1.0424E+00  1.0138E+00  7.9208E-01  4.5689E-01  1.3952E+00  9.4868E-01
             1.0828E+00
 PARAMETER:  8.3785E-02  5.7208E-01 -5.0009E-01 -3.7207E-01  1.4150E-01  1.1366E-01 -1.3309E-01 -6.8331E-01  4.3306E-01  4.7320E-02
             1.7955E-01
 GRADIENT:   8.0219E-01  9.6936E+00 -4.3801E-01  6.7492E+00 -1.2361E+00 -1.4059E+00 -6.1723E-01  3.7111E-01 -9.2814E-01 -5.9005E-01
            -1.6802E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1628.65340865484        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  9.8421E-01  1.6445E+00  5.1700E-01  5.7693E-01  1.0617E+00  1.0191E+00  7.8084E-01  2.1976E-01  1.4619E+00  9.6162E-01
             1.0834E+00
 PARAMETER:  8.4079E-02  5.9744E-01 -5.5971E-01 -4.5004E-01  1.5988E-01  1.1889E-01 -1.4738E-01 -1.4152E+00  4.7973E-01  6.0861E-02
             1.8010E-01
 GRADIENT:   1.8609E+00 -1.9813E+01  1.0346E+00 -9.6501E+00  1.1132E+00  7.4431E-01 -1.1795E-01  6.2569E-02 -6.7109E-01 -1.8663E-01
             4.2852E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1628.74267746327        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1239
 NPARAMETR:  9.8320E-01  1.6657E+00  5.0651E-01  5.7671E-01  1.0620E+00  1.0171E+00  7.7942E-01  1.2663E-01  1.4711E+00  9.6430E-01
             1.0823E+00
 PARAMETER:  8.3062E-02  6.1022E-01 -5.8021E-01 -4.5041E-01  1.6019E-01  1.1695E-01 -1.4921E-01 -1.9665E+00  4.8603E-01  6.3646E-02
             1.7907E-01
 GRADIENT:  -6.7156E-01  2.3950E+00  8.0006E-01  2.1191E+00 -1.0608E+00 -1.1286E-01  1.3497E-01  2.1482E-02 -3.3592E-02  7.7334E-02
            -1.6940E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1628.76344493323        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1420
 NPARAMETR:  9.8413E-01  1.6736E+00  4.9684E-01  5.6755E-01  1.0637E+00  1.0176E+00  7.7547E-01  3.0022E-02  1.4813E+00  9.6152E-01
             1.0824E+00
 PARAMETER:  8.3998E-02  6.1498E-01 -5.9949E-01 -4.6643E-01  1.6173E-01  1.1745E-01 -1.5429E-01 -3.4058E+00  4.9290E-01  6.0761E-02
             1.7919E-01
 GRADIENT:   1.3739E+00 -4.0607E+00  3.4470E-01 -5.7446E-01  2.6120E-01  1.0088E-01 -1.1581E-01  1.5368E-03 -1.8016E-01 -8.0353E-02
            -3.1927E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1628.76504998126        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1596
 NPARAMETR:  9.8413E-01  1.6708E+00  4.9411E-01  5.6889E-01  1.0631E+00  1.0176E+00  7.7660E-01  1.0000E-02  1.4837E+00  9.6143E-01
             1.0828E+00
 PARAMETER:  8.4012E-02  6.1451E-01 -5.9990E-01 -4.6594E-01  1.6064E-01  1.1749E-01 -1.5389E-01 -5.3642E+00  4.9319E-01  5.9667E-02
             1.7894E-01
 GRADIENT:   1.2540E-02  1.3284E+00  5.2732E-01 -5.3692E-01 -3.6369E-01  2.6847E-03 -8.0745E-02  0.0000E+00 -1.1200E-01 -6.7923E-02
            -1.1546E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1596
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1042E-04 -3.2786E-02 -2.7949E-04  2.5862E-02 -3.6479E-02
 SE:             2.9830E-02  2.3088E-02  1.0746E-04  2.3584E-02  2.2142E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9437E-01  1.5560E-01  9.2951E-03  2.7282E-01  9.9456E-02

 ETASHRINKSD(%)  6.4955E-02  2.2652E+01  9.9640E+01  2.0991E+01  2.5821E+01
 ETASHRINKVR(%)  1.2987E-01  4.0173E+01  9.9999E+01  3.7576E+01  4.4975E+01
 EBVSHRINKSD(%)  4.8540E-01  2.2197E+01  9.9692E+01  2.2203E+01  2.4605E+01
 EBVSHRINKVR(%)  9.6844E-01  3.9467E+01  9.9999E+01  3.9477E+01  4.3156E+01
 RELATIVEINF(%)  9.9005E+01  4.4606E+00  1.2034E-04  4.8308E+00  1.1069E+01
 EPSSHRINKSD(%)  4.4244E+01
 EPSSHRINKVR(%)  6.8913E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1628.7650499812596     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -893.61422341752143     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1628.765       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  1.67E+00  4.97E-01  5.68E-01  1.06E+00  1.02E+00  7.76E-01  1.00E-02  1.48E+00  9.60E-01  1.08E+00
 


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
 #CPUT: Total CPU Time in Seconds,       47.541
Stop Time:
Sun Oct 24 03:08:09 CDT 2021
