Sun Oct 24 04:19:08 CDT 2021
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
$DATA ../../../../data/SD4/D/dat43.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1672.66671151795        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4374E+02  4.2191E+00 -1.6544E+01  4.1676E+01  2.1056E+01  2.1712E+01  1.0976E+01  1.0451E+01  3.3519E+01  4.5904E+00
             3.8875E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1678.14024214818        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.7918E-01  1.0222E+00  1.0343E+00  1.0090E+00  1.0152E+00  1.0803E+00  9.6655E-01  9.5657E-01  8.8854E-01  9.9057E-01
             8.8118E-01
 PARAMETER:  7.8958E-02  1.2195E-01  1.3371E-01  1.0900E-01  1.1507E-01  1.7724E-01  6.5973E-02  5.5596E-02 -1.8172E-02  9.0529E-02
            -2.6494E-02
 GRADIENT:  -1.5474E+00  6.2127E+00 -5.5459E-01  4.4156E+00  7.5563E+00  8.3576E+00 -4.2281E-01  4.8072E+00  3.2134E+00 -6.4168E+00
            -1.4638E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1680.54347540541        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.8733E-01  8.2843E-01  1.0972E+00  1.1514E+00  9.6578E-01  1.0327E+00  1.1495E+00  4.9174E-01  7.7752E-01  1.0429E+00
             9.0403E-01
 PARAMETER:  8.7249E-02 -8.8223E-02  1.9276E-01  2.4095E-01  6.5180E-02  1.3215E-01  2.3934E-01 -6.0981E-01 -1.5165E-01  1.4199E-01
            -8.8879E-04
 GRADIENT:   1.8345E+01  3.0288E+01  1.3628E+01  4.4978E+01  4.5488E-02 -9.6291E+00 -2.6999E+00 -1.0018E+00 -5.7598E+00 -7.6979E+00
            -1.1385E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1684.02196435366        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.7768E-01  6.6783E-01  8.2964E-01  1.2158E+00  7.4629E-01  1.0622E+00  1.5333E+00  2.3530E-01  6.9464E-01  8.6577E-01
             9.1027E-01
 PARAMETER:  7.7426E-02 -3.0372E-01 -8.6761E-02  2.9539E-01 -1.9265E-01  1.6038E-01  5.2742E-01 -1.3469E+00 -2.6436E-01 -4.4137E-02
             5.9844E-03
 GRADIENT:  -4.5984E+00  2.5631E+01  1.8394E+01  2.3000E+01 -3.1117E+01  1.3481E+00  2.2386E+00  2.7635E-01 -3.6937E+00  1.0471E+00
            -2.1319E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1686.70652688754        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.7692E-01  3.8576E-01  8.7903E-01  1.3745E+00  6.9564E-01  1.0456E+00  2.0176E+00  2.3245E-02  6.6859E-01  8.9554E-01
             9.2160E-01
 PARAMETER:  7.6650E-02 -8.5253E-01 -2.8931E-02  4.1806E-01 -2.6292E-01  1.4461E-01  8.0191E-01 -3.6617E+00 -3.0258E-01 -1.0330E-02
             1.8351E-02
 GRADIENT:   2.6042E+00  8.7658E+00  1.1275E+01  2.0209E+01 -1.5915E+01 -3.4886E+00 -1.6160E+00  3.5379E-04 -1.6172E+00 -2.1492E+00
             6.9818E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1688.20834710466        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.7038E-01  1.6513E-01  9.3348E-01  1.5001E+00  6.7329E-01  1.0557E+00  3.3902E+00  1.0000E-02  6.4293E-01  9.4397E-01
             9.1986E-01
 PARAMETER:  6.9935E-02 -1.7010E+00  3.1160E-02  5.0551E-01 -2.9558E-01  1.5422E-01  1.3209E+00 -8.6110E+00 -3.4171E-01  4.2341E-02
             1.6467E-02
 GRADIENT:  -2.4119E-01  3.6083E+00  1.0146E+01  1.9425E+01 -1.4911E+01  2.0780E+00  2.2506E-01  0.0000E+00  2.5225E+00  5.8399E-01
            -1.1335E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1688.99366099056        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1055
 NPARAMETR:  9.6821E-01  7.3564E-02  9.3475E-01  1.5400E+00  6.5922E-01  1.0485E+00  4.9076E+00  1.0000E-02  6.2051E-01  9.4717E-01
             9.2466E-01
 PARAMETER:  6.7695E-02 -2.5096E+00  3.2520E-02  5.3178E-01 -3.1669E-01  1.4732E-01  1.6908E+00 -1.4038E+01 -3.7721E-01  4.5720E-02
             2.1675E-02
 GRADIENT:  -2.2705E-01  6.2889E-01  1.0148E+00 -2.8374E+00 -2.0448E+00 -9.5834E-02  1.3066E+00  0.0000E+00 -1.6245E+00  7.5024E-01
             1.5902E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1689.13747903480        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  9.6799E-01  5.1097E-02  9.3684E-01  1.5530E+00  6.5857E-01  1.0485E+00  5.6127E+00  1.0000E-02  6.2305E-01  9.2939E-01
             9.2559E-01
 PARAMETER:  6.7464E-02 -2.8740E+00  3.4758E-02  5.4018E-01 -3.1768E-01  1.4735E-01  1.8250E+00 -1.6532E+01 -3.7314E-01  2.6769E-02
             2.2681E-02
 GRADIENT:   2.4068E-01 -1.8461E+00 -1.4792E+00  1.5220E+01  2.6852E+00  6.9255E-02 -4.2459E+00  0.0000E+00  4.6931E+00 -1.6478E+00
            -4.8628E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1689.14841764060        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1299
 NPARAMETR:  9.6806E-01  5.3000E-02  9.3703E-01  1.5514E+00  6.5833E-01  1.0486E+00  5.5369E+00  1.0000E-02  6.2190E-01  9.3429E-01
             9.2569E-01
 PARAMETER:  6.7455E-02 -2.8365E+00  3.4680E-02  5.3911E-01 -3.1787E-01  1.4736E-01  1.8117E+00 -1.6275E+01 -3.7489E-01  3.1028E-02
             2.2610E-02
 GRADIENT:  -6.8784E-01  3.5383E+00 -1.7259E+00 -1.0736E+01  4.8999E+00 -2.5328E-01  3.6397E+00  0.0000E+00  9.1086E-01 -1.9254E+00
            -3.0509E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1299
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.2341E-04  2.5745E-02 -3.0606E-04 -2.1941E-02 -8.8586E-03
 SE:             2.9898E-02  1.2273E-02  2.2897E-04  2.7529E-02  2.5824E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7803E-01  3.5940E-02  1.8132E-01  4.2544E-01  7.3157E-01

 ETASHRINKSD(%)  1.0000E-10  5.8883E+01  9.9233E+01  7.7736E+00  1.3487E+01
 ETASHRINKVR(%)  1.0000E-10  8.3094E+01  9.9994E+01  1.4943E+01  2.5156E+01
 EBVSHRINKSD(%)  3.1048E-01  7.2387E+01  9.9271E+01  5.2626E+00  9.9193E+00
 EBVSHRINKVR(%)  6.2000E-01  9.2375E+01  9.9995E+01  1.0248E+01  1.8855E+01
 RELATIVEINF(%)  9.9201E+01  3.7323E+00  4.4400E-04  4.0127E+01  6.9123E+00
 EPSSHRINKSD(%)  4.3333E+01
 EPSSHRINKVR(%)  6.7888E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1689.1484176406029     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -953.99759107686475     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1689.148       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  5.31E-02  9.37E-01  1.55E+00  6.58E-01  1.05E+00  5.54E+00  1.00E-02  6.22E-01  9.33E-01  9.26E-01
 


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
 #CPUT: Total CPU Time in Seconds,       39.088
Stop Time:
Sun Oct 24 04:19:17 CDT 2021
