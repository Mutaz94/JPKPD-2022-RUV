Sun Oct 24 03:34:40 CDT 2021
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
$DATA ../../../../data/SD4/SL3/dat63.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m63.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1651.37134276002        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.1466E+02 -1.0844E+02 -6.7187E+01 -7.1893E+01  9.6365E+01  1.4573E+01 -1.7510E+01  1.3249E+01 -1.6381E+01  6.5154E+00
            -3.1866E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1672.43438997602        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  1.0536E+00  1.0930E+00  1.1187E+00  1.0547E+00  1.0076E+00  1.0913E+00  1.0476E+00  9.5228E-01  1.0537E+00  9.5584E-01
             1.1022E+00
 PARAMETER:  1.5225E-01  1.8892E-01  2.1215E-01  1.5327E-01  1.0755E-01  1.8740E-01  1.4649E-01  5.1105E-02  1.5227E-01  5.4837E-02
             1.9734E-01
 GRADIENT:   3.8431E+00 -8.6884E+00 -1.1987E+01 -1.6195E+00  3.2548E+00  7.5081E+00 -5.0528E+00  5.5304E+00 -1.9012E+00 -1.5765E-01
             9.7967E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1673.27296474869        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0488E+00  1.0972E+00  1.2031E+00  1.0558E+00  1.0500E+00  1.0661E+00  1.1121E+00  8.5163E-01  1.0399E+00  1.0257E+00
             1.0855E+00
 PARAMETER:  1.4766E-01  1.9280E-01  2.8494E-01  1.5431E-01  1.4876E-01  1.6405E-01  2.0629E-01 -6.0601E-02  1.3916E-01  1.2536E-01
             1.8201E-01
 GRADIENT:  -4.0704E+00 -2.1903E+00 -3.2646E-01 -4.7192E+00  4.4639E+00 -1.6035E+00 -7.2295E-01  6.2063E-02 -3.1284E+00  2.2854E+00
             1.9534E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1673.87999517643        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  1.0525E+00  9.5194E-01  9.9830E-01  1.1373E+00  8.9793E-01  1.0721E+00  1.3676E+00  5.1945E-01  9.4866E-01  8.7516E-01
             1.0798E+00
 PARAMETER:  1.5112E-01  5.0743E-02  9.8295E-02  2.2868E-01 -7.6641E-03  1.6960E-01  4.1309E-01 -5.5499E-01  4.7293E-02 -3.3345E-02
             1.7677E-01
 GRADIENT:   1.9456E+00  2.9990E+00 -1.2210E+00  3.6057E+00  1.5658E-01  8.9460E-02  7.1312E-01  8.4164E-01 -2.3404E-01  3.0608E-01
            -2.0205E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1674.17494905511        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  1.0511E+00  7.5418E-01  8.5701E-01  1.2324E+00  7.5550E-01  1.0720E+00  1.7050E+00  2.4789E-01  8.6216E-01  7.4285E-01
             1.0833E+00
 PARAMETER:  1.4984E-01 -1.8212E-01 -5.4310E-02  3.0895E-01 -1.8038E-01  1.6953E-01  6.3357E-01 -1.2948E+00 -4.8318E-02 -1.9726E-01
             1.8004E-01
 GRADIENT:  -5.5010E-01  9.2160E-01  1.6393E-01  1.3835E+00 -2.0509E+00 -1.2612E-01  1.5873E-01  4.1286E-01 -4.5206E-01 -1.3911E-01
             7.1723E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1674.34860430684        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  1.0520E+00  8.0964E-01  8.4750E-01  1.2016E+00  7.7022E-01  1.0725E+00  1.6099E+00  1.0039E-01  8.7825E-01  7.6917E-01
             1.0826E+00
 PARAMETER:  1.5072E-01 -1.1117E-01 -6.5460E-02  2.8364E-01 -1.6109E-01  1.7002E-01  5.7618E-01 -2.1987E+00 -2.9828E-02 -1.6245E-01
             1.7940E-01
 GRADIENT:   1.5337E-01  2.0221E+00  6.0682E-01  2.4259E+00 -2.7181E+00 -2.9462E-02  3.6333E-01  4.8638E-02 -2.5515E-01  1.6121E+00
             6.2946E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1674.39662830313        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1056
 NPARAMETR:  1.0521E+00  7.9427E-01  8.2575E-01  1.2043E+00  7.5436E-01  1.0729E+00  1.6408E+00  3.3036E-02  8.7307E-01  7.3811E-01
             1.0835E+00
 PARAMETER:  1.5078E-01 -1.3033E-01 -9.1468E-02  2.8594E-01 -1.8189E-01  1.7035E-01  5.9521E-01 -3.3102E+00 -3.5744E-02 -2.0366E-01
             1.8016E-01
 GRADIENT:   1.0216E-01  5.6018E-02 -1.1011E-01  1.2114E-02 -8.8325E-02  1.6901E-02  6.2573E-02  6.6320E-03  2.0260E-02  1.0576E-01
             5.8376E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1674.40270647376        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1241             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0534E+00  7.9616E-01  8.2429E-01  1.2027E+00  7.5468E-01  1.0735E+00  1.6412E+00  1.0000E-02  8.7374E-01  7.3755E-01
             1.0837E+00
 PARAMETER:  1.5200E-01 -1.2795E-01 -9.3233E-02  2.8455E-01 -1.8146E-01  1.7089E-01  5.9543E-01 -4.8318E+00 -3.4971E-02 -2.0442E-01
             1.8039E-01
 GRADIENT:   5.3277E+02  2.2111E+01  2.0788E+00  2.2711E+02  1.7261E+01  5.9043E+01  2.5374E+01  0.0000E+00  6.2486E+00  1.1688E+00
             1.6593E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1674.40291197472        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1410
 NPARAMETR:  1.0530E+00  7.9645E-01  8.2421E-01  1.2025E+00  7.5461E-01  1.0732E+00  1.6385E+00  1.0000E-02  8.7367E-01  7.3679E-01
             1.0834E+00
 PARAMETER:  1.5160E-01 -1.2759E-01 -9.3331E-02  2.8444E-01 -1.8155E-01  1.7067E-01  5.9377E-01 -4.8318E+00 -3.5052E-02 -2.0546E-01
             1.8011E-01
 GRADIENT:   1.6504E+00 -2.6772E-01 -5.1838E-01 -3.2880E-01  6.9463E-01  1.3128E-01  1.2848E-01  0.0000E+00  4.5353E-02  4.0573E-02
             3.8939E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1410
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.2061E-04  1.3852E-02 -4.7733E-04 -1.4514E-02 -3.0186E-03
 SE:             2.9836E-02  2.1711E-02  2.0443E-04  2.4657E-02  2.1637E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8073E-01  5.2347E-01  1.9546E-02  5.5611E-01  8.8905E-01

 ETASHRINKSD(%)  4.5046E-02  2.7265E+01  9.9315E+01  1.7395E+01  2.7513E+01
 ETASHRINKVR(%)  9.0071E-02  4.7097E+01  9.9995E+01  3.1765E+01  4.7456E+01
 EBVSHRINKSD(%)  4.4236E-01  2.7528E+01  9.9341E+01  1.6977E+01  2.6505E+01
 EBVSHRINKVR(%)  8.8276E-01  4.7478E+01  9.9996E+01  3.1072E+01  4.5985E+01
 RELATIVEINF(%)  9.8782E+01  5.7472E+00  4.7472E-04  9.6322E+00  4.3659E+00
 EPSSHRINKSD(%)  4.2856E+01
 EPSSHRINKVR(%)  6.7346E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1674.4029119747213     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -939.25208541098311     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1674.403       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  7.96E-01  8.24E-01  1.20E+00  7.55E-01  1.07E+00  1.64E+00  1.00E-02  8.74E-01  7.37E-01  1.08E+00
 


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
 #CPUT: Total CPU Time in Seconds,       40.491
Stop Time:
Sun Oct 24 03:34:49 CDT 2021
