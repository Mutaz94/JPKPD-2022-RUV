Sun Oct 24 04:14:41 CDT 2021
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
$DATA ../../../../data/SD4/D/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1609.95607116174        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2218E+02 -5.2817E+01  2.8744E+01 -8.2862E+01 -1.6069E+01  1.5924E+01 -2.7782E+01 -5.8745E+00 -1.6502E+01 -2.7182E+01
            -4.3170E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1625.99884426794        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0837E+00  1.1892E+00  1.0004E+00  1.0951E+00  1.1077E+00  1.1594E+00  1.3271E+00  1.0098E+00  1.0555E+00  1.2381E+00
             1.1040E+00
 PARAMETER:  1.8041E-01  2.7326E-01  1.0044E-01  1.9086E-01  2.0228E-01  2.4790E-01  3.8301E-01  1.0975E-01  1.5399E-01  3.1354E-01
             1.9893E-01
 GRADIENT:   5.4759E+01  6.5464E+01 -1.3530E+01  1.1426E+02  1.9099E+00  1.4247E+01 -5.1376E-02  1.7106E+00  1.2475E+01  1.0968E+00
             4.6022E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1627.36892101153        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0933E+00  9.6100E-01  1.2030E+00  1.2331E+00  1.1200E+00  1.1267E+00  1.6704E+00  1.0755E+00  8.9198E-01  1.3242E+00
             1.1226E+00
 PARAMETER:  1.8922E-01  6.0218E-02  2.8479E-01  3.0955E-01  2.1337E-01  2.1933E-01  6.1307E-01  1.7276E-01 -1.4311E-02  3.8084E-01
             2.1566E-01
 GRADIENT:   7.8794E+01  5.4823E+01 -7.9594E+00  1.1494E+02  3.8708E+00  3.4581E+00  4.4790E+00 -1.2245E+00  4.5060E+00  3.5276E+00
             9.2078E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1634.39228010084        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  1.0430E+00  9.2576E-01  1.1718E+00  1.1583E+00  1.1037E+00  1.1076E+00  1.6325E+00  1.0948E+00  8.2649E-01  1.2723E+00
             1.0933E+00
 PARAMETER:  1.4212E-01  2.2856E-02  2.5856E-01  2.4695E-01  1.9869E-01  2.0222E-01  5.9014E-01  1.9061E-01 -9.0565E-02  3.4081E-01
             1.8923E-01
 GRADIENT:  -4.2909E+00  3.2016E+00  1.1538E-01  1.1532E+00  6.3403E-01  9.7912E-02  1.0774E+00 -4.1210E-01  6.3568E-01 -7.5476E-03
             4.6385E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1634.82319907169        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0419E+00  6.1995E-01  1.4305E+00  1.3566E+00  1.0872E+00  1.0954E+00  2.1589E+00  1.3178E+00  7.4877E-01  1.3014E+00
             1.0916E+00
 PARAMETER:  1.4104E-01 -3.7811E-01  4.5805E-01  4.0495E-01  1.8364E-01  1.9116E-01  8.6958E-01  3.7597E-01 -1.8933E-01  3.6345E-01
             1.8765E-01
 GRADIENT:   2.1794E+00  5.9915E+00  4.0755E+00  7.7106E+00 -7.2458E+00 -1.7832E+00  4.9019E-01  6.6528E-02 -1.3826E+00  1.9085E-01
            -1.3403E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1634.96387828442        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  1.0385E+00  4.9962E-01  1.5349E+00  1.4301E+00  1.0964E+00  1.0968E+00  2.4632E+00  1.4069E+00  7.3499E-01  1.3212E+00
             1.0958E+00
 PARAMETER:  1.3773E-01 -5.9391E-01  5.2849E-01  4.5771E-01  1.9205E-01  1.9243E-01  1.0014E+00  4.4139E-01 -2.0790E-01  3.7851E-01
             1.9147E-01
 GRADIENT:  -3.0509E-01  3.2711E+00  1.7817E+00  5.0375E+00 -1.7833E+00 -7.7325E-02  5.8544E-01 -8.5003E-02 -5.4805E-01 -2.4140E-01
            -1.7262E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1634.97218219899        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1075
 NPARAMETR:  1.0374E+00  4.4534E-01  1.5617E+00  1.4634E+00  1.0903E+00  1.0961E+00  2.6290E+00  1.4270E+00  7.2936E-01  1.3242E+00
             1.0963E+00
 PARAMETER:  1.3676E-01 -7.0891E-01  5.4576E-01  4.8079E-01  1.8642E-01  1.9174E-01  1.0666E+00  4.5557E-01 -2.1559E-01  3.8079E-01
             1.9193E-01
 GRADIENT:  -5.1650E-01  2.6317E+00  1.0611E+00  5.9108E+00 -9.5541E-01  1.4944E-01  3.0587E-01 -1.6871E-01 -5.2713E-01 -2.0817E-01
            -4.3976E-02

0ITERATION NO.:   34    OBJECTIVE VALUE:  -1635.00737173787        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1210
 NPARAMETR:  1.0381E+00  4.3997E-01  1.5581E+00  1.4612E+00  1.0902E+00  1.0961E+00  2.6343E+00  1.4275E+00  7.2999E-01  1.3257E+00
             1.0961E+00
 PARAMETER:  1.3742E-01 -7.2104E-01  5.4349E-01  4.7929E-01  1.8637E-01  1.9179E-01  1.0686E+00  4.5594E-01 -2.1473E-01  3.8194E-01
             1.9174E-01
 GRADIENT:  -9.2268E-01  3.3702E-01  2.7389E-01  4.3822E+00  4.3669E-01 -2.1629E-01 -4.1539E-01  1.5581E-01 -4.5363E-03 -8.4447E-02
             1.4782E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1210
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1972E-03  2.8088E-02 -4.5097E-02 -3.3446E-02 -4.3115E-02
 SE:             2.9936E-02  1.9085E-02  1.7011E-02  2.2736E-02  2.0123E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4149E-01  1.4109E-01  8.0222E-03  1.4127E-01  3.2150E-02

 ETASHRINKSD(%)  1.0000E-10  3.6064E+01  4.3013E+01  2.3831E+01  3.2584E+01
 ETASHRINKVR(%)  1.0000E-10  5.9122E+01  6.7524E+01  4.1983E+01  5.4551E+01
 EBVSHRINKSD(%)  4.6211E-01  4.0260E+01  4.7826E+01  2.0657E+01  2.7329E+01
 EBVSHRINKVR(%)  9.2208E-01  6.4311E+01  7.2778E+01  3.7047E+01  4.7190E+01
 RELATIVEINF(%)  9.8471E+01  5.5824E+00  8.4947E+00  1.0225E+01  1.5817E+01
 EPSSHRINKSD(%)  4.4886E+01
 EPSSHRINKVR(%)  6.9625E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1635.0073717378702     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -899.85654517413207     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.74
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1635.007       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  4.40E-01  1.56E+00  1.46E+00  1.09E+00  1.10E+00  2.63E+00  1.43E+00  7.30E-01  1.33E+00  1.10E+00
 


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
 #CPUT: Total CPU Time in Seconds,       36.863
Stop Time:
Sun Oct 24 04:14:49 CDT 2021
