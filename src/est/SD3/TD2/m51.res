Sun Oct 24 00:45:43 CDT 2021
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
$DATA ../../../../data/SD3/TD2/dat51.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m51.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1470.71052548614        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1717E+02  9.0573E+00  1.7732E+02 -1.6139E+00  9.6769E+01  5.3336E+01  1.8088E+00 -8.8681E+02 -2.0693E+02 -1.7585E+01
            -8.4205E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1859.80430091819        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.2315E+00  1.1492E+00  8.7160E-01  9.4421E-01  9.2096E-01  1.5411E+00  9.8606E-01  3.8101E+00  1.1243E+00  9.5842E-01
             1.0187E+00
 PARAMETER:  3.0821E-01  2.3907E-01 -3.7427E-02  4.2589E-02  1.7666E-02  5.3251E-01  8.5967E-02  1.4377E+00  2.1715E-01  5.7535E-02
             1.1857E-01
 GRADIENT:   2.5077E+02 -4.3653E+01  1.4799E+00  5.2429E+01 -8.2262E+01  6.1224E+01  1.5270E-01 -7.8198E+01  9.2059E+00  2.7373E+01
            -1.3127E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1971.53794061607        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      278
 NPARAMETR:  9.9984E-01  1.1488E+00  9.4474E-01  8.8155E-01  9.7588E-01  9.2197E-01  8.8103E-01  3.8021E+00  1.0674E+00  7.8602E-01
             1.0083E+00
 PARAMETER:  9.9837E-02  2.3874E-01  4.3149E-02 -2.6078E-02  7.5586E-02  1.8755E-02 -2.6668E-02  1.4355E+00  1.6525E-01 -1.4078E-01
             1.0831E-01
 GRADIENT:   5.3156E+02  8.8068E+01 -7.2227E+00  2.9933E+01 -2.0261E+00  2.6665E+01  1.0022E+01  7.7336E+01  2.2694E+00  2.0762E+00
             5.1114E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1972.20797489373        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:      427
 NPARAMETR:  9.9606E-01  1.1483E+00  9.4975E-01  8.8021E-01  9.7537E-01  9.5436E-01  8.7426E-01  3.7915E+00  1.0670E+00  7.8613E-01
             1.0080E+00
 PARAMETER:  9.6049E-02  2.3827E-01  4.8441E-02 -2.7591E-02  7.5059E-02  5.3289E-02 -3.4372E-02  1.4328E+00  1.6481E-01 -1.4063E-01
             1.0802E-01
 GRADIENT:   9.3265E+01 -8.9271E+01 -8.1806E+00 -2.6683E+01 -1.0929E+01 -5.1418E+00  1.2913E+00 -4.9465E+01 -8.1608E+00  1.3607E+00
             4.1742E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1972.61665796926        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      614
 NPARAMETR:  9.9626E-01  1.1481E+00  9.5005E-01  8.7997E-01  9.7509E-01  9.8228E-01  8.7453E-01  3.8178E+00  1.1222E+00  7.8580E-01
             1.0084E+00
 PARAMETER:  9.6257E-02  2.3811E-01  4.8760E-02 -2.7864E-02  7.4771E-02  8.2125E-02 -3.4073E-02  1.4397E+00  2.1529E-01 -1.4106E-01
             1.0833E-01
 GRADIENT:   8.8756E+01 -9.0438E+01 -8.3270E+00 -2.3437E+01 -1.2315E+01  6.0966E+00  4.5402E+00 -4.7350E+01  1.5697E+00  1.9195E+00
             5.4882E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1976.55401622281        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      793
 NPARAMETR:  9.9543E-01  1.1592E+00  9.5161E-01  8.7940E-01  9.7416E-01  9.6661E-01  8.7553E-01  4.8105E+00  1.1183E+00  7.8429E-01
             1.0096E+00
 PARAMETER:  9.5423E-02  2.4770E-01  5.0398E-02 -2.8517E-02  7.3818E-02  6.6035E-02 -3.2927E-02  1.6708E+00  2.1179E-01 -1.4298E-01
             1.0953E-01
 GRADIENT:   8.9063E+01 -8.1369E+01 -1.4549E+01 -1.6769E+01 -5.7242E+01 -7.6574E-02  6.7727E+00  3.9114E+00  4.3015E-01  8.5042E+00
             5.7335E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1979.32435951516        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      934
 NPARAMETR:  9.6973E-01  1.1593E+00  9.9139E-01  8.7880E-01  9.8090E-01  9.5412E-01  8.5694E-01  4.4838E+00  1.1109E+00  7.7435E-01
             1.0077E+00
 PARAMETER:  6.9263E-02  2.4784E-01  9.1351E-02 -2.9198E-02  8.0716E-02  5.3035E-02 -5.4392E-02  1.6005E+00  2.0516E-01 -1.5573E-01
             1.0768E-01
 GRADIENT:   2.8065E+01 -8.3688E+01 -1.0554E+01 -1.8886E+01 -3.8294E+01 -2.1874E+00  3.6079E+00 -1.2644E+01 -1.5927E+00  3.7264E+00
             4.9377E+00

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1979.32435951516        NO. OF FUNC. EVALS.:  31
 CUMULATIVE NO. OF FUNC. EVALS.:      965
 NPARAMETR:  9.6977E-01  1.1592E+00  9.9135E-01  8.7876E-01  9.8094E-01  9.5504E-01  8.5697E-01  4.4809E+00  1.1110E+00  7.7440E-01
             1.0066E+00
 PARAMETER:  6.9263E-02  2.4784E-01  9.1351E-02 -2.9198E-02  8.0716E-02  5.3035E-02 -5.4392E-02  1.6005E+00  2.0516E-01 -1.5573E-01
             1.0768E-01
 GRADIENT:  -8.2653E+03  3.2538E+03  8.2842E+03  8.2673E+03 -4.1857E+03 -2.2806E+00 -8.2872E+03  4.9823E+02 -2.0214E+03 -5.3219E+03
             4.7754E+00
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.3         0.9         2.3         2.3         2.3         2.3
                    0.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      965
 NO. OF SIG. DIGITS IN FINAL EST.:  0.9
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.0794E-02  4.0739E-03 -2.7311E-02  2.0342E-02 -6.3474E-02
 SE:             3.0045E-02  1.7070E-02  2.3558E-02  2.5181E-02  1.7118E-02
 N:                     100         100         100         100         100

 P VAL.:         7.1939E-01  8.1137E-01  2.4634E-01  4.1919E-01  2.0886E-04

 ETASHRINKSD(%)  1.0000E-10  4.2813E+01  2.1077E+01  1.5639E+01  4.2654E+01
 ETASHRINKVR(%)  1.0000E-10  6.7296E+01  3.7711E+01  2.8832E+01  6.7114E+01
 EBVSHRINKSD(%)  3.7492E-01  4.0914E+01  1.8267E+01  1.7603E+01  4.4912E+01
 EBVSHRINKVR(%)  7.4843E-01  6.5088E+01  3.3197E+01  3.2107E+01  6.9653E+01
 RELATIVEINF(%)  9.9172E+01  5.6218E+00  3.6626E+01  1.2148E+01  1.4617E+01
 EPSSHRINKSD(%)  3.7277E+01
 EPSSHRINKVR(%)  6.0659E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1979.3243595151605     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1060.3858263104878     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.51
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1979.324       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.70E-01  1.16E+00  9.91E-01  8.79E-01  9.81E-01  9.54E-01  8.57E-01  4.48E+00  1.11E+00  7.74E-01  1.01E+00
 


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
 #CPUT: Total CPU Time in Seconds,       36.303
Stop Time:
Sun Oct 24 00:45:51 CDT 2021
