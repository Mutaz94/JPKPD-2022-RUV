Sun Oct 24 03:45:26 CDT 2021
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
$DATA ../../../../data/SD4/TD1/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1700.39352057083        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4545E+02  4.8316E+00 -5.1923E+01  8.6504E+01  1.0041E+02  6.6086E+01  1.1175E+00  1.0681E+01  4.3521E-01 -1.1512E+01
             8.2686E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1708.79381114348        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0529E+00  9.9613E-01  1.1072E+00  1.0300E+00  9.6764E-01  9.4020E-01  9.9011E-01  8.8636E-01  1.0562E+00  1.0523E+00
             9.7504E-01
 PARAMETER:  1.5154E-01  9.6119E-02  2.0179E-01  1.2953E-01  6.7106E-02  3.8337E-02  9.0061E-02 -2.0629E-02  1.5470E-01  1.5095E-01
             7.4718E-02
 GRADIENT:  -7.5101E+00  2.3704E+01  3.3969E+00  2.5361E+01 -1.3825E+01 -1.0229E+01  2.7904E+00  5.2430E+00  5.1623E+00 -4.2236E+00
            -4.7633E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1710.83343366410        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0572E+00  9.0819E-01  1.0723E+00  1.0641E+00  9.3519E-01  9.5826E-01  7.7158E-01  4.3899E-01  9.9537E-01  1.1455E+00
             9.7060E-01
 PARAMETER:  1.5564E-01  3.7014E-03  1.6977E-01  1.6209E-01  3.2995E-02  5.7369E-02 -1.5931E-01 -7.2327E-01  9.5364E-02  2.3583E-01
             7.0159E-02
 GRADIENT:   3.7852E+00  1.1536E+01 -1.7193E+00  3.9418E+00 -7.2543E+00 -2.1066E+00 -4.4938E+00  9.2175E-01 -1.4754E+01  1.9778E+00
            -7.3450E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1712.77445423266        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0542E+00  6.9342E-01  9.8625E-01  1.2036E+00  8.1531E-01  9.6403E-01  1.5295E+00  2.6718E-01  8.6858E-01  9.9131E-01
             9.8355E-01
 PARAMETER:  1.5275E-01 -2.6612E-01  8.6151E-02  2.8536E-01 -1.0418E-01  6.3367E-02  5.2493E-01 -1.2198E+00 -4.0894E-02  9.1274E-02
             8.3409E-02
 GRADIENT:  -1.6178E+00  1.4852E+01  4.4733E+00  2.6730E+01 -7.5323E+00  3.9002E-01  9.9494E-01  4.6503E-01  5.4282E+00 -1.6474E+00
             1.2655E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1714.36450317005        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0530E+00  4.4582E-01  9.4864E-01  1.3299E+00  7.2247E-01  9.5839E-01  2.0959E+00  5.9533E-02  7.4000E-01  9.7067E-01
             9.7500E-01
 PARAMETER:  1.5169E-01 -7.0783E-01  4.7269E-02  3.8512E-01 -2.2508E-01  5.7495E-02  8.3996E-01 -2.7212E+00 -2.0111E-01  7.0234E-02
             7.4685E-02
 GRADIENT:   2.1025E+00  7.3599E+00  7.6804E+00  1.2167E+01 -1.5372E+01 -5.7215E-01 -2.3732E+00  1.6817E-02 -8.4178E+00 -2.7684E+00
            -2.6898E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1715.36308200461        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      904
 NPARAMETR:  1.0468E+00  2.3768E-01  1.1222E+00  1.4700E+00  7.5415E-01  9.5562E-01  3.0291E+00  1.0000E-02  7.1711E-01  1.1096E+00
             9.8415E-01
 PARAMETER:  1.4573E-01 -1.3368E+00  2.1530E-01  4.8529E-01 -1.8216E-01  5.4606E-02  1.2083E+00 -5.1544E+00 -2.3253E-01  2.0402E-01
             8.4020E-02
 GRADIENT:  -6.7841E-01  5.5901E+00  5.5034E+00  1.9585E+01 -7.2612E+00  1.8015E-01  8.8015E-01  0.0000E+00 -2.8564E+00  6.9791E-01
             1.4058E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1716.21629268185        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1082
 NPARAMETR:  1.0437E+00  8.3630E-02  1.1389E+00  1.5547E+00  7.2882E-01  9.5336E-01  4.7882E+00  1.0000E-02  7.0632E-01  1.1414E+00
             9.8192E-01
 PARAMETER:  1.4277E-01 -2.3813E+00  2.3002E-01  5.4126E-01 -2.1632E-01  5.2236E-02  1.6662E+00 -1.0376E+01 -2.4769E-01  2.3228E-01
             8.1756E-02
 GRADIENT:  -4.8177E-01  1.6622E+00 -2.6890E+00  9.4531E+00 -5.9945E-01  2.5912E-01  3.0130E+00  0.0000E+00  2.1032E+00  3.6719E-02
             6.3652E-01

0ITERATION NO.:   32    OBJECTIVE VALUE:  -1716.23343482618        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1147
 NPARAMETR:  1.0434E+00  8.0915E-02  1.1476E+00  1.5572E+00  7.3121E-01  9.5293E-01  4.8278E+00  1.0000E-02  7.0630E-01  1.1448E+00
             9.8112E-01
 PARAMETER:  1.4267E-01 -2.4179E+00  2.3699E-01  5.4208E-01 -2.1274E-01  5.1796E-02  1.6693E+00 -1.0564E+01 -2.4735E-01  2.3758E-01
             8.1284E-02
 GRADIENT:   1.3740E+04 -8.1012E+02 -8.4982E-01 -3.5892E+03  9.1987E+03  1.3407E-02 -1.7210E+00  0.0000E+00  3.9661E+03  1.1711E+00
             2.5461E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1147
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6729E-03  1.9458E-02 -2.7058E-04 -2.4934E-02 -1.1394E-02
 SE:             2.9845E-02  1.2035E-02  1.7353E-04  2.7576E-02  2.5214E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5530E-01  1.0593E-01  1.1892E-01  3.6589E-01  6.5135E-01

 ETASHRINKSD(%)  1.4602E-02  5.9681E+01  9.9419E+01  7.6184E+00  1.5529E+01
 ETASHRINKVR(%)  2.9202E-02  8.3744E+01  9.9997E+01  1.4656E+01  2.8646E+01
 EBVSHRINKSD(%)  4.2687E-01  7.1918E+01  9.9449E+01  5.2251E+00  9.3676E+00
 EBVSHRINKVR(%)  8.5192E-01  9.2114E+01  9.9997E+01  1.0177E+01  1.7858E+01
 RELATIVEINF(%)  9.8850E+01  2.7777E+00  3.7904E-04  3.1279E+01  1.0500E+01
 EPSSHRINKSD(%)  4.3253E+01
 EPSSHRINKVR(%)  6.7798E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1716.2334348261782     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -981.08260826243998     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1716.233       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  8.06E-02  1.15E+00  1.56E+00  7.31E-01  9.53E-01  4.80E+00  1.00E-02  7.07E-01  1.15E+00  9.81E-01
 


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
 #CPUT: Total CPU Time in Seconds,       34.303
Stop Time:
Sun Oct 24 03:45:34 CDT 2021
