Sun Oct 24 00:56:29 CDT 2021
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
$DATA ../../../../data/SD3/D/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2085.79981953812        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7691E+02 -5.9457E+01  2.5984E+00 -8.7544E+01  2.9243E+00  3.3923E+01 -2.0401E+01  9.4623E+00 -2.0768E+01  6.4488E+00
            -3.6448E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2097.91654523766        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0218E+00  1.0593E+00  1.0588E+00  1.0639E+00  1.0277E+00  1.0253E+00  1.1299E+00  9.4079E-01  1.0662E+00  9.6732E-01
             1.0623E+00
 PARAMETER:  1.2160E-01  1.5766E-01  1.5716E-01  1.6195E-01  1.2731E-01  1.2499E-01  2.2214E-01  3.8963E-02  1.6409E-01  6.6771E-02
             1.6045E-01
 GRADIENT:   1.7429E+01  1.1026E+00  1.1261E+01 -1.5085E+01 -2.4774E+01  3.3198E+00 -8.5056E+00  4.7796E+00  5.8638E+00 -5.5100E-01
             5.2853E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2100.36912659825        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  1.0115E+00  8.6222E-01  1.0777E+00  1.2113E+00  9.6994E-01  1.0376E+00  1.5817E+00  6.6121E-01  9.1430E-01  1.0669E+00
             1.0469E+00
 PARAMETER:  1.1142E-01 -4.8248E-02  1.7484E-01  2.9170E-01  6.9480E-02  1.3695E-01  5.5852E-01 -3.1368E-01  1.0404E-02  1.6477E-01
             1.4580E-01
 GRADIENT:   9.1109E-01  2.0590E+01  6.8234E-04  1.9567E+01 -9.5871E+00  8.5675E+00  5.9242E+00  2.3263E+00 -3.1907E+00  1.2499E+01
            -3.4154E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2102.70039832789        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  1.0112E+00  7.2110E-01  9.6602E-01  1.2699E+00  8.5374E-01  1.0050E+00  1.6600E+00  3.1889E-01  9.0179E-01  9.5539E-01
             1.0532E+00
 PARAMETER:  1.1117E-01 -2.2698E-01  6.5434E-02  3.3891E-01 -5.8130E-02  1.0494E-01  6.0679E-01 -1.0429E+00 -3.3703E-03  5.4361E-02
             1.5182E-01
 GRADIENT:   1.6931E+00  7.3495E+00  6.5148E+00  6.3579E+00 -1.1625E+01 -3.4978E+00 -8.0615E-01  2.7900E-01  1.0024E+00  9.6655E-01
             1.2778E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2102.94250454261        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      710
 NPARAMETR:  1.0098E+00  6.3947E-01  9.5707E-01  1.3069E+00  8.2902E-01  1.0133E+00  1.8098E+00  2.0741E-01  8.7812E-01  9.5683E-01
             1.0529E+00
 PARAMETER:  1.0979E-01 -3.4711E-01  5.6126E-02  3.6766E-01 -8.7509E-02  1.1319E-01  6.9324E-01 -1.4730E+00 -2.9967E-02  5.5872E-02
             1.5157E-01
 GRADIENT:   1.0040E+00 -7.4274E-01 -1.0060E+00 -2.1310E+00  1.5156E+00  2.2185E-01 -5.8615E-02  1.8169E-02  9.1515E-02  7.2861E-02
             2.3389E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2102.94441111142        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      892
 NPARAMETR:  1.0099E+00  6.4069E-01  9.5684E-01  1.3066E+00  8.2880E-01  1.0130E+00  1.8121E+00  2.0274E-01  8.7794E-01  9.5673E-01
             1.0528E+00
 PARAMETER:  1.0989E-01 -3.4521E-01  5.5881E-02  3.6742E-01 -8.7780E-02  1.1291E-01  6.9448E-01 -1.4958E+00 -3.0180E-02  5.5761E-02
             1.5146E-01
 GRADIENT:   1.1928E+00 -2.1250E-01 -1.6505E-01 -1.7278E+00  4.7834E-01  1.0530E-01  1.4743E-01  8.9623E-03  4.8878E-02  5.9346E-02
             6.7879E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2102.94542502490        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1078             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0101E+00  6.4225E-01  9.5638E-01  1.3054E+00  8.2847E-01  1.0130E+00  1.8123E+00  1.9234E-01  8.7781E-01  9.5568E-01
             1.0527E+00
 PARAMETER:  1.1006E-01 -3.4278E-01  5.5397E-02  3.6649E-01 -8.8169E-02  1.1290E-01  6.9459E-01 -1.5485E+00 -3.0323E-02  5.4670E-02
             1.5134E-01
 GRADIENT:   3.9488E+02  4.4933E+01  4.3560E+00  3.8051E+02  6.9000E+00  3.8656E+01  2.6312E+01  1.3700E-01  9.7902E+00  5.9173E-01
             1.2762E+00

0ITERATION NO.:   33    OBJECTIVE VALUE:  -2102.94575064786        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     1173
 NPARAMETR:  1.0099E+00  6.4221E-01  9.5603E-01  1.3058E+00  8.2861E-01  1.0129E+00  1.8110E+00  1.9534E-01  8.7784E-01  9.5626E-01
             1.0528E+00
 PARAMETER:  1.0988E-01 -3.4284E-01  5.5036E-02  3.6678E-01 -8.8007E-02  1.1284E-01  6.9388E-01 -1.5330E+00 -3.0295E-02  5.5271E-02
             1.5141E-01
 GRADIENT:  -4.2582E-01 -2.2565E-02  4.4558E-01  1.2892E+00 -1.1810E-01 -2.8679E-02 -1.2703E-02 -2.3014E-03 -5.3048E-02 -3.6014E-02
            -6.7039E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1173
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.6466E-04  1.5233E-02 -9.2894E-03 -1.3224E-02 -7.9216E-03
 SE:             2.9879E-02  1.9230E-02  4.1243E-03  2.5607E-02  2.3377E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7424E-01  4.2829E-01  2.4299E-02  6.0555E-01  7.3471E-01

 ETASHRINKSD(%)  1.0000E-10  3.5576E+01  8.6183E+01  1.4214E+01  2.1685E+01
 ETASHRINKVR(%)  1.0000E-10  5.8496E+01  9.8091E+01  2.6408E+01  3.8668E+01
 EBVSHRINKSD(%)  3.7488E-01  3.8066E+01  8.7159E+01  1.3353E+01  1.8903E+01
 EBVSHRINKVR(%)  7.4836E-01  6.1641E+01  9.8351E+01  2.4923E+01  3.4233E+01
 RELATIVEINF(%)  9.8541E+01  4.6094E+00  2.8377E-01  1.1292E+01  8.9948E+00
 EPSSHRINKSD(%)  3.2834E+01
 EPSSHRINKVR(%)  5.4887E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2102.9457506478589     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1184.0072174431862     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2102.946       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  6.42E-01  9.56E-01  1.31E+00  8.29E-01  1.01E+00  1.81E+00  1.95E-01  8.78E-01  9.56E-01  1.05E+00
 


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
 #CPUT: Total CPU Time in Seconds,       39.602
Stop Time:
Sun Oct 24 00:56:38 CDT 2021
