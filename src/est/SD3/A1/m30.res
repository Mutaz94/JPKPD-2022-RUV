Sat Oct 23 21:44:48 CDT 2021
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
$DATA ../../../../data/SD3/A1/dat30.csv ignore=@
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
Current Date:       23 OCT 2021
Days until program expires : 176
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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1715.13017381967        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5323E+02  4.3738E+01  3.4872E+01  7.9314E+01  2.0258E+01  5.7577E+01 -1.2699E+01 -1.1442E+01  2.2561E+01 -6.3250E+00
            -7.9397E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1895.11113639888        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0387E+00  1.1120E+00  9.1083E-01  9.4634E-01  1.0214E+00  1.0498E+00  1.1913E+00  8.8621E-01  7.0496E-01  9.3128E-01
             1.7335E+00
 PARAMETER:  1.3801E-01  2.0615E-01  6.6054E-03  4.4843E-02  1.2118E-01  1.4856E-01  2.7504E-01 -2.0806E-02 -2.4961E-01  2.8810E-02
             6.5012E-01
 GRADIENT:   2.6608E+02  4.1281E+01  6.0620E+00  3.1675E+01  1.0016E+01  4.3000E+01 -1.7798E+00  7.2135E-01 -1.1235E+01 -5.3476E+00
             8.3827E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1897.84211896595        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0198E+00  1.0660E+00  6.3206E-01  9.4486E-01  8.2668E-01  1.0438E+00  1.2044E+00  4.4386E-01  7.2517E-01  7.4236E-01
             1.7091E+00
 PARAMETER:  1.1957E-01  1.6391E-01 -3.5877E-01  4.3286E-02 -9.0339E-02  1.4289E-01  2.8601E-01 -7.1226E-01 -2.2134E-01 -1.9793E-01
             6.3597E-01
 GRADIENT:   1.9659E+02  2.7728E+01  1.0154E+00  3.1017E+01  1.5758E+01  4.1123E+01  3.0791E-01  8.9239E-01 -1.2886E+01 -7.8628E+00
             4.6003E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1900.38957362658        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      315
 NPARAMETR:  1.0077E+00  9.9428E-01  6.9213E-01  9.9483E-01  8.3023E-01  1.0124E+00  1.2664E+00  2.0827E-01  8.1046E-01  8.5300E-01
             1.7113E+00
 PARAMETER:  1.0765E-01  9.4268E-02 -2.6798E-01  9.4815E-02 -8.6053E-02  1.1229E-01  3.3615E-01 -1.4689E+00 -1.1016E-01 -5.8990E-02
             6.3726E-01
 GRADIENT:  -9.7407E+00  4.5835E+00 -2.4118E+00  3.7271E+00  1.7030E+00  2.2143E+00  9.2083E-02  2.5636E-01  1.0938E+00  1.6883E+00
             5.4292E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1901.01883372618        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      490
 NPARAMETR:  1.0118E+00  7.9719E-01  7.3871E-01  1.1074E+00  7.7217E-01  1.0057E+00  1.5027E+00  7.7596E-02  7.5536E-01  8.3637E-01
             1.7033E+00
 PARAMETER:  1.1170E-01 -1.2666E-01 -2.0285E-01  2.0203E-01 -1.5856E-01  1.0571E-01  5.0729E-01 -2.4562E+00 -1.8056E-01 -7.8689E-02
             6.3260E-01
 GRADIENT:   3.7350E+00  1.6976E+00  8.1307E-01  1.9698E+00 -2.0762E+00  7.2691E-01 -1.7742E-01  2.2967E-02 -5.6169E-01  2.9397E-02
            -1.7922E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1901.03434773713        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      666            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0104E+00  7.8135E-01  7.4935E-01  1.1157E+00  7.7379E-01  1.0038E+00  1.5277E+00  2.5566E-02  7.5420E-01  8.4226E-01
             1.7035E+00
 PARAMETER:  1.1039E-01 -1.4673E-01 -1.8854E-01  2.0945E-01 -1.5646E-01  1.0379E-01  5.2379E-01 -3.5665E+00 -1.8210E-01 -7.1669E-02
             6.3269E-01
 GRADIENT:   1.7192E+02  1.2270E+01  3.3051E+00  8.6068E+01  6.5206E+00  2.0612E+01  9.0420E+00  4.9025E-03  4.2995E+00  8.2296E-01
             6.0435E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1901.03532292165        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      843
 NPARAMETR:  1.0101E+00  7.8116E-01  7.4932E-01  1.1165E+00  7.7392E-01  1.0037E+00  1.5263E+00  1.2984E-02  7.5477E-01  8.4367E-01
             1.7041E+00
 PARAMETER:  1.1008E-01 -1.4697E-01 -1.8858E-01  2.1016E-01 -1.5628E-01  1.0372E-01  5.2287E-01 -4.2440E+00 -1.8135E-01 -6.9990E-02
             6.3305E-01
 GRADIENT:   9.0426E-01  1.1951E-01 -3.0669E-02  6.5047E-01  1.4431E-02  1.2183E-01 -7.7212E-03  6.0719E-04  7.3474E-02  3.7416E-02
             9.2748E-02

0ITERATION NO.:   33    OBJECTIVE VALUE:  -1901.03547922388        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      935
 NPARAMETR:  1.0101E+00  7.8078E-01  7.4956E-01  1.1165E+00  7.7398E-01  1.0036E+00  1.5281E+00  1.0000E-02  7.5420E-01  8.4347E-01
             1.7040E+00
 PARAMETER:  1.1004E-01 -1.4746E-01 -1.8826E-01  2.1023E-01 -1.5620E-01  1.0364E-01  5.2400E-01 -4.6415E+00 -1.8210E-01 -7.0226E-02
             6.3297E-01
 GRADIENT:   8.3636E-01  7.0028E-02  2.1032E-02  3.7511E-01  3.8936E-02  9.6995E-02  5.3249E-02  0.0000E+00  1.6509E-03 -1.0913E-02
             2.6903E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      935
 NO. OF SIG. DIGITS IN FINAL EST.:  2.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.5799E-04  1.3683E-02 -2.9685E-04 -1.7305E-02 -6.9786E-03
 SE:             2.9680E-02  2.1293E-02  1.6843E-04  2.3018E-02  2.0666E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7694E-01  5.2049E-01  7.8000E-02  4.5218E-01  7.3560E-01

 ETASHRINKSD(%)  5.6782E-01  2.8665E+01  9.9436E+01  2.2886E+01  3.0766E+01
 ETASHRINKVR(%)  1.1324E+00  4.9114E+01  9.9997E+01  4.0535E+01  5.2066E+01
 EBVSHRINKSD(%)  9.1215E-01  2.9196E+01  9.9442E+01  2.2211E+01  2.9568E+01
 EBVSHRINKVR(%)  1.8160E+00  4.9868E+01  9.9997E+01  3.9489E+01  5.0394E+01
 RELATIVEINF(%)  9.7514E+01  4.7197E+00  2.7501E-04  6.1982E+00  4.3459E+00
 EPSSHRINKSD(%)  2.9417E+01
 EPSSHRINKVR(%)  5.0181E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1901.0354792238843     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -982.09694601921160     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1901.035       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  7.81E-01  7.50E-01  1.12E+00  7.74E-01  1.00E+00  1.53E+00  1.00E-02  7.54E-01  8.43E-01  1.70E+00
 


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
 #CPUT: Total CPU Time in Seconds,       76.144
Stop Time:
Sat Oct 23 21:45:01 CDT 2021
