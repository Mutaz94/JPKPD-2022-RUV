Sat Oct 23 16:19:53 CDT 2021
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
$DATA ../../../../data/SD1/D/dat38.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3454.88258479640        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5566E+02  7.9716E+01  1.2497E+02  3.9013E+01  6.6343E+01  7.8418E+00 -4.6770E+01 -4.4534E+02 -1.5173E+02 -3.6295E+00
            -9.3251E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3690.04543522334        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      151
 NPARAMETR:  9.7154E-01  1.0507E+00  9.1387E-01  1.0329E+00  1.0594E+00  1.1172E+00  1.1683E+00  2.5237E+00  1.2532E+00  9.3939E-01
             1.0002E+00
 PARAMETER:  7.1132E-02  1.4950E-01  9.9382E-03  1.3238E-01  1.5772E-01  2.1083E-01  2.5553E-01  1.0257E+00  3.2568E-01  3.7477E-02
             1.0016E-01
 GRADIENT:  -3.8339E+01  2.9387E+00 -3.7835E+01  5.5410E+01 -2.2208E+01  4.9532E+00 -1.9620E+01 -2.6330E+01  1.8108E+01  1.4254E+01
            -1.5574E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3706.49102612200        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      329
 NPARAMETR:  9.7996E-01  1.1049E+00  1.3018E+00  9.6112E-01  1.1857E+00  1.1395E+00  1.3005E+00  2.8850E+00  1.0836E+00  8.7769E-01
             1.0452E+00
 PARAMETER:  7.9754E-02  1.9972E-01  3.6377E-01  6.0344E-02  2.7037E-01  2.3057E-01  3.6277E-01  1.1595E+00  1.8029E-01 -3.0466E-02
             1.4425E-01
 GRADIENT:  -2.2974E+01 -6.3278E+00  6.6785E+00 -7.6431E+00 -9.1603E+00  1.3246E+01 -3.9867E+00 -2.3588E+01 -4.1505E+00  6.4283E+00
             5.5319E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3707.71034269712        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  9.8664E-01  1.0492E+00  1.3262E+00  9.9320E-01  1.1552E+00  1.1161E+00  1.3284E+00  2.9843E+00  1.0522E+00  7.9230E-01
             1.0329E+00
 PARAMETER:  8.6553E-02  1.4805E-01  3.8231E-01  9.3173E-02  2.4426E-01  2.0980E-01  3.8400E-01  1.1934E+00  1.5085E-01 -1.3282E-01
             1.3238E-01
 GRADIENT:  -1.1690E+01 -1.4265E+01  6.2061E+00 -2.1034E+00 -4.7286E+00  5.3917E+00 -6.3336E+00 -1.7514E+01 -7.7508E+00  2.3772E+00
             3.6837E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3708.03114144842        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      654
 NPARAMETR:  9.8670E-01  1.0498E+00  1.3078E+00  9.9154E-01  1.1525E+00  1.1096E+00  1.3280E+00  2.9822E+00  1.0606E+00  7.8628E-01
             1.0271E+00
 PARAMETER:  8.6611E-02  1.4858E-01  3.6831E-01  9.1508E-02  2.4197E-01  2.0402E-01  3.8365E-01  1.1927E+00  1.5884E-01 -1.4044E-01
             1.2675E-01
 GRADIENT:  -1.1651E+01 -1.3524E+01  4.3200E+00 -3.8977E+00 -5.0401E+00  3.0819E+00 -6.2574E+00 -1.6769E+01 -6.2710E+00  1.4787E+00
             2.5905E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3708.16292254625        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      833
 NPARAMETR:  9.8695E-01  1.0494E+00  1.3066E+00  9.9129E-01  1.1532E+00  1.1021E+00  1.3654E+00  2.9738E+00  1.0855E+00  7.7295E-01
             1.0274E+00
 PARAMETER:  8.6867E-02  1.4820E-01  3.6744E-01  9.1255E-02  2.4258E-01  1.9725E-01  4.1147E-01  1.1899E+00  1.8208E-01 -1.5754E-01
             1.2707E-01
 GRADIENT:  -1.1359E+01 -1.2352E+01  4.1893E+00 -2.8416E+00 -1.2382E+00  3.4164E-01 -1.7842E-02 -1.7327E+01  3.0293E-01  3.0460E-01
             2.6468E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3708.30282886760        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      998
 NPARAMETR:  9.8879E-01  1.0531E+00  1.3075E+00  9.9140E-01  1.1538E+00  1.1003E+00  1.3624E+00  2.9671E+00  1.0847E+00  7.7220E-01
             1.0213E+00
 PARAMETER:  8.8722E-02  1.5178E-01  3.6813E-01  9.1365E-02  2.4305E-01  1.9562E-01  4.0921E-01  1.1876E+00  1.8133E-01 -1.5852E-01
             1.2108E-01
 GRADIENT:  -7.9360E+00 -9.2807E+00  4.6018E+00 -1.3302E+00 -1.2316E+00 -2.9449E-01 -3.8739E-01 -1.8130E+01 -1.6802E-01 -3.2521E-01
             1.4011E+01

0ITERATION NO.:   31    OBJECTIVE VALUE:  -3708.30282886760        NO. OF FUNC. EVALS.:  27
 CUMULATIVE NO. OF FUNC. EVALS.:     1025
 NPARAMETR:  9.8893E-01  1.0529E+00  1.3069E+00  9.9126E-01  1.1542E+00  1.1017E+00  1.3634E+00  2.9621E+00  1.0850E+00  7.7342E-01
             1.0215E+00
 PARAMETER:  8.8722E-02  1.5178E-01  3.6813E-01  9.1365E-02  2.4305E-01  1.9562E-01  4.0921E-01  1.1876E+00  1.8133E-01 -1.5852E-01
             1.2108E-01
 GRADIENT:  -1.1846E+04  1.5585E+04  6.2170E+03  2.3667E+04 -1.0092E+04 -8.4022E-01 -2.6825E-01  2.0142E+03 -2.0245E-01 -2.9074E-01
            -1.9545E+04
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.3         1.7         2.2         2.3         2.3         1.5
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1025
 NO. OF SIG. DIGITS IN FINAL EST.:  1.5
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.1392E-03 -1.6667E-02 -3.1127E-02  1.6375E-02 -5.9678E-02
 SE:             2.9996E-02  2.5079E-02  2.7047E-02  2.6278E-02  1.9462E-02
 N:                     100         100         100         100         100

 P VAL.:         8.6396E-01  5.0632E-01  2.4979E-01  5.3317E-01  2.1671E-03

 ETASHRINKSD(%)  1.0000E-10  1.5982E+01  9.3902E+00  1.1967E+01  3.4799E+01
 ETASHRINKVR(%)  1.0000E-10  2.9410E+01  1.7899E+01  2.2501E+01  5.7489E+01
 EBVSHRINKSD(%)  2.3250E-01  1.6144E+01  1.3694E+01  1.3054E+01  3.5114E+01
 EBVSHRINKVR(%)  4.6447E-01  2.9682E+01  2.5513E+01  2.4405E+01  5.7898E+01
 RELATIVEINF(%)  9.9534E+01  4.5622E+01  6.4785E+01  5.5283E+01  2.6416E+01
 EPSSHRINKSD(%)  2.3360E+01
 EPSSHRINKVR(%)  4.1263E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3708.3028288675973     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2054.2134690991866     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3708.303       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.89E-01  1.05E+00  1.31E+00  9.91E-01  1.15E+00  1.10E+00  1.36E+00  2.97E+00  1.08E+00  7.72E-01  1.02E+00
 


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
 #CPUT: Total CPU Time in Seconds,       76.038
Stop Time:
Sat Oct 23 16:20:06 CDT 2021
