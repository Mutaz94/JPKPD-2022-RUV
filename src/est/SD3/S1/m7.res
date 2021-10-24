Sat Oct 23 22:53:15 CDT 2021
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
$DATA ../../../../data/SD3/S1/dat7.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1527.36638425031        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7835E+02 -8.0336E+00  1.5750E+02 -9.8074E+00  8.6745E+01  7.7336E+00  1.5937E+01 -9.2786E+02 -1.9578E+02  2.4429E+01
            -1.0186E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1947.94728575556        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      107
 NPARAMETR:  1.0758E+00  1.0803E+00  1.0131E+00  8.7122E-01  1.0107E+00  1.1436E+00  9.0057E-01  2.1368E+00  7.3719E-01  8.1048E-01
             7.6495E-01
 PARAMETER:  1.7309E-01  1.7728E-01  1.1305E-01 -3.7855E-02  1.1063E-01  2.3421E-01 -4.7274E-03  8.5930E-01 -2.0490E-01 -1.1013E-01
            -1.6794E-01
 GRADIENT:   1.5307E+02 -2.0651E+02  3.5817E+01 -2.0707E+02  4.3373E+01 -6.4890E+00 -2.6802E+01 -3.3257E+02 -4.8604E+01  3.8078E+00
            -1.3936E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2015.21903265390        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      283
 NPARAMETR:  1.0995E+00  9.8364E-01  2.3894E+00  1.0392E+00  1.2455E+00  8.5495E-01  1.4323E+00  3.4028E+00  2.3510E-01  1.0499E+00
             9.8923E-01
 PARAMETER:  1.9486E-01  8.3504E-02  9.7104E-01  1.3846E-01  3.1951E-01 -5.6714E-02  4.5931E-01  1.3246E+00 -1.3477E+00  1.4871E-01
             8.9174E-02
 GRADIENT:   3.2705E+02  4.0797E+00  2.2486E+01  3.6844E+01  1.1093E+01 -1.7768E+02 -1.8817E+01 -1.8897E+02 -8.5404E+00  3.0778E+00
             5.0350E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2039.30228013695        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      459
 NPARAMETR:  1.0600E+00  8.8119E-01  5.6400E+00  1.1439E+00  1.5592E+00  9.4013E-01  1.6246E+00  4.2118E+00  2.8427E-01  1.4580E+00
             9.4999E-01
 PARAMETER:  1.5830E-01 -2.6481E-02  1.8299E+00  2.3443E-01  5.4420E-01  3.8261E-02  5.8529E-01  1.5379E+00 -1.1578E+00  4.7704E-01
             4.8696E-02
 GRADIENT:   1.7973E+02  2.8869E+01  3.0425E+01  8.3443E+01  8.4971E+01 -9.7860E+01 -1.2978E+01 -1.8316E+02 -1.0053E+01  1.0455E+01
             2.1685E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2077.78588210987        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  1.0265E+00  7.4456E-01  7.5421E+00  1.2347E+00  1.4266E+00  1.0030E+00  1.8942E+00  5.1310E+00  3.1710E-01  1.2957E+00
             9.3240E-01
 PARAMETER:  1.2615E-01 -1.9496E-01  2.1205E+00  3.1079E-01  4.5531E-01  1.0299E-01  7.3879E-01  1.7353E+00 -1.0485E+00  3.5906E-01
             3.0004E-02
 GRADIENT:   8.6372E+01  2.3911E+01  2.8860E+01  7.4449E+01  5.4861E+01 -5.4617E+01 -9.3115E+00 -1.4940E+02 -1.0040E+01  2.3158E+01
             2.7695E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2091.29585819648        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      822
 NPARAMETR:  1.0161E+00  6.7600E-01  8.8926E+00  1.2613E+00  1.3947E+00  1.1167E+00  1.9956E+00  5.6516E+00  3.8043E-01  1.1558E+00
             9.0307E-01
 PARAMETER:  1.1596E-01 -2.9156E-01  2.2852E+00  3.3210E-01  4.3267E-01  2.1041E-01  7.9096E-01  1.8319E+00 -8.6645E-01  2.4480E-01
            -1.9599E-03
 GRADIENT:   5.1659E+01 -2.9971E+00  1.9904E+01 -1.9234E+01  2.5623E+01 -3.8220E+00 -5.2393E+00 -1.1100E+02 -4.0829E-01 -2.6480E+00
            -7.4167E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2093.61853201294        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1000
 NPARAMETR:  1.0157E+00  6.2013E-01  8.7871E+00  1.3153E+00  1.3929E+00  1.1253E+00  2.1389E+00  5.8382E+00  4.0654E-01  1.3147E+00
             9.1484E-01
 PARAMETER:  1.1562E-01 -3.7783E-01  2.2733E+00  3.7403E-01  4.3139E-01  2.1809E-01  8.6028E-01  1.8644E+00 -8.0008E-01  3.7363E-01
             1.0999E-02
 GRADIENT:   4.9763E+01  4.5714E+00  1.2009E+01  1.3923E+01 -1.6315E+00 -9.6243E-01 -8.7574E+00 -9.2637E+01  2.3722E-01  1.5892E+01
             1.0355E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2100.53561447480        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:     1120
 NPARAMETR:  9.9724E-01  5.7392E-01  7.0420E+00  1.3337E+00  1.3912E+00  1.1215E+00  2.2482E+00  6.6304E+00  3.9851E-01  1.3648E+00
             9.0617E-01
 PARAMETER:  9.7241E-02 -4.5527E-01  2.0519E+00  3.8799E-01  4.3014E-01  2.1470E-01  9.1012E-01  1.9917E+00 -8.2003E-01  4.1104E-01
             1.4667E-03
 GRADIENT:   6.0858E+02  2.6005E+02  1.4366E+01  7.1232E+02  1.7429E+01  1.9083E+02  5.0054E+02  1.0747E+02  2.1389E+01  3.3458E+01
             8.8866E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2101.14436171202        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1289             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9133E-01  5.7212E-01  7.0651E+00  1.3325E+00  1.3990E+00  1.1232E+00  2.2455E+00  6.6096E+00  3.9986E-01  1.3368E+00
             9.0082E-01
 PARAMETER:  9.1293E-02 -4.5841E-01  2.0552E+00  3.8707E-01  4.3579E-01  2.1618E-01  9.0894E-01  1.9885E+00 -8.1664E-01  3.9025E-01
            -4.4519E-03
 GRADIENT:   6.0502E+02  2.6313E+02  1.4836E+01  7.1199E+02  2.7673E+01  1.9855E+02  5.0854E+02  1.0687E+02  2.1918E+01  2.9301E+01
             4.4697E+00

0ITERATION NO.:   44    OBJECTIVE VALUE:  -2102.21962188456        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1389
 NPARAMETR:  9.8960E-01  5.7147E-01  7.0747E+00  1.3304E+00  1.3976E+00  1.1233E+00  2.2346E+00  6.5907E+00  4.0055E-01  1.2707E+00
             8.9888E-01
 PARAMETER:  8.9570E-02 -4.5920E-01  2.0566E+00  3.8575E-01  4.3446E-01  2.1532E-01  9.0432E-01  1.9871E+00 -8.1633E-01  3.3928E-01
            -5.6004E-03
 GRADIENT:   1.5198E+04  2.6615E+03  1.3371E+02  1.5985E+03 -2.8073E+03 -1.1846E+00  1.6736E+03  1.0693E+03 -2.8335E-01 -3.5917E+03
             2.5123E+00
 NUMSIGDIG:         2.7         2.3         3.5         2.3         2.3         1.6         2.7         2.3         1.9         2.3
                    1.2

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1389
 NO. OF SIG. DIGITS IN FINAL EST.:  1.2
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.8838E-03  1.8881E-02 -7.4987E-02 -5.1607E-02 -9.2334E-02
 SE:             3.0040E-02  2.4144E-02  2.2626E-02  1.6379E-02  1.7697E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5000E-01  4.3420E-01  9.1915E-04  1.6280E-03  1.8167E-07

 ETASHRINKSD(%)  1.0000E-10  1.9114E+01  2.4200E+01  4.5130E+01  4.0713E+01
 ETASHRINKVR(%)  1.0000E-10  3.4575E+01  4.2544E+01  6.9893E+01  6.4851E+01
 EBVSHRINKSD(%)  2.4094E-01  1.7667E+01  4.0473E+01  4.6775E+01  2.2486E+01
 EBVSHRINKVR(%)  4.8131E-01  3.2213E+01  6.4565E+01  7.1671E+01  3.9916E+01
 RELATIVEINF(%)  9.9449E+01  1.3441E+01  2.5733E+01  5.5057E+00  4.5392E+01
 EPSSHRINKSD(%)  3.6258E+01
 EPSSHRINKVR(%)  5.9370E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2102.2196218845556     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1183.2810886798829     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2102.220       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.90E-01  5.72E-01  7.08E+00  1.33E+00  1.40E+00  1.12E+00  2.24E+00  6.60E+00  4.00E-01  1.27E+00  9.00E-01
 


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
 #CPUT: Total CPU Time in Seconds,      117.608
Stop Time:
Sat Oct 23 22:53:33 CDT 2021
