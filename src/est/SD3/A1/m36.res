Sat Oct 23 21:45:54 CDT 2021
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
$DATA ../../../../data/SD3/A1/dat36.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m36.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1611.93474510653        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0106E+02 -4.0609E+01  6.6222E+01 -4.3820E+00  1.0603E+02  2.3308E+01 -1.4921E+01 -1.9674E+02 -5.3734E+00 -2.2280E+01
            -7.2733E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1837.56986052189        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.8711E-01  1.1901E+00  7.9906E-01  8.8292E-01  9.2763E-01  1.0740E+00  1.0460E+00  1.6045E+00  7.2676E-01  8.0714E-01
             1.5755E+00
 PARAMETER:  8.7023E-02  2.7402E-01 -1.2432E-01 -2.4515E-02  2.4883E-02  1.7142E-01  1.4498E-01  5.7284E-01 -2.1915E-01 -1.1426E-01
             5.5459E-01
 GRADIENT:   1.1424E+02 -1.6057E+01  2.8903E+01 -6.4155E+01  1.5073E+01  3.8565E+01 -1.6155E+01 -5.0060E+01 -2.0133E+01  1.3333E+00
            -5.3027E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1841.90749075510        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.8347E-01  1.0717E+00  4.9235E-01  9.4129E-01  6.9363E-01  1.0623E+00  1.2665E+00  1.4367E+00  6.2606E-01  4.4372E-01
             1.5589E+00
 PARAMETER:  8.3331E-02  1.6923E-01 -6.0856E-01  3.9495E-02 -2.6582E-01  1.6042E-01  3.3629E-01  4.6232E-01 -3.6832E-01 -7.1256E-01
             5.4398E-01
 GRADIENT:   9.6261E+01  2.7158E+00  4.3417E+00  1.9678E+00  4.0308E+01  2.8722E+01  8.4727E+00 -3.3137E+01 -2.8567E+01 -1.7078E-01
            -5.7346E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1853.98294455792        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      312
 NPARAMETR:  1.0276E+00  1.1744E+00  5.6858E-01  9.0421E-01  7.8120E-01  1.0308E+00  1.1392E+00  1.7725E+00  7.8264E-01  3.2302E-01
             1.7471E+00
 PARAMETER:  1.2725E-01  2.6073E-01 -4.6461E-01 -6.9287E-04 -1.4692E-01  1.3033E-01  2.3037E-01  6.7238E-01 -1.4509E-01 -1.0300E+00
             6.5795E-01
 GRADIENT:   1.7564E+01 -3.0763E+01 -1.9648E+00 -1.0276E+01  2.9744E+01 -1.2464E+01  1.1013E+00 -2.0350E+01 -7.6391E+00 -6.6493E-01
             2.7929E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1860.47052680257        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      497
 NPARAMETR:  1.0064E+00  1.3952E+00  5.0989E-01  7.7519E-01  8.5052E-01  1.0845E+00  9.5394E-01  2.1182E+00  1.0180E+00  3.9599E-01
             1.6656E+00
 PARAMETER:  1.0642E-01  4.3305E-01 -5.7356E-01 -1.5465E-01 -6.1909E-02  1.8112E-01  5.2842E-02  8.5057E-01  1.1786E-01 -8.2637E-01
             6.1017E-01
 GRADIENT:  -2.2859E+01 -3.0524E+01 -2.4628E-01 -2.1434E+00  1.2734E+01  6.2769E+00  1.3483E+00 -1.3761E+01  2.8092E+00 -2.7987E+00
             9.8746E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1861.17589331320        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      683            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0098E+00  1.4004E+00  5.1048E-01  7.7571E-01  8.4860E-01  1.0666E+00  9.5404E-01  2.1117E+00  9.9135E-01  5.1077E-01
             1.6693E+00
 PARAMETER:  1.0979E-01  4.3677E-01 -5.7241E-01 -1.5397E-01 -6.4172E-02  1.6447E-01  5.2949E-02  8.4751E-01  9.1316E-02 -5.7184E-01
             6.1242E-01
 GRADIENT:   1.4280E+02  8.2469E+01  7.3601E+00  2.3156E+01 -3.2572E-01  2.8337E+01  4.4265E+00  1.9176E+00  1.8468E+00 -1.0072E-01
             2.3970E+01

0ITERATION NO.:   29    OBJECTIVE VALUE:  -1861.21981486406        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      825
 NPARAMETR:  1.0101E+00  1.4017E+00  5.0984E-01  7.7559E-01  8.4896E-01  1.0658E+00  9.5115E-01  2.1115E+00  9.9045E-01  5.1763E-01
             1.6695E+00
 PARAMETER:  1.1006E-01  4.3777E-01 -5.7376E-01 -1.5410E-01 -6.3728E-02  1.6371E-01  5.0919E-02  8.4755E-01  9.0628E-02 -5.5839E-01
             6.1239E-01
 GRADIENT:  -5.4394E+03  2.7007E+03 -2.0701E+03  7.7410E+03  1.1932E+04 -1.9416E-01  2.3236E+00  1.3651E+03  3.5930E-01  2.1373E+03
            -1.9199E+03
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.3         2.5         0.6         2.3         1.2         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      825
 NO. OF SIG. DIGITS IN FINAL EST.:  0.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.2995E-03 -3.3986E-03 -2.4093E-02  1.0582E-02 -2.6001E-02
 SE:             2.9729E-02  2.3557E-02  2.0950E-02  2.1746E-02  1.2620E-02
 N:                     100         100         100         100         100

 P VAL.:         7.5443E-01  8.8529E-01  2.5015E-01  6.2653E-01  3.9366E-02

 ETASHRINKSD(%)  4.0295E-01  2.1080E+01  2.9813E+01  2.7148E+01  5.7722E+01
 ETASHRINKVR(%)  8.0427E-01  3.7717E+01  5.0738E+01  4.6926E+01  8.2126E+01
 EBVSHRINKSD(%)  7.8608E-01  2.0549E+01  3.3790E+01  2.7883E+01  5.8424E+01
 EBVSHRINKVR(%)  1.5660E+00  3.6876E+01  5.6163E+01  4.7992E+01  8.2714E+01
 RELATIVEINF(%)  9.8386E+01  7.8711E+00  1.4251E+01  7.2112E+00  3.2999E+00
 EPSSHRINKSD(%)  3.3795E+01
 EPSSHRINKVR(%)  5.6169E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1861.2198148640584     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -942.28128165938574     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1861.220       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.40E+00  5.10E-01  7.76E-01  8.49E-01  1.07E+00  9.52E-01  2.11E+00  9.91E-01  5.18E-01  1.67E+00
 


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
 #CPUT: Total CPU Time in Seconds,       65.724
Stop Time:
Sat Oct 23 21:46:05 CDT 2021
