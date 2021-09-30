Thu Sep 30 06:50:59 CDT 2021
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
$DATA ../../../../data/spa2/A3/dat94.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1306.63918443466        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2785E+02  1.7462E+02  1.7882E+02  1.4691E+01  3.3868E+02  5.7701E+01 -1.7917E+02 -1.4900E+02 -1.4094E+02 -1.6055E+02
            -6.8981E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1598.83891488411        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0192E+00  1.0108E+00  9.9933E-01  1.2029E+00  8.8981E-01  8.0866E-01  1.0009E+00  9.6570E-01  1.1006E+00  8.3923E-01
             5.2775E+00
 PARAMETER:  1.1905E-01  1.1078E-01  9.9330E-02  2.8474E-01 -1.6745E-02 -1.1238E-01  1.0092E-01  6.5094E-02  1.9588E-01 -7.5275E-02
             1.7634E+00
 GRADIENT:  -1.0135E+02  3.2073E+01 -1.5407E+01  7.1151E+01 -7.5314E+00 -2.1235E+01  1.1367E+01  8.8213E+00  1.0969E+01  3.1597E+01
             3.5598E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1617.97153759592        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.0125E+00  6.8104E-01  4.4858E-01  1.3450E+00  4.7113E-01  8.6772E-01  1.0041E+00  3.6478E-01  1.0085E+00  4.1753E-01
             4.9668E+00
 PARAMETER:  1.1243E-01 -2.8414E-01 -7.0166E-01  3.9637E-01 -6.5263E-01 -4.1884E-02  1.0405E-01 -9.0847E-01  1.0845E-01 -7.7339E-01
             1.7028E+00
 GRADIENT:  -1.0534E+02  8.6171E+01 -1.0836E+01  1.5688E+02 -4.7558E+01 -1.0284E+01  3.8477E+00  2.3164E+00 -1.5144E+01  7.9274E+00
             3.0133E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1682.33636062188        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0212E+00  3.9451E-01  3.1296E-01  1.3224E+00  3.1388E-01  8.8130E-01  1.0485E+00  1.0000E-02  1.2606E+00  4.1615E-01
             3.4657E+00
 PARAMETER:  1.2102E-01 -8.3012E-01 -1.0617E+00  3.7943E-01 -1.0588E+00 -2.6353E-02  1.4733E-01 -4.5721E+00  3.3157E-01 -7.7671E-01
             1.3429E+00
 GRADIENT:   2.6400E+01  5.1156E+01  1.2505E+01  1.7256E+02  5.3291E+01 -8.9865E+00 -5.0720E+00  0.0000E+00 -3.4145E+01 -1.6369E+01
            -1.6148E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1703.06925587164        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      406
 NPARAMETR:  1.0183E+00  3.2786E-01  2.8980E-01  1.1768E+00  2.8612E-01  9.1754E-01  5.1417E-01  1.0000E-02  1.4960E+00  6.4617E-01
             3.3116E+00
 PARAMETER:  1.1811E-01 -1.0152E+00 -1.1385E+00  2.6278E-01 -1.1514E+00  1.3940E-02 -5.6521E-01 -6.2769E+00  5.0283E-01 -3.3669E-01
             1.2974E+00
 GRADIENT:   4.2752E-01 -5.4287E+00  3.0403E-01  1.8898E+01  3.9262E+01  4.8470E+00  6.0238E-01  0.0000E+00  5.0529E+00 -1.6114E+00
            -3.7967E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1707.57292134381        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      581
 NPARAMETR:  1.0214E+00  2.9137E-01  2.2615E-01  1.0864E+00  2.4071E-01  9.0659E-01  2.6793E-01  1.0000E-02  1.5307E+00  6.6088E-01
             3.3985E+00
 PARAMETER:  1.2115E-01 -1.1332E+00 -1.3866E+00  1.8291E-01 -1.3241E+00  1.9318E-03 -1.2170E+00 -6.3041E+00  5.2572E-01 -3.1418E-01
             1.3233E+00
 GRADIENT:   4.7454E-01 -1.7006E-01  4.1602E-01 -1.5542E+00 -7.7857E-01 -3.4064E-01  3.3757E-01  0.0000E+00  1.0290E-01  4.8162E-01
             1.6495E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1707.68946782312        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      675
 NPARAMETR:  1.0204E+00  2.9050E-01  2.2470E-01  1.0899E+00  2.3969E-01  9.0945E-01  1.8903E-02  1.0000E-02  1.5314E+00  6.5912E-01
             3.3846E+00
 PARAMETER:  1.2018E-01 -1.1362E+00 -1.3930E+00  1.8611E-01 -1.3284E+00  5.0827E-03 -3.8684E+00 -6.3041E+00  5.2619E-01 -3.1684E-01
             1.3192E+00
 GRADIENT:   3.2784E+01  1.3461E+01  2.3433E+01  1.6899E+01  1.2370E+02  2.7005E+00  2.4741E-03  0.0000E+00  9.5850E+00 -8.8591E-01
             1.0725E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1707.71746068632        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      835
 NPARAMETR:  1.0213E+00  2.9043E-01  2.2504E-01  1.0867E+00  2.4004E-01  9.0762E-01  1.0527E-02  1.0000E-02  1.5343E+00  6.6559E-01
             3.3958E+00
 PARAMETER:  1.2110E-01 -1.1364E+00 -1.3915E+00  1.8316E-01 -1.3270E+00  3.0678E-03 -4.4538E+00 -6.3041E+00  5.2804E-01 -3.0709E-01
             1.3225E+00
 GRADIENT:   3.4829E+01  1.1692E+01  2.2492E+01  1.3745E+01  1.2684E+02  2.1746E+00  9.2004E-04  0.0000E+00  1.0395E+01  1.2166E+00
             1.5815E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1707.71850289154        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1017
 NPARAMETR:  1.0210E+00  2.9113E-01  2.2605E-01  1.0879E+00  2.4074E-01  9.0749E-01  1.0000E-02  1.0000E-02  1.5305E+00  6.6571E-01
             3.3958E+00
 PARAMETER:  1.2082E-01 -1.1340E+00 -1.3870E+00  1.8426E-01 -1.3240E+00  2.9247E-03 -7.2236E+00 -6.3041E+00  5.2561E-01 -3.0689E-01
             1.3225E+00
 GRADIENT:  -2.8197E-01  3.7352E-02 -2.6600E-02 -6.3949E-02 -5.7439E-02 -1.0554E-02  0.0000E+00  0.0000E+00 -1.2978E-01  1.4089E-02
            -5.7951E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1707.71860856226        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1117
 NPARAMETR:  1.0215E+00  2.9172E-01  2.2633E-01  1.0873E+00  2.4062E-01  9.0755E-01  1.0000E-02  1.0000E-02  1.5393E+00  6.6584E-01
             3.3940E+00
 PARAMETER:  1.2094E-01 -1.1338E+00 -1.3865E+00  1.8449E-01 -1.3237E+00  2.9441E-03 -7.5416E+00 -6.3041E+00  5.2605E-01 -3.0704E-01
             1.3226E+00
 GRADIENT:  -6.1752E-02 -9.8562E-02 -6.0212E-02  5.4170E-02  2.4653E-01 -1.4425E-03  0.0000E+00  0.0000E+00 -1.2267E-01 -5.8477E-03
             7.0429E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1117
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5357E-03 -5.2320E-05  1.8874E-04 -8.6091E-03  3.3301E-03
 SE:             2.8661E-02  1.5349E-04  1.9700E-04  2.6897E-02  2.5049E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5727E-01  7.3321E-01  3.3803E-01  7.4891E-01  8.9424E-01

 ETASHRINKSD(%)  3.9823E+00  9.9486E+01  9.9340E+01  9.8924E+00  1.6082E+01
 ETASHRINKVR(%)  7.8061E+00  9.9997E+01  9.9996E+01  1.8806E+01  2.9577E+01
 EBVSHRINKSD(%)  3.9298E+00  9.9473E+01  9.9315E+01  7.0340E+00  1.6703E+01
 EBVSHRINKVR(%)  7.7052E+00  9.9997E+01  9.9995E+01  1.3573E+01  3.0616E+01
 RELATIVEINF(%)  9.2128E+01  4.7548E-04  3.4364E-04  4.3644E+01  2.7733E+00
 EPSSHRINKSD(%)  2.1378E+01
 EPSSHRINKVR(%)  3.8186E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1707.7186085622573     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -604.99236871665016     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.98
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1707.719       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.91E-01  2.26E-01  1.09E+00  2.41E-01  9.08E-01  1.00E-02  1.00E-02  1.53E+00  6.66E-01  3.40E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.20E+03
 
 TH 2
+        2.14E+01  3.62E+03
 
 TH 3
+        7.50E+00  3.27E+03  9.48E+03
 
 TH 4
+       -2.06E+01  7.08E+01 -3.01E+02  3.44E+02
 
 TH 5
+        6.72E+01 -8.03E+03 -1.51E+04 -2.42E+02  3.03E+04
 
 TH 6
+        2.43E+00 -7.58E+00  3.91E+01 -7.89E+00 -6.58E+00  2.03E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.12E+01 -4.20E+01  8.17E+01 -2.07E+00  3.16E+01  1.32E+00  0.00E+00  0.00E+00  5.66E+01
 
 TH10
+       -5.17E+00 -6.76E+01 -1.49E+01  7.34E+00  1.03E+02  1.10E+00  0.00E+00  0.00E+00  4.96E+00  2.18E+02
 
 TH11
+       -2.26E+01 -1.05E+01 -3.72E+01 -5.83E+00  7.46E+01  3.55E+00  0.00E+00  0.00E+00  2.89E+00  1.49E+01  5.70E+01
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       32.480
Stop Time:
Thu Sep 30 06:51:33 CDT 2021
