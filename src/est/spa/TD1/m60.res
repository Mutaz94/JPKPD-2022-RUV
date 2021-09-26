Sat Sep 25 12:59:58 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat60.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1663.58951315805        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8176E+01 -7.9579E+01 -3.7809E+01 -5.3181E+01  3.7168E+01 -2.5377E+00 -1.7771E+01  6.0906E+00  4.8375E+00  1.6664E+01
            -2.2543E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1670.39674143981        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      111
 NPARAMETR:  1.0071E+00  1.2051E+00  1.1418E+00  9.1744E-01  1.0832E+00  1.0066E+00  1.2221E+00  9.7221E-01  9.2862E-01  8.0458E-01
             1.1180E+00
 PARAMETER:  1.0710E-01  2.8653E-01  2.3260E-01  1.3827E-02  1.7992E-01  1.0657E-01  3.0056E-01  7.1816E-02  2.5947E-02 -1.1744E-01
             2.1150E-01
 GRADIENT:   7.7173E+00  1.2644E+01  3.5976E+01 -3.0798E+01 -2.1004E+01 -8.5012E+00  5.4794E+00 -1.0292E+01 -5.8130E+00 -9.3903E+00
             9.1698E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1670.93756244557        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0194E+00  1.2087E+00  1.1112E+00  9.2377E-01  1.0731E+00  1.0372E+00  1.2396E+00  9.9158E-01  9.2900E-01  7.9193E-01
             1.1268E+00
 PARAMETER:  1.1920E-01  2.8955E-01  2.0544E-01  2.0708E-02  1.7058E-01  1.3648E-01  3.1475E-01  9.1544E-02  2.6355E-02 -1.3329E-01
             2.1938E-01
 GRADIENT:   3.1368E+01  1.9518E+01  2.9751E+01 -1.7293E+01 -1.6789E+01  2.7722E+00  7.6294E+00 -8.4234E+00 -3.9683E+00 -8.0124E+00
             1.4079E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1672.87535143946        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      490
 NPARAMETR:  1.0019E+00  1.2082E+00  1.0222E+00  9.2396E-01  1.0731E+00  1.0268E+00  1.1799E+00  9.9224E-01  9.5456E-01  8.5892E-01
             1.1268E+00
 PARAMETER:  1.0191E-01  2.8915E-01  1.2194E-01  2.0914E-02  1.7056E-01  1.2649E-01  2.6547E-01  9.2213E-02  5.3500E-02 -5.2076E-02
             2.1937E-01
 GRADIENT:  -5.6204E+00  1.6707E-01 -2.3655E+00  6.7702E+00  1.9967E+01 -9.5146E-01  2.1015E+00 -1.4388E+00  9.7813E-01  1.2158E+00
             2.1194E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1672.91032006816        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      666
 NPARAMETR:  1.0046E+00  1.2082E+00  1.0282E+00  9.2396E-01  1.0731E+00  1.0293E+00  1.1658E+00  9.9224E-01  9.5227E-01  8.5115E-01
             1.1268E+00
 PARAMETER:  1.0458E-01  2.8915E-01  1.2781E-01  2.0914E-02  1.7056E-01  1.2884E-01  2.5341E-01  9.2213E-02  5.1094E-02 -6.1165E-02
             2.1937E-01
 GRADIENT:   2.0200E-03  4.8786E-01  1.4466E-03  5.7587E+00  1.8070E+01 -2.6612E-03  1.6045E-03 -2.0350E+00 -4.7801E-04 -7.2116E-04
             2.0438E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1673.55374679664        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      767
 NPARAMETR:  1.0036E+00  1.2069E+00  1.0271E+00  9.2148E-01  1.0649E+00  1.0279E+00  1.1655E+00  1.0250E+00  9.5214E-01  8.5035E-01
             1.0911E+00
 PARAMETER:  1.0359E-01  2.8808E-01  1.2677E-01  1.8223E-02  1.6291E-01  1.2756E-01  2.5313E-01  1.2466E-01  5.0960E-02 -6.2103E-02
             1.8723E-01
 GRADIENT:   3.7730E+01  1.3096E+01  3.8269E+00  5.7642E+00  9.2996E+00  9.4394E+00  1.1373E+00 -1.8726E+00  2.2900E-01  6.6491E-01
             8.7786E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1673.57593145285        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:      958
 NPARAMETR:  1.0044E+00  1.2068E+00  1.0163E+00  9.2145E-01  1.0650E+00  1.0298E+00  1.1682E+00  1.0249E+00  9.5225E-01  8.4364E-01
             1.0912E+00
 PARAMETER:  1.0442E-01  2.8801E-01  1.1616E-01  1.8198E-02  1.6296E-01  1.2938E-01  2.5551E-01  1.2463E-01  5.1071E-02 -7.0034E-02
             1.8727E-01
 GRADIENT:   3.8073E-02 -1.4554E+00 -6.3176E-02  2.6803E+00  1.5713E+01  7.9555E-03  6.6291E-02 -1.4303E+00  1.2328E-02 -1.0903E-02
             8.7317E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1673.59619853818        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1115
 NPARAMETR:  1.0044E+00  1.2069E+00  1.0164E+00  9.2137E-01  1.0645E+00  1.0298E+00  1.1678E+00  1.0274E+00  9.5226E-01  8.4379E-01
             1.0901E+00
 PARAMETER:  1.0440E-01  2.8808E-01  1.1631E-01  1.8103E-02  1.6246E-01  1.2935E-01  2.5509E-01  1.2707E-01  5.1084E-02 -6.9848E-02
             1.8623E-01
 GRADIENT:   3.5928E-03 -1.2847E+00  2.0692E-01  2.4407E+00  1.4842E+01 -7.7308E-03  1.6225E-02 -1.4043E+00 -1.3295E-05  9.4096E-02
             8.3634E+00

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1673.59619853818        NO. OF FUNC. EVALS.:  33
 CUMULATIVE NO. OF FUNC. EVALS.:     1148
 NPARAMETR:  1.0044E+00  1.2069E+00  1.0164E+00  9.2136E-01  1.0645E+00  1.0298E+00  1.1678E+00  1.0274E+00  9.5226E-01  8.4371E-01
             1.0901E+00
 PARAMETER:  1.0440E-01  2.8808E-01  1.1631E-01  1.8103E-02  1.6246E-01  1.2935E-01  2.5509E-01  1.2707E-01  5.1084E-02 -6.9848E-02
             1.8623E-01
 GRADIENT:  -1.4911E-02  8.3272E+05  1.9608E-01  1.1994E+06 -1.4766E+06 -1.7260E-02  1.5704E-02  1.8878E+06 -5.7898E-04  9.0962E-02
            -1.2885E+06
 NUMSIGDIG:         4.2         3.3         2.0         3.3         3.3         3.5         3.4         3.3         4.4         1.7
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1148
 NO. OF SIG. DIGITS IN FINAL EST.:  1.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5825E-04 -4.9695E-03 -2.9047E-02 -6.2921E-04 -3.6651E-02
 SE:             2.9841E-02  2.2633E-02  1.2999E-02  2.1837E-02  2.0078E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9309E-01  8.2621E-01  2.5441E-02  9.7701E-01  6.7940E-02

 ETASHRINKSD(%)  2.7762E-02  2.4175E+01  5.6453E+01  2.6844E+01  3.2736E+01
 ETASHRINKVR(%)  5.5517E-02  4.2506E+01  8.1036E+01  4.6482E+01  5.4755E+01
 EBVSHRINKSD(%)  4.7465E-01  2.4356E+01  6.0683E+01  2.7734E+01  3.0916E+01
 EBVSHRINKVR(%)  9.4705E-01  4.2780E+01  8.4541E+01  4.7777E+01  5.2274E+01
 RELATIVEINF(%)  9.8619E+01  2.8023E+00  1.6386E+00  2.6116E+00  8.8477E+00
 EPSSHRINKSD(%)  4.4448E+01
 EPSSHRINKVR(%)  6.9140E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1673.5961985381821     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -938.44537197444390     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.29
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.00
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1673.596       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.21E+00  1.02E+00  9.21E-01  1.06E+00  1.03E+00  1.17E+00  1.03E+00  9.52E-01  8.44E-01  1.09E+00
 


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
+        1.03E+03
 
 TH 2
+        1.97E+03  4.96E+08
 
 TH 3
+        7.57E+00  5.01E+04  4.29E+09
 
 TH 4
+        7.44E+03 -3.14E+04  1.89E+05  7.06E+09
 
 TH 5
+       -3.97E+03  1.38E+04 -2.93E+09  5.29E+04  2.01E+09
 
 TH 6
+       -1.38E-01  3.62E+03  1.60E+00  1.37E+04 -7.27E+03  1.84E+02
 
 TH 7
+        8.54E-01  2.35E+03 -1.70E+09  8.76E+03 -4.68E+03  3.91E-01  5.77E+01
 
 TH 8
+        5.25E+03 -1.60E+05  3.89E+09 -6.05E+05 -2.66E+09  9.63E+03  6.20E+03  3.52E+09
 
 TH 9
+       -3.17E-01 -2.85E+04 -1.32E+01 -1.07E+05  5.73E+04  3.68E-02  1.72E+01 -7.58E+04  7.40E+01
 
 TH10
+       -1.30E+00 -7.73E+03 -1.04E+01 -2.92E+04  1.55E+04 -3.18E-01  7.88E+00 -2.06E+04  3.61E+00  7.39E+01
 
 TH11
+       -2.82E+09 -4.43E+05 -2.50E+09 -1.67E+06  1.71E+09 -2.22E+09  9.92E+08 -2.26E+09  3.10E+09  3.50E+09  1.46E+09
 
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
 #CPUT: Total CPU Time in Seconds,       21.359
Stop Time:
Sat Sep 25 13:00:20 CDT 2021
