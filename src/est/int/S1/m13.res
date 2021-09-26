Fri Sep 24 23:04:20 CDT 2021
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
$DATA ../../../../data/int/S1/dat13.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 RAW OUTPUT FILE (FILE): m13.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3708.25146503808        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2504E+01 -9.0272E+00 -4.7360E+01 -4.7859E+01  1.6235E+01 -1.1499E+00 -1.4041E+01 -2.0525E+01  1.2254E+00  2.2335E+01
            -3.2176E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3737.27814105396        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0007E+00  1.0222E+00  1.1646E+00  1.0230E+00  1.0713E+00  1.0052E+00  1.2056E+00  1.2905E+00  9.2666E-01  7.5215E-01
             1.1519E+00
 PARAMETER:  1.0068E-01  1.2193E-01  2.5235E-01  1.2276E-01  1.6891E-01  1.0518E-01  2.8701E-01  3.5506E-01  2.3833E-02 -1.8481E-01
             2.4143E-01
 GRADIENT:   4.1270E+01 -1.2395E+00 -2.1060E+01  1.6059E+01  5.9939E+01  1.0692E-01  3.4939E+00 -3.9625E+00 -7.2856E+00 -8.9535E+00
             2.4092E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3737.86951144180        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      184
 NPARAMETR:  9.9799E-01  1.0073E+00  1.1676E+00  1.0280E+00  1.0554E+00  1.0036E+00  1.2356E+00  1.3513E+00  9.3135E-01  7.3772E-01
             1.1471E+00
 PARAMETER:  9.7986E-02  1.0730E-01  2.5491E-01  1.2758E-01  1.5394E-01  1.0359E-01  3.1158E-01  4.0110E-01  2.8879E-02 -2.0419E-01
             2.3725E-01
 GRADIENT:   3.5640E+01 -1.7157E+00 -1.8112E+01  1.0139E+01  4.1119E+01 -5.0165E-01  6.1493E+00  7.9678E-01 -5.2098E+00 -6.3577E+00
             2.0494E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3738.23038443277        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  9.9329E-01  1.0047E+00  1.1695E+00  1.0286E+00  1.0489E+00  1.0032E+00  1.2269E+00  1.3388E+00  9.4166E-01  7.4492E-01
             1.1439E+00
 PARAMETER:  9.3264E-02  1.0465E-01  2.5660E-01  1.2824E-01  1.4771E-01  1.0320E-01  3.0449E-01  3.9175E-01  3.9890E-02 -1.9448E-01
             2.3448E-01
 GRADIENT:   2.5545E+01 -9.4630E-01 -1.2483E+01  7.6588E+00  2.9373E+01 -6.2197E-01  4.9380E+00  2.3214E-01 -2.4221E+00 -4.6013E+00
             1.5355E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3738.38531095193        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  9.9382E-01  1.0067E+00  1.1762E+00  1.0286E+00  1.0503E+00  1.0049E+00  1.2200E+00  1.3461E+00  9.4242E-01  7.5375E-01
             1.1424E+00
 PARAMETER:  9.3797E-02  1.0667E-01  2.6232E-01  1.2817E-01  1.4905E-01  1.0490E-01  2.9887E-01  3.9721E-01  4.0693E-02 -1.8270E-01
             2.3317E-01
 GRADIENT:  -6.3888E+00 -5.7621E+00 -1.1726E+01 -1.8275E+00  1.9826E+01 -3.4685E+00  2.3084E+00 -6.6497E-02 -2.5607E+00 -4.0610E+00
             1.2355E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3738.40006283010        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      638
 NPARAMETR:  9.9380E-01  1.0067E+00  1.1763E+00  1.0285E+00  1.0503E+00  1.0138E+00  1.2200E+00  1.3461E+00  9.4240E-01  7.5392E-01
             1.1425E+00
 PARAMETER:  9.3782E-02  1.0668E-01  2.6236E-01  1.2814E-01  1.4912E-01  1.1369E-01  2.9885E-01  3.9721E-01  4.0678E-02 -1.8247E-01
             2.3320E-01
 GRADIENT:  -6.2852E+00 -5.8319E+00 -1.1723E+01 -1.8845E+00  1.9842E+01  2.0210E-02  2.3111E+00 -6.8556E-02 -2.5649E+00 -4.0365E+00
             1.2465E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3738.59894313106        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      782
 NPARAMETR:  9.9634E-01  1.0121E+00  1.2031E+00  1.0285E+00  1.0400E+00  1.0135E+00  1.2042E+00  1.3476E+00  9.5244E-01  7.8248E-01
             1.1356E+00
 PARAMETER:  9.6335E-02  1.1198E-01  2.8489E-01  1.2808E-01  1.3920E-01  1.1339E-01  2.8582E-01  3.9835E-01  5.1274E-02 -1.4529E-01
             2.2714E-01
 GRADIENT:  -6.6203E-01  6.0119E+00  6.6451E+00 -7.8447E+00 -1.9024E+01 -2.8426E-02  1.5379E+00 -1.3738E+00  8.9040E-01  6.8991E-01
             3.0734E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3738.60972232141        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  9.9634E-01  1.0121E+00  1.2031E+00  1.0285E+00  1.0400E+00  1.0102E+00  1.1961E+00  1.3541E+00  9.5017E-01  7.8245E-01
             1.1355E+00
 PARAMETER:  9.6335E-02  1.1198E-01  2.8490E-01  1.2808E-01  1.3919E-01  1.1013E-01  2.7904E-01  4.0312E-01  4.8884E-02 -1.4533E-01
             2.2708E-01
 GRADIENT:   3.2984E+01  1.2156E+01  7.7756E+00  1.6533E+00 -1.4192E+01  2.4913E+00  2.2561E+00 -5.0694E-01  9.9888E-01  6.4395E-01
             5.6855E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -3738.60972232141        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      952
 NPARAMETR:  9.9634E-01  1.0121E+00  1.2031E+00  1.0285E+00  1.0400E+00  1.0103E+00  1.1961E+00  1.3541E+00  9.5017E-01  7.8245E-01
             1.1355E+00
 PARAMETER:  9.6335E-02  1.1198E-01  2.8490E-01  1.2808E-01  1.3919E-01  1.1013E-01  2.7904E-01  4.0312E-01  4.8884E-02 -1.4533E-01
             2.2708E-01
 GRADIENT:  -4.7754E+05 -8.5291E+05  1.6764E+05 -5.0823E+01 -2.6797E+01 -1.3008E+00 -3.4228E+05 -1.1854E+05  4.7753E+05 -6.5719E+05
            -4.2060E+05
 NUMSIGDIG:         3.3         3.3         3.3         7.8         8.0         1.2         3.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      952
 NO. OF SIG. DIGITS IN FINAL EST.:  1.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1068E-04 -2.0850E-02 -3.0498E-02  1.4906E-02 -2.8275E-02
 SE:             2.9975E-02  2.4316E-02  2.0127E-02  2.6846E-02  2.1227E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7842E-01  3.9119E-01  1.2971E-01  5.7873E-01  1.8286E-01

 ETASHRINKSD(%)  1.0000E-10  1.8540E+01  3.2571E+01  1.0063E+01  2.8885E+01
 ETASHRINKVR(%)  1.0000E-10  3.3642E+01  5.4533E+01  1.9114E+01  4.9427E+01
 EBVSHRINKSD(%)  3.2473E-01  1.7917E+01  3.4154E+01  1.0905E+01  2.8604E+01
 EBVSHRINKVR(%)  6.4841E-01  3.2625E+01  5.6643E+01  2.0620E+01  4.9027E+01
 RELATIVEINF(%)  9.9349E+01  3.6599E+01  2.8049E+01  5.4148E+01  2.0548E+01
 EPSSHRINKSD(%)  2.1150E+01
 EPSSHRINKVR(%)  3.7826E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3738.6097223214088     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2084.5203625529980     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.24
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3738.610       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  1.01E+00  1.20E+00  1.03E+00  1.04E+00  1.01E+00  1.20E+00  1.35E+00  9.50E-01  7.82E-01  1.14E+00
 


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
+        2.41E+09
 
 TH 2
+        2.11E+09  3.72E+09
 
 TH 3
+        5.45E-01  3.62E+01  4.06E+08
 
 TH 4
+       -1.82E+09  2.10E+02  5.29E+08  1.38E+09
 
 TH 5
+       -2.65E+00 -3.15E+02 -4.81E+08 -1.25E+09  2.28E+09
 
 TH 6
+       -2.42E+03 -2.13E+03  7.04E+02 -4.80E+00 -2.06E+00  1.95E+02
 
 TH 7
+        7.18E+08 -2.15E+03  7.15E+02  4.01E+00  6.73E+00 -7.22E+02  2.14E+08
 
 TH 8
+        5.84E+05  5.13E+05 -1.70E+05  1.06E+01 -2.08E+01 -4.43E+02  1.74E+05  8.02E+07
 
 TH 9
+        8.29E+04  7.28E+04 -2.41E+04  5.96E+01  7.36E+00  2.54E+03  2.47E+04 -4.61E+08  2.64E+09
 
 TH10
+        2.11E+09  2.32E+03 -7.72E+02  2.21E+00 -4.59E+01 -2.12E+03  6.29E+08  5.11E+05  7.26E+04  1.85E+09
 
 TH11
+        7.26E+04  6.38E+04 -2.11E+04 -3.38E+00 -1.77E+01 -9.34E+02  2.17E+04  1.70E+08 -9.75E+08  6.36E+04  3.59E+08
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       38.739
Stop Time:
Fri Sep 24 23:05:01 CDT 2021
