Wed Sep 29 23:00:10 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat4.csv ignore=@
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
Current Date:       29 SEP 2021
Days until program expires : 200
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
 RAW OUTPUT FILE (FILE): m4.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1125.17143842067        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8379E+02  4.3121E+01  1.2898E+02  4.1941E+01  1.7962E+02  2.2544E+01 -2.4734E+01 -2.3172E+02 -2.9627E+01 -8.7243E+01
            -1.4940E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1660.10575856884        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0103E+00  7.8231E-01  5.3301E-01  1.1816E+00  5.9938E-01  1.0461E+00  8.6919E-01  1.3632E+00  1.2047E+00  7.8358E-01
             2.5436E+00
 PARAMETER:  1.1022E-01 -1.4550E-01 -5.2921E-01  2.6686E-01 -4.1186E-01  1.4508E-01 -4.0190E-02  4.0983E-01  2.8626E-01 -1.4388E-01
             1.0336E+00
 GRADIENT:   1.1641E+02 -7.2919E+00 -2.6430E+01  1.1335E+02  6.7632E+01  8.7445E+00  6.4046E+00  1.7050E+01  3.6135E+01  2.2557E+01
             1.1722E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1686.95790825588        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      243             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8069E-01  6.3671E-01  3.4215E-01  1.0435E+00  4.2056E-01  1.0707E+00  4.6245E-01  1.5648E+00  1.1236E+00  3.2970E-01
             2.3412E+00
 PARAMETER:  8.0503E-02 -3.5144E-01 -9.7252E-01  1.4258E-01 -7.6618E-01  1.6828E-01 -6.7121E-01  5.4775E-01  2.1654E-01 -1.0096E+00
             9.5067E-01
 GRADIENT:   6.2113E+01 -5.3772E+01 -5.3756E+00 -7.7804E+01  1.3810E+02  1.8868E+01  1.2792E-02  3.5334E+00 -3.0335E+00  4.3824E+00
             1.3376E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1692.71003652769        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      399
 NPARAMETR:  9.8140E-01  6.8816E-01  3.4226E-01  1.0846E+00  4.2063E-01  1.0629E+00  8.7146E-01  1.5645E+00  1.0900E+00  1.8949E-01
             2.3406E+00
 PARAMETER:  8.1227E-02 -2.7374E-01 -9.7218E-01  1.8119E-01 -7.6599E-01  1.6096E-01 -3.7580E-02  5.4755E-01  1.8622E-01 -1.5634E+00
             9.5039E-01
 GRADIENT:  -7.1297E+00  6.7635E-01  6.1814E+00 -5.3416E+00  1.4303E+01  1.2054E+00  3.1835E+00  6.9538E+00 -1.8543E+00  1.4953E+00
             1.2060E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1693.95093661434        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      574
 NPARAMETR:  9.8711E-01  6.8232E-01  3.4222E-01  1.0860E+00  4.2051E-01  1.0460E+00  8.1285E-01  1.5643E+00  1.1133E+00  3.2486E-02
             2.3316E+00
 PARAMETER:  8.7023E-02 -2.8226E-01 -9.7231E-01  1.8248E-01 -7.6630E-01  1.4493E-01 -1.0721E-01  5.4742E-01  2.0735E-01 -3.3269E+00
             9.4656E-01
 GRADIENT:   4.6636E+00 -4.4778E+00  3.6503E+00 -6.2984E+00  2.4393E+01 -4.8631E+00  1.2861E+00  2.9284E+00  1.9564E+00  3.5738E-02
             1.1731E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1706.69111959952        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  9.7706E-01  6.7799E-01  3.4014E-01  1.1033E+00  4.1419E-01  1.1231E+00  8.3898E-01  1.5539E+00  1.1121E+00  1.0000E-02
             1.9222E+00
 PARAMETER:  7.6794E-02 -2.8862E-01 -9.7839E-01  1.9829E-01 -7.8144E-01  2.1611E-01 -7.5563E-02  5.4074E-01  2.0625E-01 -8.3833E+01
             7.5348E-01
 GRADIENT:  -5.0398E+00  2.2279E+01  1.8064E+01  2.5377E+01  2.4709E+01  1.9611E+01 -1.5049E+00 -2.7234E+01 -6.2076E+00  0.0000E+00
             8.7013E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1709.59211796011        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  9.7719E-01  6.0677E-01  3.4012E-01  1.1029E+00  3.8385E-01  1.0198E+00  8.6322E-01  1.5535E+00  1.1127E+00  1.0000E-02
             1.9212E+00
 PARAMETER:  7.6923E-02 -3.9961E-01 -9.7845E-01  1.9794E-01 -8.5750E-01  1.1958E-01 -4.7083E-02  5.4051E-01  2.0679E-01 -8.4176E+01
             7.5294E-01
 GRADIENT:   8.8494E+01  4.4639E+01  7.6388E+01  3.8244E+01  5.1999E+01 -2.8380E+00 -7.2254E-01 -2.2025E+01  3.5299E+00  0.0000E+00
             3.0834E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1710.51503538795        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:     1002
 NPARAMETR:  9.7770E-01  5.4950E-01  3.3834E-01  1.1029E+00  3.7594E-01  1.0545E+00  8.6629E-01  1.5534E+00  1.1127E+00  1.0000E-02
             1.9214E+00
 PARAMETER:  7.7452E-02 -4.9875E-01 -9.8370E-01  1.9796E-01 -8.7832E-01  1.5305E-01 -4.3534E-02  5.4046E-01  2.0677E-01 -8.4176E+01
             7.5305E-01
 GRADIENT:   9.1713E+01 -2.2104E-02  5.7704E+01  7.0609E+00  1.2118E+02  1.4812E+01  4.2752E-01 -1.9959E+01  5.2458E+00  0.0000E+00
             3.7333E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1710.85170590746        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1179
 NPARAMETR:  9.7770E-01  5.6803E-01  3.3825E-01  1.1029E+00  3.7592E-01  1.0637E+00  8.7222E-01  1.5534E+00  1.1127E+00  1.0000E-02
             1.9213E+00
 PARAMETER:  7.7450E-02 -4.6559E-01 -9.8396E-01  1.9796E-01 -8.7837E-01  1.6179E-01 -3.6713E-02  5.4048E-01  2.0677E-01 -8.4176E+01
             7.5302E-01
 GRADIENT:  -3.2269E+00  9.7217E-01  4.5131E+01 -2.7024E+01  3.1062E+00  1.2965E-01 -6.2841E-03 -2.4469E+01 -2.0775E+00  0.0000E+00
             2.9272E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1711.12403535670        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1364
 NPARAMETR:  9.7901E-01  5.6591E-01  3.3795E-01  1.1211E+00  3.7631E-01  1.0632E+00  8.7328E-01  1.5525E+00  1.1130E+00  1.0000E-02
             1.9212E+00
 PARAMETER:  7.8787E-02 -4.6931E-01 -9.8484E-01  2.1429E-01 -8.7735E-01  1.6128E-01 -3.5499E-02  5.3985E-01  2.0703E-01 -8.4176E+01
             7.5294E-01
 GRADIENT:  -8.9641E-01  3.2919E+00  4.0543E+01 -3.0180E+00  6.5832E+00 -1.7553E-01  5.9205E-02 -2.4259E+01 -1.2176E+00  0.0000E+00
             2.8679E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1713.23759629223        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1530
 NPARAMETR:  9.7371E-01  5.6592E-01  3.2241E-01  1.1172E+00  3.7674E-01  1.0596E+00  8.4337E-01  1.6744E+00  1.1158E+00  1.0000E-02
             1.9231E+00
 PARAMETER:  7.3354E-02 -4.6931E-01 -1.0319E+00  2.1082E-01 -8.7621E-01  1.5789E-01 -7.0349E-02  6.1546E-01  2.0958E-01 -8.4176E+01
             7.5393E-01
 GRADIENT:  -1.0761E+01 -1.7250E+01  6.5662E+00 -4.4197E+00  5.9150E+01 -1.1798E+00  8.8721E-01 -5.5941E+00 -2.2945E+00  0.0000E+00
             3.6657E+01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1713.24373861556        NO. OF FUNC. EVALS.:  69
 CUMULATIVE NO. OF FUNC. EVALS.:     1599
 NPARAMETR:  9.7392E-01  5.6628E-01  3.2190E-01  1.1176E+00  3.7719E-01  1.0613E+00  8.4264E-01  1.6737E+00  1.1158E+00  1.0000E-02
             1.9251E+00
 PARAMETER:  7.3432E-02 -4.6931E-01 -1.0321E+00  2.1086E-01 -8.7621E-01  1.5793E-01 -7.1349E-02  6.1591E-01  2.0983E-01 -8.4176E+01
             7.5393E-01
 GRADIENT:  -3.1944E+04 -6.8264E+03  3.1034E+03 -1.5152E+04 -7.2194E+03 -1.0935E+00 -6.3870E+04  1.0350E+04  3.0438E+04  0.0000E+00
            -8.5435E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1599
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.1905E-03 -1.1064E-02 -4.8276E-03  4.0100E-03 -3.1687E-04
 SE:             2.9704E-02  1.2922E-02  2.6606E-02  2.8229E-02  2.8971E-04
 N:                     100         100         100         100         100

 P VAL.:         8.3491E-01  3.9187E-01  8.5602E-01  8.8704E-01  2.7407E-01

 ETASHRINKSD(%)  4.8883E-01  5.6710E+01  1.0867E+01  5.4297E+00  9.9029E+01
 ETASHRINKVR(%)  9.7528E-01  8.1260E+01  2.0553E+01  1.0565E+01  9.9991E+01
 EBVSHRINKSD(%)  1.0807E+00  5.4870E+01  1.2745E+01  6.2446E+00  9.8959E+01
 EBVSHRINKVR(%)  2.1498E+00  7.9633E+01  2.3866E+01  1.2099E+01  9.9989E+01
 RELATIVEINF(%)  9.7787E+01  3.5663E+00  2.0627E+01  6.2771E+01  1.1305E-03
 EPSSHRINKSD(%)  3.5630E+01
 EPSSHRINKVR(%)  5.8565E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1713.2437386155605     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -794.30520541088777     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.02
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1713.244       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.74E-01  5.66E-01  3.22E-01  1.12E+00  3.77E-01  1.06E+00  8.43E-01  1.68E+00  1.12E+00  1.00E-02  1.92E+00
 


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
+        1.68E+07
 
 TH 2
+        6.17E+06  4.53E+06
 
 TH 3
+        9.86E+06 -4.04E+03  2.89E+06
 
 TH 4
+        3.43E+03  1.57E+03  2.04E+06  2.88E+06
 
 TH 5
+        4.97E+06 -1.82E+06 -1.45E+06  2.06E+06  1.47E+06
 
 TH 6
+       -8.46E+02 -3.22E+02  2.51E+02 -3.59E+02 -2.37E+02  1.71E+02
 
 TH 7
+        2.25E+03  8.32E+02 -6.58E+02  9.29E+02  6.59E+02 -9.91E+02  2.25E+07
 
 TH 8
+        4.33E+03  1.59E+03  4.64E+05  1.79E+03  4.70E+05  8.21E+01 -1.83E+06  1.50E+05
 
 TH 9
+       -1.60E+03 -6.14E+02 -4.10E+06 -1.42E+03 -2.07E+06  3.57E+02 -9.22E+02 -1.80E+03  2.91E+06
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        2.89E+04  1.06E+04  6.62E+05  1.19E+04  3.38E+05 -5.65E+01  1.32E+06 -1.08E+05 -1.20E+04  0.00E+00  7.80E+04
 
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
 #CPUT: Total CPU Time in Seconds,       34.180
Stop Time:
Wed Sep 29 23:00:46 CDT 2021
