Wed Sep 29 23:07:16 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat19.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1138.57560843397        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7300E+02  3.0367E+01  1.3097E+02  2.1527E+01  9.8788E+01  6.1708E+01 -2.2340E+01 -1.4068E+02  3.3849E+01 -9.2279E+01
            -1.7623E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1761.21571588869        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0238E+00  1.0166E+00  9.2630E-01  1.0293E+00  9.5513E-01  8.2550E-01  9.2162E-01  9.6893E-01  7.7703E-01  1.1100E+00
             2.3723E+00
 PARAMETER:  1.2350E-01  1.1642E-01  2.3440E-02  1.2889E-01  5.4094E-02 -9.1771E-02  1.8382E-02  6.8434E-02 -1.5227E-01  2.0435E-01
             9.6387E-01
 GRADIENT:   4.1740E+01 -7.7436E+00 -5.2720E+00 -6.5518E+00  8.9686E+00 -3.3880E+01 -1.6674E+00  8.6003E+00  9.1862E-01 -7.9353E-01
            -3.7009E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1767.49577352227        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0170E+00  7.9741E-01  8.0090E-01  1.1742E+00  7.8489E-01  8.8267E-01  1.2759E+00  4.1490E-01  6.0657E-01  1.0338E+00
             2.3201E+00
 PARAMETER:  1.1687E-01 -1.2639E-01 -1.2202E-01  2.6062E-01 -1.4221E-01 -2.4799E-02  3.4368E-01 -7.7973E-01 -3.9993E-01  1.3327E-01
             9.4159E-01
 GRADIENT:   2.5941E+01  2.7871E+01 -6.7535E+00  1.0070E+02  1.0107E+01 -8.0251E+00 -6.2150E+00  2.0810E+00 -1.2600E+01  4.1795E+00
            -1.5722E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1770.44215326071        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.0307E+00  7.4645E-01  6.9159E-01  1.1866E+00  6.9707E-01  9.0850E-01  1.5340E+00  3.1696E-01  6.0368E-01  8.9175E-01
             2.3496E+00
 PARAMETER:  1.3024E-01 -1.9243E-01 -2.6876E-01  2.7110E-01 -2.6086E-01  4.0376E-03  5.2791E-01 -1.0490E+00 -4.0471E-01 -1.4572E-02
             9.5425E-01
 GRADIENT:  -2.8680E+01  2.5353E+01 -9.5912E+00  4.0064E+01  1.6712E+01 -2.3647E+00  7.0366E+00  1.1211E+00 -1.0260E+01 -2.4476E+00
            -7.4158E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1777.82533953273        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  1.0335E+00  3.7719E-01  5.9425E-01  1.3155E+00  5.2257E-01  9.0796E-01  2.2806E+00  1.0000E-02  6.4200E-01  8.1165E-01
             2.3191E+00
 PARAMETER:  1.3294E-01 -8.7501E-01 -4.2046E-01  3.7423E-01 -5.4900E-01  3.4470E-03  9.2443E-01 -4.7226E+00 -3.4317E-01 -1.0868E-01
             9.4117E-01
 GRADIENT:  -4.6197E+00 -7.9733E-01  1.1490E+01 -3.7282E+01 -2.5995E+00  1.1339E-01  3.1127E+00  0.0000E+00 -3.4830E+00  1.0260E+00
            -7.4326E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1778.83013381255        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      647
 NPARAMETR:  1.0370E+00  3.2978E-01  5.0149E-01  1.3160E+00  4.5676E-01  9.1080E-01  2.3557E+00  1.0000E-02  6.7580E-01  7.4307E-01
             2.3116E+00
 PARAMETER:  1.3630E-01 -1.0093E+00 -5.9016E-01  3.7460E-01 -6.8360E-01  6.5710E-03  9.5683E-01 -5.9049E+00 -2.9186E-01 -1.9697E-01
             9.3793E-01
 GRADIENT:   2.5183E+00 -3.0858E+00 -6.9877E+00 -4.8911E-01  1.2596E+01  5.9901E-01  1.2447E+00  0.0000E+00 -3.6395E-01  5.3171E-01
             1.7718E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1784.25202380919        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      830
 NPARAMETR:  1.0473E+00  5.0890E-01  2.7928E-01  1.1229E+00  3.3966E-01  9.2857E-01  1.3279E+00  1.0000E-02  8.2909E-01  5.7924E-01
             2.1933E+00
 PARAMETER:  1.4624E-01 -5.7550E-01 -1.1755E+00  2.1593E-01 -9.7980E-01  2.5893E-02  3.8357E-01 -7.1200E+00 -8.7425E-02 -4.4604E-01
             8.8539E-01
 GRADIENT:   3.1030E+00  6.7281E+00  1.0056E+01 -2.1391E+00 -2.0549E+01 -7.2111E-01  2.1664E+00  0.0000E+00 -4.0775E+00 -2.4908E+00
            -2.5141E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1784.56202482486        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1007
 NPARAMETR:  1.0460E+00  5.2107E-01  2.7592E-01  1.1170E+00  3.4157E-01  9.3145E-01  1.1498E+00  1.0000E-02  8.6043E-01  6.2618E-01
             2.1896E+00
 PARAMETER:  1.4501E-01 -5.5187E-01 -1.1876E+00  2.1062E-01 -9.7419E-01  2.8984E-02  2.3962E-01 -7.3547E+00 -5.0319E-02 -3.6811E-01
             8.8371E-01
 GRADIENT:  -6.3316E-01 -5.2847E-01 -8.5701E-01 -2.0026E+00  2.3437E+00 -2.8252E-02  2.0647E-01  0.0000E+00  3.7689E-01  1.5506E-01
             3.5748E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1784.85615344386        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1187
 NPARAMETR:  1.0474E+00  4.7787E-01  2.3869E-01  1.0936E+00  3.0384E-01  9.3233E-01  7.1457E-01  1.0000E-02  9.1066E-01  7.0942E-01
             2.1571E+00
 PARAMETER:  1.4634E-01 -6.3843E-01 -1.3326E+00  1.8947E-01 -1.0912E+00  2.9936E-02 -2.3608E-01 -9.7448E+00  6.4106E-03 -2.4330E-01
             8.6878E-01
 GRADIENT:   8.9480E-01 -8.3263E-01 -1.5585E+00 -3.2037E+00  2.2841E+00 -4.6591E-01  4.8532E-01  0.0000E+00 -1.0682E+00 -4.4746E-01
            -3.2280E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1784.99540817296        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1362
 NPARAMETR:  1.0478E+00  4.6745E-01  2.2683E-01  1.0845E+00  2.9285E-01  9.3362E-01  3.6129E-01  1.0000E-02  9.3266E-01  7.4464E-01
             2.1569E+00
 PARAMETER:  1.4672E-01 -6.6046E-01 -1.3835E+00  1.8110E-01 -1.1281E+00  3.1316E-02 -9.1809E-01 -1.1505E+01  3.0284E-02 -1.9486E-01
             8.6865E-01
 GRADIENT:   8.1640E-01  1.6549E+00 -5.1096E-01  1.9274E+00 -1.7569E+00 -7.7700E-02  5.8460E-02  0.0000E+00  3.7297E-01  6.0954E-01
             2.5197E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1785.00474822321        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1531
 NPARAMETR:  1.0478E+00  4.6481E-01  2.2642E-01  1.0834E+00  2.9202E-01  9.3368E-01  2.6942E-01  1.0000E-02  9.3211E-01  7.4728E-01
             2.1575E+00
 PARAMETER:  1.4640E-01 -6.6611E-01 -1.3854E+00  1.7999E-01 -1.1305E+00  3.1286E-02 -1.1995E+00 -1.2111E+01  2.9756E-02 -1.9090E-01
             8.6920E-01
 GRADIENT:  -4.3668E-01  7.4946E-03 -1.4076E-02 -9.8729E-02  6.0299E-01 -1.8228E-02  4.4930E-04  0.0000E+00  9.0717E-03  4.9827E-02
             1.1099E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1531
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.1185E-04 -3.5227E-03  5.3473E-05 -6.7301E-03 -2.1577E-03
 SE:             2.9391E-02  4.2821E-03  2.3118E-04  2.7277E-02  2.5898E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7525E-01  4.1070E-01  8.1708E-01  8.0511E-01  9.3360E-01

 ETASHRINKSD(%)  1.5359E+00  8.5654E+01  9.9226E+01  8.6199E+00  1.3239E+01
 ETASHRINKVR(%)  3.0481E+00  9.7942E+01  9.9994E+01  1.6497E+01  2.4726E+01
 EBVSHRINKSD(%)  1.6207E+00  8.5907E+01  9.9290E+01  7.6874E+00  1.3398E+01
 EBVSHRINKVR(%)  3.2150E+00  9.8014E+01  9.9995E+01  1.4784E+01  2.5001E+01
 RELATIVEINF(%)  9.6692E+01  2.8052E-01  4.4044E-04  4.9590E+01  3.2077E+00
 EPSSHRINKSD(%)  2.9183E+01
 EPSSHRINKVR(%)  4.9850E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1785.0047482232080     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -866.06621501853533     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.21
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1785.005       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  4.65E-01  2.26E-01  1.08E+00  2.92E-01  9.34E-01  2.73E-01  1.00E-02  9.32E-01  7.48E-01  2.16E+00
 


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
+        1.13E+03
 
 TH 2
+       -3.44E+01  1.38E+03
 
 TH 3
+       -1.22E+01  1.93E+03  1.06E+04
 
 TH 4
+       -2.33E+01  2.52E+02 -9.16E+02  9.33E+02
 
 TH 5
+        9.39E+01 -4.07E+03 -1.18E+04 -1.89E+02  1.76E+04
 
 TH 6
+        2.81E+00 -1.47E+01  2.94E+01 -1.04E+01  2.45E+00  2.13E+02
 
 TH 7
+        4.43E-01 -1.53E+00 -1.35E+01 -2.11E+00  3.79E+00  6.33E-02  9.63E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.33E+00 -4.44E+01  1.89E+02 -1.52E+01  9.37E+00 -4.72E-01  1.86E+00  0.00E+00  1.61E+02
 
 TH10
+       -5.79E-01 -2.87E+00 -3.30E+01  1.48E+01 -1.83E+01  3.80E+00  9.82E+00  0.00E+00  5.27E+00  2.04E+02
 
 TH11
+       -1.55E+01 -1.56E+01 -7.85E+01 -1.00E+01  3.96E+01  2.14E+00  2.56E+00  0.00E+00  1.01E+01  1.56E+01  9.53E+01
 
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
 #CPUT: Total CPU Time in Seconds,       28.833
Stop Time:
Wed Sep 29 23:07:47 CDT 2021
