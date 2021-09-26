Sat Sep 25 13:55:34 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat95.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m95.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1641.76956797303        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.8983E+01 -8.8646E+01 -4.7917E+01 -1.0810E+02  5.7724E+01  1.2530E+01 -1.8495E+01  1.7454E+01 -6.0143E+01  1.4724E+01
             4.4806E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1657.45192859224        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9761E-01  1.0072E+00  1.1657E+00  1.0603E+00  1.0121E+00  9.8425E-01  1.0552E+00  8.1567E-01  1.2556E+00  9.0310E-01
             9.9859E-01
 PARAMETER:  9.7604E-02  1.0713E-01  2.5330E-01  1.5856E-01  1.1204E-01  8.4122E-02  1.5370E-01 -1.0375E-01  3.2763E-01 -1.9185E-03
             9.8584E-02
 GRADIENT:   2.5963E+01  1.3226E+00  1.5382E+01 -8.3643E+00 -4.8571E+00  6.5846E+00  2.4908E+00  4.6794E+00  7.8667E+00 -1.0579E+01
             2.5387E-02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1660.43038035968        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.9388E-01  9.0594E-01  8.8140E-01  1.1269E+00  8.4846E-01  9.7933E-01  1.2246E+00  4.2253E-01  1.1561E+00  7.8465E-01
             9.7222E-01
 PARAMETER:  9.3862E-02  1.2218E-03 -2.6239E-02  2.1949E-01 -6.4332E-02  7.9110E-02  3.0261E-01 -7.6149E-01  2.4507E-01 -1.4251E-01
             7.1830E-02
 GRADIENT:   1.5072E+01 -3.3532E+00 -2.0371E+01  2.8073E+01  1.8367E+01  3.6986E+00 -3.7247E+00  3.3067E+00  9.0237E+00 -6.2880E-01
            -6.1935E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1661.36997976401        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.8714E-01  9.3186E-01  8.9429E-01  1.0912E+00  8.6308E-01  9.6816E-01  1.2544E+00  2.3750E-01  1.1232E+00  8.1662E-01
             9.8975E-01
 PARAMETER:  8.7061E-02  2.9430E-02 -1.1729E-02  1.8725E-01 -4.7247E-02  6.7641E-02  3.2666E-01 -1.3376E+00  2.1616E-01 -1.0258E-01
             8.9693E-02
 GRADIENT:  -2.6448E+00 -2.0503E+00  8.2086E-01 -6.8717E+00 -2.8432E+00 -1.0227E+00 -3.7820E-01  7.5359E-01 -1.5402E+00  1.0926E+00
             9.7615E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1662.40384904427        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  9.9946E-01  8.7938E-01  8.8687E-01  1.1343E+00  8.3749E-01  9.7509E-01  1.3389E+00  7.1998E-02  1.0958E+00  8.0243E-01
             9.9159E-01
 PARAMETER:  9.9456E-02 -2.8533E-02 -2.0051E-02  2.2600E-01 -7.7342E-02  7.4774E-02  3.9185E-01 -2.5311E+00  1.9148E-01 -1.2011E-01
             9.1552E-02
 GRADIENT:  -1.2776E+01  2.1813E+00  2.2699E+00 -5.8061E+00 -3.2691E+00 -1.0435E+00 -6.1500E-01  4.8420E-02 -1.3829E+00  4.8422E-02
             9.2361E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1663.06909138129        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      523
 NPARAMETR:  1.0019E+00  6.2728E-01  8.9544E-01  1.2867E+00  7.5245E-01  9.7398E-01  1.7458E+00  2.7620E-02  9.9222E-01  7.8085E-01
             9.9153E-01
 PARAMETER:  1.0186E-01 -3.6635E-01 -1.0438E-02  3.5211E-01 -1.8442E-01  7.3633E-02  6.5721E-01 -3.4892E+00  9.2185E-02 -1.4737E-01
             9.1493E-02
 GRADIENT:  -2.1734E+00  1.8504E+00  3.0787E+00  1.3073E+00 -4.2796E+00 -9.0133E-02 -9.9254E-02  2.4620E-03 -6.1102E-03 -7.5775E-01
            -2.6928E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1663.09244520803        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      698
 NPARAMETR:  1.0023E+00  5.8404E-01  9.0102E-01  1.3107E+00  7.4378E-01  9.7346E-01  1.8331E+00  2.1124E-02  9.7813E-01  7.8923E-01
             9.9215E-01
 PARAMETER:  1.0225E-01 -4.3778E-01 -4.2276E-03  3.7055E-01 -1.9601E-01  7.3106E-02  7.0602E-01 -3.7574E+00  7.7892E-02 -1.3669E-01
             9.2118E-02
 GRADIENT:   2.3029E-02  6.1397E-03  1.0456E-02  1.3525E-02 -1.2794E-02  1.3033E-02  6.3148E-03  1.0444E-03  2.1245E-02  2.4895E-03
            -6.1945E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1663.09258023670        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  1.0022E+00  5.8232E-01  9.0058E-01  1.3115E+00  7.4310E-01  9.7326E-01  1.8366E+00  1.5756E-02  9.7721E-01  7.8908E-01
             9.9225E-01
 PARAMETER:  1.0218E-01 -4.4073E-01 -4.7167E-03  3.7118E-01 -1.9692E-01  7.2899E-02  7.0791E-01 -4.0505E+00  7.6942E-02 -1.3689E-01
             9.2218E-02
 GRADIENT:  -1.1689E-01 -8.7554E-02 -2.0463E-01 -7.7401E-02  2.0631E-01 -5.7330E-02 -1.8876E-02  5.9239E-04 -6.9287E-02  2.1669E-02
             3.9821E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1663.09286022915        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  1.0023E+00  5.8463E-01  9.0073E-01  1.3102E+00  7.4387E-01  9.7339E-01  1.8318E+00  1.0000E-02  9.7818E-01  7.8917E-01
             9.9219E-01
 PARAMETER:  1.0226E-01 -4.3666E-01 -4.4917E-03  3.7025E-01 -1.9591E-01  7.3110E-02  7.0525E-01 -4.6904E+00  7.8042E-02 -1.3683E-01
             9.2131E-02
 GRADIENT:   9.2715E-04  1.0300E-02  2.4981E-02  5.4561E-02 -2.1130E-02  1.0017E-02 -1.3810E-03  0.0000E+00  8.9974E-03 -3.3375E-03
            -4.1258E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1013
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3526E-04  1.6033E-02 -4.7654E-04 -1.2160E-02 -4.3733E-03
 SE:             2.9839E-02  1.8648E-02  2.2562E-04  2.6345E-02  2.2975E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9638E-01  3.8993E-01  3.4672E-02  6.4439E-01  8.4903E-01

 ETASHRINKSD(%)  3.5013E-02  3.7527E+01  9.9244E+01  1.1740E+01  2.3032E+01
 ETASHRINKVR(%)  7.0013E-02  6.0971E+01  9.9994E+01  2.2101E+01  4.0759E+01
 EBVSHRINKSD(%)  4.4356E-01  4.0453E+01  9.9265E+01  1.0775E+01  2.0461E+01
 EBVSHRINKVR(%)  8.8515E-01  6.4541E+01  9.9995E+01  2.0389E+01  3.6736E+01
 RELATIVEINF(%)  9.8529E+01  3.9917E+00  7.3497E-04  1.4038E+01  5.3189E+00
 EPSSHRINKSD(%)  4.3478E+01
 EPSSHRINKVR(%)  6.8053E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1663.0928602291463     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -927.94203366540808     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.31
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1663.093       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  5.85E-01  9.01E-01  1.31E+00  7.44E-01  9.73E-01  1.83E+00  1.00E-02  9.78E-01  7.89E-01  9.92E-01
 


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
+        1.16E+03
 
 TH 2
+       -1.19E+01  3.70E+02
 
 TH 3
+        2.13E+01  2.67E+02  7.82E+02
 
 TH 4
+       -4.67E+00  2.36E+02 -1.63E+02  5.45E+02
 
 TH 5
+       -6.27E+00 -5.22E+02 -1.13E+03  1.93E+02  2.03E+03
 
 TH 6
+       -3.27E-01 -2.08E+00  4.77E+00 -1.05E+00 -2.72E+00  2.06E+02
 
 TH 7
+        7.55E-01  2.91E+01  2.58E+00 -3.91E+00 -7.08E+00  4.99E-01  1.46E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.36E+00 -2.12E+01 -1.62E+01  1.52E+01  4.25E-01 -2.58E+00  7.15E+00  0.00E+00  1.36E+02
 
 TH10
+        2.17E+00  6.78E+00 -8.20E+01 -2.70E+01 -2.19E+01  9.92E-02  7.91E+00  0.00E+00  9.01E+00  1.25E+02
 
 TH11
+       -5.62E+00 -5.95E+00 -4.67E+01 -9.07E+00  2.22E+01 -1.51E-01  2.16E+00  0.00E+00  1.11E+01  3.05E+01  2.23E+02
 
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
 #CPUT: Total CPU Time in Seconds,       17.354
Stop Time:
Sat Sep 25 13:55:54 CDT 2021
