Wed Sep 29 11:29:43 CDT 2021
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
$DATA ../../../../data/spa/B/dat73.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m73.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1714.65756071673        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9637E+02  4.2516E+00  8.9506E+00  6.0702E+00 -3.7924E+01  7.9202E+01 -2.1475E+00  3.1657E+00 -1.6166E+00  2.6118E+01
             2.6796E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1722.94188849096        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0044E+00  1.0698E+00  1.0385E+00  1.0151E+00  1.0542E+00  8.6185E-01  1.0451E+00  1.0016E+00  1.0350E+00  8.8174E-01
             9.4485E-01
 PARAMETER:  1.0438E-01  1.6743E-01  1.3776E-01  1.1502E-01  1.5277E-01 -4.8674E-02  1.4413E-01  1.0157E-01  1.3440E-01 -2.5858E-02
             4.3269E-02
 GRADIENT:  -2.8598E+01  2.0316E+01  1.1753E+01  2.4212E+01 -1.5943E+01 -6.3184E+00 -3.7971E+00 -5.0960E+00 -3.7344E-01  6.0108E+00
            -1.5350E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1724.27453744353        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0135E+00  1.1309E+00  9.5679E-01  9.6840E-01  1.0365E+00  8.7311E-01  1.0981E+00  1.0829E+00  1.0513E+00  7.3853E-01
             9.5714E-01
 PARAMETER:  1.1337E-01  2.2305E-01  5.5825E-02  6.7887E-02  1.3587E-01 -3.5699E-02  1.9354E-01  1.7961E-01  1.5002E-01 -2.0309E-01
             5.6191E-02
 GRADIENT:  -3.8688E+00  1.7768E+01  7.3289E+00  1.9080E+01 -1.0724E+01 -1.1983E+00  5.2224E-01 -5.9794E-01  2.0488E+00 -2.0461E+00
             2.8378E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1725.95894468426        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  1.0181E+00  1.4503E+00  6.2806E-01  7.4049E-01  1.0446E+00  8.7986E-01  9.0746E-01  6.3876E-01  1.2026E+00  7.3308E-01
             9.4094E-01
 PARAMETER:  1.1795E-01  4.7175E-01 -3.6512E-01 -2.0045E-01  1.4364E-01 -2.7998E-02  2.8940E-03 -3.4823E-01  2.8448E-01 -2.1050E-01
             3.9123E-02
 GRADIENT:   2.6342E+00  1.2875E+01  4.5309E+00  8.8482E+00 -5.7358E+00  7.4170E-01 -5.0919E+00 -2.5249E-01 -5.9358E+00 -8.1466E-01
            -3.9972E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1726.93471706074        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0147E+00  1.7576E+00  4.3931E-01  5.4605E-01  1.1260E+00  8.8046E-01  8.0645E-01  3.1718E-01  1.4952E+00  7.3578E-01
             9.5291E-01
 PARAMETER:  1.1464E-01  6.6396E-01 -7.2254E-01 -5.0504E-01  2.1869E-01 -2.7311E-02 -1.1512E-01 -1.0483E+00  5.0224E-01 -2.0682E-01
             5.1763E-02
 GRADIENT:  -8.4352E+00  2.9630E+01 -7.1340E-01  1.7366E+01 -3.0161E+00  7.7658E-01 -2.2241E+00  1.6830E-01 -2.5429E+00 -2.7501E+00
             7.0317E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1727.41400961678        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.0175E+00  1.9269E+00  3.8528E-01  4.3191E-01  1.2294E+00  8.7669E-01  7.6047E-01  1.8794E-01  1.8150E+00  8.2661E-01
             9.5548E-01
 PARAMETER:  1.1736E-01  7.5591E-01 -8.5379E-01 -7.3953E-01  3.0652E-01 -3.1598E-02 -1.7382E-01 -1.5716E+00  6.9608E-01 -9.0417E-02
             5.4455E-02
 GRADIENT:  -4.7662E-02  2.1168E+01  6.7330E-01  1.0876E+01 -6.7036E+00 -6.9025E-01  2.7742E-01  4.9017E-02  2.4300E+00 -1.4142E-01
             5.8131E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1727.73573127311        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1089
 NPARAMETR:  1.0183E+00  1.9097E+00  3.8312E-01  4.2356E-01  1.2340E+00  8.7829E-01  7.5948E-01  4.0422E-02  1.8063E+00  8.3000E-01
             9.5329E-01
 PARAMETER:  1.1811E-01  7.4692E-01 -8.5942E-01 -7.5906E-01  3.1027E-01 -2.9784E-02 -1.7512E-01 -3.1084E+00  6.9126E-01 -8.6334E-02
             5.2162E-02
 GRADIENT:   2.5767E+00 -1.6580E+01 -8.7542E-01  5.6173E-01  1.8355E+00  8.3364E-02  1.5693E-01  2.7577E-03  7.3524E-01  1.5334E-01
             1.3435E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1727.73918689112        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1269
 NPARAMETR:  1.0183E+00  1.9099E+00  3.8370E-01  4.2342E-01  1.2334E+00  8.7829E-01  7.5938E-01  1.0000E-02  1.8040E+00  8.2876E-01
             9.5313E-01
 PARAMETER:  1.1810E-01  7.4704E-01 -8.5790E-01 -7.5939E-01  3.0976E-01 -2.9776E-02 -1.7526E-01 -5.2328E+00  6.9000E-01 -8.7826E-02
             5.1997E-02
 GRADIENT:   2.5653E+00 -1.5191E+01 -2.3849E-01 -5.0933E-02  6.0903E-01  8.8106E-02  2.2128E-02  0.0000E+00  3.5439E-01 -7.4141E-03
            -1.1874E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1727.73918689112        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1297
 NPARAMETR:  1.0183E+00  1.9074E+00  3.8700E-01  4.2265E-01  1.2297E+00  8.7839E-01  7.5978E-01  1.0000E-02  1.7934E+00  8.2914E-01
             9.5330E-01
 PARAMETER:  1.1810E-01  7.4704E-01 -8.5790E-01 -7.5939E-01  3.0976E-01 -2.9776E-02 -1.7526E-01 -5.2328E+00  6.9000E-01 -8.7826E-02
             5.1997E-02
 GRADIENT:  -2.1456E-02  5.7226E-01 -2.2475E-01  1.1352E-01  6.9281E-01 -6.5297E-03 -1.8756E-02  0.0000E+00  1.6732E-01 -7.8812E-03
            -1.1479E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1297
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.0735E-04 -3.3562E-02 -2.3768E-04  3.5679E-02 -4.5368E-02
 SE:             2.9845E-02  2.5528E-02  8.8348E-05  2.2877E-02  2.0186E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8911E-01  1.8860E-01  7.1382E-03  1.1886E-01  2.4609E-02

 ETASHRINKSD(%)  1.4445E-02  1.4480E+01  9.9704E+01  2.3359E+01  3.2374E+01
 ETASHRINKVR(%)  2.8888E-02  2.6862E+01  9.9999E+01  4.1261E+01  5.4267E+01
 EBVSHRINKSD(%)  5.0119E-01  1.4578E+01  9.9737E+01  2.3592E+01  3.1837E+01
 EBVSHRINKVR(%)  9.9986E-01  2.7031E+01  9.9999E+01  4.1619E+01  5.3538E+01
 RELATIVEINF(%)  9.8969E+01  8.1270E+00  8.7596E-05  6.1286E+00  1.2214E+01
 EPSSHRINKSD(%)  4.4843E+01
 EPSSHRINKVR(%)  6.9577E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1727.7391868911238     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -992.58836032738566     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1727.739       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.91E+00  3.84E-01  4.23E-01  1.23E+00  8.78E-01  7.59E-01  1.00E-02  1.80E+00  8.29E-01  9.53E-01
 


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
+        1.38E+03
 
 TH 2
+       -6.67E+00  4.04E+02
 
 TH 3
+        4.73E+00  2.12E+02  5.91E+02
 
 TH 4
+       -1.72E+01  2.52E+02 -4.66E+02  1.16E+03
 
 TH 5
+       -6.37E+00 -1.80E+02 -3.60E+02  3.15E+02  5.11E+02
 
 TH 6
+        3.16E-02 -1.06E+00  1.64E+00 -4.02E+00 -1.48E+00  2.53E+02
 
 TH 7
+        3.17E-01  9.42E+00 -4.17E+01 -1.32E+01 -8.52E-01 -7.61E-01  2.04E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.31E+00 -1.81E+01 -3.55E+01  6.93E+01 -9.91E+00  7.17E-03  6.35E+00  0.00E+00  2.91E+01
 
 TH10
+        6.83E-02 -1.55E+01 -3.19E+01 -1.28E+01 -7.68E+01  5.63E-01  2.14E+01  0.00E+00  6.44E+00  8.26E+01
 
 TH11
+       -8.33E+00 -1.71E+01 -1.97E+01 -7.87E-01 -1.40E+01  3.02E+00  1.14E+01  0.00E+00  3.60E+00  1.89E+01  2.32E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.215
Stop Time:
Wed Sep 29 11:30:09 CDT 2021
