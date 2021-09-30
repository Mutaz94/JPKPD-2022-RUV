Wed Sep 29 13:10:38 CDT 2021
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
$DATA ../../../../data/spa/A2/dat95.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1099.38180815477        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0994E+02 -4.0009E+00 -1.4540E+01 -2.5122E+01  1.6842E+02  4.7090E+01 -5.1190E+01  4.9252E+00 -1.3229E+02 -5.0790E+01
            -8.7576E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1359.89700332445        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.2599E+00  8.7200E-01  1.0106E+00  1.1304E+00  8.6186E-01  1.2879E+00  1.0799E+00  9.0828E-01  1.4978E+00  8.5810E-01
             2.3467E+00
 PARAMETER:  3.3107E-01 -3.6971E-02  1.1052E-01  2.2255E-01 -4.8668E-02  3.5303E-01  1.7691E-01  3.7972E-03  5.0401E-01 -5.3037E-02
             9.5301E-01
 GRADIENT:   5.3774E+02 -9.9054E+00 -2.0897E+00  9.5523E+00  2.1099E+01  4.9140E+01  8.9023E+00  7.5053E+00  3.7114E+01  8.4143E+00
             2.6796E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1376.33787381059        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.2245E+00  5.4573E-01  4.5685E-01  1.3140E+00  4.4321E-01  1.1995E+00  9.2369E-01  4.2235E-01  1.1436E+00  2.2719E-01
             2.4849E+00
 PARAMETER:  3.0252E-01 -5.0563E-01 -6.8339E-01  3.7307E-01 -7.1372E-01  2.8193E-01  2.0617E-02 -7.6191E-01  2.3415E-01 -1.3820E+00
             1.0102E+00
 GRADIENT:   4.5657E+02  2.9339E+01  1.2622E+00  1.2541E+02  2.0734E+01  3.7337E+01 -1.4975E+01 -3.9076E+00 -3.1412E+00 -5.1272E+00
             1.0974E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1415.56057166157        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      344
 NPARAMETR:  1.0474E+00  6.4295E-01  4.5886E-01  1.2127E+00  4.6594E-01  8.9020E-01  1.7641E+00  2.3547E-01  1.3281E+00  4.8392E-01
             1.9093E+00
 PARAMETER:  1.4634E-01 -3.4170E-01 -6.7901E-01  2.9281E-01 -6.6370E-01 -1.6313E-02  6.6763E-01 -1.3462E+00  3.8376E-01 -6.2583E-01
             7.4674E-01
 GRADIENT:   9.8029E+01  3.8110E+01  8.4205E+00  1.6035E+01 -1.1391E+01 -2.6761E+01  2.0839E+01  2.7539E-01  4.1978E+01 -3.8858E-01
            -5.5957E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1433.58131974866        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  1.0013E+00  3.7670E-01  3.8500E-01  1.2817E+00  3.6568E-01  9.4216E-01  1.8031E+00  2.9558E-02  1.1425E+00  5.4689E-01
             2.0450E+00
 PARAMETER:  1.0128E-01 -8.7630E-01 -8.5451E-01  3.4817E-01 -9.0598E-01  4.0417E-02  6.8950E-01 -3.4214E+00  2.3319E-01 -5.0351E-01
             8.1540E-01
 GRADIENT:  -2.7793E+01 -7.2720E-01 -2.4487E+01  8.3773E+00  3.9510E+01  3.9529E-02 -6.3234E+00  7.7000E-03  1.3268E+01  2.4782E+00
             1.7528E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1439.46488951940        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      697
 NPARAMETR:  1.0121E+00  2.0304E-01  4.2327E-01  1.3644E+00  3.5445E-01  9.2260E-01  2.8982E+00  1.0000E-02  1.0555E+00  5.8898E-01
             2.0873E+00
 PARAMETER:  1.1207E-01 -1.4944E+00 -7.5975E-01  4.1074E-01 -9.3718E-01  1.9436E-02  1.1641E+00 -5.7479E+00  1.5400E-01 -4.2937E-01
             8.3587E-01
 GRADIENT:   1.2055E+01  9.9046E+00  1.1662E+01 -1.2403E+01 -1.5695E+01 -4.5974E+00  6.9393E+00  0.0000E+00 -1.1680E+00 -2.3817E+00
             6.4283E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1444.02838389787        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      872
 NPARAMETR:  9.9698E-01  1.0850E-01  4.2272E-01  1.4245E+00  3.4723E-01  9.4484E-01  4.1733E+00  1.0000E-02  1.0086E+00  5.8940E-01
             2.0284E+00
 PARAMETER:  9.6975E-02 -2.1210E+00 -7.6104E-01  4.5383E-01 -9.5777E-01  4.3256E-02  1.5287E+00 -6.8008E+00  1.0860E-01 -4.2864E-01
             8.0726E-01
 GRADIENT:  -1.2330E+01  7.2516E+00 -1.0149E+01  1.7941E+01  1.0747E+01  4.4791E+00  8.2787E+00  0.0000E+00 -6.6061E+00 -7.8326E+00
            -8.5582E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1448.83402620714        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1050
 NPARAMETR:  1.0010E+00  6.7449E-02  3.7109E-01  1.3864E+00  3.1287E-01  9.2899E-01  5.1965E+00  1.0000E-02  1.0254E+00  6.0537E-01
             2.0154E+00
 PARAMETER:  1.0101E-01 -2.5964E+00 -8.9130E-01  4.2669E-01 -1.0620E+00  2.6346E-02  1.7480E+00 -7.5187E+00  1.2507E-01 -4.0192E-01
             8.0082E-01
 GRADIENT:   2.7835E+00  1.3296E+01 -1.2139E+01 -6.5172E+00  1.3223E+01 -2.8471E+00  1.8705E+01  0.0000E+00 -7.7070E+00 -6.3123E+00
            -4.9843E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1450.06228584114        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1228
 NPARAMETR:  9.9958E-01  4.5232E-02  3.4997E-01  1.3703E+00  2.9956E-01  9.3107E-01  6.2157E+00  1.0000E-02  1.0471E+00  5.8878E-01
             2.0111E+00
 PARAMETER:  9.9578E-02 -2.9959E+00 -9.4990E-01  4.1501E-01 -1.1054E+00  2.8576E-02  1.9271E+00 -7.8274E+00  1.4607E-01 -4.2970E-01
             7.9870E-01
 GRADIENT:   4.0533E+00  1.5830E+01 -2.4934E+01 -1.8393E+01  3.4576E+01 -2.3796E+00  3.0325E+01  0.0000E+00 -6.4607E-01 -1.5017E+01
            -9.4264E+00

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1450.06228584114        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1257
 NPARAMETR:  9.9938E-01  4.5369E-02  3.5054E-01  1.3714E+00  2.9912E-01  9.3200E-01  6.2282E+00  1.0000E-02  1.0481E+00  5.8928E-01
             2.0096E+00
 PARAMETER:  9.9578E-02 -2.9959E+00 -9.4990E-01  4.1501E-01 -1.1054E+00  2.8576E-02  1.9271E+00 -7.8274E+00  1.4607E-01 -4.2970E-01
             7.9870E-01
 GRADIENT:   4.1264E+02 -5.4524E+01 -8.1286E+01 -1.9329E+02  6.1956E+01 -2.3193E+00 -8.8063E+01  0.0000E+00 -5.3418E-01 -1.8940E+02
             2.0047E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1257
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0041E-03  1.8712E-02 -8.6493E-06 -1.1933E-02  3.9480E-03
 SE:             2.9367E-02  1.0075E-02  2.7309E-04  2.8110E-02  2.1816E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7273E-01  6.3280E-02  9.7473E-01  6.7119E-01  8.5639E-01

 ETASHRINKSD(%)  1.6181E+00  6.6246E+01  9.9085E+01  5.8266E+00  2.6914E+01
 ETASHRINKVR(%)  3.2100E+00  8.8607E+01  9.9992E+01  1.1314E+01  4.6584E+01
 EBVSHRINKSD(%)  1.7756E+00  7.7442E+01  9.9059E+01  5.2607E+00  2.4382E+01
 EBVSHRINKVR(%)  3.5196E+00  9.4912E+01  9.9991E+01  1.0245E+01  4.2819E+01
 RELATIVEINF(%)  9.5231E+01  3.2241E+00  3.1445E-04  3.4780E+01  2.0312E+00
 EPSSHRINKSD(%)  3.7395E+01
 EPSSHRINKVR(%)  6.0806E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1450.0622858411396     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -714.91145927740138     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.51
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1450.062       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  4.52E-02  3.50E-01  1.37E+00  3.00E-01  9.31E-01  6.22E+00  1.00E-02  1.05E+00  5.89E-01  2.01E+00
 


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
+        2.07E+05
 
 TH 2
+       -1.40E+02  6.65E+05
 
 TH 3
+       -1.12E+01 -3.70E+04  8.51E+04
 
 TH 4
+        1.84E+01  6.50E+02 -3.35E+04  2.54E+04
 
 TH 5
+        1.44E+02  2.76E+04  2.26E+04 -3.05E+04  1.24E+05
 
 TH 6
+        3.16E+02 -6.38E+01 -9.58E-01 -1.90E+01 -4.64E+00  2.37E+05
 
 TH 7
+        3.37E+00  3.56E+03 -4.44E+02 -2.97E+02  1.36E+03 -1.89E+00  9.76E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.17E+02  6.00E+02  5.39E+01  6.96E+04 -1.13E+05  1.44E+05  7.00E+00  0.00E+00  2.26E+02
 
 TH10
+        8.40E+01  2.49E+03 -1.40E+02  4.16E+04 -6.79E+04 -3.30E+01 -6.53E+02  0.00E+00 -1.21E+01  3.19E+04
 
 TH11
+        6.00E+00  6.75E+02 -7.79E+02 -1.06E+02  1.31E+03  1.88E+00 -9.94E+01  0.00E+00 -5.53E+01 -2.70E+02  3.14E+03
 
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
 #CPUT: Total CPU Time in Seconds,       23.434
Stop Time:
Wed Sep 29 13:11:03 CDT 2021
