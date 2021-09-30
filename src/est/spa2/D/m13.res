Thu Sep 30 08:36:35 CDT 2021
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
$DATA ../../../../data/spa2/D/dat13.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2018.77957381144        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3414E+02 -1.3180E+02 -4.7743E+01 -1.7609E+02  1.3374E+02 -2.0011E+02 -2.6205E+02 -3.8315E+01 -3.5946E+02 -7.7377E+01
            -4.1452E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2221.98810502395        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  9.2876E-01  1.1075E+00  1.0778E+00  9.9248E-01  1.0871E+00  1.0705E+00  1.7390E+00  1.2051E+00  1.0820E+00  1.0370E+00
             9.7892E-01
 PARAMETER:  2.6096E-02  2.0206E-01  1.7489E-01  9.2451E-02  1.8356E-01  1.6809E-01  6.5329E-01  2.8653E-01  1.7883E-01  1.3630E-01
             7.8697E-02
 GRADIENT:   1.3865E+02  4.4162E+01 -1.1706E+01  7.4050E+01 -1.9055E+01 -2.8737E+01  7.5590E+01 -9.7598E+00 -2.8614E+01  1.9742E+01
            -3.4575E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2261.80702394695        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      250
 NPARAMETR:  1.0849E+00  1.5313E+00  1.3000E+00  8.2724E-01  1.5812E+00  1.1858E+00  1.5947E+00  1.6842E+00  1.4147E+00  8.1710E-01
             1.0593E+00
 PARAMETER:  1.8151E-01  5.2612E-01  3.6234E-01 -8.9659E-02  5.5821E-01  2.7039E-01  5.6666E-01  6.2130E-01  4.4690E-01 -1.0199E-01
             1.5760E-01
 GRADIENT:  -1.3451E+00 -6.6791E+01 -3.3614E+01  7.9420E+01  2.8786E+02 -1.6339E+02 -1.5150E+02 -1.7127E+01 -1.4468E+01 -9.0441E+01
             2.3077E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2336.13319114327        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  1.0814E+00  1.5954E+00  1.1266E+00  7.3266E-01  1.4475E+00  1.4549E+00  2.0651E+00  1.6432E+00  1.1238E+00  1.0713E+00
             1.0170E+00
 PARAMETER:  1.7821E-01  5.6710E-01  2.1918E-01 -2.1107E-01  4.6984E-01  4.7491E-01  8.2519E-01  5.9666E-01  2.1671E-01  1.6885E-01
             1.1681E-01
 GRADIENT:  -1.3493E+00 -2.2372E+01 -4.3195E+00 -8.0288E+00  4.0212E+01 -3.2159E+01 -4.1786E+01 -3.4190E+00 -8.5914E+00 -1.7774E+01
             5.4352E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2344.91397790560        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      604
 NPARAMETR:  1.0762E+00  1.4235E+00  1.4005E+00  8.8529E-01  1.3288E+00  1.6033E+00  2.5570E+00  1.7725E+00  1.0126E+00  1.0892E+00
             1.0107E+00
 PARAMETER:  1.7341E-01  4.5310E-01  4.3685E-01 -2.1835E-02  3.8429E-01  5.7208E-01  1.0388E+00  6.7239E-01  1.1249E-01  1.8540E-01
             1.1065E-01
 GRADIENT:  -3.7517E+00  6.1973E+00  2.2702E+00 -7.6531E+00 -1.9707E+01  1.6555E+01 -4.6662E+00  3.2043E+00 -5.4090E+00 -2.2453E-01
             1.1579E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2345.24625227156        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      780
 NPARAMETR:  1.0833E+00  1.3734E+00  1.4601E+00  9.1500E-01  1.3347E+00  1.5575E+00  2.6232E+00  1.7815E+00  1.0221E+00  1.0852E+00
             1.0068E+00
 PARAMETER:  1.8001E-01  4.1725E-01  4.7850E-01  1.1168E-02  3.8872E-01  5.4305E-01  1.0644E+00  6.7745E-01  1.2188E-01  1.8173E-01
             1.0679E-01
 GRADIENT:   2.2561E+00  4.1786E+00  1.5105E+00  8.4456E-01 -3.4995E+00  2.4618E+00 -1.9108E+00 -9.5444E-01 -3.6844E+00  4.5324E-01
             7.3669E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2345.34658707726        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      964
 NPARAMETR:  1.0826E+00  1.3606E+00  1.4681E+00  9.2033E-01  1.3331E+00  1.5579E+00  2.6343E+00  1.7887E+00  1.0376E+00  1.0770E+00
             1.0062E+00
 PARAMETER:  1.7937E-01  4.0789E-01  4.8397E-01  1.6977E-02  3.8750E-01  5.4335E-01  1.0686E+00  6.8148E-01  1.3688E-01  1.7416E-01
             1.0619E-01
 GRADIENT:   1.6610E+00  2.9851E+00  5.7524E-01  9.7151E-01 -1.6033E+00  2.6095E+00 -1.0232E+00 -3.1079E-01 -2.4314E+00  1.1158E-01
             2.0584E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2345.59470211009        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  1.0844E+00  1.3596E+00  1.4681E+00  9.2044E-01  1.3337E+00  1.5871E+00  2.6482E+00  1.7875E+00  1.0724E+00  1.0767E+00
             1.0063E+00
 PARAMETER:  1.8107E-01  4.0721E-01  4.8398E-01  1.7101E-02  3.8798E-01  5.6193E-01  1.0739E+00  6.8081E-01  1.6990E-01  1.7387E-01
             1.0624E-01
 GRADIENT:   3.3959E+00  2.4438E+00 -3.4682E-01  3.8994E-01  2.2636E-01  1.1410E+01  2.0343E+00  4.7097E-01  5.2695E-01  9.1099E-01
             3.9997E-01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2345.59470211009        NO. OF FUNC. EVALS.:  25
 CUMULATIVE NO. OF FUNC. EVALS.:     1175
 NPARAMETR:  1.0844E+00  1.3596E+00  1.4681E+00  9.2044E-01  1.3337E+00  1.5871E+00  2.6482E+00  1.7875E+00  1.0724E+00  1.0767E+00
             1.0063E+00
 PARAMETER:  1.8107E-01  4.0721E-01  4.8398E-01  1.7101E-02  3.8798E-01  5.6193E-01  1.0739E+00  6.8081E-01  1.6990E-01  1.7387E-01
             1.0624E-01
 GRADIENT:  -1.7647E+00  4.4384E+03  7.4784E+03 -4.3444E+00 -9.3129E+03 -3.2203E+03  3.3611E+03  4.6991E-01 -4.6094E+00 -2.3090E+00
            -3.4007E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1175
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7802E-03  3.1549E-04 -4.4655E-02 -5.3959E-03 -3.0121E-02
 SE:             3.0426E-02  2.7074E-02  1.7263E-02  2.0278E-02  2.1791E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5334E-01  9.9070E-01  9.6889E-03  7.9016E-01  1.6689E-01

 ETASHRINKSD(%)  1.0000E-10  9.2994E+00  4.2166E+01  3.2067E+01  2.6996E+01
 ETASHRINKVR(%)  1.0000E-10  1.7734E+01  6.6553E+01  5.3851E+01  4.6705E+01
 EBVSHRINKSD(%)  1.3164E-01  7.7981E+00  4.6856E+01  3.7248E+01  2.2866E+01
 EBVSHRINKVR(%)  2.6311E-01  1.4988E+01  7.1758E+01  6.0622E+01  4.0504E+01
 RELATIVEINF(%)  9.9729E+01  4.5093E+01  1.5434E+01  1.4842E+01  3.1456E+01
 EPSSHRINKSD(%)  3.1003E+01
 EPSSHRINKVR(%)  5.2394E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2345.5947021100901     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1242.8684622644830     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.49
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2345.595       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.08E+00  1.36E+00  1.47E+00  9.20E-01  1.33E+00  1.59E+00  2.65E+00  1.79E+00  1.07E+00  1.08E+00  1.01E+00
 


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
+        3.74E+02
 
 TH 2
+       -8.31E+05  2.95E+05
 
 TH 3
+        6.46E+05  2.29E+05  1.78E+05
 
 TH 4
+       -6.67E+02  1.77E+06 -1.38E+06  1.07E+07
 
 TH 5
+       -9.84E+00 -6.24E+01  4.91E+05  4.39E+02  6.75E+05
 
 TH 6
+        1.05E+00  1.83E+05 -1.43E+05  1.72E-01 -1.63E+01  2.27E+05
 
 TH 7
+        1.65E+00  5.58E+00 -8.94E+04 -6.42E+01  8.47E+00  2.92E+00  2.24E+04
 
 TH 8
+        4.37E-01 -1.36E+05 -2.56E+03  2.15E+00  3.44E+03 -1.66E+05  2.58E+04  5.96E+04
 
 TH 9
+        8.10E-02  8.95E+05 -6.96E+05  7.14E+02  2.36E+02 -7.49E-02  1.74E+05  1.56E+01  3.68E+01
 
 TH10
+       -3.80E+02 -8.72E+05  6.78E+05 -1.05E+07  1.10E+02 -5.44E-02 -2.64E+01 -3.98E+00  4.20E+02  5.15E+06
 
 TH11
+       -2.89E+00 -1.54E+02  4.53E+03  2.88E+01  1.67E+02  9.47E+05 -2.45E+01  7.03E+05 -6.12E-01  2.56E+01  7.90E+06
 
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
 #CPUT: Total CPU Time in Seconds,       34.704
Stop Time:
Thu Sep 30 08:37:11 CDT 2021
