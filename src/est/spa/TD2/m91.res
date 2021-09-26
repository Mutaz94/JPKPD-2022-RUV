Sat Sep 25 13:53:16 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat91.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1606.51434016766        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8509E+02 -8.2588E+01 -5.7665E+01 -6.5949E+01  4.5303E+01 -1.8966E+01 -4.7061E+00  1.8640E+01 -8.9884E+00  2.7230E+01
            -4.7870E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1616.57304467527        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9269E-01  1.0548E+00  1.2766E+00  9.8531E-01  1.1173E+00  1.0886E+00  1.0083E+00  8.4437E-01  1.0777E+00  7.9721E-01
             1.2118E+00
 PARAMETER:  9.2667E-02  1.5339E-01  3.4418E-01  8.5202E-02  2.1087E-01  1.8485E-01  1.0831E-01 -6.9160E-02  1.7480E-01 -1.2663E-01
             2.9207E-01
 GRADIENT:   1.3584E+02 -5.0668E+01 -7.4415E+00 -5.2660E+01  6.3514E+01  1.8774E+01  2.7364E+00  1.8801E+00  2.1877E+00 -1.8040E+01
             2.4859E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1621.31799189577        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.8153E-01  8.9205E-01  1.1419E+00  1.1253E+00  1.0050E+00  1.0836E+00  7.9449E-01  2.3501E-01  9.8895E-01  9.7669E-01
             1.1642E+00
 PARAMETER:  8.1360E-02 -1.4232E-02  2.3274E-01  2.1804E-01  1.0501E-01  1.8033E-01 -1.3006E-01 -1.3481E+00  8.8893E-02  7.6411E-02
             2.5200E-01
 GRADIENT:   1.1728E+02 -1.8621E+01 -3.2196E+01  6.6940E+00  5.4152E+01  1.9797E+01 -4.5039E+00 -2.1048E-03 -1.3928E+01  5.4134E+00
             1.8795E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1627.03090000120        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.2859E-01  1.0297E+00  9.1238E-01  1.0087E+00  9.3327E-01  1.0212E+00  1.1121E+00  2.4548E-01  1.0373E+00  7.9899E-01
             1.0992E+00
 PARAMETER:  2.5909E-02  1.2930E-01  8.3059E-03  1.0867E-01  3.0943E-02  1.2097E-01  2.0625E-01 -1.3045E+00  1.3662E-01 -1.2440E-01
             1.9459E-01
 GRADIENT:   1.2763E+01 -2.1485E+01 -1.7406E+01 -1.7341E+01  2.1981E+01 -3.8372E-01  3.0095E+00  8.8073E-01  6.5750E+00  6.0958E+00
            -1.3337E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1628.50476120900        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.2400E-01  8.6172E-01  9.2678E-01  1.1336E+00  8.5264E-01  1.0280E+00  1.2693E+00  1.5513E-01  9.1617E-01  7.1680E-01
             1.1093E+00
 PARAMETER:  2.0953E-02 -4.8820E-02  2.3960E-02  2.2543E-01 -5.9418E-02  1.2765E-01  3.3843E-01 -1.7635E+00  1.2450E-02 -2.3295E-01
             2.0369E-01
 GRADIENT:   2.9766E+00  5.3645E+00  2.3496E+00  1.6267E+01 -1.6689E-01  2.5107E+00 -2.6046E+00  2.4247E-01 -1.7754E+00 -2.1087E+00
            -1.1930E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1629.38221048439        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      481
 NPARAMETR:  9.4240E-01  6.8639E-01  9.0885E-01  1.2339E+00  7.8056E-01  1.0309E+00  1.6177E+00  7.5148E-02  8.4351E-01  7.0842E-01
             1.1102E+00
 PARAMETER:  4.0674E-02 -2.7630E-01  4.4232E-03  3.1019E-01 -1.4774E-01  1.3041E-01  5.8098E-01 -2.4883E+00 -7.0188E-02 -2.4471E-01
             2.0453E-01
 GRADIENT:   1.5969E+01  4.0517E+00 -1.4112E+00  1.5364E-02 -1.5564E+00  1.7692E-01  2.1850E+00  8.2003E-02 -5.8252E-01  3.0340E+00
             1.8446E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1629.77060809256        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      657
 NPARAMETR:  9.3346E-01  5.3493E-01  9.0197E-01  1.3164E+00  7.3028E-01  1.0286E+00  1.9350E+00  3.3869E-02  8.0313E-01  6.8322E-01
             1.1062E+00
 PARAMETER:  3.1145E-02 -5.2562E-01 -3.1736E-03  3.7493E-01 -2.1433E-01  1.2816E-01  7.6011E-01 -3.2853E+00 -1.1924E-01 -2.8094E-01
             2.0094E-01
 GRADIENT:  -4.3735E-01 -2.1336E-01 -3.6879E-01 -1.0668E-01  4.7945E-01  5.8273E-02  1.4059E-02  1.3805E-02  2.8871E-01 -1.3594E-01
            -1.1905E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1629.77266590476        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      836
 NPARAMETR:  9.3341E-01  5.3636E-01  9.0635E-01  1.3166E+00  7.3305E-01  1.0287E+00  1.9299E+00  2.8781E-02  8.0242E-01  6.8609E-01
             1.1066E+00
 PARAMETER:  3.1089E-02 -5.2296E-01  1.6703E-03  3.7502E-01 -2.1055E-01  1.2825E-01  7.5745E-01 -3.4480E+00 -1.2012E-01 -2.7674E-01
             2.0132E-01
 GRADIENT:  -5.2514E-01 -9.1983E-02 -2.2049E-01  5.5635E-01  7.6835E-01  1.0121E-01 -6.5868E-02  9.5750E-03 -7.0373E-02 -2.0141E-01
            -6.2476E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1629.77761790593        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      999
 NPARAMETR:  9.3369E-01  5.3802E-01  9.0361E-01  1.3151E+00  7.3191E-01  1.0284E+00  1.9264E+00  1.0000E-02  8.0294E-01  6.8542E-01
             1.1065E+00
 PARAMETER:  3.1389E-02 -5.1986E-01 -1.3589E-03  3.7393E-01 -2.1210E-01  1.2802E-01  7.5566E-01 -4.6179E+00 -1.1947E-01 -2.7773E-01
             2.0124E-01
 GRADIENT:   4.9253E-03  9.5912E-03  1.4932E-02  2.6138E-03 -4.1292E-02 -3.0547E-03 -5.8588E-04  0.0000E+00 -5.0521E-03  5.3170E-03
             8.8744E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      999
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.4657E-04  2.2784E-02 -3.9907E-04 -2.0118E-02 -1.4229E-03
 SE:             2.9764E-02  1.9796E-02  2.1220E-04  2.5126E-02  2.1561E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8803E-01  2.4975E-01  6.0020E-02  4.2331E-01  9.4738E-01

 ETASHRINKSD(%)  2.8649E-01  3.3682E+01  9.9289E+01  1.5825E+01  2.7767E+01
 ETASHRINKVR(%)  5.7217E-01  5.6019E+01  9.9995E+01  2.9146E+01  4.7824E+01
 EBVSHRINKSD(%)  5.3892E-01  3.6105E+01  9.9254E+01  1.4629E+01  2.5461E+01
 EBVSHRINKVR(%)  1.0749E+00  5.9175E+01  9.9994E+01  2.7118E+01  4.4439E+01
 RELATIVEINF(%)  9.8436E+01  5.0754E+00  4.7041E-04  1.2518E+01  3.8799E+00
 EPSSHRINKSD(%)  4.1445E+01
 EPSSHRINKVR(%)  6.5713E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1629.7776179059283     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -894.62679134219013     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.97
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1629.778       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.34E-01  5.38E-01  9.04E-01  1.32E+00  7.32E-01  1.03E+00  1.93E+00  1.00E-02  8.03E-01  6.85E-01  1.11E+00
 


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
+        1.19E+03
 
 TH 2
+       -1.30E+01  4.58E+02
 
 TH 3
+        1.50E+01  2.69E+02  8.03E+02
 
 TH 4
+       -6.55E+00  3.09E+02 -2.19E+02  7.41E+02
 
 TH 5
+       -6.02E+00 -5.74E+02 -1.30E+03  2.72E+02  2.44E+03
 
 TH 6
+        2.33E+00 -2.60E+00  3.66E+00 -4.07E+00 -5.31E+00  1.83E+02
 
 TH 7
+        1.22E+00  3.76E+01  2.77E+00 -7.30E+00 -5.43E+00 -3.04E-01  1.59E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.01E-01 -2.28E+01 -2.43E+01  4.21E+00  1.79E+01  2.15E+00  1.04E+01  0.00E+00  1.70E+02
 
 TH10
+       -1.98E+00  4.81E+00 -7.05E+01 -3.52E+01 -3.52E+01  1.31E+00  6.02E+00  0.00E+00  7.27E+00  1.33E+02
 
 TH11
+       -1.02E+01 -9.13E+00 -4.05E+01 -8.67E+00  2.16E+01  4.02E+00  2.31E+00  0.00E+00  1.27E+01  3.61E+01  1.86E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.906
Stop Time:
Sat Sep 25 13:53:36 CDT 2021
