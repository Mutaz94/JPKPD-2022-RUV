Sat Sep 25 13:56:14 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat97.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1645.83087502219        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.9771E+02 -6.5067E+01 -5.2712E+01 -2.0995E+01  5.5524E+01 -2.0748E+01 -4.1455E-01  1.3468E+01  2.0477E+01  1.7723E+01
             2.1287E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1650.62947837762        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9562E-01  1.0897E+00  1.1477E+00  9.6065E-01  1.0884E+00  1.0794E+00  1.0110E+00  8.9823E-01  8.7365E-01  8.6491E-01
             9.1941E-01
 PARAMETER:  9.5612E-02  1.8590E-01  2.3778E-01  5.9857E-02  1.8474E-01  1.7641E-01  1.1093E-01 -7.3240E-03 -3.5078E-02 -4.5134E-02
             1.5975E-02
 GRADIENT:   1.7646E+02 -2.3987E+01 -5.8671E+00 -1.4159E+01  7.3793E+01  1.8299E+01 -8.5630E+00 -1.7187E+00 -1.0165E+01 -2.4233E+01
            -2.7359E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1654.13126351437        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      263
 NPARAMETR:  9.8656E-01  1.1125E+00  9.0176E-01  9.3883E-01  9.9716E-01  1.1063E+00  1.1343E+00  5.2777E-01  7.9839E-01  8.7136E-01
             9.2555E-01
 PARAMETER:  8.6471E-02  2.0661E-01 -3.4014E-03  3.6876E-02  9.7157E-02  2.0098E-01  2.2604E-01 -5.3909E-01 -1.2516E-01 -3.7695E-02
             2.2632E-02
 GRADIENT:   1.0383E+02 -2.3817E+01 -1.9890E+01 -1.0994E+01  6.3737E+01  1.6684E+01  1.9372E+00  8.1476E-01 -1.5583E+01 -7.0482E+00
            -1.8752E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1661.26081861072        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  9.3700E-01  8.8080E-01  8.1846E-01  1.0946E+00  8.0665E-01  1.0369E+00  1.3296E+00  3.6369E-01  7.9431E-01  6.5851E-01
             9.5710E-01
 PARAMETER:  3.4930E-02 -2.6921E-02 -1.0033E-01  1.9041E-01 -1.1486E-01  1.3625E-01  3.8486E-01 -9.1145E-01 -1.3028E-01 -3.1778E-01
             5.6157E-02
 GRADIENT:   1.1920E+01  1.1419E+01  1.5039E+00  2.9634E+01  1.5418E+01 -1.7575E+00 -4.1294E+00  2.2197E-01 -1.3343E+00 -1.3575E+01
            -2.7731E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1662.89544453481        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  9.3145E-01  7.9894E-01  7.7416E-01  1.1323E+00  7.4890E-01  1.0417E+00  1.4609E+00  1.5810E-01  7.5987E-01  7.1121E-01
             9.4676E-01
 PARAMETER:  2.8987E-02 -1.2447E-01 -1.5597E-01  2.2423E-01 -1.8915E-01  1.4088E-01  4.7904E-01 -1.7445E+00 -1.7460E-01 -2.4079E-01
             4.5293E-02
 GRADIENT:  -5.0473E-01  1.0381E+01  1.6805E+00  1.2535E+01 -7.2419E+00 -1.1135E-01  3.3787E-01  2.6040E-01 -8.5149E-01  1.0112E+00
            -1.3618E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1663.21654781391        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      795
 NPARAMETR:  9.3087E-01  6.8518E-01  7.8193E-01  1.1889E+00  7.1721E-01  1.0414E+00  1.6366E+00  7.5798E-02  7.3459E-01  7.0780E-01
             9.4778E-01
 PARAMETER:  2.8369E-02 -2.7807E-01 -1.4600E-01  2.7301E-01 -2.3238E-01  1.4057E-01  5.9264E-01 -2.4797E+00 -2.0845E-01 -2.4560E-01
             4.6364E-02
 GRADIENT:  -3.5770E-02  5.1376E-02 -1.8936E-01  1.8665E-01 -2.9804E-01  1.7187E-03  1.5424E-02  4.5719E-02  1.4614E-02  1.7576E-01
            -1.8918E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1663.23353161536        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      974
 NPARAMETR:  9.3105E-01  6.8683E-01  7.7721E-01  1.1866E+00  7.1521E-01  1.0413E+00  1.6336E+00  3.2378E-02  7.3568E-01  7.0514E-01
             9.4759E-01
 PARAMETER:  2.8562E-02 -2.7567E-01 -1.5205E-01  2.7108E-01 -2.3518E-01  1.4049E-01  5.9077E-01 -3.3303E+00 -2.0697E-01 -2.4936E-01
             4.6168E-02
 GRADIENT:   2.2576E-01 -3.9439E-01  7.3385E-02 -1.3579E+00 -6.7969E-01 -6.1315E-02  8.2854E-03  8.5166E-03  2.3231E-01  1.9264E-01
            -3.6965E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1663.23885210424        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1149
 NPARAMETR:  9.3096E-01  6.9179E-01  7.8015E-01  1.1848E+00  7.1867E-01  1.0415E+00  1.6251E+00  1.0000E-02  7.3621E-01  7.0746E-01
             9.4798E-01
 PARAMETER:  2.8456E-02 -2.6848E-01 -1.4827E-01  2.6960E-01 -2.3035E-01  1.4065E-01  5.8559E-01 -4.8152E+00 -2.0624E-01 -2.4608E-01
             4.6583E-02
 GRADIENT:  -1.5048E-03  1.2803E-02 -7.8903E-03  2.1306E-02 -9.2646E-03  3.5437E-03 -1.4523E-04  0.0000E+00  2.8225E-04  6.3918E-03
            -2.6638E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1663.23885210424        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1171
 NPARAMETR:  9.3096E-01  6.9179E-01  7.8015E-01  1.1848E+00  7.1867E-01  1.0415E+00  1.6251E+00  1.0000E-02  7.3621E-01  7.0746E-01
             9.4798E-01
 PARAMETER:  2.8456E-02 -2.6848E-01 -1.4827E-01  2.6960E-01 -2.3035E-01  1.4065E-01  5.8559E-01 -4.8152E+00 -2.0624E-01 -2.4608E-01
             4.6583E-02
 GRADIENT:  -1.5048E-03  1.2803E-02 -7.8903E-03  2.1306E-02 -9.2646E-03  3.5437E-03 -1.4523E-04  0.0000E+00  2.8225E-04  6.3918E-03
            -2.6638E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1171
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3274E-04  1.6286E-02 -5.0044E-04 -1.6914E-02  1.0003E-03
 SE:             2.9860E-02  2.1736E-02  2.2455E-04  2.4316E-02  2.2450E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9378E-01  4.5370E-01  2.5841E-02  4.8670E-01  9.6446E-01

 ETASHRINKSD(%)  1.0000E-10  2.7182E+01  9.9248E+01  1.8537E+01  2.4790E+01
 ETASHRINKVR(%)  1.0000E-10  4.6975E+01  9.9994E+01  3.3638E+01  4.3435E+01
 EBVSHRINKSD(%)  3.6112E-01  2.7711E+01  9.9273E+01  1.8041E+01  2.3481E+01
 EBVSHRINKVR(%)  7.2094E-01  4.7744E+01  9.9995E+01  3.2828E+01  4.1449E+01
 RELATIVEINF(%)  9.8901E+01  5.7319E+00  4.5872E-04  8.8945E+00  4.1411E+00
 EPSSHRINKSD(%)  4.3510E+01
 EPSSHRINKVR(%)  6.8089E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1663.2388521042403     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -928.08802554050214     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.52
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1663.239       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.31E-01  6.92E-01  7.80E-01  1.18E+00  7.19E-01  1.04E+00  1.63E+00  1.00E-02  7.36E-01  7.07E-01  9.48E-01
 


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
+        1.17E+03
 
 TH 2
+       -8.87E+00  4.81E+02
 
 TH 3
+        1.84E+01  2.95E+02  1.04E+03
 
 TH 4
+       -8.27E+00  3.59E+02 -3.98E+02  9.95E+02
 
 TH 5
+       -5.56E+00 -5.75E+02 -1.48E+03  4.98E+02  2.49E+03
 
 TH 6
+        7.31E-01 -1.08E+00  5.27E+00 -2.96E+00 -1.67E+00  1.82E+02
 
 TH 7
+        1.22E+00  3.76E+01 -4.84E+00 -1.12E+01 -8.39E-01 -2.85E-02  2.69E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.35E+00 -2.43E+01 -3.70E+01  1.17E+01  2.18E+01  1.82E-01  1.36E+01  0.00E+00  1.85E+02
 
 TH10
+       -2.35E+00 -3.36E+00 -1.01E+02 -4.44E+01 -2.86E+01 -3.27E-02  1.23E+01  0.00E+00  1.41E+01  1.41E+02
 
 TH11
+       -1.15E+01 -1.24E+01 -4.34E+01 -4.81E+00  1.67E+01  1.74E+00  2.69E+00  0.00E+00  1.21E+01  3.16E+01  2.45E+02
 
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
 #CPUT: Total CPU Time in Seconds,       20.653
Stop Time:
Sat Sep 25 13:56:37 CDT 2021
