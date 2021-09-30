Wed Sep 29 12:06:28 CDT 2021
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
$DATA ../../../../data/spa/A1/dat39.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m39.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1323.47639664066        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6930E+02  6.7309E+01  1.9956E+01  6.1927E+01  3.2164E+01  5.5402E+01 -1.4184E+01  1.4959E+01 -2.2578E+01 -3.4892E+01
            -6.1488E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1462.32751158300        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1198E+00  9.5259E-01  1.0277E+00  1.0875E+00  9.4961E-01  9.9001E-01  1.0098E+00  8.8903E-01  1.0706E+00  1.0157E+00
             2.4631E+00
 PARAMETER:  2.1316E-01  5.1425E-02  1.2731E-01  1.8385E-01  4.8298E-02  8.9962E-02  1.0973E-01 -1.7620E-02  1.6823E-01  1.1561E-01
             1.0014E+00
 GRADIENT:   3.9797E+02  3.5658E+01  3.1133E+00  5.5411E+01 -2.7098E+01 -4.4248E-01  5.9548E+00  9.5829E+00  1.1250E+01  1.4603E+01
             1.0613E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1479.78427336157        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0986E+00  8.6127E-01  7.1981E-01  1.1214E+00  8.0252E-01  1.0310E+00  9.0981E-01  1.5773E-01  1.1008E+00  8.4899E-01
             2.1815E+00
 PARAMETER:  1.9406E-01 -4.9349E-02 -2.2877E-01  2.1461E-01 -1.1999E-01  1.3051E-01  5.4754E-03 -1.7469E+00  1.9599E-01 -6.3712E-02
             8.8000E-01
 GRADIENT:   3.6332E+02  2.4878E-01 -6.2478E+01  8.8202E+01  8.1859E+01  2.2559E+01  7.6590E-01  6.1341E-01  1.5301E+01  9.9772E+00
             6.3993E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1501.00641039262        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  9.8866E-01  6.3627E-01  5.9135E-01  1.2209E+00  6.0062E-01  8.9507E-01  1.2058E+00  1.2928E-01  9.8023E-01  7.3602E-01
             1.7957E+00
 PARAMETER:  8.8597E-02 -3.5213E-01 -4.2534E-01  2.9957E-01 -4.0978E-01 -1.0851E-02  2.8716E-01 -1.9458E+00  8.0032E-02 -2.0650E-01
             6.8539E-01
 GRADIENT:  -2.2860E+01  2.2422E+01 -4.4483E+01  7.0484E+01  4.0641E+01 -2.1412E+01 -7.9571E-01  6.2277E-01  2.1641E+00  1.1693E+01
            -1.7861E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1508.08592832762        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  9.7709E-01  3.9718E-01  8.2423E-01  1.3804E+00  6.7018E-01  9.2796E-01  1.6855E+00  8.0996E-02  9.0651E-01  7.8683E-01
             1.8569E+00
 PARAMETER:  7.6820E-02 -8.2336E-01 -9.3300E-02  4.2240E-01 -3.0022E-01  2.5237E-02  6.2205E-01 -2.4134E+00  1.8486E-03 -1.3975E-01
             7.1890E-01
 GRADIENT:  -3.2714E+01  1.5238E+01 -2.6606E+00  5.3198E+01  1.3770E+01 -4.2046E+00  1.0637E+00  1.5580E-01  3.9650E-01 -4.1136E+00
            -3.6848E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1515.12597415860        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      608
 NPARAMETR:  9.9193E-01  1.1788E-01  6.0357E-01  1.4312E+00  4.8714E-01  9.8526E-01  1.4716E+00  1.0000E-02  8.9038E-01  6.3857E-01
             1.8454E+00
 PARAMETER:  9.1896E-02 -2.0381E+00 -4.0490E-01  4.5852E-01 -6.1921E-01  8.5148E-02  4.8637E-01 -6.5151E+00 -1.6110E-02 -3.4853E-01
             7.1271E-01
 GRADIENT:   1.8111E+01 -5.9315E-01 -2.8482E+01 -1.6851E+01  5.2618E+01  1.8067E+01 -4.3827E-01  0.0000E+00  8.7544E+00 -9.9421E+00
            -2.7498E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1517.91396110813        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      783
 NPARAMETR:  9.8598E-01  5.4300E-02  5.2363E-01  1.4348E+00  4.2116E-01  9.3711E-01  1.4249E+00  1.0000E-02  8.6170E-01  6.4927E-01
             1.8229E+00
 PARAMETER:  8.5878E-02 -2.8132E+00 -5.4697E-01  4.6104E-01 -7.6473E-01  3.5043E-02  4.5409E-01 -9.3499E+00 -4.8844E-02 -3.3191E-01
             7.0043E-01
 GRADIENT:   7.5261E+00  1.3268E+00  5.9082E+00 -1.0203E+00 -1.1783E+01 -2.0524E-01 -1.1478E-01  0.0000E+00 -7.0451E-02  3.3449E-01
            -8.1087E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1518.83747839642        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      959
 NPARAMETR:  9.8052E-01  1.0000E-02  5.5447E-01  1.4688E+00  4.3484E-01  9.3628E-01  1.4896E+00  1.0000E-02  8.4835E-01  6.5926E-01
             1.8214E+00
 PARAMETER:  8.0327E-02 -4.7531E+00 -4.8974E-01  4.8442E-01 -7.3278E-01  3.4159E-02  4.9854E-01 -1.6390E+01 -6.4461E-02 -3.1664E-01
             6.9960E-01
 GRADIENT:   1.4526E-01  0.0000E+00  1.1460E+00 -4.7665E-01 -1.1798E+00  1.7139E-01 -3.6530E-03  0.0000E+00  5.2209E-01 -5.8178E-01
            -4.1839E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1518.83795264807        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1016
 NPARAMETR:  9.8047E-01  1.0000E-02  5.5521E-01  1.4692E+00  4.3544E-01  9.3593E-01  1.4919E+00  1.0000E-02  8.4715E-01  6.6189E-01
             1.8218E+00
 PARAMETER:  8.0279E-02 -4.7357E+00 -4.8840E-01  4.8472E-01 -7.3141E-01  3.3781E-02  5.0006E-01 -1.6327E+01 -6.5875E-02 -3.1266E-01
             6.9983E-01
 GRADIENT:   9.0484E-03  0.0000E+00  2.0021E-01 -3.0404E-01 -1.6886E-01  3.0348E-02 -3.6671E-03  0.0000E+00  1.0039E-01 -7.6536E-02
            -4.0725E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1016
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1736E-04 -4.5042E-04 -1.9912E-05 -8.0004E-03 -9.8232E-03
 SE:             2.9470E-02  2.8282E-04  2.4427E-04  2.8257E-02  2.2882E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9682E-01  1.1125E-01  9.3503E-01  7.7708E-01  6.6771E-01

 ETASHRINKSD(%)  1.2704E+00  9.9053E+01  9.9182E+01  5.3341E+00  2.3341E+01
 ETASHRINKVR(%)  2.5246E+00  9.9991E+01  9.9993E+01  1.0384E+01  4.1234E+01
 EBVSHRINKSD(%)  1.4329E+00  9.9196E+01  9.9120E+01  5.0035E+00  2.2780E+01
 EBVSHRINKVR(%)  2.8452E+00  9.9994E+01  9.9992E+01  9.7566E+00  4.0370E+01
 RELATIVEINF(%)  8.6052E+01  4.0449E-04  3.0325E-04  1.2710E+01  1.5080E+00
 EPSSHRINKSD(%)  3.6648E+01
 EPSSHRINKVR(%)  5.9866E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1518.8379526480692     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -783.68712608433100     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1518.838       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.00E-02  5.55E-01  1.47E+00  4.35E-01  9.36E-01  1.49E+00  1.00E-02  8.47E-01  6.62E-01  1.82E+00
 


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
+        1.28E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.03E+01  0.00E+00  2.81E+03
 
 TH 4
+       -2.58E+01  0.00E+00 -1.96E+02  6.48E+02
 
 TH 5
+        6.11E+01  0.00E+00 -4.52E+03 -1.43E+02  7.96E+03
 
 TH 6
+       -2.16E-01  0.00E+00  7.11E+00 -6.49E+00 -1.37E-01  2.14E+02
 
 TH 7
+        8.78E-03  0.00E+00 -3.29E-03 -4.01E-03  8.15E-03  2.05E-03 -2.39E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.19E+00  0.00E+00  3.47E+01 -1.04E+01 -5.74E-02  2.59E+00  6.20E-03  0.00E+00  2.28E+02
 
 TH10
+       -5.33E+00  0.00E+00 -3.73E+01  1.22E-01 -6.75E+01  2.27E-01  4.05E-04  0.00E+00  3.81E-01  1.65E+02
 
 TH11
+       -1.19E+01  0.00E+00 -2.15E+01 -7.06E+00  7.33E+00  1.85E+00  1.15E-03  0.00E+00  7.23E+00  3.55E+01  7.88E+01
 
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
 #CPUT: Total CPU Time in Seconds,       16.937
Stop Time:
Wed Sep 29 12:06:47 CDT 2021
