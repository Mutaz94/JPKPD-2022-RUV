Wed Sep 29 16:56:02 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat86.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1667.63365793835        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6865E+02 -5.5821E+01 -6.2928E+01  2.7704E+01  1.0675E+02  5.2906E+01 -2.1374E+00  8.6048E+00  1.7417E+01  5.9583E+00
            -6.6900E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1679.49461123477        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0123E+00  1.1250E+00  1.1472E+00  9.7259E-01  1.0188E+00  9.7678E-01  1.0005E+00  9.4072E-01  9.2640E-01  8.9309E-01
             1.2901E+00
 PARAMETER:  1.1223E-01  2.1779E-01  2.3734E-01  7.2210E-02  1.1859E-01  7.6509E-02  1.0054E-01  3.8887E-02  2.3555E-02 -1.3063E-02
             3.5469E-01
 GRADIENT:  -5.7700E+01 -2.1641E+00  5.7539E+00 -1.0228E+01 -3.9390E+00 -5.0034E+00 -8.1579E-01  7.5395E-02 -5.7336E+00 -4.0891E+00
             3.6433E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1680.72680049696        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.0338E+00  1.3680E+00  1.0359E+00  8.3353E-01  1.0664E+00  9.7385E-01  8.6899E-01  6.2000E-01  1.1454E+00  9.5326E-01
             1.2570E+00
 PARAMETER:  1.3326E-01  4.1334E-01  1.3532E-01 -8.2083E-02  1.6432E-01  7.3502E-02 -4.0429E-02 -3.7804E-01  2.3574E-01  5.2136E-02
             3.2870E-01
 GRADIENT:  -1.0842E+01  3.4171E+01  2.2006E+01  1.0888E+01 -3.1153E+01 -5.5548E+00  3.6821E+00 -1.4708E+00  4.0686E+00 -4.8751E+00
             2.3112E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1683.56636454749        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      543
 NPARAMETR:  1.0397E+00  1.3619E+00  8.3730E-01  8.2297E-01  1.0069E+00  9.8972E-01  8.9278E-01  4.2614E-01  1.0890E+00  9.2655E-01
             1.1630E+00
 PARAMETER:  1.3897E-01  4.0889E-01 -7.7573E-02 -9.4838E-02  1.0689E-01  8.9665E-02 -1.3418E-02 -7.5299E-01  1.8527E-01  2.3711E-02
             2.5101E-01
 GRADIENT:   1.3920E+00  1.3031E+01  4.3491E+00  9.8739E+00 -6.0992E+00 -6.9488E-02  4.3022E-01  9.8786E-02 -3.3705E-01 -6.6327E-01
            -9.6881E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1684.25623816602        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      718
 NPARAMETR:  1.0397E+00  1.6401E+00  6.2871E-01  6.3798E-01  1.0584E+00  9.9156E-01  7.7964E-01  1.3393E-01  1.2868E+00  9.1612E-01
             1.1654E+00
 PARAMETER:  1.3897E-01  5.9474E-01 -3.6408E-01 -3.4945E-01  1.5674E-01  9.1523E-02 -1.4892E-01 -1.9105E+00  3.5216E-01  1.2393E-02
             2.5305E-01
 GRADIENT:  -1.0421E+00  1.1869E+01 -5.5280E-01  9.8569E+00 -1.4125E+00  7.4792E-02 -1.2696E+00  4.6243E-02 -1.2447E+00 -8.6239E-01
             4.1098E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1684.37446005138        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      894
 NPARAMETR:  1.0399E+00  1.7724E+00  5.7930E-01  5.5166E-01  1.1150E+00  9.8902E-01  7.3708E-01  6.9492E-02  1.4424E+00  9.4665E-01
             1.1598E+00
 PARAMETER:  1.3917E-01  6.7231E-01 -4.4593E-01 -4.9482E-01  2.0883E-01  8.8964E-02 -2.0506E-01 -2.5665E+00  4.6628E-01  4.5177E-02
             2.4828E-01
 GRADIENT:  -4.7194E-01  1.2546E+01  5.4894E-01  8.0455E+00  3.5352E-01 -9.7484E-01 -1.0109E+00  1.0502E-02 -1.2078E+00 -1.3188E+00
            -2.9295E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1684.48520141845        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1079             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0413E+00  1.7867E+00  5.6690E-01  5.3101E-01  1.1238E+00  9.9151E-01  7.3170E-01  1.2466E-02  1.4848E+00  9.5421E-01
             1.1660E+00
 PARAMETER:  1.4048E-01  6.8035E-01 -4.6757E-01 -5.3297E-01  2.1668E-01  9.1477E-02 -2.1239E-01 -4.2847E+00  4.9530E-01  5.3129E-02
             2.5359E-01
 GRADIENT:   4.7674E+02  5.6648E+02  1.7780E+00  8.2147E+01  1.2257E+01  3.5311E+01  8.6634E+00  6.3310E-04  1.4555E+01  3.6274E-01
             2.0470E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1684.48616794090        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1263
 NPARAMETR:  1.0413E+00  1.7871E+00  5.6602E-01  5.3122E-01  1.1229E+00  9.9150E-01  7.3202E-01  1.0000E-02  1.4854E+00  9.5425E-01
             1.1660E+00
 PARAMETER:  1.4044E-01  6.8057E-01 -4.6912E-01 -5.3258E-01  2.1593E-01  9.1465E-02 -2.1195E-01 -5.2907E+00  4.9566E-01  5.3170E-02
             2.5357E-01
 GRADIENT:   2.6321E+00 -8.1362E+00 -2.3650E-01  6.5648E-02  1.0791E+00  9.8986E-02 -1.4523E-03  0.0000E+00  2.0710E-01  3.4840E-02
             6.7343E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1684.48616794090        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1285
 NPARAMETR:  1.0413E+00  1.7871E+00  5.6602E-01  5.3122E-01  1.1229E+00  9.9150E-01  7.3202E-01  1.0000E-02  1.4854E+00  9.5425E-01
             1.1660E+00
 PARAMETER:  1.4044E-01  6.8057E-01 -4.6912E-01 -5.3258E-01  2.1593E-01  9.1465E-02 -2.1195E-01 -5.2907E+00  4.9566E-01  5.3170E-02
             2.5357E-01
 GRADIENT:   2.6321E+00 -8.1362E+00 -2.3650E-01  6.5648E-02  1.0791E+00  9.8986E-02 -1.4523E-03  0.0000E+00  2.0710E-01  3.4840E-02
             6.7343E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1285
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0776E-04 -3.0083E-02 -2.2496E-04  2.4074E-02 -3.5511E-02
 SE:             2.9791E-02  2.2880E-02  8.7594E-05  2.2356E-02  2.2316E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9711E-01  1.8857E-01  1.0224E-02  2.8154E-01  1.1154E-01

 ETASHRINKSD(%)  1.9483E-01  2.3350E+01  9.9707E+01  2.5105E+01  2.5240E+01
 ETASHRINKVR(%)  3.8929E-01  4.1247E+01  9.9999E+01  4.3908E+01  4.4109E+01
 EBVSHRINKSD(%)  5.8710E-01  2.2686E+01  9.9734E+01  2.6944E+01  2.4174E+01
 EBVSHRINKVR(%)  1.1708E+00  4.0225E+01  9.9999E+01  4.6628E+01  4.2504E+01
 RELATIVEINF(%)  9.8746E+01  3.6779E+00  8.2154E-05  3.4066E+00  1.2622E+01
 EPSSHRINKSD(%)  4.2688E+01
 EPSSHRINKVR(%)  6.7153E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1684.4861679409014     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -949.33534137716322     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.46
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1684.486       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.79E+00  5.66E-01  5.31E-01  1.12E+00  9.92E-01  7.32E-01  1.00E-02  1.49E+00  9.54E-01  1.17E+00
 


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
+        1.03E+03
 
 TH 2
+       -7.97E+00  4.09E+02
 
 TH 3
+        9.09E+00  1.30E+02  2.60E+02
 
 TH 4
+       -2.03E+01  3.61E+02 -2.30E+02  9.79E+02
 
 TH 5
+       -5.05E+00 -1.92E+02 -2.81E+02  2.29E+02  5.58E+02
 
 TH 6
+        9.93E-02 -1.71E+00  1.66E+00 -4.87E+00 -1.66E+00  1.98E+02
 
 TH 7
+        1.24E+00  8.49E+00  6.42E+00 -2.25E+01 -1.42E+01 -5.25E-01  1.45E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.31E-01 -1.84E+01 -2.79E+01  5.80E+01  6.46E-01 -3.23E-01  1.73E+01  0.00E+00  3.52E+01
 
 TH10
+       -9.21E-01 -1.19E+01 -2.81E+01 -5.88E+00 -5.66E+01  7.94E-01  1.77E+01  0.00E+00  4.86E+00  8.01E+01
 
 TH11
+       -7.03E+00 -2.01E+01 -2.40E+01  1.29E+00 -4.95E+00  2.52E+00  1.15E+01  0.00E+00  5.52E+00  1.84E+01  1.61E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.751
Stop Time:
Wed Sep 29 16:56:27 CDT 2021
