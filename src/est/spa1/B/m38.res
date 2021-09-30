Wed Sep 29 21:05:29 CDT 2021
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
$DATA ../../../../data/spa1/B/dat38.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1593.56142652478        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5181E+02  1.0644E+01  5.1324E+01  3.5982E+01 -1.6412E+00  2.1902E+01 -1.9852E+01 -1.7653E+02 -3.7895E+01  6.6561E+00
            -7.6969E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2070.44710720867        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  1.1120E+00  1.0025E+00  1.0709E+00  9.7605E-01  1.0502E+00  1.0994E+00  1.0819E+00  1.0597E+00  1.0830E+00  9.3125E-01
             1.0293E+00
 PARAMETER:  2.0616E-01  1.0250E-01  1.6846E-01  7.5753E-02  1.4896E-01  1.9475E-01  1.7875E-01  1.5802E-01  1.7973E-01  2.8767E-02
             1.2888E-01
 GRADIENT:   1.1059E+03 -3.4644E+00 -1.1877E+01  3.4783E+01  9.7789E+00  7.7466E+01 -8.0770E+00  1.2726E+01  2.5527E+00  1.4693E+00
            -2.0345E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2077.43154760915        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0805E+00  8.5497E-01  1.1198E+00  1.0825E+00  1.0578E+00  1.0852E+00  1.3628E+00  5.9878E-01  9.8542E-01  1.0171E+00
             9.7928E-01
 PARAMETER:  1.7745E-01 -5.6690E-02  2.1314E-01  1.7927E-01  1.5617E-01  1.8177E-01  4.0951E-01 -4.1286E-01  8.5309E-02  1.1694E-01
             7.9065E-02
 GRADIENT:   1.0106E+03  9.9840E-02 -1.8645E+01  1.8011E+02  5.9444E+01  8.7005E+01  1.1612E+01  1.5697E+00 -3.4015E+00  4.4391E-01
            -4.5210E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2081.58058630552        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  1.0069E+00  8.9004E-01  9.7413E-01  1.0363E+00  9.8112E-01  1.0079E+00  1.3376E+00  5.1221E-01  1.0210E+00  9.2855E-01
             9.3701E-01
 PARAMETER:  1.0685E-01 -1.6491E-02  7.3787E-02  1.3566E-01  8.0944E-02  1.0791E-01  3.9086E-01 -5.6902E-01  1.2080E-01  2.5866E-02
             3.4934E-02
 GRADIENT:   5.8570E+02 -1.9782E+00 -1.7360E+01  1.1951E+02  3.6258E+01  4.5751E+01  1.2392E+01  3.2583E+00  5.0143E+00  1.4385E+00
            -8.2820E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2088.74257127194        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.8875E-01  1.0170E+00  9.3755E-01  9.9969E-01  9.9359E-01  1.0579E+00  1.3533E+00  4.0117E-01  1.0627E+00  9.0398E-01
             1.0243E+00
 PARAMETER:  8.8681E-02  1.1689E-01  3.5517E-02  9.9685E-02  9.3572E-02  1.5630E-01  4.0256E-01 -8.1337E-01  1.6078E-01 -9.4781E-04
             1.2397E-01
 GRADIENT:  -7.3363E+00  3.5408E+00 -3.2610E+00  9.6125E+00  5.6420E-01 -3.0766E+00  1.7329E+00  1.4463E+00 -3.5353E+00  5.5885E-01
            -3.7771E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2090.14228661581        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      546
 NPARAMETR:  9.9545E-01  1.2131E+00  8.4282E-01  8.7119E-01  1.0466E+00  1.0625E+00  1.1646E+00  5.5282E-02  1.1877E+00  9.2263E-01
             1.0242E+00
 PARAMETER:  9.5437E-02  2.9321E-01 -7.1004E-02 -3.7889E-02  1.4554E-01  1.6066E-01  2.5237E-01 -2.7953E+00  2.7202E-01  1.9468E-02
             1.2394E-01
 GRADIENT:   2.6917E+00  7.4636E-01 -1.8192E+00  4.5514E+00  7.5210E-01 -2.0032E+00  5.0764E-01  3.2436E-02 -3.4781E-01  9.4332E-01
            -2.5809E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2090.31637961112        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  9.9414E-01  1.4023E+00  7.9310E-01  7.5306E-01  1.1274E+00  1.0702E+00  1.0313E+00  1.0000E-02  1.3292E+00  9.5759E-01
             1.0295E+00
 PARAMETER:  9.4124E-02  4.3815E-01 -1.3181E-01 -1.8361E-01  2.1996E-01  1.6786E-01  1.3080E-01 -4.7141E+00  3.8454E-01  5.6660E-02
             1.2904E-01
 GRADIENT:  -1.7473E+00  3.5347E+00  1.1750E+00  2.6789E+00 -1.4808E+00  5.3989E-01 -1.0355E-01  0.0000E+00 -1.4431E-01 -4.2110E-01
             4.7036E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2090.33109027995        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      911
 NPARAMETR:  9.9583E-01  1.4013E+00  7.9090E-01  7.4946E-01  1.1286E+00  1.0704E+00  1.0305E+00  1.0000E-02  1.3318E+00  9.5968E-01
             1.0292E+00
 PARAMETER:  9.5818E-02  4.3741E-01 -1.3458E-01 -1.8841E-01  2.2096E-01  1.6800E-01  1.3003E-01 -4.7432E+00  3.8653E-01  5.8848E-02
             1.2876E-01
 GRADIENT:   1.5341E+00 -1.7141E+00  2.9428E-01 -5.0168E-01 -1.7706E-01  6.0018E-01 -1.8007E-02  0.0000E+00  5.4651E-02  2.2768E-03
             7.7525E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2090.33109027995        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:      939
 NPARAMETR:  9.9599E-01  1.4012E+00  7.8984E-01  7.4990E-01  1.1287E+00  1.0704E+00  1.0308E+00  1.0000E-02  1.3323E+00  9.5967E-01
             1.0290E+00
 PARAMETER:  9.5818E-02  4.3741E-01 -1.3458E-01 -1.8841E-01  2.2096E-01  1.6800E-01  1.3003E-01 -4.7432E+00  3.8653E-01  5.8848E-02
             1.2876E-01
 GRADIENT:  -2.2247E-01  8.5082E-02  1.7106E-01 -3.3178E-01 -1.3736E-01 -1.2885E-02 -3.2311E-02  0.0000E+00 -5.6437E-02  1.6112E-03
             7.8214E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      939
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.1736E-04 -2.0796E-02 -3.5498E-04  1.6339E-02 -2.9984E-02
 SE:             2.9878E-02  2.2555E-02  1.2740E-04  2.4196E-02  2.2260E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8885E-01  3.5653E-01  5.3303E-03  4.9951E-01  1.7799E-01

 ETASHRINKSD(%)  1.0000E-10  2.4437E+01  9.9573E+01  1.8939E+01  2.5426E+01
 ETASHRINKVR(%)  1.0000E-10  4.2902E+01  9.9998E+01  3.4291E+01  4.4387E+01
 EBVSHRINKSD(%)  3.1986E-01  2.3892E+01  9.9595E+01  1.9710E+01  2.4220E+01
 EBVSHRINKVR(%)  6.3871E-01  4.2076E+01  9.9998E+01  3.5536E+01  4.2575E+01
 RELATIVEINF(%)  9.9202E+01  4.6176E+00  2.9981E-04  5.7539E+00  1.1092E+01
 EPSSHRINKSD(%)  3.3050E+01
 EPSSHRINKVR(%)  5.5176E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2090.3310902799471     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1171.3925570752745     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.77
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2090.331       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  1.40E+00  7.91E-01  7.49E-01  1.13E+00  1.07E+00  1.03E+00  1.00E-02  1.33E+00  9.60E-01  1.03E+00
 


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
+        9.72E+02
 
 TH 2
+       -4.39E+00  3.26E+02
 
 TH 3
+        5.44E+00  1.36E+02  2.96E+02
 
 TH 4
+       -2.32E+00  2.58E+02 -1.92E+02  7.28E+02
 
 TH 5
+        8.94E-01 -1.88E+02 -3.03E+02  2.20E+02  5.42E+02
 
 TH 6
+        8.91E-01 -7.45E-01  1.87E+00 -1.37E+00 -1.91E-01  1.72E+02
 
 TH 7
+        8.27E-01  1.97E+01 -2.28E+00 -1.39E+01 -6.44E+00 -3.85E-02  7.00E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.06E-01 -2.11E+01 -2.43E+01  4.31E+01 -1.40E-01 -3.83E-01  1.42E+01  0.00E+00  5.68E+01
 
 TH10
+        9.20E-01 -9.67E+00 -3.37E+01 -9.02E+00 -5.34E+01  2.01E-01  1.29E+01  0.00E+00  5.95E+00  7.74E+01
 
 TH11
+       -7.23E+00 -1.75E+01 -3.48E+01 -6.17E-01 -1.41E+00  1.86E+00  6.98E+00  0.00E+00  2.54E+00  2.05E+01  3.90E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.184
Stop Time:
Wed Sep 29 21:05:52 CDT 2021
