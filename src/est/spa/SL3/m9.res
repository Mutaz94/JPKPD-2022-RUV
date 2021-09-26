Sat Sep 25 11:31:11 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat9.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m9.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1598.40765317160        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7634E+02 -7.0264E+01 -3.0817E+01 -6.0593E+01  7.8880E+01 -1.3452E+00  4.8496E+00  2.8686E+00  4.7338E+00 -1.9194E+01
            -1.4950E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1608.23743708081        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.3155E-01  9.7458E-01  1.0297E+00  1.0576E+00  9.5041E-01  9.8187E-01  8.7951E-01  9.5556E-01  9.4335E-01  1.0973E+00
             1.0149E+00
 PARAMETER:  2.9092E-02  7.4251E-02  1.2925E-01  1.5599E-01  4.9135E-02  8.1702E-02 -2.8391E-02  5.4539E-02  4.1679E-02  1.9288E-01
             1.1483E-01
 GRADIENT:   1.7847E+01 -5.3000E+00 -6.7878E+00 -1.0071E+01  5.8384E+00  8.1879E-01  1.9212E+00  1.6164E+00 -9.7743E+00  2.1703E+00
            -6.6212E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1609.87341886657        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.3267E-01  9.2428E-01  1.1958E+00  1.1031E+00  9.7829E-01  1.0134E+00  4.6382E-01  1.0644E+00  1.1084E+00  1.0952E+00
             1.0353E+00
 PARAMETER:  3.0292E-02  2.1254E-02  2.7884E-01  1.9814E-01  7.8047E-02  1.1332E-01 -6.6825E-01  1.6238E-01  2.0291E-01  1.9092E-01
             1.3473E-01
 GRADIENT:   2.2892E+01  9.4887E+00  1.2858E+01  2.1953E+01 -1.8841E+01  1.3620E+01  2.9466E+00 -3.5512E+00  2.1515E+01 -1.0128E+01
             1.9608E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1612.08818212875        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.2077E-01  1.0561E+00  1.2807E+00  1.0042E+00  1.0829E+00  9.7465E-01  2.3908E-01  1.3267E+00  1.1419E+00  1.2365E+00
             1.0331E+00
 PARAMETER:  1.7454E-02  1.5460E-01  3.4738E-01  1.0416E-01  1.7964E-01  7.4326E-02 -1.3310E+00  3.8270E-01  2.3266E-01  3.1228E-01
             1.3256E-01
 GRADIENT:  -7.9305E+00  2.2694E+00  8.1121E-01 -3.0525E+00 -3.8232E+00 -2.3186E+00  1.2532E+00  1.1153E+00 -3.9227E-01  8.3276E-01
            -5.8075E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1612.78462033243        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.2412E-01  1.0102E+00  1.2921E+00  1.0353E+00  1.0743E+00  9.8225E-01  5.6670E-02  1.2723E+00  1.1200E+00  1.2302E+00
             1.0333E+00
 PARAMETER:  2.1088E-02  1.1018E-01  3.5624E-01  1.3466E-01  1.7167E-01  8.2086E-02 -2.7705E+00  3.4082E-01  2.1332E-01  3.0717E-01
             1.3272E-01
 GRADIENT:   1.4939E+00  9.4185E-01 -1.8428E-01  1.8030E+00  1.0049E+00  1.0596E+00  8.2233E-02 -2.0948E-01 -5.6405E-01 -7.4368E-01
            -9.3592E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1612.78754793821        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  9.2367E-01  1.0150E+00  1.2943E+00  1.0314E+00  1.0766E+00  9.8009E-01  3.6033E-02  1.2813E+00  1.1254E+00  1.2355E+00
             1.0332E+00
 PARAMETER:  2.0601E-02  1.1490E-01  3.5793E-01  1.3090E-01  1.7381E-01  7.9886E-02 -3.2233E+00  3.4785E-01  2.1817E-01  3.1145E-01
             1.3267E-01
 GRADIENT:   2.3861E-01  1.6967E-01 -8.7073E-02  5.7165E-01  3.3451E-01  1.6155E-01  3.5837E-02 -3.6153E-02 -4.7800E-02 -1.4079E-01
            -1.3982E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1612.79383173737        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  9.2357E-01  1.0177E+00  1.2943E+00  1.0292E+00  1.0773E+00  9.7960E-01  1.0000E-02  1.2845E+00  1.1279E+00  1.2369E+00
             1.0332E+00
 PARAMETER:  2.0489E-02  1.1751E-01  3.5798E-01  1.2880E-01  1.7449E-01  7.9390E-02 -4.8711E+00  3.5039E-01  2.2036E-01  3.1264E-01
             1.3267E-01
 GRADIENT:  -6.6535E-02 -3.6812E-02  2.3628E-02 -5.7056E-02 -5.0814E-02 -4.3167E-02  0.0000E+00  5.6733E-03  3.8501E-02  2.4050E-02
             2.8721E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1613.10317824989        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      621
 NPARAMETR:  9.3759E-01  1.0646E+00  1.3022E+00  1.0024E+00  1.1013E+00  9.8733E-01  1.0000E-02  1.3350E+00  1.1655E+00  1.2532E+00
             1.0351E+00
 PARAMETER:  3.5561E-02  1.6262E-01  3.6407E-01  1.0235E-01  1.9645E-01  8.7250E-02 -5.1313E+00  3.8896E-01  2.5317E-01  3.2571E-01
             1.3448E-01
 GRADIENT:  -1.0498E-01  2.2799E-01 -1.1575E-01  3.4534E-01  2.4390E-01  4.1499E-01  0.0000E+00 -6.3837E-02 -3.1315E-01 -2.5536E-01
            -1.0977E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1613.10582753642        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      783
 NPARAMETR:  9.3772E-01  1.0828E+00  1.2999E+00  9.8978E-01  1.1078E+00  9.8637E-01  1.0000E-02  1.3515E+00  1.1811E+00  1.2588E+00
             1.0353E+00
 PARAMETER:  3.5698E-02  1.7953E-01  3.6226E-01  8.9732E-02  2.0234E-01  8.6280E-02 -5.0362E+00  4.0122E-01  2.6642E-01  3.3019E-01
             1.3470E-01
 GRADIENT:   1.7188E-02  7.9172E-03  1.8106E-02 -6.0736E-03 -2.4083E-02  2.2878E-03  0.0000E+00 -2.9386E-03  1.4060E-03 -2.5855E-03
            -1.8362E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      783
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.5402E-04 -8.6960E-04 -3.5035E-02 -3.0503E-03 -3.5009E-02
 SE:             2.9821E-02  2.3448E-04  1.4554E-02  2.9106E-02  2.3136E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7448E-01  2.0844E-04  1.6075E-02  9.1653E-01  1.3023E-01

 ETASHRINKSD(%)  9.5626E-02  9.9214E+01  5.1241E+01  2.4928E+00  2.2491E+01
 ETASHRINKVR(%)  1.9116E-01  9.9994E+01  7.6226E+01  4.9235E+00  3.9924E+01
 EBVSHRINKSD(%)  4.7343E-01  9.9282E+01  5.5957E+01  2.9710E+00  1.8883E+01
 EBVSHRINKVR(%)  9.4463E-01  9.9995E+01  8.0602E+01  5.8537E+00  3.4200E+01
 RELATIVEINF(%)  9.8869E+01  5.2397E-04  6.2865E+00  1.1256E+01  2.2254E+01
 EPSSHRINKSD(%)  4.4441E+01
 EPSSHRINKVR(%)  6.9132E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1613.1058275364187     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -877.95500097268052     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.12
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1613.106       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.38E-01  1.08E+00  1.30E+00  9.90E-01  1.11E+00  9.86E-01  1.00E-02  1.35E+00  1.18E+00  1.26E+00  1.04E+00
 


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
+       -1.63E+01  5.36E+02
 
 TH 3
+        3.84E+00  6.91E+01  7.60E+01
 
 TH 4
+       -1.48E+01  5.30E+02 -6.57E+00  7.66E+02
 
 TH 5
+        1.36E+00 -1.69E+02 -1.20E+02 -2.48E+01  3.83E+02
 
 TH 6
+       -1.48E+00 -1.86E+00  3.63E-01 -4.01E+00 -2.43E+00  1.93E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -1.08E+00 -1.94E+01 -1.71E+01 -2.38E+00  2.63E+00  1.31E+00  0.00E+00  1.43E+01
 
 TH 9
+        3.27E+00 -1.01E+02  1.90E+00  5.08E+00  2.10E+00 -3.39E-01  0.00E+00  2.30E-01  1.27E+02
 
 TH10
+       -1.34E-01 -5.35E+00 -8.75E+00 -3.67E-01 -4.70E+01  8.04E-02  0.00E+00  7.96E+00  5.29E-01  5.87E+01
 
 TH11
+       -8.56E+00 -2.95E+01 -1.38E+01 -1.06E+01 -3.38E-01  5.34E+00  0.00E+00  6.95E+00  8.89E+00  1.01E+01  1.92E+02
 
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
 #CPUT: Total CPU Time in Seconds,       12.655
Stop Time:
Sat Sep 25 11:31:25 CDT 2021
