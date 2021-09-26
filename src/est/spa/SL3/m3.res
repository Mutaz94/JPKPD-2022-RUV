Sat Sep 25 11:28:28 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat3.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1675.67867300047        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.3411E+01 -8.8917E+01 -3.4793E+01 -7.7954E+01  6.9957E+01  6.4268E+00  8.1484E-02 -3.0324E-01  8.4927E-02  1.7675E+01
            -4.4241E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1686.53512069452        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0379E+00  1.0950E+00  1.0202E+00  1.0223E+00  9.6426E-01  9.8461E-01  9.3350E-01  1.0560E+00  9.8141E-01  7.2426E-01
             1.1461E+00
 PARAMETER:  1.3723E-01  1.9073E-01  1.2001E-01  1.2202E-01  6.3606E-02  8.4493E-02  3.1187E-02  1.5448E-01  8.1239E-02 -2.2260E-01
             2.3633E-01
 GRADIENT:   6.6431E+01  2.3018E+01  1.2868E+01  2.5672E+01  2.7042E+00 -7.2539E-01 -2.7086E+00 -9.8771E+00 -5.4126E+00 -6.3510E+00
             6.1479E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1687.49484025621        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0368E+00  1.1609E+00  9.0131E-01  9.8085E-01  9.3231E-01  9.9059E-01  9.7371E-01  1.2970E+00  1.0064E+00  6.3402E-01
             1.1279E+00
 PARAMETER:  1.3614E-01  2.4921E-01 -3.9100E-03  8.0665E-02  2.9906E-02  9.0546E-02  7.3363E-02  3.6002E-01  1.0638E-01 -3.5568E-01
             2.2032E-01
 GRADIENT:   6.2783E+01  2.6802E+01 -3.7317E+00  3.5894E+01 -9.5012E-01  1.2854E+00  3.5452E+00  2.5095E+00  2.3795E-01 -2.6629E+00
             4.1004E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1688.22332696791        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0185E+00  1.1998E+00  9.5161E-01  9.4881E-01  9.7617E-01  9.8782E-01  8.6088E-01  1.3155E+00  1.0584E+00  7.0827E-01
             1.1182E+00
 PARAMETER:  1.1832E-01  2.8219E-01  5.0401E-02  4.7453E-02  7.5877E-02  8.7747E-02 -4.9801E-02  3.7419E-01  1.5679E-01 -2.4493E-01
             2.1171E-01
 GRADIENT:   1.8348E+01  1.9134E+01  3.6081E-01  2.0387E+01 -6.2472E-01  3.8790E-01 -2.2074E+00 -4.9256E-01  7.1244E-01 -1.9756E+00
             1.1088E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1688.22339594355        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  1.0182E+00  1.1984E+00  9.5204E-01  9.4923E-01  9.7586E-01  9.8779E-01  8.6223E-01  1.3144E+00  1.0575E+00  7.0870E-01
             1.1180E+00
 PARAMETER:  1.1808E-01  2.8102E-01  5.0851E-02  4.7898E-02  7.5561E-02  8.7717E-02 -4.8231E-02  3.7341E-01  1.5593E-01 -2.4432E-01
             2.1151E-01
 GRADIENT:   1.7765E+01  1.8535E+01  3.5407E-01  1.9746E+01 -6.1267E-01  3.7602E-01 -2.1387E+00 -4.7759E-01  6.8875E-01 -1.9131E+00
             1.0731E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1689.35070237543        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      488
 NPARAMETR:  1.0302E+00  1.4593E+00  7.8574E-01  7.7676E-01  1.0263E+00  1.0001E+00  7.8714E-01  1.3276E+00  1.2051E+00  7.4342E-01
             1.1150E+00
 PARAMETER:  1.2977E-01  4.7794E-01 -1.4112E-01 -1.5262E-01  1.2592E-01  1.0009E-01 -1.3935E-01  3.8336E-01  2.8654E-01 -1.9650E-01
             2.0886E-01
 GRADIENT:  -5.7040E-01  1.5937E+01  4.5931E+00  9.9416E+00 -7.7206E+00 -5.1714E-01 -1.1845E+00 -1.0223E+00 -7.8331E-01  3.2888E-02
            -1.7966E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1690.28862925808        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      664
 NPARAMETR:  1.0316E+00  1.6836E+00  5.5211E-01  6.1901E-01  1.0382E+00  1.0028E+00  7.4299E-01  1.2086E+00  1.3691E+00  7.2206E-01
             1.1154E+00
 PARAMETER:  1.3110E-01  6.2092E-01 -4.9402E-01 -3.7964E-01  1.3748E-01  1.0282E-01 -1.9707E-01  2.8945E-01  4.1418E-01 -2.2565E-01
             2.0925E-01
 GRADIENT:   8.9021E-01  1.5196E+01  2.3011E+00  7.1981E+00 -5.9158E+00  2.0618E-01 -5.2133E-01 -2.1402E-01 -6.9328E-01 -6.3850E-02
            -8.0773E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1690.37290910544        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  1.0321E+00  1.7083E+00  5.2780E-01  6.0203E-01  1.0393E+00  1.0024E+00  7.3921E-01  1.1976E+00  1.3894E+00  7.2080E-01
             1.1152E+00
 PARAMETER:  1.3159E-01  6.3547E-01 -5.3904E-01 -4.0746E-01  1.3852E-01  1.0238E-01 -2.0217E-01  2.8028E-01  4.2887E-01 -2.2740E-01
             2.0908E-01
 GRADIENT:   4.9720E+01  6.8934E+01  2.7282E+00  1.6919E+01 -6.0519E+00  5.5660E+00  8.8241E-01 -3.1132E-02  1.0480E+00  2.5059E-01
             1.0844E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1690.37332586315        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      992             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0321E+00  1.7083E+00  5.2780E-01  6.0203E-01  1.0393E+00  1.0023E+00  7.3943E-01  1.2008E+00  1.3894E+00  7.1998E-01
             1.1154E+00
 PARAMETER:  1.3159E-01  6.3547E-01 -5.3904E-01 -4.0746E-01  1.3852E-01  1.0229E-01 -2.0188E-01  2.8297E-01  4.2887E-01 -2.2853E-01
             2.0919E-01
 GRADIENT:   4.9711E+01  6.8815E+01  2.6665E+00  1.6965E+01 -5.8969E+00  5.5262E+00  9.3406E-01 -7.3711E-03  1.0802E+00  2.0217E-01
             1.4412E-01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1690.37332586315        NO. OF FUNC. EVALS.:  90
 CUMULATIVE NO. OF FUNC. EVALS.:     1082
 NPARAMETR:  1.0321E+00  1.7083E+00  5.2780E-01  6.0203E-01  1.0393E+00  1.0023E+00  7.3953E-01  1.2008E+00  1.3894E+00  7.1981E-01
             1.1154E+00
 PARAMETER:  1.3159E-01  6.3547E-01 -5.3904E-01 -4.0746E-01  1.3852E-01  1.0229E-01 -2.0188E-01  2.8297E-01  4.2887E-01 -2.2853E-01
             2.0919E-01
 GRADIENT:   4.2174E+00  7.2544E+00  5.9809E+01 -2.8053E+01  3.5219E+00 -6.6498E-03 -1.4530E-01  2.5475E+05 -6.7393E+01  1.2814E-01
            -3.4475E+05
 NUMSIGDIG:         8.7         7.8         7.0         7.4         8.8         3.8         2.1         3.3         7.0         1.9
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1082
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.9540E-04 -3.1635E-02 -2.3565E-02  1.7612E-02 -3.6364E-02
 SE:             2.9843E-02  2.3291E-02  1.1247E-02  2.3181E-02  1.8909E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9210E-01  1.7439E-01  3.6157E-02  4.4740E-01  5.4462E-02

 ETASHRINKSD(%)  2.1635E-02  2.1971E+01  6.2320E+01  2.2341E+01  3.6653E+01
 ETASHRINKVR(%)  4.3266E-02  3.9116E+01  8.5802E+01  3.9692E+01  5.9872E+01
 EBVSHRINKSD(%)  5.1720E-01  2.2078E+01  6.3990E+01  2.2726E+01  3.6218E+01
 EBVSHRINKVR(%)  1.0317E+00  3.9281E+01  8.7033E+01  4.0287E+01  5.9319E+01
 RELATIVEINF(%)  9.8935E+01  2.9353E+00  7.7261E-01  3.0375E+00  7.5341E+00
 EPSSHRINKSD(%)  4.3863E+01
 EPSSHRINKVR(%)  6.8486E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1690.3733258631473     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -955.22249929940915     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.75
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
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
 





 #OBJV:********************************************    -1690.373       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.71E+00  5.28E-01  6.02E-01  1.04E+00  1.00E+00  7.39E-01  1.20E+00  1.39E+00  7.20E-01  1.12E+00
 


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
+        1.95E+09
 
 TH 2
+       -7.81E+00  3.06E+07
 
 TH 3
+        5.01E+00  1.84E+02  4.45E+08
 
 TH 4
+        5.41E+08  3.24E+02  2.58E+08  5.99E+08
 
 TH 5
+       -5.30E+00 -2.86E+02 -4.40E+08 -5.11E+08  1.74E+09
 
 TH 6
+       -3.70E+00 -1.13E+00  1.00E+00 -3.55E+00  3.96E-02  1.95E+02
 
 TH 7
+        1.19E+00  4.49E+00 -7.72E+00 -1.09E+01 -9.53E+00  1.35E+00  1.62E+09
 
 TH 8
+       -6.98E-01 -1.13E+01 -8.73E+01  4.02E+01  1.66E+00  1.29E+03  6.07E+03  1.56E+08
 
 TH 9
+        2.39E+00 -1.72E+01 -3.91E+01  5.65E+01 -1.26E+01  5.22E-01  1.41E+01  6.71E+01  1.02E+08
 
 TH10
+        1.60E+00 -7.49E+00 -1.59E+01 -1.46E+01 -9.14E+01  1.23E+00 -7.34E+08  1.90E+03  4.94E+00  6.66E+08
 
 TH11
+       -8.50E+00 -1.39E+01 -1.52E+01 -5.42E+00 -2.18E+01 -1.87E+03  5.18E+08 -6.45E+03  8.41E+00 -4.70E+08  3.31E+08
 
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
 #CPUT: Total CPU Time in Seconds,       18.792
Stop Time:
Sat Sep 25 11:28:48 CDT 2021
