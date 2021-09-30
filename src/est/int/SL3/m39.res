Wed Sep 29 04:15:34 CDT 2021
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
$DATA ../../../../data/int/SL3/dat39.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      979
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E19.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      879
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -633.816225554604        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5651E+02  7.3047E+01  1.5479E+02  9.4619E+01  2.4792E+01  5.5556E+01 -6.4252E+01 -2.4403E+02 -7.3249E+01 -4.5286E+01
            -5.8297E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2700.14657459699        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0783E+00  1.3738E+00  1.0453E+00  8.8412E-01  1.3172E+00  8.5629E-01  1.0057E+00  8.9771E-01  7.6580E-01  1.1449E+00
             2.6528E+00
 PARAMETER:  1.7534E-01  4.1755E-01  1.4429E-01 -2.3162E-02  3.7550E-01 -5.5150E-02  1.0566E-01 -7.9065E-03 -1.6683E-01  2.3535E-01
             1.0756E+00
 GRADIENT:   3.0743E+02  1.0678E+02 -2.0591E+01  9.9750E+01  2.6691E+01 -4.7874E+01  1.7640E+01  3.7208E+00 -9.6297E+00 -2.9196E+01
             1.5026E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2709.78418654444        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      222
 NPARAMETR:  1.0789E+00  1.5976E+00  1.4587E+00  8.2181E-01  1.6788E+00  9.6193E-01  8.7749E-01  7.8145E-01  9.9102E-01  1.4746E+00
             2.6569E+00
 PARAMETER:  1.7591E-01  5.6848E-01  4.7754E-01 -9.6252E-02  6.1806E-01  6.1190E-02 -3.0685E-02 -1.4661E-01  9.0975E-02  4.8836E-01
             1.0771E+00
 GRADIENT:   1.5843E+02  1.2725E+02 -1.8829E+01  1.6391E+02  3.4255E+01 -5.3372E+00  1.1241E+01  2.0843E-01 -3.3369E+00 -8.9571E+00
             2.7734E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2742.78408321218        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      400
 NPARAMETR:  9.8553E-01  1.4828E+00  4.1265E+00  7.5914E-01  2.0804E+00  9.1776E-01  8.7513E-01  2.8530E+00  7.1681E-01  1.7907E+00
             2.5631E+00
 PARAMETER:  8.5422E-02  4.9393E-01  1.5174E+00 -1.7557E-01  8.3256E-01  1.4178E-02 -3.3387E-02  1.1484E+00 -2.3294E-01  6.8258E-01
             1.0412E+00
 GRADIENT:  -4.6903E+01 -7.5340E+01 -1.1990E+01 -2.8264E+01  2.0605E+01 -1.3124E+01 -3.3446E+00 -1.9260E+00  1.9489E+00 -1.2239E+00
             1.2051E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2752.06502607918        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  1.0133E+00  1.5165E+00  2.2843E+01  7.5512E-01  2.3194E+00  9.3677E-01  9.9644E-01  6.7460E+00  3.2245E-01  2.0757E+00
             2.5399E+00
 PARAMETER:  1.1318E-01  5.1643E-01  3.2286E+00 -1.8088E-01  9.4131E-01  3.4679E-02  9.6432E-02  2.0090E+00 -1.0318E+00  8.3029E-01
             1.0321E+00
 GRADIENT:   2.1761E+01 -3.2520E+01  7.3926E-01 -2.6919E+01  3.3190E+00 -5.3219E+00  5.2634E+00 -2.1536E+00  5.9603E-01  1.8674E+01
             2.2092E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2753.30219014196        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      756
 NPARAMETR:  1.0071E+00  1.5246E+00  2.8820E+01  7.6651E-01  2.3157E+00  9.4686E-01  9.7690E-01  7.4727E+00  3.0665E-01  1.9588E+00
             2.5379E+00
 PARAMETER:  1.0711E-01  5.2171E-01  3.4611E+00 -1.6591E-01  9.3972E-01  4.5398E-02  7.6630E-02  2.1113E+00 -1.0821E+00  7.7231E-01
             1.0313E+00
 GRADIENT:   5.6164E+00 -3.4220E+00 -1.4809E-01 -3.4886E-01  3.2917E+00 -1.4439E+00 -2.2927E+00 -1.2317E+00 -2.0601E-01  3.4487E+00
            -2.7187E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2753.44230575684        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      935
 NPARAMETR:  1.0047E+00  1.5206E+00  3.3621E+01  7.7165E-01  2.3027E+00  9.5182E-01  9.8825E-01  8.0157E+00  2.8547E-01  1.9156E+00
             2.5424E+00
 PARAMETER:  1.0465E-01  5.1911E-01  3.6151E+00 -1.5922E-01  9.3410E-01  5.0617E-02  8.8176E-02  2.1814E+00 -1.1536E+00  7.5004E-01
             1.0331E+00
 GRADIENT:  -5.9036E-01  1.1256E+00 -3.8802E+00  7.7590E+00 -2.0024E+00  5.1618E-01 -7.6915E-01  5.5679E+00 -3.3023E-01 -5.0051E+00
            -2.3261E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2753.46440722118        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1112
 NPARAMETR:  1.0055E+00  1.5405E+00  3.3859E+01  7.5765E-01  2.3019E+00  9.4994E-01  9.7552E-01  8.1725E+00  2.9663E-01  1.9131E+00
             2.5439E+00
 PARAMETER:  1.0548E-01  5.3212E-01  3.6222E+00 -1.7753E-01  9.3371E-01  4.8640E-02  7.5211E-02  2.2008E+00 -1.1153E+00  7.4873E-01
             1.0337E+00
 GRADIENT:   1.4857E+00 -1.5697E-01 -4.6013E-01  9.1250E-01 -3.4136E-01 -3.2304E-01 -2.9758E-01  8.4661E-01 -1.7111E-01  9.4466E-01
            -2.2672E+00

0ITERATION NO.:   39    OBJECTIVE VALUE:  -2753.47727054412        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:     1267
 NPARAMETR:  1.0048E+00  1.5774E+00  3.3860E+01  7.3322E-01  2.3050E+00  9.5060E-01  9.5163E-01  8.3813E+00  3.2841E-01  1.9081E+00
             2.5456E+00
 PARAMETER:  1.0476E-01  5.5578E-01  3.6199E+00 -2.1032E-01  9.3460E-01  4.9313E-02  5.0421E-02  2.2274E+00 -1.0237E+00  7.4656E-01
             1.0341E+00
 GRADIENT:  -5.6010E-01  7.8704E-02 -1.8644E+01 -1.0266E-01 -1.1125E+01 -7.1533E-02 -6.3438E-03  3.0614E+01 -7.9539E-02  1.5255E+01
            -1.1251E+01
 NUMSIGDIG:         2.6         4.2         2.3         3.4         2.4         2.7         3.7         2.4         1.2         2.4
                    2.7

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1267
 NO. OF SIG. DIGITS IN FINAL EST.:  1.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6298E-03 -4.3182E-03 -2.2698E-02 -1.2878E-02 -2.2651E-02
 SE:             2.9374E-02  2.7260E-02  9.7816E-03  7.0446E-03  2.6659E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5575E-01  8.7414E-01  2.0312E-02  6.7539E-02  3.9552E-01

 ETASHRINKSD(%)  1.5938E+00  8.6753E+00  6.7230E+01  7.6400E+01  1.0690E+01
 ETASHRINKVR(%)  3.1621E+00  1.6598E+01  8.9262E+01  9.4430E+01  2.0237E+01
 EBVSHRINKSD(%)  1.6105E+00  8.7037E+00  7.7388E+01  7.8113E+01  9.1051E+00
 EBVSHRINKVR(%)  3.1950E+00  1.6650E+01  9.4887E+01  9.5210E+01  1.7381E+01
 RELATIVEINF(%)  9.6725E+01  5.1075E+00  4.3603E+00  2.9135E-01  7.4669E+01
 EPSSHRINKSD(%)  1.6383E+01
 EPSSHRINKVR(%)  3.0081E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          879
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1615.4939413738146     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2753.4772705441178     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1137.9833291703033     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.50
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2753.477       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.58E+00  3.38E+01  7.33E-01  2.30E+00  9.51E-01  9.52E-01  8.39E+00  3.25E-01  1.91E+00  2.54E+00
 


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
+        1.18E+03
 
 TH 2
+       -1.92E+01  4.35E+02
 
 TH 3
+       -4.05E-02  3.55E-01  2.12E-01
 
 TH 4
+       -1.92E+01  5.14E+02 -3.46E-01  1.27E+03
 
 TH 5
+       -1.08E+01  1.55E+02  1.63E-01 -8.41E+02  1.32E+02
 
 TH 6
+        6.12E+00 -2.38E+00  4.94E-02 -1.87E+01  9.58E+00  2.07E+02
 
 TH 7
+        6.82E+00 -5.67E+00 -2.90E-01 -3.39E+01 -5.80E+01 -1.83E+00  1.71E+02
 
 TH 8
+        8.95E-02 -1.73E+00 -1.25E+00  5.37E-01  1.79E-01 -1.86E-01  1.11E+00  7.96E+00
 
 TH 9
+        3.30E-01 -2.47E+00  1.91E-02 -7.17E+01  1.39E+01  1.62E+00  1.98E+01  1.25E-01  1.04E+01
 
 TH10
+        1.05E+01 -2.54E+02 -1.18E+00  1.25E+03  3.38E+02 -1.22E+01  6.30E+01  5.41E+00 -1.17E+01  9.45E+02
 
 TH11
+       -1.77E+01 -8.18E+00  1.45E+00 -3.48E+01  1.28E+01  3.49E+00  4.79E-01 -5.99E+00  2.69E+00 -1.46E+01  2.55E+02
 
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
 #CPUT: Total CPU Time in Seconds,       50.203
Stop Time:
Wed Sep 29 04:16:26 CDT 2021
