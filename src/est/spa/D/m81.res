Sat Sep 18 15:39:28 CDT 2021
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
$DATA ../../../../data/spa/D/dat81.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   23321.0151035634        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2743E+02  3.5971E+02 -5.7475E+01  4.5114E+02  4.1022E+01 -2.3242E+03 -1.0795E+03 -2.6039E+01 -1.3426E+03 -2.4924E+02
            -4.4656E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -526.735953180546        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1776E+00  1.0771E+00  9.6682E-01  1.2815E+00  1.2232E+00  1.5808E+00  1.0907E+00  9.7377E-01  9.9483E-01  9.6360E-01
             1.4911E+01
 PARAMETER:  2.6351E-01  1.7425E-01  6.6259E-02  3.4802E-01  3.0149E-01  5.5793E-01  1.8684E-01  7.3415E-02  9.4815E-02  6.2917E-02
             2.8021E+00
 GRADIENT:  -6.1027E+01  8.8753E+00 -1.5066E+00  7.3766E+00 -8.2483E+00  7.3826E+00  7.7520E-01  2.5541E+00  7.1627E+00  2.5255E+00
             4.1320E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -534.311889613102        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.2612E+00  9.5027E-01  1.0813E+00  1.4318E+00  2.0941E+00  1.4911E+00  1.6956E+00  2.8853E-01  6.8432E-01  1.1306E+00
             1.4977E+01
 PARAMETER:  3.3209E-01  4.8996E-02  1.7820E-01  4.5893E-01  8.3911E-01  4.9949E-01  6.2806E-01 -1.1430E+00 -2.7933E-01  2.2277E-01
             2.8065E+00
 GRADIENT:  -9.4246E+00  1.3743E+01 -1.2795E+00  1.1884E+01 -4.7804E+00 -7.1125E+00  2.8378E+00  1.9641E-01  3.8556E-01 -6.4831E-02
             1.2523E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -548.935982027533        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1588E+00  2.7569E-01  8.6132E-01  1.7291E+00  7.9914E+00  1.5961E+00  1.0294E+00  1.0000E-02  3.8069E-01  7.1650E+00
             1.4558E+01
 PARAMETER:  2.4741E-01 -1.1885E+00 -4.9285E-02  6.4758E-01  2.1784E+00  5.6756E-01  1.2900E-01 -6.4718E+00 -8.6577E-01  2.0692E+00
             2.7781E+00
 GRADIENT:  -2.0972E+01  9.4978E+00  1.1122E+01  8.4695E+01  3.1338E+00  2.1747E+00 -6.7636E-02  0.0000E+00 -1.2923E+00 -5.8221E+00
            -3.1662E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -549.791305904057        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  1.1366E+00  2.3297E-01  7.8456E-01  1.7412E+00  1.0028E+01  1.6047E+00  9.4594E-01  1.0000E-02  3.3200E-01  1.0813E+01
             1.4506E+01
 PARAMETER:  2.2808E-01 -1.3569E+00 -1.4264E-01  6.5460E-01  2.4054E+00  5.7297E-01  4.4425E-02 -7.4210E+00 -1.0026E+00  2.4808E+00
             2.7745E+00
 GRADIENT:  -4.3782E+01  4.0276E+00  1.1733E+01  3.4620E+01  1.4708E+01  3.3590E+00  6.4788E-02  0.0000E+00  1.9638E+00 -1.8384E+01
             1.5443E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -551.007244115141        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      388
 NPARAMETR:  1.1330E+00  2.1524E-01  7.5206E-01  1.7358E+00  1.0742E+01  1.5926E+00  9.1797E-01  1.0000E-02  3.0750E-01  1.1681E+01
             1.4510E+01
 PARAMETER:  2.2486E-01 -1.4360E+00 -1.8494E-01  6.5148E-01  2.4741E+00  5.6536E-01  1.4414E-02 -7.7938E+00 -1.0793E+00  2.5580E+00
             2.7749E+00
 GRADIENT:  -3.5991E+01  2.8299E+00  4.5557E+00  2.5077E+01  9.4455E+00 -8.5336E+00  8.0802E-02  0.0000E+00  2.2527E+00  3.8603E+00
             1.7259E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -555.241321144669        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      458
 NPARAMETR:  1.1618E+00  9.4886E-02  6.3318E-01  1.6634E+00  8.8985E+00  1.4972E+00  1.0483E+00  1.0000E-02  1.3587E-01  9.3487E+00
             1.4398E+01
 PARAMETER:  2.5001E-01 -2.2551E+00 -3.5700E-01  6.0888E-01  2.2859E+00  5.0361E-01  1.4715E-01 -9.6949E+00 -1.8960E+00  2.3352E+00
             2.7671E+00
 GRADIENT:   3.3362E+00  1.0811E+00  3.3980E+00  1.2774E+01  1.6087E+01 -1.2909E+01  2.1365E-02  0.0000E+00  5.2934E-01 -1.5556E+01
             1.6126E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -555.462692838113        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:      539
 NPARAMETR:  1.1601E+00  8.9278E-02  6.3019E-01  1.6603E+00  8.8125E+00  1.5007E+00  1.0563E+00  1.0000E-02  1.2794E-01  9.2948E+00
             1.4383E+01
 PARAMETER:  2.4850E-01 -2.3160E+00 -3.6174E-01  6.0697E-01  2.2762E+00  5.0595E-01  1.5482E-01 -9.8352E+00 -1.9562E+00  2.3295E+00
             2.7660E+00
 GRADIENT:   9.7258E+00  2.8040E+00 -3.9786E+00  2.6536E+01 -8.1347E+00 -1.8054E+01  1.9405E-02  0.0000E+00  4.4294E-01  8.7388E-01
            -2.0276E+01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -555.544622928340        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  1.1592E+00  8.7221E-02  6.2881E-01  1.6587E+00  8.7930E+00  1.5030E+00  1.0592E+00  1.0000E-02  1.2504E-01  9.2725E+00
             1.4387E+01
 PARAMETER:  2.4783E-01 -2.3387E+00 -3.6401E-01  6.0618E-01  2.2734E+00  5.0734E-01  1.5771E-01 -9.8914E+00 -1.9796E+00  2.3276E+00
             2.7656E+00
 GRADIENT:   1.1040E+03  1.1777E+02 -7.4694E+02  4.6106E+02 -2.2348E+02 -5.5037E+02  1.1485E-01  0.0000E+00 -2.7452E+02  2.1878E+02
            -9.9173E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      648
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0751E-03 -1.0670E-03 -9.4741E-05 -5.8829E-03 -4.8732E-02
 SE:             2.8537E-02  4.7341E-04  6.0924E-05  3.0729E-03  8.7571E-03
 N:                     100         100         100         100         100

 P VAL.:         9.6995E-01  2.4200E-02  1.1993E-01  5.5560E-02  2.6296E-08

 ETASHRINKSD(%)  4.3986E+00  9.8414E+01  9.9796E+01  8.9706E+01  7.0663E+01
 ETASHRINKVR(%)  8.6038E+00  9.9975E+01  1.0000E+02  9.8940E+01  9.1393E+01
 EBVSHRINKSD(%)  7.1222E+00  9.8145E+01  9.9747E+01  8.9351E+01  6.0390E+01
 EBVSHRINKVR(%)  1.3737E+01  9.9966E+01  9.9999E+01  9.8866E+01  8.4310E+01
 RELATIVEINF(%)  3.5739E+01  1.2488E-03  8.0016E-05  2.4835E-02  6.9003E+00
 EPSSHRINKSD(%)  2.0286E+00
 EPSSHRINKVR(%)  4.0160E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -555.54462292833955     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       179.60620363539863     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -555.545       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.16E+00  8.73E-02  6.29E-01  1.66E+00  8.79E+00  1.50E+00  1.06E+00  1.00E-02  1.25E-01  9.28E+00  1.44E+01
 


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
+        1.65E+06
 
 TH 2
+       -1.69E+03  3.26E+06
 
 TH 3
+        1.77E+03  1.02E+03  2.59E+06
 
 TH 4
+       -3.52E+02 -7.59E+02 -3.73E+01  1.34E+05
 
 TH 5
+        8.21E+00  1.63E+01  1.37E+01  1.52E+01  3.38E+02
 
 TH 6
+        5.12E+02  1.43E+03 -5.14E+02  1.74E+02 -6.15E+00  2.34E+05
 
 TH 7
+       -2.83E+00 -1.53E+00  3.12E+00  4.37E-01 -9.97E-02 -4.33E-01 -8.88E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.79E+00 -2.31E+01  6.23E+00 -2.57E+00  5.12E-01  4.65E+00  3.92E+00  0.00E+00  2.22E+06
 
 TH10
+       -1.26E+01 -1.98E+01 -8.70E+00 -2.10E+01  1.40E+00  7.65E+00 -7.25E-02  0.00E+00 -7.90E-03  2.89E+02
 
 TH11
+       -2.41E+00  8.77E+00  3.17E+00  6.87E+00 -1.72E+00 -2.04E+00 -9.33E-03  0.00E+00  4.65E-01  1.47E+00  8.58E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.254
Stop Time:
Sat Sep 18 15:39:52 CDT 2021
