Sat Sep 18 13:26:44 CDT 2021
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
$DATA ../../../../data/spa/S2/dat42.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1714.20467440892        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -6.2469E+01  2.2859E+01 -2.6186E+01  6.9403E+01  6.4694E+01  1.2211E+00  9.1420E+00  1.4611E+00  5.1769E+00 -5.6114E+00
             8.0564E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1720.47296714680        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  1.0603E+00  9.3071E-01  9.6320E-01  9.9983E-01  9.0171E-01  9.9107E-01  8.9056E-01  9.9196E-01  9.6009E-01  9.7857E-01
             9.5000E-01
 PARAMETER:  1.5856E-01  2.8189E-02  6.2502E-02  9.9826E-02 -3.4613E-03  9.1033E-02 -1.5905E-02  9.1932E-02  5.9272E-02  7.8341E-02
             4.8705E-02
 GRADIENT:   1.7147E+01  3.7960E+00  6.8643E+00 -1.2225E+01 -1.6234E+01 -2.3851E+00  2.5197E+00  2.7451E+00 -5.4540E+00  6.4388E+00
            -1.2257E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1721.24166817996        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  1.0627E+00  8.6191E-01  9.7511E-01  1.0431E+00  8.8219E-01  9.9350E-01  7.9285E-01  8.6327E-01  9.7905E-01  9.6249E-01
             9.7410E-01
 PARAMETER:  1.6085E-01 -4.8600E-02  7.4792E-02  1.4217E-01 -2.5347E-02  9.3479E-02 -1.3212E-01 -4.7023E-02  7.8823E-02  6.1771E-02
             7.3755E-02
 GRADIENT:   2.2937E+01  4.2563E+00  1.0371E+01 -5.7298E+00 -1.2731E+01 -1.2383E+00  6.2660E-01 -1.7326E+00  8.2959E-01  7.7084E-01
            -3.2417E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1721.53072181968        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  1.0514E+00  8.0991E-01  9.8404E-01  1.0810E+00  8.7014E-01  9.9641E-01  8.0589E-01  8.9386E-01  9.5664E-01  9.4293E-01
             9.8296E-01
 PARAMETER:  1.5016E-01 -1.1083E-01  8.3915E-02  1.7790E-01 -3.9107E-02  9.6406E-02 -1.1581E-01 -1.2208E-02  5.5671E-02  4.1236E-02
             8.2814E-02
 GRADIENT:  -3.1710E-01  4.6031E+00  1.7793E+00  6.2747E+00 -2.6416E+00  3.4041E-01  5.3427E-01 -1.9875E-01  1.0417E+00 -2.3843E-01
             1.0596E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1721.77595747038        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  1.0469E+00  5.9053E-01  1.1357E+00  1.2170E+00  8.6118E-01  9.9196E-01  5.4213E-01  9.9198E-01  8.9704E-01  9.7472E-01
             9.8151E-01
 PARAMETER:  1.4583E-01 -4.2674E-01  2.2721E-01  2.9635E-01 -4.9453E-02  9.1930E-02 -5.1226E-01  9.1945E-02 -8.6578E-03  7.4396E-02
             8.1339E-02
 GRADIENT:  -2.5118E+00  1.0918E+00  1.3119E+00  2.8576E+00 -1.6578E+00 -1.3225E-01  1.6705E-01 -2.1562E-01  1.0509E+00  7.0070E-02
             7.2279E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1721.79908453244        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  1.0471E+00  5.3970E-01  1.1535E+00  1.2457E+00  8.5400E-01  9.9156E-01  4.3321E-01  1.0040E+00  8.7944E-01  9.7488E-01
             9.7998E-01
 PARAMETER:  1.4600E-01 -5.1674E-01  2.4281E-01  3.1973E-01 -5.7829E-02  9.1524E-02 -7.3654E-01  1.0400E-01 -2.8473E-02  7.4554E-02
             7.9777E-02
 GRADIENT:  -3.4325E-01 -9.0730E-01 -3.9498E-01 -1.8284E+00  8.1448E-01 -6.7866E-03  9.5301E-02  5.4992E-02  2.6609E-01  2.0763E-02
             1.6755E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1721.83185859545        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  1.0473E+00  5.4963E-01  1.1647E+00  1.2400E+00  8.6152E-01  9.9176E-01  1.5663E-01  1.0191E+00  8.9276E-01  9.8325E-01
             9.7952E-01
 PARAMETER:  1.4624E-01 -4.9851E-01  2.5243E-01  3.1515E-01 -4.9053E-02  9.1724E-02 -1.7539E+00  1.1892E-01 -1.3440E-02  8.3108E-02
             7.9303E-02
 GRADIENT:  -1.4361E-02 -7.7281E-02  6.3419E-02 -7.9281E-02 -2.0786E-02  7.7967E-03  8.7083E-03 -5.9919E-03  3.5186E-01  6.1040E-02
             2.0654E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1721.83582541645        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  1.0472E+00  5.4223E-01  1.1673E+00  1.2446E+00  8.6026E-01  9.9163E-01  2.3644E-02  1.0210E+00  8.8976E-01  9.8275E-01
             9.7949E-01
 PARAMETER:  1.4614E-01 -5.1206E-01  2.5472E-01  3.1881E-01 -5.0519E-02  9.1590E-02 -3.6446E+00  1.2078E-01 -1.6799E-02  8.2595E-02
             7.9281E-02
 GRADIENT:   6.4632E-04  3.9649E-02  1.9126E-02  8.9707E-02 -7.2898E-02 -5.8729E-03  1.8506E-04  1.0229E-02 -1.7968E-02  1.1522E-03
            -4.6269E-03

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1721.83590898905        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1327
 NPARAMETR:  1.0472E+00  5.4254E-01  1.1672E+00  1.2444E+00  8.6038E-01  9.9165E-01  1.0000E-02  1.0207E+00  8.8995E-01  9.8286E-01
             9.7951E-01
 PARAMETER:  1.4614E-01 -5.1149E-01  2.5463E-01  3.1862E-01 -5.0382E-02  9.1611E-02 -4.7035E+00  1.2048E-01 -1.6590E-02  8.2708E-02
             7.9300E-02
 GRADIENT:  -5.5538E-03  3.6973E-03 -1.0775E-02  6.4260E-03  1.8997E-02  6.3028E-04  0.0000E+00 -9.3612E-04 -1.3180E-02 -6.2656E-03
             1.5538E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1327
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.4624E-04 -2.8696E-04 -3.0455E-02 -3.8804E-03 -2.9259E-02
 SE:             2.9874E-02  1.1655E-04  1.5997E-02  2.9235E-02  2.2802E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8541E-01  1.3811E-02  5.6937E-02  8.9441E-01  1.9943E-01

 ETASHRINKSD(%)  1.0000E-10  9.9610E+01  4.6409E+01  2.0594E+00  2.3609E+01
 ETASHRINKVR(%)  1.0000E-10  9.9998E+01  7.1280E+01  4.0764E+00  4.1645E+01
 EBVSHRINKSD(%)  4.0525E-01  9.9650E+01  5.0312E+01  2.3834E+00  2.0825E+01
 EBVSHRINKVR(%)  8.0887E-01  9.9999E+01  7.5311E+01  4.7100E+00  3.7313E+01
 RELATIVEINF(%)  9.8248E+01  8.3186E-05  5.2967E+00  8.2280E+00  9.9631E+00
 EPSSHRINKSD(%)  4.5374E+01
 EPSSHRINKVR(%)  7.0160E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1721.8359089890519     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -986.68508242531368     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.06
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1721.836       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  5.43E-01  1.17E+00  1.24E+00  8.60E-01  9.92E-01  1.00E-02  1.02E+00  8.90E-01  9.83E-01  9.80E-01
 


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
+        1.02E+03
 
 TH 2
+       -2.49E+01  4.95E+02
 
 TH 3
+        5.14E+00  1.10E+02  2.06E+02
 
 TH 4
+       -1.14E+01  5.43E+02 -2.94E+01  8.64E+02
 
 TH 5
+        3.95E+00 -3.25E+02 -3.68E+02 -5.18E+01  9.80E+02
 
 TH 6
+       -6.19E-01 -4.15E+00  1.19E+00 -2.92E+00  4.01E+00  1.98E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -2.22E-01 -1.42E+01 -3.83E+01 -3.62E+00  9.32E+00 -2.07E+00  0.00E+00  2.75E+01
 
 TH 9
+        3.00E+00 -1.06E+02  4.78E+00  2.65E-01 -5.60E+00  4.43E+00  0.00E+00  4.79E-01  2.24E+02
 
 TH10
+        1.36E+00  2.83E+00 -1.36E+01 -1.16E+00 -7.92E+01 -8.42E-01  0.00E+00  2.05E+01 -4.25E-01  8.47E+01
 
 TH11
+       -7.87E+00 -1.37E+01 -1.23E+01 -8.07E+00 -5.53E-02  3.44E+00  0.00E+00  1.17E+01  9.18E+00  1.52E+01  2.13E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.678
Stop Time:
Sat Sep 18 13:27:07 CDT 2021
