Wed Sep 29 03:57:51 CDT 2021
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
$DATA ../../../../data/int/SL3/dat14.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      981
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      881
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
 RAW OUTPUT FILE (FILE): m14.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   234.508038607567        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1781E+02 -5.3036E+01  9.4368E+01  2.0487E+02  2.5633E+02  3.3730E+01 -1.4876E+02 -2.4299E+02 -1.0231E+02 -3.0207E+01
            -7.6033E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2366.54769636943        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0127E+00  1.4696E+00  8.6525E-01  9.0737E-01  1.0311E+00  9.3743E-01  1.1455E+00  1.0539E+00  7.8188E-01  9.6213E-01
             5.2161E+00
 PARAMETER:  1.1265E-01  4.8496E-01 -4.4736E-02  2.7947E-03  1.3059E-01  3.5384E-02  2.3585E-01  1.5254E-01 -1.4605E-01  6.1395E-02
             1.7517E+00
 GRADIENT:  -1.3665E+02  8.6451E+01 -9.5096E+00  4.9384E+01 -1.7043E+01 -2.9874E+01  3.0554E+01  5.7989E+00  1.3303E+01  9.0429E+00
             8.1051E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2602.89071542049        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.5592E-01  8.0399E-01  5.0317E+00  1.2336E+00  1.1673E+00  1.1813E+00  2.1735E+00  3.7553E+00  1.9290E-01  4.2106E-01
             3.2339E+00
 PARAMETER:  5.4924E-02 -1.1817E-01  1.7158E+00  3.0994E-01  2.5469E-01  2.6660E-01  8.7634E-01  1.4232E+00 -1.5456E+00 -7.6498E-01
             1.2737E+00
 GRADIENT:  -8.1329E+01  1.8896E+01  3.6205E+01  2.3775E+01 -7.0731E+01  5.4368E+01  4.2333E+01 -1.4091E+01 -9.9769E+00 -1.1422E+01
             2.1105E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2641.44401939334        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.9495E-01  1.0254E+00  2.3095E+00  1.1039E+00  1.1075E+00  1.0501E+00  1.2724E+00  3.3013E+00  5.9824E-01  8.0538E-01
             2.8342E+00
 PARAMETER:  9.4942E-02  1.2510E-01  9.3703E-01  1.9887E-01  2.0215E-01  1.4893E-01  3.4093E-01  1.2943E+00 -4.1376E-01 -1.1645E-01
             1.1418E+00
 GRADIENT:  -1.6554E+01  1.9479E+01 -8.6475E-01  6.3806E+01 -9.9422E+00  1.6451E+01  1.3612E+00 -2.0232E+00  2.1477E-01 -4.1316E+00
            -1.3537E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2655.07671905879        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0031E+00  1.1123E+00  3.3017E+00  1.0354E+00  1.2696E+00  1.0105E+00  1.1616E+00  3.7019E+00  5.7882E-01  1.0520E+00
             2.8342E+00
 PARAMETER:  1.0314E-01  2.0641E-01  1.2944E+00  1.3475E-01  3.3872E-01  1.1043E-01  2.4978E-01  1.4088E+00 -4.4676E-01  1.5068E-01
             1.1418E+00
 GRADIENT:  -2.7135E+00 -1.6775E+00  1.0518E+00  8.3606E-01 -4.4759E-01  5.7344E-01 -8.5694E-01  5.8543E-02  2.0432E-01 -4.9321E-01
             2.8419E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2656.43803582310        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      462
 NPARAMETR:  1.0456E+00  1.1465E+00  3.4037E+00  1.0322E+00  1.2835E+00  1.0364E+00  1.1591E+00  3.7114E+00  6.0908E-01  1.1525E+00
             2.8557E+00
 PARAMETER:  1.4454E-01  2.3667E-01  1.3249E+00  1.3168E-01  3.4959E-01  1.3580E-01  2.4768E-01  1.4114E+00 -3.9580E-01  2.4191E-01
             1.1493E+00
 GRADIENT:   1.0628E+02  2.3083E+01 -4.6976E+00  2.6928E+01 -2.2527E+01  1.2608E+01  5.4768E+00  4.8516E+00  3.0064E+00  4.7456E+00
             1.3781E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2656.67361140763        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:      596
 NPARAMETR:  1.0332E+00  1.1417E+00  3.4053E+00  1.0322E+00  1.2837E+00  1.0241E+00  1.1592E+00  3.7103E+00  5.9803E-01  1.1187E+00
             2.8564E+00
 PARAMETER:  1.3263E-01  2.3255E-01  1.3253E+00  1.3165E-01  3.4971E-01  1.2382E-01  2.4772E-01  1.4111E+00 -4.1412E-01  2.1217E-01
             1.1496E+00
 GRADIENT:  -3.7207E-01 -1.8707E-01 -1.0418E+01  4.3398E+00 -3.6019E+01  3.1403E-02  1.5529E+00 -2.1147E+00  1.3163E-01 -2.0309E-01
            -6.3138E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2662.40652182523        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      759             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0337E+00  1.1495E+00  4.5618E+00  1.0319E+00  1.3896E+00  1.0256E+00  1.1346E+00  3.9703E+00  6.0146E-01  1.1945E+00
             2.8661E+00
 PARAMETER:  1.3317E-01  2.3936E-01  1.6177E+00  1.3144E-01  4.2903E-01  1.2523E-01  2.2631E-01  1.4789E+00 -4.0839E-01  2.7774E-01
             1.1530E+00
 GRADIENT:   7.5064E+01  1.3028E+01  2.6940E+00  1.5228E+01  1.1840E+01  8.1927E+00  3.1423E+00  4.0523E+00  2.0735E+00  1.5200E+00
             2.9193E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2662.50814808314        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      936
 NPARAMETR:  1.0341E+00  1.1530E+00  4.6761E+00  1.0333E+00  1.3944E+00  1.0248E+00  1.1278E+00  3.9559E+00  6.0290E-01  1.1946E+00
             2.8604E+00
 PARAMETER:  1.3354E-01  2.4236E-01  1.6425E+00  1.3278E-01  4.3248E-01  1.2449E-01  2.2025E-01  1.4752E+00 -4.0600E-01  2.7783E-01
             1.1509E+00
 GRADIENT:   3.5968E-01 -3.2167E+00 -1.8893E+00  2.4265E+00 -6.4738E+00 -1.1110E-01  1.9122E-01 -4.7274E+00 -1.0607E-01 -8.3111E-01
             3.7361E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2662.63494173613        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1120
 NPARAMETR:  1.0343E+00  1.1663E+00  4.6361E+00  1.0246E+00  1.4059E+00  1.0255E+00  1.1202E+00  3.9849E+00  6.0364E-01  1.2166E+00
             2.8556E+00
 PARAMETER:  1.3371E-01  2.5385E-01  1.6339E+00  1.2426E-01  4.4065E-01  1.2522E-01  2.1353E-01  1.4825E+00 -4.0478E-01  2.9605E-01
             1.1493E+00
 GRADIENT:   7.6131E-01 -4.1411E+00 -2.9694E+00  3.6647E+00 -2.8704E+00  1.1144E-01  4.5662E-01 -3.6473E+00 -3.5946E-02  3.9675E-01
             3.5531E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -2662.64740931716        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1186
 NPARAMETR:  1.0345E+00  1.1676E+00  4.6232E+00  1.0220E+00  1.4075E+00  1.0255E+00  1.1169E+00  3.9740E+00  6.0401E-01  1.2171E+00
             2.8610E+00
 PARAMETER:  1.3365E-01  2.5542E-01  1.6340E+00  1.2301E-01  4.4099E-01  1.2521E-01  2.1218E-01  1.4825E+00 -4.0469E-01  2.9705E-01
             1.1493E+00
 GRADIENT:  -7.3663E+02  3.8202E+02  5.5603E+01  2.9520E+00 -1.1514E+02  8.1250E-02  3.3773E-01  9.6279E+02 -6.9316E-02  1.6545E+02
            -9.0126E+01
 NUMSIGDIG:         2.3         2.3         2.3         1.6         2.3         2.8         1.7         2.3         2.5         2.3
                    2.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1186
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0586E-03 -3.4515E-03 -3.6425E-02 -5.4580E-03 -3.4573E-02
 SE:             2.9352E-02  2.3031E-02  1.9526E-02  1.6876E-02  2.1456E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4409E-01  8.8087E-01  6.2117E-02  7.4638E-01  1.0710E-01

 ETASHRINKSD(%)  1.6660E+00  2.2843E+01  3.4586E+01  4.3462E+01  2.8120E+01
 ETASHRINKVR(%)  3.3043E+00  4.0469E+01  5.7210E+01  6.8034E+01  4.8332E+01
 EBVSHRINKSD(%)  1.7415E+00  2.3262E+01  4.0050E+01  4.5115E+01  2.4216E+01
 EBVSHRINKVR(%)  3.4528E+00  4.1113E+01  6.4060E+01  6.9877E+01  4.2568E+01
 RELATIVEINF(%)  9.6457E+01  5.8688E+00  1.7353E+01  2.9813E+00  3.2128E+01
 EPSSHRINKSD(%)  1.6688E+01
 EPSSHRINKVR(%)  3.0591E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          881
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1619.1696955066332     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2662.6474093171560     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1043.4777138105228     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    33.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2662.647       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.17E+00  4.64E+00  1.02E+00  1.41E+00  1.03E+00  1.12E+00  3.98E+00  6.04E-01  1.22E+00  2.86E+00
 


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
+        1.30E+05
 
 TH 2
+       -5.98E+04  2.81E+04
 
 TH 3
+        1.39E-01  5.44E+00  4.18E+01
 
 TH 4
+       -7.52E+02  8.52E+02  5.48E+00  8.58E+02
 
 TH 5
+       -3.15E+00 -1.36E+04 -2.14E+01 -1.41E+02  6.78E+03
 
 TH 6
+       -6.69E-01 -2.72E+00  1.49E-01 -8.99E+00 -3.54E+00  1.75E+02
 
 TH 7
+       -2.61E+02  1.45E+02  4.86E+00 -2.05E+01 -5.82E+01  4.19E-02  6.05E+01
 
 TH 8
+        2.94E+03  9.05E+00  9.11E+00 -3.33E+03 -3.93E+01  1.77E+00  3.11E+01  1.07E+03
 
 TH 9
+       -3.42E+02  1.45E+02  5.98E+00 -3.77E+01 -7.08E+01 -1.73E-01  3.55E+01  3.83E+01  6.92E+01
 
 TH10
+       -1.03E+00  2.28E+04 -7.89E+00  2.97E+02 -1.11E+04  2.67E+00  1.03E+02  5.61E+01  1.33E+02  1.87E+04
 
 TH11
+       -1.38E+01 -2.23E+01 -1.10E+01 -4.98E+01  2.52E+01  2.03E+00 -7.12E+00 -1.51E+01 -8.39E+00  2.09E+01  3.91E+02
 
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
 #CPUT: Total CPU Time in Seconds,       48.630
Stop Time:
Wed Sep 29 03:58:42 CDT 2021
