Sat Sep 18 11:03:33 CDT 2021
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
$DATA ../../../../data/spa/S1/dat37.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1636.67121017797        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.7119E+01 -1.9884E+01  7.0650E+00 -1.8502E+01  1.2165E+01  1.0573E+01 -2.9259E+01 -6.2641E+00 -1.7889E+01 -1.5561E+01
            -1.5040E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1641.53131572075        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  1.0166E+00  1.1332E+00  9.1479E-01  9.3121E-01  1.0322E+00  9.6491E-01  1.5460E+00  1.0486E+00  1.0659E+00  1.1107E+00
             1.0220E+00
 PARAMETER:  1.1646E-01  2.2508E-01  1.0941E-02  2.8733E-02  1.3174E-01  6.4282E-02  5.3565E-01  1.4741E-01  1.6385E-01  2.0495E-01
             1.2171E-01
 GRADIENT:   6.8298E+01  1.8986E+01 -9.5689E+00  5.5791E+00  1.3312E+00 -4.4142E+00  3.0252E+01  5.3499E+00  1.9987E+01  9.6336E+00
             1.4922E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1644.65425922847        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0135E+00  7.8249E-01  1.0719E+00  1.1861E+00  9.2661E-01  9.6442E-01  2.0035E+00  1.0207E+00  8.5685E-01  1.0690E+00
             1.0019E+00
 PARAMETER:  1.1341E-01 -1.4527E-01  1.6944E-01  2.7069E-01  2.3782E-02  6.3772E-02  7.9489E-01  1.2049E-01 -5.4496E-02  1.6674E-01
             1.0190E-01
 GRADIENT:   6.8946E+01  3.3630E+01  3.5305E+00  7.5064E+01 -1.7479E+01 -2.9764E+00  1.7350E+01  1.8198E+00  1.1640E+01  8.4983E+00
            -8.5454E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1645.43106992483        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      238
 NPARAMETR:  1.0092E+00  8.0497E-01  1.0320E+00  1.1597E+00  9.2352E-01  9.6266E-01  1.9131E+00  9.5236E-01  8.6323E-01  1.0600E+00
             1.0058E+00
 PARAMETER:  1.0917E-01 -1.1695E-01  1.3153E-01  2.4814E-01  2.0437E-02  6.1948E-02  7.4875E-01  5.1186E-02 -4.7073E-02  1.5823E-01
             1.0576E-01
 GRADIENT:   5.5325E+01  2.7216E+01  3.0583E+00  5.9369E+01 -1.3034E+01 -3.7611E+00  1.4883E+01  5.4266E-01  1.0783E+01  8.0946E+00
            -6.8436E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1646.12002936043        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  1.0092E+00  8.0507E-01  1.0318E+00  1.1596E+00  9.2344E-01  9.7154E-01  1.9130E+00  9.9251E-01  8.0527E-01  1.0032E+00
             1.0273E+00
 PARAMETER:  1.0918E-01 -1.1682E-01  1.3134E-01  2.4809E-01  2.0355E-02  7.1130E-02  7.4869E-01  9.2477E-02 -1.1658E-01  1.0321E-01
             1.2696E-01
 GRADIENT:   5.2500E+01  2.9102E+01  1.5264E+00  6.2345E+01 -6.6594E+00 -1.5648E-01  1.2718E+01  1.0380E+00  6.3553E-01  1.5393E+00
             5.4308E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1646.16639950879        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      517
 NPARAMETR:  1.0092E+00  8.0507E-01  1.0318E+00  1.1596E+00  9.2344E-01  9.8455E-01  1.9130E+00  9.7710E-01  8.0693E-01  9.9577E-01
             1.0260E+00
 PARAMETER:  1.0918E-01 -1.1682E-01  1.3134E-01  2.4809E-01  2.0355E-02  8.4430E-02  7.4869E-01  7.6832E-02 -1.1452E-01  9.5759E-02
             1.2571E-01
 GRADIENT:   6.4683E+00  2.6852E+01  2.9767E+00  3.5494E+01 -6.5785E+00  1.9936E-01  8.0509E+00 -1.2051E-02 -2.2408E-01  8.5494E-02
            -3.1032E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1646.79134700218        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      672
 NPARAMETR:  9.9839E-01  7.5190E-01  1.0292E+00  1.1384E+00  9.2490E-01  9.7872E-01  1.8068E+00  9.6863E-01  8.0570E-01  9.9397E-01
             1.0259E+00
 PARAMETER:  9.8387E-02 -1.8515E-01  1.2879E-01  2.2959E-01  2.1929E-02  7.8492E-02  6.9158E-01  6.8125E-02 -1.1605E-01  9.3955E-02
             1.2561E-01
 GRADIENT:   2.4883E+01 -9.0201E+00  7.3303E-01 -6.8387E+00  5.4796E+00  3.0564E+00 -5.1463E-01  3.9569E-01  2.1801E-01  2.9825E-01
             3.9208E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1647.24340795096        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      828
 NPARAMETR:  1.0080E+00  7.5741E-01  1.0263E+00  1.1562E+00  9.1644E-01  9.8431E-01  1.8640E+00  9.5496E-01  8.0628E-01  9.9074E-01
             1.0260E+00
 PARAMETER:  1.0797E-01 -1.7785E-01  1.2596E-01  2.4516E-01  1.2744E-02  8.4185E-02  7.2270E-01  5.3919E-02 -1.1532E-01  9.0695E-02
             1.2566E-01
 GRADIENT:   4.4914E+00  3.2529E+00  1.8514E+00 -1.5619E+00 -2.7704E-01  3.0049E-01  1.6421E-01 -1.3087E-01  4.2325E-01  2.9928E-01
             2.5925E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1647.73288706169        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1010
 NPARAMETR:  1.0046E+00  5.8557E-01  9.9508E-01  1.2488E+00  8.3905E-01  9.8177E-01  2.2337E+00  8.8901E-01  7.6146E-01  9.3454E-01
             1.0263E+00
 PARAMETER:  1.0462E-01 -4.3517E-01  9.5064E-02  3.2219E-01 -7.5488E-02  8.1602E-02  9.0364E-01 -1.7645E-02 -1.7252E-01  3.2304E-02
             1.2598E-01
 GRADIENT:   7.5181E-01  3.9535E-01 -7.7325E-01  7.1938E-02 -1.4193E+00  1.1780E-01 -1.5339E-01 -1.4062E-01 -5.5330E-02 -8.5646E-02
             2.8343E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1647.78434372207        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1189
 NPARAMETR:  1.0041E+00  5.4110E-01  1.0308E+00  1.2740E+00  8.4338E-01  9.8119E-01  2.3538E+00  9.3197E-01  7.5163E-01  9.4571E-01
             1.0262E+00
 PARAMETER:  1.0414E-01 -5.1415E-01  1.3030E-01  3.4216E-01 -7.0342E-02  8.1010E-02  9.5604E-01  2.9542E-02 -1.8551E-01  4.4180E-02
             1.2583E-01
 GRADIENT:   1.7034E+00 -1.2480E+00 -1.2018E+00 -3.8942E+00 -6.3251E-01  2.4543E-01 -1.2622E-01  8.3203E-02  1.8994E-01  7.0981E-02
             5.5389E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1647.78611061987        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1380
 NPARAMETR:  1.0041E+00  5.4056E-01  1.0320E+00  1.2744E+00  8.4380E-01  9.8109E-01  2.3553E+00  9.3325E-01  7.5137E-01  9.4616E-01
             1.0261E+00
 PARAMETER:  1.0411E-01 -5.1512E-01  1.3146E-01  3.4248E-01 -6.9843E-02  8.0988E-02  9.5672E-01  3.1019E-02 -1.8572E-01  4.4711E-02
             1.2582E-01
 GRADIENT:   1.1860E+06  4.7941E+05  1.8786E+06  7.2106E+05 -1.2348E+06  2.2873E-01  2.5813E+05  7.8140E-02  1.7844E-01  6.5391E-02
             4.1509E-02
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         2.0         3.3         1.9         2.0         2.2
                    2.8

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1380
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2681E-04  2.8116E-02 -3.5313E-02 -2.7930E-02 -1.8129E-02
 SE:             2.9832E-02  2.0861E-02  1.5352E-02  2.2586E-02  2.0573E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9126E-01  1.7773E-01  2.1436E-02  2.1623E-01  3.7820E-01

 ETASHRINKSD(%)  5.8277E-02  3.0114E+01  4.8569E+01  2.4334E+01  3.1079E+01
 ETASHRINKVR(%)  1.1652E-01  5.1159E+01  7.3548E+01  4.2747E+01  5.2499E+01
 EBVSHRINKSD(%)  4.9725E-01  3.1737E+01  5.1860E+01  2.2779E+01  2.7556E+01
 EBVSHRINKVR(%)  9.9202E-01  5.3402E+01  7.6825E+01  4.0369E+01  4.7519E+01
 RELATIVEINF(%)  9.8297E+01  7.4708E+00  4.0779E+00  1.0052E+01  9.0637E+00
 EPSSHRINKSD(%)  4.5007E+01
 EPSSHRINKVR(%)  6.9758E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1647.7861106198748     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -912.63528405613658     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.60
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1647.786       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  5.41E-01  1.03E+00  1.27E+00  8.44E-01  9.81E-01  2.36E+00  9.33E-01  7.51E-01  9.46E-01  1.03E+00
 


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
+        5.65E+09
 
 TH 2
+        3.41E+03  7.96E+08
 
 TH 3
+        7.04E+03  6.65E+04  3.35E+09
 
 TH 4
+        2.18E+03 -5.32E+03  4.22E+04  3.24E+08
 
 TH 5
+       -1.13E+04  2.01E+04 -2.20E+05  1.31E+04  8.67E+09
 
 TH 6
+        1.11E+04  4.15E+03  8.52E+03  2.65E+03 -1.37E+04  2.07E+02
 
 TH 7
+        3.00E+03  1.16E+03  2.31E+03  7.07E+02 -3.71E+03  5.13E+02  1.22E+07
 
 TH 8
+       -5.97E+05 -2.24E+05 -4.60E+05 -1.43E+05  7.39E+05  1.10E+00 -2.77E+04  3.88E+01
 
 TH 9
+       -4.96E+04 -1.86E+04 -3.82E+04 -1.19E+04  6.14E+04 -3.50E-01 -2.29E+03 -4.74E+09  3.17E+09
 
 TH10
+       -2.75E+03 -1.03E+03 -2.13E+03 -6.76E+02  3.31E+03  2.10E+00 -1.27E+02 -6.99E+09  4.68E+09  6.90E+09
 
 TH11
+       -1.69E+06 -6.36E+05 -1.30E+06 -4.06E+05  2.10E+06  3.43E+00 -7.85E+04  2.44E+00  6.49E+00  9.64E+00  2.05E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.104
Stop Time:
Sat Sep 18 11:03:58 CDT 2021
