Sat Sep 18 15:31:53 CDT 2021
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
$DATA ../../../../data/spa/D/dat65.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m65.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13639.6606314099        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.4285E+02  3.2157E+02 -2.0148E+01  1.6021E+01  3.8330E+02 -2.6560E+03 -8.8848E+02 -1.1354E+02 -1.7521E+03 -1.0503E+03
            -2.4281E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -609.723778717470        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3432E+00  9.6588E-01  9.2592E-01  1.5872E+00  1.3532E+00  1.8364E+00  1.2097E+00  9.6843E-01  1.3348E+00  1.1739E+00
             1.4803E+01
 PARAMETER:  3.9502E-01  6.5286E-02  2.3031E-02  5.6198E-01  4.0247E-01  7.0782E-01  2.9037E-01  6.7925E-02  3.8875E-01  2.6031E-01
             2.7948E+00
 GRADIENT:  -3.2191E+01  1.3855E+01 -8.8748E+00  1.8510E+01 -6.7638E+00  1.4207E+01  9.0355E-02  4.9682E+00  8.2814E+00  1.8884E+00
             9.3077E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -621.587723365675        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3188E+00  9.2383E-01  1.2839E+00  1.8865E+00  3.2646E+00  2.3072E+00  3.7512E+00  4.5029E-01  3.5227E+00  4.9791E+00
             1.0909E+01
 PARAMETER:  3.7675E-01  2.0768E-02  3.4987E-01  7.3470E-01  1.2831E+00  9.3604E-01  1.4221E+00 -6.9786E-01  1.3592E+00  1.7053E+00
             2.4896E+00
 GRADIENT:  -3.3032E+01  9.6659E+00  1.1556E+01  1.3108E+01 -1.6501E+01 -9.6112E+01  1.5011E+01 -6.9476E-01  5.8388E+01  4.2238E+00
             6.1277E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -666.403156102744        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.3974E+00  3.9328E-01  7.9029E-01  1.6884E+00  3.9252E+00  2.7808E+00  2.9745E+00  1.0000E-02  2.0766E+00  4.2626E+00
             1.0402E+01
 PARAMETER:  4.3464E-01 -8.3323E-01 -1.3535E-01  6.2377E-01  1.4674E+00  1.1227E+00  1.1901E+00 -6.9280E+00  8.3073E-01  1.5499E+00
             2.4420E+00
 GRADIENT:   1.9211E+01  1.0112E+01  2.7252E+01 -3.0132E+01 -2.5438E+01  8.5012E+00  4.3990E+00  0.0000E+00 -1.3772E+01  1.3478E+00
             3.6666E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -710.860468506364        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.3702E+00  1.8798E-01  3.0009E-01  1.3301E+00  1.1355E+01  2.7699E+00  3.8540E-01  1.0000E-02  1.5053E+00  6.4604E+00
             9.3101E+00
 PARAMETER:  4.1499E-01 -1.5714E+00 -1.1037E+00  3.8525E-01  2.5297E+00  1.1188E+00 -8.5347E-01 -1.5556E+01  5.0902E-01  1.9657E+00
             2.3311E+00
 GRADIENT:   6.8616E+01  2.7980E+00  3.5328E+01 -1.4452E+01 -6.1339E-01  2.0084E+01  2.7684E-01  0.0000E+00 -2.4458E+01  4.2840E+00
            -3.8732E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -750.317631017631        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  8.0897E-01  8.3741E-02  9.9937E-02  8.3025E-01  2.0218E+01  2.2325E+00  1.0000E-02  1.0000E-02  1.2540E+00  9.6601E+00
             8.9561E+00
 PARAMETER: -1.1200E-01 -2.3800E+00 -2.2032E+00 -8.6031E-02  3.1066E+00  9.0312E-01 -5.7588E+00 -2.2262E+01  3.2631E-01  2.3680E+00
             2.2923E+00
 GRADIENT:  -4.4865E+00  3.1743E+01 -3.4638E+01  5.4799E+01 -9.9955E+00 -8.1732E+00  0.0000E+00  0.0000E+00  3.2689E+00  1.3756E+01
            -4.6603E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -793.900130014956        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  4.1698E-01  1.0320E-02  1.6242E-02  2.0909E-01  4.5597E+01  1.9509E+00  1.0000E-02  1.0000E-02  8.1673E-01  9.1228E+00
             8.9113E+00
 PARAMETER: -7.7472E-01 -4.4737E+00 -4.0201E+00 -1.4650E+00  3.9198E+00  7.6830E-01 -1.7432E+01 -3.9418E+01 -1.0245E-01  2.3108E+00
             2.2873E+00
 GRADIENT:   1.3703E+01  1.3827E+00 -4.7865E+01  7.1252E+01  2.8733E-02 -1.5127E+01  0.0000E+00  0.0000E+00 -2.3709E+01 -9.4185E-04
            -3.5311E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -796.782180669150        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      582
 NPARAMETR:  4.0959E-01  1.0000E-02  1.6731E-02  2.0327E-01  3.8721E+01  2.0235E+00  1.0000E-02  1.0000E-02  8.7921E-01  8.7706E+00
             9.1704E+00
 PARAMETER: -7.9261E-01 -4.6078E+00 -3.9905E+00 -1.4932E+00  3.7564E+00  8.0483E-01 -1.8252E+01 -4.0134E+01 -2.8731E-02  2.2714E+00
             2.3160E+00
 GRADIENT:  -1.0640E+01  0.0000E+00  2.9693E+01 -2.9776E+01 -1.7744E-02 -2.4926E+00  0.0000E+00  0.0000E+00 -7.5666E+00 -1.1351E-03
            -3.1035E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -797.447324444543        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      757
 NPARAMETR:  4.1164E-01  1.0000E-02  1.6019E-02  2.0031E-01  3.8828E+01  2.0358E+00  1.0000E-02  1.0000E-02  9.2951E-01  9.4871E+00
             9.1729E+00
 PARAMETER: -7.8761E-01 -4.6632E+00 -4.0340E+00 -1.5079E+00  3.7592E+00  8.1091E-01 -1.8720E+01 -4.0767E+01  2.6898E-02  2.3499E+00
             2.3163E+00
 GRADIENT:   1.9839E-01  0.0000E+00 -3.1656E-01  2.1472E-01 -4.8177E-03  2.3470E-01  0.0000E+00  0.0000E+00 -1.4545E-01  1.9974E-04
             2.2123E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -797.447926580549        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      942
 NPARAMETR:  4.1197E-01  1.0000E-02  1.6070E-02  2.0080E-01  4.7925E+01  2.0350E+00  1.0000E-02  1.0000E-02  9.3061E-01  9.4996E+00
             9.1689E+00
 PARAMETER: -7.8680E-01 -4.6596E+00 -4.0308E+00 -1.5055E+00  3.9696E+00  8.1049E-01 -1.8699E+01 -4.0728E+01  2.8084E-02  2.3513E+00
             2.3158E+00
 GRADIENT:   1.2425E-01  0.0000E+00 -2.6147E-02  3.3066E-02 -1.1183E-04  9.4080E-03  0.0000E+00  0.0000E+00  5.0593E-03 -3.3338E-04
            -1.9967E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -797.447938098430        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1106
 NPARAMETR:  4.1190E-01  1.0000E-02  1.6071E-02  2.0080E-01  4.8746E+01  2.0350E+00  1.0000E-02  1.0000E-02  9.3054E-01  9.5060E+00
             9.1691E+00
 PARAMETER: -7.8698E-01 -4.6596E+00 -4.0308E+00 -1.5055E+00  3.9866E+00  8.1050E-01 -1.8699E+01 -4.0728E+01  2.8009E-02  2.3519E+00
             2.3158E+00
 GRADIENT:   3.4865E-02  0.0000E+00  1.4677E-02  3.2266E-02  6.6552E-06  1.1578E-02  0.0000E+00  0.0000E+00 -5.3140E-03 -3.4782E-04
             8.8763E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1106
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.0895E-04 -2.1591E-07  1.2962E-04 -2.1734E-02 -1.5439E-04
 SE:             2.9311E-02  1.7197E-06  2.2743E-04  2.4067E-02  2.0023E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8887E-01  9.0008E-01  5.6873E-01  3.6649E-01  4.4067E-01

 ETASHRINKSD(%)  1.8054E+00  9.9994E+01  9.9238E+01  1.9373E+01  9.9329E+01
 ETASHRINKVR(%)  3.5782E+00  1.0000E+02  9.9994E+01  3.4993E+01  9.9996E+01
 EBVSHRINKSD(%)  1.7329E+00  9.9987E+01  9.9269E+01  2.0048E+01  9.9381E+01
 EBVSHRINKVR(%)  3.4357E+00  1.0000E+02  9.9995E+01  3.6077E+01  9.9996E+01
 RELATIVEINF(%)  2.7587E+00  4.2125E-07  2.6094E-05  2.8746E-01  1.9029E-04
 EPSSHRINKSD(%)  1.4624E+01
 EPSSHRINKVR(%)  2.7110E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -797.44793809843009     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -62.297111534691908     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.77
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -797.448       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.12E-01  1.00E-02  1.61E-02  2.01E-01  4.87E+01  2.04E+00  1.00E-02  1.00E-02  9.31E-01  9.51E+00  9.17E+00
 


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
+        1.52E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.39E+04  0.00E+00  2.17E+06
 
 TH 4
+       -1.79E+02  0.00E+00 -1.93E+05  1.89E+04
 
 TH 5
+        1.15E-02  0.00E+00 -2.53E-01  2.51E-02  1.28E-06
 
 TH 6
+        1.79E+00  0.00E+00 -9.01E+01 -2.03E+01  6.62E-05  4.40E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.66E+00  0.00E+00  1.68E+03 -1.72E+02 -3.86E-03 -9.08E-02  0.00E+00  0.00E+00  8.63E+01
 
 TH10
+       -1.50E-03  0.00E+00 -7.34E-02 -2.82E-03 -2.76E-06  7.85E-04  0.00E+00  0.00E+00  2.84E-02 -9.78E-05
 
 TH11
+       -1.39E+01  0.00E+00  3.39E+02 -2.17E+01 -1.57E-04  9.80E-01  0.00E+00  0.00E+00  3.59E+00 -9.35E-06  4.13E+00
 
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
 #CPUT: Total CPU Time in Seconds,       24.276
Stop Time:
Sat Sep 18 15:32:19 CDT 2021
