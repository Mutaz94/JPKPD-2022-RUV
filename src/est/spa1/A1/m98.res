Wed Sep 29 22:57:38 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat98.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m98.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1763.31024637469        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8692E+02  2.8867E+01  5.7333E+00  6.7818E+01  9.6636E+01  5.6839E+01 -3.0633E+01 -6.9090E+00 -1.2844E+00 -3.5543E+01
            -5.5603E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1859.23397327709        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.9189E-01  9.4963E-01  1.0202E+00  1.0727E+00  8.9379E-01  9.2086E-01  1.2363E+00  7.8000E-01  1.0274E+00  8.6843E-01
             2.0324E+00
 PARAMETER:  9.1853E-02  4.8315E-02  1.2003E-01  1.7017E-01 -1.2285E-02  1.7550E-02  3.1214E-01 -1.4846E-01  1.2707E-01 -4.1069E-02
             8.0921E-01
 GRADIENT:   1.4599E+02  3.7838E+01  1.6571E+01  4.6114E+01 -4.2539E+01  2.1743E+00  1.1746E+01  5.2092E+00  1.7952E+01  4.0251E+00
             1.7260E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1870.70476829440        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      178
 NPARAMETR:  9.5747E-01  7.5101E-01  5.5068E-01  1.1598E+00  5.9680E-01  9.6354E-01  1.5193E+00  1.9394E-01  8.7008E-01  6.8078E-01
             1.8751E+00
 PARAMETER:  5.6543E-02 -1.8633E-01 -4.9661E-01  2.4822E-01 -4.1617E-01  6.2858E-02  5.1828E-01 -1.5402E+00 -3.9173E-02 -2.8451E-01
             7.2867E-01
 GRADIENT:  -4.1624E+01  2.6146E+01 -4.8275E+01  9.6224E+01  3.6728E+01  7.3533E+00  1.0353E+01  6.6290E-01 -5.8729E-01  4.4470E+00
             1.3895E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1891.66582139963        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      355
 NPARAMETR:  9.7550E-01  5.7682E-01  8.8069E-01  1.2489E+00  7.1754E-01  9.4442E-01  1.6981E+00  2.8946E-01  8.4553E-01  8.5657E-01
             1.5235E+00
 PARAMETER:  7.5197E-02 -4.5022E-01 -2.7048E-02  3.2224E-01 -2.3193E-01  4.2814E-02  6.2950E-01 -1.1398E+00 -6.7796E-02 -5.4813E-02
             5.2101E-01
 GRADIENT:   2.6600E+01  1.1974E+01  2.3224E+01 -1.1069E+01 -2.1369E+01  2.0173E+00 -6.8826E-01  3.9098E-01 -5.5130E+00 -1.0322E+01
            -1.7045E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1897.07603369645        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  9.5728E-01  2.6813E-01  8.6315E-01  1.4086E+00  6.3640E-01  9.2884E-01  2.6579E+00  1.2444E-02  8.0638E-01  9.1601E-01
             1.5656E+00
 PARAMETER:  5.6340E-02 -1.2163E+00 -4.7165E-02  4.4260E-01 -3.5192E-01  2.6178E-02  1.0775E+00 -4.2865E+00 -1.1520E-01  1.2273E-02
             5.4829E-01
 GRADIENT:  -5.1242E+00  4.4615E+00 -2.6914E+00 -2.2721E+01  2.4454E+00 -1.6393E+00  4.9859E+00  1.5240E-03 -4.7541E+00  3.1701E-01
             1.4224E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1900.12740908623        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      710
 NPARAMETR:  9.5108E-01  6.7785E-02  9.5420E-01  1.5511E+00  6.4013E-01  9.3145E-01  4.9748E+00  1.0000E-02  7.7526E-01  9.8049E-01
             1.5385E+00
 PARAMETER:  4.9848E-02 -2.5914E+00  5.3116E-02  5.3894E-01 -3.4608E-01  2.8988E-02  1.7044E+00 -1.0544E+01 -1.5455E-01  8.0301E-02
             5.3079E-01
 GRADIENT:  -7.6476E+00  2.0694E+00 -2.9707E+00  2.9446E+01  1.4539E+00  8.2682E-01 -2.1506E-02  0.0000E+00 -3.8167E+00  3.3789E+00
            -2.1093E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1901.67335110615        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      893             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5308E-01  2.4844E-02  9.1117E-01  1.5361E+00  6.1074E-01  9.2758E-01  7.6161E+00  1.0000E-02  7.7628E-01  9.3215E-01
             1.5406E+00
 PARAMETER:  5.1946E-02 -3.5951E+00  6.9688E-03  5.2928E-01 -3.9308E-01  2.4829E-02  2.1303E+00 -1.5685E+01 -1.5324E-01  2.9735E-02
             5.3216E-01
 GRADIENT:   1.6454E+02  4.0848E+00  6.1328E+00  3.8773E+02  2.1233E+01  1.5506E+01  3.7432E+01  0.0000E+00  1.0033E+01  1.4574E-03
             5.8468E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1901.75261020731        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  9.5254E-01  2.4927E-02  9.0728E-01  1.5489E+00  6.1012E-01  9.2820E-01  7.6249E+00  1.0000E-02  7.7585E-01  9.3489E-01
             1.5402E+00
 PARAMETER:  5.1374E-02 -3.5918E+00  2.6959E-03  5.3755E-01 -3.9409E-01  2.5493E-02  2.1314E+00 -1.5685E+01 -1.5379E-01  3.2669E-02
             5.3194E-01
 GRADIENT:  -3.4708E-01 -9.1704E-01 -1.2424E+00  1.2278E+00  1.8204E+00 -5.2893E-02 -8.9634E-01  0.0000E+00  2.8367E-01  2.4638E-01
             2.8597E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1901.79479228646        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1199
 NPARAMETR:  9.5248E-01  2.5307E-02  9.0425E-01  1.5478E+00  6.0800E-01  9.2817E-01  7.6646E+00  1.0000E-02  7.7644E-01  9.3372E-01
             1.5402E+00
 PARAMETER:  5.1309E-02 -3.5767E+00 -6.4426E-04  5.3681E-01 -3.9759E-01  2.5458E-02  2.1366E+00 -1.5685E+01 -1.5304E-01  3.1426E-02
             5.3192E-01
 GRADIENT:  -4.7593E-01 -3.8351E-01 -9.7764E-02 -4.0079E-01  3.5059E-03 -7.2778E-02  2.2759E-01  0.0000E+00  1.9593E-01  1.5228E-02
             3.2850E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1199
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1146E-03  7.8946E-03 -1.1087E-04 -8.4012E-03 -1.5386E-02
 SE:             2.9706E-02  6.5774E-03  1.9173E-04  2.8757E-02  2.4585E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7007E-01  2.3004E-01  5.6310E-01  7.7018E-01  5.3143E-01

 ETASHRINKSD(%)  4.8094E-01  7.7965E+01  9.9358E+01  3.6588E+00  1.7637E+01
 ETASHRINKVR(%)  9.5957E-01  9.5145E+01  9.9996E+01  7.1837E+00  3.2164E+01
 EBVSHRINKSD(%)  8.4260E-01  8.5318E+01  9.9341E+01  3.5856E+00  1.6287E+01
 EBVSHRINKVR(%)  1.6781E+00  9.7844E+01  9.9996E+01  7.0426E+00  2.9922E+01
 RELATIVEINF(%)  9.7907E+01  1.0366E+00  4.4872E-04  4.1822E+01  7.4370E+00
 EPSSHRINKSD(%)  3.0338E+01
 EPSSHRINKVR(%)  5.1472E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1901.7947922864585     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -982.85625908178577     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1901.795       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.52E-01  2.53E-02  9.04E-01  1.55E+00  6.08E-01  9.28E-01  7.66E+00  1.00E-02  7.76E-01  9.34E-01  1.54E+00
 


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
+        1.40E+03
 
 TH 2
+       -1.24E+02  2.24E+05
 
 TH 3
+       -1.45E+01  5.14E+02  6.02E+02
 
 TH 4
+       -3.61E+00  1.20E+03 -1.62E+02  3.38E+03
 
 TH 5
+        2.82E+01 -1.61E+03 -1.16E+03 -9.28E+03  3.47E+04
 
 TH 6
+        1.74E+00  1.22E+00  3.17E-01 -4.10E+00 -5.26E-01  2.24E+02
 
 TH 7
+       -4.09E-01 -4.48E+01  2.17E+00  4.82E+00 -6.60E+00  4.13E-02  6.97E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.02E+00  8.61E+01  3.30E+01 -2.80E+01 -5.94E+01  2.01E-01  6.98E-01  0.00E+00  3.08E+02
 
 TH10
+       -1.65E+00  3.48E+02  1.66E+00 -3.90E+01 -1.36E+02  9.89E-01  1.87E+00  0.00E+00  1.86E+01  1.36E+02
 
 TH11
+       -1.46E+01  2.39E+02 -4.72E+00 -4.35E+01 -9.52E+03  1.89E+00  1.42E+00  0.00E+00  2.27E+01  4.06E+01  3.00E+03
 
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
 #CPUT: Total CPU Time in Seconds,       26.769
Stop Time:
Wed Sep 29 22:58:06 CDT 2021
