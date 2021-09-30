Wed Sep 29 23:02:25 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat9.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1388.80949833333        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2250E+02  2.2729E+01  3.0882E+01  1.3574E+01  1.5186E+02  2.2961E+01 -3.8243E+01 -1.1611E+01 -2.7543E+01 -8.1365E+01
            -1.1380E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1714.73266198192        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8463E-01  1.0027E+00  9.6682E-01  1.0748E+00  9.1703E-01  1.0015E+00  9.7926E-01  9.4341E-01  9.1508E-01  9.3097E-01
             2.2242E+00
 PARAMETER:  8.4514E-02  1.0268E-01  6.6261E-02  1.7218E-01  1.3380E-02  1.0146E-01  7.9046E-02  4.1741E-02  1.1254E-02  2.8468E-02
             8.9939E-01
 GRADIENT:   1.6546E+02  2.1587E+01  1.0715E+00  3.7376E+01  2.2098E+00  1.6939E+01 -8.0811E-01  4.4291E+00 -1.1455E+00  7.0132E+00
             3.0094E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1721.29847321401        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.7244E-01  8.4794E-01  4.5148E-01  1.1256E+00  5.4693E-01  1.0203E+00  1.1599E+00  4.9693E-01  9.3009E-01  4.4011E-01
             2.1551E+00
 PARAMETER:  7.2057E-02 -6.4950E-02 -6.9521E-01  2.1831E-01 -5.0344E-01  1.2009E-01  2.4834E-01 -5.9931E-01  2.7525E-02 -7.2074E-01
             8.6783E-01
 GRADIENT:   1.2936E+02  3.8333E+01 -1.3970E+01  1.1439E+02  4.0179E+01  2.4474E+01 -1.2278E+00  2.8183E+00  1.1305E+01 -1.0115E+00
             1.7953E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1731.02156491382        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      339
 NPARAMETR:  9.4102E-01  7.6459E-01  4.5956E-01  1.1387E+00  5.1973E-01  9.5850E-01  1.3898E+00  2.1327E-01  8.3494E-01  5.7296E-01
             2.0859E+00
 PARAMETER:  3.9214E-02 -1.6842E-01 -6.7748E-01  2.2989E-01 -5.5445E-01  5.7616E-02  4.2913E-01 -1.4452E+00 -8.0394E-02 -4.5694E-01
             8.3522E-01
 GRADIENT:  -6.1233E+00  3.9402E+01 -8.9401E+00  4.8593E+01  2.4677E+01 -2.1874E+00  3.3837E-01  5.8500E-01 -2.6588E+01 -1.1025E+01
             4.5508E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1747.81155472379        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  9.4905E-01  4.3119E-01  3.5156E-01  1.1889E+00  3.5994E-01  9.8550E-01  1.5822E+00  1.4779E-01  9.1071E-01  5.6431E-01
             1.9270E+00
 PARAMETER:  4.7710E-02 -7.4122E-01 -9.4537E-01  2.7299E-01 -9.2182E-01  8.5392E-02  5.5880E-01 -1.8120E+00  6.4665E-03 -4.7215E-01
             7.5597E-01
 GRADIENT:   2.2980E+01 -8.0435E+00 -2.1955E+00 -2.9399E+01  1.7390E+01  8.3700E+00 -2.4486E+00  1.1058E-01 -8.8102E+00 -7.6506E+00
             1.0367E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1749.18828314073        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      690
 NPARAMETR:  9.3882E-01  4.0107E-01  3.5586E-01  1.2181E+00  3.5316E-01  9.6237E-01  1.5739E+00  1.0450E-01  9.4013E-01  6.5086E-01
             1.8910E+00
 PARAMETER:  3.6864E-02 -8.1361E-01 -9.3322E-01  2.9727E-01 -9.4084E-01  6.1644E-02  5.5359E-01 -2.1586E+00  3.8265E-02 -3.2946E-01
             7.3709E-01
 GRADIENT:   3.1505E-01 -4.5233E-01 -4.6462E-01  1.2339E-01  8.4325E-01 -1.4436E-01 -8.0832E-02  1.6104E-01 -2.4060E-01  4.2475E-01
             1.0965E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1749.25614415879        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      865
 NPARAMETR:  9.3907E-01  4.1000E-01  3.4874E-01  1.2111E+00  3.5014E-01  9.6367E-01  1.5308E+00  4.7540E-02  9.4461E-01  6.4602E-01
             1.8900E+00
 PARAMETER:  3.7135E-02 -7.9161E-01 -9.5344E-01  2.9149E-01 -9.4941E-01  6.2991E-02  5.2578E-01 -2.9462E+00  4.3016E-02 -3.3692E-01
             7.3660E-01
 GRADIENT:   1.0958E-01  7.6993E-01  7.3558E-01  1.1422E+00 -1.6566E+00  1.3536E-01 -4.0764E-01  3.2172E-02 -2.6020E-01  1.5023E-03
             2.1902E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1749.27313353567        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1041
 NPARAMETR:  9.3882E-01  4.0511E-01  3.5234E-01  1.2143E+00  3.5175E-01  9.6301E-01  1.5626E+00  1.0000E-02  9.4265E-01  6.4860E-01
             1.8902E+00
 PARAMETER:  3.6872E-02 -8.0360E-01 -9.4315E-01  2.9421E-01 -9.4483E-01  6.2306E-02  5.4633E-01 -4.8855E+00  4.0941E-02 -3.3295E-01
             7.3671E-01
 GRADIENT:   2.1486E-03 -2.3831E-02 -1.8972E-02 -1.6179E-02  4.7829E-02 -2.6290E-03  1.9154E-03  0.0000E+00  2.6874E-03  6.2752E-03
             3.4861E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1749.27313353567        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  9.3882E-01  4.0511E-01  3.5234E-01  1.2143E+00  3.5175E-01  9.6301E-01  1.5626E+00  1.0000E-02  9.4265E-01  6.4860E-01
             1.8902E+00
 PARAMETER:  3.6872E-02 -8.0360E-01 -9.4315E-01  2.9421E-01 -9.4483E-01  6.2306E-02  5.4633E-01 -4.8855E+00  4.0941E-02 -3.3295E-01
             7.3671E-01
 GRADIENT:   2.1486E-03 -2.3831E-02 -1.8972E-02 -1.6179E-02  4.7829E-02 -2.6290E-03  1.9154E-03  0.0000E+00  2.6874E-03  6.2752E-03
             3.4861E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1063
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2094E-03  2.3104E-02 -2.6334E-04 -1.1568E-02  7.2708E-03
 SE:             2.9519E-02  1.6765E-02  2.1871E-04  2.7566E-02  2.1404E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6732E-01  1.6815E-01  2.2858E-01  6.7474E-01  7.3408E-01

 ETASHRINKSD(%)  1.1083E+00  4.3836E+01  9.9267E+01  7.6492E+00  2.8295E+01
 ETASHRINKVR(%)  2.2044E+00  6.8456E+01  9.9995E+01  1.4713E+01  4.8584E+01
 EBVSHRINKSD(%)  1.2735E+00  4.7049E+01  9.9210E+01  7.0389E+00  2.6325E+01
 EBVSHRINKVR(%)  2.5307E+00  7.1962E+01  9.9994E+01  1.3582E+01  4.5720E+01
 RELATIVEINF(%)  9.7214E+01  6.3691E+00  4.1099E-04  4.4903E+01  2.6646E+00
 EPSSHRINKSD(%)  3.0594E+01
 EPSSHRINKVR(%)  5.1828E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1749.2731335356668     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -830.33460033099414     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.40
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1749.273       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.39E-01  4.05E-01  3.52E-01  1.21E+00  3.52E-01  9.63E-01  1.56E+00  1.00E-02  9.43E-01  6.49E-01  1.89E+00
 


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
+        1.33E+03
 
 TH 2
+       -3.27E+01  7.91E+02
 
 TH 3
+        1.91E+01  1.33E+03  5.48E+03
 
 TH 4
+       -1.20E+01  2.28E+02 -5.35E+02  7.36E+02
 
 TH 5
+        3.38E+01 -2.42E+03 -7.55E+03  1.19E+02  1.17E+04
 
 TH 6
+        3.34E+00 -7.26E+00  1.34E+01 -6.18E+00 -6.55E+00  2.02E+02
 
 TH 7
+        1.53E+00  3.48E+01 -1.32E+00 -1.24E+00 -4.57E+01  1.75E-01  1.13E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.14E+00 -1.81E+01  4.43E+01  3.52E-02  1.24E+01 -8.60E-01  5.88E+00  0.00E+00  1.62E+02
 
 TH10
+       -2.28E+00  4.44E+01 -2.08E+02 -4.02E+00  1.18E+02  7.40E-01  1.54E+01  0.00E+00 -4.49E+00  1.51E+02
 
 TH11
+       -1.49E+01  1.71E+00 -5.97E+01 -1.34E+01  6.31E+01  3.32E+00  3.62E+00  0.00E+00  1.06E+01  2.04E+01  1.17E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       24.463
Stop Time:
Wed Sep 29 23:02:51 CDT 2021
