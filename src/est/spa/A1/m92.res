Sat Sep 25 08:20:07 CDT 2021
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
$DATA ../../../../data/spa/A1/dat92.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1456.23369270618        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.1563E+01 -5.8298E+01 -2.9154E+01 -4.2357E+01  9.9309E+01  1.9790E+01 -2.6822E+01 -3.7338E+00 -2.6249E+01 -1.7904E+01
            -3.2083E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1525.99105377117        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.7245E-01  1.0473E+00  1.0625E+00  1.0252E+00  9.6877E-01  9.2471E-01  1.0903E+00  9.7644E-01  1.0643E+00  9.6042E-01
             1.4650E+00
 PARAMETER:  7.2065E-02  1.4624E-01  1.6065E-01  1.2488E-01  6.8268E-02  2.1725E-02  1.8646E-01  7.6162E-02  1.6231E-01  5.9616E-02
             4.8187E-01
 GRADIENT:  -1.2834E+01  1.0826E+01  2.1065E+00  1.3158E+01 -2.6693E+00 -8.7062E+00 -7.5736E+00  2.2552E-01 -1.4676E+00 -2.3178E+00
            -4.1652E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1528.93285926534        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.7656E-01  8.9993E-01  8.2062E-01  1.1095E+00  8.1004E-01  9.5381E-01  1.5394E+00  4.9090E-01  8.7852E-01  8.3537E-01
             1.5035E+00
 PARAMETER:  7.6286E-02 -5.4416E-03 -9.7690E-02  2.0387E-01 -1.1067E-01  5.2714E-02  5.3136E-01 -6.1151E-01 -2.9516E-02 -7.9875E-02
             5.0782E-01
 GRADIENT:  -7.0012E+00  1.7305E+01 -1.3657E+01  3.8233E+01  2.3489E+01  2.7578E+00  5.7571E+00  1.0775E+00 -4.5900E+00 -2.9323E+00
            -2.1995E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1532.09571599018        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.7866E-01  6.1569E-01  6.4385E-01  1.2336E+00  5.9657E-01  9.4617E-01  1.9329E+00  1.2793E-01  7.9556E-01  7.2392E-01
             1.5520E+00
 PARAMETER:  7.8431E-02 -3.8501E-01 -3.4028E-01  3.0991E-01 -4.1657E-01  4.4666E-02  7.5904E-01 -1.9563E+00 -1.2870E-01 -2.2307E-01
             5.3952E-01
 GRADIENT:  -6.6627E+00  1.5288E+01  2.1099E+00  3.0232E+01 -8.9508E+00 -7.9597E-01  2.2089E+00  2.3191E-01 -4.1825E+00  1.2132E-01
             5.0940E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1532.37952941342        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.7934E-01  5.3591E-01  5.8457E-01  1.2466E+00  5.4115E-01  9.4902E-01  2.0510E+00  6.9928E-02  7.9343E-01  6.9589E-01
             1.5059E+00
 PARAMETER:  7.9119E-02 -5.2379E-01 -4.3688E-01  3.2044E-01 -5.1406E-01  4.7674E-02  8.1832E-01 -2.5603E+00 -1.3139E-01 -2.6256E-01
             5.0942E-01
 GRADIENT:  -3.5669E+00  5.7387E+00  1.8573E+00  9.5603E+00 -4.2544E+00 -3.6133E-02  5.0605E-01  8.7258E-02 -1.3151E+00 -5.9695E-01
            -3.1387E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1532.39041036013        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  9.7996E-01  5.0056E-01  5.6962E-01  1.2576E+00  5.2290E-01  9.4903E-01  2.1230E+00  5.2471E-02  7.9061E-01  6.9371E-01
             1.5037E+00
 PARAMETER:  7.9754E-02 -5.9203E-01 -4.6279E-01  3.2917E-01 -5.4836E-01  4.7681E-02  8.5281E-01 -2.8475E+00 -1.3496E-01 -2.6570E-01
             5.0795E-01
 GRADIENT:  -1.7956E+00  3.5712E+00  1.8507E+00  5.3642E+00 -3.7758E+00 -3.2131E-02  1.8073E-01  5.3045E-02 -8.3475E-01 -3.7053E-01
            -1.8552E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1532.95637608309        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      487
 NPARAMETR:  9.8593E-01  5.0335E-01  6.5315E-01  1.2728E+00  5.7086E-01  9.4870E-01  2.0985E+00  6.1902E-02  8.1218E-01  7.6314E-01
             1.5197E+00
 PARAMETER:  8.5826E-02 -5.8647E-01 -3.2595E-01  3.4124E-01 -4.6061E-01  4.7338E-02  8.4121E-01 -2.6822E+00 -1.0803E-01 -1.7031E-01
             5.1851E-01
 GRADIENT:   1.8711E-01  2.0089E-01  1.0392E+01 -2.2453E+01 -1.7734E+01 -7.6218E-01 -2.5411E+00  5.4610E-02  1.7474E+00  1.6082E+00
            -6.1430E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1533.76454813605        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      662
 NPARAMETR:  9.8239E-01  3.6173E-01  7.3718E-01  1.3747E+00  5.8746E-01  9.4592E-01  2.6355E+00  2.2288E-02  7.9088E-01  8.3433E-01
             1.5257E+00
 PARAMETER:  8.2229E-02 -9.1686E-01 -2.0492E-01  4.1822E-01 -4.3196E-01  4.4399E-02  1.0691E+00 -3.7037E+00 -1.3461E-01 -8.1132E-02
             5.2245E-01
 GRADIENT:   2.0124E-01  9.2514E-01  5.9171E-01 -7.5439E-02 -1.1891E+00 -3.6572E-01  5.0751E-01  4.5923E-03  1.2447E-01  3.6985E-01
             4.5586E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1533.79374255395        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      824
 NPARAMETR:  9.8151E-01  3.2269E-01  7.3800E-01  1.3948E+00  5.7937E-01  9.4612E-01  2.7934E+00  1.4317E-02  7.8707E-01  8.4144E-01
             1.5230E+00
 PARAMETER:  8.1332E-02 -1.0311E+00 -2.0381E-01  4.3277E-01 -4.4581E-01  4.4611E-02  1.1273E+00 -4.1463E+00 -1.3944E-01 -7.2638E-02
             5.2067E-01
 GRADIENT:   7.4465E-03  6.4488E-04  4.2506E-03  5.4621E-02 -2.2548E-02  7.1784E-04 -1.0470E-03  1.8939E-03  5.4931E-03 -3.6415E-04
            -4.4548E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      824
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.0535E-04  3.7174E-02 -5.0420E-04 -2.9906E-02  5.8347E-03
 SE:             2.9644E-02  1.8763E-02  2.7266E-04  2.4848E-02  2.1902E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7564E-01  4.7564E-02  6.4435E-02  2.2875E-01  7.8993E-01

 ETASHRINKSD(%)  6.8724E-01  3.7142E+01  9.9087E+01  1.6757E+01  2.6627E+01
 ETASHRINKVR(%)  1.3698E+00  6.0489E+01  9.9992E+01  3.0707E+01  4.6164E+01
 EBVSHRINKSD(%)  1.0539E+00  4.4740E+01  9.9022E+01  1.3762E+01  2.1845E+01
 EBVSHRINKVR(%)  2.0967E+00  6.9463E+01  9.9990E+01  2.5630E+01  3.8918E+01
 RELATIVEINF(%)  9.6988E+01  6.9744E+00  6.2444E-04  2.1391E+01  3.7577E+00
 EPSSHRINKSD(%)  4.0315E+01
 EPSSHRINKVR(%)  6.4377E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1533.7937425539531     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -798.64291599021487     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1533.794       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  3.23E-01  7.38E-01  1.39E+00  5.79E-01  9.46E-01  2.79E+00  1.43E-02  7.87E-01  8.41E-01  1.52E+00
 


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
+        1.26E+03
 
 TH 2
+       -2.99E+01  4.68E+02
 
 TH 3
+        2.46E+01  2.32E+02  1.01E+03
 
 TH 4
+       -1.65E+01  2.98E+02 -2.20E+02  7.04E+02
 
 TH 5
+        9.95E-01 -5.69E+02 -1.63E+03  1.63E+02  2.87E+03
 
 TH 6
+       -1.85E+00 -5.38E+00  6.20E+00 -4.49E+00 -9.24E-01  2.11E+02
 
 TH 7
+        1.28E+00  3.42E+01 -1.65E+00 -4.81E+00 -2.30E+00  4.41E-02  7.71E+00
 
 TH 8
+       -9.19E+00 -1.58E+00  1.85E+00 -7.78E-01  3.43E+00 -1.23E+00 -1.95E-01 -1.40E+00
 
 TH 9
+        5.51E+00 -1.68E+01 -1.38E+01 -1.59E+01  3.77E+01 -5.03E+00  4.28E+00  2.06E+00  1.97E+02
 
 TH10
+       -2.39E+00  1.96E+01 -5.14E+01 -2.94E+01 -6.17E+00  1.69E+00  2.73E+00 -3.60E+00 -1.91E+00  1.02E+02
 
 TH11
+       -1.05E+01 -7.40E+00 -3.80E+01 -1.03E+01  1.52E+01  3.23E+00  1.14E+00 -1.84E-01  8.22E+00  2.31E+01  1.02E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.972
Stop Time:
Sat Sep 25 08:20:23 CDT 2021
