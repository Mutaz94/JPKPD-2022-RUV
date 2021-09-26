Sat Sep 25 13:48:43 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat81.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1688.79320229922        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.0877E+00 -6.1160E+01 -3.4926E+01 -2.9306E+01  6.0203E+01 -3.3169E+01 -5.6586E-01 -2.4834E+00  2.0743E+01  8.4599E+00
             1.7205E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1691.57244599665        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.9189E-01  1.0485E+00  1.0101E+00  9.6071E-01  1.0292E+00  1.0950E+00  1.0244E+00  1.0635E+00  8.8017E-01  9.2583E-01
             1.0045E+00
 PARAMETER:  9.1862E-02  1.4734E-01  1.1007E-01  5.9914E-02  1.2876E-01  1.9080E-01  1.2409E-01  1.6154E-01 -2.7644E-02  2.2937E-02
             1.0452E-01
 GRADIENT:  -1.2249E+01 -7.0264E+01 -2.5640E+01 -5.3580E+01  7.3660E+01  1.4476E+01 -7.2022E+00 -4.3418E+00 -3.4345E+00 -4.5828E+00
             1.9132E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1691.99397162290        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      269             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9785E-01  1.0533E+00  1.0037E+00  9.5750E-01  1.0284E+00  1.0927E+00  1.0334E+00  1.0769E+00  8.7588E-01  9.2150E-01
             1.0051E+00
 PARAMETER:  9.7849E-02  1.5190E-01  1.0374E-01  5.6566E-02  1.2804E-01  1.8861E-01  1.3289E-01  1.7407E-01 -3.2526E-02  1.8251E-02
             1.0504E-01
 GRADIENT:  -1.6127E+00 -6.9913E+01 -2.6153E+01 -5.2811E+01  7.3475E+01  1.3719E+01 -6.5562E+00 -3.6816E+00 -3.5105E+00 -4.4578E+00
             3.9840E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1692.03455173299        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:      464
 NPARAMETR:  9.9786E-01  1.0533E+00  1.0038E+00  9.5760E-01  1.0284E+00  1.0927E+00  1.0334E+00  1.0769E+00  8.7589E-01  9.2975E-01
             1.0051E+00
 PARAMETER:  9.7855E-02  1.5189E-01  1.0375E-01  5.6679E-02  1.2805E-01  1.8863E-01  1.3288E-01  1.7406E-01 -3.2520E-02  2.7164E-02
             1.0505E-01
 GRADIENT:  -1.5937E+00 -6.9604E+01 -2.6119E+01 -5.2977E+01  7.1493E+01  1.3730E+01 -6.4000E+00 -3.4620E+00 -3.4458E+00 -3.1553E+00
             6.4972E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1692.88711587166        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  1.0225E+00  1.0533E+00  1.0037E+00  9.5763E-01  1.0284E+00  1.1005E+00  1.1022E+00  1.0769E+00  8.9496E-01  9.3097E-01
             1.0036E+00
 PARAMETER:  1.2229E-01  1.5192E-01  1.0372E-01  5.6704E-02  1.2802E-01  1.9578E-01  1.9732E-01  1.7410E-01 -1.0971E-02  2.8467E-02
             1.0360E-01
 GRADIENT:   5.4315E+01 -6.2318E+01 -2.4927E+01 -4.7235E+01  6.8454E+01  1.7381E+01  1.8662E+00 -2.8638E+00  3.6385E+00 -1.8427E+00
             6.2900E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1692.92188315678        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      729
 NPARAMETR:  1.0239E+00  1.0533E+00  1.0037E+00  9.5763E-01  1.0284E+00  1.0995E+00  1.0975E+00  1.0769E+00  8.8002E-01  9.3097E-01
             1.0028E+00
 PARAMETER:  1.2360E-01  1.5192E-01  1.0372E-01  5.6704E-02  1.2802E-01  1.9490E-01  1.9302E-01  1.7410E-01 -2.7805E-02  2.8467E-02
             1.0279E-01
 GRADIENT:  -7.5378E-01 -6.9960E+01 -2.4553E+01 -5.7142E+01  6.7387E+01 -1.3039E-01 -4.4820E-01 -3.1013E+00  5.3589E-01 -2.1327E+00
            -1.3135E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1693.05831181439        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  1.0238E+00  1.0533E+00  1.0037E+00  9.5993E-01  1.0284E+00  1.0997E+00  1.1004E+00  1.0769E+00  8.7795E-01  9.3191E-01
             1.0030E+00
 PARAMETER:  1.2355E-01  1.5192E-01  1.0372E-01  5.9101E-02  1.2802E-01  1.9502E-01  1.9568E-01  1.7410E-01 -3.0163E-02  2.9477E-02
             1.0296E-01
 GRADIENT:   5.7188E+01 -6.0390E+01 -2.4904E+01 -4.5279E+01  6.8707E+01  1.6986E+01  7.4917E-01 -3.0269E+00  1.0869E+00 -1.9922E+00
             8.9716E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1693.46367165072        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  1.0240E+00  1.0533E+00  1.0037E+00  9.6804E-01  1.0284E+00  1.0997E+00  1.0933E+00  1.0769E+00  8.7795E-01  9.3854E-01
             1.0030E+00
 PARAMETER:  1.2373E-01  1.5192E-01  1.0372E-01  6.7517E-02  1.2802E-01  1.9502E-01  1.8916E-01  1.7410E-01 -3.0163E-02  3.6569E-02
             1.0296E-01
 GRADIENT:  -6.8436E-01 -5.9985E+01 -2.7551E+01 -4.0292E+01  6.9236E+01 -1.2466E-01 -3.3759E-01 -2.8499E+00  9.5971E-01 -1.4441E+00
             1.1010E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1693.48300212348        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:     1167
 NPARAMETR:  1.0240E+00  1.0533E+00  1.0037E+00  9.6850E-01  1.0284E+00  1.0997E+00  1.0930E+00  1.0769E+00  8.7795E-01  9.3886E-01
             1.0030E+00
 PARAMETER:  1.2374E-01  1.5192E-01  1.0372E-01  6.7991E-02  1.2802E-01  1.9502E-01  1.8894E-01  1.7410E-01 -3.0163E-02  3.6907E-02
             1.0296E-01
 GRADIENT:  -9.8517E+05  1.6048E+06 -1.1754E+06  2.4381E+06 -1.9045E+06 -6.2512E+05  1.2905E+06  1.4004E+06 -2.4382E+06 -1.2191E+06
            -2.3685E+06

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1167
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.5339E-04  2.0677E-02 -1.8238E-02  1.7313E-02 -6.2143E-02
 SE:             2.9885E-02  2.0604E-02  1.4171E-02  2.2847E-02  2.1228E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8523E-01  3.1558E-01  1.9810E-01  4.4857E-01  3.4183E-03

 ETASHRINKSD(%)  1.0000E-10  3.0975E+01  5.2526E+01  2.3461E+01  2.8883E+01
 ETASHRINKVR(%)  1.0000E-10  5.2355E+01  7.7462E+01  4.1417E+01  4.9424E+01
 EBVSHRINKSD(%)  3.6418E-01  3.0627E+01  5.6842E+01  2.3799E+01  2.6089E+01
 EBVSHRINKVR(%)  7.2703E-01  5.1873E+01  8.1374E+01  4.1935E+01  4.5372E+01
 RELATIVEINF(%)  9.8641E+01  1.4753E+00  1.6358E+00  1.8870E+00  9.4556E+00
 EPSSHRINKSD(%)  4.4950E+01
 EPSSHRINKVR(%)  6.9695E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1693.4830021234777     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -958.33217555973954     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.30
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1693.483       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.05E+00  1.00E+00  9.68E-01  1.03E+00  1.10E+00  1.09E+00  1.08E+00  8.78E-01  9.39E-01  1.00E+00
 


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
+        3.80E+09
 
 TH 2
+        3.54E+03  2.38E+09
 
 TH 3
+       -5.44E+03  9.33E+04  5.62E+09
 
 TH 4
+       -4.97E+09 -3.48E+04  5.40E+04  6.50E+09
 
 TH 5
+       -4.31E+03 -2.89E+09 -1.14E+05  4.30E+04  3.52E+09
 
 TH 6
+       -2.65E+03  5.92E+03 -9.10E+03  2.63E+04 -7.20E+03  1.33E+09
 
 TH 7
+       -2.45E+03  1.98E+03 -2.97E+03  3.05E+09 -2.37E+03 -1.45E+03  1.43E+09
 
 TH 8
+        3.03E+03 -1.07E+04  7.95E+04 -3.01E+04  2.53E+04  5.05E+03  1.65E+03  1.73E+09
 
 TH 9
+       -6.46E+03 -1.77E+04  2.71E+04  6.43E+04  2.15E+04 -1.08E+04 -3.50E+03 -1.51E+04  7.91E+09
 
 TH10
+        8.69E+03 -6.87E+03  1.06E+04 -1.14E+04  8.24E+03  5.13E+03 -5.32E+03 -5.86E+03  1.25E+04  6.92E+09
 
 TH11
+       -5.50E+03  1.94E+04 -1.44E+05  5.46E+04 -4.59E+04 -9.17E+03 -3.00E+03 -1.00E+06  2.74E+04  1.07E+04  5.72E+09
 
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
 
 Elapsed finaloutput time in seconds:     0.15
 #CPUT: Total CPU Time in Seconds,       22.742
Stop Time:
Sat Sep 25 13:49:09 CDT 2021
