Sat Sep 18 14:17:16 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat80.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1717.43041807304        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.2906E+01 -2.7151E+01  9.1046E+00 -4.8634E+01 -3.7586E+01  3.6930E+01  1.2959E+01  7.0380E+00  3.4135E+01  1.9735E+01
            -3.2579E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1725.89217370901        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0374E+00  1.0623E+00  1.0400E+00  9.9922E-01  1.0559E+00  8.5098E-01  9.1066E-01  9.4921E-01  8.0005E-01  8.8436E-01
             1.1404E+00
 PARAMETER:  1.3672E-01  1.6045E-01  1.3925E-01  9.9223E-02  1.5436E-01 -6.1371E-02  6.4121E-03  4.7871E-02 -1.2308E-01 -2.2895E-02
             2.3136E-01
 GRADIENT:   7.5409E+01  9.9314E+00  1.4932E+01 -8.1694E+00 -8.5075E+00 -2.4171E+01 -1.5595E+00 -1.6788E+00 -9.8881E+00 -9.1780E-01
             1.2045E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1727.01056397776        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0309E+00  1.0628E+00  7.6914E-01  9.7983E-01  9.0201E-01  8.8273E-01  1.0056E+00  7.4715E-01  7.5975E-01  6.6988E-01
             1.1357E+00
 PARAMETER:  1.3040E-01  1.6089E-01 -1.6248E-01  7.9626E-02 -3.1315E-03 -2.4742E-02  1.0562E-01 -1.9149E-01 -1.7476E-01 -3.0066E-01
             2.2723E-01
 GRADIENT:   4.6377E+01  7.4284E+00  7.3601E+00 -5.8784E+00 -6.8843E+00 -9.3858E+00  6.5649E-01  1.1930E+00 -1.1304E+01 -4.4351E+00
             1.3031E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1728.54324965940        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0154E+00  1.1072E+00  6.6287E-01  9.4253E-01  8.6644E-01  9.0380E-01  9.4511E-01  4.4927E-01  8.3432E-01  7.1971E-01
             1.0884E+00
 PARAMETER:  1.1531E-01  2.0179E-01 -3.1118E-01  4.0808E-02 -4.3363E-02 -1.1513E-03  4.3548E-02 -7.0013E-01 -8.1143E-02 -2.2890E-01
             1.8469E-01
 GRADIENT:   1.7054E+00  8.5608E-01 -2.7766E+00  2.1417E+00  2.1794E-01 -8.4577E-01  1.3873E-01  1.1670E+00  1.6612E-01  1.2181E+00
             1.7103E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1728.58097270851        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0136E+00  1.0768E+00  6.5762E-01  9.5792E-01  8.5005E-01  9.0713E-01  9.7470E-01  3.2621E-01  8.2505E-01  7.2147E-01
             1.0842E+00
 PARAMETER:  1.1352E-01  1.7396E-01 -3.1913E-01  5.7013E-02 -6.2465E-02  2.5352E-03  7.4378E-02 -1.0202E+00 -9.2315E-02 -2.2646E-01
             1.8081E-01
 GRADIENT:  -3.0999E+00 -2.8579E+00 -2.9752E+00  5.3758E-01  3.1992E+00  5.4843E-01  7.6606E-01  5.2781E-01  8.6089E-01  9.8341E-01
             4.0114E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1728.78608423207        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0159E+00  1.1087E+00  5.9690E-01  9.3133E-01  8.2204E-01  9.0522E-01  9.5321E-01  7.1595E-02  8.2399E-01  6.7246E-01
             1.0845E+00
 PARAMETER:  1.1581E-01  2.0322E-01 -4.1601E-01  2.8861E-02 -9.5962E-02  4.2704E-04  5.2077E-02 -2.5367E+00 -9.3598E-02 -2.9681E-01
             1.8110E-01
 GRADIENT:   2.1117E+00  2.4083E+00  3.0708E+00 -4.6574E-01 -3.2350E+00 -4.8407E-01 -1.4211E+00  2.1697E-02 -1.0072E+00 -1.3751E+00
            -4.0606E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1729.41545598799        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      552
 NPARAMETR:  1.0321E+00  1.2507E+00  5.9202E-01  8.5537E-01  8.9701E-01  9.1155E-01  8.6578E-01  2.7225E-02  8.9508E-01  7.3160E-01
             1.0906E+00
 PARAMETER:  1.3156E-01  3.2373E-01 -4.2421E-01 -5.6218E-02 -8.6858E-03  7.3865E-03 -4.4123E-02 -3.5036E+00 -1.0843E-02 -2.1252E-01
             1.8670E-01
 GRADIENT:   1.5368E+00  4.9673E+00  1.4794E+00  3.8165E+00 -2.9065E+00 -1.1234E+00 -2.3209E-01  1.9860E-03 -3.7156E-01 -9.5482E-03
             7.6360E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1729.46488488323        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0317E+00  1.3322E+00  5.5621E-01  8.0065E-01  9.2166E-01  9.1453E-01  8.2429E-01  2.0336E-02  9.3813E-01  7.2930E-01
             1.0891E+00
 PARAMETER:  1.3116E-01  3.8685E-01 -4.8661E-01 -1.2233E-01  1.8421E-02  1.0659E-02 -9.3239E-02 -3.7954E+00  3.6134E-02 -2.1567E-01
             1.8533E-01
 GRADIENT:   1.8878E-01  4.7114E-01  1.3224E-01  2.2913E-01 -2.8271E-01  4.1876E-03 -2.4333E-02  1.4140E-03 -2.3789E-02 -1.5695E-02
            -1.2401E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1729.46530352653        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  1.0315E+00  1.3367E+00  5.5435E-01  7.9769E-01  9.2327E-01  9.1457E-01  8.2258E-01  1.1421E-02  9.4088E-01  7.3025E-01
             1.0892E+00
 PARAMETER:  1.3098E-01  3.9020E-01 -4.8996E-01 -1.2603E-01  2.0171E-02  1.0699E-02 -9.5305E-02 -4.3723E+00  3.9057E-02 -2.1437E-01
             1.8543E-01
 GRADIENT:  -2.9953E-01  1.0969E-01 -1.3129E-01  1.7276E-01  1.0967E-02  1.3961E-02  1.1768E-01  4.5837E-04  7.3119E-02  1.1255E-01
             8.0077E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1729.46551066897        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      999
 NPARAMETR:  1.0316E+00  1.3344E+00  5.5520E-01  7.9905E-01  9.2249E-01  9.1453E-01  8.2327E-01  1.0000E-02  9.3949E-01  7.2946E-01
             1.0891E+00
 PARAMETER:  1.3109E-01  3.8849E-01 -4.8842E-01 -1.2434E-01  1.9322E-02  1.0660E-02 -9.4472E-02 -4.7301E+00  3.7578E-02 -2.1545E-01
             1.8533E-01
 GRADIENT:   6.7277E-05  3.0748E-03 -8.3883E-03  1.4820E-02  1.2359E-02  2.2870E-03  1.9472E-03  0.0000E+00  3.0541E-05  1.0326E-03
            -1.6966E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      999
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0902E-05 -1.3561E-02 -3.4383E-04  7.9810E-03 -2.1425E-02
 SE:             2.9808E-02  2.3057E-02  1.5794E-04  2.3873E-02  2.0924E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9917E-01  5.5643E-01  2.9480E-02  7.3815E-01  3.0586E-01

 ETASHRINKSD(%)  1.4002E-01  2.2756E+01  9.9471E+01  2.0022E+01  2.9902E+01
 ETASHRINKVR(%)  2.7984E-01  4.0334E+01  9.9997E+01  3.6035E+01  5.0863E+01
 EBVSHRINKSD(%)  5.7327E-01  2.2619E+01  9.9507E+01  2.0466E+01  2.9868E+01
 EBVSHRINKVR(%)  1.1432E+00  4.0123E+01  9.9998E+01  3.6744E+01  5.0815E+01
 RELATIVEINF(%)  9.8777E+01  2.5844E+00  1.6207E-04  3.0297E+00  4.2144E+00
 EPSSHRINKSD(%)  4.2703E+01
 EPSSHRINKVR(%)  6.7171E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1729.4655106689718     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -994.31468410523360     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.10
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1729.466       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.33E+00  5.55E-01  7.99E-01  9.22E-01  9.15E-01  8.23E-01  1.00E-02  9.39E-01  7.29E-01  1.09E+00
 


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
+        1.23E+03
 
 TH 2
+       -1.05E+01  5.67E+02
 
 TH 3
+        1.38E+01  3.40E+02  9.32E+02
 
 TH 4
+       -2.45E+01  4.01E+02 -5.39E+02  1.26E+03
 
 TH 5
+       -6.05E+00 -4.57E+02 -9.24E+02  4.65E+02  1.21E+03
 
 TH 6
+        1.82E-01 -1.73E+00  3.23E+00 -3.50E+00 -2.92E+00  2.33E+02
 
 TH 7
+       -1.27E-01  2.47E+01 -2.59E+01 -1.27E+01 -2.00E+00 -1.01E+00  1.14E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.19E+00 -2.46E+01 -4.36E+01  5.35E+01 -7.32E+00 -1.43E+00  2.07E+01  0.00E+00  1.02E+02
 
 TH10
+       -1.37E+00 -1.20E+01 -5.84E+01 -2.22E+01 -6.74E+01  9.70E-02  3.01E+01  0.00E+00  1.28E+01  9.63E+01
 
 TH11
+       -8.64E+00 -1.95E+01 -3.46E+01 -3.93E+00 -6.52E+00  3.14E+00  9.39E+00  0.00E+00  1.34E+01  2.61E+01  1.83E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.670
Stop Time:
Sat Sep 18 14:17:33 CDT 2021
