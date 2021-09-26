Sat Sep 25 08:02:24 CDT 2021
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
$DATA ../../../../data/spa/A1/dat40.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1333.66219875379        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4565E+01  1.9401E+01  1.0240E+01  1.4513E+01  4.9588E+01  6.0051E+01 -2.4734E+01 -9.0633E-01 -3.2843E+01 -3.5097E+01
            -6.3007E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1507.31493418258        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0631E+00  9.9330E-01  1.0216E+00  1.0237E+00  9.5748E-01  8.3417E-01  1.0558E+00  9.6738E-01  1.0849E+00  9.7322E-01
             1.8637E+00
 PARAMETER:  1.6123E-01  9.3276E-02  1.2136E-01  1.2339E-01  5.6554E-02 -8.1319E-02  1.5434E-01  6.6837E-02  1.8145E-01  7.2854E-02
             7.2256E-01
 GRADIENT:   1.7204E+02  2.0207E+01  9.0514E+00  1.6956E+01 -1.5663E+01 -6.4009E+00 -2.1801E+00  2.9978E+00  1.2654E-01  3.2568E+00
            -3.0181E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1513.42120334965        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  1.0417E+00  7.4439E-01  5.6508E-01  1.1567E+00  6.2900E-01  8.5721E-01  1.5684E+00  3.2893E-01  8.6061E-01  5.1263E-01
             1.8825E+00
 PARAMETER:  1.4090E-01 -1.9519E-01 -4.7078E-01  2.4559E-01 -3.6363E-01 -5.4068E-02  5.5006E-01 -1.0119E+00 -5.0115E-02 -5.6820E-01
             7.3261E-01
 GRADIENT:   8.4841E+01  1.7276E+01 -5.4846E+01  1.2044E+02  8.2466E+01  6.4478E+00  1.5568E+00  5.9203E-01 -1.1221E+01 -4.6072E+00
            -2.0910E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1520.35137003311        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  1.0126E+00  5.2782E-01  4.1952E-01  1.1695E+00  4.4835E-01  8.4738E-01  1.8281E+00  9.5535E-02  8.3847E-01  4.5695E-01
             1.8795E+00
 PARAMETER:  1.1248E-01 -5.3900E-01 -7.6864E-01  2.5654E-01 -7.0218E-01 -6.5606E-02  7.0327E-01 -2.2483E+00 -7.6178E-02 -6.8318E-01
             7.3101E-01
 GRADIENT:  -3.0517E+00  1.3947E+01  1.1048E+01  2.5815E+01 -1.9580E+01  2.4030E+00 -6.2367E-01  3.5769E-02 -2.8382E+00 -2.8073E+00
             3.4484E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1520.95657017316        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0154E+00  5.0678E-01  4.4430E-01  1.1808E+00  4.5953E-01  8.4519E-01  1.8744E+00  7.7505E-02  8.4602E-01  5.1554E-01
             1.8619E+00
 PARAMETER:  1.1524E-01 -5.7968E-01 -7.1126E-01  2.6617E-01 -6.7756E-01 -6.8198E-02  7.2830E-01 -2.4574E+00 -6.7217E-02 -5.6255E-01
             7.2160E-01
 GRADIENT:  -4.5126E+00  6.9128E+00  1.0197E+01 -1.6055E+00 -2.0651E+01  1.3114E+00 -1.8253E+00  5.6016E-02 -4.8556E-01 -9.5863E-02
             8.1087E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1521.74278179505        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  1.0126E+00  3.8555E-01  5.4180E-01  1.2729E+00  4.9583E-01  8.3707E-01  2.3393E+00  3.6666E-02  8.3398E-01  6.1834E-01
             1.8779E+00
 PARAMETER:  1.1257E-01 -8.5307E-01 -5.1287E-01  3.4127E-01 -6.0152E-01 -7.7843E-02  9.4983E-01 -3.2059E+00 -8.1540E-02 -3.8071E-01
             7.3015E-01
 GRADIENT:  -5.6553E-01  6.8715E-01 -3.2196E-01 -1.6131E-02 -3.5212E-01  1.6065E-01  6.0910E-01  1.3685E-02 -1.9227E-01  2.5424E-01
             1.0664E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1521.83435313193        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  1.0088E+00  2.8807E-01  6.0072E-01  1.3448E+00  5.0839E-01  8.3807E-01  2.7090E+00  1.3761E-02  8.3368E-01  6.9351E-01
             1.8822E+00
 PARAMETER:  1.0881E-01 -1.1446E+00 -4.0962E-01  3.9626E-01 -5.7650E-01 -7.6649E-02  1.0966E+00 -4.1859E+00 -8.1900E-02 -2.6600E-01
             7.3243E-01
 GRADIENT:  -1.3523E+00  2.4252E+00  3.1926E+00  8.2040E+00 -6.9142E+00  1.8118E+00 -2.2246E-02  2.1143E-03 -1.3921E+00  1.1407E+00
             1.4242E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1522.03884263324        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  1.0063E+00  2.0168E-01  5.9794E-01  1.3830E+00  4.9219E-01  8.3368E-01  3.2600E+00  1.0000E-02  8.2608E-01  7.0776E-01
             1.8670E+00
 PARAMETER:  1.0633E-01 -1.5011E+00 -4.1427E-01  4.2423E-01 -6.0889E-01 -8.1906E-02  1.2817E+00 -5.8929E+00 -9.1058E-02 -2.4565E-01
             7.2431E-01
 GRADIENT:  -1.4081E-01  1.0233E-01  7.9549E-01  6.5854E+00 -2.0746E+00  3.5120E-01 -1.1044E+00  0.0000E+00  3.3111E-02  1.4067E-01
            -1.4326E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1522.05520347907        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  1.0059E+00  1.9117E-01  5.9923E-01  1.3851E+00  4.9157E-01  8.3249E-01  3.3728E+00  1.0000E-02  8.2356E-01  7.0866E-01
             1.8725E+00
 PARAMETER:  1.0593E-01 -1.5546E+00 -4.1212E-01  4.2575E-01 -6.1016E-01 -8.3331E-02  1.3158E+00 -6.1914E+00 -9.4116E-02 -2.4438E-01
             7.2728E-01
 GRADIENT:   6.5518E-03  4.8431E-03  2.6501E-02  1.2776E-02 -4.5759E-02 -5.3707E-03  1.7122E-03  0.0000E+00 -4.7784E-03  6.4887E-03
             2.6292E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1522.05520347907        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1082
 NPARAMETR:  1.0059E+00  1.9117E-01  5.9923E-01  1.3851E+00  4.9157E-01  8.3249E-01  3.3728E+00  1.0000E-02  8.2356E-01  7.0866E-01
             1.8725E+00
 PARAMETER:  1.0593E-01 -1.5546E+00 -4.1212E-01  4.2575E-01 -6.1016E-01 -8.3331E-02  1.3158E+00 -6.1914E+00 -9.4116E-02 -2.4438E-01
             7.2728E-01
 GRADIENT:   6.5518E-03  4.8431E-03  2.6501E-02  1.2776E-02 -4.5759E-02 -5.3707E-03  1.7122E-03  0.0000E+00 -4.7784E-03  6.4887E-03
             2.6292E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1082
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5033E-03  4.2544E-02 -2.8306E-04 -2.8034E-02  8.1288E-03
 SE:             2.9346E-02  1.7430E-02  2.0548E-04  2.5645E-02  2.0094E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5914E-01  1.4652E-02  1.6835E-01  2.7432E-01  6.8582E-01

 ETASHRINKSD(%)  1.6884E+00  4.1608E+01  9.9312E+01  1.4088E+01  3.2682E+01
 ETASHRINKVR(%)  3.3483E+00  6.5904E+01  9.9995E+01  2.6191E+01  5.4683E+01
 EBVSHRINKSD(%)  1.9502E+00  5.5054E+01  9.9229E+01  1.0681E+01  2.6988E+01
 EBVSHRINKVR(%)  3.8623E+00  7.9799E+01  9.9994E+01  2.0220E+01  4.6692E+01
 RELATIVEINF(%)  9.4801E+01  6.9163E+00  3.0816E-04  3.1728E+01  2.7075E+00
 EPSSHRINKSD(%)  3.7775E+01
 EPSSHRINKVR(%)  6.1281E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1522.0552034790719     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -786.90437691533373     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.53
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1522.055       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.91E-01  5.99E-01  1.39E+00  4.92E-01  8.32E-01  3.37E+00  1.00E-02  8.24E-01  7.09E-01  1.87E+00
 


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
+        1.53E+03
 
 TH 2
+       -4.16E+01  1.21E+03
 
 TH 3
+        6.51E+00  2.15E+01  2.12E+03
 
 TH 4
+       -4.04E+01  1.47E+02 -1.52E+02  7.36E+02
 
 TH 5
+        4.99E+01 -4.00E+02 -3.23E+03 -6.45E+01  5.30E+03
 
 TH 6
+       -2.58E+00 -1.01E+01  5.33E+00 -7.50E+00  2.05E+00  2.65E+02
 
 TH 7
+        4.00E+00  9.95E+01 -3.89E+01 -1.90E+01  4.68E+01  3.47E-01  1.13E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.41E+00 -1.90E+02  6.88E+01  1.56E+01 -4.58E+01 -4.26E+00 -1.32E+01  0.00E+00  2.67E+02
 
 TH10
+       -6.38E+00 -2.77E+01 -9.47E-02  1.39E+00 -6.63E+01  4.98E+00 -6.74E+00  0.00E+00 -4.62E+00  1.33E+02
 
 TH11
+       -1.50E+01 -7.18E+00 -2.72E+01 -7.05E+00  6.60E+00  3.47E+00 -5.55E-01  0.00E+00  5.11E+00  3.11E+01  7.32E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.744
Stop Time:
Sat Sep 25 08:02:46 CDT 2021
