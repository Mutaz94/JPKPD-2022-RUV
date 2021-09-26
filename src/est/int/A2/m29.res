Fri Sep 24 21:18:45 CDT 2021
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
$DATA ../../../../data/int/A2/dat29.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 NO. OF DATA RECS IN DATA SET:     1000
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

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2506.24757947040        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.5638E+01  4.8759E+01  1.0718E+02  4.4868E+00  1.4973E+02 -7.6177E+00 -1.8886E+02 -7.4775E+01 -3.0180E+01 -1.1280E+02
            -2.3106E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3139.52749157984        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9925E-01  7.2024E-01  7.3571E-01  1.1331E+00  6.9331E-01  9.8594E-01  1.7348E+00  6.8165E-01  9.2684E-01  9.9987E-01
             1.8116E+00
 PARAMETER:  9.9252E-02 -2.2817E-01 -2.0692E-01  2.2499E-01 -2.6628E-01  8.5839E-02  6.5088E-01 -2.8323E-01  2.4026E-02  9.9871E-02
             6.9420E-01
 GRADIENT:   3.9912E+01 -9.1792E+00  2.8147E+01  5.3147E+00  7.6924E+00 -9.9442E+00  1.9050E+01  1.4889E+01 -2.8247E+01  2.1034E+01
             1.0964E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3151.39711952219        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0006E+00  5.0680E-01  4.7082E-01  1.2975E+00  4.6309E-01  9.5351E-01  1.6968E+00  4.0708E-01  9.5595E-01  8.4305E-01
             1.7712E+00
 PARAMETER:  1.0065E-01 -5.7964E-01 -6.5327E-01  3.6044E-01 -6.6984E-01  5.2393E-02  6.2874E-01 -7.9874E-01  5.4951E-02 -7.0726E-02
             6.7168E-01
 GRADIENT:   4.3856E+01  1.4834E+01  9.4648E+00  2.1365E+02  3.2865E+01 -2.4773E+01  2.1237E+01  8.5861E+00 -5.4623E+01  2.8470E+01
             1.1818E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3167.18650257683        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.9983E-01  3.9653E-01  3.2065E-01  1.2896E+00  3.3037E-01  1.0023E+00  1.5507E+00  1.2306E-01  1.0620E+00  7.2454E-01
             1.7422E+00
 PARAMETER:  9.9829E-02 -8.2501E-01 -1.0374E+00  3.5435E-01 -1.0075E+00  1.0226E-01  5.3872E-01 -1.9951E+00  1.6016E-01 -2.2222E-01
             6.5518E-01
 GRADIENT:   3.9936E+01  5.7282E+01  1.7031E+01  2.0473E+02 -2.8129E+01 -4.6939E+00  1.4127E+00  3.5565E-01 -3.3767E+01  3.8549E+00
             2.2637E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3169.99680793516        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0002E+00  2.6445E-01  1.8405E-01  1.1610E+00  2.1556E-01  1.0248E+00  1.5661E+00  1.0771E-02  1.2099E+00  7.3720E-01
             1.6606E+00
 PARAMETER:  1.0024E-01 -1.2301E+00 -1.5926E+00  2.4928E-01 -1.4345E+00  1.2449E-01  5.4856E-01 -4.4309E+00  2.9051E-01 -2.0489E-01
             6.0719E-01
 GRADIENT:   3.8633E+01  4.6363E+01  2.7448E+01  1.7300E+02 -6.0374E+01  4.3509E+00 -2.9531E+00 -6.7364E-03 -2.3158E+01 -1.4992E+01
            -2.6127E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3183.50006214196        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      478
 NPARAMETR:  9.8673E-01  2.7655E-01  2.0886E-01  1.0727E+00  2.4308E-01  1.0167E+00  1.5729E+00  1.0000E-02  1.1956E+00  7.4940E-01
             1.6674E+00
 PARAMETER:  8.6645E-02 -1.1854E+00 -1.4661E+00  1.7022E-01 -1.3144E+00  1.1656E-01  5.5293E-01 -5.7729E+00  2.7867E-01 -1.8848E-01
             6.1125E-01
 GRADIENT:  -8.3571E-01 -8.0528E-01  2.9928E+00 -6.1125E-01 -4.4810E+00  1.0501E+00 -8.1568E-01  0.0000E+00  9.0698E-02 -1.6940E+00
            -1.8703E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3183.50416719018        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      647
 NPARAMETR:  9.8670E-01  2.7669E-01  2.0897E-01  1.0729E+00  2.4319E-01  1.0160E+00  1.5738E+00  1.0000E-02  1.1953E+00  7.5041E-01
             1.6675E+00
 PARAMETER:  8.6612E-02 -1.1849E+00 -1.4656E+00  1.7036E-01 -1.3139E+00  1.1585E-01  5.5346E-01 -5.7699E+00  2.7841E-01 -1.8714E-01
             6.1135E-01
 GRADIENT:  -9.1914E-01 -7.1086E-01  2.8420E+00 -5.7635E-01 -4.2482E+00  7.7911E-01 -6.6041E-01  0.0000E+00  4.7424E-02 -1.3171E+00
            -1.6428E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3183.50594460573        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:      750
 NPARAMETR:  9.8670E-01  2.7669E-01  2.0896E-01  1.0729E+00  2.4318E-01  1.0140E+00  1.5737E+00  1.0000E-02  1.1929E+00  7.5146E-01
             1.6675E+00
 PARAMETER:  8.6611E-02 -1.1848E+00 -1.4656E+00  1.7036E-01 -1.3140E+00  1.1393E-01  5.5346E-01 -5.7699E+00  2.7642E-01 -1.8574E-01
             6.1135E-01
 GRADIENT:   1.1527E+01  1.2718E+01  1.6021E+01  5.6224E+00  6.4728E+01  1.3273E+00  7.8151E-01  0.0000E+00  1.3460E+00 -5.0605E-01
            -4.8102E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -3183.50594460573        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  9.8667E-01  2.7679E-01  2.0887E-01  1.0728E+00  2.4309E-01  1.0139E+00  1.5735E+00  1.0000E-02  1.1928E+00  7.5142E-01
             1.6678E+00
 PARAMETER:  8.6611E-02 -1.1848E+00 -1.4656E+00  1.7036E-01 -1.3140E+00  1.1393E-01  5.5346E-01 -5.7699E+00  2.7642E-01 -1.8574E-01
             6.1135E-01
 GRADIENT:   4.5635E+04 -1.9244E+03  1.5585E+03  2.6792E+04  3.4791E+03  4.2923E-02  8.2448E+03  0.0000E+00  1.6497E+04  2.4572E+04
            -3.7349E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      815
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3792E-03  4.4840E-03 -1.5837E-05 -2.4541E-03  3.3145E-03
 SE:             2.9674E-02  2.6556E-02  3.1910E-04  2.8858E-02  2.7252E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6293E-01  8.6591E-01  9.6042E-01  9.3223E-01  9.0320E-01

 ETASHRINKSD(%)  5.8865E-01  1.1035E+01  9.8931E+01  3.3231E+00  8.7012E+00
 ETASHRINKVR(%)  1.1738E+00  2.0852E+01  9.9989E+01  6.5359E+00  1.6645E+01
 EBVSHRINKSD(%)  7.5238E-01  9.4269E+00  9.9048E+01  2.8107E+00  9.5201E+00
 EBVSHRINKVR(%)  1.4991E+00  1.7965E+01  9.9991E+01  5.5423E+00  1.8134E+01
 RELATIVEINF(%)  9.8500E+01  2.1993E+01  6.2458E-04  5.2969E+01  6.0963E+00
 EPSSHRINKSD(%)  2.1448E+01
 EPSSHRINKVR(%)  3.8295E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3183.5059446057271     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1529.4165848373163     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.81
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3183.506       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.87E-01  2.77E-01  2.09E-01  1.07E+00  2.43E-01  1.01E+00  1.57E+00  1.00E-02  1.19E+00  7.51E-01  1.67E+00
 


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
+        1.17E+08
 
 TH 2
+       -1.25E+02  1.06E+07
 
 TH 3
+        1.57E+02  2.03E+04  1.22E+07
 
 TH 4
+        2.40E+04 -7.26E+03  7.05E+03  3.42E+07
 
 TH 5
+        1.23E+02  2.00E+04 -2.81E+04  7.25E+03  1.13E+07
 
 TH 6
+        3.90E+02 -1.23E+02  1.39E+02  2.07E+02  1.06E+02  8.55E+07
 
 TH 7
+       -8.16E+02  2.67E+02 -2.37E+02  7.17E+06 -2.87E+02  4.46E+01  1.50E+06
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.99E+04  1.50E+04 -1.59E+04 -2.69E+04 -1.54E+04  1.16E+02 -5.65E+03  0.00E+00  1.05E+07
 
 TH10
+        2.50E+04 -7.54E+03  8.22E+03  1.35E+04  7.67E+03  2.73E+02  2.84E+03  0.00E+00 -3.53E+04  5.86E+07
 
 TH11
+       -5.02E+01  2.25E+03 -2.65E+03 -2.33E+03 -2.21E+03 -3.55E+01  8.86E+01  0.00E+00  4.84E+03 -2.40E+03  1.10E+06
 
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
 #CPUT: Total CPU Time in Seconds,       29.622
Stop Time:
Fri Sep 24 21:19:16 CDT 2021
