Sat Sep 25 07:50:28 CDT 2021
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
$DATA ../../../../data/spa/A1/dat5.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m5.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1588.07027269619        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6152E+02 -3.3829E+01  7.5567E+00 -4.8170E+01 -3.2053E+01  1.2722E+01 -1.0873E+01  5.9153E+00  1.4534E+01  1.8761E+00
            -9.5675E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1602.09906595187        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9527E-01  1.1286E+00  1.0733E+00  9.3742E-01  1.1279E+00  9.6857E-01  1.1392E+00  9.3774E-01  8.7427E-01  9.8369E-01
             1.2347E+00
 PARAMETER:  9.5262E-02  2.2102E-01  1.7076E-01  3.5371E-02  2.2036E-01  6.8069E-02  2.3032E-01  3.5713E-02 -3.4370E-02  8.3558E-02
             3.1083E-01
 GRADIENT:   1.4113E+02 -1.7222E+01  1.0515E+01 -3.1437E+01  2.8231E-01  2.4157E+00 -2.9862E+00 -7.3975E-01  1.8913E-01 -6.8766E+00
             5.3997E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1606.07966473404        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.5720E-01  1.2505E+00  6.9005E-01  8.6066E-01  9.8099E-01  9.5185E-01  1.1207E+00  3.8584E-01  8.8170E-01  9.4863E-01
             1.2087E+00
 PARAMETER:  5.6254E-02  3.2352E-01 -2.7100E-01 -5.0054E-02  8.0809E-02  5.0648E-02  2.1391E-01 -8.5233E-01 -2.5903E-02  4.7259E-02
             2.8958E-01
 GRADIENT:   4.3420E+01  4.6772E+00 -1.0716E+01  3.4294E+00 -5.4535E-01  1.2542E-01  4.3884E+00  1.2569E+00  2.6775E+00  1.0907E+01
             6.6817E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1606.78440648933        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.3889E-01  1.0647E+00  7.6228E-01  9.7873E-01  9.1788E-01  9.4922E-01  1.2570E+00  3.9989E-01  8.0932E-01  8.4904E-01
             1.1970E+00
 PARAMETER:  3.6939E-02  1.6265E-01 -1.7145E-01  7.8500E-02  1.4308E-02  4.7882E-02  3.2873E-01 -8.1656E-01 -1.1156E-01 -6.3643E-02
             2.7980E-01
 GRADIENT:  -6.9585E-01  4.7120E+00 -1.1234E+00  7.1676E+00  1.7046E-02 -8.2970E-01 -4.4585E-01  7.7905E-01  5.3343E-01  6.1263E-01
            -3.5887E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1607.10840990077        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  9.4219E-01  9.6329E-01  6.8982E-01  1.0255E+00  8.2141E-01  9.4866E-01  1.3661E+00  1.3141E-01  7.6352E-01  7.5599E-01
             1.1974E+00
 PARAMETER:  4.0456E-02  6.2602E-02 -2.7132E-01  1.2523E-01 -9.6728E-02  4.7296E-02  4.1198E-01 -1.9294E+00 -1.6982E-01 -1.7973E-01
             2.8012E-01
 GRADIENT:   6.4281E+00 -1.3220E+00 -1.9145E+00  5.5411E+00  2.9131E+00 -1.1721E+00 -6.9356E-01  1.2031E-01 -4.8296E-01 -4.1114E-01
             2.8444E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1607.19960859700        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      425
 NPARAMETR:  9.4458E-01  9.7328E-01  7.0573E-01  1.0205E+00  8.3637E-01  9.5336E-01  1.3606E+00  1.3008E-01  7.7122E-01  7.7305E-01
             1.1969E+00
 PARAMETER:  4.2985E-02  7.2921E-02 -2.4852E-01  1.2028E-01 -7.8682E-02  5.2233E-02  4.0791E-01 -1.9396E+00 -1.5979E-01 -1.5741E-01
             2.7974E-01
 GRADIENT:  -1.2569E+01 -2.9760E+00  1.1585E+00 -4.4322E+00 -1.0814E+00 -1.2305E+00 -1.1635E+00  8.2168E-02 -3.3128E-01 -7.5910E-01
            -3.6238E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1607.29154215965        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      600
 NPARAMETR:  9.4938E-01  1.0219E+00  7.0240E-01  9.9477E-01  8.5838E-01  9.5588E-01  1.3139E+00  6.7143E-02  7.8844E-01  7.9283E-01
             1.1975E+00
 PARAMETER:  4.8056E-02  1.2164E-01 -2.5325E-01  9.4756E-02 -5.2713E-02  5.4880E-02  3.7300E-01 -2.6009E+00 -1.3770E-01 -1.3215E-01
             2.8021E-01
 GRADIENT:  -8.1576E-01 -7.6509E-01 -3.3240E-01 -7.6449E-01  5.0134E-01 -2.0549E-01 -1.7135E-01  2.1764E-02  5.9650E-02 -1.4710E-01
            -1.0609E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1607.30336658071        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      776
 NPARAMETR:  9.4970E-01  1.0135E+00  7.0452E-01  1.0005E+00  8.5534E-01  9.5634E-01  1.3248E+00  1.0648E-02  7.8482E-01  7.9399E-01
             1.1977E+00
 PARAMETER:  4.8388E-02  1.1344E-01 -2.5024E-01  1.0050E-01 -5.6253E-02  5.5354E-02  3.8128E-01 -4.4424E+00 -1.4230E-01 -1.3068E-01
             2.8038E-01
 GRADIENT:   6.1824E-02  9.5228E-02  1.7733E-02  1.7130E-01 -2.8010E-02 -1.2581E-02 -1.4150E-02  5.1930E-04 -1.0926E-02 -1.4232E-02
             7.5734E-03

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1607.30342112025        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  9.4968E-01  1.0147E+00  7.0443E-01  9.9971E-01  8.5593E-01  9.5637E-01  1.3235E+00  1.0000E-02  7.8531E-01  7.9442E-01
             1.1976E+00
 PARAMETER:  4.8365E-02  1.1462E-01 -2.5036E-01  9.9709E-02 -5.5562E-02  5.5389E-02  3.8031E-01 -4.7404E+00 -1.4167E-01 -1.3014E-01
             2.8035E-01
 GRADIENT:   6.5098E-04  7.6964E-04 -6.1002E-04  1.6996E-03  1.3093E-04  4.2579E-04  1.0971E-04  0.0000E+00  3.7263E-05  1.2407E-04
             1.1381E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      903
 NO. OF SIG. DIGITS IN FINAL EST.:  4.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.5733E-05  3.9552E-03 -4.0795E-04 -9.8838E-03 -1.1519E-02
 SE:             2.9738E-02  2.3495E-02  1.7665E-04  2.2345E-02  2.1033E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9824E-01  8.6631E-01  2.0921E-02  6.5826E-01  5.8394E-01

 ETASHRINKSD(%)  3.7265E-01  2.1288E+01  9.9408E+01  2.5141E+01  2.9535E+01
 ETASHRINKVR(%)  7.4391E-01  3.8044E+01  9.9996E+01  4.3961E+01  5.0347E+01
 EBVSHRINKSD(%)  6.6446E-01  2.0593E+01  9.9443E+01  2.5554E+01  2.8788E+01
 EBVSHRINKVR(%)  1.3245E+00  3.6946E+01  9.9997E+01  4.4579E+01  4.9289E+01
 RELATIVEINF(%)  9.8350E+01  5.1081E+00  2.9580E-04  4.5542E+00  4.1450E+00
 EPSSHRINKSD(%)  4.1315E+01
 EPSSHRINKVR(%)  6.5561E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1607.3034211202480     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -872.15259455650983     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.46
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1607.303       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.50E-01  1.01E+00  7.04E-01  1.00E+00  8.56E-01  9.56E-01  1.32E+00  1.00E-02  7.85E-01  7.94E-01  1.20E+00
 


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
+       -8.92E+00  3.98E+02
 
 TH 3
+        2.18E+01  2.06E+02  7.74E+02
 
 TH 4
+       -2.04E+01  3.29E+02 -4.21E+02  1.01E+03
 
 TH 5
+       -8.69E+00 -3.26E+02 -8.55E+02  4.51E+02  1.21E+03
 
 TH 6
+        3.75E-01 -2.66E-01  6.01E+00 -6.19E+00 -3.39E+00  2.11E+02
 
 TH 7
+        1.56E+00  3.06E+01 -1.52E+01 -1.78E+01  2.74E+00 -4.02E-01  4.80E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.37E+00 -2.08E+01 -3.60E+01  1.91E+01  9.76E+00 -2.57E+00  1.87E+01  0.00E+00  1.12E+02
 
 TH10
+       -1.41E+00 -8.71E+00 -5.38E+01 -2.47E+01 -6.40E+01 -8.36E-03  9.62E+00  0.00E+00  1.29E+01  8.71E+01
 
 TH11
+       -9.62E+00 -1.34E+01 -3.67E+01 -1.55E+00  4.81E+00  3.12E+00  5.53E+00  0.00E+00  1.42E+01  2.41E+01  1.57E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.249
Stop Time:
Sat Sep 25 07:50:44 CDT 2021
