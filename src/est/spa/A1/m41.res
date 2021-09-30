Wed Sep 29 12:07:11 CDT 2021
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
$DATA ../../../../data/spa/A1/dat41.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1380.53691913981        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6761E+02 -3.0914E+01 -3.3872E+01  1.0248E+01  1.2932E+02  5.0110E+01 -9.9493E+00  9.4015E-02 -2.3202E+01 -2.0536E+01
            -5.3903E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1518.61219370171        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0999E+00  1.0192E+00  1.0163E+00  1.0497E+00  9.1586E-01  9.8297E-01  1.0164E+00  9.7067E-01  1.0380E+00  9.7390E-01
             1.8744E+00
 PARAMETER:  1.9522E-01  1.1903E-01  1.1620E-01  1.4850E-01  1.2108E-02  8.2827E-02  1.1627E-01  7.0233E-02  1.3726E-01  7.3558E-02
             7.2827E-01
 GRADIENT:   3.7623E+02  1.9625E+01  5.7404E+00  2.7856E+01 -8.1120E+00  1.1980E+01 -1.0187E-01  3.8156E+00  2.4130E+00  2.4305E+00
             1.2765E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1522.88187305980        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      206
 NPARAMETR:  1.0791E+00  8.7549E-01  8.6958E-01  1.1304E+00  8.0353E-01  9.7415E-01  1.1324E+00  4.8859E-01  9.8281E-01  9.4574E-01
             1.8263E+00
 PARAMETER:  1.7614E-01 -3.2968E-02 -3.9746E-02  2.2253E-01 -1.1874E-01  7.3815E-02  2.2435E-01 -6.1623E-01  8.2661E-02  4.4209E-02
             7.0230E-01
 GRADIENT:   9.6305E+01  6.2745E+00 -2.2692E+00  6.1547E+00  1.4146E+00  4.1799E+00  6.4927E-02  7.6776E-01 -1.1779E+00  6.0201E+00
            -4.7889E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1526.17364760769        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      383
 NPARAMETR:  1.0213E+00  6.9102E-01  7.4686E-01  1.2247E+00  6.6879E-01  9.4872E-01  1.3009E+00  2.6717E-01  9.1786E-01  7.9272E-01
             1.8680E+00
 PARAMETER:  1.2103E-01 -2.6959E-01 -1.9188E-01  3.0271E-01 -3.0228E-01  4.7356E-02  3.6306E-01 -1.2199E+00  1.4288E-02 -1.3228E-01
             7.2488E-01
 GRADIENT:  -3.3893E+01  7.3029E+00 -5.7303E+00  1.5425E+01  4.7384E+00 -2.1680E+00 -7.0398E-01  3.0654E-01 -9.3977E-01 -8.2781E-01
             5.3194E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1527.39428887170        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      561
 NPARAMETR:  1.0339E+00  4.8964E-01  7.3569E-01  1.3280E+00  6.0020E-01  9.4929E-01  1.5935E+00  1.1969E-01  8.6610E-01  7.8769E-01
             1.8590E+00
 PARAMETER:  1.3331E-01 -6.1408E-01 -2.0694E-01  3.8370E-01 -4.1049E-01  4.7956E-02  5.6593E-01 -2.0228E+00 -4.3750E-02 -1.3865E-01
             7.2006E-01
 GRADIENT:   8.3919E-01  4.4923E+00  4.6635E+00  5.5036E+00 -8.2918E+00 -7.2793E-01 -9.8933E-02 -3.2224E-03 -1.0002E+00 -1.2573E+00
            -8.2771E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1527.75219413681        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      736
 NPARAMETR:  1.0279E+00  2.8611E-01  7.9037E-01  1.4487E+00  5.8514E-01  9.4708E-01  1.9728E+00  2.0451E-02  8.2666E-01  8.4858E-01
             1.8688E+00
 PARAMETER:  1.2756E-01 -1.1514E+00 -1.3525E-01  4.7067E-01 -4.3590E-01  4.5624E-02  7.7944E-01 -3.7897E+00 -9.0360E-02 -6.4186E-02
             7.2530E-01
 GRADIENT:  -1.7121E+00  1.7931E+00  1.2066E+00  1.0248E+01 -2.7721E+00 -1.1488E-01 -4.1563E-01 -1.0920E-03 -6.3527E-01  3.6859E-01
            -1.1052E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1527.80243025952        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  1.0267E+00  2.1077E-01  8.0221E-01  1.4873E+00  5.7668E-01  9.4553E-01  2.3234E+00  1.0000E-02  8.1136E-01  8.5467E-01
             1.8771E+00
 PARAMETER:  1.2632E-01 -1.4570E+00 -1.2038E-01  4.9699E-01 -4.5048E-01  4.3986E-02  9.4302E-01 -4.8704E+00 -1.0904E-01 -5.7042E-02
             7.2971E-01
 GRADIENT:   2.4793E-01  6.4363E-01  1.1799E+00  4.7529E+00 -1.5305E+00 -8.8991E-02 -2.3808E-01  0.0000E+00 -7.1779E-02  4.4160E-01
            -2.0422E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1527.81617508620        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1095             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0260E+00  1.8023E-01  7.9605E-01  1.4984E+00  5.6858E-01  9.4522E-01  2.5520E+00  1.0000E-02  8.0551E-01  8.4991E-01
             1.8793E+00
 PARAMETER:  1.2567E-01 -1.6135E+00 -1.2810E-01  5.0437E-01 -4.6462E-01  4.3662E-02  1.0369E+00 -5.4022E+00 -1.1628E-01 -6.2622E-02
             7.3092E-01
 GRADIENT:   1.4089E+02  6.8816E+00  1.2901E+00  2.1172E+02  2.6194E+01  8.7186E+00  1.5052E+00  0.0000E+00  4.6143E+00  8.5819E-01
             5.8571E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1527.81691301708        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1264
 NPARAMETR:  1.0259E+00  1.8056E-01  7.9579E-01  1.4983E+00  5.6833E-01  9.4516E-01  2.5517E+00  1.0000E-02  8.0511E-01  8.4954E-01
             1.8793E+00
 PARAMETER:  1.2553E-01 -1.6117E+00 -1.2843E-01  5.0434E-01 -4.6505E-01  4.3595E-02  1.0368E+00 -5.4022E+00 -1.1678E-01 -6.3064E-02
             7.3092E-01
 GRADIENT:   2.8398E-01 -1.4744E-02 -2.8806E-01 -1.7724E+00  1.1366E+00  1.7406E-03 -6.9859E-02  0.0000E+00 -3.4549E-02 -2.2748E-02
            -2.5307E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1264
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.2086E-05  2.1746E-03 -6.6991E-05 -9.8268E-03 -1.6810E-02
 SE:             2.9444E-02  7.6266E-03  2.0325E-04  2.7261E-02  2.3025E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9750E-01  7.7554E-01  7.4170E-01  7.1849E-01  4.6534E-01

 ETASHRINKSD(%)  1.3592E+00  7.4450E+01  9.9319E+01  8.6738E+00  2.2864E+01
 ETASHRINKVR(%)  2.6999E+00  9.3472E+01  9.9995E+01  1.6595E+01  4.0500E+01
 EBVSHRINKSD(%)  1.4975E+00  7.5887E+01  9.9332E+01  8.2763E+00  2.1658E+01
 EBVSHRINKVR(%)  2.9727E+00  9.4186E+01  9.9996E+01  1.5868E+01  3.8626E+01
 RELATIVEINF(%)  8.9178E+01  1.9167E-01  2.5061E-04  4.6039E+00  2.0455E+00
 EPSSHRINKSD(%)  3.6767E+01
 EPSSHRINKVR(%)  6.0015E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1527.8169130170807     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -792.66608645334247     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.77
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1527.817       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.81E-01  7.96E-01  1.50E+00  5.68E-01  9.45E-01  2.55E+00  1.00E-02  8.05E-01  8.50E-01  1.88E+00
 


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
+        1.15E+03
 
 TH 2
+       -5.07E+01  3.40E+02
 
 TH 3
+        1.30E+01  2.23E+02  7.95E+02
 
 TH 4
+       -2.90E+01  3.38E+02 -1.02E+02  6.45E+02
 
 TH 5
+        1.50E+01 -5.88E+02 -1.45E+03 -5.97E+01  2.93E+03
 
 TH 6
+       -1.65E-01 -6.73E+00  5.05E+00 -7.34E+00 -2.41E+00  2.09E+02
 
 TH 7
+        8.49E-02  5.87E+00  4.65E-01 -5.83E-01 -8.98E-01  7.62E-02  8.95E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.13E+00 -2.82E+01  1.35E+01 -3.65E+00 -9.05E+00  1.74E+00  2.05E+00  0.00E+00  2.17E+02
 
 TH10
+       -3.25E+00  2.04E+01 -1.20E+01 -4.33E+00 -6.92E+01 -2.20E-01  8.51E-01  0.00E+00  1.23E+00  1.07E+02
 
 TH11
+       -1.17E+01 -2.54E+00 -1.80E+01 -9.99E+00  8.86E+00  2.38E+00  1.79E-01  0.00E+00  1.16E+01  2.35E+01  7.30E+01
 
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
 #CPUT: Total CPU Time in Seconds,       23.533
Stop Time:
Wed Sep 29 12:07:37 CDT 2021
