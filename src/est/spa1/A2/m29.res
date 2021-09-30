Wed Sep 29 23:13:24 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat29.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -996.186129839564        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7486E+02  4.7610E+01  2.1653E+01  5.9181E+01  1.8407E+02  2.7474E+01 -8.0303E+01 -7.5102E+00 -5.0013E+01 -8.7085E+01
            -1.9466E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1637.25085748585        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.7430E-01  7.5139E-01  1.1768E+00  1.2216E+00  8.8131E-01  9.9740E-01  1.6327E+00  5.2238E-01  1.1667E+00  9.4015E-01
             1.9200E+00
 PARAMETER:  7.3961E-02 -1.8583E-01  2.6277E-01  3.0017E-01 -2.6344E-02  9.7399E-02  5.9026E-01 -5.4937E-01  2.5418E-01  3.8282E-02
             7.5231E-01
 GRADIENT:   9.2209E+01  4.3048E+01  2.2206E+01  1.2070E+02  1.8673E+01  6.4866E+00  1.1919E+01  1.8845E+00  2.8590E+01 -2.8620E+01
            -2.4513E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1658.25650949487        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.9382E-01  6.2296E-01  3.6803E-01  1.1753E+00  4.4699E-01  8.7838E-01  1.5117E+00  1.7570E-02  1.0662E+00  7.3125E-01
             1.9243E+00
 PARAMETER:  9.3804E-02 -3.7326E-01 -8.9958E-01  2.6151E-01 -7.0522E-01 -2.9672E-02  5.1323E-01 -3.9415E+00  1.6409E-01 -2.1300E-01
             7.5458E-01
 GRADIENT:   1.0322E+02  5.5261E+00 -8.2529E+01  1.7041E+02  2.2725E+02 -5.4108E+01  5.4734E+00  9.7986E-03 -3.2800E-01  9.8894E+00
            -1.4805E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1684.05325553680        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      340
 NPARAMETR:  9.9561E-01  4.8151E-01  4.7041E-01  1.2778E+00  4.4411E-01  1.0502E+00  1.6127E+00  1.0000E-02  1.0071E+00  6.2535E-01
             2.3765E+00
 PARAMETER:  9.5598E-02 -6.3083E-01 -6.5416E-01  3.4512E-01 -7.1169E-01  1.4899E-01  5.7793E-01 -5.0100E+00  1.0704E-01 -3.6945E-01
             9.6563E-01
 GRADIENT:  -6.4587E-02  1.6386E+01 -2.0009E+01  7.4429E+01  3.3548E+01  1.6039E+01 -5.9809E+00  0.0000E+00 -1.2611E+01 -2.4484E+00
             2.6604E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1694.44912505118        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      517
 NPARAMETR:  9.8215E-01  1.5841E-01  4.9097E-01  1.3851E+00  3.7943E-01  9.8192E-01  2.7470E+00  1.0000E-02  9.9144E-01  7.4982E-01
             2.2879E+00
 PARAMETER:  8.1992E-02 -1.7426E+00 -6.1138E-01  4.2576E-01 -8.6908E-01  8.1758E-02  1.1105E+00 -7.1278E+00  9.1399E-02 -1.8792E-01
             9.2763E-01
 GRADIENT:  -4.0619E+00  8.6467E+00  5.9129E+01  2.2830E+01 -8.5843E+01 -4.0998E+00 -2.9053E+00  0.0000E+00 -3.4907E+00  5.7270E+00
             5.9038E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1698.61160043728        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      694
 NPARAMETR:  9.7987E-01  1.0382E-01  4.1572E-01  1.3668E+00  3.4303E-01  9.9402E-01  3.8464E+00  1.0000E-02  1.0150E+00  6.8254E-01
             2.2381E+00
 PARAMETER:  7.9666E-02 -2.1651E+00 -7.7775E-01  4.1248E-01 -9.6992E-01  9.4005E-02  1.4471E+00 -8.3817E+00  1.1486E-01 -2.8194E-01
             9.0564E-01
 GRADIENT:  -1.9963E+00  1.1864E+00  1.3513E+01  3.0850E+01 -2.6478E+01  5.2575E-01 -5.1763E+00  0.0000E+00  2.8426E-01  2.6285E+00
            -6.1354E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1698.96807399513        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  9.8065E-01  1.0234E-01  4.1494E-01  1.3485E+00  3.4485E-01  9.9218E-01  3.8339E+00  1.0000E-02  1.0170E+00  6.7921E-01
             2.2604E+00
 PARAMETER:  8.0458E-02 -2.1794E+00 -7.7963E-01  3.9902E-01 -9.6464E-01  9.2153E-02  1.4439E+00 -8.3817E+00  1.1689E-01 -2.8683E-01
             9.1554E-01
 GRADIENT:   2.1936E-01  1.5349E+00  1.5825E+00  1.4836E+00 -2.0541E+00 -3.3263E-02 -2.3372E+00  0.0000E+00  8.2924E-01  1.6007E+00
             2.9396E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1699.09365071583        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  9.8108E-01  9.5096E-02  4.1450E-01  1.3469E+00  3.4447E-01  9.8890E-01  3.9409E+00  1.0000E-02  1.0079E+00  6.6261E-01
             2.2707E+00
 PARAMETER:  8.0903E-02 -2.2529E+00 -7.8067E-01  3.9783E-01 -9.6574E-01  8.8833E-02  1.4714E+00 -8.3817E+00  1.0787E-01 -3.1157E-01
             9.2008E-01
 GRADIENT:   1.6641E+00 -2.0244E+00  4.3808E+00 -1.4367E+00 -2.5279E+00 -8.7247E-01 -7.7388E+00  0.0000E+00 -6.4734E-01  8.8086E-01
             6.7951E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1700.66853214012        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  9.7802E-01  3.3526E-02  4.0607E-01  1.3631E+00  3.3355E-01  9.8643E-01  6.8751E+00  1.0000E-02  1.0031E+00  6.5147E-01
             2.2721E+00
 PARAMETER:  7.7774E-02 -3.2954E+00 -8.0122E-01  4.0979E-01 -9.9796E-01  8.6336E-02  2.0279E+00 -8.3817E+00  1.0311E-01 -3.2853E-01
             9.2072E-01
 GRADIENT:   5.3973E+00  1.2555E+01 -4.7415E+00 -1.7162E+01  8.4306E+00 -3.5293E+00  2.5610E+01  0.0000E+00 -9.7484E+00 -8.3851E+00
             1.4596E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1700.82035517580        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1404
 NPARAMETR:  9.7581E-01  3.0250E-02  4.0484E-01  1.3650E+00  3.3322E-01  9.9021E-01  7.1617E+00  1.0000E-02  1.0108E+00  6.5284E-01
             2.2534E+00
 PARAMETER:  7.5333E-02 -3.3923E+00 -8.0569E-01  4.1038E-01 -9.9737E-01  8.9159E-02  2.0724E+00 -8.3817E+00  1.1051E-01 -3.2705E-01
             9.1076E-01
 GRADIENT:  -1.0837E+00  2.2825E+01 -1.0199E+02 -2.0369E+02  8.1253E+01 -5.9436E-01  7.8705E+01  0.0000E+00 -7.7751E+02 -2.6193E+02
            -9.3503E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1404
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4564E-03  1.5495E-02  3.0702E-05 -1.0553E-02 -2.2407E-03
 SE:             2.9382E-02  8.7434E-03  2.2816E-04  2.8072E-02  2.1894E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6047E-01  7.6353E-02  8.9295E-01  7.0698E-01  9.1848E-01

 ETASHRINKSD(%)  1.5679E+00  7.0709E+01  9.9236E+01  5.9543E+00  2.6653E+01
 ETASHRINKVR(%)  3.1113E+00  9.1420E+01  9.9994E+01  1.1554E+01  4.6201E+01
 EBVSHRINKSD(%)  1.6330E+00  7.9734E+01  9.9164E+01  4.9670E+00  2.5429E+01
 EBVSHRINKVR(%)  3.2394E+00  9.5893E+01  9.9993E+01  9.6872E+00  4.4391E+01
 RELATIVEINF(%)  9.6147E+01  2.9823E+00  3.8856E-04  3.9369E+01  3.1488E+00
 EPSSHRINKSD(%)  2.7682E+01
 EPSSHRINKVR(%)  4.7701E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1700.8203551757963     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -781.88182197112360     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.43
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1700.820       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  3.04E-02  4.04E-01  1.36E+00  3.34E-01  9.89E-01  7.19E+00  1.00E-02  1.01E+00  6.52E-01  2.25E+00
 


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
+        1.22E+03
 
 TH 2
+       -5.24E+01  4.04E+05
 
 TH 3
+       -3.72E+00 -6.42E+04  2.39E+04
 
 TH 4
+        4.78E+00 -3.82E+04 -4.68E+02  7.28E+03
 
 TH 5
+        6.46E+01  6.01E+04 -6.54E+03 -6.43E+01  5.13E+04
 
 TH 6
+       -5.33E+01 -3.18E+00 -1.63E+01 -2.12E+01  1.13E+01  2.37E+02
 
 TH 7
+       -8.89E-01  1.43E+03 -4.50E+02 -2.67E+02  4.26E+02  6.81E-01  1.99E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.97E+05  3.35E+02 -4.31E+01 -7.34E+01  1.01E+02  1.94E+05  9.40E-01  0.00E+00  1.72E+05
 
 TH10
+        5.81E+01 -2.51E+02 -2.12E+02 -9.70E+01  1.29E+02 -5.18E+01  4.27E+00  0.00E+00 -2.04E+02  4.71E+04
 
 TH11
+       -9.40E+00 -1.06E+04 -4.84E+01  1.84E+03  2.27E+01 -1.07E+00 -7.36E+01  0.00E+00 -1.44E+01 -1.82E-01  5.99E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.843
Stop Time:
Wed Sep 29 23:14:01 CDT 2021
