Sat Sep 25 10:01:44 CDT 2021
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
$DATA ../../../../data/spa/S1/dat66.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1626.32669588902        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4666E+02 -7.4006E+01 -1.8114E+00 -8.3307E+01  9.4714E+00  1.0296E+01 -9.1618E+00  3.6679E+00  2.2347E+01 -1.2128E+01
            -1.6119E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1636.13152687923        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.6717E-01  1.1704E+00  9.8016E-01  9.5034E-01  1.0942E+00  9.4537E-01  1.1472E+00  9.4416E-01  7.3431E-01  1.1526E+00
             1.0661E+00
 PARAMETER:  6.6619E-02  2.5733E-01  7.9957E-02  4.9062E-02  1.9007E-01  4.3820E-02  2.3732E-01  4.2536E-02 -2.0883E-01  2.4201E-01
             1.6400E-01
 GRADIENT:   6.8943E+01  1.1949E+00 -1.6255E+00 -9.6565E+00  1.2020E+01 -6.7140E+00 -1.1753E+00  5.6747E-01 -7.4596E+00  5.3948E-01
             6.0790E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1637.15721617660        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.6786E-01  9.3161E-01  9.7126E-01  1.1096E+00  9.5234E-01  9.5815E-01  1.5146E+00  6.7563E-01  6.4125E-01  1.0357E+00
             1.1046E+00
 PARAMETER:  6.7333E-02  2.9158E-02  7.0843E-02  2.0403E-01  5.1171E-02  5.7249E-02  5.1515E-01 -2.9211E-01 -3.4434E-01  1.3512E-01
             1.9947E-01
 GRADIENT:   6.8692E+01  2.1299E+01  3.9800E+00  1.7906E+01  1.7693E+00 -1.1143E+00  1.0533E+01 -3.3313E-01 -6.0336E+00 -2.0112E+00
             1.8899E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1639.06010430636        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.3942E-01  9.4422E-01  7.8846E-01  1.0753E+00  8.5270E-01  9.5838E-01  1.3393E+00  4.0187E-01  7.3787E-01  9.4201E-01
             1.0437E+00
 PARAMETER:  3.7510E-02  4.2609E-02 -1.3767E-01  1.7261E-01 -5.9348E-02  5.7493E-02  3.9215E-01 -8.1162E-01 -2.0398E-01  4.0259E-02
             1.4276E-01
 GRADIENT:  -2.9849E+00  4.9457E-01 -3.9544E+00  2.8902E+00  3.4661E+00 -1.3086E+00 -1.0541E-01  9.0613E-01  4.3444E-01  1.2482E+00
             1.2733E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1640.37696181787        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.4677E-01  7.4826E-01  7.0656E-01  1.1774E+00  7.1474E-01  9.5453E-01  1.6346E+00  1.1636E-01  6.6536E-01  8.1317E-01
             1.0611E+00
 PARAMETER:  4.5303E-02 -1.9000E-01 -2.4735E-01  2.6331E-01 -2.3584E-01  5.3466E-02  5.9139E-01 -2.0511E+00 -3.0743E-01 -1.0682E-01
             1.5931E-01
 GRADIENT:   1.3714E+01  6.5084E+00  7.6399E-01  1.4495E+01 -1.4087E+00 -2.5271E+00  4.9232E-01  8.8859E-02 -3.3822E+00 -1.0706E+00
             6.3171E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1641.16481576700        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      458
 NPARAMETR:  9.5411E-01  6.3877E-01  7.8295E-01  1.2505E+00  7.2463E-01  9.7223E-01  1.8901E+00  5.6451E-02  6.5490E-01  8.8401E-01
             1.0406E+00
 PARAMETER:  5.3027E-02 -3.4821E-01 -1.4469E-01  3.2352E-01 -2.2209E-01  7.1832E-02  7.3663E-01 -2.7744E+00 -3.2328E-01 -2.3291E-02
             1.3976E-01
 GRADIENT:   6.1154E+00  5.8984E+00  9.1755E+00  8.4134E-01 -1.0207E+01  2.3886E+00  7.8057E-01 -5.4719E-03 -7.8506E-01 -4.8349E-01
            -3.3036E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1641.40334578663        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      633
 NPARAMETR:  9.5081E-01  5.4780E-01  7.6346E-01  1.2925E+00  6.8838E-01  9.6519E-01  2.1062E+00  1.8866E-02  6.4003E-01  8.5876E-01
             1.0476E+00
 PARAMETER:  4.9559E-02 -5.0184E-01 -1.6990E-01  3.5659E-01 -2.7341E-01  6.4568E-02  8.4487E-01 -3.8704E+00 -3.4624E-01 -5.2263E-02
             1.4646E-01
 GRADIENT:  -1.1372E-01  4.2464E-02 -5.1845E-02  1.6192E-01  5.3927E-02 -3.9615E-03  6.9512E-03  2.5892E-04 -6.4187E-02 -8.2601E-03
            -7.1911E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1641.40345986431        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      811
 NPARAMETR:  9.5085E-01  5.4753E-01  7.6357E-01  1.2926E+00  6.8836E-01  9.6519E-01  2.1066E+00  1.0307E-02  6.4026E-01  8.5892E-01
             1.0478E+00
 PARAMETER:  4.9599E-02 -5.0233E-01 -1.6976E-01  3.5667E-01 -2.7345E-01  6.4568E-02  8.4506E-01 -4.4750E+00 -3.4588E-01 -5.2082E-02
             1.4665E-01
 GRADIENT:  -7.2143E-03  8.4332E-04  1.7611E-02 -1.3103E-03 -1.6963E-02 -5.0587E-04  9.3714E-07  7.6649E-05  3.5377E-03  6.7641E-03
             1.0543E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1641.40346180509        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      868
 NPARAMETR:  9.5085E-01  5.4753E-01  7.6356E-01  1.2926E+00  6.8836E-01  9.6519E-01  2.1065E+00  1.0000E-02  6.4027E-01  8.5890E-01
             1.0478E+00
 PARAMETER:  4.9601E-02 -5.0233E-01 -1.6977E-01  3.5667E-01 -2.7344E-01  6.4568E-02  8.4503E-01 -4.5324E+00 -3.4586E-01 -5.2104E-02
             1.4665E-01
 GRADIENT:  -4.6648E-03 -3.3459E-03 -2.0027E-04  2.5743E-02  8.5898E-03 -6.4316E-04 -5.2485E-03  0.0000E+00  5.3606E-03  3.3154E-03
             7.5174E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      868
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.4699E-04  2.7172E-02 -5.0599E-04 -2.8553E-02  4.4087E-03
 SE:             2.9823E-02  2.1301E-02  2.3043E-04  2.3095E-02  2.3467E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8269E-01  2.0208E-01  2.8101E-02  2.1633E-01  8.5098E-01

 ETASHRINKSD(%)  8.9916E-02  2.8639E+01  9.9228E+01  2.2630E+01  2.1383E+01
 ETASHRINKVR(%)  1.7975E-01  4.9076E+01  9.9994E+01  4.0139E+01  3.8193E+01
 EBVSHRINKSD(%)  5.0776E-01  2.9577E+01  9.9281E+01  2.1379E+01  1.9132E+01
 EBVSHRINKVR(%)  1.0129E+00  5.0407E+01  9.9995E+01  3.8187E+01  3.4603E+01
 RELATIVEINF(%)  9.8369E+01  8.0823E+00  4.5261E-04  1.1026E+01  5.1371E+00
 EPSSHRINKSD(%)  4.2878E+01
 EPSSHRINKVR(%)  6.7371E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1641.4034618050944     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -906.25263524135619     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1641.403       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.51E-01  5.48E-01  7.64E-01  1.29E+00  6.88E-01  9.65E-01  2.11E+00  1.00E-02  6.40E-01  8.59E-01  1.05E+00
 


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
+        1.31E+03
 
 TH 2
+       -1.60E+01  4.44E+02
 
 TH 3
+        2.32E+01  2.03E+02  1.07E+03
 
 TH 4
+       -8.66E+00  3.78E+02 -3.88E+02  1.04E+03
 
 TH 5
+       -4.51E+00 -4.03E+02 -1.37E+03  3.92E+02  2.09E+03
 
 TH 6
+       -9.72E-01 -3.49E+00  1.70E+00 -3.47E+00 -3.94E+00  2.10E+02
 
 TH 7
+        1.91E+00  3.80E+01 -4.66E+00 -1.49E+01  3.63E+00  3.61E-01  1.72E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.96E+00 -2.33E+01 -4.02E+01 -1.88E+01  4.85E+01  6.00E-02  1.27E+01  0.00E+00  1.98E+02
 
 TH10
+       -2.21E+00 -5.73E+00 -8.88E+01 -3.15E+01 -4.22E+01 -2.58E-01  2.89E+00  0.00E+00  1.91E+01  1.15E+02
 
 TH11
+       -7.57E+00 -8.65E+00 -4.66E+01 -9.75E+00  1.88E+01  2.19E+00  1.48E+00  0.00E+00  1.92E+01  2.83E+01  2.01E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.034
Stop Time:
Sat Sep 25 10:02:05 CDT 2021
