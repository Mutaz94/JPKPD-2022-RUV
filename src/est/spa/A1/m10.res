Wed Sep 29 11:54:41 CDT 2021
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
$DATA ../../../../data/spa/A1/dat10.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m10.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1458.35143757829        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7053E+02  3.3614E+01 -1.7723E+01  6.8906E+01  5.4372E+01  2.3014E+01  2.9455E+00  6.7035E+00 -1.5108E+01  8.8958E+00
            -3.4509E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1523.84886927799        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.6548E-01  9.6189E-01  1.0723E+00  1.0268E+00  1.0008E+00  1.0427E+00  9.1765E-01  9.1413E-01  1.0510E+00  7.8954E-01
             2.0201E+00
 PARAMETER:  6.4875E-02  6.1148E-02  1.6983E-01  1.2642E-01  1.0076E-01  1.4180E-01  1.4061E-02  1.0212E-02  1.4978E-01 -1.3630E-01
             8.0316E-01
 GRADIENT:   5.6551E+01 -1.0189E+01 -1.7676E+01  9.1233E+00  3.0317E+01  1.5734E+01  7.2687E+00  5.1126E+00  8.1175E+00  8.3822E+00
             9.2389E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1527.05224975216        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.5838E-01  9.5425E-01  9.0995E-01  1.0302E+00  8.8161E-01  1.0534E+00  7.6107E-01  7.1090E-01  1.1001E+00  6.5210E-01
             1.9770E+00
 PARAMETER:  5.7490E-02  5.3168E-02  5.6369E-03  1.2979E-01 -2.6010E-02  1.5200E-01 -1.7303E-01 -2.4123E-01  1.9544E-01 -3.2756E-01
             7.8156E-01
 GRADIENT:   4.5996E+01  1.7192E+01  1.6091E+00  2.7820E+01 -1.6935E+01  2.0568E+01  1.5718E+00  2.9408E+00  1.1265E+01  2.8665E+00
             8.1490E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1536.50796271867        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      271
 NPARAMETR:  9.7255E-01  7.2741E-01  7.9824E-01  1.1595E+00  7.3356E-01  1.0445E+00  1.2563E+00  2.5864E-01  9.0247E-01  6.3274E-01
             1.7170E+00
 PARAMETER:  7.2170E-02 -2.1827E-01 -1.2535E-01  2.4797E-01 -2.0984E-01  1.4352E-01  3.2821E-01 -1.2523E+00 -2.6220E-03 -3.5770E-01
             6.4058E-01
 GRADIENT:  -3.1674E+01  1.7777E+01  2.2729E+01  4.4625E+00 -3.2354E+01 -7.5580E+00 -9.4898E-02  9.9587E-03 -4.4738E-01 -2.4708E+00
             2.3654E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1539.91018932218        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  9.9099E-01  4.7784E-01  6.4728E-01  1.2754E+00  5.7427E-01  1.0919E+00  1.7964E+00  7.5305E-02  8.0785E-01  6.2058E-01
             1.5637E+00
 PARAMETER:  9.0948E-02 -6.3848E-01 -3.3498E-01  3.4325E-01 -4.5465E-01  1.8795E-01  6.8580E-01 -2.4862E+00 -1.1338E-01 -3.7710E-01
             5.4707E-01
 GRADIENT:   7.8986E+00  1.0859E+01  1.2026E+01  1.3657E+01 -2.0309E+01  8.4233E+00  6.4325E-01  4.0058E-02 -1.2501E+00 -1.6960E+00
            -4.2520E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1541.23483021551        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      623
 NPARAMETR:  9.8373E-01  2.6246E-01  6.6130E-01  1.3891E+00  5.3359E-01  1.0472E+00  2.4317E+00  1.5956E-02  7.7850E-01  6.8991E-01
             1.5542E+00
 PARAMETER:  8.3592E-02 -1.2376E+00 -3.1355E-01  4.2868E-01 -5.2812E-01  1.4608E-01  9.8859E-01 -4.0379E+00 -1.5039E-01 -2.7120E-01
             5.4096E-01
 GRADIENT:   3.6368E+00  4.2567E+00  5.0472E+00  1.4303E+01 -1.1829E+01 -6.1000E+00 -1.9774E-01  1.6931E-03 -1.2211E+00 -1.0407E+00
            -4.1425E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1541.64062760577        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  9.7648E-01  1.7138E-01  7.5041E-01  1.4566E+00  5.6792E-01  1.0645E+00  3.0029E+00  1.0000E-02  7.6997E-01  7.4296E-01
             1.5766E+00
 PARAMETER:  7.6196E-02 -1.6639E+00 -1.8714E-01  4.7610E-01 -4.6577E-01  1.6252E-01  1.1996E+00 -5.0121E+00 -1.6140E-01 -1.9712E-01
             5.5528E-01
 GRADIENT:  -3.5599E+00  9.4341E-01  5.8017E+00  1.6972E+01 -5.4879E+00  1.4498E+00 -2.9660E+00  0.0000E+00  2.1109E+00  8.9157E-01
             4.0778E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1543.41754629847        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      980
 NPARAMETR:  9.7209E-01  4.7839E-02  6.5763E-01  1.4900E+00  4.9607E-01  1.0798E+00  5.5834E+00  1.0000E-02  7.5130E-01  6.9677E-01
             1.5778E+00
 PARAMETER:  7.1693E-02 -2.9399E+00 -3.1911E-01  4.9878E-01 -6.0103E-01  1.7677E-01  1.8198E+00 -8.5060E+00 -1.8595E-01 -2.6130E-01
             5.5605E-01
 GRADIENT:  -5.5040E+00  2.4633E+00  6.5005E+00  1.4684E+01 -9.6044E+00  7.2688E+00  3.1884E+00  0.0000E+00  4.6068E+00 -4.6899E+00
             1.7876E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1544.43084841607        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1169
 NPARAMETR:  9.7526E-01  2.7059E-02  6.0984E-01  1.4699E+00  4.7033E-01  1.0572E+00  7.0162E+00  1.0000E-02  7.3301E-01  6.8644E-01
             1.5454E+00
 PARAMETER:  7.4947E-02 -3.5097E+00 -3.9456E-01  4.8523E-01 -6.5433E-01  1.5562E-01  2.0482E+00 -1.0167E+01 -2.1059E-01 -2.7624E-01
             5.3531E-01
 GRADIENT:   1.7261E+00 -2.8376E+00 -8.5899E-01 -1.0774E+00 -1.2105E-01 -7.5508E-01 -5.6907E+00  0.0000E+00 -9.1517E-01  1.2914E+00
            -1.0259E+00

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1544.46646413204        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  9.7428E-01  2.7022E-02  6.1062E-01  1.4727E+00  4.6963E-01  1.0589E+00  7.0139E+00  1.0000E-02  7.3552E-01  6.8485E-01
             1.5567E+00
 PARAMETER:  7.2947E-02 -3.5048E+00 -3.9411E-01  4.8613E-01 -6.5454E-01  1.5615E-01  2.0511E+00 -1.0167E+01 -2.0704E-01 -2.7802E-01
             5.4142E-01
 GRADIENT:  -2.3237E+00  2.0266E+01 -2.1572E+02 -1.6520E+02  1.2655E+02 -4.9111E-01  3.0163E+01  0.0000E+00  1.1171E-01  6.5500E-01
            -1.5481E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1233
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5557E-03  1.9263E-02 -9.4637E-05 -1.1263E-02 -7.7867E-03
 SE:             2.9725E-02  9.8337E-03  2.5326E-04  2.8040E-02  2.3472E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5826E-01  5.0127E-02  7.0865E-01  6.8794E-01  7.4008E-01

 ETASHRINKSD(%)  4.1643E-01  6.7056E+01  9.9152E+01  6.0608E+00  2.1366E+01
 ETASHRINKVR(%)  8.3113E-01  8.9147E+01  9.9993E+01  1.1754E+01  3.8167E+01
 EBVSHRINKSD(%)  8.6378E-01  7.9194E+01  9.9122E+01  5.7991E+00  1.9734E+01
 EBVSHRINKVR(%)  1.7201E+00  9.5671E+01  9.9992E+01  1.1262E+01  3.5574E+01
 RELATIVEINF(%)  9.8000E+01  2.9479E+00  3.4792E-04  4.0349E+01  2.9539E+00
 EPSSHRINKSD(%)  3.9085E+01
 EPSSHRINKVR(%)  6.2894E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1544.4664641320389     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -809.31563756830076     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.85
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
 





 #OBJV:********************************************    -1544.466       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.73E-01  2.72E-02  6.10E-01  1.47E+00  4.70E-01  1.06E+00  7.04E+00  1.00E-02  7.36E-01  6.85E-01  1.55E+00
 


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
+        1.03E+03
 
 TH 2
+       -1.07E+02  4.19E+05
 
 TH 3
+        3.63E+01  8.34E+02  3.88E+04
 
 TH 4
+       -8.70E+00  9.14E+02 -3.45E+02  4.89E+03
 
 TH 5
+        5.09E+00 -2.43E+03 -3.57E+03 -9.48E+03  2.87E+04
 
 TH 6
+        2.96E+00 -8.26E+00  1.70E+00 -5.29E+00  2.26E+00  1.72E+02
 
 TH 7
+       -2.79E-01  1.23E+03  2.90E+00 -1.63E+02  3.73E+02 -4.28E-02  2.80E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.86E+00 -1.91E+02  5.78E+04 -2.52E+01  2.76E+01 -4.57E-01 -2.04E+00  0.00E+00  3.14E+02
 
 TH10
+       -4.58E+01  3.69E+02  4.61E+04 -8.51E+01  1.27E+02  2.16E+00  2.01E+00  0.00E+00  1.37E+02  5.85E+04
 
 TH11
+        1.25E+00  1.82E+02  1.04E+04 -4.42E+01  9.05E+01  2.63E+00  1.21E+00  0.00E+00  1.65E+04  1.32E+04  3.07E+03
 
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
 #CPUT: Total CPU Time in Seconds,       22.614
Stop Time:
Wed Sep 29 11:55:06 CDT 2021
