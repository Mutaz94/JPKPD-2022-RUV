Wed Sep 29 04:50:22 CDT 2021
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
$DATA ../../../../data/int/SL3/dat87.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      986
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      886
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1203.44108986238        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1645E+02  1.7232E+01  3.0701E+02  1.3153E+02  2.8546E+02  3.9598E+01 -1.4370E+02 -7.5877E+02 -2.4135E+02 -5.0604E+01
            -8.8842E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2349.93403040945        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1124E+00  1.4009E+00  9.4502E-01  9.5309E-01  1.0233E+00  9.4863E-01  1.0640E+00  1.0119E+00  9.8215E-01  8.4414E-01
             5.2792E+00
 PARAMETER:  2.0655E-01  4.3710E-01  4.3447E-02  5.1950E-02  1.2299E-01  4.7268E-02  1.6200E-01  1.1185E-01  8.1989E-02 -6.9439E-02
             1.7638E+00
 GRADIENT:   9.0494E+01  8.9822E+01 -6.1918E-01  4.8744E+01 -4.8927E+01 -5.1563E+00  1.8433E+01  6.4391E+00  1.0761E+01  9.2919E+00
             7.5904E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2430.83447780529        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0270E+00  6.8564E-01  8.5414E-01  1.3825E+00  6.5299E-01  9.7988E-01  2.1556E+00  5.7399E-01  8.0418E-01  1.1296E-01
             4.4830E+00
 PARAMETER:  1.2666E-01 -2.7741E-01 -5.7658E-02  4.2391E-01 -3.2619E-01  7.9674E-02  8.6806E-01 -4.5514E-01 -1.1793E-01 -2.0807E+00
             1.6003E+00
 GRADIENT:  -5.5618E+01  3.5585E+01  1.5563E+01  2.8061E+02 -2.0521E+01 -1.4575E+00  4.8084E+01  4.1291E+00 -9.0863E+00  1.2501E-01
             5.6022E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2582.25008351817        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0129E+00  6.9517E-01  5.3020E-01  1.1235E+00  5.5346E-01  9.9637E-01  1.5381E+00  1.5190E-01  9.0038E-01  2.2761E-01
             3.0713E+00
 PARAMETER:  1.1282E-01 -2.6360E-01 -5.3450E-01  2.1643E-01 -4.9156E-01  9.6367E-02  5.3056E-01 -1.7845E+00 -4.9420E-03 -1.3801E+00
             1.2221E+00
 GRADIENT:   5.9745E+00  1.4062E+01  3.9216E+01 -2.5040E+00 -3.1236E+01  5.9026E+00 -1.1596E+01  2.4572E-01 -5.4174E+00 -4.5083E-01
            -2.0246E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2587.36348136524        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      403
 NPARAMETR:  1.0400E+00  8.9843E-01  6.7315E-01  1.0746E+00  7.1605E-01  9.9600E-01  1.4474E+00  7.6442E-02  9.8023E-01  4.1458E-01
             3.1246E+00
 PARAMETER:  1.3925E-01 -7.1027E-03 -2.9579E-01  1.7199E-01 -2.3401E-01  9.5995E-02  4.6977E-01 -2.4712E+00  8.0035E-02 -7.8050E-01
             1.2393E+00
 GRADIENT:   5.8277E+00  6.1254E+00 -3.3376E+00  6.5000E+00 -1.4708E+00  9.3356E-01  1.3148E+00  1.1093E-01  1.4202E+00 -1.6825E-01
             2.9645E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2591.28774971695        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      583
 NPARAMETR:  1.0380E+00  1.3055E+00  1.0492E+00  9.3862E-01  1.0923E+00  9.8852E-01  1.0009E+00  1.4373E-02  1.0275E+00  1.0257E+00
             3.1094E+00
 PARAMETER:  1.3731E-01  3.6657E-01  1.4801E-01  3.6656E-02  1.8830E-01  8.8452E-02  1.0093E-01 -4.1424E+00  1.2713E-01  1.2540E-01
             1.2344E+00
 GRADIENT:  -1.7907E+00  3.3314E+01 -6.0390E+00  7.9955E+01  7.7589E+00 -2.2380E+00  6.8530E+00  1.4924E-03 -6.6953E+00  4.2360E+00
            -3.0818E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2605.01005533505        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      759
 NPARAMETR:  1.0331E+00  1.8044E+00  1.5212E+00  5.8389E-01  1.5504E+00  9.9504E-01  6.6343E-01  1.0000E-02  1.5396E+00  1.4383E+00
             3.0499E+00
 PARAMETER:  1.3259E-01  6.9023E-01  5.1947E-01 -4.3805E-01  5.3849E-01  9.5029E-02 -3.1034E-01 -5.0347E+00  5.3150E-01  4.6343E-01
             1.2151E+00
 GRADIENT:  -7.2621E+00 -5.4952E+00 -7.0722E-01  1.8626E+01 -2.0249E+00  4.0138E-01 -7.8304E+00  0.0000E+00  3.4192E+00  5.2842E+00
            -5.1073E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2610.14170660058        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      934
 NPARAMETR:  1.0351E+00  2.2094E+00  1.8302E+00  3.1204E-01  1.8361E+00  9.9522E-01  6.6836E-01  1.0000E-02  2.0103E+00  1.6138E+00
             3.0334E+00
 PARAMETER:  1.3452E-01  8.9271E-01  7.0444E-01 -1.0646E+00  7.0763E-01  9.5209E-02 -3.0293E-01 -5.2705E+00  7.9830E-01  5.7859E-01
             1.2097E+00
 GRADIENT:  -2.9765E+00 -4.7073E-01  2.3130E+00  2.2902E+00  1.7777E+00  1.2183E-01 -3.0889E+00  0.0000E+00 -2.2469E+00 -2.0428E+00
            -4.1729E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2611.60006500258        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1112
 NPARAMETR:  1.0380E+00  2.5532E+00  1.6796E+00  9.8696E-02  1.9097E+00  9.9562E-01  6.6186E-01  1.0000E-02  3.7828E+00  1.7250E+00
             3.0324E+00
 PARAMETER:  1.3732E-01  1.0374E+00  6.1855E-01 -2.2157E+00  7.4696E-01  9.5606E-02 -3.1270E-01 -5.3742E+00  1.4305E+00  6.4521E-01
             1.2093E+00
 GRADIENT:   1.9894E+00  4.2394E+01  4.6650E+00  5.0524E+00 -1.0361E+01 -1.8370E-01 -4.1827E+00  0.0000E+00  4.0365E+00  4.3810E+00
            -9.2803E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2613.25123775272        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1288
 NPARAMETR:  1.0358E+00  2.5911E+00  1.4169E+00  4.6368E-02  1.9283E+00  9.9567E-01  6.5201E-01  1.0000E-02  4.7082E+00  1.6789E+00
             3.0247E+00
 PARAMETER:  1.3518E-01  1.0521E+00  4.4845E-01 -2.9711E+00  7.5665E-01  9.5658E-02 -3.2770E-01 -5.7560E+00  1.6493E+00  6.1814E-01
             1.2068E+00
 GRADIENT:  -1.4924E+00  6.7013E-02  4.5148E+00 -2.9417E+00 -2.5217E+00 -2.4778E-03  1.3270E+00  0.0000E+00 -9.2159E+00 -2.6128E+00
            -2.3720E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2615.42620244181        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1474
 NPARAMETR:  1.0422E+00  2.6004E+00  6.9814E-01  3.9085E-02  1.9562E+00  1.0000E+00  6.5657E-01  1.0000E-02  5.0856E+00  1.7101E+00
             3.0154E+00
 PARAMETER:  1.4133E-01  1.0556E+00 -2.5933E-01 -3.1420E+00  7.7099E-01  1.0004E-01 -3.2072E-01 -5.8079E+00  1.7264E+00  6.3653E-01
             1.2037E+00
 GRADIENT:   1.2319E+01 -8.1814E+00  1.1703E+00 -2.6763E+00  6.5857E+00  1.2822E+00  4.4793E+00  0.0000E+00 -8.6138E+00  1.9349E+00
            -5.3818E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2615.93171051358        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1653
 NPARAMETR:  1.0302E+00  2.6092E+00  5.7917E-01  4.0521E-02  1.9325E+00  9.9380E-01  6.5629E-01  1.0000E-02  5.2813E+00  1.6777E+00
             3.0333E+00
 PARAMETER:  1.2974E-01  1.0590E+00 -4.4615E-01 -3.1059E+00  7.5884E-01  9.3777E-02 -3.2116E-01 -5.8079E+00  1.7642E+00  6.1743E-01
             1.2097E+00
 GRADIENT:  -1.3386E+01  8.4417E+00  1.9964E-01 -9.9498E-01 -3.8651E-01 -1.1789E+00  3.5574E+00  0.0000E+00 -5.6622E+00 -1.8153E+00
             7.2989E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2616.08340785194        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1828
 NPARAMETR:  1.0364E+00  2.6049E+00  5.3570E-01  4.0858E-02  1.9354E+00  9.9647E-01  6.5617E-01  1.0000E-02  5.3831E+00  1.6913E+00
             3.0233E+00
 PARAMETER:  1.3573E-01  1.0574E+00 -5.2418E-01 -3.0977E+00  7.6034E-01  9.6465E-02 -3.2133E-01 -5.8079E+00  1.7833E+00  6.2552E-01
             1.2064E+00
 GRADIENT:   2.5818E-01 -6.7124E+00  6.5461E-02  1.3407E+00 -3.8853E-01 -9.3479E-03 -7.0668E-01  0.0000E+00  5.0906E-01 -1.4305E-01
            -1.0420E+00

0ITERATION NO.:   61    OBJECTIVE VALUE:  -2616.08340785194        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:     1879
 NPARAMETR:  1.0364E+00  2.6083E+00  5.3290E-01  4.0376E-02  1.9357E+00  9.9653E-01  6.5663E-01  1.0000E-02  5.3633E+00  1.6920E+00
             3.0240E+00
 PARAMETER:  1.3573E-01  1.0574E+00 -5.2418E-01 -3.0977E+00  7.6034E-01  9.6465E-02 -3.2133E-01 -5.8079E+00  1.7833E+00  6.2552E-01
             1.2064E+00
 GRADIENT:  -1.4706E-01 -2.8735E+01  5.2804E-02  1.0646E+00 -1.8335E-01 -4.5620E-02 -7.4908E-01  0.0000E+00  1.6084E+01 -1.8482E-01
            -9.7897E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1879
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2363E-03 -1.5169E-02  1.0236E-06  1.2430E-02 -1.6907E-02
 SE:             2.9178E-02  2.6741E-02  2.7437E-06  8.9289E-03  2.6554E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6620E-01  5.7054E-01  7.0909E-01  1.6390E-01  5.2432E-01

 ETASHRINKSD(%)  2.2488E+00  1.0415E+01  9.9991E+01  7.0087E+01  1.1040E+01
 ETASHRINKVR(%)  4.4470E+00  1.9745E+01  1.0000E+02  9.1052E+01  2.0861E+01
 EBVSHRINKSD(%)  2.0795E+00  8.6623E+00  9.9952E+01  7.8976E+01  9.7714E+00
 EBVSHRINKVR(%)  4.1158E+00  1.6574E+01  1.0000E+02  9.5580E+01  1.8588E+01
 RELATIVEINF(%)  9.5773E+01  2.1700E+01  2.2316E-05  1.1620E+00  7.7061E+01
 EPSSHRINKSD(%)  1.5187E+01
 EPSSHRINKVR(%)  2.8067E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          886
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1628.3590808386800     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2616.0834078519420     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -987.72432701326193     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    49.64
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2616.083       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.60E+00  5.36E-01  4.09E-02  1.94E+00  9.96E-01  6.56E-01  1.00E-02  5.38E+00  1.69E+00  3.02E+00
 


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
+        1.00E+03
 
 TH 2
+       -2.80E+01  4.94E+02
 
 TH 3
+        1.78E-01 -2.52E+01  1.72E+01
 
 TH 4
+        3.17E+02 -2.59E+03 -3.11E+02  5.23E+04
 
 TH 5
+       -3.91E+00  1.36E+01 -3.39E+00 -3.02E+03  8.55E+01
 
 TH 6
+        7.30E+00 -1.93E+01  1.94E-01  4.09E+02 -1.82E+00  1.83E+02
 
 TH 7
+        5.55E-01  1.41E+03 -4.15E+01 -2.13E+04  2.09E+01 -3.43E+00  6.13E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.56E+00 -4.99E+01  8.09E+00  8.61E+02 -8.70E+00  4.07E+00 -5.54E+01  0.00E+00  1.45E+01
 
 TH10
+       -8.65E-01  1.77E+01  1.60E+00 -3.86E+03 -5.55E+00 -5.06E-01 -2.53E+00  0.00E+00 -5.54E+00  4.32E+01
 
 TH11
+       -1.69E+01  6.74E+01 -1.52E+00 -1.25E+03  2.22E+00  2.29E+00  2.18E+01  0.00E+00 -2.79E+00  6.54E+00  1.29E+02
 
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
 #CPUT: Total CPU Time in Seconds,       64.546
Stop Time:
Wed Sep 29 04:51:28 CDT 2021
