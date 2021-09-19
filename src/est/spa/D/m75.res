Sat Sep 18 15:36:51 CDT 2021
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
$DATA ../../../../data/spa/D/dat75.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m75.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   14875.2884872341        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4147E+01  1.6931E+02 -9.4055E+01 -9.1416E+01  2.8815E+02 -2.2115E+03 -9.0361E+02 -5.2137E+01 -1.7545E+03 -5.7378E+02
            -2.7569E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -584.426298758842        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4772E+00  1.0806E+00  9.2852E-01  1.9141E+00  1.4087E+00  2.2375E+00  1.2770E+00  9.5658E-01  1.4760E+00  1.0380E+00
             1.3856E+01
 PARAMETER:  4.9018E-01  1.7753E-01  2.5837E-02  7.4925E-01  4.4267E-01  9.0537E-01  3.4447E-01  5.5613E-02  4.8936E-01  1.3725E-01
             2.7287E+00
 GRADIENT:   1.7863E+01  2.6955E+01 -1.0593E+01  3.9816E+01 -6.0513E+00  4.3784E+01 -2.9703E+00  6.2963E+00 -1.2487E+01  9.3683E-01
             4.1241E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -600.101807994811        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.4290E+00  8.4359E-01  1.2760E+00  2.2257E+00  2.4771E+00  2.3466E+00  2.6041E+00  2.3598E-01  2.5380E+00  2.4946E+00
             1.1763E+01
 PARAMETER:  4.5700E-01 -7.0085E-02  3.4374E-01  9.0008E-01  1.0071E+00  9.5296E-01  1.0571E+00 -1.3440E+00  1.0314E+00  1.0141E+00
             2.5650E+00
 GRADIENT:   1.2449E+01  1.8785E+01 -7.4285E+00  3.2849E+01 -1.5723E+01  1.9829E+00  7.1730E+00  8.8892E-02  1.5693E+01 -8.6456E-01
             4.0727E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -622.246237992301        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.2985E+00  8.0583E-01  9.9017E-01  1.8527E+00  7.1454E+00  2.2901E+00  8.8882E-01  1.1544E-01  2.5634E+00  1.0376E+01
             1.1229E+01
 PARAMETER:  3.6117E-01 -1.1588E-01  9.0116E-02  7.1666E-01  2.0665E+00  9.2858E-01 -1.7864E-02 -2.0590E+00  1.0414E+00  2.4395E+00
             2.5185E+00
 GRADIENT:   9.2141E-01  7.6362E+00  6.4353E+00  2.0720E+01 -8.5775E+00 -2.3425E+00  2.7014E+00  6.3931E-03  1.4741E+01  9.4795E+00
             4.4646E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -667.216622091029        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0788E+00  2.7652E-01  3.4033E-01  1.4169E+00  2.9495E+01  2.1692E+00  1.8217E-01  2.7954E-02  1.7344E+00  1.2209E+01
             1.0239E+01
 PARAMETER:  1.7588E-01 -1.1855E+00 -9.7785E-01  4.4847E-01  3.4842E+00  8.7435E-01 -1.6028E+00 -3.4772E+00  6.5068E-01  2.6021E+00
             2.4262E+00
 GRADIENT:   3.2588E+00  9.0521E+00  2.6795E+01 -2.0722E+01 -1.8025E+00  6.8768E+00  1.6046E-01 -1.6600E-03  1.5437E+01  9.0272E+00
             2.8568E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -713.473727522946        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  6.8198E-01  2.1779E-02  6.1387E-02  6.2272E-01  2.2465E+02  1.6562E+00  1.0000E-02  1.8851E-02  8.0848E-01  9.0297E+00
             9.6550E+00
 PARAMETER: -2.8275E-01 -3.7268E+00 -2.6906E+00 -3.7366E-01  5.5146E+00  6.0451E-01 -6.3509E+00 -3.8712E+00 -1.1260E-01  2.3005E+00
             2.3675E+00
 GRADIENT:   3.6036E+01 -6.1612E-02 -7.1820E+01  1.3109E+02  3.3632E-02 -3.3032E+01  0.0000E+00 -1.8564E-02 -3.6418E+01 -3.5986E-04
            -6.9530E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -733.440835625308        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  3.7591E-01  1.0000E-02  1.4304E-02  1.9834E-01  3.0393E+03  1.6393E+00  1.0000E-02  1.0000E-02  9.4490E-01  1.1455E+01
             9.5454E+00
 PARAMETER: -8.7840E-01 -5.7417E+00 -4.1472E+00 -1.5178E+00  8.1194E+00  5.9424E-01 -1.3896E+01 -4.8076E+00  4.3320E-02  2.5385E+00
             2.3561E+00
 GRADIENT:   6.9910E+00  0.0000E+00 -3.3932E+01  5.4170E+01  1.9354E-04 -7.3346E-01  0.0000E+00  0.0000E+00 -4.6637E+00  1.3149E-06
            -1.5636E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -734.562386803930        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      619
 NPARAMETR:  3.9225E-01  1.0000E-02  1.5758E-02  2.0741E-01  2.6287E+03  1.6474E+00  1.0000E-02  1.0000E-02  9.4469E-01  1.1606E+01
             9.7560E+00
 PARAMETER: -8.3586E-01 -5.5931E+00 -4.0504E+00 -1.4731E+00  7.9742E+00  5.9921E-01 -1.3683E+01 -4.8906E+00  4.3107E-02  2.5515E+00
             2.3779E+00
 GRADIENT:  -2.1784E-02  0.0000E+00  4.5038E-01 -5.8335E-01  8.9980E-05 -3.6923E-01  0.0000E+00  0.0000E+00  1.4218E-01  7.1433E-07
            -4.8921E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -734.564253464168        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      799
 NPARAMETR:  3.9321E-01  1.0000E-02  1.5844E-02  2.0838E-01  6.4363E+01  1.6471E+00  1.0000E-02  1.0000E-02  9.4432E-01  1.1591E+01
             9.7532E+00
 PARAMETER: -8.3341E-01 -5.5853E+00 -4.0449E+00 -1.4684E+00  4.2645E+00  5.9900E-01 -1.3654E+01 -4.8949E+00  4.2709E-02  2.5502E+00
             2.3776E+00
 GRADIENT:  -3.3313E-01  0.0000E+00  6.0435E-01 -6.3473E-01 -4.8044E-04 -5.7916E-01  0.0000E+00  0.0000E+00  9.2675E-02  1.9892E-03
            -1.8292E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -734.565100147598        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      978
 NPARAMETR:  3.9344E-01  1.0000E-02  1.5847E-02  2.0850E-01  6.0063E+01  1.6500E+00  1.0000E-02  1.0000E-02  9.4380E-01  1.1324E+01
             9.7570E+00
 PARAMETER: -8.3284E-01 -5.5853E+00 -4.0447E+00 -1.4678E+00  4.1954E+00  6.0078E-01 -1.3654E+01 -4.8949E+00  4.2159E-02  2.5269E+00
             2.3780E+00
 GRADIENT:  -7.8077E-02  0.0000E+00  3.7726E-02 -1.2401E-01 -4.7935E-04  9.2289E-03  0.0000E+00  0.0000E+00 -2.6235E-03  2.2556E-03
             2.9510E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -734.566540029027        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1157
 NPARAMETR:  3.9720E-01  1.0000E-02  1.6197E-02  2.1233E-01  2.9815E+01  1.6521E+00  1.0000E-02  1.0000E-02  9.4364E-01  4.1440E+00
             9.7583E+00
 PARAMETER: -8.2332E-01 -5.5853E+00 -4.0229E+00 -1.4496E+00  3.4950E+00  6.0202E-01 -1.3654E+01 -4.8949E+00  4.1995E-02  1.5217E+00
             2.3781E+00
 GRADIENT:  -1.5789E-01  0.0000E+00  2.1040E-01 -1.3304E-01 -1.0215E-03 -4.9201E-02  0.0000E+00  0.0000E+00 -2.0901E-03  1.4013E-03
             2.7038E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -734.566765189587        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1345
 NPARAMETR:  3.9724E-01  1.0000E-02  1.6201E-02  2.1237E-01  2.9536E+01  1.6522E+00  1.0000E-02  1.0000E-02  9.4379E-01  3.3094E+00
             9.7590E+00
 PARAMETER: -8.2323E-01 -5.5853E+00 -4.0227E+00 -1.4494E+00  3.4856E+00  6.0212E-01 -1.3654E+01 -4.8949E+00  4.2144E-02  1.2968E+00
             2.3782E+00
 GRADIENT:  -1.8502E-01  0.0000E+00  2.3826E-01 -1.5286E-01 -4.3782E-04 -1.9682E-02  0.0000E+00  0.0000E+00  2.3644E-02  8.9978E-04
             9.0881E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -734.567488901481        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1522
 NPARAMETR:  3.9663E-01  1.0000E-02  1.6131E-02  2.1167E-01  2.8257E+01  1.6519E+00  1.0000E-02  1.0000E-02  9.4370E-01  3.8948E-01
             9.7584E+00
 PARAMETER: -8.2475E-01 -5.5853E+00 -4.0270E+00 -1.4527E+00  3.4413E+00  6.0195E-01 -1.3654E+01 -4.8949E+00  4.2051E-02 -8.4295E-01
             2.3781E+00
 GRADIENT:  -7.3509E-02  0.0000E+00 -7.7943E-02  1.1493E-01 -4.9917E-04  3.2716E-04  0.0000E+00  0.0000E+00 -7.8393E-03  1.4382E-05
             3.2438E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -734.567518596203        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1700
 NPARAMETR:  3.9658E-01  1.0000E-02  1.6125E-02  2.1158E-01  2.8911E+01  1.6519E+00  1.0000E-02  1.0000E-02  9.4374E-01  3.0807E-01
             9.7580E+00
 PARAMETER: -8.2487E-01 -5.5853E+00 -4.0274E+00 -1.4531E+00  3.4642E+00  6.0192E-01 -1.3654E+01 -4.8949E+00  4.2091E-02 -1.0774E+00
             2.3781E+00
 GRADIENT:   1.3599E-04  0.0000E+00 -4.9487E-02  4.5301E-02 -4.2954E-05  1.1138E-03  0.0000E+00  0.0000E+00 -2.0641E-03  8.4053E-06
            -2.4953E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1700
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.6447E-04  1.8924E-06  1.3669E-04 -2.2827E-02 -5.5843E-06
 SE:             2.8948E-02  2.0609E-06  2.2646E-04  2.3534E-02  1.0407E-05
 N:                     100         100         100         100         100

 P VAL.:         9.7893E-01  3.5849E-01  5.4610E-01  3.3206E-01  5.9155E-01

 ETASHRINKSD(%)  3.0212E+00  9.9993E+01  9.9241E+01  2.1159E+01  9.9965E+01
 ETASHRINKVR(%)  5.9511E+00  1.0000E+02  9.9994E+01  3.7841E+01  1.0000E+02
 EBVSHRINKSD(%)  2.8131E+00  9.9984E+01  9.9286E+01  2.1981E+01  9.9961E+01
 EBVSHRINKVR(%)  5.5471E+00  1.0000E+02  9.9995E+01  3.9131E+01  1.0000E+02
 RELATIVEINF(%)  2.1752E+00  3.7212E-07  2.3720E-05  2.8610E-01  6.2348E-07
 EPSSHRINKSD(%)  1.3495E+01
 EPSSHRINKVR(%)  2.5168E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -734.56751859620272     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      0.58330796753546110     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.24
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.32
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -734.568       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.97E-01  1.00E-02  1.61E-02  2.12E-01  2.89E+01  1.65E+00  1.00E-02  1.00E-02  9.44E-01  3.08E-01  9.76E+00
 


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
+        2.42E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.15E+04  0.00E+00  2.04E+06
 
 TH 4
+       -2.35E+02  0.00E+00 -1.68E+05  1.58E+04
 
 TH 5
+        5.25E-02  0.00E+00 -1.03E+00  6.80E-02  1.57E-05
 
 TH 6
+       -9.86E-01  0.00E+00  1.16E+01 -3.58E+01  6.19E-04  6.00E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.06E+00  0.00E+00  1.92E+03 -1.71E+02 -2.21E-03 -1.59E+00  0.00E+00  0.00E+00  8.56E+01
 
 TH10
+        3.19E-01  0.00E+00  7.55E-01  6.40E-02  8.55E-04 -7.16E-02  0.00E+00  0.00E+00 -6.84E-01  6.53E-01
 
 TH11
+       -2.13E+01  0.00E+00  3.87E+02 -1.88E+01 -6.41E-04  1.29E+00  0.00E+00  0.00E+00  3.26E+00 -3.45E-04  3.70E+00
 
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
 #CPUT: Total CPU Time in Seconds,       34.636
Stop Time:
Sat Sep 18 15:37:27 CDT 2021
