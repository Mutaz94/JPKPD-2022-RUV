Wed Sep 29 13:48:40 CDT 2021
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
$DATA ../../../../data/spa/A3/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   10.1272602380245        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2739E+02  7.4637E+01  4.6702E+01  4.2056E+01  2.6445E+02  2.8521E+01 -8.5874E+01 -1.8759E+01 -1.5968E+02 -1.9677E+02
            -2.8045E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1213.97762160052        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8064E-01  9.3172E-01  9.8739E-01  1.0998E+00  9.7173E-01  8.6320E-01  9.9694E-01  9.3580E-01  1.1744E+00  8.6862E-01
             4.0742E+00
 PARAMETER:  8.0452E-02  2.9279E-02  8.7311E-02  1.9511E-01  7.1321E-02 -4.7109E-02  9.6931E-02  3.3645E-02  2.6079E-01 -4.0851E-02
             1.5047E+00
 GRADIENT:   2.5048E+00 -3.6870E+01 -3.8262E+01 -1.4957E+01  5.2080E+01 -2.2440E+01  6.3266E+00  5.9539E+00  9.4035E+00  1.5239E+01
             4.1844E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1224.24280290617        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.8250E-01  8.6201E-01  6.9995E-01  1.1624E+00  6.9229E-01  9.1305E-01  8.4307E-01  4.0565E-01  1.2820E+00  4.3838E-01
             4.1849E+00
 PARAMETER:  8.2344E-02 -4.8486E-02 -2.5674E-01  2.5052E-01 -2.6775E-01  9.0336E-03 -7.0699E-02 -8.0226E-01  3.4843E-01 -7.2468E-01
             1.5315E+00
 GRADIENT:  -8.4138E+00  1.1418E+01 -1.0408E+01  3.6886E+01  3.3371E+00 -5.1849E+00  4.5002E+00  1.9748E+00  1.6714E+01  5.7955E+00
             7.1060E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1229.17086525196        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.7331E-01  6.2577E-01  6.0507E-01  1.2470E+00  5.4461E-01  9.3278E-01  9.1513E-01  1.5134E-01  1.1172E+00  3.7171E-01
             3.7840E+00
 PARAMETER:  7.2949E-02 -3.6877E-01 -4.0241E-01  3.2075E-01 -5.0768E-01  3.0413E-02  1.1314E-02 -1.7882E+00  2.1080E-01 -8.8965E-01
             1.4308E+00
 GRADIENT:  -7.8033E+00  2.5467E+01  1.5657E+01  3.8176E+01 -2.4156E+01 -9.3343E-01 -4.9350E-01  2.3055E-01 -1.5746E+00  1.7463E+00
            -9.4371E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1231.19668849527        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      386
 NPARAMETR:  9.9087E-01  5.0887E-01  6.3243E-01  1.3138E+00  5.3489E-01  9.4163E-01  9.6142E-01  9.3674E-02  1.0718E+00  3.2357E-01
             3.8497E+00
 PARAMETER:  9.0826E-02 -5.7556E-01 -3.5818E-01  3.7291E-01 -5.2569E-01  3.9862E-02  6.0654E-02 -2.2679E+00  1.6930E-01 -1.0284E+00
             1.4480E+00
 GRADIENT:   9.9704E+00  1.0889E+01  8.3112E+00  1.3297E+01 -1.7058E+01  1.3389E+00 -3.4661E-01  1.0287E-01 -2.0859E+00  1.1349E+00
            -5.1064E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1234.56176343853        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      561
 NPARAMETR:  9.7097E-01  1.7155E-01  4.6549E-01  1.4047E+00  3.7196E-01  9.2524E-01  2.1177E+00  1.0000E-02  1.0392E+00  6.2648E-02
             3.7981E+00
 PARAMETER:  7.0536E-02 -1.6629E+00 -6.6467E-01  4.3985E-01 -8.8896E-01  2.2297E-02  8.5032E-01 -5.3499E+00  1.3849E-01 -2.6702E+00
             1.4345E+00
 GRADIENT:  -1.7780E+01  5.1657E+00  2.9467E+01  2.2546E+01 -4.1921E+01 -3.2148E+00 -1.2596E+00  0.0000E+00 -7.7453E+00 -2.0502E-01
            -6.8553E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1236.77168119187        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      738
 NPARAMETR:  9.7505E-01  1.1878E-01  3.7921E-01  1.3594E+00  3.1983E-01  9.3579E-01  2.9856E+00  1.0000E-02  1.1248E+00  2.5128E-02
             3.7662E+00
 PARAMETER:  7.4738E-02 -2.0305E+00 -8.6967E-01  4.0706E-01 -1.0400E+00  3.3631E-02  1.1938E+00 -6.6828E+00  2.1764E-01 -3.5838E+00
             1.4261E+00
 GRADIENT:   3.3528E-01  2.5466E+00  3.3285E+00  7.8236E+00 -5.4292E+00 -3.3640E-01  6.5466E-01  0.0000E+00  4.4503E-01 -4.6521E-02
            -1.1863E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1237.27903371355        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      913
 NPARAMETR:  9.7167E-01  6.8831E-02  3.4785E-01  1.3358E+00  2.9614E-01  9.3593E-01  4.2194E+00  1.0000E-02  1.1473E+00  1.0955E-02
             3.7866E+00
 PARAMETER:  7.1259E-02 -2.5761E+00 -9.5597E-01  3.8950E-01 -1.1169E+00  3.3790E-02  1.5397E+00 -8.5431E+00  2.3738E-01 -4.4140E+00
             1.4315E+00
 GRADIENT:  -1.2392E+00  2.5611E+00  3.6490E+00 -6.7168E+00 -8.4283E+00 -5.3185E-01  2.8517E+00  0.0000E+00  1.8013E+00 -8.7867E-03
             6.4705E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1237.34336826624        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     1032
 NPARAMETR:  9.6814E-01  6.5587E-02  3.4886E-01  1.3379E+00  2.9691E-01  9.3615E-01  4.2400E+00  1.0000E-02  1.1298E+00  1.3133E-02
             3.7716E+00
 PARAMETER:  6.7626E-02 -2.6244E+00 -9.5309E-01  3.9112E-01 -1.1143E+00  3.4020E-02  1.5446E+00 -8.6187E+00  2.2206E-01 -4.2327E+00
             1.4275E+00
 GRADIENT:   1.4501E+01  1.1763E+00  1.0094E+01  2.9740E+01  2.6362E+01  1.7532E+00 -2.5571E-01  0.0000E+00  1.4862E+00 -1.1301E-02
             1.6770E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1241.26889359630        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1213
 NPARAMETR:  9.7266E-01  6.3252E-02  3.4845E-01  1.3409E+00  2.9740E-01  9.3719E-01  4.2617E+00  1.0000E-02  1.1457E+00  4.6241E-01
             3.6654E+00
 PARAMETER:  7.2275E-02 -2.6606E+00 -9.5425E-01  3.9336E-01 -1.1127E+00  3.5129E-02  1.5497E+00 -8.6187E+00  2.3603E-01 -6.7131E-01
             1.3989E+00
 GRADIENT:  -3.9612E+00 -4.5192E-01 -3.0386E+01 -1.0290E+01  4.6686E+01 -1.8776E+00  7.3326E-02  0.0000E+00  1.6712E+00 -9.7985E-02
             3.5276E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1243.41717548251        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1390
 NPARAMETR:  9.7268E-01  5.7753E-02  3.4663E-01  1.3390E+00  2.8857E-01  9.4791E-01  4.2205E+00  1.0000E-02  1.1573E+00  5.9419E-01
             3.3732E+00
 PARAMETER:  7.2304E-02 -2.7516E+00 -9.5949E-01  3.9193E-01 -1.1428E+00  4.6508E-02  1.5399E+00 -8.6187E+00  2.4609E-01 -4.2055E-01
             1.3159E+00
 GRADIENT:   2.9273E-01  2.1287E-01 -1.7250E-01 -5.2790E+00  3.1083E+00  5.9350E-01 -1.6583E-01  0.0000E+00  1.3240E-01  5.0160E-01
            -8.1915E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1243.55443410424        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1566
 NPARAMETR:  9.7069E-01  3.7557E-02  3.4664E-01  1.3546E+00  2.8695E-01  9.4653E-01  5.3234E+00  1.0000E-02  1.1498E+00  5.8669E-01
             3.3796E+00
 PARAMETER:  7.0255E-02 -3.1819E+00 -9.5947E-01  4.0348E-01 -1.1485E+00  4.5050E-02  1.7721E+00 -8.6187E+00  2.3958E-01 -4.3326E-01
             1.3178E+00
 GRADIENT:  -2.0411E+00  2.2631E-01 -1.3227E+00  3.5235E+00  7.9474E-01  8.0635E-02 -9.7456E-02  0.0000E+00 -7.4459E-01 -1.5032E-02
            -3.5742E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1243.70775698841        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1730
 NPARAMETR:  9.6969E-01  1.0000E-02  3.4746E-01  1.3619E+00  2.8485E-01  9.4468E-01  1.0657E+01  1.0000E-02  1.1464E+00  5.8846E-01
             3.3816E+00
 PARAMETER:  6.9218E-02 -4.5775E+00 -9.5710E-01  4.0888E-01 -1.1558E+00  4.3087E-02  2.4662E+00 -8.6187E+00  2.3661E-01 -4.3024E-01
             1.3183E+00
 GRADIENT:  -1.1547E-01  0.0000E+00  7.3281E-01  6.2408E-01 -1.4616E+00 -1.6701E-01  5.9242E-02  0.0000E+00 -4.8053E-02 -1.5093E-02
            -8.6613E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1730
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.1044E-04 -2.2656E-04  1.1253E-04 -1.3152E-02  5.3974E-04
 SE:             2.8370E-02  1.5103E-03  2.1383E-04  2.6589E-02  1.7712E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9408E-01  8.8076E-01  5.9871E-01  6.2085E-01  9.7569E-01

 ETASHRINKSD(%)  4.9555E+00  9.4940E+01  9.9284E+01  1.0924E+01  4.0661E+01
 ETASHRINKVR(%)  9.6653E+00  9.9744E+01  9.9995E+01  2.0655E+01  6.4789E+01
 EBVSHRINKSD(%)  4.4696E+00  9.5539E+01  9.9274E+01  9.4163E+00  4.0199E+01
 EBVSHRINKVR(%)  8.7394E+00  9.9801E+01  9.9995E+01  1.7946E+01  6.4239E+01
 RELATIVEINF(%)  7.5851E+01  2.2367E-02  2.1461E-04  1.9206E+01  1.0828E+00
 EPSSHRINKSD(%)  2.7957E+01
 EPSSHRINKVR(%)  4.8097E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1243.7077569884061     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -508.55693042466794     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.52
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1243.708       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.70E-01  1.00E-02  3.47E-01  1.36E+00  2.85E-01  9.45E-01  1.07E+01  1.00E-02  1.15E+00  5.88E-01  3.38E+00
 


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
+        1.21E+03
 
 TH 2
+        0.00E+00  1.60E+03
 
 TH 3
+       -6.82E+01  0.00E+00  4.78E+03
 
 TH 4
+       -3.82E+01  0.00E+00 -2.10E+02  3.75E+02
 
 TH 5
+        2.59E+02  0.00E+00 -7.59E+03 -2.95E+02  1.41E+04
 
 TH 6
+       -9.82E-01  0.00E+00  3.36E+01 -1.39E+01 -1.51E+01  1.84E+02
 
 TH 7
+       -5.29E-01  0.00E+00  1.97E+01  3.35E+00 -2.96E+01  7.88E-01  8.24E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.53E+01  0.00E+00  5.58E+01 -1.00E+01  4.65E+01  3.30E+00 -1.08E+00  0.00E+00  9.98E+01
 
 TH10
+       -1.56E+01  0.00E+00 -1.10E+02 -1.81E+00  1.83E+02  1.70E+00  3.68E+00  0.00E+00  3.93E-01  6.71E+01
 
 TH11
+       -1.96E+01  0.00E+00 -8.40E+00 -3.95E+00  4.84E+00  3.05E+00  1.42E+00  0.00E+00  4.84E+00  2.13E+01  5.42E+01
 
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
 #CPUT: Total CPU Time in Seconds,       29.247
Stop Time:
Wed Sep 29 13:49:10 CDT 2021
