Sat Sep 25 07:56:47 CDT 2021
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
$DATA ../../../../data/spa/A1/dat22.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m22.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1129.68316942528        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3499E+02  2.1090E+01 -1.3117E+01  4.7807E+01  1.0896E+02  5.8753E+01 -4.0333E+01 -7.4535E+00 -7.1968E+01 -4.7783E+01
            -8.9836E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1420.86197725736        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.1591E-01  1.0044E+00  1.0779E+00  9.9110E-01  9.4577E-01  7.8086E-01  1.0063E+00  9.5429E-01  1.0095E+00  9.2992E-01
             2.2513E+00
 PARAMETER:  1.2162E-02  1.0442E-01  1.7505E-01  9.1062E-02  4.4240E-02 -1.4736E-01  1.0632E-01  5.3213E-02  1.0942E-01  2.7347E-02
             9.1150E-01
 GRADIENT:  -1.8712E+02  1.4287E+01  1.0047E+01 -6.2538E+00 -2.8033E+01 -3.0044E+01  1.0911E-01  3.4720E+00 -3.5611E+00  3.7985E+00
            -3.2321E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1427.75061000889        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.3401E-01  7.9984E-01  1.1029E+00  1.1469E+00  9.0752E-01  8.0470E-01  8.8674E-01  3.0396E-01  9.5681E-01  9.3314E-01
             2.3086E+00
 PARAMETER:  3.1737E-02 -1.2334E-01  1.9792E-01  2.3709E-01  2.9593E-03 -1.1729E-01 -2.0201E-02 -1.0909E+00  5.5852E-02  3.0800E-02
             9.3664E-01
 GRADIENT:  -1.1597E+02  2.0360E+01 -1.3503E+00  3.4528E+01  4.4446E-01 -1.3404E+01 -4.2670E+00  3.0057E-01 -8.6718E+00  3.0428E+00
            -2.4397E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1432.93465347628        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.7457E-01  6.0639E-01  7.9776E-01  1.2403E+00  6.7827E-01  8.3019E-01  1.3202E+00  1.1453E-01  8.9016E-01  6.9777E-01
             2.4026E+00
 PARAMETER:  7.4240E-02 -4.0023E-01 -1.2595E-01  3.1536E-01 -2.8820E-01 -8.6105E-02  3.7780E-01 -2.0669E+00 -1.6358E-02 -2.5987E-01
             9.7656E-01
 GRADIENT:   1.3274E+01  1.7374E+01  1.0152E+01  3.2711E+01 -2.0569E+01  2.4326E+00 -1.4301E+00  9.2791E-02 -6.8067E-01  1.0980E-02
             9.1532E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1433.81544288547        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.6905E-01  4.5205E-01  6.4602E-01  1.2754E+00  5.5757E-01  8.2296E-01  1.8298E+00  2.7148E-02  8.3868E-01  5.9014E-01
             2.3737E+00
 PARAMETER:  6.8561E-02 -6.9397E-01 -3.3693E-01  3.4326E-01 -4.8417E-01 -9.4850E-02  7.0418E-01 -3.5064E+00 -7.5926E-02 -4.2739E-01
             9.6445E-01
 GRADIENT:  -1.5344E+00  2.3643E+00  5.5487E-01  2.1242E+00 -1.8876E+00 -5.0822E-01  7.3647E-01  7.1626E-03 -6.3503E-02  1.8834E-01
            -2.2687E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1434.70553074503        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  9.6614E-01  2.8986E-01  7.4425E-01  1.3870E+00  5.7157E-01  8.2394E-01  2.1300E+00  1.0000E-02  8.3542E-01  6.9232E-01
             2.3970E+00
 PARAMETER:  6.5556E-02 -1.1383E+00 -1.9538E-01  4.2711E-01 -4.5938E-01 -9.3657E-02  8.5612E-01 -6.1189E+00 -7.9815E-02 -2.6770E-01
             9.7424E-01
 GRADIENT:   2.9284E+00  4.9761E+00  3.6793E-01  8.5417E+00 -7.1457E+00  1.2210E+00  3.9123E+00  0.0000E+00  7.6415E-02  2.8827E+00
             3.3084E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1435.94581515780        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  9.6324E-01  1.1939E-01  7.8898E-01  1.4811E+00  5.6068E-01  8.1901E-01  2.8426E+00  1.0000E-02  8.0553E-01  6.7184E-01
             2.4239E+00
 PARAMETER:  6.2542E-02 -2.0254E+00 -1.3702E-01  4.9281E-01 -4.7860E-01 -9.9656E-02  1.1447E+00 -1.1712E+01 -1.1626E-01 -2.9774E-01
             9.8538E-01
 GRADIENT:   9.3051E+00  7.8849E-01  4.8397E+00  8.2472E+00 -7.9118E+00  2.5423E-01 -3.7017E-01  0.0000E+00  5.0706E-01  3.7748E-01
             3.9066E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1436.01026566996        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  9.5750E-01  6.6463E-02  7.4445E-01  1.4948E+00  5.3003E-01  8.1765E-01  3.8692E+00  1.0000E-02  7.9913E-01  6.5777E-01
             2.3998E+00
 PARAMETER:  5.6568E-02 -2.6111E+00 -1.9511E-01  5.0202E-01 -5.3482E-01 -1.0132E-01  1.4531E+00 -1.5877E+01 -1.2424E-01 -3.1890E-01
             9.7538E-01
 GRADIENT:  -3.8887E+00  2.0401E-01  1.4560E+00  3.6536E+00 -3.4751E+00 -2.7033E-01 -2.6184E-01  0.0000E+00 -6.4616E-01  1.1633E-01
            -1.3355E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1436.25647269874        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  9.5785E-01  2.1203E-02  7.5758E-01  1.5217E+00  5.2864E-01  8.1780E-01  7.3524E+00  1.0000E-02  7.9157E-01  6.5511E-01
             2.4139E+00
 PARAMETER:  5.6932E-02 -3.7536E+00 -1.7762E-01  5.1980E-01 -5.3745E-01 -1.0114E-01  2.0950E+00 -2.4189E+01 -1.3373E-01 -3.2295E-01
             9.8126E-01
 GRADIENT:  -4.7371E+00  1.9096E-01  1.7988E+00 -8.1478E+00 -5.1128E+00 -2.3282E-01  2.7151E-01  0.0000E+00 -4.1396E-01 -1.4481E-01
             1.5002E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1436.47895349191        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      796
 NPARAMETR:  9.5886E-01  1.0256E-02  8.5771E-01  1.5542E+00  5.7541E-01  8.1799E-01  1.0835E+01  1.0000E-02  7.8026E-01  6.9408E-01
             2.4143E+00
 PARAMETER:  5.7989E-02 -4.4799E+00 -5.3489E-02  5.4096E-01 -4.5267E-01 -1.0091E-01  2.4828E+00 -2.9463E+01 -1.4813E-01 -2.6517E-01
             9.8142E-01
 GRADIENT:  -2.4297E-01  1.2303E-01  9.6665E-01  1.0233E+00 -1.5253E+00  1.7505E-02  1.5851E-01  0.0000E+00 -6.3217E-02 -3.8400E-02
            -1.5274E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1436.48215567700        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      972
 NPARAMETR:  9.5892E-01  1.0000E-02  8.5576E-01  1.5536E+00  5.7496E-01  8.1788E-01  1.0741E+01  1.0000E-02  7.8036E-01  6.9315E-01
             2.4152E+00
 PARAMETER:  5.8057E-02 -4.5108E+00 -5.5761E-02  5.4058E-01 -4.5345E-01 -1.0104E-01  2.4741E+00 -2.9690E+01 -1.4800E-01 -2.6652E-01
             9.8179E-01
 GRADIENT:   3.3312E-02  0.0000E+00 -1.4439E-02  1.7639E-01 -6.5116E-05 -1.4701E-02  2.6433E-02  0.0000E+00 -2.5273E-02 -2.8345E-02
             7.5559E-02

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1436.48238759184        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1099
 NPARAMETR:  9.5890E-01  1.0000E-02  8.5592E-01  1.5535E+00  5.7503E-01  8.1792E-01  1.0546E+01  1.0000E-02  7.8045E-01  6.9367E-01
             2.4146E+00
 PARAMETER:  5.8036E-02 -4.5135E+00 -5.5576E-02  5.4053E-01 -4.5333E-01 -1.0099E-01  2.4558E+00 -2.9704E+01 -1.4789E-01 -2.6576E-01
             9.8154E-01
 GRADIENT:  -5.3775E-03  0.0000E+00 -1.5104E-03 -4.0457E-03  2.2093E-03  1.2026E-03 -2.1927E-05  0.0000E+00  2.8489E-04  3.8232E-04
            -4.1906E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1099
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.2457E-04  3.0273E-04  7.7368E-05 -1.0765E-02 -1.5664E-02
 SE:             2.8960E-02  1.8100E-03  1.7654E-04  2.7138E-02  1.8717E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8555E-01  8.6717E-01  6.6122E-01  6.9160E-01  4.0266E-01

 ETASHRINKSD(%)  2.9817E+00  9.3936E+01  9.9409E+01  9.0834E+00  3.7297E+01
 ETASHRINKVR(%)  5.8745E+00  9.9632E+01  9.9997E+01  1.7342E+01  6.0683E+01
 EBVSHRINKSD(%)  2.9875E+00  9.4358E+01  9.9382E+01  9.0018E+00  3.7208E+01
 EBVSHRINKVR(%)  5.8857E+00  9.9682E+01  9.9996E+01  1.7193E+01  6.0572E+01
 RELATIVEINF(%)  8.1099E+01  8.0112E-03  1.7417E-04  3.3281E+00  1.2942E+00
 EPSSHRINKSD(%)  3.0775E+01
 EPSSHRINKVR(%)  5.2079E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1436.4823875918410     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -701.33156102810278     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.64
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1436.482       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.59E-01  1.00E-02  8.56E-01  1.55E+00  5.75E-01  8.18E-01  1.05E+01  1.00E-02  7.80E-01  6.94E-01  2.41E+00
 


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
+        1.70E+03
 
 TH 2
+        0.00E+00  1.47E+03
 
 TH 3
+       -7.40E+00  0.00E+00  5.61E+02
 
 TH 4
+       -6.67E+01  0.00E+00 -4.60E+01  6.26E+02
 
 TH 5
+        4.18E+01  0.00E+00 -1.18E+03 -1.63E+02  2.66E+03
 
 TH 6
+       -6.74E+00  0.00E+00  9.02E+00 -1.22E+01 -1.95E+00  2.64E+02
 
 TH 7
+       -2.56E-02  0.00E+00 -2.75E-02 -1.97E-02  1.25E-02  1.63E-02  3.95E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.18E+01  0.00E+00  2.59E+01 -8.03E+00  7.20E+00  4.40E+00  9.90E-03  0.00E+00  2.23E+02
 
 TH10
+       -1.32E+01  0.00E+00 -6.85E+00 -6.92E+00 -3.81E+01  1.26E+00 -1.67E-02  0.00E+00  4.59E+00  7.23E+01
 
 TH11
+       -1.89E+01  0.00E+00 -3.16E+00 -8.66E+00 -1.03E+01  4.00E+00  2.85E-03  0.00E+00  1.12E+01  2.54E+01  5.32E+01
 
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
 #CPUT: Total CPU Time in Seconds,       18.114
Stop Time:
Sat Sep 25 07:57:07 CDT 2021
