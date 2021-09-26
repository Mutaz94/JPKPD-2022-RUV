Sat Sep 25 14:51:34 CDT 2021
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
$DATA ../../../../data/spa/D/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26847.0731335074        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3096E+02  6.2107E+02  2.2853E+01  6.0521E+02  1.5495E+02 -3.8548E+03 -1.3273E+03 -6.7732E+01 -1.8852E+03 -9.1326E+02
            -4.8653E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -433.279023381708        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1600E+00  9.8963E-01  8.0896E-01  1.3514E+00  1.7826E+00  1.5407E+00  8.9567E-01  9.4869E-01  6.5619E-01  8.8443E-01
             1.5708E+01
 PARAMETER:  2.4842E-01  8.9580E-02 -1.1201E-01  4.0111E-01  6.7806E-01  5.3225E-01 -1.0183E-02  4.7328E-02 -3.2130E-01 -2.2811E-02
             2.8542E+00
 GRADIENT:  -4.1526E+01  1.7556E+01 -5.2388E+00  2.9432E+01 -4.3952E+00  5.2505E+00  1.6087E+00  3.5305E+00  6.3115E+00  8.5255E-02
            -4.0195E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -443.190185542478        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.2181E+00  7.6288E-01  8.7014E-01  1.4798E+00  1.8583E+00  1.4137E+00  4.7156E-01  3.1216E-01  3.2259E-01  1.4109E+00
             1.7044E+01
 PARAMETER:  2.9726E-01 -1.7065E-01 -3.9103E-02  4.9191E-01  7.1968E-01  4.4623E-01 -6.5170E-01 -1.0642E+00 -1.0314E+00  4.4421E-01
             2.9358E+00
 GRADIENT:  -2.3978E+01  1.0872E+01 -2.1119E+00  1.6603E+01 -4.3351E+00 -9.1447E+00  2.6547E-01  3.9010E-01  2.2097E+00 -4.4523E-03
             7.7311E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -453.980804888927        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  1.2655E+00  3.9666E-01  8.5264E-01  1.6705E+00  6.9173E+00  1.4438E+00  4.7075E-01  5.8497E-02  8.0982E-02  1.1195E+01
             1.7477E+01
 PARAMETER:  3.3547E-01 -8.2468E-01 -5.9413E-02  6.1314E-01  2.0340E+00  4.6731E-01 -6.5342E-01 -2.7388E+00 -2.4135E+00  2.5155E+00
             2.9609E+00
 GRADIENT:   2.1328E+01  2.4071E+00  5.7234E+00 -7.0987E-01  5.6544E+00 -1.3009E+00  1.0577E-01  4.8485E-03  2.4519E-01  1.0077E+00
             1.7521E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -457.918745920022        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  1.1355E+00  2.2349E-01  6.9852E-01  1.6174E+00  4.8747E+00  1.4575E+00  1.1208E+00  7.8045E-02  6.2291E-02  9.5057E+00
             1.6384E+01
 PARAMETER:  2.2708E-01 -1.3984E+00 -2.5879E-01  5.8081E-01  1.6841E+00  4.7670E-01  2.1407E-01 -2.4505E+00 -2.6759E+00  2.3519E+00
             2.8963E+00
 GRADIENT:  -2.5383E+01  4.3919E+00 -4.6888E+00  1.0309E+01  2.5620E+00  2.1216E+00  1.9865E-01  1.2974E-02  1.8445E-01 -6.6320E+00
             1.3766E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -460.664955502177        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.1544E+00  4.1645E-02  6.6317E-01  1.6980E+00  3.9068E+00  1.4926E+00  6.9066E-02  1.4941E-01  7.0767E-02  8.7598E+00
             1.6501E+01
 PARAMETER:  2.4357E-01 -3.0786E+00 -3.1072E-01  6.2944E-01  1.4627E+00  5.0049E-01 -2.5727E+00 -1.8011E+00 -2.5484E+00  2.2702E+00
             2.9034E+00
 GRADIENT:  -3.1921E+00  6.7715E-01 -1.7256E+00 -3.2216E-01  5.0224E-01  5.9447E+00  2.4935E-05  5.8984E-02  2.6556E-01 -4.9975E+00
             1.0954E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -461.455925286825        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.1293E+00  1.0000E-02  6.6365E-01  1.6634E+00  4.5419E+00  1.4462E+00  1.3185E-02  2.4627E-01  5.1670E-02  9.3780E+00
             1.6234E+01
 PARAMETER:  2.2164E-01 -4.7426E+00 -3.1000E-01  6.0888E-01  1.6133E+00  4.6891E-01 -4.2287E+00 -1.3013E+00 -2.8629E+00  2.3384E+00
             2.8871E+00
 GRADIENT:   6.0897E+00  0.0000E+00 -1.5368E+00  4.2005E+00 -3.7856E+00 -2.7729E+00  5.6001E-08  1.4838E-01  1.3653E-01  1.0563E+01
            -6.8369E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -461.552362747152        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      516
 NPARAMETR:  1.1250E+00  1.0000E-02  6.5846E-01  1.6542E+00  4.6354E+00  1.4406E+00  4.0091E-02  8.3065E-02  3.6780E-02  9.4432E+00
             1.6189E+01
 PARAMETER:  2.1779E-01 -4.5212E+00 -3.1785E-01  6.0335E-01  1.6337E+00  4.6507E-01 -3.1166E+00 -2.3881E+00 -3.2028E+00  2.3453E+00
             2.8843E+00
 GRADIENT:   1.5101E+01  0.0000E+00 -3.2451E+00  6.5072E+00 -8.1740E+00 -6.3513E+00  5.1507E-07  1.6062E-02  6.8199E-02  2.3359E+01
            -1.6322E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -461.590109949175        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      589
 NPARAMETR:  1.1298E+00  1.0000E-02  6.6082E-01  1.6592E+00  4.4840E+00  1.4430E+00  2.0151E-01  1.0404E-02  1.5357E-02  9.3142E+00
             1.6237E+01
 PARAMETER:  2.2201E-01 -4.5958E+00 -3.1427E-01  6.0634E-01  1.6005E+00  4.6675E-01 -1.5019E+00 -4.4656E+00 -4.0762E+00  2.3315E+00
             2.8873E+00
 GRADIENT:  -3.7572E+00  0.0000E+00  5.5874E-01 -1.9525E+00  2.0069E+00  1.6044E+00  1.2467E-05  3.0035E-04  1.2611E-02 -5.5094E+00
             4.3932E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -461.617772265542        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  1.1348E+00  1.0000E-02  6.6830E-01  1.6694E+00  4.6417E+00  1.4470E+00  4.8348E-01  1.0000E-02  1.0000E-02  9.4416E+00
             1.6320E+01
 PARAMETER:  2.2643E-01 -4.6123E+00 -3.0302E-01  6.1246E-01  1.6351E+00  4.6948E-01 -6.2674E-01 -5.6346E+00 -4.5719E+00  2.3451E+00
             2.8924E+00
 GRADIENT:  -3.5306E+00  0.0000E+00  6.2416E-01 -1.1620E+00  1.7356E+00  9.4753E-01  7.0345E-05  0.0000E+00  0.0000E+00 -4.7724E+00
             3.6398E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -461.619374963581        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  1.1343E+00  1.0000E-02  6.6732E-01  1.6684E+00  4.6525E+00  1.4472E+00  4.6957E-01  1.0000E-02  1.0000E-02  9.4547E+00
             1.6314E+01
 PARAMETER:  2.2605E-01 -4.5584E+00 -3.0449E-01  6.1186E-01  1.6374E+00  4.6962E-01 -6.5593E-01 -5.5079E+00 -4.5168E+00  2.3465E+00
             2.8920E+00
 GRADIENT:  -2.1837E-01  0.0000E+00  4.2072E-02 -1.8874E-01  1.4493E-01  2.4246E-01  6.6612E-05  0.0000E+00  0.0000E+00 -4.0453E-01
             2.9551E-01

0ITERATION NO.:   54    OBJECTIVE VALUE:  -461.619400526821        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1034
 NPARAMETR:  1.1343E+00  1.0000E-02  6.6733E-01  1.6683E+00  4.6548E+00  1.4470E+00  4.6385E-01  1.0000E-02  1.0000E-02  9.4564E+00
             1.6312E+01
 PARAMETER:  2.2600E-01 -4.5562E+00 -3.0446E-01  6.1182E-01  1.6379E+00  4.6950E-01 -6.6820E-01 -5.4878E+00 -4.5084E+00  2.3467E+00
             2.8919E+00
 GRADIENT:   2.1424E-02  0.0000E+00  7.7669E-03 -2.1991E-02  8.5146E-03  5.8341E-02  6.5011E-05  0.0000E+00  7.2337E-04 -2.3151E-02
            -1.9294E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1034
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0154E-02 -3.6239E-05 -1.1879E-04 -2.2879E-04 -2.8239E-02
 SE:             2.7863E-02  2.1527E-05  4.7248E-05  2.0252E-04  7.5353E-03
 N:                     100         100         100         100         100

 P VAL.:         7.1554E-01  9.2290E-02  1.1932E-02  2.5860E-01  1.7860E-04

 ETASHRINKSD(%)  6.6549E+00  9.9928E+01  9.9842E+01  9.9322E+01  7.4756E+01
 ETASHRINKVR(%)  1.2867E+01  1.0000E+02  1.0000E+02  9.9995E+01  9.3627E+01
 EBVSHRINKSD(%)  8.2272E+00  9.9908E+01  9.9767E+01  9.9178E+01  5.7526E+01
 EBVSHRINKVR(%)  1.5778E+01  1.0000E+02  9.9999E+01  9.9993E+01  8.1959E+01
 RELATIVEINF(%)  2.7401E+01  2.8743E-06  1.0436E-04  1.7011E-04  9.3134E+00
 EPSSHRINKSD(%)  2.5461E+00
 EPSSHRINKVR(%)  5.0275E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -461.61940052682058     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       273.53142603691759     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.93
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.66
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -461.619       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.13E+00  1.00E-02  6.67E-01  1.67E+00  4.65E+00  1.45E+00  4.64E-01  1.00E-02  1.00E-02  9.46E+00  1.63E+01
 


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
+        2.96E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -7.59E+02  0.00E+00  3.95E+02
 
 TH 4
+        4.18E+02  0.00E+00 -2.59E+02  4.39E+02
 
 TH 5
+       -3.20E+02  0.00E+00  9.81E+01 -7.34E+01  4.10E+01
 
 TH 6
+       -5.94E+02  0.00E+00  1.79E+02 -3.15E+02  9.17E+01  4.23E+02
 
 TH 7
+       -1.89E+01  0.00E+00 -8.16E+00 -1.89E-01 -1.88E-01  1.67E+00  1.14E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.80E+00  0.00E+00 -2.38E+01  3.52E+00 -2.76E+00 -1.97E+01  3.47E+01  0.00E+00  2.30E+03
 
 TH10
+        4.14E+02  0.00E+00 -1.84E+02  5.03E+01 -5.12E+01 -1.27E+02  1.23E-01  0.00E+00  2.89E-01  6.00E+01
 
 TH11
+       -1.52E+02  0.00E+00  5.24E+01 -3.37E+01  1.78E+01  4.27E+01 -6.83E-02  0.00E+00 -1.26E-01 -2.13E+01  9.96E+00
 
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
 #CPUT: Total CPU Time in Seconds,       26.122
Stop Time:
Sat Sep 25 14:52:09 CDT 2021
