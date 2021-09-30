Wed Sep 29 15:26:01 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat88.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1640.60826789921        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1602E+02 -4.9440E+01 -4.5331E+01  1.2123E+01  8.2483E+01  4.7296E+01 -8.0302E+00 -3.9851E+00  6.2794E+00  2.3454E+01
            -6.4615E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1650.98628831476        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  1.0078E+00  1.0855E+00  1.0734E+00  9.9547E-01  1.0161E+00  9.7766E-01  1.0504E+00  1.0445E+00  9.9095E-01  8.2841E-01
             1.3104E+00
 PARAMETER:  1.0778E-01  1.8207E-01  1.7085E-01  9.5462E-02  1.1596E-01  7.7405E-02  1.4918E-01  1.4354E-01  9.0909E-02 -8.8244E-02
             3.7030E-01
 GRADIENT:   2.6552E+00 -2.2674E+01 -1.2828E+01 -2.3837E+00  3.1307E+01 -2.2010E+00 -4.8651E+00 -7.8733E+00  2.9166E+00  1.1117E+01
             4.8157E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1654.48354582249        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.9774E-01  1.1659E+00  9.3429E-01  9.4372E-01  9.1397E-01  1.0026E+00  1.1857E+00  1.7608E+00  9.6018E-01  3.9875E-01
             1.2124E+00
 PARAMETER:  9.7742E-02  2.5349E-01  3.2031E-02  4.2074E-02  1.0041E-02  1.0264E-01  2.7035E-01  6.6578E-01  5.9368E-02 -8.1943E-01
             2.9260E-01
 GRADIENT:  -1.9352E+01  1.4683E+01  9.8018E+00  2.5150E-01 -5.4165E+01  6.3892E+00  7.7526E+00  1.4695E+01  6.9932E+00  1.9334E+00
             1.9785E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1657.79544773083        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  1.0093E+00  1.2028E+00  9.3311E-01  9.1919E-01  9.6992E-01  9.8687E-01  1.1099E+00  1.3489E+00  9.7359E-01  5.7040E-01
             1.1577E+00
 PARAMETER:  1.0923E-01  2.8469E-01  3.0770E-02  1.5741E-02  6.9463E-02  8.6782E-02  2.0425E-01  3.9928E-01  7.3236E-02 -4.6142E-01
             2.4641E-01
 GRADIENT:   6.4524E+00 -2.5497E-01 -1.9149E+00  3.5693E+00  1.0153E+01  2.5217E-01  8.2351E-01 -7.8073E-01 -7.9077E-01 -1.5264E+00
            -4.0014E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1657.90842495479        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  1.0090E+00  1.3497E+00  8.0576E-01  8.2331E-01  9.7674E-01  9.8819E-01  1.0063E+00  1.3513E+00  1.0590E+00  5.7003E-01
             1.1616E+00
 PARAMETER:  1.0900E-01  3.9990E-01 -1.1597E-01 -9.4423E-02  7.6463E-02  8.8119E-02  1.0626E-01  4.0109E-01  1.5730E-01 -4.6206E-01
             2.4979E-01
 GRADIENT:   3.9314E+00  5.4624E+00  1.5179E+00  4.6699E+00 -3.5394E+00  4.2335E-01 -8.4790E-01  9.3492E-02  7.8122E-01 -2.7424E-01
            -1.4298E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1658.03314055402        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  1.0086E+00  1.5090E+00  5.7698E-01  7.0614E-01  9.4040E-01  9.8805E-01  9.5869E-01  1.1915E+00  1.1156E+00  4.8914E-01
             1.1650E+00
 PARAMETER:  1.0854E-01  5.1142E-01 -4.4994E-01 -2.4795E-01  3.8551E-02  8.7979E-02  5.7812E-02  2.7525E-01  2.0939E-01 -6.1510E-01
             2.5270E-01
 GRADIENT:   8.1241E-01  1.2334E+01  6.0187E+00 -1.4310E+00 -1.3978E+01 -5.7300E-02  1.4046E+00 -2.0311E-01 -7.2962E-01 -6.1877E-02
             3.2401E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1660.26142754587        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  9.9591E-01  1.7624E+00  1.7109E-01  4.6516E-01  8.1232E-01  9.8247E-01  8.4364E-01  4.9339E-01  1.2608E+00  3.3054E-01
             1.1603E+00
 PARAMETER:  9.5897E-02  6.6665E-01 -1.6656E+00 -6.6538E-01 -1.0786E-01  8.2310E-02 -7.0025E-02 -6.0647E-01  3.3176E-01 -1.0070E+00
             2.4870E-01
 GRADIENT:  -8.4856E+00  3.8356E+00  1.1433E+01 -2.4618E+01 -4.3612E+01 -1.0373E+00  8.8328E+00 -3.7374E-01 -4.2555E+00  4.1984E+00
             1.9956E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1662.10545372705        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  9.9765E-01  2.0157E+00  1.1031E-01  3.3086E-01  9.4034E-01  9.8788E-01  7.4906E-01  3.6634E-01  1.7047E+00  3.9863E-01
             1.0894E+00
 PARAMETER:  9.7651E-02  8.0098E-01 -2.1044E+00 -1.0061E+00  3.8484E-02  8.7807E-02 -1.8894E-01 -9.0420E-01  6.3337E-01 -8.1972E-01
             1.8561E-01
 GRADIENT:  -1.3280E+00  7.3776E+00 -4.3361E+00  1.4676E+01  1.4341E+01 -1.8210E-01 -3.7816E+00 -1.1216E+00 -7.6206E-01 -1.7345E+00
            -5.6361E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1662.35273331052        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1417
 NPARAMETR:  1.0013E+00  2.0610E+00  1.0618E-01  2.9897E-01  9.7161E-01  9.8951E-01  7.4664E-01  4.1902E-01  1.8418E+00  4.0011E-01
             1.1012E+00
 PARAMETER:  1.0130E-01  8.2317E-01 -2.1426E+00 -1.1074E+00  7.1196E-02  8.9452E-02 -1.9217E-01 -7.6983E-01  7.1077E-01 -8.1602E-01
             1.9638E-01
 GRADIENT:   5.3660E+00  3.8995E+00  1.1072E+00  1.9625E+00  1.4042E+00  2.3424E-01 -1.4567E+00 -1.2674E+00  6.9009E-01 -2.7967E+00
            -1.9796E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1662.72830138027        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1603
 NPARAMETR:  9.9877E-01  2.0567E+00  1.0678E-01  2.9720E-01  9.7892E-01  9.8837E-01  7.4071E-01  4.7815E-01  1.7773E+00  4.8342E-01
             1.0947E+00
 PARAMETER:  9.8772E-02  8.2109E-01 -2.1370E+00 -1.1133E+00  7.8690E-02  8.8301E-02 -2.0015E-01 -6.3783E-01  6.7510E-01 -6.2686E-01
             1.9045E-01
 GRADIENT:  -1.3048E+00 -1.6400E+01  1.5537E+00 -2.3942E+00  1.8108E+00 -1.6094E-01 -1.1496E+00 -1.4411E+00 -3.9358E+00  4.3339E-01
             4.8757E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1662.73544994227        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1669
 NPARAMETR:  9.9884E-01  2.0542E+00  1.0712E-01  2.9669E-01  9.7868E-01  9.8858E-01  7.4219E-01  4.8076E-01  1.7750E+00  4.8422E-01
             1.0948E+00
 PARAMETER:  9.8687E-02  8.2109E-01 -2.1370E+00 -1.1134E+00  7.8781E-02  8.8253E-02 -2.0015E-01 -6.3145E-01  6.7480E-01 -6.2614E-01
             1.9031E-01
 GRADIENT:  -3.7398E+04  9.1090E+03 -3.4585E+03  6.6929E+03  1.3441E+00 -1.6969E-01 -1.0220E+00  1.1812E+04  1.0977E+04 -1.1948E+04
            -3.9316E+04
 NUMSIGDIG:         2.3         2.3         2.3         2.3         1.9         2.1         1.5         2.3         2.3         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1669
 NO. OF SIG. DIGITS IN FINAL EST.:  1.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0994E-03 -1.9468E-02 -2.8779E-03  3.0961E-02 -3.6609E-02
 SE:             2.9836E-02  2.8108E-02  4.4422E-03  2.3752E-02  1.4898E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7061E-01  4.8856E-01  5.1709E-01  1.9240E-01  1.3995E-02

 ETASHRINKSD(%)  4.5939E-02  5.8350E+00  8.5118E+01  2.0428E+01  5.0091E+01
 ETASHRINKVR(%)  9.1857E-02  1.1330E+01  9.7785E+01  3.6682E+01  7.5091E+01
 EBVSHRINKSD(%)  5.1426E-01  7.2148E+00  8.7662E+01  2.0180E+01  5.0550E+01
 EBVSHRINKVR(%)  1.0259E+00  1.3909E+01  9.8478E+01  3.6288E+01  7.5547E+01
 RELATIVEINF(%)  9.8304E+01  2.4046E+01  7.0032E-01  1.1736E+01  5.4906E+00
 EPSSHRINKSD(%)  4.3976E+01
 EPSSHRINKVR(%)  6.8614E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1662.7354499422702     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -927.58462337853200     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1662.735       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  2.06E+00  1.07E-01  2.97E-01  9.79E-01  9.88E-01  7.41E-01  4.81E-01  1.78E+00  4.84E-01  1.09E+00
 


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
+        1.87E+07
 
 TH 2
+        1.19E+02  6.60E+04
 
 TH 3
+       -1.05E+03  1.15E+04  3.51E+06
 
 TH 4
+        6.61E+02 -2.04E+03  1.41E+04  1.70E+06
 
 TH 5
+       -3.46E+03 -1.20E+02 -2.83E+03  2.06E+03  1.95E+07
 
 TH 6
+       -7.36E+02  4.20E+01 -3.09E+02  2.15E+02 -2.83E+00  2.00E+02
 
 TH 7
+       -3.60E+02  2.58E+01 -2.32E+02  1.10E+02 -1.29E+07 -6.33E-01  2.76E+02
 
 TH 8
+       -6.15E+06 -1.92E+03  1.40E+04 -9.78E+03  1.12E+03  2.41E+02  1.20E+02  2.01E+06
 
 TH 9
+        1.87E+02 -1.73E+03  1.26E+04  4.65E+05  2.59E+02  6.00E+01  3.66E+01 -2.66E+03  1.28E+05
 
 TH10
+        2.82E+03 -1.80E+02  1.32E+03 -8.88E+02  6.31E+06 -2.42E+02 -4.16E+06 -9.26E+02 -2.23E+02  2.04E+06
 
 TH11
+        6.29E+03 -3.88E+02  2.64E+03 -1.89E+03 -1.71E+03 -3.48E+02 -1.63E+02 -2.05E+03 -5.07E+02  1.39E+03  4.31E+06
 
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
 #CPUT: Total CPU Time in Seconds,       28.106
Stop Time:
Wed Sep 29 15:26:31 CDT 2021
