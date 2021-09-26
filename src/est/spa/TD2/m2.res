Sat Sep 25 13:15:46 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat2.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1690.91597628148        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.5286E+01 -4.9214E+01 -1.9966E+01 -5.2971E+01 -3.7515E+01  2.6373E+00  8.7939E+00  1.7407E+01  3.5410E+01  1.7061E+01
             5.4612E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1697.89076370676        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.7682E-01  1.1478E+00  1.2668E+00  9.6014E-01  1.2484E+00  9.8366E-01  9.2794E-01  8.4855E-01  7.3469E-01  9.4486E-01
             1.0165E+00
 PARAMETER:  7.6547E-02  2.3786E-01  3.3647E-01  5.9320E-02  3.2190E-01  8.3526E-02  2.5210E-02 -6.4228E-02 -2.0831E-01  4.3278E-02
             1.1641E-01
 GRADIENT:   3.2050E+01  5.1006E+00  6.4631E-01  3.7724E+00  5.3025E+01 -2.5858E+00 -1.0122E+01 -3.3458E+00 -2.3128E+01 -3.4462E+01
            -9.2429E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1701.27805499387        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.8780E-01  9.4947E-01  1.2342E+00  1.0975E+00  1.1376E+00  9.5516E-01  8.9181E-01  3.4907E-01  8.3849E-01  1.0370E+00
             9.9322E-01
 PARAMETER:  8.7723E-02  4.8146E-02  3.1045E-01  1.9302E-01  2.2891E-01  5.4129E-02 -1.4499E-02 -9.5248E-01 -7.6151E-02  1.3631E-01
             9.3200E-02
 GRADIENT:   6.3314E+01  8.1654E+00 -1.7472E+00  4.5340E+01  3.2300E+01 -1.5273E+01 -5.1060E+00 -1.3978E+00 -2.7475E-01 -1.6890E+01
            -1.4906E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1703.41651993846        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.6862E-01  9.9605E-01  1.0667E+00  1.0495E+00  1.0652E+00  9.8594E-01  1.0225E+00  3.0841E-01  8.0371E-01  1.0119E+00
             1.0011E+00
 PARAMETER:  6.8119E-02  9.6043E-02  1.6456E-01  1.4829E-01  1.6318E-01  8.5840E-02  1.2228E-01 -1.0763E+00 -1.1851E-01  1.1185E-01
             1.0108E-01
 GRADIENT:   1.3968E+01  2.8073E-01 -5.4934E+00  1.4400E+01  1.0086E+01 -2.2890E+00 -1.1512E+00  2.8384E-02 -1.1183E+00 -1.8785E+00
            -1.6227E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1703.69046936484        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      367
 NPARAMETR:  9.7926E-01  1.0214E+00  1.0835E+00  1.0332E+00  1.0804E+00  1.0011E+00  1.0051E+00  3.2075E-01  8.2097E-01  1.0341E+00
             1.0057E+00
 PARAMETER:  7.9045E-02  1.2114E-01  1.8016E-01  1.3266E-01  1.7735E-01  1.0111E-01  1.0513E-01 -1.0371E+00 -9.7270E-02  1.3349E-01
             1.0565E-01
 GRADIENT:  -9.6942E-01 -1.1098E+00  8.5246E-01 -1.3488E+00 -8.0926E-01  2.0237E-01  3.4742E-02 -9.3163E-02 -1.5536E-01 -3.4265E-02
            -1.4835E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1703.82270081939        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      543
 NPARAMETR:  9.7946E-01  1.2503E+00  9.9713E-01  8.9357E-01  1.1556E+00  1.0009E+00  8.4759E-01  3.6049E-01  9.2653E-01  1.0587E+00
             1.0066E+00
 PARAMETER:  7.9246E-02  3.2339E-01  9.7125E-02 -1.2530E-02  2.4459E-01  1.0086E-01 -6.5364E-02 -9.2028E-01  2.3695E-02  1.5704E-01
             1.0656E-01
 GRADIENT:  -3.6640E+00  8.8679E+00 -1.1485E-02  1.0018E+01 -1.1514E+00 -2.8142E-01 -6.8083E-01  4.8575E-02 -8.3016E-01 -5.8351E-01
             2.3463E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1703.99829223663        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      719
 NPARAMETR:  9.8327E-01  1.4685E+00  8.8462E-01  7.4577E-01  1.2318E+00  1.0024E+00  7.5068E-01  2.7788E-01  1.0753E+00  1.0975E+00
             1.0053E+00
 PARAMETER:  8.3127E-02  4.8424E-01 -2.2593E-02 -1.9333E-01  3.0845E-01  1.0237E-01 -1.8677E-01 -1.1806E+00  1.7259E-01  1.9301E-01
             1.0532E-01
 GRADIENT:   3.0067E+00  2.5473E+00  1.0011E+00  1.3384E+00 -1.9571E+00  9.7417E-02  1.5766E-01  5.4552E-02  3.3912E-01  3.7954E-01
            -4.1510E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1704.05330032971        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      896
 NPARAMETR:  9.8192E-01  1.5931E+00  7.9851E-01  6.6313E-01  1.2709E+00  1.0021E+00  7.1218E-01  1.8948E-01  1.1672E+00  1.1072E+00
             1.0066E+00
 PARAMETER:  8.1755E-02  5.6571E-01 -1.2500E-01 -3.1079E-01  3.3972E-01  1.0209E-01 -2.3942E-01 -1.5635E+00  2.5460E-01  2.0181E-01
             1.0658E-01
 GRADIENT:  -1.1007E+00  2.2425E+00 -4.4113E-01  1.6862E+00 -1.0121E-02 -1.0890E-01 -1.2413E-01  5.6139E-02 -1.0681E-01 -3.0341E-02
             2.2534E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1704.06013315565        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  9.8229E-01  1.6282E+00  7.8054E-01  6.3820E-01  1.2866E+00  1.0022E+00  7.0076E-01  1.5854E-01  1.2041E+00  1.1163E+00
             1.0068E+00
 PARAMETER:  8.2128E-02  5.8745E-01 -1.4776E-01 -3.4910E-01  3.5203E-01  1.0218E-01 -2.5559E-01 -1.7418E+00  2.8572E-01  2.1004E-01
             1.0674E-01
 GRADIENT:  -3.9367E-01 -6.8813E-01 -4.7948E-02 -4.7976E-01 -1.1146E-01 -7.4937E-02  1.6133E-02  3.8816E-02  1.2396E-01  1.0003E-01
             1.6212E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1704.07947635063        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1250
 NPARAMETR:  9.8243E-01  1.6008E+00  7.9093E-01  6.5636E-01  1.2736E+00  1.0024E+00  7.0962E-01  5.0120E-02  1.1770E+00  1.1089E+00
             1.0063E+00
 PARAMETER:  8.2272E-02  5.7053E-01 -1.3454E-01 -3.2104E-01  3.4184E-01  1.0236E-01 -2.4302E-01 -2.8933E+00  2.6295E-01  2.0336E-01
             1.0626E-01
 GRADIENT:   1.4255E-02 -6.7053E-01 -2.4018E-01 -6.0927E-02  4.5626E-01 -1.3047E-02 -7.1066E-02  3.6181E-03 -2.7657E-03 -1.9753E-02
             2.7318E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1704.08130568733        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1425
 NPARAMETR:  9.8242E-01  1.6050E+00  7.8907E-01  6.5381E-01  1.2750E+00  1.0024E+00  7.0885E-01  1.0000E-02  1.1804E+00  1.1098E+00
             1.0063E+00
 PARAMETER:  8.2266E-02  5.7313E-01 -1.3691E-01 -3.2494E-01  3.4295E-01  1.0239E-01 -2.4411E-01 -4.5171E+00  2.6589E-01  2.0415E-01
             1.0630E-01
 GRADIENT:  -1.8090E-02  4.7658E-02  4.2029E-03  3.0370E-02  4.0169E-02 -2.3713E-03  9.8012E-04  0.0000E+00 -1.4391E-02 -9.0584E-03
            -4.5298E-03

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1704.08130568733        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1454
 NPARAMETR:  9.8246E-01  1.6049E+00  7.8903E-01  6.5375E-01  1.2749E+00  1.0024E+00  7.0884E-01  1.0000E-02  1.1808E+00  1.1099E+00
             1.0063E+00
 PARAMETER:  8.2266E-02  5.7313E-01 -1.3691E-01 -3.2494E-01  3.4295E-01  1.0239E-01 -2.4411E-01 -4.5171E+00  2.6589E-01  2.0415E-01
             1.0630E-01
 GRADIENT:  -3.2717E-02  8.2400E-02  4.3750E-03  3.2407E-02  4.0938E-02 -3.8213E-03  5.8869E-04  0.0000E+00 -1.4559E-02 -8.3201E-03
            -4.3020E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1454
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.4147E-04 -2.7159E-02 -2.6109E-04  1.6344E-02 -3.6251E-02
 SE:             2.9821E-02  2.1149E-02  1.0619E-04  2.2709E-02  2.3160E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9354E-01  1.9909E-01  1.3946E-02  4.7169E-01  1.1753E-01

 ETASHRINKSD(%)  9.4859E-02  2.9148E+01  9.9644E+01  2.3923E+01  2.2411E+01
 ETASHRINKVR(%)  1.8963E-01  4.9799E+01  9.9999E+01  4.2123E+01  3.9799E+01
 EBVSHRINKSD(%)  4.2741E-01  2.7740E+01  9.9682E+01  2.6010E+01  2.0162E+01
 EBVSHRINKVR(%)  8.5299E-01  4.7785E+01  9.9999E+01  4.5255E+01  3.6259E+01
 RELATIVEINF(%)  9.8949E+01  1.8566E+00  1.1695E-04  2.1091E+00  1.3108E+01
 EPSSHRINKSD(%)  4.2195E+01
 EPSSHRINKVR(%)  6.6585E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1704.0813056873330     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -968.93047912359486     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.68
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.83
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1704.081       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  1.61E+00  7.89E-01  6.54E-01  1.28E+00  1.00E+00  7.09E-01  1.00E-02  1.18E+00  1.11E+00  1.01E+00
 


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
+        1.13E+03
 
 TH 2
+       -7.86E+00  4.71E+02
 
 TH 3
+        6.68E+00  1.16E+02  2.06E+02
 
 TH 4
+       -1.66E+01  4.90E+02 -1.55E+02  1.03E+03
 
 TH 5
+       -2.33E+00 -1.52E+02 -1.91E+02  1.44E+02  3.55E+02
 
 TH 6
+        2.45E+00 -1.22E+00  7.48E-01 -5.75E+00 -2.09E+00  1.98E+02
 
 TH 7
+       -2.27E-02  3.46E+00  1.52E+01 -1.90E+01 -1.48E+01 -1.21E+00  1.14E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.23E-01 -1.86E+01 -2.45E+01  4.79E+01  5.17E+00 -4.62E-01  3.26E+01  0.00E+00  5.00E+01
 
 TH10
+       -1.73E+00 -1.00E+01 -2.26E+01 -3.52E+00 -4.78E+01 -1.65E-02  9.92E+00  0.00E+00  4.82E+00  7.19E+01
 
 TH11
+       -8.92E+00 -2.68E+01 -3.56E+01  5.08E+00  2.18E+00  3.20E+00  1.15E+01  0.00E+00  9.34E+00  1.63E+01  2.26E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.572
Stop Time:
Sat Sep 25 13:16:11 CDT 2021
