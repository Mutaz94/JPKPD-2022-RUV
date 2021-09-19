Sat Sep 18 09:22:08 CDT 2021
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
$DATA ../../../../data/spa/A1/dat66.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1295.63412895042        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5671E+02 -4.1607E+01  1.7420E+01 -7.9130E+01  4.6013E+01  1.4084E+01 -3.1489E+01  2.5766E+00 -3.4772E+01 -7.4532E+01
            -5.4266E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1432.36060883883        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0358E+00  9.9853E-01  1.1978E+00  1.0705E+00  1.0989E+00  1.0076E+00  1.4347E+00  7.3679E-01  1.0624E+00  1.7451E+00
             2.1159E+00
 PARAMETER:  1.3522E-01  9.8527E-02  2.8050E-01  1.6811E-01  1.9434E-01  1.0758E-01  4.6098E-01 -2.0546E-01  1.6057E-01  6.5678E-01
             8.4950E-01
 GRADIENT:   1.8261E+02 -7.9031E+00 -8.7213E+00 -2.2714E+01 -2.1343E+01  1.1888E+01  1.5710E+01  3.0964E+00  1.5516E+01  3.7037E+01
             7.6507E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1446.87398602542        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.8457E-01  6.9117E-01  1.0018E+00  1.3236E+00  8.3863E-01  9.8697E-01  1.6984E+00  1.3484E-01  8.5788E-01  1.0994E+00
             2.0377E+00
 PARAMETER:  8.4450E-02 -2.6936E-01  1.0179E-01  3.8038E-01 -7.5991E-02  8.6888E-02  6.2967E-01 -1.9037E+00 -5.3293E-02  1.9480E-01
             8.1180E-01
 GRADIENT:   7.4385E+01  3.7554E+01  2.6527E+00  9.7148E+01 -1.0280E+01  1.5072E+01  1.3059E+00  1.5730E-01 -8.2652E+00 -6.4383E+00
             4.5408E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1454.12538039755        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.4641E-01  7.1688E-01  1.0525E+00  1.2420E+00  8.9965E-01  9.3942E-01  1.5188E+00  1.4269E-01  9.2547E-01  1.2402E+00
             1.8112E+00
 PARAMETER:  4.4924E-02 -2.3284E-01  1.5117E-01  3.1670E-01 -5.7461E-03  3.7510E-02  5.1793E-01 -1.8471E+00  2.2548E-02  3.1526E-01
             6.9399E-01
 GRADIENT:  -2.3864E+00  3.8436E+00 -1.7478E+00  6.6126E+00  2.9721E-01 -4.6222E-01 -1.1324E+00  1.5574E-01 -2.9936E-01 -3.6106E-01
            -6.1285E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1455.42702049271        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.4416E-01  4.6391E-01  1.0525E+00  1.3841E+00  8.1154E-01  9.3784E-01  2.0050E+00  3.8516E-02  8.6962E-01  1.1864E+00
             1.8172E+00
 PARAMETER:  4.2535E-02 -6.6807E-01  1.5119E-01  4.2503E-01 -1.0882E-01  3.5827E-02  7.9563E-01 -3.1567E+00 -3.9695E-02  2.7089E-01
             6.9729E-01
 GRADIENT:  -1.1260E+00  9.1601E-01 -9.2438E-01  3.3813E+00  1.5159E+00 -2.5151E-01 -1.5728E-01  1.1566E-02  9.6668E-01 -8.2801E-01
             2.4359E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1456.02596724602        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.4219E-01  3.1850E-01  1.0482E+00  1.4647E+00  7.6395E-01  9.3757E-01  2.4639E+00  1.0000E-02  8.3767E-01  1.1751E+00
             1.8191E+00
 PARAMETER:  4.0453E-02 -1.0441E+00  1.4709E-01  4.8164E-01 -1.6926E-01  3.5532E-02  1.0018E+00 -4.5986E+00 -7.7136E-02  2.6133E-01
             6.9834E-01
 GRADIENT:  -6.0934E-01  1.1461E+00  1.2070E+00  2.9182E+00 -2.3584E+00  3.1760E-01  1.4316E-01  0.0000E+00 -1.0285E+00  4.3013E-01
            -3.5902E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1456.60373998326        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  9.4001E-01  1.9612E-01  1.0366E+00  1.5342E+00  7.2350E-01  9.3507E-01  3.1840E+00  1.0000E-02  8.1728E-01  1.1565E+00
             1.8233E+00
 PARAMETER:  3.8135E-02 -1.5290E+00  1.3598E-01  5.2802E-01 -2.2365E-01  3.2866E-02  1.2581E+00 -6.7072E+00 -1.0177E-01  2.4539E-01
             7.0063E-01
 GRADIENT:  -9.8600E-01  1.5355E+00  3.7167E+00  1.1040E+01 -6.6921E+00 -1.8353E-01  4.2867E-02  0.0000E+00 -1.2302E+00  6.9291E-01
            -6.2009E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1457.48441382159        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  9.3779E-01  7.1619E-02  1.0338E+00  1.6067E+00  6.9929E-01  9.3364E-01  4.9266E+00  1.0000E-02  8.1409E-01  1.1469E+00
             1.8320E+00
 PARAMETER:  3.5774E-02 -2.5364E+00  1.3322E-01  5.7419E-01 -2.5768E-01  3.1337E-02  1.6946E+00 -1.1644E+01 -1.0569E-01  2.3707E-01
             7.0541E-01
 GRADIENT:  -1.0959E+01 -1.6877E+00 -2.0656E+00  3.7677E+00 -1.1890E+00 -1.1871E+00 -5.9345E+00  0.0000E+00  5.8842E+00  1.4870E+00
             8.1045E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1458.77165532940        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  9.3928E-01  1.3520E-02  1.1350E+00  1.6644E+00  7.3146E-01  9.2964E-01  1.0760E+01  1.0000E-02  7.9235E-01  1.1752E+00
             1.8402E+00
 PARAMETER:  3.7362E-02 -4.2036E+00  2.2667E-01  6.0946E-01 -2.1271E-01  2.7047E-02  2.4759E+00 -1.9866E+01 -1.3275E-01  2.6146E-01
             7.0989E-01
 GRADIENT:  -4.2535E+00 -6.8053E-01  2.7040E+00  3.1956E+01 -5.6017E+00 -2.5119E+00 -2.5942E+00  0.0000E+00  2.7010E+00 -7.8332E-01
            -2.9549E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1459.16106840453        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  9.3992E-01  1.0000E-02  1.1066E+00  1.6512E+00  7.2079E-01  9.3275E-01  1.3206E+01  1.0000E-02  7.8926E-01  1.1630E+00
             1.8389E+00
 PARAMETER:  3.8042E-02 -4.6150E+00  2.0128E-01  6.0152E-01 -2.2741E-01  3.0379E-02  2.6807E+00 -2.2096E+01 -1.3666E-01  2.5098E-01
             7.0918E-01
 GRADIENT:  -1.7462E+00  0.0000E+00 -2.0179E-01  9.2261E+00 -2.0252E-01 -1.1376E+00 -9.0709E-01  0.0000E+00  1.4948E+00 -7.8523E-01
            -7.3541E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1459.18226637498        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  9.4053E-01  1.0000E-02  1.1149E+00  1.6477E+00  7.2483E-01  9.3561E-01  1.3229E+01  1.0000E-02  7.8607E-01  1.1754E+00
             1.8379E+00
 PARAMETER:  3.8687E-02 -4.6072E+00  2.0877E-01  5.9940E-01 -2.2181E-01  3.3442E-02  2.6824E+00 -2.2047E+01 -1.4071E-01  2.6161E-01
             7.0863E-01
 GRADIENT:   7.9693E-02  0.0000E+00  9.2923E-02 -4.2414E-02 -1.6493E-01  6.9682E-02 -4.0791E-01  0.0000E+00  1.3858E-01  2.3168E-01
             3.6503E-02

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1459.18238724760        NO. OF FUNC. EVALS.:  88
 CUMULATIVE NO. OF FUNC. EVALS.:     1167
 NPARAMETR:  9.4050E-01  1.0000E-02  1.1145E+00  1.6479E+00  7.2468E-01  9.3546E-01  1.3225E+01  1.0000E-02  7.8611E-01  1.1745E+00
             1.8380E+00
 PARAMETER:  3.8655E-02 -4.6082E+00  2.0844E-01  5.9951E-01 -2.2206E-01  3.3274E-02  2.6826E+00 -2.2052E+01 -1.4052E-01  2.6108E-01
             7.0867E-01
 GRADIENT:  -1.6349E-03  0.0000E+00  5.2682E-02 -6.9251E-02 -1.0477E-01 -7.0096E-03  2.9391E+02  0.0000E+00  1.0804E-01  9.5100E-02
             2.6039E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1167
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.8168E-05  5.4215E-03 -5.9519E-05 -1.2494E-02 -2.9079E-02
 SE:             2.9429E-02  4.5016E-03  1.3496E-04  2.7676E-02  2.3011E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9788E-01  2.2846E-01  6.5921E-01  6.5167E-01  2.0634E-01

 ETASHRINKSD(%)  1.4091E+00  8.4919E+01  9.9548E+01  7.2818E+00  2.2910E+01
 ETASHRINKVR(%)  2.7983E+00  9.7726E+01  9.9998E+01  1.4033E+01  4.0572E+01
 EBVSHRINKSD(%)  1.4650E+00  8.8038E+01  9.9517E+01  6.4682E+00  2.0100E+01
 EBVSHRINKVR(%)  2.9086E+00  9.8569E+01  9.9998E+01  1.2518E+01  3.6160E+01
 RELATIVEINF(%)  9.6912E+01  1.0066E+00  2.0521E-04  4.8541E+01  5.6355E+00
 EPSSHRINKSD(%)  3.6588E+01
 EPSSHRINKVR(%)  5.9789E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1459.1823872476000     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -724.03156068386181     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.93
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1459.182       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.40E-01  1.00E-02  1.11E+00  1.65E+00  7.25E-01  9.35E-01  1.32E+01  1.00E-02  7.86E-01  1.17E+00  1.84E+00
 


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
+        1.38E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.94E+01  0.00E+00  2.28E+02
 
 TH 4
+       -1.48E+02  0.00E+00  8.59E+05  2.02E+05
 
 TH 5
+        2.17E+01  0.00E+00 -4.61E+02 -1.24E+06  1.11E+03
 
 TH 6
+       -1.84E+01  0.00E+00 -1.62E+01  3.08E+00  1.43E+01  2.08E+02
 
 TH 7
+       -1.94E+00  0.00E+00  1.03E+01  3.46E+01 -3.38E+01  1.81E-01  1.57E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.62E+00  0.00E+00  6.00E+01  1.81E+06 -1.78E+02  7.17E-02  1.27E+01  0.00E+00  4.32E+02
 
 TH10
+       -1.50E+01  0.00E+00  2.43E+01  6.51E+05 -1.56E+02 -1.50E+01  1.43E+01  0.00E+00  9.69E+01  1.27E+02
 
 TH11
+       -1.29E+01  0.00E+00 -6.58E+00  1.53E+05 -5.82E+00  3.64E-01  1.29E+00  0.00E+00  2.65E+01  2.25E+01  7.97E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.091
Stop Time:
Sat Sep 18 09:22:31 CDT 2021
