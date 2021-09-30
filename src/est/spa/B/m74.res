Wed Sep 29 11:30:09 CDT 2021
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
$DATA ../../../../data/spa/B/dat74.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1691.60143315052        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3202E+02 -2.1934E+01 -1.4754E+01 -8.2179E+00 -1.0401E+01  5.6969E+01 -7.5844E+00  1.4405E+01  4.0563E-01  1.6740E+01
             2.3375E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1697.84390642895        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9265E-01  1.0580E+00  1.0854E+00  1.0196E+00  1.0623E+00  9.3713E-01  1.0463E+00  9.4113E-01  1.0151E+00  9.4986E-01
             9.3623E-01
 PARAMETER:  9.2621E-02  1.5638E-01  1.8191E-01  1.1945E-01  1.6042E-01  3.5068E-02  1.4525E-01  3.9324E-02  1.1499E-01  4.8561E-02
             3.4101E-02
 GRADIENT:  -4.4854E-01  9.2970E-01  1.9665E-01 -9.4958E-01 -8.2475E+00 -4.8905E+00 -6.6126E+00  6.2866E+00 -1.6658E+00 -4.2600E+00
            -9.8804E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1699.21225896334        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.8769E-01  1.0140E+00  1.1581E+00  1.0524E+00  1.1043E+00  9.5690E-01  1.1623E+00  7.1741E-01  9.7786E-01  1.0682E+00
             9.6046E-01
 PARAMETER:  8.7612E-02  1.1388E-01  2.4675E-01  1.5103E-01  1.9924E-01  5.5940E-02  2.5044E-01 -2.3211E-01  7.7608E-02  1.6600E-01
             5.9658E-02
 GRADIENT:  -1.0602E+01  1.5505E+00  3.9576E-01  1.2439E+00  7.6131E+00  3.9750E+00 -6.9246E-01 -9.2861E-02 -2.2506E+00  2.2936E+00
             5.6822E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.12985145474        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  9.9238E-01  8.3742E-01  1.0321E+00  1.1561E+00  9.4630E-01  9.4520E-01  1.3772E+00  5.0545E-01  9.0605E-01  9.2423E-01
             9.5576E-01
 PARAMETER:  9.2355E-02 -7.7435E-02  1.3164E-01  2.4504E-01  4.4805E-02  4.3639E-02  4.2005E-01 -5.8231E-01  1.3420E-03  2.1207E-02
             5.4747E-02
 GRADIENT:   1.4865E+00  5.9134E+00  3.6322E+00  5.0649E+00 -4.7095E+00 -6.2184E-01 -2.8580E-02 -2.6265E-01 -4.6023E-01 -9.1174E-01
            -9.3648E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1700.29948747660        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      710
 NPARAMETR:  9.8954E-01  6.2965E-01  1.0859E+00  1.2841E+00  8.9330E-01  9.4360E-01  1.6674E+00  4.9644E-01  8.4883E-01  9.2930E-01
             9.5740E-01
 PARAMETER:  8.9480E-02 -3.6259E-01  1.8241E-01  3.5006E-01 -1.2833E-02  4.1952E-02  6.1126E-01 -6.0029E-01 -6.3893E-02  2.6679E-02
             5.6470E-02
 GRADIENT:   3.8759E-01  4.6711E+00  5.2589E+00  5.1038E+00 -6.9316E+00 -1.6313E-01 -3.5062E-01 -1.0331E+00 -5.9045E-01  4.8776E-01
            -7.8876E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1700.52278679432        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  9.8502E-01  4.0323E-01  1.1884E+00  1.4300E+00  8.6617E-01  9.4105E-01  2.1986E+00  6.4990E-01  8.0462E-01  9.4639E-01
             9.5728E-01
 PARAMETER:  8.4906E-02 -8.0826E-01  2.7261E-01  4.5771E-01 -4.3674E-02  3.9242E-02  8.8783E-01 -3.3094E-01 -1.1738E-01  4.4896E-02
             5.6344E-02
 GRADIENT:  -1.4972E+00  5.2776E+00  4.2273E+00  1.2175E+01 -1.1161E+01  3.0096E-01 -1.7672E-01 -9.4610E-02 -5.2196E-01  2.7398E+00
             4.6905E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1700.67608814474        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  9.8251E-01  2.6919E-01  1.2669E+00  1.5171E+00  8.6294E-01  9.3856E-01  2.7720E+00  7.6367E-01  7.8299E-01  9.4626E-01
             9.5660E-01
 PARAMETER:  8.2351E-02 -1.2123E+00  3.3655E-01  5.1677E-01 -4.7414E-02  3.6590E-02  1.1196E+00 -1.6962E-01 -1.4464E-01  4.4762E-02
             5.5625E-02
 GRADIENT:  -1.4440E+00  4.5618E+00  2.0352E+00  1.6367E+01 -8.5273E+00  2.0705E-01  3.6093E-01  4.8068E-01 -8.5535E-01  1.6491E+00
             3.6586E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1700.88733316498        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  9.8057E-01  1.5563E-01  1.3691E+00  1.5934E+00  8.7478E-01  9.3606E-01  3.6801E+00  8.8576E-01  7.6548E-01  9.4674E-01
             9.5637E-01
 PARAMETER:  8.0378E-02 -1.7602E+00  4.1418E-01  5.6587E-01 -3.3783E-02  3.3921E-02  1.4029E+00 -2.1312E-02 -1.6725E-01  4.5271E-02
             5.5386E-02
 GRADIENT:  -8.8449E-01  3.4698E+00  9.0728E-01  2.0219E+01 -4.9797E+00 -7.6261E-02  8.3287E-01  3.9984E-01 -1.2015E+00 -3.4835E-01
            -2.2244E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1701.47141210368        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1413
 NPARAMETR:  9.7944E-01  5.7476E-02  1.5712E+00  1.6584E+00  9.2542E-01  9.3438E-01  5.6884E+00  1.0749E+00  7.4082E-01  9.8464E-01
             9.5754E-01
 PARAMETER:  7.9227E-02 -2.7564E+00  5.5182E-01  6.0588E-01  2.2490E-02  3.2127E-02  1.8384E+00  1.7224E-01 -1.9999E-01  8.4520E-02
             5.6614E-02
 GRADIENT:   3.0063E-01 -2.1358E-02  2.0833E+00  1.2897E+01 -2.1766E+00 -6.6408E-02 -1.6387E+00  2.5683E-01 -5.6429E-01  2.5871E-01
            -2.8422E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1701.73774943957        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1589
 NPARAMETR:  9.7905E-01  4.9537E-02  1.4898E+00  1.6530E+00  8.9670E-01  9.3439E-01  5.9418E+00  1.0000E+00  7.4329E-01  9.7034E-01
             9.5737E-01
 PARAMETER:  7.8831E-02 -2.9050E+00  4.9862E-01  6.0258E-01 -9.0318E-03  3.2138E-02  1.8820E+00  1.0003E-01 -1.9667E-01  6.9893E-02
             5.6439E-02
 GRADIENT:  -2.9263E-02 -1.2243E-02 -5.0753E-02  1.5148E+00 -2.1269E-01 -3.1798E-02 -2.2288E-01  1.2625E-01  1.2252E-01  2.4837E-01
             1.0375E-02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1701.91060095558        NO. OF FUNC. EVALS.:  67
 CUMULATIVE NO. OF FUNC. EVALS.:     1656
 NPARAMETR:  9.7888E-01  3.9551E-02  1.4186E+00  1.6483E+00  8.6911E-01  9.3535E-01  6.3355E+00  9.1677E-01  7.4211E-01  9.4804E-01
             9.5901E-01
 PARAMETER:  7.8617E-02 -3.1295E+00  4.4981E-01  5.9964E-01 -4.0291E-02  3.3250E-02  1.9465E+00  1.2103E-02 -1.9806E-01  4.7017E-02
             5.8288E-02
 GRADIENT:  -9.3774E-01  1.1955E+01  1.2428E+00 -5.9663E+01 -6.2535E-02  3.0363E-01  1.5625E+01 -4.1095E-01  6.6752E-01  5.5916E-01
             5.1293E-01
 NUMSIGDIG:         2.1         2.4         2.4         2.5         3.4         1.8         2.5         0.8         1.7         1.2
                    1.6

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1656
 NO. OF SIG. DIGITS IN FINAL EST.:  0.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.9919E-04  1.2870E-02 -2.3261E-02 -1.4633E-02 -3.0654E-02
 SE:             2.9803E-02  9.1103E-03  1.5099E-02  2.8695E-02  2.1646E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7861E-01  1.5776E-01  1.2342E-01  6.1009E-01  1.5673E-01

 ETASHRINKSD(%)  1.5596E-01  6.9480E+01  4.9417E+01  3.8687E+00  2.7483E+01
 ETASHRINKVR(%)  3.1169E-01  9.0685E+01  7.4414E+01  7.5878E+00  4.7414E+01
 EBVSHRINKSD(%)  4.3500E-01  8.0475E+01  5.1622E+01  3.6682E+00  2.4223E+01
 EBVSHRINKVR(%)  8.6810E-01  9.6188E+01  7.6596E+01  7.2018E+00  4.2578E+01
 RELATIVEINF(%)  9.8971E+01  1.9250E+00  4.1437E+00  4.8932E+01  9.9674E+00
 EPSSHRINKSD(%)  4.3936E+01
 EPSSHRINKVR(%)  6.8568E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1701.9106009555842     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -966.75977439184601     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.53
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1701.911       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  3.96E-02  1.42E+00  1.65E+00  8.69E-01  9.35E-01  6.34E+00  9.16E-01  7.42E-01  9.48E-01  9.59E-01
 


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
+       -2.69E+01  7.52E+04
 
 TH 3
+       -1.61E+00 -7.50E+01  1.64E+02
 
 TH 4
+       -8.96E+00 -9.11E+03  1.85E+03  1.93E+03
 
 TH 5
+        2.15E+00  8.70E+01 -3.40E+02 -1.37E+04  9.76E+02
 
 TH 6
+        1.99E+00 -6.83E+00  5.88E-01 -2.26E+00 -1.05E+00  2.23E+02
 
 TH 7
+        3.08E-02  7.94E+02 -9.48E-01 -4.03E+00  2.25E+00 -9.76E-03  8.51E+00
 
 TH 8
+        2.15E-02 -6.50E+01 -2.67E+01  2.63E+00 -5.69E+00  3.24E-01 -4.15E-01  2.77E+01
 
 TH 9
+        3.01E+00 -6.30E+01  1.91E+01  8.08E+03 -2.14E+01 -2.96E-01  1.17E-01  4.98E+00  3.40E+02
 
 TH10
+        7.53E-01 -2.47E+02  2.48E+00  1.25E+04 -9.81E+01  8.09E-01 -1.71E+00  2.35E+01  1.71E+01  9.41E+01
 
 TH11
+       -1.01E+01 -2.66E+01 -1.75E+01 -9.50E+00 -3.38E+00  2.00E+00 -7.47E-02  1.64E+01  1.35E+01  2.23E+01  2.30E+02
 
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
 #CPUT: Total CPU Time in Seconds,       29.408
Stop Time:
Wed Sep 29 11:30:40 CDT 2021
