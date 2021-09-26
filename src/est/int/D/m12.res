Sat Sep 25 05:23:23 CDT 2021
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
$DATA ../../../../data/int/D/dat12.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m12.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3333.82713074705        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.3186E+02 -1.2384E+02 -9.1506E+01 -2.2100E+02  1.2426E+02 -4.0524E+02 -2.4942E+02 -2.1029E+01 -4.4943E+02 -9.3820E+01
            -1.4645E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3644.21974340909        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  1.1303E+00  1.3137E+00  1.3356E+00  1.0958E+00  1.2649E+00  1.4880E+00  1.9546E+00  1.1722E+00  1.6086E+00  1.1464E+00
             1.1713E+00
 PARAMETER:  2.2246E-01  3.7282E-01  3.8942E-01  1.9150E-01  3.3496E-01  4.9741E-01  7.7018E-01  2.5888E-01  5.7535E-01  2.3662E-01
             2.5813E-01
 GRADIENT:   1.1407E+02  2.7753E+01 -3.2526E+01  1.2278E+02  5.6700E+01 -3.7746E+01 -2.5445E+01 -2.0966E+00  9.0974E+00  9.5985E+00
             1.5909E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3653.45492139248        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  1.1035E+00  1.4472E+00  1.8271E+00  1.0395E+00  1.4827E+00  1.7507E+00  2.4384E+00  1.4823E+00  1.3586E+00  1.0715E+00
             1.1744E+00
 PARAMETER:  1.9853E-01  4.6964E-01  7.0272E-01  1.3876E-01  4.9388E-01  6.6000E-01  9.9136E-01  4.9361E-01  4.0642E-01  1.6906E-01
             2.6079E-01
 GRADIENT:   7.7307E+01  6.1438E+01 -1.0318E+01  1.2607E+02  1.3677E+02  5.3370E+01  4.1579E+01 -4.7726E+00 -8.2690E+00 -3.0839E+01
             1.4128E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3672.20340726302        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      359
 NPARAMETR:  1.0951E+00  1.4438E+00  1.5710E+00  9.1632E-01  1.3698E+00  1.7528E+00  2.5721E+00  1.2680E+00  1.3414E+00  1.1624E+00
             1.1038E+00
 PARAMETER:  1.9080E-01  4.6727E-01  5.5171E-01  1.2616E-02  4.1464E-01  6.6123E-01  1.0447E+00  3.3744E-01  3.9372E-01  2.5049E-01
             1.9878E-01
 GRADIENT:   7.8352E+00  1.4678E+01  8.7777E+00  1.6978E+01  1.9410E+01  1.6037E+00  2.1355E+01 -6.8565E+00 -8.2412E+00 -3.2868E+00
             2.4745E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3675.01772480295        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  1.0797E+00  1.3817E+00  1.4505E+00  9.0184E-01  1.3220E+00  1.7477E+00  2.3536E+00  1.2964E+00  1.4105E+00  1.1567E+00
             1.0875E+00
 PARAMETER:  1.7666E-01  4.2329E-01  4.7189E-01 -3.3180E-03  3.7917E-01  6.5828E-01  9.5595E-01  3.5962E-01  4.4397E-01  2.4561E-01
             1.8384E-01
 GRADIENT:  -2.3286E+00 -1.8842E+00  1.3103E+00  1.6886E+00 -9.6059E-01  5.5971E-01 -2.5816E+00 -3.4892E+00 -5.5648E-01  2.2595E-01
             3.9356E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3675.08126865127        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      712
 NPARAMETR:  1.0805E+00  1.3853E+00  1.4567E+00  9.0238E-01  1.3241E+00  1.7476E+00  2.3567E+00  1.3222E+00  1.4077E+00  1.1575E+00
             1.0866E+00
 PARAMETER:  1.7739E-01  4.2593E-01  4.7618E-01 -2.7194E-03  3.8073E-01  6.5822E-01  9.5725E-01  3.7928E-01  4.4195E-01  2.4625E-01
             1.8303E-01
 GRADIENT:  -1.8009E+00 -9.8854E-01  5.6726E-01  2.4597E+00 -1.2422E+00  5.4363E-01 -2.2556E+00 -2.2719E+00 -7.6731E-01 -1.3835E-02
             2.8118E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3675.11584369952        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      876
 NPARAMETR:  1.0810E+00  1.3873E+00  1.4600E+00  9.0127E-01  1.3262E+00  1.7453E+00  2.3590E+00  1.3426E+00  1.4079E+00  1.1590E+00
             1.0863E+00
 PARAMETER:  1.7786E-01  4.2737E-01  4.7845E-01 -3.9548E-03  3.8235E-01  6.5692E-01  9.5824E-01  3.9460E-01  4.4211E-01  2.4754E-01
             1.8275E-01
 GRADIENT:   6.5807E+01  3.8111E+01  1.0279E+00  9.1470E+00  1.1642E+01  6.2350E+01  5.0973E+01 -4.8459E-01  4.6671E+00  9.1214E-01
             2.6573E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3675.12366796521        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  1.0810E+00  1.3873E+00  1.4604E+00  9.0127E-01  1.3305E+00  1.7417E+00  2.3753E+00  1.3479E+00  1.4079E+00  1.1591E+00
             1.0842E+00
 PARAMETER:  1.7786E-01  4.2737E-01  4.7869E-01 -3.9548E-03  3.8552E-01  6.5488E-01  9.6512E-01  3.9852E-01  4.4211E-01  2.4766E-01
             1.8083E-01
 GRADIENT:   6.6062E+01  3.8084E+01 -3.4421E-01  1.0155E+01  1.5100E+01  6.1577E+01  5.3958E+01 -9.4631E-02  4.9255E+00  3.7834E-01
            -1.5412E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3675.12712704840        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1244             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0815E+00  1.3873E+00  1.4603E+00  9.0127E-01  1.3293E+00  1.7442E+00  2.3738E+00  1.3488E+00  1.4080E+00  1.1591E+00
             1.0848E+00
 PARAMETER:  1.7838E-01  4.2737E-01  4.7866E-01 -3.9548E-03  3.8465E-01  6.5630E-01  9.6449E-01  3.9921E-01  4.4213E-01  2.4768E-01
             1.8139E-01
 GRADIENT:   6.6563E+01  3.8099E+01 -3.3099E-01  9.9188E+00  1.3959E+01  6.2176E+01  5.3626E+01  1.3718E-01  4.8997E+00  4.6990E-01
            -3.7566E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3675.12877677793        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1443
 NPARAMETR:  1.0815E+00  1.3872E+00  1.4612E+00  9.0126E-01  1.3293E+00  1.7447E+00  2.3736E+00  1.3494E+00  1.4082E+00  1.1596E+00
             1.0848E+00
 PARAMETER:  1.7836E-01  4.2731E-01  4.7929E-01 -3.9573E-03  3.8464E-01  6.5660E-01  9.6439E-01  3.9966E-01  4.4233E-01  2.4805E-01
             1.8143E-01
 GRADIENT:  -1.1050E+00 -1.0492E+00 -1.8984E+00  2.9165E+00  3.5744E-01 -9.2867E-02  1.0098E-01 -1.6304E-01 -3.2385E-01 -5.5084E-01
            -5.7217E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -3675.12884577242        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1500
 NPARAMETR:  1.0815E+00  1.3872E+00  1.4612E+00  9.0126E-01  1.3293E+00  1.7449E+00  2.3736E+00  1.3499E+00  1.4082E+00  1.1596E+00
             1.0848E+00
 PARAMETER:  1.7836E-01  4.2731E-01  4.7929E-01 -3.9573E-03  3.8464E-01  6.5670E-01  9.6439E-01  4.0006E-01  4.4233E-01  2.4805E-01
             1.8143E-01
 GRADIENT:  -1.1069E+00 -1.0553E+00 -1.9836E+00  2.9326E+00  3.0464E-01 -5.0700E-02  1.0191E-01 -7.8782E-02 -3.2542E-01 -5.6409E-01
            -5.9290E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1500
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0796E-03 -1.5818E-02 -2.9224E-02  2.1797E-02 -4.0433E-02
 SE:             2.9974E-02  2.7666E-02  1.6489E-02  2.5532E-02  2.3103E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7127E-01  5.6749E-01  7.6337E-02  3.9326E-01  8.0099E-02

 ETASHRINKSD(%)  1.0000E-10  7.3144E+00  4.4760E+01  1.4464E+01  2.2602E+01
 ETASHRINKVR(%)  1.0000E-10  1.4094E+01  6.9485E+01  2.6835E+01  4.0095E+01
 EBVSHRINKSD(%)  1.0273E-01  7.4152E+00  5.0235E+01  1.8639E+01  2.0149E+01
 EBVSHRINKVR(%)  2.0535E-01  1.4280E+01  7.5234E+01  3.3804E+01  3.6238E+01
 RELATIVEINF(%)  9.9794E+01  6.9223E+01  2.0142E+01  4.8729E+01  4.2176E+01
 EPSSHRINKSD(%)  2.1641E+01
 EPSSHRINKVR(%)  3.8598E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3675.1288457724245     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2021.0394860040137     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    47.23
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.60
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3675.129       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.08E+00  1.39E+00  1.46E+00  9.01E-01  1.33E+00  1.74E+00  2.37E+00  1.35E+00  1.41E+00  1.16E+00  1.08E+00
 


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
+        1.42E+09
 
 TH 2
+       -2.47E+08  1.39E+08
 
 TH 3
+        5.35E+04  1.74E+04  4.98E+07
 
 TH 4
+       -1.63E+09  5.29E+08 -6.05E+07  3.01E+09
 
 TH 5
+       -3.32E+03 -1.12E+03  7.89E+07 -6.99E+03  1.08E+08
 
 TH 6
+        3.31E+02  1.07E+02 -2.85E+00  7.07E+02  1.28E+02  6.56E+01
 
 TH 7
+        9.64E+02  3.19E+02 -2.38E+06 -1.37E+08 -2.09E+07 -2.86E+01  1.01E+07
 
 TH 8
+       -2.73E+05 -8.88E+04  7.47E+07 -5.84E+05 -1.03E+05  1.21E+02  2.29E+04  9.68E+07
 
 TH 9
+       -6.93E+04 -2.26E+04 -8.74E+06 -6.81E+07  8.87E+07 -2.84E+00 -2.68E+06 -2.67E+02  6.30E+07
 
 TH10
+        1.65E+02  4.74E+01 -1.89E+07 -1.47E+08  1.92E+08 -6.19E+00 -5.80E+06 -7.22E+02 -2.13E+07  2.95E+08
 
 TH11
+       -3.56E+04 -2.42E+08  2.77E+07  6.89E+08 -2.81E+08 -3.24E+02  6.27E+07  2.67E+05  3.11E+07  6.74E+07  6.31E+08
 
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
 #CPUT: Total CPU Time in Seconds,       62.948
Stop Time:
Sat Sep 25 05:24:27 CDT 2021
