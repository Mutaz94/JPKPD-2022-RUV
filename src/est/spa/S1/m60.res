Sat Sep 25 09:59:42 CDT 2021
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
$DATA ../../../../data/spa/S1/dat60.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1662.35579681627        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5217E+01 -7.8104E+01 -3.5831E+01 -5.6310E+01  4.0779E+01  2.1886E+00 -2.2006E+01  5.5775E+00 -3.8502E+00 -1.2888E+01
             6.2979E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1669.76618046476        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      113
 NPARAMETR:  1.0087E+00  1.1953E+00  1.1868E+00  9.1805E-01  1.1204E+00  9.9127E-01  1.2792E+00  9.5087E-01  1.0084E+00  1.1817E+00
             9.8297E-01
 PARAMETER:  1.0869E-01  2.7838E-01  2.7122E-01  1.4502E-02  2.1368E-01  9.1232E-02  3.4624E-01  4.9620E-02  1.0834E-01  2.6698E-01
             8.2820E-02
 GRADIENT:   1.0929E+01 -1.2587E+00  2.4470E+01 -3.3821E+01 -3.0684E+01 -1.0374E+01  1.0175E+01 -5.5094E+00  3.4253E+00  1.7463E+00
            -4.3698E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1670.57048383242        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      292
 NPARAMETR:  1.0100E+00  1.3227E+00  1.1814E+00  8.6290E-01  1.1660E+00  1.0306E+00  1.1461E+00  1.1208E+00  1.0493E+00  1.1969E+00
             9.6993E-01
 PARAMETER:  1.0995E-01  3.7966E-01  2.6670E-01 -4.7459E-02  2.5357E-01  1.3011E-01  2.3640E-01  2.1407E-01  1.4811E-01  2.7970E-01
             6.9472E-02
 GRADIENT:   1.1338E+01  2.1386E+01  2.1556E+01 -5.9233E+00 -2.9277E+01  4.7104E+00  4.1101E+00 -5.9655E+00 -2.1014E-01  9.0930E-01
            -9.7987E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1670.65614789711        NO. OF FUNC. EVALS.: 214
 CUMULATIVE NO. OF FUNC. EVALS.:      506             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0098E+00  1.3239E+00  1.1889E+00  8.6234E-01  1.1697E+00  1.0067E+00  1.0828E+00  1.1432E+00  1.0505E+00  1.1987E+00
             9.7041E-01
 PARAMETER:  1.0974E-01  3.8061E-01  2.7307E-01 -4.8101E-02  2.5676E-01  1.0666E-01  1.7952E-01  2.3385E-01  1.4924E-01  2.8126E-01
             6.9958E-02
 GRADIENT:   6.4927E+01  4.5337E+01  2.0917E+01  2.9154E+00 -2.6427E+01  5.7529E+00 -1.8522E+00 -6.1186E+00 -2.0699E+00  9.2385E-01
            -1.0155E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1670.72279309332        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      684
 NPARAMETR:  1.0098E+00  1.3239E+00  1.1889E+00  8.6234E-01  1.1697E+00  1.0181E+00  1.1103E+00  1.1432E+00  1.0505E+00  1.1987E+00
             9.7041E-01
 PARAMETER:  1.0974E-01  3.8061E-01  2.7307E-01 -4.8101E-02  2.5676E-01  1.1799E-01  2.0464E-01  2.3385E-01  1.4924E-01  2.8126E-01
             6.9958E-02
 GRADIENT:   1.1113E+01  2.0133E+01  2.0787E+01 -4.3510E+00 -2.8208E+01  7.7917E-03 -2.9348E-03 -6.0394E+00 -1.6746E+00  8.1703E-01
            -9.8511E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1671.29980258787        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:      787
 NPARAMETR:  1.0023E+00  1.3125E+00  1.1800E+00  8.6286E-01  1.1731E+00  1.0104E+00  1.1072E+00  1.1907E+00  1.0445E+00  1.1896E+00
             9.7367E-01
 PARAMETER:  1.0231E-01  3.7195E-01  2.6548E-01 -4.7498E-02  2.5965E-01  1.1035E-01  2.0187E-01  2.7456E-01  1.4351E-01  2.7363E-01
             7.3313E-02
 GRADIENT:   4.4753E+01  3.3317E+01  1.4677E+01 -1.4222E+00 -1.7560E+01  7.7417E+00  6.7866E-01 -4.2881E+00 -4.7923E-01  6.8791E-01
            -7.1081E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1671.36129907358        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      963             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0026E+00  1.3119E+00  1.1794E+00  8.6404E-01  1.1733E+00  1.0114E+00  1.1086E+00  1.1947E+00  1.0455E+00  1.1887E+00
             9.7610E-01
 PARAMETER:  1.0258E-01  3.7147E-01  2.6497E-01 -4.6142E-02  2.5982E-01  1.1136E-01  2.0309E-01  2.7793E-01  1.4448E-01  2.7284E-01
             7.5807E-02
 GRADIENT:   4.5167E+01  3.3408E+01  1.3923E+01 -2.1176E-01 -1.6704E+01  8.1866E+00  8.8175E-01 -4.0637E+00 -1.2049E-01  6.8667E-01
            -5.9344E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1671.38892073615        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1149             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0029E+00  1.3117E+00  1.1792E+00  8.6498E-01  1.1731E+00  1.0232E+00  1.1092E+00  1.1946E+00  1.0465E+00  1.1880E+00
             9.7870E-01
 PARAMETER:  1.0291E-01  3.7129E-01  2.6484E-01 -4.5050E-02  2.5969E-01  1.2291E-01  2.0363E-01  2.7779E-01  1.4546E-01  2.7225E-01
             7.8474E-02
 GRADIENT:   4.5802E+01  3.3807E+01  1.3662E+01  7.3028E-01 -1.6478E+01  1.3783E+01  9.9798E-01 -3.9739E+00  1.1309E-01  6.7899E-01
            -4.7509E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1671.41699434555        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1310            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0029E+00  1.3115E+00  1.1782E+00  8.6505E-01  1.1731E+00  1.0179E+00  1.1093E+00  1.1971E+00  1.0466E+00  1.1879E+00
             9.7888E-01
 PARAMETER:  1.0293E-01  3.7114E-01  2.6396E-01 -4.4968E-02  2.5969E-01  1.1775E-01  2.0369E-01  2.7990E-01  1.4557E-01  2.7222E-01
             7.8658E-02
 GRADIENT:   4.5800E+01  3.3522E+01  1.3272E+01  9.0920E-01 -1.6110E+01  1.1266E+01  1.0118E+00 -3.8486E+00  1.9599E-01  7.4104E-01
            -4.6132E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1671.42744311208        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1497
 NPARAMETR:  1.0029E+00  1.3113E+00  1.1778E+00  8.6509E-01  1.1731E+00  1.0178E+00  1.1093E+00  1.1986E+00  1.0466E+00  1.1879E+00
             9.7899E-01
 PARAMETER:  1.0288E-01  3.7103E-01  2.6365E-01 -4.4919E-02  2.5969E-01  1.1768E-01  2.0371E-01  2.8117E-01  1.4558E-01  2.7216E-01
             7.8763E-02
 GRADIENT:  -3.7048E+00  9.0058E+00  1.2844E+01 -5.4136E+00 -1.7507E+01 -3.0787E-02 -5.6720E-01 -3.8277E+00 -8.0956E-01  4.8686E-01
            -4.6113E+00

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1671.42869034098        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1567
 NPARAMETR:  1.0029E+00  1.3113E+00  1.1778E+00  8.6510E-01  1.1731E+00  1.0178E+00  1.1093E+00  1.1990E+00  1.0466E+00  1.1879E+00
             9.7901E-01
 PARAMETER:  1.0288E-01  3.7103E-01  2.6365E-01 -4.4909E-02  2.5969E-01  1.1768E-01  2.0372E-01  2.8145E-01  1.4559E-01  2.7215E-01
             7.8784E-02
 GRADIENT:  -1.2971E+06 -7.1931E+05 -1.0122E+06 -3.9518E+01 -5.1389E+05 -1.1340E+06 -1.1355E+06  9.4816E+05 -9.1658E+05  4.9034E+05
             2.6685E+06

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1567
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7601E-03 -1.8489E-02 -4.0242E-02  9.6630E-03 -2.9103E-02
 SE:             2.9853E-02  2.1824E-02  1.2823E-02  2.2053E-02  2.2267E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5298E-01  3.9689E-01  1.6992E-03  6.6126E-01  1.9122E-01

 ETASHRINKSD(%)  1.0000E-10  2.6886E+01  5.7042E+01  2.6119E+01  2.5402E+01
 ETASHRINKVR(%)  1.0000E-10  4.6543E+01  8.1546E+01  4.5416E+01  4.4352E+01
 EBVSHRINKSD(%)  3.9716E-01  2.6569E+01  6.4424E+01  2.8403E+01  2.1385E+01
 EBVSHRINKVR(%)  7.9275E-01  4.6079E+01  8.7344E+01  4.8739E+01  3.8197E+01
 RELATIVEINF(%)  9.8680E+01  1.8050E+00  1.5275E+00  1.7231E+00  1.6214E+01
 EPSSHRINKSD(%)  4.4410E+01
 EPSSHRINKVR(%)  6.9097E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1671.4286903409825     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -936.27786377724431     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.99
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1671.429       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.31E+00  1.18E+00  8.65E-01  1.17E+00  1.02E+00  1.11E+00  1.20E+00  1.05E+00  1.19E+00  9.79E-01
 


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
+        6.27E+09
 
 TH 2
+       -1.49E+03  2.82E+08
 
 TH 3
+       -2.32E+03  2.70E+03  6.92E+08
 
 TH 4
+        9.56E+04  2.06E+04  3.17E+04  8.92E+09
 
 TH 5
+       -2.37E+03  2.62E+03 -3.92E+04  3.25E+04  1.44E+09
 
 TH 6
+       -1.20E+04 -2.54E+03 -3.98E+03  1.60E+00 -4.06E+03  4.65E+09
 
 TH 7
+       -8.49E+03 -1.79E+03 -2.82E+03 -1.25E+01 -2.88E+03  2.47E+09  2.61E+09
 
 TH 8
+        3.95E+05  8.38E+04  1.31E+05  2.29E+09  6.49E+08  3.66E+03  2.60E+03  5.86E+08
 
 TH 9
+        6.09E+04  1.29E+04  2.02E+04 -5.06E+09  2.06E+04 -8.11E+03 -5.72E+03  1.30E+09  2.87E+09
 
 TH10
+        2.67E+08  2.73E+03  4.28E+03  2.39E+09  4.33E+03  3.82E+03  2.71E+03 -3.94E+03  1.35E+09  6.38E+08
 
 TH11
+        2.04E+06  4.33E+05  6.78E+05 -1.35E+00  6.92E+05 -5.69E+09 -3.02E+09 -6.24E+05  1.38E+06 -6.52E+05  6.96E+09
 
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
 #CPUT: Total CPU Time in Seconds,       28.473
Stop Time:
Sat Sep 25 10:00:12 CDT 2021
