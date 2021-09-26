Sat Sep 25 01:04:56 CDT 2021
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
$DATA ../../../../data/int/SL2/dat26.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      997
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      897
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -567.636502666011        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8150E+01 -3.0853E+01  3.6979E+02  3.4838E+01  1.1343E+02  2.2299E+01 -1.1914E+02 -7.9592E+02 -1.9446E+02 -6.3272E+01
            -5.4233E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2794.47151076188        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0140E+00  1.3140E+00  8.4057E-01  8.1726E-01  1.1377E+00  8.8588E-01  1.1160E+00  1.0035E+00  9.1096E-01  1.0555E+00
             2.3741E+00
 PARAMETER:  1.1388E-01  3.7308E-01 -7.3678E-02 -1.0180E-01  2.2899E-01 -2.1169E-02  2.0976E-01  1.0354E-01  6.7452E-03  1.5401E-01
             9.6462E-01
 GRADIENT:   1.6419E+01 -1.2383E+01 -1.2396E+01  3.0730E+00  1.0990E+01 -2.0495E+01  1.8125E+01  5.4036E+00 -1.1518E+01 -1.9764E+01
            -1.6442E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2807.63589553962        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0204E+00  1.7186E+00  8.2448E-01  6.3500E-01  1.3600E+00  9.3204E-01  8.9023E-01  4.5907E-01  1.2902E+00  1.4054E+00
             2.3928E+00
 PARAMETER:  1.2024E-01  6.4150E-01 -9.3003E-02 -3.5413E-01  4.0746E-01  2.9625E-02 -1.6278E-02 -6.7856E-01  3.5482E-01  4.4034E-01
             9.7245E-01
 GRADIENT:   3.0575E+01  9.6671E+01  4.3421E+00  5.3498E+01 -3.5049E+01 -6.8800E-01  1.7925E+01  1.6738E-01  9.1259E+00  6.1867E+00
            -1.3066E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2821.51552168204        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0112E+00  1.9831E+00  7.2443E-01  4.1474E-01  1.6671E+00  9.1321E-01  7.0823E-01  1.0813E-01  1.6771E+00  1.4632E+00
             2.4891E+00
 PARAMETER:  1.1118E-01  7.8468E-01 -2.2237E-01 -7.8011E-01  6.1109E-01  9.2075E-03 -2.4498E-01 -2.1245E+00  6.1709E-01  4.8065E-01
             1.0119E+00
 GRADIENT:   4.0551E+00  3.0569E+01  3.3551E+00  1.5465E+01  9.7478E+00 -7.2030E+00  4.3031E-01  1.6605E-03  3.7703E+00 -9.4596E+00
            -1.5383E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2825.73051946242        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  1.0091E+00  2.2500E+00  4.2662E-01  2.4850E-01  1.8187E+00  9.3889E-01  6.6660E-01  2.8522E-02  2.0547E+00  1.6500E+00
             2.4957E+00
 PARAMETER:  1.0902E-01  9.1094E-01 -7.5186E-01 -1.2923E+00  6.9813E-01  3.6942E-02 -3.0556E-01 -3.4571E+00  8.2013E-01  6.0075E-01
             1.0146E+00
 GRADIENT:  -9.9393E+00  3.7522E+01 -1.3006E+00  9.8393E+00 -1.1995E+00  2.1705E+00 -2.7899E+00  1.1150E-03 -2.7673E+00  3.1703E+00
            -1.5032E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2827.06096613691        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      523
 NPARAMETR:  1.0119E+00  2.3621E+00  3.5355E-01  1.5592E-01  1.9401E+00  9.3410E-01  6.6152E-01  1.0000E-02  2.6572E+00  1.6859E+00
             2.4938E+00
 PARAMETER:  1.1180E-01  9.5954E-01 -9.3974E-01 -1.7584E+00  7.6274E-01  3.1833E-02 -3.1322E-01 -4.6790E+00  1.0773E+00  6.2228E-01
             1.0138E+00
 GRADIENT:  -2.0400E+00 -1.8642E+00 -1.5346E+00  6.7016E-01  2.6996E+00  5.8351E-01  1.9052E+00  0.0000E+00 -8.7973E-01 -5.8349E-01
            -7.3606E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2827.37599833559        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  1.0128E+00  2.3862E+00  4.2525E-01  1.4487E-01  1.9732E+00  9.3279E-01  6.5590E-01  1.0000E-02  2.8343E+00  1.6964E+00
             2.4956E+00
 PARAMETER:  1.1267E-01  9.6970E-01 -7.5509E-01 -1.8319E+00  7.7964E-01  3.0420E-02 -3.2174E-01 -4.5328E+00  1.1418E+00  6.2849E-01
             1.0145E+00
 GRADIENT:   7.5950E-01 -1.2671E+01 -7.6390E+00  1.1078E+01 -5.9720E-01  1.7419E-02 -6.9998E+00  0.0000E+00  1.4107E+01  3.7174E-01
            -6.6194E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2827.82423481450        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  1.0125E+00  2.4455E+00  4.0383E-01  9.8540E-02  2.0147E+00  9.3253E-01  6.4994E-01  1.0000E-02  3.5066E+00  1.7279E+00
             2.4907E+00
 PARAMETER:  1.1241E-01  9.9427E-01 -8.0676E-01 -2.2173E+00  8.0047E-01  3.0150E-02 -3.3088E-01 -5.3535E+00  1.3546E+00  6.4690E-01
             1.0126E+00
 GRADIENT:  -6.7419E-02 -6.6515E-01 -1.9251E-01  4.2049E-01 -2.2805E-01  6.7589E-02 -4.4180E-01  0.0000E+00  7.9324E-01  1.9278E-02
            -2.8515E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2827.83615240199        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1055
 NPARAMETR:  1.0127E+00  2.4583E+00  3.9824E-01  8.8666E-02  2.0263E+00  9.3156E-01  6.5057E-01  1.0000E-02  3.6679E+00  1.7332E+00
             2.4898E+00
 PARAMETER:  1.1265E-01  9.9949E-01 -8.2070E-01 -2.3229E+00  8.0622E-01  2.9105E-02 -3.2991E-01 -5.5715E+00  1.3996E+00  6.4995E-01
             1.0122E+00
 GRADIENT:   4.4950E-01  1.9431E+00  2.2088E-01 -1.1483E+00  1.2833E+00 -3.6770E-01  1.5346E+00  0.0000E+00 -2.6784E+00  1.3832E-01
             7.1800E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2827.87594061587        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1231
 NPARAMETR:  1.0124E+00  2.4825E+00  3.8771E-01  7.1764E-02  2.0349E+00  9.3299E-01  6.5008E-01  1.0000E-02  4.1527E+00  1.7345E+00
             2.4873E+00
 PARAMETER:  1.1234E-01  1.0093E+00 -8.4751E-01 -2.5344E+00  8.1043E-01  3.0640E-02 -3.3066E-01 -6.0447E+00  1.5238E+00  6.5073E-01
             1.0112E+00
 GRADIENT:  -1.4244E-01 -1.2890E+00 -4.0987E-01  1.1888E+00 -9.0718E-01  2.2489E-01 -1.1977E+00  0.0000E+00  1.9291E+00 -6.8824E-02
            -2.9118E-01

0ITERATION NO.:   49    OBJECTIVE VALUE:  -2827.87809621871        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1358
 NPARAMETR:  1.0125E+00  2.4844E+00  3.8610E-01  6.9523E-02  2.0374E+00  9.3251E-01  6.4915E-01  1.0000E-02  4.2228E+00  1.7349E+00
             2.4869E+00
 PARAMETER:  1.1243E-01  1.0100E+00 -8.5165E-01 -2.5661E+00  8.1167E-01  3.0123E-02 -3.3208E-01 -6.1184E+00  1.5405E+00  6.5093E-01
             1.0110E+00
 GRADIENT:   7.7924E-03 -3.1973E-02 -1.0950E-03  1.3403E-03  7.1688E-04 -6.5272E-04 -8.6540E-03  0.0000E+00  2.0905E-03  9.7241E-04
             4.4056E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1358
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1693E-03 -2.7290E-02 -9.4576E-06  3.8485E-02 -2.0157E-02
 SE:             2.9446E-02  2.5776E-02  5.6147E-06  1.5277E-02  2.6711E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6832E-01  2.8972E-01  9.2100E-02  1.1765E-02  4.5047E-01

 ETASHRINKSD(%)  1.3527E+00  1.3646E+01  9.9981E+01  4.8820E+01  1.0515E+01
 ETASHRINKVR(%)  2.6870E+00  2.5430E+01  1.0000E+02  7.3806E+01  1.9924E+01
 EBVSHRINKSD(%)  1.5579E+00  9.4280E+00  9.9921E+01  6.3481E+01  7.2260E+00
 EBVSHRINKVR(%)  3.0915E+00  1.7967E+01  1.0000E+02  8.6663E+01  1.3930E+01
 RELATIVEINF(%)  9.6844E+01  3.4967E+01  5.3545E-05  5.5632E+00  7.6407E+01
 EPSSHRINKSD(%)  1.6314E+01
 EPSSHRINKVR(%)  2.9967E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          897
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1648.5757285691827     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2827.8780962187052     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1179.3023676495225     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    33.85
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2827.878       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  2.48E+00  3.86E-01  6.95E-02  2.04E+00  9.33E-01  6.49E-01  1.00E-02  4.22E+00  1.73E+00  2.49E+00
 


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
+       -3.26E+01  1.34E+04
 
 TH 3
+       -3.49E+01  9.64E+04  4.56E+02
 
 TH 4
+        1.41E+03 -1.80E+05 -2.29E+04  2.51E+06
 
 TH 5
+       -1.09E+01  1.92E+04  9.69E+01 -2.70E+05  1.01E+02
 
 TH 6
+        7.17E+00 -1.34E+01 -1.49E+01  7.39E+02 -4.99E+00  2.18E+02
 
 TH 7
+       -7.90E+01  1.48E+05  1.02E+03 -2.07E+06  2.55E+02 -4.64E+01  3.05E+03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.19E+01 -5.01E+03 -3.74E+04  6.97E+04 -7.46E+03  4.79E+01 -5.73E+04  0.00E+00  1.93E+03
 
 TH10
+       -7.85E-01 -2.05E+01  1.83E+01 -6.41E+02 -2.23E+00 -1.74E-01  3.75E+01  0.00E+00 -4.19E+01  4.54E+01
 
 TH11
+       -1.97E+01  1.26E+04  5.29E+01 -2.87E+03  1.30E+01  1.75E+00  1.32E+02  0.00E+00 -4.91E+03  6.62E+00  1.97E+02
 
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
 #CPUT: Total CPU Time in Seconds,       48.329
Stop Time:
Sat Sep 25 01:05:46 CDT 2021
