Sat Sep 18 13:39:37 CDT 2021
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
$DATA ../../../../data/spa/S2/dat80.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1728.36408606716        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.2954E-01 -1.0852E+01  8.7376E+00 -1.5935E+01 -4.3769E+01  4.1502E+01  6.2091E+00  3.7855E+00  3.8022E+01  1.4822E+01
            -5.4234E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1735.58606112484        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  1.0083E+00  1.0650E+00  1.0275E+00  9.7788E-01  1.0735E+00  8.5130E-01  9.7553E-01  9.8389E-01  8.1956E-01  9.3866E-01
             1.0120E+00
 PARAMETER:  1.0825E-01  1.6297E-01  1.2708E-01  7.7636E-02  1.7094E-01 -6.0992E-02  7.5228E-02  8.3759E-02 -9.8989E-02  3.6699E-02
             1.1196E-01
 GRADIENT:  -3.7608E+01  1.8411E-01  3.4818E+00  7.3458E-01  4.5806E+00 -2.2778E+01 -4.1036E+00 -2.9252E+00 -1.5663E+00 -2.3296E+00
            -2.8147E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1736.15551924754        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      290
 NPARAMETR:  1.0202E+00  1.0886E+00  9.7461E-01  9.7365E-01  1.0404E+00  8.7798E-01  1.0257E+00  1.0235E+00  7.8376E-01  8.7422E-01
             1.0233E+00
 PARAMETER:  1.2002E-01  1.8486E-01  7.4283E-02  7.3300E-02  1.3960E-01 -3.0135E-02  1.2534E-01  1.2324E-01 -1.4366E-01 -3.4420E-02
             1.2299E-01
 GRADIENT:  -3.9572E+00  2.2896E+01  5.9759E+00  2.1334E+01 -8.6494E+00 -9.4956E+00 -6.9926E-01 -8.1591E-01 -5.6666E+00 -3.0077E+00
             1.7717E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1736.70125071635        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      468
 NPARAMETR:  1.0226E+00  1.1463E+00  8.6770E-01  9.1967E-01  1.0273E+00  9.0066E-01  9.8122E-01  8.8865E-01  8.3526E-01  8.8891E-01
             1.0135E+00
 PARAMETER:  1.2233E-01  2.3655E-01 -4.1904E-02  1.6264E-02  1.2696E-01 -4.6228E-03  8.1043E-02 -1.8051E-02 -8.0008E-02 -1.7762E-02
             1.1343E-01
 GRADIENT:   3.6559E-01  4.9363E-01 -1.2139E-01  1.2107E+00 -9.7076E-01  2.2825E-01  9.3531E-02  1.1187E-01  8.5982E-02  4.6365E-01
             1.0839E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1736.76057593246        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      643
 NPARAMETR:  1.0241E+00  1.3232E+00  6.8652E-01  8.0010E-01  1.0204E+00  9.0168E-01  8.9575E-01  7.0232E-01  9.0084E-01  8.5965E-01
             1.0131E+00
 PARAMETER:  1.2382E-01  3.8003E-01 -2.7612E-01 -1.2302E-01  1.2022E-01 -3.4992E-03 -1.0089E-02 -2.5337E-01 -4.4303E-03 -5.1234E-02
             1.1306E-01
 GRADIENT:  -2.1978E-01  2.7129E+00  3.9356E-01  1.5328E+00 -2.1729E+00 -1.4784E-01  1.7682E-01  1.6527E-01  1.3300E-02  3.3012E-01
             2.8915E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1736.78370505286        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      819
 NPARAMETR:  1.0247E+00  1.4211E+00  5.9010E-01  7.3097E-01  1.0218E+00  9.0231E-01  8.5354E-01  5.3898E-01  9.4787E-01  8.4678E-01
             1.0132E+00
 PARAMETER:  1.2438E-01  4.5141E-01 -4.2747E-01 -2.1339E-01  1.2155E-01 -2.7998E-03 -5.8365E-02 -5.1807E-01  4.6464E-02 -6.6309E-02
             1.1310E-01
 GRADIENT:  -3.9615E-01 -1.0036E+00 -4.7419E-01 -4.6710E-01  5.3758E-01 -2.1386E-01  2.3561E-01  1.5145E-01  4.7171E-02  1.0079E-01
             6.0490E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1736.79917157173        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      994
 NPARAMETR:  1.0252E+00  1.4891E+00  5.2276E-01  6.8335E-01  1.0203E+00  9.0329E-01  8.2662E-01  3.3407E-01  9.8438E-01  8.3670E-01
             1.0132E+00
 PARAMETER:  1.2493E-01  4.9817E-01 -5.4863E-01 -2.8075E-01  1.2012E-01 -1.7145E-03 -9.0407E-02 -9.9642E-01  8.4259E-02 -7.8290E-02
             1.1309E-01
 GRADIENT:   2.4278E-01 -1.0061E+00 -4.1694E-01 -6.8879E-01  6.2740E-01  2.0873E-02  2.2197E-01  7.6807E-02  7.1541E-03  4.8274E-02
             6.3633E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1736.82639477487        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1170
 NPARAMETR:  1.0254E+00  1.4602E+00  5.1628E-01  7.0027E-01  9.9446E-01  9.0381E-01  8.3973E-01  1.2895E-01  9.6526E-01  8.2057E-01
             1.0125E+00
 PARAMETER:  1.2508E-01  4.7857E-01 -5.6112E-01 -2.5628E-01  9.4447E-02 -1.1336E-03 -7.4672E-02 -1.9483E+00  6.4643E-02 -9.7756E-02
             1.1242E-01
 GRADIENT:   4.8071E-01  4.6630E-01 -5.5609E-02  3.4519E-01 -2.9046E-01  1.9207E-01  4.1982E-02  7.5036E-03 -3.8470E-02  8.4981E-02
             5.3839E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1736.83028831061        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1346
 NPARAMETR:  1.0252E+00  1.4721E+00  5.0894E-01  6.9213E-01  9.9754E-01  9.0338E-01  8.3448E-01  2.2426E-02  9.7297E-01  8.2124E-01
             1.0125E+00
 PARAMETER:  1.2491E-01  4.8671E-01 -5.7542E-01 -2.6798E-01  9.7533E-02 -1.6109E-03 -8.0950E-02 -3.6975E+00  7.2599E-02 -9.6945E-02
             1.1240E-01
 GRADIENT:  -2.1226E-02  3.7306E-02 -4.2269E-02  4.9697E-02 -5.9726E-02 -5.9865E-03 -1.8939E-04  2.3521E-04  1.6449E-02  2.5248E-02
             7.8738E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1736.83040038597        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1521
 NPARAMETR:  1.0252E+00  1.4719E+00  5.0943E-01  6.9229E-01  9.9785E-01  9.0338E-01  8.3463E-01  1.0000E-02  9.7292E-01  8.2157E-01
             1.0125E+00
 PARAMETER:  1.2491E-01  4.8656E-01 -5.7446E-01 -2.6775E-01  9.7843E-02 -1.6068E-03 -8.0770E-02 -4.6266E+00  7.2547E-02 -9.6541E-02
             1.1242E-01
 GRADIENT:  -7.0208E-03 -9.7622E-04 -3.8780E-03  1.1155E-02  1.7805E-02 -3.2435E-03 -7.0432E-04  0.0000E+00 -4.2555E-04 -1.6666E-03
            -6.8740E-04

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1736.83040038597        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1549
 NPARAMETR:  1.0253E+00  1.4718E+00  5.0947E-01  6.9222E-01  9.9775E-01  9.0348E-01  8.3470E-01  1.0000E-02  9.7295E-01  8.2163E-01
             1.0125E+00
 PARAMETER:  1.2491E-01  4.8656E-01 -5.7446E-01 -2.6775E-01  9.7843E-02 -1.6068E-03 -8.0770E-02 -4.6266E+00  7.2547E-02 -9.6541E-02
             1.1242E-01
 GRADIENT:  -3.3559E-02  2.1584E-02 -3.6556E-03  1.4338E-02  1.8349E-02 -4.3930E-03 -1.6902E-03  0.0000E+00 -5.0166E-04 -1.1559E-03
            -2.2542E-04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1549
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.0284E-05 -1.6602E-02 -3.3751E-04  1.3289E-02 -2.5398E-02
 SE:             2.9826E-02  2.4328E-02  1.4194E-04  2.2731E-02  2.1513E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9839E-01  4.9499E-01  1.7414E-02  5.5881E-01  2.3777E-01

 ETASHRINKSD(%)  7.8341E-02  1.8497E+01  9.9524E+01  2.3848E+01  2.7929E+01
 ETASHRINKVR(%)  1.5662E-01  3.3572E+01  9.9998E+01  4.2009E+01  4.8058E+01
 EBVSHRINKSD(%)  5.1333E-01  1.8407E+01  9.9581E+01  2.4923E+01  2.7156E+01
 EBVSHRINKVR(%)  1.0240E+00  3.3426E+01  9.9998E+01  4.3635E+01  4.6938E+01
 RELATIVEINF(%)  9.8916E+01  3.2262E+00  1.1893E-04  2.4506E+00  6.0573E+00
 EPSSHRINKSD(%)  4.3660E+01
 EPSSHRINKVR(%)  6.8258E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1736.8304003859703     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1001.6795738222321     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.87
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.72
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1736.830       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.47E+00  5.09E-01  6.92E-01  9.98E-01  9.03E-01  8.35E-01  1.00E-02  9.73E-01  8.22E-01  1.01E+00
 


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
+        1.28E+03
 
 TH 2
+       -9.09E+00  5.04E+02
 
 TH 3
+        1.29E+01  2.48E+02  8.32E+02
 
 TH 4
+       -2.33E+01  4.07E+02 -6.03E+02  1.40E+03
 
 TH 5
+       -4.63E+00 -3.05E+02 -7.10E+02  4.90E+02  8.80E+02
 
 TH 6
+       -2.86E-01 -1.60E+00  2.85E+00 -5.78E+00 -9.85E-01  2.41E+02
 
 TH 7
+        4.15E-01  1.93E+01 -3.59E+01 -1.29E+01  5.25E-01  7.03E-01  1.35E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.67E-01 -2.08E+01 -4.72E+01  6.17E+01 -2.19E+00 -2.15E+00  2.22E+01  0.00E+00  8.02E+01
 
 TH10
+       -5.33E-01 -1.51E+01 -4.91E+01 -1.36E+01 -7.50E+01 -6.24E-01  1.91E+01  0.00E+00  1.35E+01  9.37E+01
 
 TH11
+       -1.00E+01 -1.80E+01 -3.49E+01  4.25E+00 -4.73E+00  2.06E+00  1.10E+01  0.00E+00  1.29E+01  2.03E+01  2.09E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.656
Stop Time:
Sat Sep 18 13:40:04 CDT 2021
