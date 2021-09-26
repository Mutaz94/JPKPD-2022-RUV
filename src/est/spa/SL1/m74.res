Sat Sep 25 10:41:51 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat74.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1671.61454849599        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0653E+01 -4.8665E+01 -1.2186E+01 -6.2645E+01  1.4163E+01  2.3409E+01 -6.1233E+00  8.1295E+00 -9.4933E+00  1.1500E+01
             1.2347E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1674.98371850942        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0022E+00  1.0180E+00  1.0050E+00  1.0349E+00  9.7683E-01  9.6583E-01  1.0790E+00  8.8617E-01  1.0520E+00  8.8376E-01
             1.0177E+00
 PARAMETER:  1.0223E-01  1.1779E-01  1.0501E-01  1.3426E-01  7.6562E-02  6.5235E-02  1.7602E-01 -2.0844E-02  1.5073E-01 -2.3572E-02
             1.1750E-01
 GRADIENT:   5.5126E+01  1.1440E+01  1.5040E+01  5.4343E+00 -9.5876E+00  1.0304E+01  1.0657E+00  1.3344E+00  7.1672E+00 -3.2431E+00
             4.6086E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1676.70576824875        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.9790E-01  9.0806E-01  7.0908E-01  1.0796E+00  7.6918E-01  9.5392E-01  1.2057E+00  4.5611E-01  9.6301E-01  7.1091E-01
             1.0235E+00
 PARAMETER:  9.7893E-02  3.5504E-03 -2.4379E-01  1.7663E-01 -1.6243E-01  5.2826E-02  2.8707E-01 -6.8502E-01  6.2305E-02 -2.4121E-01
             1.2322E-01
 GRADIENT:   3.6728E+01  5.9844E+00 -6.5890E+00  2.1865E+01 -2.2952E+00  4.6566E+00 -5.9872E-01  1.5469E+00  4.8542E+00  2.0246E+00
             1.0535E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1677.34569067479        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      276
 NPARAMETR:  9.9603E-01  8.6473E-01  7.9382E-01  1.1116E+00  8.0414E-01  9.4949E-01  1.2662E+00  4.3899E-01  9.4059E-01  7.7674E-01
             1.0058E+00
 PARAMETER:  9.6020E-02 -4.5343E-02 -1.3089E-01  2.0578E-01 -1.1798E-01  4.8175E-02  3.3601E-01 -7.2328E-01  3.8748E-02 -1.5264E-01
             1.0578E-01
 GRADIENT:  -3.1028E+00  4.6575E+00  2.9399E+00  1.1808E+00 -5.3177E+00 -2.2172E-01 -9.5582E-02 -1.1512E-01 -3.7075E-02  5.2657E-01
             1.4584E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1677.86602415367        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      453
 NPARAMETR:  9.9352E-01  5.8342E-01  8.9362E-01  1.2878E+00  7.5403E-01  9.4748E-01  1.6720E+00  5.7171E-01  8.5498E-01  7.7711E-01
             1.0001E+00
 PARAMETER:  9.3494E-02 -4.3884E-01 -1.2477E-02  3.5297E-01 -1.8232E-01  4.6050E-02  6.1402E-01 -4.5912E-01 -5.6680E-02 -1.5217E-01
             1.0012E-01
 GRADIENT:  -4.1483E-01  5.4750E+00  7.6289E+00  6.4175E+00 -1.0851E+01  2.7428E-01 -2.0964E-01 -1.5763E-01 -7.5865E-01  1.6903E-01
            -8.7275E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1678.16364964818        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      635
 NPARAMETR:  9.9174E-01  4.1791E-01  8.9537E-01  1.3788E+00  7.1009E-01  9.4457E-01  2.0875E+00  5.7108E-01  8.2012E-01  7.6771E-01
             1.0019E+00
 PARAMETER:  9.1702E-02 -7.7249E-01 -1.0517E-02  4.2124E-01 -2.4236E-01  4.2979E-02  8.3597E-01 -4.6023E-01 -9.8301E-02 -1.6434E-01
             1.0187E-01
 GRADIENT:   1.6774E+00  1.7091E+00  5.9750E-01  5.5474E+00 -1.6735E+00  7.9713E-03  9.3883E-02 -2.0272E-01 -8.1247E-02 -5.2355E-01
            -9.8941E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1678.28990704700        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      810
 NPARAMETR:  9.8792E-01  2.8469E-01  9.3602E-01  1.4559E+00  6.9749E-01  9.4206E-01  2.6109E+00  6.2638E-01  7.9519E-01  7.8526E-01
             1.0017E+00
 PARAMETER:  8.7851E-02 -1.1564E+00  3.3884E-02  4.7559E-01 -2.6027E-01  4.0312E-02  1.0597E+00 -3.6780E-01 -1.2918E-01 -1.4174E-01
             1.0172E-01
 GRADIENT:  -8.7025E-02  7.7374E-01  1.6890E+00  2.6336E+00 -2.3532E+00  2.4074E-02  1.7107E-01 -2.5603E-01 -1.4184E-01 -3.2473E-02
            -1.4474E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1678.39559657833        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      987
 NPARAMETR:  9.8399E-01  1.4332E-01  1.0082E+00  1.5431E+00  6.9859E-01  9.3907E-01  3.7211E+00  7.5444E-01  7.7140E-01  7.9328E-01
             1.0009E+00
 PARAMETER:  8.3864E-02 -1.8427E+00  1.0813E-01  5.3381E-01 -2.5869E-01  3.7131E-02  1.4140E+00 -1.8177E-01 -1.5955E-01 -1.3158E-01
             1.0090E-01
 GRADIENT:  -6.6313E-01  7.3478E-01  2.2288E+00  5.0173E+00 -3.6580E+00  6.8149E-02  1.5149E-01  2.8756E-01 -5.2388E-01  2.3425E-01
            -1.3096E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1678.67045516902        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  9.8346E-01  6.5426E-02  9.6786E-01  1.5790E+00  6.6453E-01  9.3662E-01  5.3225E+00  7.2088E-01  7.6905E-01  7.7206E-01
             1.0024E+00
 PARAMETER:  8.3324E-02 -2.6268E+00  6.7334E-02  5.5681E-01 -3.0868E-01  3.4518E-02  1.7719E+00 -2.2728E-01 -1.6261E-01 -1.5870E-01
             1.0241E-01
 GRADIENT:   2.4718E+00 -4.0257E-01  6.1275E-01  1.3321E+01 -3.6862E+00 -4.9298E-01 -2.8473E+00  8.5524E-02  3.2268E+00  6.9246E-01
             5.6852E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1679.77995770948        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1343
 NPARAMETR:  9.7996E-01  1.0000E-02  9.2400E-01  1.5846E+00  6.3560E-01  9.3730E-01  1.2428E+01  6.8739E-01  7.4830E-01  7.3254E-01
             1.0008E+00
 PARAMETER:  7.9756E-02 -4.6007E+00  2.0960E-02  5.6031E-01 -3.5319E-01  3.5251E-02  2.6200E+00 -2.7486E-01 -1.8995E-01 -2.1124E-01
             1.0078E-01
 GRADIENT:  -2.9339E+00  0.0000E+00  1.2327E+00 -5.7408E+00 -3.3971E-01  3.4347E-01 -2.9296E+00 -2.8168E-02 -5.9481E-01  2.8375E-01
            -6.2001E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1679.81518734644        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1520
 NPARAMETR:  9.8108E-01  1.0000E-02  9.2865E-01  1.5891E+00  6.3738E-01  9.3674E-01  1.2590E+01  6.9832E-01  7.5276E-01  7.3270E-01
             1.0022E+00
 PARAMETER:  8.0897E-02 -4.6140E+00  2.5974E-02  5.6317E-01 -3.5039E-01  3.4651E-02  2.6329E+00 -2.5908E-01 -1.8401E-01 -2.1102E-01
             1.0221E-01
 GRADIENT:  -2.9828E-02  0.0000E+00 -4.7495E-02 -1.4685E-01 -1.3712E-01 -2.6492E-03  1.7103E-01  1.6592E-02 -5.1560E-02  7.5282E-02
             3.1902E-02

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1679.81529175455        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1647
 NPARAMETR:  9.8109E-01  1.0000E-02  9.3003E-01  1.5893E+00  6.3807E-01  9.3671E-01  1.2591E+01  7.0002E-01  7.5259E-01  7.3221E-01
             1.0022E+00
 PARAMETER:  8.0913E-02 -4.6140E+00  2.7464E-02  5.6331E-01 -3.4931E-01  3.4623E-02  2.6330E+00 -2.5664E-01 -1.8423E-01 -2.1169E-01
             1.0220E-01
 GRADIENT:  -8.1276E-03  0.0000E+00 -2.2228E-03  3.2361E-02  5.8695E-03 -5.6286E-03 -2.0244E-02  1.1128E-04  1.0496E-02  1.3764E-03
             3.6153E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1647
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0891E-04  1.4392E-02 -1.6609E-02 -7.3144E-03 -1.7978E-02
 SE:             2.9819E-02  6.9775E-03  1.5709E-02  2.8855E-02  2.1670E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9441E-01  3.9150E-02  2.9038E-01  7.9989E-01  4.0674E-01

 ETASHRINKSD(%)  1.0421E-01  7.6625E+01  4.7372E+01  3.3328E+00  2.7402E+01
 ETASHRINKVR(%)  2.0832E-01  9.4536E+01  7.2303E+01  6.5545E+00  4.7296E+01
 EBVSHRINKSD(%)  4.6557E-01  8.1193E+01  4.8387E+01  3.4193E+00  2.5559E+01
 EBVSHRINKVR(%)  9.2896E-01  9.6463E+01  7.3361E+01  6.7216E+00  4.4586E+01
 RELATIVEINF(%)  9.9003E+01  2.9788E+00  2.6244E+00  6.7343E+01  5.4622E+00
 EPSSHRINKSD(%)  4.4223E+01
 EPSSHRINKVR(%)  6.8889E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1679.8152917545522     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -944.66446519081398     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.90
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.76
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1679.815       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.00E-02  9.30E-01  1.59E+00  6.38E-01  9.37E-01  1.26E+01  7.00E-01  7.53E-01  7.32E-01  1.00E+00
 


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
+        0.00E+00  0.00E+00
 
 TH 3
+       -8.76E+01  0.00E+00  1.87E+03
 
 TH 4
+        3.61E+01  0.00E+00 -7.69E+02  3.25E+05
 
 TH 5
+       -2.13E+02  0.00E+00 -1.26E+07  1.57E+03  5.25E+06
 
 TH 6
+       -3.93E+00  0.00E+00  8.34E+01 -6.40E+01  2.80E+02  2.26E+02
 
 TH 7
+       -1.08E+00  0.00E+00  1.57E+01  2.58E+01 -3.85E+01  1.44E+00  2.37E+02
 
 TH 8
+       -2.58E+01  0.00E+00  2.39E+02 -1.40E+02  6.62E+02  2.20E+01  2.93E+00  1.11E+02
 
 TH 9
+       -1.41E+02  0.00E+00  2.31E+03 -8.86E+02 -8.44E+06  1.55E+02  1.92E+01  5.55E+02  1.36E+07
 
 TH10
+       -6.71E+01  0.00E+00  5.92E+02 -2.27E+02 -7.55E+06  5.23E+01  4.17E+00  2.06E+02  1.26E+03  4.50E+02
 
 TH11
+       -3.60E+01  0.00E+00  6.92E+01 -4.71E+01  1.71E+02  2.74E+01  7.30E-01  5.19E+01  1.77E+02  5.78E+01  1.94E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.728
Stop Time:
Sat Sep 25 10:42:21 CDT 2021
