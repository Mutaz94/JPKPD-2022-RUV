Sat Sep 25 13:02:52 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1708.85724646088        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9358E+00  1.4153E+01 -2.4061E+01  6.5214E+01  4.1560E+01  5.6796E+00  1.9050E+01  3.8296E+00  4.0725E+01  1.3361E+00
             2.0105E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1716.54630759645        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0022E+00  9.7799E-01  9.7690E-01  9.8837E-01  9.5858E-01  9.9143E-01  8.6831E-01  9.8834E-01  8.1137E-01  9.8786E-01
             9.9367E-01
 PARAMETER:  1.0218E-01  7.7741E-02  7.6629E-02  8.8299E-02  5.7700E-02  9.1396E-02 -4.1202E-02  8.8269E-02 -1.0904E-01  8.7789E-02
             9.3654E-02
 GRADIENT:   9.3293E+00  9.3866E+00 -8.7047E+00  2.1430E+01  1.1068E+01  3.0069E+00  2.9368E+00  3.2174E+00  1.1125E+00  4.4490E+00
             1.4969E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1717.30442758427        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      270
 NPARAMETR:  1.0171E+00  9.2187E-01  9.8185E-01  1.0275E+00  9.2807E-01  9.9880E-01  8.2262E-01  9.2023E-01  8.2271E-01  9.5409E-01
             9.8018E-01
 PARAMETER:  1.1696E-01  1.8652E-02  8.1682E-02  1.2715E-01  2.5355E-02  9.8798E-02 -9.5265E-02  1.6873E-02 -9.5153E-02  5.3006E-02
             7.9986E-02
 GRADIENT:  -3.0781E+00  1.3962E+01 -1.8091E+00  2.6006E+01  1.6656E+00  1.6379E+00  8.4357E-02  8.7255E-01  3.6524E+00 -3.3686E-01
             8.5758E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1717.83041625965        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  1.0190E+00  8.4586E-01  8.6875E-01  1.0577E+00  8.3958E-01  9.9530E-01  9.7022E-01  6.9556E-01  7.5460E-01  8.8840E-01
             9.5950E-01
 PARAMETER:  1.1879E-01 -6.7403E-02 -4.0697E-02  1.5610E-01 -7.4855E-02  9.5287E-02  6.9768E-02 -2.6304E-01 -1.8157E-01 -1.8334E-02
             5.8655E-02
 GRADIENT:  -1.2451E-01  5.6241E+00  1.8267E+00  6.7438E+00 -3.1498E+00 -1.0238E-01  8.2317E-03 -1.6372E-01 -8.8657E-01 -5.0750E-01
            -1.1419E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1717.91561836225        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      623
 NPARAMETR:  1.0166E+00  7.2213E-01  1.0355E+00  1.1453E+00  8.6869E-01  9.9428E-01  1.0121E+00  8.8247E-01  7.2825E-01  9.4113E-01
             9.6387E-01
 PARAMETER:  1.1649E-01 -2.2555E-01  1.3485E-01  2.3570E-01 -4.0765E-02  9.4265E-02  1.1203E-01 -2.5026E-02 -2.1711E-01  3.9326E-02
             6.3201E-02
 GRADIENT:   3.2072E-01  8.2164E+00  3.1850E+00  1.2224E+01 -9.0022E+00  5.5820E-01 -1.5085E-01  2.4520E-01 -3.9316E-01  1.9911E+00
             1.2252E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1718.15726614647        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  1.0140E+00  5.8645E-01  1.1527E+00  1.2340E+00  8.7216E-01  9.8988E-01  1.1349E+00  9.8346E-01  6.9025E-01  9.4106E-01
             9.6031E-01
 PARAMETER:  1.1391E-01 -4.3366E-01  2.4215E-01  3.1026E-01 -3.6779E-02  8.9826E-02  2.2657E-01  8.3318E-02 -2.7070E-01  3.9256E-02
             5.9503E-02
 GRADIENT:  -1.9836E-01  7.0928E+00  3.2572E+00  1.5814E+01 -5.5966E+00 -2.5424E-01 -2.7649E-01 -5.6477E-01 -7.9355E-01 -7.6351E-01
            -8.5229E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1718.82532205967        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      977
 NPARAMETR:  1.0092E+00  3.2274E-01  1.4081E+00  1.4051E+00  8.8494E-01  9.8534E-01  1.6408E+00  1.2534E+00  6.2505E-01  9.7738E-01
             9.6152E-01
 PARAMETER:  1.0916E-01 -1.0309E+00  4.4225E-01  4.4008E-01 -2.2241E-02  8.5233E-02  5.9516E-01  3.2588E-01 -3.6992E-01  7.7120E-02
             6.0759E-02
 GRADIENT:   5.3774E-02  4.4053E+00  5.0030E-01  1.8164E+01 -6.0351E+00  4.8421E-02  1.8217E-01  1.1725E+00 -9.3153E-01  1.6451E+00
             2.3811E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1719.65520616881        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1152
 NPARAMETR:  1.0051E+00  1.2256E-01  1.6641E+00  1.5374E+00  9.0806E-01  9.8017E-01  2.7517E+00  1.5227E+00  5.8748E-01  9.9357E-01
             9.5964E-01
 PARAMETER:  1.0512E-01 -1.9991E+00  6.0927E-01  5.3012E-01  3.5537E-03  7.9966E-02  1.1122E+00  5.2046E-01 -4.3191E-01  9.3551E-02
             5.8808E-02
 GRADIENT:  -1.0363E+00  1.7011E+00 -2.6467E+00  2.0846E+01 -2.1115E+00 -3.4102E-01  2.7002E-02  2.3122E+00 -1.7504E+00  8.8002E-01
            -1.3328E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1720.02970301873        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1327
 NPARAMETR:  1.0030E+00  4.4988E-02  1.8459E+00  1.5900E+00  9.2962E-01  9.7756E-01  4.3165E+00  1.6558E+00  5.7829E-01  9.9418E-01
             9.6073E-01
 PARAMETER:  1.0298E-01 -3.0014E+00  7.1295E-01  5.6372E-01  2.7017E-02  7.7308E-02  1.5624E+00  6.0429E-01 -4.4768E-01  9.4161E-02
             5.9941E-02
 GRADIENT:  -3.2815E+00  5.0067E-01  1.3635E+00  1.4382E+01 -2.9528E+00 -7.3897E-01 -4.0144E-02  1.7105E-01 -2.7451E-01 -1.6001E+00
             1.3466E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1720.19945109258        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1502
 NPARAMETR:  1.0038E+00  1.1682E-02  1.7973E+00  1.6073E+00  9.1346E-01  9.7924E-01  7.5961E+00  1.6075E+00  5.7397E-01  1.0073E+00
             9.5984E-01
 PARAMETER:  1.0379E-01 -4.3497E+00  6.8628E-01  5.7455E-01  9.4853E-03  7.9020E-02  2.1276E+00  5.7471E-01 -4.5518E-01  1.0727E-01
             5.9008E-02
 GRADIENT:   2.6544E-01  1.0112E-01 -3.1733E-01  1.2419E+01 -2.6001E+00  1.9778E-01 -2.2345E-02  1.4047E-01 -6.9258E-01  1.1240E+00
             1.5270E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1720.23076323103        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1679
 NPARAMETR:  1.0035E+00  1.0000E-02  1.8512E+00  1.6075E+00  9.2623E-01  9.7867E-01  8.7826E+00  1.6524E+00  5.7408E-01  1.0071E+00
             9.5938E-01
 PARAMETER:  1.0354E-01 -4.6895E+00  7.1585E-01  5.7467E-01  2.3370E-02  7.8443E-02  2.2728E+00  6.0223E-01 -4.5499E-01  1.0708E-01
             5.8534E-02
 GRADIENT:  -2.4881E-01  0.0000E+00  6.1504E-02 -8.2826E-01  2.1980E-01  4.1742E-03 -1.3860E-02 -4.0700E-02  7.7105E-03 -7.6668E-02
            -3.0016E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1720.23178300582        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1841
 NPARAMETR:  1.0037E+00  1.0000E-02  1.8438E+00  1.6073E+00  9.2445E-01  9.7872E-01  9.5895E+00  1.6469E+00  5.7398E-01  1.0063E+00
             9.5945E-01
 PARAMETER:  1.0369E-01 -4.8883E+00  7.1185E-01  5.7459E-01  2.1439E-02  7.8494E-02  2.3607E+00  5.9887E-01 -4.5517E-01  1.0628E-01
             5.8603E-02
 GRADIENT:   4.6263E-02  0.0000E+00  1.2208E-02  8.5691E-03 -2.3296E-02  2.1590E-02 -2.8973E-05 -3.3790E-04 -1.1226E-05  3.3490E-03
             7.5161E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1841
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.1262E-04  6.4901E-05 -3.5511E-02 -9.2557E-03 -4.3507E-02
 SE:             2.9868E-02  1.9864E-03  1.8897E-02  2.8842E-02  2.0330E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9699E-01  9.7394E-01  6.0218E-02  7.4828E-01  3.2354E-02

 ETASHRINKSD(%)  1.0000E-10  9.3345E+01  3.6692E+01  3.3771E+00  3.1890E+01
 ETASHRINKVR(%)  1.0000E-10  9.9557E+01  5.9921E+01  6.6401E+00  5.3611E+01
 EBVSHRINKSD(%)  3.8599E-01  9.3534E+01  4.0125E+01  3.8140E+00  2.8098E+01
 EBVSHRINKVR(%)  7.7049E-01  9.9582E+01  6.4150E+01  7.4825E+00  4.8301E+01
 RELATIVEINF(%)  9.5559E+01  9.2482E-03  8.6922E+00  2.2759E+00  8.4794E+00
 EPSSHRINKSD(%)  4.5163E+01
 EPSSHRINKVR(%)  6.9929E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1720.2317830058162     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -985.08095644207799     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.29
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1720.232       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  1.84E+00  1.61E+00  9.24E-01  9.79E-01  9.59E+00  1.65E+00  5.74E-01  1.01E+00  9.59E-01
 


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
+        1.15E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.46E+00  0.00E+00  4.65E+01
 
 TH 4
+       -1.52E+01  0.00E+00 -2.83E+01  1.21E+03
 
 TH 5
+       -1.31E+01  0.00E+00 -1.33E+02 -1.03E+02  6.86E+02
 
 TH 6
+       -1.82E+01  0.00E+00  2.86E+00 -1.03E+00  7.59E+00  2.09E+02
 
 TH 7
+       -3.46E-02  0.00E+00  2.69E-03  1.61E-03 -6.13E-02  9.01E-02  1.97E-03
 
 TH 8
+       -4.21E-01  0.00E+00 -1.60E+01 -4.73E+00 -8.19E+00  5.67E-01 -1.48E-02  2.11E+01
 
 TH 9
+       -3.31E+00  0.00E+00  6.17E+00 -8.38E-01 -2.55E+00  7.74E+00  6.75E-02  6.01E-01  5.24E+02
 
 TH10
+        9.02E+00  0.00E+00 -1.06E+00 -4.38E+00 -9.28E+01 -1.07E+01  1.23E-02  7.67E+00  1.91E+00  6.16E+01
 
 TH11
+       -3.58E+01  0.00E+00 -1.57E+00 -1.47E+01  6.50E+00  7.33E-02  4.79E-02  5.62E+00  2.46E+01  2.67E+01  2.46E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.413
Stop Time:
Sat Sep 25 13:03:22 CDT 2021
