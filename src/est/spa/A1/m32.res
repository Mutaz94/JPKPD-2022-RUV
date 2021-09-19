Sat Sep 18 09:11:12 CDT 2021
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
$DATA ../../../../data/spa/A1/dat32.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m32.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1519.41554505378        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.7211E+01 -3.5652E+01  1.4322E+01 -5.5523E+01  5.9083E+01  3.8397E+01 -6.1970E-01 -1.5046E+01 -4.5866E+00  9.5275E-01
            -3.2656E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1589.18795255060        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0389E+00  9.7489E-01  9.1329E-01  1.0297E+00  9.0189E-01  8.4096E-01  9.5430E-01  1.0486E+00  9.3926E-01  8.3106E-01
             1.3837E+00
 PARAMETER:  1.3819E-01  7.4565E-02  9.2987E-03  1.2928E-01 -3.2607E-03 -7.3206E-02  5.3224E-02  1.4748E-01  3.7333E-02 -8.5052E-02
             4.2473E-01
 GRADIENT:   2.5945E+01 -1.8269E+01  1.0758E+01 -4.0545E+01  1.5166E+01 -2.3567E+01 -2.3142E+00 -4.5730E+00 -7.0516E+00  6.4313E+00
            -7.3057E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1597.07742808294        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0338E+00  7.4554E-01  5.6278E-01  1.1751E+00  6.0184E-01  8.9314E-01  1.6425E+00  4.5714E-01  8.0190E-01  4.0755E-01
             1.4209E+00
 PARAMETER:  1.3325E-01 -1.9364E-01 -4.7486E-01  2.6136E-01 -4.0776E-01 -1.3017E-02  5.9621E-01 -6.8278E-01 -1.2078E-01 -7.9760E-01
             4.5126E-01
 GRADIENT:  -1.8820E+00  3.7388E+01 -5.3701E+00  9.0950E+01  3.2827E+01 -1.9092E-01  2.1679E+01 -6.2795E+00 -9.1968E-01 -6.8917E+00
            -5.8331E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1605.19553000340        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0367E+00  8.8191E-01  3.3349E-01  9.8824E-01  5.0029E-01  8.9665E-01  1.0812E+00  2.2125E-01  8.3413E-01  3.1324E-01
             1.5767E+00
 PARAMETER:  1.3602E-01 -2.5670E-02 -9.9813E-01  8.8167E-02 -5.9257E-01 -9.0894E-03  1.7803E-01 -1.4085E+00 -8.1368E-02 -1.0608E+00
             5.5533E-01
 GRADIENT:   1.7892E+00  1.2349E+01  7.0302E+00  7.5355E+00 -1.4613E+01  9.0360E-01  2.0632E-01 -1.8980E-01  6.8730E-01  1.2946E+00
             7.5687E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1605.93483291417        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      411
 NPARAMETR:  1.0464E+00  9.6459E-01  3.5069E-01  9.5133E-01  5.4766E-01  8.9708E-01  1.0420E+00  4.9357E-01  8.5702E-01  2.0276E-01
             1.5352E+00
 PARAMETER:  1.4535E-01  6.3945E-02 -9.4784E-01  5.0107E-02 -5.0209E-01 -8.6115E-03  1.4117E-01 -6.0609E-01 -5.4298E-02 -1.4957E+00
             5.2865E-01
 GRADIENT:   7.6797E+00 -6.7981E-01 -1.7452E+00  6.1825E+00  7.7546E+00  1.2756E-01 -3.4928E-01 -5.7309E-01 -1.6211E+00 -1.2544E-02
            -9.5317E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1606.31428387399        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      587
 NPARAMETR:  1.0392E+00  1.0729E+00  2.9845E-01  8.7092E-01  5.5460E-01  8.9408E-01  9.4153E-01  5.7552E-01  8.9753E-01  1.0911E-01
             1.5762E+00
 PARAMETER:  1.3841E-01  1.7034E-01 -1.1091E+00 -3.8209E-02 -4.8950E-01 -1.1962E-02  3.9753E-02 -4.5249E-01 -8.1102E-03 -2.1154E+00
             5.5503E-01
 GRADIENT:  -2.9556E+00  4.4581E+00  2.6712E+00  2.4546E-02 -7.8821E+00  3.4334E-02 -1.4416E-02  1.8161E-01  1.1832E+00  1.8398E-01
             4.9593E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1607.54932463941        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      767
 NPARAMETR:  1.0328E+00  1.3781E+00  2.1322E-01  6.7835E-01  6.3044E-01  8.8822E-01  7.7550E-01  1.1523E+00  1.0180E+00  1.6689E-02
             1.5064E+00
 PARAMETER:  1.3231E-01  4.2072E-01 -1.4454E+00 -2.8809E-01 -3.6134E-01 -1.8533E-02 -1.5424E-01  2.4179E-01  1.1779E-01 -3.9930E+00
             5.0974E-01
 GRADIENT:   4.9388E+00  2.7495E+01  5.8257E+00  4.7916E+00 -3.4247E+01  1.7215E+00 -3.4986E+00 -2.1809E+00  1.3167E+00  6.4748E-03
            -1.1084E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1607.69076844905        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      961
 NPARAMETR:  1.0326E+00  1.3879E+00  2.1357E-01  6.7519E-01  6.3518E-01  8.8766E-01  7.7299E-01  1.1976E+00  1.0199E+00  1.6103E-02
             1.5053E+00
 PARAMETER:  1.3210E-01  4.2778E-01 -1.4438E+00 -2.9276E-01 -3.5385E-01 -1.9165E-02 -1.5749E-01  2.8033E-01  1.1966E-01 -4.0287E+00
             5.0901E-01
 GRADIENT:   4.3796E+00  3.0128E+01  5.3525E+00  8.3113E+00 -3.3484E+01  1.6117E+00 -3.5108E+00 -2.0496E+00  1.4330E+00  5.7246E-03
            -2.8425E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1607.71640173009        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     1104
 NPARAMETR:  1.0327E+00  1.3876E+00  2.1372E-01  6.7509E-01  6.3506E-01  8.8771E-01  7.8450E-01  1.1974E+00  1.0198E+00  1.0000E-02
             1.5057E+00
 PARAMETER:  1.3216E-01  4.2757E-01 -1.4431E+00 -2.9290E-01 -3.5403E-01 -1.9115E-02 -1.4271E-01  2.8019E-01  1.1960E-01 -4.6506E+00
             5.0926E-01
 GRADIENT:   4.3681E+00  2.9766E+01  5.5439E+00  6.2143E+00 -3.4130E+01  1.6320E+00 -6.8035E-04 -1.9179E+00  1.6513E+00  0.0000E+00
             2.9611E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1607.73356954646        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     1304
 NPARAMETR:  1.0326E+00  1.3870E+00  2.1386E-01  6.7492E-01  6.3504E-01  8.8746E-01  7.8438E-01  1.1992E+00  1.0195E+00  1.0000E-02
             1.5062E+00
 PARAMETER:  1.3204E-01  4.2740E-01 -1.4433E+00 -2.9297E-01 -3.5385E-01 -1.9454E-02 -1.4272E-01  2.8187E-01  1.1943E-01 -4.6506E+00
             5.0926E-01
 GRADIENT:  -3.6115E+05  5.5814E+04 -3.3010E+04  8.1384E+04  1.3473E+05 -4.7688E+05  1.4631E-02  1.6909E+05  1.9964E+05  0.0000E+00
            -9.3638E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1304
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.3386E-03 -1.5910E-02 -1.2438E-02  6.8478E-04 -3.0236E-04
 SE:             2.9511E-02  2.6418E-02  1.4614E-02  2.4268E-02  3.2821E-04
 N:                     100         100         100         100         100

 P VAL.:         9.6382E-01  5.4702E-01  3.9471E-01  9.7749E-01  3.5692E-01

 ETASHRINKSD(%)  1.1357E+00  1.1496E+01  5.1041E+01  1.8698E+01  9.8900E+01
 ETASHRINKVR(%)  2.2585E+00  2.1670E+01  7.6031E+01  3.3900E+01  9.9988E+01
 EBVSHRINKSD(%)  1.0864E+00  1.1324E+01  5.2955E+01  1.7772E+01  9.8872E+01
 EBVSHRINKVR(%)  2.1609E+00  2.1365E+01  7.7867E+01  3.2386E+01  9.9987E+01
 RELATIVEINF(%)  9.5927E+01  5.2614E+00  4.7175E+00  6.4789E+00  9.2836E-04
 EPSSHRINKSD(%)  4.1195E+01
 EPSSHRINKVR(%)  6.5419E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1607.7335695464617     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -872.58274298272352     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.62
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0T MATRIX UNOBTAINABLE
 Elapsed covariance  time in seconds:     5.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1607.734       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.39E+00  2.14E-01  6.75E-01  6.35E-01  8.87E-01  7.84E-01  1.20E+00  1.02E+00  1.00E-02  1.51E+00
 


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
+        6.41E+08
 
 TH 2
+        1.34E+03  3.39E+07
 
 TH 3
+       -2.71E+03  9.28E+03  1.25E+08
 
 TH 4
+        4.02E+03 -1.31E+04  2.98E+05  3.05E+08
 
 TH 5
+        3.54E+03 -1.29E+04  2.61E+05 -5.33E+04  2.36E+08
 
 TH 6
+        9.85E+08  2.15E+03 -4.15E+03  6.45E+03  5.68E+03  1.51E+09
 
 TH 7
+       -1.44E+04  3.32E+03 -6.48E+03  9.93E+03  8.80E+03 -2.22E+04  2.02E+02
 
 TH 8
+        2.38E+03 -7.81E+03  1.22E+05 -1.90E+05 -1.67E+05  3.78E+03  5.82E+03  1.04E+08
 
 TH 9
+        6.61E+03  1.65E+08  2.20E+04 -3.42E+04 -3.02E+04  1.05E+04  1.62E+04 -2.00E+04  8.04E+08
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.06E+03  3.43E+03 -7.74E+04  6.30E+03  5.53E+03 -1.67E+03 -2.55E+03  4.90E+04  8.85E+03  0.00E+00  2.03E+07
 
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
 #CPUT: Total CPU Time in Seconds,       21.478
Stop Time:
Sat Sep 18 09:11:35 CDT 2021
