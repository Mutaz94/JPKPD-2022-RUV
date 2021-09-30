Wed Sep 29 17:57:40 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1650.64023845834        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5073E+02 -3.3033E+01  2.4154E+01 -4.4919E+01 -2.4638E+01  4.4747E+01 -1.4222E+01  1.0477E+00  5.1462E+00 -3.2888E+01
            -3.0261E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1663.46573579492        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0159E+00  1.1410E+00  1.0082E+00  1.0119E+00  1.1234E+00  1.0102E+00  1.1583E+00  9.6182E-01  9.2618E-01  1.3549E+00
             1.1442E+00
 PARAMETER:  1.1574E-01  2.3190E-01  1.0821E-01  1.1178E-01  2.1633E-01  1.1019E-01  2.4696E-01  6.1075E-02  2.3314E-02  4.0373E-01
             2.3471E-01
 GRADIENT:  -4.8306E+01  9.4647E+00 -1.3892E+01  2.3847E+01 -7.2657E+00  1.0491E+00 -3.0465E-01  5.3272E+00  7.1991E-01  8.6756E+00
             3.0045E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1664.73332467817        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0332E+00  1.2396E+00  1.0273E+00  9.4154E-01  1.2587E+00  1.0097E+00  1.1552E+00  8.8726E-01  1.0048E+00  1.4672E+00
             1.0822E+00
 PARAMETER:  1.3265E-01  3.1479E-01  1.2694E-01  3.9767E-02  3.3008E-01  1.0969E-01  2.4429E-01 -1.9619E-02  1.0474E-01  4.8337E-01
             1.7900E-01
 GRADIENT:  -9.3921E+00 -2.1153E+00 -1.2934E+01  1.3151E+01  2.2976E+01  1.6417E+00  8.8185E+00  5.5202E-01  5.4047E+00  5.8305E+00
             5.6604E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1666.06451620698        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0398E+00  1.3042E+00  9.0010E-01  8.9546E-01  1.1639E+00  1.0064E+00  1.0476E+00  6.9916E-01  1.0333E+00  1.3323E+00
             1.0725E+00
 PARAMETER:  1.3907E-01  3.6558E-01 -5.2503E-03 -1.0415E-02  2.5175E-01  1.0636E-01  1.4650E-01 -2.5788E-01  1.3273E-01  3.8691E-01
             1.7003E-01
 GRADIENT:   1.9620E+00  5.3404E+00 -7.3785E-01  8.1555E+00 -1.6516E+00 -2.6815E-01  9.7446E-01  3.0814E-01  1.2777E+00  8.2460E-01
             1.3069E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1666.57852117824        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0391E+00  1.6691E+00  6.2889E-01  6.5255E-01  1.2509E+00  1.0113E+00  8.6517E-01  2.4161E-01  1.2659E+00  1.3516E+00
             1.0721E+00
 PARAMETER:  1.3831E-01  6.1228E-01 -3.6380E-01 -3.2687E-01  3.2385E-01  1.1123E-01 -4.4831E-02 -1.3204E+00  3.3576E-01  4.0127E-01
             1.6963E-01
 GRADIENT:  -4.0091E+00 -2.2716E-01 -3.9023E+00  7.0235E+00  2.9035E+00  5.8169E-01 -8.0706E-01  8.3649E-02 -1.1472E+00  8.6557E-01
             8.4706E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1666.60415947601        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  1.0395E+00  1.7872E+00  5.6907E-01  5.8105E-01  1.2975E+00  1.0110E+00  8.2454E-01  1.6277E-01  1.3702E+00  1.3677E+00
             1.0700E+00
 PARAMETER:  1.3876E-01  6.8063E-01 -4.6376E-01 -4.4292E-01  3.6045E-01  1.1096E-01 -9.2935E-02 -1.7154E+00  4.1498E-01  4.1311E-01
             1.6770E-01
 GRADIENT:  -3.4602E+00  1.1037E+01 -1.4426E+00  9.4086E+00  1.1883E+00  3.4268E-01 -1.5062E+00  3.1792E-02 -1.8758E+00 -5.6282E-01
            -8.5674E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1666.74142366103        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1082
 NPARAMETR:  1.0406E+00  1.7895E+00  5.6999E-01  5.7074E-01  1.3029E+00  1.0097E+00  8.2535E-01  5.6655E-02  1.4046E+00  1.3740E+00
             1.0725E+00
 PARAMETER:  1.3977E-01  6.8191E-01 -4.6214E-01 -4.6082E-01  3.6462E-01  1.0963E-01 -9.1942E-02 -2.7708E+00  4.3972E-01  4.1773E-01
             1.7000E-01
 GRADIENT:  -8.8585E-01 -1.5238E+00  5.0087E-02  2.0495E+00 -1.1844E+00 -1.1082E-01  1.9419E-02  3.5958E-03  5.1997E-02 -9.4209E-02
             1.7075E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1666.74193772725        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1263
 NPARAMETR:  1.0405E+00  1.8004E+00  5.6635E-01  5.6405E-01  1.3101E+00  1.0094E+00  8.2144E-01  3.3507E-02  1.4179E+00  1.3793E+00
             1.0729E+00
 PARAMETER:  1.3970E-01  6.8800E-01 -4.6854E-01 -4.7261E-01  3.7010E-01  1.0936E-01 -9.6700E-02 -3.2960E+00  4.4918E-01  4.2157E-01
             1.7034E-01
 GRADIENT:  -1.0369E+00 -8.4888E-01  2.8110E-01  2.2657E+00 -1.0046E+00 -2.1290E-01 -7.3182E-02  1.1989E-03  6.0935E-02 -1.0381E-01
             1.3229E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1666.75592948314        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1451             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0423E+00  1.7991E+00  5.6510E-01  5.6154E-01  1.3122E+00  1.0103E+00  8.2133E-01  1.0000E-02  1.4203E+00  1.3810E+00
             1.0721E+00
 PARAMETER:  1.4147E-01  6.8731E-01 -4.7076E-01 -4.7708E-01  3.7168E-01  1.1023E-01 -9.6826E-02 -4.9416E+00  4.5088E-01  4.2278E-01
             1.6966E-01
 GRADIENT:   5.4395E+02  6.1095E+02  2.5668E+00  9.2279E+01  1.6813E+01  4.5100E+01  5.2270E+00  0.0000E+00  1.5091E+01  5.4134E+00
             1.3387E+00

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1666.75592948314        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:     1512
 NPARAMETR:  1.0424E+00  1.7957E+00  5.6621E-01  5.5905E-01  1.3126E+00  1.0104E+00  8.2051E-01  1.0000E-02  1.4179E+00  1.3817E+00
             1.0725E+00
 PARAMETER:  1.4147E-01  6.8731E-01 -4.7076E-01 -4.7708E-01  3.7168E-01  1.1023E-01 -9.6826E-02 -4.9416E+00  4.5088E-01  4.2278E-01
             1.6966E-01
 GRADIENT:  -1.8295E-02  5.5960E-01 -5.1855E-02  3.2701E-01 -2.8155E-02 -8.1822E-03  2.1692E-02  0.0000E+00  3.2528E-02 -1.2164E-02
            -1.7576E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1512
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.4100E-05 -3.4919E-02 -3.2541E-04  2.9031E-02 -4.5375E-02
 SE:             2.9827E-02  2.3399E-02  1.0377E-04  2.2178E-02  2.2673E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9909E-01  1.3562E-01  1.7130E-03  1.9055E-01  4.5358E-02

 ETASHRINKSD(%)  7.5952E-02  2.1610E+01  9.9652E+01  2.5700E+01  2.4044E+01
 ETASHRINKVR(%)  1.5185E-01  3.8550E+01  9.9999E+01  4.4795E+01  4.2307E+01
 EBVSHRINKSD(%)  4.8530E-01  2.0391E+01  9.9717E+01  2.8713E+01  2.0600E+01
 EBVSHRINKVR(%)  9.6825E-01  3.6623E+01  9.9999E+01  4.9181E+01  3.6957E+01
 RELATIVEINF(%)  9.8908E+01  4.1499E+00  1.0508E-04  3.2395E+00  1.8209E+01
 EPSSHRINKSD(%)  4.3903E+01
 EPSSHRINKVR(%)  6.8531E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1666.7559294831392     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -931.60510291940102     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.20
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.61
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1666.756       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.80E+00  5.65E-01  5.62E-01  1.31E+00  1.01E+00  8.21E-01  1.00E-02  1.42E+00  1.38E+00  1.07E+00
 


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
+        9.92E+02
 
 TH 2
+       -5.28E+00  3.46E+02
 
 TH 3
+        1.18E+01  1.25E+02  3.21E+02
 
 TH 4
+       -1.56E+01  3.16E+02 -2.57E+02  9.08E+02
 
 TH 5
+       -1.07E+00 -8.79E+01 -1.60E+02  1.45E+02  2.17E+02
 
 TH 6
+       -3.33E-02 -9.25E-01  3.02E+00 -4.24E+00 -3.53E-01  1.92E+02
 
 TH 7
+        4.50E-01  6.64E+00 -9.09E+00 -1.49E+01 -6.74E+00 -6.67E-01  1.25E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.24E-01 -1.80E+01 -3.47E+01  5.79E+01  9.75E-01 -6.55E-01  1.73E+01  0.00E+00  3.70E+01
 
 TH10
+        5.07E-01 -9.63E+00 -2.43E+01  1.81E+00 -3.50E+01  2.77E-01  7.59E+00  0.00E+00  4.12E+00  4.74E+01
 
 TH11
+       -7.43E+00 -2.14E+01 -3.90E+01  6.75E+00 -1.65E+00  2.86E+00  1.34E+01  0.00E+00  6.09E+00  7.66E+00  1.88E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.865
Stop Time:
Wed Sep 29 17:58:10 CDT 2021
