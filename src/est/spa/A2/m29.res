Sat Sep 18 09:46:14 CDT 2021
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
$DATA ../../../../data/spa/A2/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -816.018534139288        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.0778E+01  6.1192E+00 -9.3839E+00 -1.1377E+01  1.5380E+02 -1.1540E+01 -7.2821E+01  2.1137E+00 -1.1629E+02 -8.3449E+01
            -1.3920E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1304.82525676663        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0356E+00  9.1779E-01  1.1742E+00  1.1758E+00  8.9888E-01  9.2075E-01  1.3054E+00  8.4759E-01  1.5679E+00  8.2382E-01
             3.7295E+00
 PARAMETER:  1.3499E-01  1.4216E-02  2.6056E-01  2.6191E-01 -6.6093E-03  1.7435E-02  3.6647E-01 -6.5352E-02  5.4971E-01 -9.3802E-02
             1.4163E+00
 GRADIENT:   6.0947E+01  1.9356E+01  1.0607E+01  3.7971E+01 -3.8155E+01 -1.8506E+01  1.2509E+01  5.1852E+00  5.1940E+01  1.5745E+01
             1.4317E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1321.52366877107        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0129E+00  5.4606E-01  8.8811E-01  1.4204E+00  6.4997E-01  9.3940E-01  2.8612E+00  1.7170E-01  1.0420E+00  5.3711E-01
             3.3822E+00
 PARAMETER:  1.1277E-01 -5.0502E-01 -1.8657E-02  4.5096E-01 -3.3083E-01  3.7486E-02  1.1512E+00 -1.6620E+00  1.4111E-01 -5.2156E-01
             1.3185E+00
 GRADIENT:   1.9901E+01  5.2652E+01  1.8814E+01  7.4797E+01 -4.9393E+01 -1.4865E+01  3.4063E+01  2.7710E-01  3.3470E+00  4.2150E+00
             9.6058E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1352.69107007373        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.8352E-01  3.5140E-01  5.5256E-01  1.3675E+00  4.5327E-01  1.0112E+00  1.9958E+00  1.0000E-02  1.0417E+00  7.0849E-01
             2.4145E+00
 PARAMETER:  8.3380E-02 -9.4584E-01 -4.9319E-01  4.1300E-01 -6.9126E-01  1.1114E-01  7.9105E-01 -4.6848E+00  1.4086E-01 -2.4462E-01
             9.8150E-01
 GRADIENT:  -1.8047E+01  2.1958E+01  1.4716E+01  5.8522E+01 -2.4244E+01  1.4096E+00 -8.2243E-01  0.0000E+00 -1.0050E+01  2.9508E+00
            -2.8575E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1356.37942972968        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      363
 NPARAMETR:  1.0012E+00  2.8253E-01  3.6945E-01  1.2706E+00  3.3646E-01  1.0111E+00  1.9431E+00  1.0000E-02  1.0611E+00  6.2711E-01
             2.4562E+00
 PARAMETER:  1.0116E-01 -1.1640E+00 -8.9574E-01  3.3946E-01 -9.8927E-01  1.1100E-01  7.6427E-01 -5.3891E+00  1.5931E-01 -3.6663E-01
             9.9863E-01
 GRADIENT:   3.1108E+00  6.5079E+00 -1.3956E+00  3.6874E+00 -1.7041E+00 -1.9406E-02 -2.2952E+00  0.0000E+00 -5.7584E+00 -7.5662E-01
            -5.4055E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1361.03281624064        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      539
 NPARAMETR:  9.8934E-01  1.1953E-01  4.1822E-01  1.3717E+00  3.4593E-01  9.9829E-01  3.7166E+00  1.0000E-02  1.0352E+00  6.3711E-01
             2.5063E+00
 PARAMETER:  8.9283E-02 -2.0242E+00 -7.7174E-01  4.1608E-01 -9.6151E-01  9.8287E-02  1.4128E+00 -6.8632E+00  1.3458E-01 -3.5082E-01
             1.0188E+00
 GRADIENT:  -3.8003E+00  8.3663E-01  4.8475E+00  1.2143E+01 -9.7313E+00 -1.6929E+00 -1.9304E+00  0.0000E+00 -1.5546E+00  1.0320E+00
             5.4759E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1361.45963250419        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  9.9017E-01  9.1905E-02  3.9895E-01  1.3528E+00  3.3261E-01  1.0074E+00  4.4011E+00  1.0000E-02  1.0356E+00  6.4320E-01
             2.4510E+00
 PARAMETER:  9.0124E-02 -2.2870E+00 -8.1891E-01  4.0214E-01 -1.0008E+00  1.0741E-01  1.5819E+00 -7.5173E+00  1.3498E-01 -3.4129E-01
             9.9648E-01
 GRADIENT:   3.6812E+00  6.6242E+00 -3.0447E+00 -1.6170E+01  5.9024E+00  9.4381E-01  1.0077E+01  0.0000E+00  3.9537E-02 -3.3476E+00
            -9.0082E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1361.78671937335        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      891
 NPARAMETR:  9.8726E-01  6.0008E-02  3.7987E-01  1.3593E+00  3.1978E-01  1.0025E+00  5.4353E+00  1.0000E-02  1.0335E+00  6.2138E-01
             2.4852E+00
 PARAMETER:  8.7183E-02 -2.7133E+00 -8.6792E-01  4.0699E-01 -1.0401E+00  1.0247E-01  1.7929E+00 -8.5161E+00  1.3298E-01 -3.7581E-01
             1.0104E+00
 GRADIENT:  -1.4997E+00  7.3832E-01 -2.7347E-01  4.3830E+00 -1.7989E+00 -4.5416E-01  7.1989E-01  0.0000E+00 -3.1799E-01 -3.3991E-01
             2.7894E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1361.83830805527        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1067
 NPARAMETR:  9.8693E-01  4.8229E-02  3.7568E-01  1.3555E+00  3.1703E-01  1.0036E+00  5.9928E+00  1.0000E-02  1.0360E+00  6.1955E-01
             2.4687E+00
 PARAMETER:  8.6840E-02 -2.9318E+00 -8.7901E-01  4.0416E-01 -1.0487E+00  1.0355E-01  1.8906E+00 -9.0941E+00  1.3537E-01 -3.7877E-01
             1.0037E+00
 GRADIENT:   1.9039E-01  1.1610E+00 -3.0124E+00 -2.1765E+00  4.4575E+00 -1.8213E-01  2.4154E+00  0.0000E+00  5.9734E-01 -1.2468E+00
            -1.1520E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1361.84223118194        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1242
 NPARAMETR:  9.8674E-01  4.4111E-02  3.7485E-01  1.3560E+00  3.1612E-01  1.0034E+00  6.2519E+00  1.0000E-02  1.0349E+00  6.1992E-01
             2.4663E+00
 PARAMETER:  8.6656E-02 -3.0210E+00 -8.8123E-01  4.0453E-01 -1.0516E+00  1.0344E-01  1.9329E+00 -9.3193E+00  1.3431E-01 -3.7817E-01
             1.0027E+00
 GRADIENT:  -4.4801E-03 -2.7460E-02  2.6599E-02  2.0403E-02 -5.3060E-02 -7.9965E-04 -5.8178E-02  0.0000E+00 -3.8368E-02  2.1357E-02
             4.1848E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1361.84223118194        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1264
 NPARAMETR:  9.8674E-01  4.4111E-02  3.7485E-01  1.3560E+00  3.1612E-01  1.0034E+00  6.2519E+00  1.0000E-02  1.0349E+00  6.1992E-01
             2.4663E+00
 PARAMETER:  8.6656E-02 -3.0210E+00 -8.8123E-01  4.0453E-01 -1.0516E+00  1.0344E-01  1.9329E+00 -9.3193E+00  1.3431E-01 -3.7817E-01
             1.0027E+00
 GRADIENT:  -4.4801E-03 -2.7460E-02  2.6599E-02  2.0403E-02 -5.3060E-02 -7.9965E-04 -5.8178E-02  0.0000E+00 -3.8368E-02  2.1357E-02
             4.1848E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1264
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.8208E-04  1.8771E-02 -7.4530E-06 -1.3763E-02  2.6372E-03
 SE:             2.9105E-02  9.9436E-03  2.2404E-04  2.7563E-02  2.0251E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8678E-01  5.9060E-02  9.7346E-01  6.1755E-01  8.9639E-01

 ETASHRINKSD(%)  2.4931E+00  6.6688E+01  9.9249E+01  7.6597E+00  3.2157E+01
 ETASHRINKVR(%)  4.9240E+00  8.8903E+01  9.9994E+01  1.4733E+01  5.3973E+01
 EBVSHRINKSD(%)  2.2807E+00  7.6511E+01  9.9177E+01  7.2387E+00  2.9860E+01
 EBVSHRINKVR(%)  4.5094E+00  9.4482E+01  9.9993E+01  1.3953E+01  5.0804E+01
 RELATIVEINF(%)  9.4482E+01  3.8438E+00  2.5656E-04  3.1217E+01  1.8860E+00
 EPSSHRINKSD(%)  3.4502E+01
 EPSSHRINKVR(%)  5.7100E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1361.8422311819374     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -626.69140461819927     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.36
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1361.842       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.87E-01  4.41E-02  3.75E-01  1.36E+00  3.16E-01  1.00E+00  6.25E+00  1.00E-02  1.03E+00  6.20E-01  2.47E+00
 


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
+        1.12E+03
 
 TH 2
+        1.36E+04  1.06E+07
 
 TH 3
+       -5.67E+02 -4.29E+06  1.59E+04
 
 TH 4
+       -1.68E+02 -2.58E+06  2.24E+03  1.16E+03
 
 TH 5
+        1.88E+03  4.26E+06 -1.74E+06 -1.04E+06  1.73E+06
 
 TH 6
+       -3.03E+01 -8.02E+03  3.30E+02  7.64E+01 -1.02E+03  1.83E+02
 
 TH 7
+       -1.67E+01  1.17E+05 -4.71E+04  6.56E+01  4.70E+04  9.90E+00  1.29E+03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.91E+01  1.31E+04 -4.63E+02 -1.42E+02  1.66E+03 -2.10E+01 -1.12E+01  0.00E+00  1.69E+02
 
 TH10
+       -1.73E+02 -6.03E+06  2.78E+03  7.71E+02 -1.02E+04  1.07E+02  5.89E+01  0.00E+00 -1.58E+02  1.01E+03
 
 TH11
+       -7.39E+01 -5.73E+05  1.15E+03  2.82E+02 -2.30E+05  3.90E+01 -6.28E+03  0.00E+00 -5.23E+01  3.55E+02  1.74E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.561
Stop Time:
Sat Sep 18 09:46:39 CDT 2021
