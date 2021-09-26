Sat Sep 25 09:38:51 CDT 2021
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
$DATA ../../../../data/spa/S1/dat1.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1671.63419348557        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3500E+01 -7.0118E+00  3.0982E+01 -4.0503E+01 -4.4151E+01  8.0069E+00  5.0244E-02 -5.9678E-01  7.0307E+00 -2.4477E+01
             2.1029E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1675.35890074269        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0179E+00  1.0520E+00  9.7187E-01  9.8340E-01  1.1070E+00  9.6513E-01  9.9935E-01  9.4984E-01  9.3273E-01  1.3649E+00
             9.1289E-01
 PARAMETER:  1.1770E-01  1.5068E-01  7.1469E-02  8.3263E-02  2.0170E-01  6.4509E-02  9.9346E-02  4.8536E-02  3.0356E-02  4.1106E-01
             8.8629E-03
 GRADIENT:   1.0345E+02 -1.7833E+01 -2.1872E+01 -9.4370E+00  1.3521E+01 -6.3626E+00  1.8736E+00  4.2534E+00 -1.0188E+01  1.3449E+01
            -9.8547E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1678.17441196486        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.9974E-01  1.1923E+00  1.0627E+00  9.1299E-01  1.2159E+00  9.7855E-01  7.6454E-01  1.1803E+00  1.0880E+00  1.4129E+00
             9.1772E-01
 PARAMETER:  9.9739E-02  2.7586E-01  1.6083E-01  8.9683E-03  2.9548E-01  7.8319E-02 -1.6848E-01  2.6580E-01  1.8433E-01  4.4564E-01
             1.4140E-02
 GRADIENT:   5.1887E+01  1.5686E+01 -5.5616E+00  1.2778E+01  4.9924E+00  3.5803E-01 -7.8088E-01  8.3661E-01 -3.5806E+00  5.6061E+00
            -7.4619E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1678.90175107709        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      336
 NPARAMETR:  9.9952E-01  1.3891E+00  9.9583E-01  7.8444E-01  1.2763E+00  9.9396E-01  6.8700E-01  1.2893E+00  1.2655E+00  1.3715E+00
             9.4587E-01
 PARAMETER:  9.9518E-02  4.2863E-01  9.5824E-02 -1.4278E-01  3.4397E-01  9.3946E-02 -2.7541E-01  3.5407E-01  3.3544E-01  4.1588E-01
             4.4347E-02
 GRADIENT:  -2.8909E+00  8.2716E+00  5.1153E+00  5.8184E+00 -2.9450E+00  5.8418E-01 -1.3570E+00 -1.3959E+00  8.4179E-01 -3.8852E+00
             4.5158E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1679.35877914369        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  1.0022E+00  1.5100E+00  8.7674E-01  7.0023E-01  1.2933E+00  9.9495E-01  7.1125E-01  1.2477E+00  1.3378E+00  1.3983E+00
             9.3304E-01
 PARAMETER:  1.0216E-01  5.1208E-01 -3.1550E-02 -2.5634E-01  3.5723E-01  9.4935E-02 -2.4074E-01  3.2133E-01  3.9103E-01  4.3527E-01
             3.0690E-02
 GRADIENT:   1.7370E+00  5.7592E+00  4.7940E+00  6.0703E-02 -7.5295E+00  6.1142E-01  1.4750E+00 -4.5807E-01 -6.0185E-01  8.8124E-01
            -5.7479E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1679.42185974998        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      696            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0018E+00  1.5105E+00  8.6614E-01  7.0025E-01  1.2993E+00  9.9342E-01  7.0539E-01  1.2572E+00  1.3400E+00  1.3916E+00
             9.3376E-01
 PARAMETER:  1.0177E-01  5.1241E-01 -4.3710E-02 -2.5632E-01  3.6184E-01  9.3403E-02 -2.4901E-01  3.2886E-01  3.9269E-01  4.3046E-01
             3.1464E-02
 GRADIENT:   5.2002E+01  4.9406E+01  6.1206E-01  1.2248E+01  1.7028E+00  5.5375E+00  1.9068E+00  5.1846E-01  1.6219E+00  3.9315E-01
            -3.7553E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1679.42946061139        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      860             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0020E+00  1.5108E+00  8.6175E-01  7.0034E-01  1.2997E+00  9.9397E-01  6.9377E-01  1.2523E+00  1.3475E+00  1.3959E+00
             9.3420E-01
 PARAMETER:  1.0198E-01  5.1267E-01 -4.8784E-02 -2.5620E-01  3.6210E-01  9.3955E-02 -2.6562E-01  3.2500E-01  3.9822E-01  4.3354E-01
             3.1939E-02
 GRADIENT:   5.2592E+01  4.9059E+01 -5.6754E-01  1.3737E+01  2.7623E+00  5.6923E+00  9.9071E-01  6.1610E-01  2.0471E+00  7.7123E-01
             2.3993E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1679.43093930011        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      997
 NPARAMETR:  1.0015E+00  1.5108E+00  8.6445E-01  7.0034E-01  1.2997E+00  9.9346E-01  6.9489E-01  1.2523E+00  1.3478E+00  1.3954E+00
             9.3394E-01
 PARAMETER:  1.0149E-01  5.1267E-01 -4.5664E-02 -2.5620E-01  3.6210E-01  9.3437E-02 -2.6400E-01  3.2500E-01  3.9846E-01  4.3321E-01
             3.1661E-02
 GRADIENT:   4.2836E-02  1.2181E+00  6.4318E-02  3.6643E+00 -1.4196E-01 -2.1657E-02 -3.1535E-02  3.8798E-01 -8.5543E-02  7.1490E-03
             2.9610E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1679.43419750574        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1160
 NPARAMETR:  1.0012E+00  1.5102E+00  8.6419E-01  6.9980E-01  1.2996E+00  9.9336E-01  6.9487E-01  1.2501E+00  1.3481E+00  1.3951E+00
             9.3389E-01
 PARAMETER:  1.0118E-01  5.1228E-01 -4.5959E-02 -2.5697E-01  3.6207E-01  9.3338E-02 -2.6402E-01  3.2326E-01  3.9869E-01  4.3299E-01
             3.1603E-02
 GRADIENT:  -6.2409E-01 -2.0353E-01  8.0496E-02  2.7089E+00  1.3691E-02 -5.7675E-02 -1.4419E-02  3.6453E-01 -6.8185E-02 -3.5656E-02
            -3.7703E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1679.43457267418        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1336
 NPARAMETR:  1.0013E+00  1.5107E+00  8.6403E-01  6.9984E-01  1.2999E+00  9.9353E-01  6.9485E-01  1.2500E+00  1.3481E+00  1.3955E+00
             9.3385E-01
 PARAMETER:  1.0132E-01  5.1232E-01 -4.6198E-02 -2.5703E-01  3.6209E-01  9.3410E-02 -2.6394E-01  3.2295E-01  3.9885E-01  4.3305E-01
             3.1606E-02
 GRADIENT:   7.4107E+05 -1.4657E+05 -7.5081E+05 -2.9214E+05 -2.0736E+05 -2.0924E-02  2.8448E+05 -2.3260E+05  1.8822E+05 -1.7338E+05
             7.5083E+05
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         2.9         3.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1336
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.6949E-04 -4.0929E-02 -3.3369E-02  2.0518E-02 -4.8229E-02
 SE:             2.9898E-02  1.8268E-02  1.1797E-02  2.4674E-02  2.2622E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8480E-01  2.5058E-02  4.6774E-03  4.0566E-01  3.3011E-02

 ETASHRINKSD(%)  1.0000E-10  3.8800E+01  6.0477E+01  1.7338E+01  2.4213E+01
 ETASHRINKVR(%)  1.0000E-10  6.2546E+01  8.4379E+01  3.1669E+01  4.2563E+01
 EBVSHRINKSD(%)  3.8635E-01  3.8098E+01  6.3390E+01  1.8482E+01  1.9854E+01
 EBVSHRINKVR(%)  7.7122E-01  6.1681E+01  8.6597E+01  3.3549E+01  3.5766E+01
 RELATIVEINF(%)  9.8997E+01  1.8971E+00  2.1284E+00  3.7702E+00  2.2237E+01
 EPSSHRINKSD(%)  4.5481E+01
 EPSSHRINKVR(%)  7.0276E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1679.4345726741751     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -944.28374611043694     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.23
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1679.435       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.51E+00  8.64E-01  7.00E-01  1.30E+00  9.93E-01  6.95E-01  1.25E+00  1.35E+00  1.40E+00  9.34E-01
 


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
+        1.82E+09
 
 TH 2
+       -7.26E+02  3.14E+07
 
 TH 3
+       -6.43E+03 -1.95E+04  2.51E+09
 
 TH 4
+       -3.10E+03  1.29E+04 -8.45E+04  5.80E+08
 
 TH 5
+       -1.18E+03 -8.77E+02 -3.24E+04 -3.28E+03  8.48E+07
 
 TH 6
+        7.89E+03 -1.04E+03 -9.26E+03 -4.45E+03 -1.70E+03  1.92E+02
 
 TH 7
+        3.04E+03  1.90E+03  8.28E+04  8.19E+03  3.13E+03  4.37E+03  5.58E+08
 
 TH 8
+       -1.38E+03  5.64E+04 -3.77E+04  2.43E+05  9.27E+04 -1.99E+03  3.66E+03  1.15E+08
 
 TH 9
+        1.04E+03  1.48E+04  1.33E+05  6.40E+04  2.44E+04  1.49E+03 -6.27E+04  2.85E+04  6.49E+07
 
 TH10
+       -9.20E+02 -6.77E+02  3.60E+08 -2.90E+03 -1.15E+03 -1.33E+03 -1.69E+08 -1.29E+03  1.90E+04  5.14E+07
 
 TH11
+       -4.67E+04  6.11E+03  5.48E+04  2.63E+04  1.01E+04  8.57E+03 -2.58E+04  1.17E+04 -8.80E+03  7.85E+03  2.15E+09
 
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
 #CPUT: Total CPU Time in Seconds,       23.415
Stop Time:
Sat Sep 25 09:39:16 CDT 2021
