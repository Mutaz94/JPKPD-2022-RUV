Wed Sep 29 12:36:21 CDT 2021
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
$DATA ../../../../data/spa/A2/dat15.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -907.731025195219        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4890E+02 -1.0381E+01  2.7020E+01 -5.6625E+01  1.5636E+02  5.6862E+01 -4.3967E+01 -1.7407E+01 -1.1564E+02 -9.0002E+01
            -1.2922E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1373.10860541291        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0360E+00  9.5585E-01  9.5405E-01  1.1073E+00  8.8389E-01  9.1455E-01  1.0555E+00  9.7785E-01  1.3058E+00  1.0758E+00
             2.0257E+00
 PARAMETER:  1.3539E-01  5.4844E-02  5.2965E-02  2.0192E-01 -2.3421E-02  1.0678E-02  1.5399E-01  7.7600E-02  3.6679E-01  1.7305E-01
             8.0590E-01
 GRADIENT:   1.1203E+02  1.2053E+01  2.0217E+00  3.0189E+01  2.0366E+01  2.7317E+00 -4.8718E-01  3.4094E+00  9.5637E+00 -8.2502E+00
            -1.5016E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1382.47881698669        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      201
 NPARAMETR:  1.0391E+00  6.9920E-01  8.1165E-01  1.2696E+00  7.1541E-01  8.9559E-01  1.0230E+00  2.7014E-01  1.1871E+00  1.0345E+00
             2.0846E+00
 PARAMETER:  1.3833E-01 -2.5782E-01 -1.0868E-01  3.3873E-01 -2.3490E-01 -1.0269E-02  1.2277E-01 -1.2088E+00  2.7150E-01  1.3396E-01
             8.3458E-01
 GRADIENT:  -2.7948E+01  1.1254E+01 -5.1822E-02  1.6715E+01  1.8850E+01 -1.3524E+01 -4.0070E+00 -4.9271E-01 -4.6338E+00 -4.9857E+00
            -1.3677E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1401.64252734871        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  1.0593E+00  4.4960E-01  6.2539E-01  1.3855E+00  5.2035E-01  9.1015E-01  1.4392E+00  1.0000E-02  1.0526E+00  8.0527E-01
             2.6326E+00
 PARAMETER:  1.5762E-01 -6.9939E-01 -3.6938E-01  4.2606E-01 -5.5326E-01  5.8566E-03  4.6412E-01 -4.7893E+00  1.5125E-01 -1.1658E-01
             1.0680E+00
 GRADIENT:   1.8929E+00  1.2566E+01  6.6440E+00  2.1874E+01 -1.2478E+01 -3.3753E+00 -7.4803E-01  0.0000E+00 -1.5395E+00  9.2869E-02
             8.5332E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1406.37918424095        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      555
 NPARAMETR:  1.0516E+00  1.3947E-01  5.5375E-01  1.5058E+00  4.2908E-01  9.1648E-01  3.3194E+00  1.0000E-02  9.8528E-01  7.9139E-01
             2.5228E+00
 PARAMETER:  1.5029E-01 -1.8699E+00 -4.9104E-01  5.0933E-01 -7.4611E-01  1.2781E-02  1.2998E+00 -1.1552E+01  8.5171E-02 -1.3396E-01
             1.0254E+00
 GRADIENT:   6.9268E+00  3.2418E+00 -7.7906E-01  2.2763E+01 -2.0731E+00 -2.0516E-01  8.2598E-01  0.0000E+00 -1.6658E-01 -1.6472E+00
            -6.7548E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1407.75346300146        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  1.0445E+00  2.0139E-02  5.0804E-01  1.5105E+00  3.8769E-01  9.1529E-01  8.8541E+00  1.0000E-02  9.7248E-01  8.0367E-01
             2.5132E+00
 PARAMETER:  1.4352E-01 -3.8051E+00 -5.7719E-01  5.1244E-01 -8.4754E-01  1.1483E-02  2.2809E+00 -2.4865E+01  7.2098E-02 -1.1857E-01
             1.0216E+00
 GRADIENT:   3.1603E+00  2.8170E-01  4.6072E+00 -6.6917E+00 -8.1700E+00  2.4784E-01  1.3522E-01  0.0000E+00  1.4702E+00  2.2604E+00
            -3.5851E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1407.95418325071        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  1.0427E+00  1.0000E-02  5.3621E-01  1.5337E+00  4.0394E-01  9.1403E-01  1.2991E+01  1.0000E-02  9.6033E-01  7.9700E-01
             2.5445E+00
 PARAMETER:  1.4182E-01 -4.5847E+00 -5.2322E-01  5.2768E-01 -8.0648E-01  1.0113E-02  2.6643E+00 -2.9753E+01  5.9519E-02 -1.2690E-01
             1.0339E+00
 GRADIENT:   1.1083E-01  0.0000E+00 -1.6379E-01 -2.2779E+00  2.9720E-01  1.6043E-01  1.2224E-01  0.0000E+00  2.1646E-01  2.2061E-01
             1.4586E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1407.96013926920        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1092
 NPARAMETR:  1.0428E+00  1.0000E-02  5.3624E-01  1.5336E+00  4.0386E-01  9.1368E-01  1.1365E+01  1.0000E-02  9.5978E-01  7.9547E-01
             2.5446E+00
 PARAMETER:  1.4194E-01 -4.5847E+00 -5.2318E-01  5.2759E-01 -8.0669E-01  9.7277E-03  2.5305E+00 -2.9753E+01  5.8951E-02 -1.2882E-01
             1.0340E+00
 GRADIENT:   4.1761E-01  0.0000E+00  2.5309E-01 -2.3631E+00 -1.8168E-01  2.0828E-02  5.7201E-03  0.0000E+00  1.9026E-02 -4.3200E-03
            -5.1325E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1407.96013926920        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1114
 NPARAMETR:  1.0428E+00  1.0000E-02  5.3624E-01  1.5336E+00  4.0386E-01  9.1368E-01  1.1365E+01  1.0000E-02  9.5978E-01  7.9547E-01
             2.5446E+00
 PARAMETER:  1.4194E-01 -4.5847E+00 -5.2318E-01  5.2759E-01 -8.0669E-01  9.7277E-03  2.5305E+00 -2.9753E+01  5.8951E-02 -1.2882E-01
             1.0340E+00
 GRADIENT:   4.1761E-01  0.0000E+00  2.5309E-01 -2.3631E+00 -1.8168E-01  2.0828E-02  5.7201E-03  0.0000E+00  1.9026E-02 -4.3200E-03
            -5.1325E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1114
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.4237E-04 -2.0497E-04  1.2845E-05 -9.7259E-03 -1.0185E-02
 SE:             2.9016E-02  1.4931E-03  2.0436E-04  2.7502E-02  2.1887E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9334E-01  8.9081E-01  9.4988E-01  7.2360E-01  6.4168E-01

 ETASHRINKSD(%)  2.7922E+00  9.4998E+01  9.9315E+01  7.8651E+00  2.6677E+01
 ETASHRINKVR(%)  5.5064E+00  9.9750E+01  9.9995E+01  1.5112E+01  4.6238E+01
 EBVSHRINKSD(%)  2.7807E+00  9.5233E+01  9.9297E+01  7.2656E+00  2.5863E+01
 EBVSHRINKVR(%)  5.4841E+00  9.9773E+01  9.9995E+01  1.4003E+01  4.5037E+01
 RELATIVEINF(%)  7.8446E+01  1.5389E-02  2.1369E-04  1.1808E+01  1.6427E+00
 EPSSHRINKSD(%)  3.3632E+01
 EPSSHRINKVR(%)  5.5952E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1407.9601392692000     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -672.80931270546182     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.32
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1407.960       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  5.36E-01  1.53E+00  4.04E-01  9.14E-01  1.14E+01  1.00E-02  9.60E-01  7.95E-01  2.54E+00
 


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
+        1.16E+03
 
 TH 2
+        0.00E+00  7.39E+02
 
 TH 3
+       -1.24E+01  0.00E+00  1.90E+03
 
 TH 4
+       -3.37E+01  0.00E+00 -1.56E+02  4.43E+02
 
 TH 5
+        8.34E+01  0.00E+00 -3.16E+03 -1.25E+02  5.88E+03
 
 TH 6
+       -1.28E+00  0.00E+00  7.36E+00 -8.08E+00  9.06E-01  2.12E+02
 
 TH 7
+        1.16E-03  0.00E+00  7.77E-04 -6.80E-03  1.48E-02  1.81E-03  1.18E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.42E+00  0.00E+00  2.40E+01 -9.22E+00  7.74E+00  1.62E+00  1.64E-02  0.00E+00  1.56E+02
 
 TH10
+       -7.03E+00  0.00E+00 -3.27E+01  7.86E-01 -1.85E+01  1.75E-01  6.90E-03  0.00E+00 -6.41E-01  9.48E+01
 
 TH11
+       -1.51E+01  0.00E+00 -8.62E+00 -6.49E+00 -1.89E+00  3.60E+00  2.68E-03  0.00E+00  8.84E+00  2.21E+01  4.23E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.470
Stop Time:
Wed Sep 29 12:36:45 CDT 2021
