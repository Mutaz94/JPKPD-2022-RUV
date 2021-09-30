Wed Sep 29 22:25:29 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat39.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m39.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1633.19362332818        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5790E+02  7.1386E+01  4.6863E+01  5.6373E+01  3.8153E+01  5.5156E+01 -2.3148E+01  6.7916E+00  4.4243E+00 -3.8621E+01
            -8.9489E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1830.41210132181        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      150
 NPARAMETR:  1.0195E+00  8.8126E-01  1.0015E+00  1.1617E+00  9.4059E-01  9.3110E-01  1.0254E+00  8.0563E-01  9.6650E-01  1.0029E+00
             2.3346E+00
 PARAMETER:  1.1930E-01 -2.6405E-02  1.0154E-01  2.4992E-01  3.8753E-02  2.8611E-02  1.2505E-01 -1.1614E-01  6.5923E-02  1.0292E-01
             9.4783E-01
 GRADIENT:   3.9769E+01  4.4959E+01 -1.8693E+01  8.8505E+01  4.8206E+00 -4.4862E+00  1.2099E+00  1.1003E+01  2.4367E+00  1.2751E+01
             1.8009E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1852.39608080589        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      327
 NPARAMETR:  1.0290E+00  7.4140E-01  6.1416E-01  1.1913E+00  6.4472E-01  9.7218E-01  1.4087E+00  1.4095E-01  9.1396E-01  7.8140E-01
             2.1446E+00
 PARAMETER:  1.2862E-01 -1.9922E-01 -3.8751E-01  2.7501E-01 -3.3894E-01  7.1782E-02  4.4268E-01 -1.8594E+00  1.0028E-02 -1.4667E-01
             8.6294E-01
 GRADIENT:   5.9096E+01  4.4964E+01 -3.4763E+01  8.4822E+01  2.4795E+01  9.3233E+00  1.3632E+01  6.7801E-01  4.6868E+00  1.9155E+01
             1.5301E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1885.68171476634        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      504
 NPARAMETR:  1.0006E+00  4.9161E-01  4.6047E-01  1.2427E+00  4.5182E-01  9.4485E-01  1.4783E+00  5.3359E-02  8.8963E-01  5.5915E-01
             1.6858E+00
 PARAMETER:  1.0057E-01 -6.1007E-01 -6.7550E-01  3.1730E-01 -6.9447E-01  4.3271E-02  4.9089E-01 -2.8307E+00 -1.6952E-02 -4.8133E-01
             6.2222E-01
 GRADIENT:   1.0115E+01  4.3826E+01  4.2433E+01  6.1640E+01 -6.9199E+01 -7.1646E-01 -2.0597E+00  9.3465E-02 -3.9199E+00 -2.8584E+00
            -1.6976E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1891.66830827063        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      679
 NPARAMETR:  9.8510E-01  2.6862E-01  4.8333E-01  1.3071E+00  4.2731E-01  9.2812E-01  2.1151E+00  1.7434E-02  8.7095E-01  6.0328E-01
             1.7197E+00
 PARAMETER:  8.4987E-02 -1.2145E+00 -6.2705E-01  3.6784E-01 -7.5024E-01  2.5403E-02  8.4911E-01 -3.9493E+00 -3.8171E-02 -4.0538E-01
             6.4214E-01
 GRADIENT:  -1.0920E+01  3.9236E+00  1.4316E+01 -1.3360E+01 -2.0866E+01 -4.5115E+00  4.4782E-01  9.5665E-03  2.4744E+00  5.3800E-02
             1.3870E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1893.81626986813        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      855
 NPARAMETR:  9.8427E-01  1.3152E-01  5.5888E-01  1.4121E+00  4.5669E-01  9.4273E-01  3.1794E+00  1.0000E-02  8.2395E-01  6.6016E-01
             1.7117E+00
 PARAMETER:  8.4141E-02 -1.9286E+00 -4.8182E-01  4.4510E-01 -6.8375E-01  4.1026E-02  1.2567E+00 -5.1388E+00 -9.3648E-02 -3.1527E-01
             6.3748E-01
 GRADIENT:   4.4736E+00  3.5540E+00 -4.5017E+00  9.6804E+00  6.1797E+00  2.6910E+00  2.0704E+00  0.0000E+00 -5.6141E+00 -5.4015E-01
            -4.4009E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1894.93589233412        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1030
 NPARAMETR:  9.7875E-01  6.1183E-02  5.3193E-01  1.4266E+00  4.3027E-01  9.3384E-01  4.8158E+00  1.0000E-02  8.2607E-01  6.4296E-01
             1.7260E+00
 PARAMETER:  7.8525E-02 -2.6939E+00 -5.3124E-01  4.5530E-01 -7.4335E-01  3.1554E-02  1.6719E+00 -6.9705E+00 -9.1072E-02 -3.4168E-01
             6.4578E-01
 GRADIENT:  -1.3560E+00  3.7438E+00 -3.7237E-02 -1.1009E+00 -1.9107E+00 -4.6193E-01  6.1047E+00  0.0000E+00 -1.1602E+00 -1.1094E+00
             8.6549E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1895.67940362103        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1207
 NPARAMETR:  9.7788E-01  2.3280E-02  5.4119E-01  1.4455E+00  4.3233E-01  9.3485E-01  7.5051E+00  1.0000E-02  8.0589E-01  6.4343E-01
             1.7197E+00
 PARAMETER:  7.7631E-02 -3.6601E+00 -5.1399E-01  4.6846E-01 -7.3857E-01  3.2629E-02  2.1156E+00 -9.2386E+00 -1.1581E-01 -3.4094E-01
             6.4214E-01
 GRADIENT:   6.8414E-01  6.4533E-01 -3.3804E-01  5.3748E+00 -1.4003E+00  6.2061E-01  1.5355E+00  0.0000E+00 -5.0446E+00  8.5203E-01
            -5.2674E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1895.70237758663        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1272
 NPARAMETR:  9.7785E-01  2.3679E-02  5.3742E-01  1.4423E+00  4.3034E-01  9.3418E-01  7.3995E+00  1.0000E-02  8.1181E-01  6.3928E-01
             1.7207E+00
 PARAMETER:  7.7605E-02 -3.6421E+00 -5.2113E-01  4.6612E-01 -7.4298E-01  3.2015E-02  2.1020E+00 -9.2138E+00 -1.0959E-01 -3.4752E-01
             6.4253E-01
 GRADIENT:   2.3454E-01  4.8858E+01 -1.7613E+02 -1.8716E+02  1.2029E+02  3.2717E-01  8.2746E+01  0.0000E+00 -3.1015E+00 -2.6406E+02
            -1.4199E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1272
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.7990E-04  1.3056E-02 -2.8873E-05 -1.0142E-02 -5.7377E-03
 SE:             2.9549E-02  7.8385E-03  2.5368E-04  2.8813E-02  2.2716E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8434E-01  9.5797E-02  9.0938E-01  7.2485E-01  8.0058E-01

 ETASHRINKSD(%)  1.0082E+00  7.3740E+01  9.9150E+01  3.4714E+00  2.3900E+01
 ETASHRINKVR(%)  2.0062E+00  9.3104E+01  9.9993E+01  6.8224E+00  4.2088E+01
 EBVSHRINKSD(%)  1.0435E+00  8.2512E+01  9.9076E+01  4.0643E+00  2.2460E+01
 EBVSHRINKVR(%)  2.0760E+00  9.6942E+01  9.9991E+01  7.9634E+00  3.9875E+01
 RELATIVEINF(%)  9.7453E+01  2.0291E+00  5.2828E-04  4.4955E+01  3.7928E+00
 EPSSHRINKSD(%)  2.9023E+01
 EPSSHRINKVR(%)  4.9623E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1895.7023775866289     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -976.76384438195623     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.90
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.01
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1895.702       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.78E-01  2.37E-02  5.37E-01  1.44E+00  4.30E-01  9.34E-01  7.40E+00  1.00E-02  8.11E-01  6.39E-01  1.72E+00
 


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
+       -2.43E+02  8.55E+05
 
 TH 3
+        1.20E+01 -9.35E+04  3.22E+04
 
 TH 4
+       -2.94E-01 -3.85E+04 -4.90E+02  5.73E+03
 
 TH 5
+        2.65E+01  7.90E+04 -4.61E+03  1.37E+02  3.05E+04
 
 TH 6
+       -2.93E+00  1.16E+02 -2.49E+01 -1.53E+01  1.85E+01  2.21E+02
 
 TH 7
+       -9.20E-01  3.04E+03  1.25E+01  8.66E+00  4.44E+02  6.92E-01  2.66E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.33E+00 -3.97E+02  1.24E+02  1.22E+01 -5.87E+01  2.78E+00 -2.24E+00  0.00E+00  2.65E+02
 
 TH10
+        5.14E+01 -1.19E+05  3.65E+04 -9.56E+01  1.22E+02 -3.87E+01  5.54E+00  0.00E+00  7.23E+01  4.67E+04
 
 TH11
+       -3.57E+00 -2.39E+04  7.32E+03 -5.74E+01  1.10E+02 -5.40E+00  2.88E+00  0.00E+00  2.34E+01 -3.37E+01  2.01E+03
 
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
 #CPUT: Total CPU Time in Seconds,       27.975
Stop Time:
Wed Sep 29 22:25:59 CDT 2021
