Wed Sep 29 15:18:53 CDT 2021
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
 GRADIENT:   4.2626E+02 -1.6348E+01 -1.0569E+01 -2.0986E+00  2.0481E+01  5.6622E+01 -3.2019E+00  8.2803E+00 -3.2987E+00  1.2144E+01
             1.9531E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1674.74980587422        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  9.8910E-01  1.0244E+00  1.0264E+00  1.0293E+00  1.0015E+00  9.3935E-01  1.0255E+00  9.6402E-01  1.0309E+00  9.5024E-01
             1.0018E+00
 PARAMETER:  8.9042E-02  1.2415E-01  1.2610E-01  1.2891E-01  1.0150E-01  3.7432E-02  1.2515E-01  6.3353E-02  1.3045E-01  4.8964E-02
             1.0176E-01
 GRADIENT:  -1.6017E+01 -1.4282E+00  2.8720E+00 -6.0453E+00 -3.4052E-01 -3.9523E+00 -2.6402E+00  3.8935E+00 -5.1255E-01  2.5811E+00
            -5.8250E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1676.83806259881        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  9.9402E-01  8.7694E-01  8.1633E-01  1.1198E+00  8.2747E-01  9.5697E-01  1.2549E+00  6.0380E-01  9.3982E-01  7.4978E-01
             1.0026E+00
 PARAMETER:  9.4007E-02 -3.1316E-02 -1.0293E-01  2.1313E-01 -8.9382E-02  5.6014E-02  3.2707E-01 -4.0452E-01  3.7929E-02 -1.8798E-01
             1.0255E-01
 GRADIENT:  -6.3574E+00  5.6846E+00 -1.4689E+01  2.6177E+01  2.1517E+01  2.8834E+00 -4.2076E-01  1.7373E+00 -1.4364E+00 -2.9040E+00
             2.7815E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1677.46458942734        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.9685E-01  7.6884E-01  7.9197E-01  1.1727E+00  7.6326E-01  9.4683E-01  1.3994E+00  4.5137E-01  9.0442E-01  7.4026E-01
             1.0001E+00
 PARAMETER:  9.6845E-02 -1.6287E-01 -1.3323E-01  2.5934E-01 -1.7015E-01  4.5364E-02  4.3602E-01 -6.9546E-01 -4.5698E-04 -2.0075E-01
             1.0007E-01
 GRADIENT:   1.7277E+00  8.1466E+00  1.3771E+00  1.0519E+01 -4.6279E+00 -1.0736E+00  9.0036E-01 -2.3913E-01 -1.3416E-01 -1.3483E+00
            -1.0004E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1677.99856553225        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  9.9116E-01  5.0299E-01  9.0492E-01  1.3346E+00  7.3787E-01  9.4824E-01  1.8279E+00  5.6473E-01  8.3455E-01  7.9188E-01
             1.0041E+00
 PARAMETER:  9.1116E-02 -5.8718E-01  8.7552E-05  3.8863E-01 -2.0399E-01  4.6853E-02  7.0316E-01 -4.7141E-01 -8.0868E-02 -1.3335E-01
             1.0408E-01
 GRADIENT:  -1.8110E+00  3.3471E+00  4.0429E+00  2.8506E+00 -6.8097E+00  1.1316E+00 -8.9192E-01 -1.3909E-01 -1.2260E+00  1.4806E+00
             9.1705E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1678.13305269451        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  9.8963E-01  3.8734E-01  9.2606E-01  1.4055E+00  7.1545E-01  9.4397E-01  2.1857E+00  6.1793E-01  8.1044E-01  7.7799E-01
             1.0005E+00
 PARAMETER:  8.9573E-02 -8.4845E-01  2.3179E-02  4.4039E-01 -2.3485E-01  4.2337E-02  8.8192E-01 -3.8138E-01 -1.1018E-01 -1.5104E-01
             1.0054E-01
 GRADIENT:  -4.2686E-02  4.3143E+00  4.7025E+00  9.0566E+00 -8.2860E+00  1.1974E-01  9.6711E-02  2.4707E-02 -1.0937E+00 -2.1829E-02
            -6.7106E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1678.23705885736        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  9.8711E-01  2.6630E-01  9.6224E-01  1.4764E+00  7.0484E-01  9.4017E-01  2.7207E+00  6.7963E-01  7.9116E-01  7.8153E-01
             9.9970E-01
 PARAMETER:  8.7031E-02 -1.2231E+00  6.1511E-02  4.8962E-01 -2.4978E-01  3.8301E-02  1.1009E+00 -2.8621E-01 -1.3426E-01 -1.4650E-01
             9.9701E-02
 GRADIENT:   7.3364E-01  2.8118E+00  1.9143E+00  6.6225E+00 -3.7729E+00 -5.1166E-01  8.2775E-01  4.8842E-02 -3.0015E-01 -7.1981E-01
            -1.0344E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1678.30360933319        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  9.8721E-01  2.5913E-01  9.6223E-01  1.4717E+00  7.0470E-01  9.4137E-01  2.7490E+00  6.6849E-01  7.9070E-01  7.8794E-01
             1.0015E+00
 PARAMETER:  8.7129E-02 -1.2504E+00  6.1494E-02  4.8644E-01 -2.4999E-01  3.9578E-02  1.1113E+00 -3.0273E-01 -1.3483E-01 -1.3833E-01
             1.0154E-01
 GRADIENT:   1.4330E+00  7.2282E-01  1.4404E+00 -1.0853E+01 -7.0457E-01  7.8635E-02  7.7423E-01 -1.1416E-01 -7.1835E-02 -3.0438E-01
            -1.8345E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1678.30615107798        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     1375
 NPARAMETR:  9.8719E-01  2.5746E-01  9.6166E-01  1.4723E+00  7.0463E-01  9.4136E-01  2.7454E+00  6.6949E-01  7.9071E-01  7.9021E-01
             1.0018E+00
 PARAMETER:  8.7106E-02 -1.2569E+00  6.0909E-02  4.8680E-01 -2.5009E-01  3.9567E-02  1.1099E+00 -3.0123E-01 -1.3482E-01 -1.3546E-01
             1.0178E-01
 GRADIENT:   1.4490E+00  3.2185E-01  8.0823E-02 -1.1394E+01  6.6507E-01  8.8344E-02  3.3219E-01  8.3459E-02  6.8268E-02  7.7896E-02
             7.2968E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1375
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8922E-04  2.3249E-02 -2.1979E-02 -1.5160E-02 -1.4785E-02
 SE:             2.9823E-02  1.4659E-02  1.4202E-02  2.6900E-02  2.1427E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8424E-01  1.1276E-01  1.2172E-01  5.7304E-01  4.9018E-01

 ETASHRINKSD(%)  8.8626E-02  5.0889E+01  5.2422E+01  9.8814E+00  2.8215E+01
 ETASHRINKVR(%)  1.7717E-01  7.5881E+01  7.7363E+01  1.8786E+01  4.8469E+01
 EBVSHRINKSD(%)  4.8847E-01  5.7677E+01  5.3405E+01  8.4887E+00  2.4849E+01
 EBVSHRINKVR(%)  9.7456E-01  8.2088E+01  7.8289E+01  1.6257E+01  4.3523E+01
 RELATIVEINF(%)  9.7818E+01  2.5479E+00  2.5725E+00  1.6354E+01  5.3363E+00
 EPSSHRINKSD(%)  4.4459E+01
 EPSSHRINKVR(%)  6.9152E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1678.3061510779789     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -943.15532451424076     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.10
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1678.306       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.87E-01  2.57E-01  9.62E-01  1.47E+00  7.05E-01  9.41E-01  2.75E+00  6.69E-01  7.91E-01  7.90E-01  1.00E+00
 


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
+        1.27E+03
 
 TH 2
+       -2.94E+01  4.92E+02
 
 TH 3
+        8.49E+00  1.67E+02  6.07E+02
 
 TH 4
+       -5.84E+00  3.44E+02 -9.37E+01  7.04E+02
 
 TH 5
+        1.64E+00 -4.39E+02 -1.01E+03  2.12E+01  2.06E+03
 
 TH 6
+       -4.76E-02 -4.93E+00  1.60E+00 -2.05E+00 -1.15E+00  2.21E+02
 
 TH 7
+        5.92E-01  3.38E+01 -3.13E+00 -5.75E+00  5.60E+00 -7.89E-03  6.67E+00
 
 TH 8
+        1.09E-01  5.31E-01 -6.00E+01 -5.48E+00  1.86E+01  2.34E-01  3.71E-01  3.65E+01
 
 TH 9
+        2.74E+00 -4.00E+01  1.48E+01  1.20E+01 -1.92E+01 -7.96E-01 -2.28E+00  6.11E-01  2.51E+02
 
 TH10
+       -6.61E-02  2.01E+01 -2.35E+01 -1.70E+01 -8.57E+01  2.97E-01  9.73E-01  3.51E+01  4.41E+00  1.00E+02
 
 TH11
+       -8.58E+00 -3.65E+00 -2.39E+01 -8.60E+00  2.28E+00  2.59E+00  2.54E-01  1.49E+01  8.87E+00  2.05E+01  2.09E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.947
Stop Time:
Wed Sep 29 15:19:20 CDT 2021
