Sat Sep 18 13:12:01 CDT 2021
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
$DATA ../../../../data/spa/S2/dat2.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1677.25428482005        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.7397E+01 -5.7535E+01 -2.2153E+01 -6.0370E+01 -1.1757E+01  6.9697E+00  1.0700E+01  1.4298E+01  3.1735E+01  1.6344E+01
            -1.9818E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1682.32525538736        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8147E-01  1.1283E+00  1.2092E+00  9.8423E-01  1.1609E+00  9.6908E-01  9.0077E-01  8.6623E-01  7.4132E-01  9.2290E-01
             1.1367E+00
 PARAMETER:  8.1292E-02  2.2072E-01  2.8993E-01  8.4103E-02  2.4920E-01  6.8587E-02 -4.5053E-03 -4.3599E-02 -1.9932E-01  1.9767E-02
             2.2812E-01
 GRADIENT:   3.4709E+01  1.9771E+01  8.8923E+00  1.3142E+01  2.2926E+01 -4.1122E+00 -1.0856E+01 -3.5300E+00 -2.2642E+01 -2.1618E+01
             1.2635E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1685.90388857006        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.8575E-01  9.2642E-01  1.0231E+00  1.1110E+00  1.0038E+00  9.7036E-01  9.7392E-01  2.7780E-01  8.1932E-01  9.5493E-01
             1.0955E+00
 PARAMETER:  8.5648E-02  2.3574E-02  1.2283E-01  2.0528E-01  1.0378E-01  6.9911E-02  7.3572E-02 -1.1808E+00 -9.9275E-02  5.3883E-02
             1.9124E-01
 GRADIENT:   4.7519E+01  6.9503E-01 -2.0016E+01  5.2626E+01  4.3175E+01 -4.1540E+00 -4.7552E+00 -2.7174E-01  2.1492E+00 -8.5860E+00
             1.0175E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1687.48303698996        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      271
 NPARAMETR:  9.7767E-01  1.0093E+00  9.5165E-01  1.0386E+00  9.8347E-01  9.8748E-01  1.0246E+00  3.3992E-01  8.1335E-01  9.4751E-01
             1.0623E+00
 PARAMETER:  7.7413E-02  1.0927E-01  5.0440E-02  1.3784E-01  8.3328E-02  8.7402E-02  1.2428E-01 -9.7903E-01 -1.0659E-01  4.6078E-02
             1.6047E-01
 GRADIENT:  -5.4883E+00 -2.9808E+00 -1.7847E+00 -2.3712E+00  2.8334E+00 -4.0260E-01  2.7768E-01  1.9218E-01 -1.8803E-01  4.7287E-01
             5.8310E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1687.56793645241        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  9.8015E-01  1.0929E+00  8.8835E-01  9.8731E-01  9.8726E-01  9.8850E-01  9.5864E-01  1.8452E-01  8.4583E-01  9.3945E-01
             1.0612E+00
 PARAMETER:  7.9951E-02  1.8885E-01 -1.8388E-02  8.7229E-02  8.7175E-02  8.8430E-02  5.7760E-02 -1.5900E+00 -6.7434E-02  3.7537E-02
             1.5944E-01
 GRADIENT:  -1.6783E+00  2.1227E+00 -1.0892E+00  4.0751E+00  5.8267E-02 -3.5691E-01 -4.1553E-01  6.7934E-02 -4.8669E-01  2.8734E-01
             1.2229E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1687.60016465944        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      623
 NPARAMETR:  9.8130E-01  1.1909E+00  8.5452E-01  9.2302E-01  1.0209E+00  9.8986E-01  8.9465E-01  8.1949E-02  8.9523E-01  9.5120E-01
             1.0617E+00
 PARAMETER:  8.1124E-02  2.7470E-01 -5.7211E-02  1.9901E-02  1.2069E-01  8.9809E-02 -1.1328E-02 -2.4017E+00 -1.0671E-02  4.9968E-02
             1.5983E-01
 GRADIENT:  -3.6863E-02 -2.3158E-02 -4.0143E-03 -1.3291E-01  2.4520E-01  1.8357E-04  1.0408E-02  1.0488E-02 -5.1226E-02 -9.9153E-02
            -9.5787E-03

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1687.60673796901        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  9.8117E-01  1.1537E+00  8.6410E-01  9.4683E-01  1.0059E+00  9.8974E-01  9.1938E-01  2.5902E-02  8.7621E-01  9.4530E-01
             1.0615E+00
 PARAMETER:  8.0991E-02  2.4297E-01 -4.6069E-02  4.5369E-02  1.0587E-01  8.9684E-02  1.5946E-02 -3.5535E+00 -3.2146E-02  4.3750E-02
             1.5970E-01
 GRADIENT:  -2.5453E-02  5.8626E-01  1.3856E-01  5.5373E-01 -2.3284E-01  5.3055E-03 -2.8364E-02  9.9656E-04 -5.5522E-02 -3.4135E-02
            -2.2485E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1687.60733852286        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      976
 NPARAMETR:  9.8119E-01  1.1563E+00  8.6290E-01  9.4482E-01  1.0068E+00  9.8972E-01  9.1769E-01  1.0000E-02  8.7781E-01  9.4566E-01
             1.0615E+00
 PARAMETER:  8.1013E-02  2.4523E-01 -4.7455E-02  4.3236E-02  1.0682E-01  8.9668E-02  1.4100E-02 -4.8408E+00 -3.0323E-02  4.4127E-02
             1.5966E-01
 GRADIENT:   7.8515E-04  2.6013E-02  1.8021E-02  1.7703E-02 -9.2453E-03 -5.5020E-03 -1.5169E-03  0.0000E+00 -1.9225E-03 -1.0368E-02
            -1.2752E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1687.60733938144        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1033
 NPARAMETR:  9.8119E-01  1.1564E+00  8.6280E-01  9.4473E-01  1.0069E+00  9.8974E-01  9.1762E-01  1.0000E-02  8.7787E-01  9.4568E-01
             1.0615E+00
 PARAMETER:  8.1013E-02  2.4533E-01 -4.7569E-02  4.3146E-02  1.0683E-01  8.9683E-02  1.4030E-02 -4.8403E+00 -3.0260E-02  4.4154E-02
             1.5968E-01
 GRADIENT:  -5.7747E-04  1.6851E-03 -2.5237E-04  3.0788E-03  6.0329E-04 -1.7837E-04  5.2878E-04  0.0000E+00  7.4722E-04  1.2246E-04
             9.5088E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1033
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.1172E-04 -1.1354E-02 -3.3172E-04  1.8836E-03 -2.4006E-02
 SE:             2.9795E-02  2.0292E-02  1.5481E-04  2.3565E-02  2.3237E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9433E-01  5.7579E-01  3.2138E-02  9.3629E-01  3.0156E-01

 ETASHRINKSD(%)  1.8291E-01  3.2020E+01  9.9481E+01  2.1053E+01  2.2153E+01
 ETASHRINKVR(%)  3.6549E-01  5.3788E+01  9.9997E+01  3.7673E+01  3.9399E+01
 EBVSHRINKSD(%)  4.8738E-01  3.1516E+01  9.9518E+01  2.1593E+01  2.0624E+01
 EBVSHRINKVR(%)  9.7239E-01  5.3100E+01  9.9998E+01  3.8523E+01  3.6995E+01
 RELATIVEINF(%)  9.8574E+01  1.3426E+00  2.1525E-04  2.0450E+00  6.2679E+00
 EPSSHRINKSD(%)  4.1932E+01
 EPSSHRINKVR(%)  6.6281E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1687.6073393814377     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -952.45651281769949     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.66
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1687.607       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.16E+00  8.63E-01  9.45E-01  1.01E+00  9.90E-01  9.18E-01  1.00E-02  8.78E-01  9.46E-01  1.06E+00
 


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
+        1.17E+03
 
 TH 2
+       -9.48E+00  4.81E+02
 
 TH 3
+        1.47E+01  1.73E+02  3.81E+02
 
 TH 4
+       -1.66E+01  4.73E+02 -1.98E+02  9.97E+02
 
 TH 5
+       -1.78E+00 -2.86E+02 -4.51E+02  2.11E+02  7.67E+02
 
 TH 6
+        8.44E-01 -3.04E+00  8.07E+00 -6.70E+00  7.33E-02  1.96E+02
 
 TH 7
+        2.23E+00  2.09E+01  5.19E+00 -1.05E+01 -1.25E+01 -1.46E+00  5.31E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.21E+00 -2.22E+01 -2.16E+01  3.31E+01  8.40E+00 -5.99E-01  3.01E+01  0.00E+00  1.03E+02
 
 TH10
+       -7.04E-01 -7.73E+00 -3.88E+01 -1.45E+01 -5.87E+01 -1.02E+00  1.29E+01  0.00E+00  2.59E+00  9.29E+01
 
 TH11
+       -7.82E+00 -2.04E+01 -3.76E+01 -9.00E-01  5.59E+00  3.56E+00  6.63E+00  0.00E+00  1.50E+01  2.15E+01  1.99E+02
 
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
 #CPUT: Total CPU Time in Seconds,       17.515
Stop Time:
Sat Sep 18 13:12:20 CDT 2021
