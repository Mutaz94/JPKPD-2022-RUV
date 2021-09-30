Wed Sep 29 16:34:14 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat40.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1618.06560095993        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3782E+02  1.4716E+01 -4.0951E+01  1.0336E+02  1.0125E+02  8.9406E+01 -5.7976E+00  1.1853E+00 -3.3491E+00 -1.0262E+01
            -1.1773E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1628.47563954894        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  1.0372E+00  1.0325E+00  1.0445E+00  9.7112E-01  9.3796E-01  8.2204E-01  1.0351E+00  9.9689E-01  1.0469E+00  1.0189E+00
             1.5174E+00
 PARAMETER:  1.3650E-01  1.3203E-01  1.4356E-01  7.0699E-02  3.5948E-02 -9.5961E-02  1.3453E-01  9.6881E-02  1.4586E-01  1.1873E-01
             5.1697E-01
 GRADIENT:   8.5203E+01 -2.7844E-01  5.5105E+00 -1.0779E+01 -3.9350E+01 -1.8078E+01  2.7558E+00  4.4337E+00  5.3852E+00  1.4174E+01
             7.2306E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1631.65793381938        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  1.0237E+00  1.2628E+00  1.2033E+00  8.5719E-01  1.1038E+00  8.1732E-01  8.6607E-01  8.4702E-01  1.1878E+00  1.2370E+00
             1.4548E+00
 PARAMETER:  1.2342E-01  3.3332E-01  2.8503E-01 -5.4099E-02  1.9871E-01 -1.0172E-01 -4.3784E-02 -6.6034E-02  2.7210E-01  3.1268E-01
             4.7485E-01
 GRADIENT:   4.4473E+01  3.3773E+01  7.9276E+00  2.2368E+01 -2.5272E+01 -1.9854E+01  5.6371E+00  6.9723E-01 -6.5905E-01  1.6821E+01
             4.9173E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1638.74219084358        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  1.0104E+00  1.1653E+00  8.3275E-01  8.8117E-01  9.4195E-01  8.5796E-01  1.0455E+00  6.7117E-01  1.0695E+00  9.6234E-01
             1.2468E+00
 PARAMETER:  1.1031E-01  2.5301E-01 -8.3019E-02 -2.6501E-02  4.0192E-02 -5.3197E-02  1.4450E-01 -2.9873E-01  1.6718E-01  6.1612E-02
             3.2055E-01
 GRADIENT:   2.0734E+00 -4.5620E+00 -5.7315E+00  4.1374E+00  6.1294E+00 -1.1433E+00  1.1958E+00  8.7628E-01  1.3222E+00  3.2031E+00
            -5.2761E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1639.12618579948        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  1.0105E+00  1.3644E+00  7.0854E-01  7.5540E-01  9.7216E-01  8.6229E-01  9.0346E-01  4.1034E-01  1.1971E+00  9.5541E-01
             1.2533E+00
 PARAMETER:  1.1042E-01  4.1068E-01 -2.4455E-01 -1.8051E-01  7.1768E-02 -4.8163E-02 -1.5213E-03 -7.9076E-01  2.7989E-01  5.4389E-02
             3.2580E-01
 GRADIENT:  -2.4456E-01  7.8983E+00  1.8308E+00  6.9603E+00 -5.2430E+00  5.0223E-01 -2.6920E+00  4.4153E-03 -4.4772E-01 -6.0225E-01
            -7.5672E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1639.32326257783        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  1.0109E+00  1.5414E+00  5.9927E-01  6.3706E-01  1.0192E+00  8.6086E-01  8.5901E-01  1.6844E-01  1.3130E+00  9.7531E-01
             1.2536E+00
 PARAMETER:  1.1084E-01  5.3269E-01 -4.1204E-01 -3.5089E-01  1.1903E-01 -4.9823E-02 -5.1975E-02 -1.6812E+00  3.7230E-01  7.4997E-02
             3.2603E-01
 GRADIENT:  -1.8530E-01  7.1591E+00  7.4221E-01  3.0483E+00 -2.0139E+00 -2.8775E-01  5.5990E-01  1.9908E-02 -2.1158E+00  4.9075E-01
            -4.0306E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1639.36587013238        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  1.0111E+00  1.5842E+00  5.7517E-01  6.0454E-01  1.0335E+00  8.6146E-01  8.3318E-01  4.7588E-02  1.3840E+00  9.7659E-01
             1.2547E+00
 PARAMETER:  1.1099E-01  5.6010E-01 -4.5309E-01 -4.0328E-01  1.3293E-01 -4.9126E-02 -8.2505E-02 -2.9452E+00  4.2500E-01  7.6317E-02
             3.2693E-01
 GRADIENT:   4.4127E-01  2.6505E-03  6.0047E-01  3.5508E-01 -1.4583E+00  2.8771E-02  4.9568E-02  1.7499E-03  1.0253E-01 -1.6929E-01
            -1.8834E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1639.36950281120        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1241             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0115E+00  1.5829E+00  5.7426E-01  6.0418E-01  1.0342E+00  8.6157E-01  8.3284E-01  1.0000E-02  1.3838E+00  9.7783E-01
             1.2548E+00
 PARAMETER:  1.1143E-01  5.5926E-01 -4.5467E-01 -4.0388E-01  1.3362E-01 -4.9003E-02 -8.2916E-02 -5.6484E+00  4.2480E-01  7.7576E-02
             3.2697E-01
 GRADIENT:   3.2737E+02  3.5251E+02  2.3454E+00  7.3254E+01  6.9657E+00  2.4780E+01  3.7544E+00  0.0000E+00  1.1776E+01  5.1811E-01
             2.5740E+00

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1639.36950317496        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:     1324
 NPARAMETR:  1.0110E+00  1.5815E+00  5.7579E-01  6.0397E-01  1.0328E+00  8.6150E-01  8.3365E-01  1.0000E-02  1.3798E+00  9.7782E-01
             1.2553E+00
 PARAMETER:  1.1141E-01  5.5927E-01 -4.5467E-01 -4.0392E-01  1.3363E-01 -4.9006E-02 -8.2905E-02 -5.5719E+00  4.2478E-01  7.7590E-02
             3.2698E-01
 GRADIENT:   2.6920E-01  3.6270E-01 -1.0987E-01  4.2079E-02  3.6101E-01  5.7839E-03 -2.7304E-02  0.0000E+00  8.7885E-02  6.3475E-04
            -3.3319E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1324
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.1719E-04 -2.6752E-02 -2.4974E-04  1.9491E-02 -3.2620E-02
 SE:             2.9698E-02  2.2629E-02  1.0045E-04  2.2846E-02  2.2240E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9685E-01  2.3714E-01  1.2914E-02  3.9358E-01  1.4245E-01

 ETASHRINKSD(%)  5.0844E-01  2.4189E+01  9.9663E+01  2.3461E+01  2.5493E+01
 ETASHRINKVR(%)  1.0143E+00  4.2527E+01  9.9999E+01  4.1418E+01  4.4487E+01
 EBVSHRINKSD(%)  8.7647E-01  2.3514E+01  9.9707E+01  2.4954E+01  2.4429E+01
 EBVSHRINKVR(%)  1.7453E+00  4.1499E+01  9.9999E+01  4.3681E+01  4.2890E+01
 RELATIVEINF(%)  9.8105E+01  3.2608E+00  8.6617E-05  3.3082E+00  1.0240E+01
 EPSSHRINKSD(%)  4.2508E+01
 EPSSHRINKVR(%)  6.6946E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1639.3695031749646     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -904.21867661122644     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1639.370       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.58E+00  5.74E-01  6.04E-01  1.03E+00  8.62E-01  8.33E-01  1.00E-02  1.38E+00  9.78E-01  1.25E+00
 


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
+        1.44E+03
 
 TH 2
+       -1.16E+01  3.93E+02
 
 TH 3
+        1.48E+01  1.43E+02  3.08E+02
 
 TH 4
+       -2.86E+01  3.36E+02 -2.42E+02  9.13E+02
 
 TH 5
+       -7.25E+00 -2.17E+02 -3.40E+02  2.46E+02  6.20E+02
 
 TH 6
+       -9.86E-01 -2.34E+00  1.99E+00 -4.97E+00 -1.92E+00  2.59E+02
 
 TH 7
+        9.15E-01  9.78E+00  2.46E+00 -1.76E+01 -1.55E+01 -1.23E-01  9.97E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.17E+00 -1.71E+01 -2.56E+01  4.82E+01  2.28E+00  5.19E-01  1.99E+01  0.00E+00  4.01E+01
 
 TH10
+       -5.47E-01 -1.23E+01 -3.17E+01 -5.80E+00 -5.54E+01  4.40E-01  1.45E+01  0.00E+00  5.10E+00  7.60E+01
 
 TH11
+       -1.09E+01 -1.78E+01 -2.83E+01  1.95E+00 -2.31E+00  2.57E+00  1.09E+01  0.00E+00  5.54E+00  1.67E+01  1.39E+02
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       23.616
Stop Time:
Wed Sep 29 16:34:39 CDT 2021
