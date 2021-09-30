Wed Sep 29 19:50:23 CDT 2021
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
$DATA ../../../../data/spa/D/dat21.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1208.87325936155        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9029E+01 -9.2624E+01 -4.8665E+01 -2.0969E+02  2.3314E+02 -8.5822E+02 -4.1131E+02 -1.8371E+01 -7.2501E+02 -3.9951E+02
            -3.4018E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1067.14196207688        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1448E+00  8.0976E-01  2.2386E+00  1.6308E+00  2.0204E+00  2.2101E+00  2.0612E+00  9.6065E-01  4.2046E+00  1.9507E+00
             4.4666E+00
 PARAMETER:  2.3521E-01 -1.1101E-01  9.0583E-01  5.8906E-01  8.0331E-01  8.9302E-01  8.2327E-01  5.9856E-02  1.5362E+00  7.6821E-01
             1.5966E+00
 GRADIENT:   1.9869E+01 -2.4517E+01 -1.5559E+01  4.5259E+01 -2.4002E+01  7.9456E+01  1.5153E+01  9.1305E-01  1.5640E+02  8.7940E+00
             1.3860E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1101.94230254659        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.4365E+00  7.6424E-01  1.6558E+00  1.5615E+00  1.1396E+01  2.3428E+00  4.5744E+00  8.7901E-01  2.4555E+00  1.0950E+01
             4.1428E+00
 PARAMETER:  4.6223E-01 -1.6888E-01  6.0431E-01  5.4562E-01  2.5332E+00  9.5133E-01  1.6205E+00 -2.8964E-02  9.9832E-01  2.4933E+00
             1.5214E+00
 GRADIENT:   2.1117E+02  1.6199E+01  1.3494E+00  5.4204E+01 -5.2888E-01  1.1153E+02  3.1214E+01 -1.9364E+00  6.6335E+01  7.1453E+00
             9.2299E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1151.83527135748        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0531E+00  2.6450E-01  1.6393E+00  1.7005E+00  9.0735E+00  1.8910E+00  6.2838E+00  4.6855E+00  1.6203E+00  1.0747E+01
             3.8406E+00
 PARAMETER:  1.5176E-01 -1.2299E+00  5.9429E-01  6.3095E-01  2.3054E+00  7.3710E-01  1.9380E+00  1.6445E+00  5.8264E-01  2.4746E+00
             1.4456E+00
 GRADIENT:   1.6596E+01  2.4256E+01  9.6422E+00  1.1741E+02 -7.9441E+00  7.9782E+01  6.3177E+00  1.1668E+02  4.2191E+01  2.4782E+01
             3.8600E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1153.91272920960        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:      360
 NPARAMETR:  9.9793E-01  1.9896E-01  1.5527E+00  1.7037E+00  8.2661E+00  1.7762E+00  6.4870E+00  4.0970E+00  1.4397E+00  1.0610E+01
             3.7989E+00
 PARAMETER:  9.7923E-02 -1.5146E+00  5.4000E-01  6.3280E-01  2.2122E+00  6.7449E-01  1.9698E+00  1.5103E+00  4.6445E-01  2.4618E+00
             1.4347E+00
 GRADIENT:  -3.0625E+01  2.3502E+01  1.7253E+01  1.2851E+02 -6.0571E+00  5.9547E+01  6.0209E+00  6.8345E+01  2.1490E+01  3.4356E+01
             2.3409E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1173.42847366368        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      485
 NPARAMETR:  9.9674E-01  1.8439E-01  5.6181E-01  1.5741E+00  9.1438E+00  1.7179E+00  6.3365E+00  2.6666E+00  1.4420E+00  6.4466E+00
             3.7304E+00
 PARAMETER:  9.6737E-02 -1.5907E+00 -4.7659E-01  5.5368E-01  2.3131E+00  6.4113E-01  1.9463E+00  1.0808E+00  4.6602E-01  1.9636E+00
             1.4165E+00
 GRADIENT:   1.2261E+01  3.2298E+01 -1.5666E+01  1.9491E+02 -5.7789E+00  5.4358E+01  9.0417E+00  2.0116E+01  2.9276E+01  7.7996E+00
             1.3491E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1179.28370196535        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      602
 NPARAMETR:  9.9595E-01  1.7545E-01  5.7348E-01  1.4953E+00  9.7505E+00  1.6826E+00  6.2202E+00  2.0721E+00  1.4431E+00  4.7734E+00
             3.6901E+00
 PARAMETER:  9.5947E-02 -1.6404E+00 -4.5604E-01  5.0233E-01  2.3773E+00  6.2034E-01  1.9278E+00  8.2857E-01  4.6679E-01  1.6631E+00
             1.4056E+00
 GRADIENT:  -2.1134E+01  1.6635E+01  7.4912E+00  5.8155E+01 -2.2902E+00 -1.0179E+01  9.7868E+00 -1.8638E+01  1.7437E+01  2.5150E-01
             4.5045E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1181.66555090339        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      788
 NPARAMETR:  9.9579E-01  1.7396E-01  5.4376E-01  1.4890E+00  9.8311E+00  1.7024E+00  6.1816E+00  2.4319E+00  1.3874E+00  4.2774E+00
             3.6906E+00
 PARAMETER:  9.5784E-02 -1.6489E+00 -5.0926E-01  4.9811E-01  2.3856E+00  6.3203E-01  1.9216E+00  9.8866E-01  4.2743E-01  1.5533E+00
             1.4058E+00
 GRADIENT:  -1.4857E+01  1.7772E+01 -5.3082E+00  8.1227E+01 -2.6594E+00  3.2855E+00  3.3484E+00 -7.1738E+00  1.8223E+01  4.6350E-02
             5.7535E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1197.43901614336        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      913
 NPARAMETR:  9.9052E-01  7.0971E-02  5.5364E-01  1.3602E+00  1.1703E+01  1.6472E+00  3.2892E+00  2.2980E+00  1.3937E+00  1.0000E-02
             3.4246E+00
 PARAMETER:  9.0471E-02 -2.5455E+00 -4.9123E-01  4.0762E-01  2.5598E+00  5.9911E-01  1.2907E+00  9.3203E-01  4.3194E-01 -8.0697E+00
             1.3310E+00
 GRADIENT:  -1.1520E-01 -1.6888E-01  3.5098E+01 -1.3554E+01  1.3848E-01 -5.4424E+00  4.7385E-01 -2.8908E+01  3.4744E+01  0.0000E+00
            -3.5790E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1198.43998713309        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1084
 NPARAMETR:  9.9250E-01  8.5790E-02  5.5192E-01  1.3905E+00  1.1327E+01  1.6635E+00  3.7399E+00  2.3251E+00  1.3936E+00  1.0000E-02
             3.5028E+00
 PARAMETER:  9.1472E-02 -2.3545E+00 -4.9459E-01  4.2944E-01  2.5258E+00  6.0859E-01  1.4198E+00  9.4416E-01  4.3166E-01 -5.9918E+00
             1.3526E+00
 GRADIENT:  -3.5999E+00  2.6491E+01 -2.1634E+02 -2.7260E+02 -4.7671E+01 -1.0200E+02  4.3912E+01  1.0230E+02 -2.4798E+02  0.0000E+00
            -6.1640E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1084
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.0073E-03 -1.0527E-02 -5.8228E-02 -2.5178E-02 -2.2561E-06
 SE:             2.9041E-02  1.5959E-03  1.7884E-02  2.2183E-02  1.7766E-06
 N:                     100         100         100         100         100

 P VAL.:         9.4489E-01  4.2391E-11  1.1307E-03  2.5637E-01  2.0413E-01

 ETASHRINKSD(%)  2.7106E+00  9.4654E+01  4.0086E+01  2.5684E+01  9.9994E+01
 ETASHRINKVR(%)  5.3478E+00  9.9714E+01  6.4103E+01  4.4771E+01  1.0000E+02
 EBVSHRINKSD(%)  2.7015E+00  9.2574E+01  5.3511E+01  1.3557E+01  9.9994E+01
 EBVSHRINKVR(%)  5.3300E+00  9.9449E+01  7.8387E+01  2.5275E+01  1.0000E+02
 RELATIVEINF(%)  1.9306E+01  5.2327E-02  1.8451E+00  1.9973E+00  5.9529E-08
 EPSSHRINKSD(%)  2.3856E+01
 EPSSHRINKVR(%)  4.2021E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1198.4399871330909     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -463.28916056935270     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.64
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1198.440       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.92E-01  8.59E-02  5.52E-01  1.39E+00  1.13E+01  1.66E+00  3.74E+00  2.33E+00  1.39E+00  1.00E-02  3.50E+00
 


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
+        3.95E+02
 
 TH 2
+       -4.63E+01  7.32E+04
 
 TH 3
+       -8.32E+01 -3.98E+01  4.08E+04
 
 TH 4
+       -3.74E+01  1.22E+02  1.84E+04  8.69E+03
 
 TH 5
+        1.29E-01  3.11E+00 -2.36E-01 -1.97E-01  3.68E+00
 
 TH 6
+       -4.98E+00 -1.21E+00  1.09E+04  4.99E+03 -7.05E-02  3.01E+03
 
 TH 7
+        6.37E-01  5.49E+01 -2.10E+03 -9.61E+02 -3.62E-01 -5.67E+02  1.09E+02
 
 TH 8
+        6.40E+00  1.15E+01 -2.24E+01  1.48E+01  3.79E-01  6.68E+00  4.72E+00  6.33E+02
 
 TH 9
+       -1.27E+01 -2.14E+01  1.84E+04  8.39E+03 -2.63E-01  4.96E+03 -9.54E+02  1.02E+01  8.40E+03
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.40E+00 -5.48E+01  1.57E+00 -7.95E+00  2.24E+01  1.66E+00 -2.12E+00 -4.16E+00  1.17E+00  0.00E+00  1.70E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       23.075
Stop Time:
Wed Sep 29 19:50:47 CDT 2021
