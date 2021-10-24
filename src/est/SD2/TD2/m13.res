Sat Oct 23 19:56:35 CDT 2021
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
$DATA ../../../../data/SD2/TD2/dat13.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      800
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

 TOT. NO. OF OBS RECS:      700
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m13.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2680.61326852347        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2788E+02 -3.3204E+01  7.3946E+01 -3.6586E+01  4.4882E+01  3.8768E+01  4.1433E+00 -4.1886E+02 -8.5640E+01  4.0748E+01
            -8.4705E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2796.67028843236        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  8.6708E-01  1.0623E+00  1.0508E+00  1.0034E+00  1.0979E+00  7.9903E-01  9.6790E-01  1.7047E+00  9.4937E-01  8.7428E-01
             1.0764E+00
 PARAMETER: -4.2628E-02  1.6042E-01  1.4952E-01  1.0341E-01  1.9339E-01 -1.2435E-01  6.7374E-02  6.3336E-01  4.8044E-02 -3.4354E-02
             1.7361E-01
 GRADIENT:  -4.7799E+02 -7.7479E+01 -3.4896E+01 -2.5371E+01  6.6437E+01 -1.7990E+02 -6.7277E+00 -1.7882E+02 -9.6152E+00  1.7928E+01
             4.2691E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2830.23622534709        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      336
 NPARAMETR:  9.1758E-01  9.6711E-01  1.7822E+00  1.1206E+00  1.1249E+00  7.0084E-01  1.2325E+00  2.2222E+00  7.3504E-01  8.6450E-01
             1.0997E+00
 PARAMETER:  1.3981E-02  6.6556E-02  6.7787E-01  2.1390E-01  2.1772E-01 -2.5547E-01  3.0907E-01  8.9848E-01 -2.0784E-01 -4.5602E-02
             1.9503E-01
 GRADIENT:  -3.6446E+02  2.4335E+01  4.0172E+01  2.3813E+00 -7.6051E+01 -2.3280E+02 -5.6200E+00 -1.3419E+02 -1.2932E+01  2.5333E+01
             8.0884E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2898.10399338196        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      513
 NPARAMETR:  9.5118E-01  9.6502E-01  2.3018E+00  1.1408E+00  1.3089E+00  8.5285E-01  1.4143E+00  3.2191E+00  6.0926E-01  7.0361E-01
             1.0374E+00
 PARAMETER:  4.9951E-02  6.4391E-02  9.3368E-01  2.3173E-01  3.6918E-01 -5.9174E-02  4.4665E-01  1.2691E+00 -3.9550E-01 -2.5153E-01
             1.3669E-01
 GRADIENT:  -1.3539E+02  2.3222E+01  9.8973E-01  5.0951E+01  7.7380E+01 -8.2804E+01 -1.5210E+00 -4.7942E+01 -1.3433E+01 -1.2563E+01
            -1.6491E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2906.70915804799        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.5178E-01  9.6383E-01  2.2956E+00  1.1385E+00  1.3044E+00  1.0824E+00  1.4109E+00  3.3879E+00  6.1195E-01  7.0363E-01
             1.0380E+00
 PARAMETER:  5.0579E-02  6.3157E-02  9.3100E-01  2.2971E-01  3.6571E-01  1.7923E-01  4.4422E-01  1.3202E+00 -3.9110E-01 -2.5150E-01
             1.3728E-01
 GRADIENT:  -8.2095E+01  1.9457E+01 -1.1895E+00  4.3292E+01  6.2762E+01  2.3642E+01 -1.6617E+00 -3.2821E+01 -1.2330E+01 -9.6512E+00
            -1.4296E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2912.31898433753        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      882            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9497E-01  9.4607E-01  2.2973E+00  1.1326E+00  1.2898E+00  1.0164E+00  1.4231E+00  3.8002E+00  6.8962E-01  7.0419E-01
             1.0489E+00
 PARAMETER:  9.4953E-02  4.4556E-02  9.3172E-01  2.2448E-01  3.5447E-01  1.1623E-01  4.5280E-01  1.4350E+00 -2.7161E-01 -2.5071E-01
             1.4771E-01
 GRADIENT:   3.9245E+02  4.2999E+01  2.0311E+01  2.1051E+02  9.1662E+01  5.0752E+01  5.0739E+01  4.9374E+01  1.3071E+01  2.5198E+00
             1.9521E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2912.74205072630        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1038
 NPARAMETR:  9.9416E-01  9.6071E-01  2.2973E+00  1.1326E+00  1.2897E+00  1.0053E+00  1.2747E+00  3.7988E+00  7.4476E-01  6.9063E-01
             1.0512E+00
 PARAMETER:  9.4144E-02  5.9922E-02  9.3175E-01  2.2448E-01  3.5442E-01  1.0530E-01  3.4275E-01  1.4347E+00 -1.9469E-01 -2.7015E-01
             1.4989E-01
 GRADIENT:  -3.1421E-01  1.7856E-01 -4.4320E+00  1.5862E+01  3.2888E+01  1.1413E-02  5.5210E-01  5.0866E-01  2.0426E-02 -7.3515E-03
            -3.7715E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2913.02522243724        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1216
 NPARAMETR:  9.9416E-01  9.5997E-01  2.2976E+00  1.1325E+00  1.2896E+00  1.0054E+00  1.1742E+00  3.7526E+00  8.0255E-01  7.0943E-01
             1.0499E+00
 PARAMETER:  9.4147E-02  5.9151E-02  9.3185E-01  2.2446E-01  3.5432E-01  1.0543E-01  2.6061E-01  1.4225E+00 -1.1997E-01 -2.4330E-01
             1.4867E-01
 GRADIENT:  -1.5101E-01 -2.6630E+00 -5.0215E+00  7.6038E+00  2.3768E+01  8.3731E-02 -6.9131E-01 -1.3465E+00  2.1301E-02 -2.1088E-02
             2.4322E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -2913.02534792361        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     1275
 NPARAMETR:  9.9418E-01  9.6048E-01  2.2976E+00  1.1325E+00  1.2896E+00  1.0053E+00  1.1719E+00  3.7562E+00  8.0386E-01  7.0987E-01
             1.0494E+00
 PARAMETER:  9.4163E-02  5.9680E-02  9.3186E-01  2.2446E-01  3.5431E-01  1.0526E-01  2.5861E-01  1.4234E+00 -1.1833E-01 -2.4268E-01
             1.4826E-01
 GRADIENT:  -1.0232E+00  1.2793E+02 -1.3169E+01 -6.0774E+01 -2.2018E+01 -1.0503E-01  5.1847E+01  1.3505E+01 -2.8440E-03  3.3233E-02
            -1.0302E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1275
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4278E-03 -1.2965E-03 -4.0958E-02 -6.3349E-03 -8.0336E-02
 SE:             2.9905E-02  2.0941E-02  2.4056E-02  2.3411E-02  1.8169E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6192E-01  9.5063E-01  8.8640E-02  7.8671E-01  9.8063E-06

 ETASHRINKSD(%)  1.0000E-10  2.9844E+01  1.9410E+01  2.1569E+01  3.9131E+01
 ETASHRINKVR(%)  1.0000E-10  5.0781E+01  3.5053E+01  3.8485E+01  6.2950E+01
 EBVSHRINKSD(%)  2.9890E-01  3.1597E+01  1.6282E+01  2.2120E+01  4.0378E+01
 EBVSHRINKVR(%)  5.9691E-01  5.3210E+01  2.9913E+01  3.9346E+01  6.4453E+01
 RELATIVEINF(%)  9.9389E+01  1.0724E+01  4.3786E+01  1.4776E+01  2.0269E+01
 EPSSHRINKSD(%)  2.6148E+01
 EPSSHRINKVR(%)  4.5458E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2913.0253479236062     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1626.5114014370645     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2913.025       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  9.60E-01  2.30E+00  1.13E+00  1.29E+00  1.01E+00  1.17E+00  3.76E+00  8.04E-01  7.10E-01  1.05E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,      121.945
Stop Time:
Sat Oct 23 19:56:53 CDT 2021
