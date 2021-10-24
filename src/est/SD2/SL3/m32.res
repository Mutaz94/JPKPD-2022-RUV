Sat Oct 23 19:17:07 CDT 2021
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
$DATA ../../../../data/SD2/SL3/dat32.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      799
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

 TOT. NO. OF OBS RECS:      699
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
 RAW OUTPUT FILE (FILE): m32.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1510.30377619059        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6126E+02 -5.2115E+01 -2.6442E+01  1.2306E+02  1.2848E+02  8.4225E+01 -4.8945E+01 -4.7311E+01 -4.5736E+01 -3.5470E+01
            -2.8114E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2382.71802586800        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0649E+00  1.1034E+00  1.0249E+00  9.9085E-01  9.9016E-01  8.4289E-01  1.0202E+00  9.8758E-01  9.8738E-01  1.0227E+00
             2.1161E+00
 PARAMETER:  1.6290E-01  1.9842E-01  1.2461E-01  9.0807E-02  9.0111E-02 -7.0921E-02  1.2001E-01  8.7503E-02  8.7299E-02  1.2249E-01
             8.4958E-01
 GRADIENT:   2.0530E+02  9.7755E+00 -1.1799E+01  3.5036E+01  1.6996E+01 -1.0295E+01  6.3790E-01  3.1162E+00 -1.3513E+00 -6.8948E+00
            -6.2624E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2385.31476902499        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      224
 NPARAMETR:  1.0614E+00  1.4504E+00  9.7597E-01  8.2131E-01  1.1825E+00  8.9338E-01  8.6244E-01  9.7363E-01  1.2534E+00  1.1405E+00
             2.0978E+00
 PARAMETER:  1.5963E-01  4.7187E-01  7.5678E-02 -9.6858E-02  2.6762E-01 -1.2739E-02 -4.7986E-02  7.3272E-02  3.2587E-01  2.3142E-01
             8.4088E-01
 GRADIENT:   2.5073E+01  4.3181E+01 -7.4850E+00  5.8614E+01  2.4732E+01  2.5170E+00  4.0808E+00  2.9807E-01  4.4649E+00 -1.8402E+01
            -7.6314E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2390.52218602602        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      401
 NPARAMETR:  1.0501E+00  1.5085E+00  8.2817E-01  7.4487E-01  1.1455E+00  8.9038E-01  7.5652E-01  6.8085E-01  1.2895E+00  1.2811E+00
             2.2039E+00
 PARAMETER:  1.4893E-01  5.1113E-01 -8.8533E-02 -1.9455E-01  2.3588E-01 -1.6102E-02 -1.7902E-01 -2.8441E-01  3.5424E-01  3.4774E-01
             8.9021E-01
 GRADIENT:  -6.8736E+00  1.8780E+01 -2.2338E-01  1.5466E+01 -1.2354E+01  2.4512E+00  1.6283E+00  9.4313E-01  1.1857E+00  1.4244E+00
             8.7637E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2391.54332936920        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      577
 NPARAMETR:  1.0528E+00  1.7308E+00  6.9738E-01  5.8500E-01  1.2771E+00  8.8429E-01  6.7869E-01  4.2905E-01  1.5383E+00  1.3793E+00
             2.1808E+00
 PARAMETER:  1.5147E-01  6.4856E-01 -2.6043E-01 -4.3615E-01  3.4457E-01 -2.2968E-02 -2.8759E-01 -7.4619E-01  5.3066E-01  4.2160E-01
             8.7968E-01
 GRADIENT:   1.2469E+00  2.2786E+00 -1.0422E-01  1.9685E+00  1.0496E+00 -1.1012E-01 -6.8967E-01  9.9228E-02 -2.1813E-01 -3.2103E-01
            -2.1925E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2391.57280438128        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      758
 NPARAMETR:  1.0523E+00  1.7373E+00  6.8423E-01  5.7856E-01  1.2774E+00  8.8455E-01  6.8548E-01  2.9804E-01  1.5489E+00  1.3815E+00
             2.1835E+00
 PARAMETER:  1.5099E-01  6.5232E-01 -2.7946E-01 -4.4721E-01  3.4481E-01 -2.2681E-02 -2.7764E-01 -1.1105E+00  5.3757E-01  4.2314E-01
             8.8092E-01
 GRADIENT:  -8.5138E-02 -5.6554E-01  2.0751E-01  3.9059E-01  7.3879E-01 -5.2231E-03  1.7472E-01  2.3773E-02  1.4117E-01  1.8240E-01
             2.6629E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2391.58342085688        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      933
 NPARAMETR:  1.0517E+00  1.7397E+00  6.7044E-01  5.7543E-01  1.2727E+00  8.8435E-01  6.8358E-01  1.2746E-01  1.5535E+00  1.3784E+00
             2.1831E+00
 PARAMETER:  1.5043E-01  6.5372E-01 -2.9982E-01 -4.5263E-01  3.4112E-01 -2.2903E-02 -2.8040E-01 -1.9599E+00  5.4051E-01  4.2094E-01
             8.8075E-01
 GRADIENT:  -1.5728E+00 -1.8416E+00 -7.5007E-02 -4.0047E-01 -9.6507E-02 -1.0680E-01 -3.2826E-02  3.0616E-03 -2.3176E-03 -5.1637E-03
            -1.5969E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2391.58575940661        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1111            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0527E+00  1.7388E+00  6.6972E-01  5.7641E-01  1.2718E+00  8.8468E-01  6.8462E-01  4.8291E-02  1.5523E+00  1.3780E+00
             2.1834E+00
 PARAMETER:  1.5131E-01  6.5321E-01 -3.0090E-01 -4.5093E-01  3.4043E-01 -2.2525E-02 -2.7889E-01 -2.9305E+00  5.3975E-01  4.2063E-01
             8.8089E-01
 GRADIENT:   1.4770E+02  1.9338E+02  9.0077E-01  2.6851E+01  1.4086E+01  9.6430E+00  3.4800E+00  2.2844E-03  7.0947E+00  3.2052E+00
             1.2269E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2391.58587830624        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1290            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0527E+00  1.7387E+00  6.6953E-01  5.7645E-01  1.2717E+00  8.8468E-01  6.8464E-01  1.1608E-02  1.5521E+00  1.3779E+00
             2.1835E+00
 PARAMETER:  1.5132E-01  6.5317E-01 -3.0118E-01 -4.5087E-01  3.4033E-01 -2.2533E-02 -2.7887E-01 -4.3561E+00  5.3959E-01  4.2057E-01
             8.8091E-01
 GRADIENT:   1.4770E+02  1.9338E+02  9.3422E-01  2.6835E+01  1.4056E+01  9.6393E+00  3.4683E+00  1.3529E-04  7.0783E+00  3.1853E+00
             1.2289E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2391.58589150512        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1475
 NPARAMETR:  1.0527E+00  1.7387E+00  6.6949E-01  5.7649E-01  1.2715E+00  8.8468E-01  6.8461E-01  1.0000E-02  1.5519E+00  1.3778E+00
             2.1835E+00
 PARAMETER:  1.5132E-01  6.5311E-01 -3.0124E-01 -4.5080E-01  3.4023E-01 -2.2532E-02 -2.7890E-01 -4.5845E+00  5.3947E-01  4.2048E-01
             8.8092E-01
 GRADIENT:   8.3833E-01 -1.3504E+00 -2.5873E-02 -3.3680E-02  1.2081E-02  3.3307E-02  2.1600E-02  0.0000E+00  5.6546E-02  2.9581E-02
            -2.4215E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -2391.58589150512        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1497
 NPARAMETR:  1.0527E+00  1.7387E+00  6.6949E-01  5.7649E-01  1.2715E+00  8.8468E-01  6.8461E-01  1.0000E-02  1.5519E+00  1.3778E+00
             2.1835E+00
 PARAMETER:  1.5132E-01  6.5311E-01 -3.0124E-01 -4.5080E-01  3.4023E-01 -2.2532E-02 -2.7890E-01 -4.5845E+00  5.3947E-01  4.2048E-01
             8.8092E-01
 GRADIENT:   8.3833E-01 -1.3504E+00 -2.5873E-02 -3.3680E-02  1.2081E-02  3.3307E-02  2.1600E-02  0.0000E+00  5.6546E-02  2.9581E-02
            -2.4215E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1497
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1484E-03 -3.6211E-02 -1.2993E-04  2.3464E-02 -2.6458E-02
 SE:             2.9510E-02  2.0160E-02  8.6430E-05  2.3389E-02  2.4939E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6896E-01  7.2467E-02  1.3277E-01  3.1576E-01  2.8872E-01

 ETASHRINKSD(%)  1.1367E+00  3.2462E+01  9.9710E+01  2.1643E+01  1.6452E+01
 ETASHRINKVR(%)  2.2605E+00  5.4386E+01  9.9999E+01  3.8602E+01  3.0198E+01
 EBVSHRINKSD(%)  1.4444E+00  3.1564E+01  9.9738E+01  2.3076E+01  1.4759E+01
 EBVSHRINKVR(%)  2.8680E+00  5.3166E+01  9.9999E+01  4.0827E+01  2.7340E+01
 RELATIVEINF(%)  9.7066E+01  4.9594E+00  2.9102E-04  7.2015E+00  2.1617E+01
 EPSSHRINKSD(%)  2.1388E+01
 EPSSHRINKVR(%)  3.8202E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          699
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1284.6760694201323     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2391.5858915051176     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1106.9098220849853     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2391.586       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.74E+00  6.69E-01  5.76E-01  1.27E+00  8.85E-01  6.85E-01  1.00E-02  1.55E+00  1.38E+00  2.18E+00
 


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
 #CPUT: Total CPU Time in Seconds,      137.449
Stop Time:
Sat Oct 23 19:17:27 CDT 2021
