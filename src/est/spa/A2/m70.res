Sat Sep 25 08:49:11 CDT 2021
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
$DATA ../../../../data/spa/A2/dat70.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m70.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1156.94124576461        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.9804E+01  2.2412E+01  3.5352E+01 -2.5223E+01  4.8593E+01  1.5568E+01 -2.2666E+01 -5.9886E+00 -4.7827E+01 -2.3200E+01
            -9.5658E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1468.52927316997        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0296E+00  9.4945E-01  1.0819E+00  1.0393E+00  1.0001E+00  8.8483E-01  9.8474E-01  9.0641E-01  1.0594E+00  7.7907E-01
             2.3102E+00
 PARAMETER:  1.2918E-01  4.8129E-02  1.7869E-01  1.3852E-01  1.0009E-01 -2.2356E-02  8.4624E-02  1.7415E-03  1.5768E-01 -1.4965E-01
             9.3733E-01
 GRADIENT:  -1.7080E+01 -1.1813E+01  8.4151E-01 -2.7988E+01  1.0241E+01 -2.1170E+01  2.8779E+00  2.9566E+00 -1.4376E-01  5.7670E+00
            -1.6215E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1472.97838337547        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0320E+00  7.7751E-01  8.3631E-01  1.1523E+00  7.9344E-01  9.4328E-01  9.1935E-01  6.8579E-01  1.0751E+00  5.1831E-01
             2.2918E+00
 PARAMETER:  1.3149E-01 -1.5166E-01 -7.8754E-02  2.4177E-01 -1.3138E-01  4.1605E-02  1.5910E-02 -2.7719E-01  1.7243E-01 -5.5718E-01
             9.2935E-01
 GRADIENT:  -1.1978E+01  3.2961E+00 -4.9633E+00  1.2019E+01  6.1779E+00  1.3218E+00 -1.7986E+00  1.5657E+00  9.8217E+00  5.6067E-01
            -1.9386E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1478.12225823814        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  1.0446E+00  5.6489E-01  6.0279E-01  1.2651E+00  5.8472E-01  9.4201E-01  1.8121E+00  2.1260E-01  7.9529E-01  3.4939E-01
             2.3949E+00
 PARAMETER:  1.4362E-01 -4.7112E-01 -4.0619E-01  3.3516E-01 -4.3662E-01  4.0261E-02  6.9450E-01 -1.4483E+00 -1.2905E-01 -9.5157E-01
             9.7333E-01
 GRADIENT:   4.8092E+00  1.8659E+01 -1.1576E+01  8.0352E+01  1.2677E+01  2.1754E-01  5.3334E-01  2.9278E-01 -8.6566E+00  6.9523E-01
             1.3974E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1481.82707390094        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0414E+00  4.2929E-01  3.3894E-01  1.1750E+00  3.7487E-01  9.4189E-01  1.9827E+00  2.7217E-02  8.3280E-01  2.7549E-01
             2.1719E+00
 PARAMETER:  1.4061E-01 -7.4562E-01 -9.8193E-01  2.6128E-01 -8.8117E-01  4.0135E-02  7.8445E-01 -3.5039E+00 -8.2963E-02 -1.1892E+00
             8.7561E-01
 GRADIENT:  -1.3386E+00  3.4590E+00 -5.5861E+00  1.4722E+01  7.4089E+00 -5.7006E+00  4.5746E+00  6.0425E-03  2.4439E-01 -1.9060E+00
             3.4268E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1482.55861608777        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  1.0460E+00  4.2662E-01  3.7161E-01  1.1916E+00  3.9429E-01  9.5616E-01  1.9457E+00  2.2665E-02  8.3276E-01  3.6950E-01
             2.1378E+00
 PARAMETER:  1.4500E-01 -7.5187E-01 -8.8991E-01  2.7530E-01 -8.3068E-01  5.5171E-02  7.6562E-01 -3.6869E+00 -8.3015E-02 -8.9562E-01
             8.5977E-01
 GRADIENT:   3.6970E-02  1.1108E+00  3.8570E+00 -8.9443E+00 -4.8597E+00 -2.1967E-01 -8.8844E-01  8.5546E-03 -1.3369E-01 -1.1477E-01
            -3.6864E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1482.72437946590        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      616
 NPARAMETR:  1.0458E+00  3.9072E-01  3.9224E-01  1.2236E+00  4.0065E-01  9.5358E-01  2.1209E+00  1.6354E-02  8.2320E-01  3.9673E-01
             2.1636E+00
 PARAMETER:  1.4476E-01 -8.3976E-01 -8.3587E-01  3.0176E-01 -8.1466E-01  5.2470E-02  8.5184E-01 -4.0133E+00 -9.4560E-02 -8.2451E-01
             8.7179E-01
 GRADIENT:   7.4467E-01  5.7398E-01 -1.0830E-01  1.0158E+00 -7.9362E-01 -4.6280E-01  2.9438E-01  4.4422E-03 -1.1954E-01 -3.5713E-02
             1.3480E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1482.72965573859        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:      790
 NPARAMETR:  1.0452E+00  3.8689E-01  3.9561E-01  1.2265E+00  4.0218E-01  9.5465E-01  2.1347E+00  1.5526E-02  8.2334E-01  4.0286E-01
             2.1578E+00
 PARAMETER:  1.4421E-01 -8.4962E-01 -8.2733E-01  3.0420E-01 -8.1087E-01  5.3588E-02  8.5833E-01 -4.0653E+00 -9.4384E-02 -8.0916E-01
             8.6908E-01
 GRADIENT:  -1.3217E-02  3.9598E-02  1.0951E-01  1.3881E-01 -1.9256E-01 -1.0633E-03  6.8570E-03  4.0205E-03 -1.3793E-02 -1.3594E-02
            -5.4511E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1482.73084642306        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      934
 NPARAMETR:  1.0452E+00  3.8686E-01  3.9560E-01  1.2265E+00  4.0219E-01  9.5465E-01  2.1346E+00  1.0000E-02  8.2338E-01  4.0300E-01
             2.1580E+00
 PARAMETER:  1.4422E-01 -8.4970E-01 -8.2735E-01  3.0414E-01 -8.1083E-01  5.3587E-02  8.5826E-01 -6.5387E+00 -9.4339E-02 -8.0881E-01
             8.6917E-01
 GRADIENT:   3.7161E-03 -2.5669E-02 -1.9950E-02 -1.5461E-02  3.0537E-02  1.5623E-03  8.4226E-04  0.0000E+00  2.4705E-03  1.2032E-03
             1.2114E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1482.73085643945        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1097
 NPARAMETR:  1.0452E+00  3.8725E-01  3.9556E-01  1.2263E+00  4.0224E-01  9.5465E-01  2.1330E+00  1.0000E-02  8.2343E-01  4.0278E-01
             2.1579E+00
 PARAMETER:  1.4423E-01 -8.4869E-01 -8.2744E-01  3.0400E-01 -8.1070E-01  5.3593E-02  8.5751E-01 -6.4856E+00 -9.4282E-02 -8.0936E-01
             8.6915E-01
 GRADIENT:   5.4436E-03 -2.2366E-03  1.6603E-03  3.7961E-03 -2.5043E-03  6.7592E-04 -5.2671E-03  0.0000E+00  3.0887E-04 -3.0981E-05
            -1.1379E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1097
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4712E-03  3.9843E-02 -5.1941E-04 -2.4375E-02  1.7942E-02
 SE:             2.9364E-02  2.0305E-02  2.2811E-04  2.5642E-02  1.3760E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6004E-01  4.9730E-02  2.2788E-02  3.4181E-01  1.9228E-01

 ETASHRINKSD(%)  1.6268E+00  3.1977E+01  9.9236E+01  1.4096E+01  5.3901E+01
 ETASHRINKVR(%)  3.2271E+00  5.3729E+01  9.9994E+01  2.6205E+01  7.8749E+01
 EBVSHRINKSD(%)  1.9008E+00  3.4647E+01  9.9148E+01  1.3314E+01  5.1383E+01
 EBVSHRINKVR(%)  3.7656E+00  5.7290E+01  9.9993E+01  2.4855E+01  7.6364E+01
 RELATIVEINF(%)  9.5650E+01  9.9260E+00  2.3360E-04  2.7197E+01  6.8160E-01
 EPSSHRINKSD(%)  3.5101E+01
 EPSSHRINKVR(%)  5.7882E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1482.7308564394505     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -747.58002987571228     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.06
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1482.731       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  3.87E-01  3.96E-01  1.23E+00  4.02E-01  9.55E-01  2.13E+00  1.00E-02  8.23E-01  4.03E-01  2.16E+00
 


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
+        1.07E+03
 
 TH 2
+       -4.14E+01  7.50E+02
 
 TH 3
+        1.06E+01  1.08E+03  5.66E+03
 
 TH 4
+       -3.69E+01  2.24E+02 -7.25E+02  8.29E+02
 
 TH 5
+        6.15E+01 -1.90E+03 -7.36E+03  4.01E+02  1.04E+04
 
 TH 6
+       -3.57E-01 -6.87E+00  1.06E+01 -8.59E+00  9.49E+00  2.01E+02
 
 TH 7
+        1.74E+00  5.03E+01 -1.69E+01 -5.04E+00 -5.33E+00  4.47E-01  1.31E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.00E+00 -1.44E+01 -3.11E+01 -1.53E+01  7.83E+01  4.34E+00  7.52E+00  0.00E+00  1.76E+02
 
 TH10
+       -3.62E+00  1.05E+01 -2.49E+02 -2.59E+01  2.82E+02 -1.77E+00  7.00E+00  0.00E+00 -4.22E+00  1.00E+02
 
 TH11
+       -1.27E+01 -3.86E+00 -8.05E+01 -1.30E+01  6.25E+01  3.85E+00  3.08E+00  0.00E+00  1.09E+01  2.60E+01  5.83E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.683
Stop Time:
Sat Sep 25 08:49:32 CDT 2021
