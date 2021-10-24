Sun Oct 24 02:17:21 CDT 2021
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
$DATA ../../../../data/SD4/A2/dat64.csv ignore=@
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
Current Date:       24 OCT 2021
Days until program expires : 175
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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -720.006069378636        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1253E+02  3.8752E+01  4.3362E+01  4.0345E+01  1.7780E+02  2.1675E+01 -4.8291E+01 -3.0554E+01 -5.2679E+01 -1.1481E+02
            -1.6143E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1294.06072510276        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1552E+00  9.4937E-01  8.7308E-01  1.0962E+00  8.2251E-01  1.3142E+00  1.1116E+00  1.0300E+00  1.1932E+00  1.1379E+00
             2.1309E+00
 PARAMETER:  2.4431E-01  4.8046E-02 -3.5730E-02  1.9183E-01 -9.5395E-02  3.7324E-01  2.0578E-01  1.2959E-01  2.7666E-01  2.2915E-01
             8.5656E-01
 GRADIENT:   3.8789E+02  3.5073E+01  1.4137E+01  6.1905E+01 -4.7287E-01  8.2813E+01 -1.8567E+00  2.2273E+00  1.8680E+01  3.2977E+00
            -1.8511E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1301.85733352567        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      198
 NPARAMETR:  1.1444E+00  7.3546E-01  5.7840E-01  1.2272E+00  5.9008E-01  1.2489E+00  1.5573E+00  8.2679E-01  1.0170E+00  8.6096E-01
             2.0854E+00
 PARAMETER:  2.3492E-01 -2.0726E-01 -4.4749E-01  3.0477E-01 -4.2749E-01  3.2229E-01  5.4293E-01 -9.0204E-02  1.1684E-01 -4.9705E-02
             8.3496E-01
 GRADIENT:   1.4967E+02  3.5654E+01 -2.6817E+01  9.9353E+01  4.1796E+01  4.1696E+01  6.5013E+00  8.7781E+00  8.4990E-01  3.7792E+00
            -1.9779E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1341.67703412472        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0272E+00  6.2074E-01  5.6270E-01  1.2217E+00  5.4540E-01  9.6929E-01  1.8896E+00  4.5108E-01  7.8305E-01  7.8233E-01
             2.4890E+00
 PARAMETER:  1.2684E-01 -3.7684E-01 -4.7501E-01  3.0026E-01 -5.0624E-01  6.8805E-02  7.3635E-01 -6.9612E-01 -1.4456E-01 -1.4548E-01
             1.0119E+00
 GRADIENT:  -1.3068E+01  2.2256E+01 -1.1611E+01  2.8448E+01  2.6234E+01 -2.3690E+01  1.0857E+01  4.9225E-01 -3.0871E+01  3.6446E+00
            -7.6481E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1356.24687214843        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0373E+00  3.6029E-01  4.0703E-01  1.2674E+00  3.7799E-01  1.0338E+00  2.3095E+00  7.1782E-01  8.2318E-01  4.9854E-01
             2.7242E+00
 PARAMETER:  1.3658E-01 -9.2084E-01 -7.9886E-01  3.3697E-01 -8.7289E-01  1.3324E-01  9.3703E-01 -2.3153E-01 -9.4580E-02 -5.9607E-01
             1.1022E+00
 GRADIENT:  -1.0726E+00  5.4544E+00 -1.6435E-01  8.9949E+00 -5.8097E-01  5.4763E-03  6.4590E+00  4.0413E+00 -1.2676E+01  3.7027E+00
            -1.5211E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1357.97868576447        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0456E+00  3.9762E-01  3.5066E-01  1.2127E+00  3.5017E-01  1.0283E+00  1.9486E+00  5.6285E-01  9.0843E-01  4.8124E-01
             2.7087E+00
 PARAMETER:  1.4459E-01 -8.2227E-01 -9.4795E-01  2.9289E-01 -9.4932E-01  1.2790E-01  7.6713E-01 -4.7474E-01  3.9650E-03 -6.3138E-01
             1.0965E+00
 GRADIENT:   1.1340E+01  2.7799E+00  5.1683E+00  1.1200E+00 -8.9526E+00 -3.9919E+00  3.9179E+00 -1.5655E-01  1.3447E+00  1.4082E+00
            -8.0969E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1362.61306815134        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      904
 NPARAMETR:  1.0285E+00  6.0035E-01  2.0494E-01  1.0041E+00  3.0225E-01  1.0773E+00  1.0839E+00  1.0575E+00  1.0089E+00  2.7470E-01
             2.6209E+00
 PARAMETER:  1.2808E-01 -4.1024E-01 -1.4850E+00  1.0410E-01 -1.0965E+00  1.7449E-01  1.8059E-01  1.5592E-01  1.0888E-01 -1.1921E+00
             1.0635E+00
 GRADIENT:  -8.4887E+00  2.4855E+01  2.7452E+01  1.5247E+01 -6.4856E+01  5.4592E+00 -2.5154E+00 -6.4717E+00 -4.1023E+00  1.8516E+00
             2.2747E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1375.04428896116        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1081
 NPARAMETR:  1.0136E+00  6.4066E-01  1.8533E-01  9.5508E-01  3.0556E-01  1.0547E+00  9.2271E-01  1.6356E+00  1.1364E+00  2.4177E-01
             2.2427E+00
 PARAMETER:  1.1348E-01 -3.4526E-01 -1.5856E+00  5.4036E-02 -1.0856E+00  1.5324E-01  1.9555E-02  5.9203E-01  2.2788E-01 -1.3198E+00
             9.0766E-01
 GRADIENT:  -9.1238E+00  5.0913E+00  1.5641E+01  8.6554E-01 -3.0762E+01  4.1124E+00 -3.2365E+00 -2.7270E-01 -3.3665E+00  1.9031E+00
             1.0903E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1376.67498866921        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1256
 NPARAMETR:  1.0126E+00  7.3868E-01  1.7641E-01  9.1096E-01  3.3448E-01  1.0387E+00  9.2567E-01  1.6579E+00  1.1386E+00  1.7623E-01
             2.2411E+00
 PARAMETER:  1.1251E-01 -2.0290E-01 -1.6349E+00  6.7438E-03 -9.9518E-01  1.3800E-01  2.2767E-02  6.0552E-01  2.2983E-01 -1.6360E+00
             9.0698E-01
 GRADIENT:  -7.0488E-02 -1.9737E+00 -4.0339E+00 -2.5804E+00  5.6721E+00  2.6771E-01 -1.8781E-02  2.4481E-01 -2.6570E-01  1.4774E+00
             4.1039E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1377.42908304553        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1432
 NPARAMETR:  1.0133E+00  7.1372E-01  1.7832E-01  9.2307E-01  3.2763E-01  1.0403E+00  9.5077E-01  1.6588E+00  1.1370E+00  4.3209E-02
             2.2207E+00
 PARAMETER:  1.1325E-01 -2.3726E-01 -1.6242E+00  1.9955E-02 -1.0159E+00  1.3955E-01  4.9518E-02  6.0609E-01  2.2840E-01 -3.0417E+00
             8.9782E-01
 GRADIENT:  -4.7510E-01 -1.9930E+00 -1.5376E+00  1.5785E-03  6.9218E-01  3.5742E-01 -8.0338E-01  4.1290E-02 -5.3105E-01  8.1166E-02
             1.2622E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1377.47903387287        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1594
 NPARAMETR:  1.0135E+00  7.2047E-01  1.7996E-01  9.2244E-01  3.3024E-01  1.0390E+00  9.5669E-01  1.6643E+00  1.1362E+00  1.0000E-02
             2.2160E+00
 PARAMETER:  1.1341E-01 -2.2785E-01 -1.6150E+00  1.9263E-02 -1.0079E+00  1.3827E-01  5.5727E-02  6.0940E-01  2.2769E-01 -4.5890E+00
             8.9569E-01
 GRADIENT:  -3.0217E-02  1.9805E-01  1.9290E-01  2.1480E-01  1.4233E-01  5.2394E-02  6.7290E-02 -3.1711E-02  2.2576E-02  0.0000E+00
            -1.2061E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1594
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.1749E-04 -3.3981E-03 -6.0471E-03 -5.6864E-03  6.2440E-05
 SE:             2.9318E-02  2.1542E-02  2.0890E-02  2.5601E-02  3.0348E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8592E-01  8.7466E-01  7.7221E-01  8.2422E-01  8.3699E-01

 ETASHRINKSD(%)  1.7798E+00  2.7832E+01  3.0017E+01  1.4235E+01  9.8983E+01
 ETASHRINKVR(%)  3.5279E+00  4.7917E+01  5.1024E+01  2.6443E+01  9.9990E+01
 EBVSHRINKSD(%)  1.8573E+00  2.7211E+01  3.1121E+01  1.3932E+01  9.8931E+01
 EBVSHRINKVR(%)  3.6802E+00  4.7018E+01  5.2557E+01  2.5923E+01  9.9989E+01
 RELATIVEINF(%)  9.3537E+01  2.7914E+00  1.7498E+01  2.9589E+01  5.5373E-04
 EPSSHRINKSD(%)  4.1207E+01
 EPSSHRINKVR(%)  6.5434E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1377.4790338728678     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -642.32820730912965     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1377.479       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  7.20E-01  1.80E-01  9.22E-01  3.30E-01  1.04E+00  9.57E-01  1.66E+00  1.14E+00  1.00E-02  2.22E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       46.240
Stop Time:
Sun Oct 24 02:17:31 CDT 2021
