Sun Oct 24 01:55:30 CDT 2021
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
$DATA ../../../../data/SD4/A1/dat16.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m16.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1400.71775192189        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8091E+02 -2.9509E+01 -5.0980E+01  5.3354E+01  1.8859E+02  7.1089E+01 -2.2228E+01  6.3969E+00 -4.2372E+00 -6.4775E+01
            -4.9991E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1526.28313120118        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.1563E+00  1.0338E+00  1.0403E+00  1.0283E+00  8.9601E-01  1.0188E+00  1.0451E+00  9.4495E-01  1.0019E+00  1.0265E+00
             1.9091E+00
 PARAMETER:  2.4522E-01  1.3327E-01  1.3951E-01  1.2787E-01 -9.7992E-03  1.1865E-01  1.4409E-01  4.3382E-02  1.0188E-01  1.2619E-01
             7.4662E-01
 GRADIENT:   5.2707E+02  2.0170E+01  1.3927E+00  2.9044E+01 -8.2952E+00  2.5924E+01 -8.5957E-01  6.4561E+00  4.1513E+00  4.5712E+00
             2.2799E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1536.42443283822        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.1187E+00  8.4884E-01  7.9198E-01  1.1242E+00  7.4790E-01  9.8932E-01  1.3356E+00  2.6003E-01  9.3055E-01  1.0155E+00
             1.7285E+00
 PARAMETER:  2.1217E-01 -6.3880E-02 -1.3322E-01  2.1708E-01 -1.9048E-01  8.9265E-02  3.8939E-01 -1.2470E+00  2.8017E-02  1.1535E-01
             6.4724E-01
 GRADIENT:   4.7141E+02  7.3730E+00 -2.7231E+01  8.7676E+01  3.9073E+01  2.7066E+01  5.5482E+00  1.1839E+00  9.7513E+00  1.5522E+01
            -5.5130E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1539.72586154312        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      294
 NPARAMETR:  1.0523E+00  8.4398E-01  5.8616E-01  1.0749E+00  6.2022E-01  9.2545E-01  1.3894E+00  1.8609E-01  9.9759E-01  8.1378E-01
             1.5935E+00
 PARAMETER:  1.5099E-01 -6.9630E-02 -4.3417E-01  1.7227E-01 -3.7768E-01  2.2522E-02  4.2889E-01 -1.5816E+00  9.7589E-02 -1.0606E-01
             5.6593E-01
 GRADIENT:   2.4161E+01 -1.7900E+00 -1.9759E+01 -2.4255E+00  1.5439E+01  2.6585E+00  6.4595E+00  7.4643E-01  2.4591E+01  9.6681E+00
            -4.7407E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1547.83456074483        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      469
 NPARAMETR:  1.0434E+00  4.9414E-01  6.7021E-01  1.2895E+00  5.4838E-01  9.1083E-01  2.0005E+00  6.2854E-02  8.0133E-01  7.7233E-01
             1.7832E+00
 PARAMETER:  1.4247E-01 -6.0493E-01 -3.0017E-01  3.5424E-01 -5.0078E-01  6.5966E-03  7.9342E-01 -2.6669E+00 -1.2149E-01 -1.5835E-01
             6.7841E-01
 GRADIENT:   4.5971E+00  9.5586E+00  1.0992E+01  8.6587E+00 -2.3098E+01 -1.4965E+00  1.9480E-01  6.3316E-02 -2.1938E+00 -2.6600E-01
             8.3829E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1548.68524875982        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      644
 NPARAMETR:  1.0409E+00  3.1753E-01  6.9158E-01  1.3735E+00  5.2658E-01  9.1333E-01  2.5978E+00  1.2123E-02  7.8231E-01  8.1668E-01
             1.7441E+00
 PARAMETER:  1.4005E-01 -1.0472E+00 -2.6878E-01  4.1736E-01 -5.4136E-01  9.3449E-03  1.0547E+00 -4.3126E+00 -1.4550E-01 -1.0251E-01
             6.5621E-01
 GRADIENT:   8.6956E+00 -8.3587E-01  1.3514E+00 -6.3167E+00 -5.1316E+00  5.4970E-01 -1.6084E+00  2.5778E-03 -1.4104E+00  3.9470E-01
             1.8324E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1548.71744198547        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      820
 NPARAMETR:  1.0382E+00  2.7414E-01  7.1006E-01  1.3983E+00  5.2871E-01  9.1138E-01  2.8956E+00  1.0000E-02  7.8222E-01  8.3217E-01
             1.7378E+00
 PARAMETER:  1.3750E-01 -1.1941E+00 -2.4240E-01  4.3525E-01 -5.3731E-01  7.2100E-03  1.1632E+00 -4.8797E+00 -1.4562E-01 -8.3717E-02
             6.5261E-01
 GRADIENT:   5.8721E+00 -9.8265E-01  6.7096E-01 -7.8524E+00 -2.8546E+00  1.1951E-01  2.3054E-03  0.0000E+00  8.3916E-01  5.9263E-01
             9.1842E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1548.82471057851        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      998
 NPARAMETR:  1.0360E+00  3.0692E-01  7.6295E-01  1.3975E+00  5.6213E-01  9.1072E-01  2.7204E+00  1.1677E-02  7.8266E-01  8.5710E-01
             1.7421E+00
 PARAMETER:  1.3536E-01 -1.0812E+00 -1.7056E-01  4.3468E-01 -4.7602E-01  6.4824E-03  1.1008E+00 -4.3501E+00 -1.4505E-01 -5.4204E-02
             6.5511E-01
 GRADIENT:   3.0604E-01  8.0665E-02  6.7986E-01 -6.6207E-01 -1.0504E+00 -1.4381E-03 -1.5235E-02  2.0827E-03  4.6587E-02  1.1800E-01
             2.0020E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1548.82681521743        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1094
 NPARAMETR:  1.0362E+00  3.0703E-01  7.6241E-01  1.3963E+00  5.6241E-01  9.1078E-01  2.7256E+00  1.0000E-02  7.8258E-01  8.5612E-01
             1.7414E+00
 PARAMETER:  1.3559E-01 -1.0808E+00 -1.7127E-01  4.3383E-01 -4.7552E-01  6.5490E-03  1.1027E+00 -5.3932E+00 -1.4516E-01 -5.5342E-02
             6.5467E-01
 GRADIENT:   9.6041E-01 -1.5306E-01 -1.7467E-01 -2.5291E+00  5.6911E-01  2.6464E-02  2.1900E-01  0.0000E+00  6.3805E-02  1.6767E-02
             2.6628E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1094
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.0637E-04  3.1800E-02 -2.3646E-04 -2.6208E-02 -5.9868E-04
 SE:             2.9498E-02  1.7256E-02  1.6854E-04  2.4919E-02  2.2024E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8630E-01  6.5356E-02  1.6064E-01  2.9293E-01  9.7831E-01

 ETASHRINKSD(%)  1.1767E+00  4.2189E+01  9.9435E+01  1.6518E+01  2.6216E+01
 ETASHRINKVR(%)  2.3396E+00  6.6579E+01  9.9997E+01  3.0307E+01  4.5560E+01
 EBVSHRINKSD(%)  1.4354E+00  4.9545E+01  9.9361E+01  1.3983E+01  2.2131E+01
 EBVSHRINKVR(%)  2.8501E+00  7.4543E+01  9.9996E+01  2.6011E+01  3.9365E+01
 RELATIVEINF(%)  9.5783E+01  4.9855E+00  2.5319E-04  1.7586E+01  3.6808E+00
 EPSSHRINKSD(%)  3.8783E+01
 EPSSHRINKVR(%)  6.2525E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1548.8268152174260     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -813.67598865368780     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     4.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1548.827       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  3.07E-01  7.62E-01  1.40E+00  5.62E-01  9.11E-01  2.73E+00  1.00E-02  7.83E-01  8.56E-01  1.74E+00
 


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
 #CPUT: Total CPU Time in Seconds,       30.669
Stop Time:
Sun Oct 24 01:55:37 CDT 2021
