Sun Oct 24 03:32:19 CDT 2021
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
$DATA ../../../../data/SD4/SL3/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1538.35823561129        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0576E+02 -3.8497E+01 -1.7347E+01  1.5533E+01  7.3971E+01  3.4628E+01 -4.2812E+01 -1.8865E-03 -9.9065E+00 -4.9479E+01
            -6.7909E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1549.97399910283        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.5172E-01  1.0915E+00  9.1972E-01  9.5055E-01  1.1132E+00  9.2467E-01  1.5505E+00  9.4606E-01  1.0013E+00  1.3256E+00
             1.2269E+00
 PARAMETER:  5.0513E-02  1.8759E-01  1.6310E-02  4.9289E-02  2.0722E-01  2.1682E-02  5.3859E-01  4.4553E-02  1.0132E-01  3.8188E-01
             3.0449E-01
 GRADIENT:   2.4810E+02  3.0394E+00 -6.0499E+01  5.1966E+01  8.8032E+01 -1.5651E+01  4.3580E+01  8.6535E+00  2.0798E+01  7.2343E+00
             3.9367E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1556.65148005438        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      244
 NPARAMETR:  9.6010E-01  1.0267E+00  1.1611E+00  1.0666E+00  1.0868E+00  9.5373E-01  1.7361E+00  8.4441E-01  1.0425E+00  1.4347E+00
             1.1936E+00
 PARAMETER:  5.9282E-02  1.2634E-01  2.4934E-01  1.6449E-01  1.8328E-01  5.2629E-02  6.5165E-01 -6.9119E-02  1.4165E-01  4.6096E-01
             2.7697E-01
 GRADIENT:   3.6425E+00  1.9202E+01 -9.3119E+00  3.0131E+01  1.0123E+01 -4.5959E+01  2.0695E+01  3.7502E-01  2.2526E+01  1.3256E+01
             2.4359E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1564.40698879400        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      425
 NPARAMETR:  9.5341E-01  9.6236E-01  9.6769E-01  1.0600E+00  9.4705E-01  1.0542E+00  1.6418E+00  5.7723E-01  8.6657E-01  1.2181E+00
             1.1252E+00
 PARAMETER:  5.2291E-02  6.1634E-02  6.7151E-02  1.5827E-01  4.5601E-02  1.5283E-01  5.9580E-01 -4.4951E-01 -4.3218E-02  2.9727E-01
             2.1795E-01
 GRADIENT:  -1.1887E+01  2.6755E+00  4.1860E+00 -1.7555E+00 -4.8966E+00 -2.7045E+00  5.4571E-01  4.9202E-02 -2.0271E+00  2.0989E+00
             5.7898E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1564.82693350749        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      601
 NPARAMETR:  9.6139E-01  9.6637E-01  8.3242E-01  1.0464E+00  8.7959E-01  1.0617E+00  1.6319E+00  3.0257E-01  8.7323E-01  1.1205E+00
             1.1178E+00
 PARAMETER:  6.0624E-02  6.5792E-02 -8.3421E-02  1.4539E-01 -2.8299E-02  1.5989E-01  5.8974E-01 -1.0954E+00 -3.5555E-02  2.1379E-01
             2.1134E-01
 GRADIENT:   2.0904E+00 -4.2645E-01 -1.3340E+00  2.5163E+00  2.3482E+00 -4.1395E-01 -3.3487E-01  2.1418E-01 -1.6541E-01  6.7607E-02
            -3.3697E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1564.96090415379        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      776
 NPARAMETR:  9.6138E-01  1.0118E+00  7.7856E-01  1.0148E+00  8.6711E-01  1.0624E+00  1.5822E+00  6.5633E-02  8.8593E-01  1.0932E+00
             1.1176E+00
 PARAMETER:  6.0619E-02  1.1174E-01 -1.5031E-01  1.1471E-01 -4.2595E-02  1.6054E-01  5.5883E-01 -2.6237E+00 -2.1120E-02  1.8913E-01
             2.1118E-01
 GRADIENT:   7.3386E-01  4.7969E-01 -5.2902E-01  1.4453E+00  8.2807E-01 -4.2021E-01  1.2476E-01  1.3063E-02  1.7977E-01 -5.1856E-02
            -5.6197E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1564.96904369905        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      952
 NPARAMETR:  9.6106E-01  1.0053E+00  7.7560E-01  1.0172E+00  8.6191E-01  1.0636E+00  1.5883E+00  1.3351E-02  8.8281E-01  1.0887E+00
             1.1174E+00
 PARAMETER:  6.0286E-02  1.0533E-01 -1.5412E-01  1.1707E-01 -4.8605E-02  1.6163E-01  5.6266E-01 -4.2162E+00 -2.4641E-02  1.8496E-01
             2.1099E-01
 GRADIENT:   7.8492E-02 -6.5539E-02  1.2806E-03 -1.6772E-01 -1.8914E-01  1.3006E-02 -4.9812E-02  6.0982E-04 -1.7038E-02  1.3236E-03
             4.5759E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1564.97048959544        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1113
 NPARAMETR:  9.6164E-01  1.0065E+00  7.7655E-01  1.0166E+00  8.6318E-01  1.0656E+00  1.5895E+00  1.0000E-02  8.8338E-01  1.0902E+00
             1.1175E+00
 PARAMETER:  6.0889E-02  1.0650E-01 -1.5290E-01  1.1642E-01 -4.7132E-02  1.6358E-01  5.6339E-01 -4.5956E+00 -2.3995E-02  1.8635E-01
             2.1108E-01
 GRADIENT:   1.2734E+00 -1.7224E-02  9.1308E-02 -3.7678E-01 -9.9986E-02  7.9878E-01  2.4250E-01  0.0000E+00  5.3753E-02  3.6846E-02
             1.6810E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1113
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.4676E-05  3.7569E-03 -4.4763E-04 -1.0080E-02 -1.3057E-02
 SE:             2.9799E-02  2.3332E-02  1.6367E-04  2.1987E-02  2.3199E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9961E-01  8.7207E-01  6.2372E-03  6.4662E-01  5.7353E-01

 ETASHRINKSD(%)  1.6870E-01  2.1836E+01  9.9452E+01  2.6339E+01  2.2282E+01
 ETASHRINKVR(%)  3.3712E-01  3.8904E+01  9.9997E+01  4.5741E+01  3.9599E+01
 EBVSHRINKSD(%)  4.7800E-01  2.1315E+01  9.9526E+01  2.6907E+01  2.0055E+01
 EBVSHRINKVR(%)  9.5372E-01  3.8086E+01  9.9998E+01  4.6574E+01  3.6088E+01
 RELATIVEINF(%)  9.8768E+01  6.3506E+00  3.6542E-04  5.2730E+00  8.9962E+00
 EPSSHRINKSD(%)  4.3208E+01
 EPSSHRINKVR(%)  6.7747E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1564.9704895954421     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -829.81966303170395     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1564.970       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.62E-01  1.01E+00  7.77E-01  1.02E+00  8.63E-01  1.07E+00  1.59E+00  1.00E-02  8.83E-01  1.09E+00  1.12E+00
 


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
 #CPUT: Total CPU Time in Seconds,       32.793
Stop Time:
Sun Oct 24 03:32:27 CDT 2021
