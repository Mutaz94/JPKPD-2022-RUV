Sun Oct 24 02:46:35 CDT 2021
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
$DATA ../../../../data/SD4/S1/dat53.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m53.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1677.75975674208        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4980E+02 -5.4063E+01 -2.4839E+01 -1.5456E+01  4.0585E+01  2.7136E+01  6.5316E+00  6.4679E+00  4.0688E+01  1.2496E+01
             3.4338E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1689.23464081300        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.7493E-01  1.0839E+00  1.0236E+00  1.0278E+00  1.0246E+00  1.1268E+00  9.7369E-01  9.8360E-01  8.3415E-01  9.4589E-01
             9.0768E-01
 PARAMETER:  7.4610E-02  1.8054E-01  1.2336E-01  1.2738E-01  1.2434E-01  2.1934E-01  7.3334E-02  8.3468E-02 -8.1339E-02  4.4371E-02
             3.1397E-03
 GRADIENT:  -7.9326E+00  7.0203E+00 -8.6289E+00  3.0330E+01  2.9617E+01  1.1683E+01 -3.4194E+00 -1.0655E+00  3.5859E+00 -3.4496E+00
            -1.0171E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1690.39969380160        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.8010E-01  1.0225E+00  8.9206E-01  1.0531E+00  9.1199E-01  1.1226E+00  1.0909E+00  9.1488E-01  7.3748E-01  8.0707E-01
             9.3642E-01
 PARAMETER:  7.9903E-02  1.2224E-01 -1.4219E-02  1.5178E-01  7.8768E-03  2.1563E-01  1.8704E-01  1.1038E-02 -2.0451E-01 -1.1434E-01
             3.4310E-02
 GRADIENT:  -7.6594E-01  8.5452E+00 -6.1709E+00  2.3948E+01  1.3633E+01  9.8464E+00 -2.7541E+00  2.4248E+00 -4.5307E+00 -3.6774E+00
             3.3373E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1691.96626116297        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.8351E-01  1.0773E+00  6.1693E-01  9.8117E-01  7.7282E-01  1.0964E+00  1.0740E+00  4.4905E-01  7.5803E-01  7.1899E-01
             9.2324E-01
 PARAMETER:  8.3368E-02  1.7443E-01 -3.8299E-01  8.0987E-02 -1.5771E-01  1.9208E-01  1.7135E-01 -7.0063E-01 -1.7703E-01 -2.2991E-01
             2.0134E-02
 GRADIENT:   1.0018E+00  1.8375E+00  1.4028E+00 -3.8923E+00 -9.1256E+00 -7.1951E-01  1.1243E+00  1.4685E+00  1.4632E+00  2.4940E+00
             1.0554E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1692.46728438070        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      700
 NPARAMETR:  9.8351E-01  1.0866E+00  5.8625E-01  9.7202E-01  7.6071E-01  1.0983E+00  1.0648E+00  1.3935E-01  7.5659E-01  7.2343E-01
             9.2133E-01
 PARAMETER:  8.3372E-02  1.8310E-01 -4.3402E-01  7.1623E-02 -1.7350E-01  1.9376E-01  1.6275E-01 -1.8708E+00 -1.7893E-01 -2.2374E-01
             1.8063E-02
 GRADIENT:   6.1327E-01  2.8188E+00  2.4099E+00 -1.2993E+00 -5.4837E+00 -1.5839E-01  2.4687E-03  9.7974E-02 -2.3184E-01  9.8380E-01
             1.4652E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1692.52351971189        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      875
 NPARAMETR:  9.8321E-01  1.0750E+00  5.7871E-01  9.7633E-01  7.5233E-01  1.0990E+00  1.0760E+00  2.7121E-02  7.5293E-01  7.1212E-01
             9.2062E-01
 PARAMETER:  8.3068E-02  1.7229E-01 -4.4695E-01  7.6045E-02 -1.8458E-01  1.9436E-01  1.7321E-01 -3.5074E+00 -1.8379E-01 -2.3951E-01
             1.7292E-02
 GRADIENT:  -1.3387E-03 -4.0784E-01 -2.1218E-01 -6.3387E-01  9.3881E-02  5.9033E-02  2.5318E-01  3.9859E-03 -9.1116E-02  1.9932E-01
            -1.2361E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1692.52840205120        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1051            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8412E-01  1.0809E+00  5.7797E-01  9.7311E-01  7.5447E-01  1.1021E+00  1.0697E+00  1.0000E-02  7.5575E-01  7.1221E-01
             9.2109E-01
 PARAMETER:  8.3993E-02  1.7777E-01 -4.4824E-01  7.2742E-02 -1.8174E-01  1.9722E-01  1.6734E-01 -4.9692E+00 -1.8004E-01 -2.3938E-01
             1.7803E-02
 GRADIENT:   4.5688E+02  7.3951E+01  1.7371E+01  7.5646E+01  2.4948E+01  1.5692E+02  6.9954E+00  0.0000E+00  7.5222E+00  1.5532E+00
             8.0371E-01

0ITERATION NO.:   32    OBJECTIVE VALUE:  -1692.52840205120        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1108
 NPARAMETR:  9.8412E-01  1.0809E+00  5.7797E-01  9.7311E-01  7.5447E-01  1.1021E+00  1.0697E+00  1.0000E-02  7.5575E-01  7.1221E-01
             9.2109E-01
 PARAMETER:  8.3993E-02  1.7777E-01 -4.4824E-01  7.2742E-02 -1.8174E-01  1.9722E-01  1.6734E-01 -4.9692E+00 -1.8004E-01 -2.3938E-01
             1.7803E-02
 GRADIENT:   1.6770E+00 -6.5974E-02  9.2511E-02 -2.5929E-01 -1.3420E-01  1.2173E+00  1.2435E-02  0.0000E+00  2.1684E-02  1.0216E-02
             3.7174E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1108
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.2218E-05 -2.4479E-03 -4.5046E-04 -8.2887E-04 -9.4403E-03
 SE:             2.9895E-02  2.3510E-02  1.9937E-04  2.4103E-02  2.1905E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9807E-01  9.1707E-01  2.3855E-02  9.7257E-01  6.6649E-01

 ETASHRINKSD(%)  1.0000E-10  2.1239E+01  9.9332E+01  1.9252E+01  2.6617E+01
 ETASHRINKVR(%)  1.0000E-10  3.7967E+01  9.9996E+01  3.4798E+01  4.6149E+01
 EBVSHRINKSD(%)  2.8896E-01  2.0874E+01  9.9387E+01  1.9860E+01  2.6468E+01
 EBVSHRINKVR(%)  5.7708E-01  3.7390E+01  9.9996E+01  3.5776E+01  4.5930E+01
 RELATIVEINF(%)  9.9365E+01  3.9229E+00  2.9192E-04  4.2582E+00  3.4258E+00
 EPSSHRINKSD(%)  4.4613E+01
 EPSSHRINKVR(%)  6.9323E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1692.5284020511961     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -957.37757548745788     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     4.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1692.528       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  1.08E+00  5.78E-01  9.73E-01  7.54E-01  1.10E+00  1.07E+00  1.00E-02  7.56E-01  7.12E-01  9.21E-01
 


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
 #CPUT: Total CPU Time in Seconds,       30.244
Stop Time:
Sun Oct 24 02:46:42 CDT 2021
