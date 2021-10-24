Sat Oct 23 13:59:21 CDT 2021
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
$DATA ../../../../data/SD1/A2/dat52.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m52.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2608.79360752000        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0984E+02  8.2073E+01  9.7078E+01 -6.1682E+00  1.3678E+02  4.8434E+01 -1.0206E+02 -3.9209E+01 -1.5382E+01 -3.1728E+01
            -2.3403E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3209.94475663760        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.6944E-01  7.0465E-01  7.0599E-01  1.1723E+00  6.5720E-01  9.3307E-01  1.3777E+00  7.5289E-01  1.0384E+00  8.3193E-01
             1.7328E+00
 PARAMETER:  6.8967E-02 -2.5005E-01 -2.4815E-01  2.5897E-01 -3.1976E-01  3.0730E-02  4.2040E-01 -1.8383E-01  1.3773E-01 -8.4004E-02
             6.4975E-01
 GRADIENT:   4.3930E+01  4.3749E+01  4.8721E+01  6.5130E+01  3.8249E+01  1.9849E+00  2.4525E+01  2.0866E+01  8.6431E+00  1.4339E+01
            -4.9259E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3219.01183980603        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      187
 NPARAMETR:  1.0538E+00  4.5872E-01  4.0406E-01  1.2083E+00  3.8425E-01  1.0354E+00  1.1510E+00  5.8579E-01  1.0461E+00  7.0043E-01
             1.6789E+00
 PARAMETER:  1.5243E-01 -6.7932E-01 -8.0620E-01  2.8918E-01 -8.5647E-01  1.3477E-01  2.4065E-01 -4.3479E-01  1.4511E-01 -2.5607E-01
             6.1812E-01
 GRADIENT:   1.0166E+02  8.1927E+01  7.9620E+01 -2.8998E+01 -1.5115E+02  2.5306E+01 -3.1667E+01  2.1976E+01 -2.7403E+01  1.1284E+01
            -6.1209E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3229.17582064240        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      365
 NPARAMETR:  1.0787E+00  4.3796E-01  4.0708E-01  1.2153E+00  3.8674E-01  9.5294E-01  1.2215E+00  2.6451E-01  1.1306E+00  7.5508E-01
             1.7202E+00
 PARAMETER:  1.7578E-01 -7.2563E-01 -7.9875E-01  2.9502E-01 -8.5001E-01  5.1794E-02  3.0006E-01 -1.2299E+00  2.2274E-01 -1.8093E-01
             6.4245E-01
 GRADIENT:   1.7544E+02  1.7955E+01  1.0236E+02 -2.1200E+01 -1.1713E+02 -1.0837E+01 -1.3823E+01  3.9836E+00  6.2548E+00  2.0742E+01
            -1.1056E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3244.31924122873        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      541
 NPARAMETR:  1.0026E+00  3.7847E-01  3.3531E-01  1.2150E+00  3.4364E-01  9.5593E-01  1.3249E+00  1.4093E-01  1.1427E+00  6.3899E-01
             1.7140E+00
 PARAMETER:  1.0257E-01 -8.7163E-01 -9.9270E-01  2.9476E-01 -9.6815E-01  5.4925E-02  3.8132E-01 -1.8595E+00  2.3344E-01 -3.4787E-01
             6.3884E-01
 GRADIENT:  -1.6713E+00  1.1713E+00  5.4794E+00  3.4369E+00 -6.7447E+00  9.5030E-01 -3.0811E-01  7.5028E-01  1.6664E+00  9.9046E-01
             1.2955E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3244.60376684911        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      719
 NPARAMETR:  1.0027E+00  3.7942E-01  3.3625E-01  1.2137E+00  3.4502E-01  9.5652E-01  1.3270E+00  6.3587E-02  1.1368E+00  6.4406E-01
             1.7142E+00
 PARAMETER:  1.0274E-01 -8.6912E-01 -9.8989E-01  2.9368E-01 -9.6415E-01  5.5547E-02  3.8292E-01 -2.6554E+00  2.2824E-01 -3.3997E-01
             6.3896E-01
 GRADIENT:  -1.2679E+00  6.5663E-01  6.4582E-01  1.1190E+00 -4.1846E-02  1.2241E+00  1.7284E-01  1.4784E-01  1.1080E-01  6.8083E-01
            -2.4188E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3244.70038618764        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      895
 NPARAMETR:  1.0033E+00  3.7622E-01  3.3291E-01  1.2116E+00  3.4220E-01  9.5343E-01  1.3259E+00  1.2466E-02  1.1377E+00  6.4242E-01
             1.7142E+00
 PARAMETER:  1.0326E-01 -8.7759E-01 -9.9989E-01  2.9195E-01 -9.7237E-01  5.2313E-02  3.8205E-01 -4.2848E+00  2.2898E-01 -3.4251E-01
             6.3896E-01
 GRADIENT:  -4.0030E-02 -4.0601E-02  2.3734E-02 -9.4267E-03 -1.1488E-02 -3.1879E-02 -7.9887E-03  5.5988E-03 -9.3211E-03  3.1717E-02
            -2.5871E-03

0ITERATION NO.:   33    OBJECTIVE VALUE:  -3244.70105818398        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  1.0033E+00  3.7628E-01  3.3296E-01  1.2116E+00  3.4225E-01  9.5350E-01  1.3259E+00  1.0000E-02  1.1377E+00  6.4243E-01
             1.7142E+00
 PARAMETER:  1.0326E-01 -8.7743E-01 -9.9972E-01  2.9197E-01 -9.7223E-01  5.2381E-02  3.8208E-01 -4.5775E+00  2.2897E-01 -3.4250E-01
             6.3897E-01
 GRADIENT:  -4.3366E-02 -1.3024E-02  1.3795E-02 -7.9929E-03  2.5475E-03 -5.6732E-03 -3.1806E-03  0.0000E+00 -3.7943E-03  2.2117E-02
            -5.6619E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      988
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1659E-03  1.2602E-03 -1.1175E-04 -2.7662E-03 -2.0179E-03
 SE:             2.9639E-02  2.5157E-02  3.0818E-04  2.9133E-02  2.5220E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6862E-01  9.6005E-01  7.1691E-01  9.2435E-01  9.3623E-01

 ETASHRINKSD(%)  7.0531E-01  1.5722E+01  9.8968E+01  2.3999E+00  1.5511E+01
 ETASHRINKVR(%)  1.4056E+00  2.8972E+01  9.9989E+01  4.7421E+00  2.8616E+01
 EBVSHRINKSD(%)  8.7580E-01  1.4230E+01  9.8908E+01  2.1982E+00  1.6223E+01
 EBVSHRINKVR(%)  1.7439E+00  2.6434E+01  9.9988E+01  4.3481E+00  2.9814E+01
 RELATIVEINF(%)  9.8237E+01  1.8561E+01  9.4759E-04  8.6581E+01  5.1503E+00
 EPSSHRINKSD(%)  2.0273E+01
 EPSSHRINKVR(%)  3.6436E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3244.7010581839800     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1590.6116984155692     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3244.701       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  3.76E-01  3.33E-01  1.21E+00  3.42E-01  9.53E-01  1.33E+00  1.00E-02  1.14E+00  6.42E-01  1.71E+00
 


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
 #CPUT: Total CPU Time in Seconds,       63.936
Stop Time:
Sat Oct 23 13:59:33 CDT 2021
