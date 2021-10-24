Sat Oct 23 15:23:15 CDT 2021
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
$DATA ../../../../data/SD1/SL3/dat18.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      979
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

 TOT. NO. OF OBS RECS:      879
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
 RAW OUTPUT FILE (FILE): m18.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1233.68841498429        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4498E+02  1.0460E+02  9.5925E+01  6.0024E+01  2.0158E+02  4.3636E+01 -8.3355E+01 -2.3866E+02 -1.5699E+02 -1.9787E+01
            -4.6076E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2762.56143187818        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0896E+00  1.1155E+00  1.1246E+00  1.0164E+00  1.0566E+00  9.3000E-01  1.0427E+00  1.0215E+00  1.4370E+00  8.3484E-01
             2.1409E+00
 PARAMETER:  1.8585E-01  2.0928E-01  2.1741E-01  1.1625E-01  1.5504E-01  2.7430E-02  1.4183E-01  1.2131E-01  4.6258E-01 -8.0519E-02
             8.6120E-01
 GRADIENT:   3.8509E+02  2.4741E+01 -6.4677E+00  6.5991E+01  4.2086E+01 -2.6995E+01 -4.9979E+00 -1.2701E+01  3.4805E+01 -1.6618E+01
            -2.8744E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2795.50329325517        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0051E+00  1.0331E+00  2.3201E+00  1.0658E+00  1.2622E+00  9.0989E-01  9.5096E-01  3.9645E+00  1.2496E+00  9.8537E-01
             2.2094E+00
 PARAMETER:  1.0510E-01  1.3253E-01  9.4160E-01  1.6374E-01  3.3288E-01  5.5725E-03  4.9714E-02  1.4774E+00  3.2286E-01  8.5259E-02
             8.9273E-01
 GRADIENT:   9.6605E+01 -2.4972E+01 -5.2093E+00  2.7929E+01  3.7922E+01 -1.8396E+01  1.2550E+01  4.2115E+01  3.3972E+01  1.2231E+00
            -1.1554E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2803.21261583744        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  9.9024E-01  1.1579E+00  2.8132E+00  1.0154E+00  1.4700E+00  8.8942E-01  7.1957E-01  3.1403E+00  1.2987E+00  1.5147E+00
             2.2805E+00
 PARAMETER:  9.0188E-02  2.4660E-01  1.1343E+00  1.1528E-01  4.8523E-01 -1.7181E-02 -2.2910E-01  1.2443E+00  3.6134E-01  5.1525E-01
             9.2441E-01
 GRADIENT:  -4.1455E+01 -3.1629E+01 -2.5687E+01  2.2130E+01  8.1071E+00 -3.6786E+01  7.4521E+00 -9.0210E-01  1.2806E+01  3.8793E+01
            -7.5198E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2816.24551048116        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      478
 NPARAMETR:  1.0078E+00  1.4066E+00  4.0774E+00  8.3922E-01  1.6559E+00  9.7222E-01  4.4190E-01  3.4685E+00  1.5011E+00  1.4072E+00
             2.3770E+00
 PARAMETER:  1.0775E-01  4.4116E-01  1.5054E+00 -7.5279E-02  6.0436E-01  7.1830E-02 -7.1667E-01  1.3437E+00  5.0623E-01  4.4163E-01
             9.6586E-01
 GRADIENT:   1.9222E+00  6.5451E+00  3.5349E+00  5.9990E+00 -3.2807E+00 -1.4858E-01  2.5243E-02 -2.6580E+00 -1.1728E+00  1.9480E+00
             6.7559E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2816.59738042442        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      653
 NPARAMETR:  1.0067E+00  1.4954E+00  3.5345E+00  7.6866E-01  1.6739E+00  9.7296E-01  4.2393E-01  3.5205E+00  1.6161E+00  1.4063E+00
             2.3781E+00
 PARAMETER:  1.0666E-01  5.0237E-01  1.3626E+00 -1.6311E-01  6.1517E-01  7.2585E-02 -7.5819E-01  1.3586E+00  5.7999E-01  4.4098E-01
             9.6630E-01
 GRADIENT:  -5.3463E-01  3.8028E-02 -5.1044E-01  9.3436E-01  4.9496E-01  1.3755E-01 -7.5295E-01  5.8930E-01 -8.8516E-02 -3.8029E-01
            -6.2648E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2816.64753209692        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  1.0070E+00  1.5411E+00  3.4259E+00  7.3929E-01  1.6894E+00  9.7292E-01  4.7453E-01  3.4922E+00  1.6384E+00  1.4171E+00
             2.3801E+00
 PARAMETER:  1.0694E-01  5.3251E-01  1.3314E+00 -2.0206E-01  6.2438E-01  7.2551E-02 -6.4543E-01  1.3505E+00  5.9369E-01  4.4858E-01
             9.6715E-01
 GRADIENT:  -7.0546E-02  3.0268E-01 -5.1677E-03  2.2131E-01  1.6748E-01  4.5696E-02 -2.7181E-02 -7.6778E-02 -1.1878E-01  7.9817E-02
            -9.3187E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      815
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5547E-03 -4.7532E-02 -3.3770E-02  1.8201E-02 -3.1410E-02
 SE:             2.9478E-02  1.2890E-02  1.5964E-02  2.6390E-02  2.4026E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3094E-01  2.2647E-04  3.4404E-02  4.9039E-01  1.9110E-01

 ETASHRINKSD(%)  1.2447E+00  5.6818E+01  4.6517E+01  1.1590E+01  1.9510E+01
 ETASHRINKVR(%)  2.4738E+00  8.1353E+01  7.1396E+01  2.1836E+01  3.5214E+01
 EBVSHRINKSD(%)  1.5011E+00  5.9488E+01  5.1092E+01  1.1091E+01  1.5976E+01
 EBVSHRINKVR(%)  2.9796E+00  8.3588E+01  7.6081E+01  2.0951E+01  2.9399E+01
 RELATIVEINF(%)  9.6966E+01  2.6366E+00  1.1845E+01  1.3959E+01  3.9851E+01
 EPSSHRINKSD(%)  1.7383E+01
 EPSSHRINKVR(%)  3.1745E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          879
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1615.4939413738146     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2816.6475320969198     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1201.1535907231053     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2816.648       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.54E+00  3.43E+00  7.39E-01  1.69E+00  9.73E-01  4.75E-01  3.49E+00  1.64E+00  1.42E+00  2.38E+00
 


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
 #CPUT: Total CPU Time in Seconds,       52.725
Stop Time:
Sat Oct 23 15:23:25 CDT 2021
