Sat Oct 23 22:46:07 CDT 2021
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
$DATA ../../../../data/SD3/A3/dat79.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   183.219142539519        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3214E+02  3.8843E+01  2.8990E+02 -1.6143E+02  1.6554E+02  6.0374E+01 -7.3917E+01 -3.6320E+02 -2.0249E+02 -1.2758E+02
            -3.7893E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1459.84638442678        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0620E+00  1.0673E+00  9.1197E-01  1.1827E+00  1.0008E+00  8.0907E-01  1.0063E+00  1.0224E+00  9.3678E-01  9.8629E-01
             5.3171E+00
 PARAMETER:  1.6015E-01  1.6510E-01  7.8514E-03  2.6778E-01  1.0078E-01 -1.1187E-01  1.0626E-01  1.2215E-01  3.4692E-02  8.6199E-02
             1.7709E+00
 GRADIENT:   4.4508E+01  4.0073E+00 -1.7232E+01  3.1885E+01 -1.1491E+01 -3.1811E+00  9.8848E+00  8.0417E+00  2.5052E+01  1.9201E+01
             3.2926E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1518.40074554331        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0074E+00  7.6615E-01  6.0315E-01  1.2894E+00  6.7025E-01  8.3579E-01  1.8886E+00  4.5145E-02  8.4722E-01  3.8043E-01
             4.2535E+00
 PARAMETER:  1.0735E-01 -1.6637E-01 -4.0559E-01  3.5414E-01 -3.0011E-01 -7.9375E-02  7.3582E-01 -2.9979E+00 -6.5798E-02 -8.6645E-01
             1.5477E+00
 GRADIENT:  -4.2308E+01  1.9022E+01 -4.3677E+01  6.2430E+01  3.9720E+01 -2.5478E+00  3.1227E+01  3.2014E-02  2.6870E+01  5.5898E+00
             1.9788E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1541.47561838394        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9620E-01  7.3355E-01  2.9293E-01  1.1368E+00  4.0657E-01  8.6535E-01  1.1557E+00  2.0333E-02  9.0792E-01  2.4659E-01
             3.4432E+00
 PARAMETER:  9.6191E-02 -2.0986E-01 -1.1278E+00  2.2825E-01 -8.0000E-01 -4.4621E-02  2.4467E-01 -3.7955E+00  3.4032E-03 -1.3000E+00
             1.3364E+00
 GRADIENT:  -4.3767E+01  3.0268E+01  1.0532E+00  2.0180E+01 -3.9604E+00 -3.7013E+00 -7.8628E+00 -1.4112E-03 -1.0339E+01  9.1862E-01
             3.4650E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1545.73960644054        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      388
 NPARAMETR:  1.0219E+00  6.7445E-01  3.4763E-01  1.1876E+00  4.3624E-01  8.7270E-01  1.3761E+00  3.1092E-02  8.9863E-01  1.7308E-01
             3.3508E+00
 PARAMETER:  1.2169E-01 -2.9386E-01 -9.5661E-01  2.7190E-01 -7.2956E-01 -3.6161E-02  4.1927E-01 -3.3708E+00 -6.8847E-03 -1.6540E+00
             1.3092E+00
 GRADIENT:   8.8736E+00  2.7766E+00 -7.0590E+00  7.9015E+00  9.8226E+00  8.0836E-01 -2.4550E+00 -7.7147E-03  5.9969E-01  1.8049E-01
            -1.4850E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1547.48487033938        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      564
 NPARAMETR:  1.0181E+00  4.9243E-01  3.5731E-01  1.2749E+00  3.9299E-01  8.6733E-01  1.7971E+00  3.9237E-01  8.4794E-01  1.8415E-02
             3.3629E+00
 PARAMETER:  1.1795E-01 -6.0839E-01 -9.2916E-01  3.4283E-01 -8.3397E-01 -4.2335E-02  6.8618E-01 -8.3555E-01 -6.4943E-02 -3.8946E+00
             1.3128E+00
 GRADIENT:   2.8326E+00  2.4570E+00 -6.0316E+00  3.1543E+01  9.2583E+00  2.3468E-01 -3.7080E-03 -1.6628E+00 -4.4639E+00 -3.5044E-03
            -5.6239E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1562.89049356910        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      743
 NPARAMETR:  1.0105E+00  3.6037E-01  2.0782E-01  1.2159E+00  2.5364E-01  8.7849E-01  1.6478E+00  1.1859E+00  1.2684E+00  1.0000E-02
             2.9079E+00
 PARAMETER:  1.1049E-01 -9.2063E-01 -1.4711E+00  2.9551E-01 -1.2718E+00 -2.9548E-02  5.9946E-01  2.7049E-01  3.3775E-01 -5.9283E+00
             1.1674E+00
 GRADIENT:  -6.1037E+00  1.0975E+01 -1.5692E+00  5.0213E+01  4.9949E+00 -1.8206E+00  1.3979E+01 -6.8498E-01 -6.2036E+00  0.0000E+00
             1.1526E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1569.18472753669        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      918
 NPARAMETR:  1.0094E+00  3.5397E-01  1.8508E-01  1.1315E+00  2.3849E-01  8.8847E-01  8.1576E-01  1.3149E+00  1.4540E+00  1.0000E-02
             2.8184E+00
 PARAMETER:  1.0940E-01 -9.3855E-01 -1.5869E+00  2.2356E-01 -1.3334E+00 -1.8252E-02 -1.0364E-01  3.7374E-01  4.7430E-01 -6.0589E+00
             1.1362E+00
 GRADIENT:  -2.0058E-01 -2.6717E-01  2.0662E+00 -4.8902E+00 -1.8293E-01  8.2283E-01  7.7945E-01 -1.3484E+00  1.7384E+00  0.0000E+00
             2.2759E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1569.29590978015        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1093
 NPARAMETR:  1.0093E+00  3.5109E-01  1.8473E-01  1.1366E+00  2.3847E-01  8.8655E-01  6.3234E-01  1.3466E+00  1.4592E+00  1.0000E-02
             2.8129E+00
 PARAMETER:  1.0922E-01 -9.4671E-01 -1.5889E+00  2.2800E-01 -1.3335E+00 -2.0417E-02 -3.5833E-01  3.9761E-01  4.7792E-01 -5.8991E+00
             1.1342E+00
 GRADIENT:  -7.9663E-02 -1.9133E-01  8.3342E-02 -8.4722E-01  8.3835E-02 -1.1389E-02  2.1169E-02  1.2173E-01  9.2397E-02  0.0000E+00
             2.9272E-01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1569.29622411404        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  1.0093E+00  3.5135E-01  1.8483E-01  1.1374E+00  2.3860E-01  8.8662E-01  6.2685E-01  1.3465E+00  1.4590E+00  1.0000E-02
             2.8123E+00
 PARAMETER:  1.0924E-01 -9.4597E-01 -1.5883E+00  2.2878E-01 -1.3330E+00 -2.0342E-02 -3.6705E-01  3.9752E-01  4.7774E-01 -5.8780E+00
             1.1340E+00
 GRADIENT:  -2.5409E-02 -2.0285E-02  3.3331E-02 -2.4568E-01 -1.1531E-02 -2.9751E-03 -4.1765E-03  2.8139E-02  3.1406E-03  0.0000E+00
             2.4479E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1150
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2641E-03 -1.5104E-02  8.4317E-03 -5.3691E-03  6.1677E-04
 SE:             2.8873E-02  7.8837E-03  2.4305E-02  2.7484E-02  4.1109E-04
 N:                     100         100         100         100         100

 P VAL.:         9.3750E-01  5.5375E-02  7.2866E-01  8.4511E-01  1.3353E-01

 ETASHRINKSD(%)  3.2706E+00  7.3589E+01  1.8574E+01  7.9263E+00  9.8623E+01
 ETASHRINKVR(%)  6.4342E+00  9.3024E+01  3.3697E+01  1.5224E+01  9.9981E+01
 EBVSHRINKSD(%)  3.2682E+00  7.3409E+01  1.8651E+01  6.6553E+00  9.8797E+01
 EBVSHRINKVR(%)  6.4297E+00  9.2929E+01  3.3823E+01  1.2868E+01  9.9986E+01
 RELATIVEINF(%)  9.3394E+01  1.4106E+00  1.0801E+01  6.2808E+01  1.0278E-03
 EPSSHRINKSD(%)  2.7798E+01
 EPSSHRINKVR(%)  4.7869E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1569.2962241140449     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -650.35769090937220     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1569.296       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  3.51E-01  1.85E-01  1.14E+00  2.39E-01  8.87E-01  6.27E-01  1.35E+00  1.46E+00  1.00E-02  2.81E+00
 


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
 #CPUT: Total CPU Time in Seconds,       93.743
Stop Time:
Sat Oct 23 22:46:22 CDT 2021
