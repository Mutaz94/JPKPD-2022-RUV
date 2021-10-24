Sat Oct 23 19:09:37 CDT 2021
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
$DATA ../../../../data/SD2/SL3/dat7.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      798
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

 TOT. NO. OF OBS RECS:      698
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1505.08651944128        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0833E+02  3.9416E+01 -1.6045E+01  1.8261E+02  1.2358E+02  2.5420E+01 -3.6282E+01 -2.7953E+01  1.2217E+01 -3.5341E+01
            -2.7646E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2335.97295003144        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0279E+00  1.1172E+00  9.7599E-01  9.3604E-01  9.9164E-01  1.1212E+00  1.0674E+00  9.7950E-01  7.8219E-01  1.0624E+00
             2.2440E+00
 PARAMETER:  1.2755E-01  2.1078E-01  7.5699E-02  3.3902E-02  9.1606E-02  2.1437E-01  1.6519E-01  7.9287E-02 -1.4566E-01  1.6058E-01
             9.0825E-01
 GRADIENT:   1.7337E+02  4.0354E+01 -1.0487E+00  3.2335E+01 -5.2299E+00  3.5511E+01 -6.0349E-01  2.5648E+00 -3.2322E-01 -2.1299E+00
             5.8846E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2336.36841820095        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0193E+00  1.0780E+00  9.8339E-01  9.6006E-01  9.7139E-01  1.1011E+00  1.1133E+00  8.4901E-01  7.3627E-01  1.0333E+00
             2.2395E+00
 PARAMETER:  1.1911E-01  1.7515E-01  8.3248E-02  5.9244E-02  7.0973E-02  1.9630E-01  2.0733E-01 -6.3685E-02 -2.0616E-01  1.3271E-01
             9.0625E-01
 GRADIENT:   1.5189E+02  3.7507E+01  4.1049E+00  3.4688E+01 -6.4530E+00  2.7518E+01 -1.2588E+00  5.2586E-01 -4.3836E+00 -6.4753E+00
             5.2644E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2338.26769798828        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      275
 NPARAMETR:  9.9414E-01  1.0971E+00  9.3100E-01  9.3356E-01  9.7748E-01  1.0893E+00  1.0673E+00  5.1609E-01  8.1227E-01  1.1212E+00
             2.1818E+00
 PARAMETER:  9.4119E-02  1.9269E-01  2.8507E-02  3.1247E-02  7.7224E-02  1.8552E-01  1.6513E-01 -5.6147E-01 -1.0793E-01  2.1437E-01
             8.8015E-01
 GRADIENT:  -6.6603E+00 -2.5449E+00  8.2716E-01 -5.4164E-01 -1.9443E+00 -3.2108E+00 -1.6134E+00  3.2550E-01 -1.1474E-01  7.1447E-01
             5.5193E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2338.42895671117        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  9.9729E-01  1.1488E+00  9.0121E-01  9.0213E-01  9.9886E-01  1.0983E+00  1.0391E+00  3.0024E-01  8.3282E-01  1.1438E+00
             2.1750E+00
 PARAMETER:  9.7284E-02  2.3873E-01 -4.0218E-03 -3.0019E-03  9.8862E-02  1.9377E-01  1.3833E-01 -1.1032E+00 -8.2937E-02  2.3436E-01
             8.7704E-01
 GRADIENT:  -7.2652E-01  4.8970E-01  3.1153E-02 -3.7597E-01  4.0595E-01 -3.3526E-02 -1.3425E-01  7.3456E-02 -1.6210E-01 -9.1026E-02
            -8.2350E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2338.48771164886        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      627
 NPARAMETR:  9.9758E-01  1.1158E+00  8.8250E-01  9.2036E-01  9.7059E-01  1.0985E+00  1.0637E+00  7.0736E-02  8.2511E-01  1.1253E+00
             2.1756E+00
 PARAMETER:  9.7581E-02  2.0954E-01 -2.4998E-02  1.7013E-02  7.0150E-02  1.9397E-01  1.6172E-01 -2.5488E+00 -9.2240E-02  2.1804E-01
             8.7732E-01
 GRADIENT:  -2.7264E-01 -9.9415E-03 -2.5408E-01  8.0880E-02  2.0615E-01  4.1367E-02  7.5833E-02  6.1303E-03  3.9171E-02  1.0145E-01
             1.8901E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2338.49075106818        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      789
 NPARAMETR:  9.9773E-01  1.1134E+00  8.8303E-01  9.2177E-01  9.6922E-01  1.0984E+00  1.0649E+00  1.0000E-02  8.2448E-01  1.1242E+00
             2.1754E+00
 PARAMETER:  9.7729E-02  2.0740E-01 -2.4399E-02  1.8543E-02  6.8731E-02  1.9381E-01  1.6289E-01 -4.6707E+00 -9.3006E-02  2.1704E-01
             8.7723E-01
 GRADIENT:   3.6578E-03  2.1915E-02  1.7669E-02  1.2686E-02 -3.9640E-02 -1.2776E-02 -1.8688E-03  0.0000E+00  1.3056E-02  8.1963E-05
            -2.9201E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      789
 NO. OF SIG. DIGITS IN FINAL EST.:  2.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.7686E-04 -6.9433E-03 -1.9557E-04 -1.3278E-03 -1.1946E-02
 SE:             2.9618E-02  2.1906E-02  1.1612E-04  2.2188E-02  2.4293E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7369E-01  7.5127E-01  9.2125E-02  9.5228E-01  6.2290E-01

 ETASHRINKSD(%)  7.7478E-01  2.6613E+01  9.9611E+01  2.5667E+01  1.8615E+01
 ETASHRINKVR(%)  1.5436E+00  4.6143E+01  9.9998E+01  4.4747E+01  3.3765E+01
 EBVSHRINKSD(%)  1.0003E+00  2.6923E+01  9.9635E+01  2.6420E+01  1.7789E+01
 EBVSHRINKVR(%)  1.9905E+00  4.6598E+01  9.9999E+01  4.5860E+01  3.2414E+01
 RELATIVEINF(%)  9.7957E+01  5.8327E+00  3.3254E-04  6.5868E+00  8.6846E+00
 EPSSHRINKSD(%)  2.1169E+01
 EPSSHRINKVR(%)  3.7856E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          698
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1282.8381923537231     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2338.4907510681846     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1055.6525587144615     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2338.491       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.98E-01  1.11E+00  8.83E-01  9.22E-01  9.69E-01  1.10E+00  1.06E+00  1.00E-02  8.24E-01  1.12E+00  2.18E+00
 


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
 #CPUT: Total CPU Time in Seconds,       70.665
Stop Time:
Sat Oct 23 19:09:49 CDT 2021
