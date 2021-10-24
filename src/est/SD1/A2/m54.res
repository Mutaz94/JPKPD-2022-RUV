Sat Oct 23 13:59:41 CDT 2021
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
$DATA ../../../../data/SD1/A2/dat54.csv ignore=@
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1873.15458665954        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8867E+02  2.8379E+02  2.3090E+02  1.3649E+02  2.0709E+02  2.2868E+01 -1.7271E+02 -2.8433E+02 -3.1802E+01 -1.0290E+02
            -3.3616E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2937.41707949957        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0944E+00  8.7519E-01  8.1882E-01  1.0483E+00  8.5388E-01  1.0846E+00  9.5778E-01  1.0504E+00  9.0840E-01  7.2087E-01
             2.7519E+00
 PARAMETER:  1.9022E-01 -3.3311E-02 -9.9896E-02  1.4714E-01 -5.7970E-02  1.8125E-01  5.6865E-02  1.4921E-01  3.9259E-03 -2.2730E-01
             1.1123E+00
 GRADIENT:   1.7479E+02  1.2059E+00 -1.1943E+01  1.0745E+00  9.7421E+00  2.9510E+01 -9.9351E+00  1.2687E+01 -6.7660E+00  8.8484E+00
             4.0088E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2946.51537871149        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      203
 NPARAMETR:  1.0938E+00  6.2930E-01  5.9498E-01  1.1897E+00  6.1161E-01  9.8883E-01  1.2890E+00  9.4665E-01  9.8915E-01  5.0059E-01
             2.6649E+00
 PARAMETER:  1.8962E-01 -3.6315E-01 -4.1923E-01  2.7370E-01 -3.9166E-01  8.8766E-02  3.5389E-01  4.5171E-02  8.9092E-02 -5.9197E-01
             1.0802E+00
 GRADIENT:   6.3165E+01 -2.3708E-02 -2.7438E+01  7.6741E+01  2.5781E+01 -1.6570E+01  1.5300E+01  1.7943E+01  2.2880E+01  5.6989E+00
             3.5319E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2989.35801088960        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  1.0483E+00  6.5703E-01  6.0894E-01  1.1650E+00  6.2699E-01  1.0219E+00  1.1366E+00  7.3411E-01  9.1990E-01  6.5200E-01
             2.1463E+00
 PARAMETER:  1.4720E-01 -3.2003E-01 -3.9604E-01  2.5273E-01 -3.6682E-01  1.2167E-01  2.2801E-01 -2.0910E-01  1.6512E-02 -3.2770E-01
             8.6376E-01
 GRADIENT:  -1.2076E+01  4.1195E+01 -7.3281E+00  7.6752E+01  3.4235E+01 -3.4471E+00 -4.8524E+00  3.9708E+00 -1.9576E+01 -1.2435E+01
            -2.9734E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3008.36296205130        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      556
 NPARAMETR:  1.0573E+00  4.8695E-01  4.6933E-01  1.2102E+00  4.7340E-01  1.0286E+00  1.3210E+00  8.0134E-01  1.0621E+00  7.2359E-01
             2.0494E+00
 PARAMETER:  1.5574E-01 -6.1959E-01 -6.5644E-01  2.9077E-01 -6.4780E-01  1.2822E-01  3.7839E-01 -1.2147E-01  1.6029E-01 -2.2352E-01
             8.1752E-01
 GRADIENT:   9.1486E+00  9.0402E+00  3.9578E+01  9.4776E+01  5.5892E+00 -2.9530E+00  2.5464E+01  9.4145E+00  1.0184E+01  7.7315E+00
            -2.8646E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3045.55847710174        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  1.0564E+00  2.7909E-01  2.0596E-01  1.0380E+00  2.5784E-01  1.0492E+00  1.1462E+00  1.0429E+00  1.1707E+00  5.1405E-01
             1.9276E+00
 PARAMETER:  1.5485E-01 -1.1762E+00 -1.4801E+00  1.3731E-01 -1.2554E+00  1.4805E-01  2.3641E-01  1.4204E-01  2.5757E-01 -5.6543E-01
             7.5627E-01
 GRADIENT:   3.7747E+00  6.1347E+00  7.3052E+00  1.4140E+01 -1.4702E+01 -8.2585E-02  7.4969E-01 -1.0189E+01 -2.5538E+00 -3.4334E+00
            -3.3056E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3046.34331791340        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  1.0546E+00  2.6771E-01  1.9349E-01  1.0122E+00  2.4840E-01  1.0492E+00  1.1279E+00  1.1169E+00  1.1987E+00  5.4172E-01
             1.9219E+00
 PARAMETER:  1.5312E-01 -1.2178E+00 -1.5425E+00  1.1212E-01 -1.2927E+00  1.4800E-01  2.2035E-01  2.1052E-01  2.8126E-01 -5.1302E-01
             7.5333E-01
 GRADIENT:   2.1069E-01 -2.8732E-01  1.5134E-01 -3.2433E-01 -2.2352E-02  1.3277E-02  1.4478E-01  1.4659E-01  4.8634E-02 -8.4077E-02
             2.1511E-01

0ITERATION NO.:   32    OBJECTIVE VALUE:  -3046.34400892529        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      964
 NPARAMETR:  1.0545E+00  2.6784E-01  1.9351E-01  1.0125E+00  2.4843E-01  1.0491E+00  1.1254E+00  1.1143E+00  1.1986E+00  5.4406E-01
             1.9219E+00
 PARAMETER:  1.5302E-01 -1.2174E+00 -1.5424E+00  1.1241E-01 -1.2926E+00  1.4791E-01  2.1816E-01  2.0822E-01  2.8113E-01 -5.0870E-01
             7.5330E-01
 GRADIENT:   1.5305E-02  9.8965E-02 -2.7523E-02  4.3809E-02 -4.1686E-04 -2.1835E-02 -3.0929E-02  1.5937E-03  1.1476E-01  1.9155E-02
            -9.8064E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      964
 NO. OF SIG. DIGITS IN FINAL EST.:  2.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.5934E-04  8.9738E-03  2.2974E-03 -1.4157E-03  8.7193E-03
 SE:             2.9676E-02  2.4541E-02  2.2722E-02  2.8437E-02  1.9492E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7421E-01  7.1461E-01  9.1946E-01  9.6029E-01  6.5463E-01

 ETASHRINKSD(%)  5.8084E-01  1.7786E+01  2.3879E+01  4.7331E+00  3.4701E+01
 ETASHRINKVR(%)  1.1583E+00  3.2408E+01  4.2055E+01  9.2423E+00  5.7360E+01
 EBVSHRINKSD(%)  8.5458E-01  1.7088E+01  2.3527E+01  4.3607E+00  3.5520E+01
 EBVSHRINKVR(%)  1.7019E+00  3.1256E+01  4.1519E+01  8.5312E+00  5.8423E+01
 RELATIVEINF(%)  9.8293E+01  1.6616E+01  9.2768E+00  5.3560E+01  4.6877E+00
 EPSSHRINKSD(%)  2.1696E+01
 EPSSHRINKVR(%)  3.8685E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3046.3440089252927     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1392.2546491568819     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.85
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3046.344       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  2.68E-01  1.94E-01  1.01E+00  2.48E-01  1.05E+00  1.13E+00  1.11E+00  1.20E+00  5.44E-01  1.92E+00
 


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
 #CPUT: Total CPU Time in Seconds,       65.031
Stop Time:
Sat Oct 23 13:59:52 CDT 2021
