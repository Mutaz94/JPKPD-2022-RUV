Sat Oct 23 22:23:16 CDT 2021
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
$DATA ../../../../data/SD3/A2/dat90.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m90.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1388.20979876345        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6573E+02  1.4433E+01  4.5067E+01  2.3559E+01  3.9462E+01  2.5456E+01  9.5804E+00 -3.4624E+01  2.0138E+01 -3.8051E+01
            -1.3281E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1748.92612904739        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1224E+00  1.0571E+00  9.6889E-01  1.0443E+00  9.8091E-01  1.1493E+00  9.1721E-01  1.0245E+00  8.7635E-01  1.0082E+00
             2.2624E+00
 PARAMETER:  2.1543E-01  1.5554E-01  6.8391E-02  1.4331E-01  8.0726E-02  2.3912E-01  1.3586E-02  1.2425E-01 -3.1992E-02  1.0820E-01
             9.1640E-01
 GRADIENT:   3.8784E+02  3.5435E+01  2.9479E+00  5.1905E+01 -1.4333E+01  3.3319E+01  4.6039E+00 -3.8146E-02  8.5083E-01  7.0976E+00
             7.3754E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1768.16192074734        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0076E+00  1.0747E+00  1.2993E+00  1.0240E+00  1.1858E+00  9.8698E-01  5.1136E-01  1.2486E+00  1.0629E+00  1.2731E+00
             2.1206E+00
 PARAMETER:  1.0758E-01  1.7205E-01  3.6184E-01  1.2374E-01  2.7038E-01  8.6898E-02 -5.7068E-01  3.2201E-01  1.6104E-01  3.4146E-01
             8.5170E-01
 GRADIENT:   1.5835E+02  1.7768E+01  1.5724E+00  3.5653E+01  9.5385E+00  2.4112E+00 -1.8426E-01 -5.6033E+00  1.7188E+01  6.1256E+00
             2.7482E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1770.63355305153        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.7419E-01  1.0655E+00  1.1577E+00  1.0134E+00  1.0959E+00  9.6110E-01  6.3398E-01  1.3898E+00  9.8006E-01  1.0927E+00
             2.0680E+00
 PARAMETER:  7.3854E-02  1.6349E-01  2.4643E-01  1.1333E-01  1.9158E-01  6.0323E-02 -3.5574E-01  4.2913E-01  7.9856E-02  1.8862E-01
             8.2656E-01
 GRADIENT:   7.9692E+01  7.9549E+00 -4.4126E+00  1.6737E+01  4.4611E+00 -6.1587E+00 -9.0691E-01  1.8336E+00  4.0040E+00 -3.2337E-01
             7.3673E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1771.46190475638        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      405
 NPARAMETR:  9.7989E-01  1.0883E+00  1.3896E+00  1.0119E+00  1.1955E+00  9.9630E-01  7.1859E-01  1.6504E+00  9.4976E-01  1.1743E+00
             2.0731E+00
 PARAMETER:  7.9682E-02  1.8458E-01  4.2898E-01  1.1181E-01  2.7853E-01  9.6289E-02 -2.3046E-01  6.0102E-01  4.8454E-02  2.6067E-01
             8.2904E-01
 GRADIENT:  -4.3805E-01 -2.7332E-01 -7.9344E-02  1.5209E+00  1.7215E-01 -1.9654E-01  1.4289E-01  1.3986E-01  4.9679E-01  2.8278E-01
            -2.2580E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1771.74733294119        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      583
 NPARAMETR:  9.8292E-01  1.3558E+00  1.3073E+00  8.4791E-01  1.2617E+00  1.0015E+00  6.2941E-01  1.9463E+00  1.0756E+00  1.1675E+00
             2.0725E+00
 PARAMETER:  8.2769E-02  4.0443E-01  3.6798E-01 -6.4981E-02  3.3246E-01  1.0147E-01 -3.6296E-01  7.6592E-01  1.7289E-01  2.5484E-01
             8.2877E-01
 GRADIENT:   1.3651E+00  2.1415E+01  5.1369E+00  1.5283E+01 -9.9486E+00  7.9768E-01 -2.1616E+00 -2.2243E+00 -4.3296E+00 -1.3492E+00
            -2.2620E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1774.88876075677        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      761
 NPARAMETR:  9.8477E-01  1.8286E+00  8.2210E-01  5.6033E-01  1.3868E+00  1.0098E+00  5.9558E-01  2.3530E+00  1.4288E+00  1.2492E+00
             2.0445E+00
 PARAMETER:  8.4649E-02  7.0357E-01 -9.5890E-02 -4.7924E-01  4.2702E-01  1.0978E-01 -4.1822E-01  9.5568E-01  4.5687E-01  3.2247E-01
             8.1514E-01
 GRADIENT:  -1.4462E-01  6.8196E+01 -1.3468E-01  3.8396E+01  2.1161E+00  2.2423E+00 -3.4205E+00 -1.1444E+01 -6.3603E+00  8.2881E-02
             7.5696E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1779.03681031605        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  9.8083E-01  2.0087E+00  6.9979E-01  4.1305E-01  1.4898E+00  1.0007E+00  5.9690E-01  2.7081E+00  1.8698E+00  1.3801E+00
             2.0087E+00
 PARAMETER:  8.0645E-02  7.9748E-01 -2.5697E-01 -7.8418E-01  4.9867E-01  1.0073E-01 -4.1601E-01  1.0963E+00  7.2582E-01  4.2219E-01
             7.9749E-01
 GRADIENT:  -6.4185E+00  1.2268E+01 -8.4469E-01  1.2978E+01  3.0867E+00 -8.2476E-01  4.8626E+00 -3.0821E+00  4.3854E+00  3.7084E+00
            -1.8443E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1780.16136969498        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1113
 NPARAMETR:  9.8354E-01  2.2308E+00  4.4367E-01  2.5037E-01  1.5175E+00  1.0001E+00  5.6300E-01  2.4787E+00  2.4713E+00  1.3377E+00
             2.0386E+00
 PARAMETER:  8.3401E-02  9.0236E-01 -7.1267E-01 -1.2848E+00  5.1705E-01  1.0009E-01 -4.7448E-01  1.0077E+00  1.0047E+00  3.9092E-01
             8.1229E-01
 GRADIENT:  -2.3755E+00  1.7376E+01  3.6590E+00 -6.4563E-01 -6.3387E+00 -6.5131E-01 -8.5536E-01 -1.1340E+00 -2.2862E+00  4.4121E-01
             2.6163E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1780.83824684791        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1290
 NPARAMETR:  9.8416E-01  2.3965E+00  2.1515E-01  1.3798E-01  1.5883E+00  1.0013E+00  5.7584E-01  2.1367E+00  3.3740E+00  1.3612E+00
             2.0381E+00
 PARAMETER:  8.4037E-02  9.7399E-01 -1.4364E+00 -1.8807E+00  5.6265E-01  1.0128E-01 -4.5192E-01  8.5925E-01  1.3161E+00  4.0838E-01
             8.1203E-01
 GRADIENT:  -9.9422E-01  7.4467E-01 -1.0455E-01  9.4858E-01 -3.0913E-01 -1.5940E-02  1.0849E+00 -1.9515E+00  1.8965E+00 -1.8382E-02
             4.3055E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1781.21431974874        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1468
 NPARAMETR:  9.8573E-01  2.4037E+00  2.3816E-01  1.3285E-01  1.6571E+00  1.0032E+00  5.8307E-01  2.7728E+00  3.2334E+00  1.4351E+00
             1.9940E+00
 PARAMETER:  8.5631E-02  9.7699E-01 -1.3348E+00 -1.9186E+00  6.0510E-01  1.0322E-01 -4.3946E-01  1.1198E+00  1.2735E+00  4.6126E-01
             7.9013E-01
 GRADIENT:   3.7326E+00 -2.3268E+00  1.7278E+00 -1.9326E+00  1.2694E+01  3.6435E-01  7.9613E-01 -4.0306E+00 -6.0686E+00 -1.8174E+00
            -4.4519E+00

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1781.40389025040        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1566
 NPARAMETR:  9.8517E-01  2.4078E+00  2.4543E-01  1.3111E-01  1.6618E+00  1.0019E+00  5.8178E-01  2.9482E+00  3.2578E+00  1.4524E+00
             1.9980E+00
 PARAMETER:  8.5062E-02  9.7865E-01 -1.3080E+00 -1.9322E+00  6.0943E-01  1.0290E-01 -4.4097E-01  1.1841E+00  1.2842E+00  4.7208E-01
             7.9021E-01
 GRADIENT:   2.8987E+00 -7.8618E+01 -8.9189E+02 -2.1750E-01  9.7028E+02  4.0193E-01  3.1007E-01  9.7036E+02  9.0501E+02 -1.2421E+03
            -1.4892E+03
 NUMSIGDIG:         6.2         3.8         2.3         3.3         2.3         1.7         2.5         2.3         2.3         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1566
 NO. OF SIG. DIGITS IN FINAL EST.:  1.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4777E-03 -4.9145E-02 -1.8746E-02  5.4654E-02 -6.8731E-02
 SE:             2.9685E-02  2.3495E-02  9.2465E-03  1.9721E-02  2.0020E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6030E-01  3.6460E-02  4.2625E-02  5.5817E-03  5.9662E-04

 ETASHRINKSD(%)  5.5042E-01  2.1290E+01  6.9023E+01  3.3933E+01  3.2932E+01
 ETASHRINKVR(%)  1.0978E+00  3.8047E+01  9.0404E+01  5.6352E+01  5.5019E+01
 EBVSHRINKSD(%)  1.2299E+00  1.6536E+01  7.4238E+01  4.8795E+01  3.0750E+01
 EBVSHRINKVR(%)  2.4447E+00  3.0338E+01  9.3363E+01  7.3780E+01  5.2044E+01
 RELATIVEINF(%)  9.7473E+01  2.5905E+01  5.0902E+00  9.7665E+00  3.7176E+01
 EPSSHRINKSD(%)  3.0011E+01
 EPSSHRINKVR(%)  5.1016E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1781.4038902503994     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -862.46535704572671     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1781.404       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  2.41E+00  2.45E-01  1.31E-01  1.66E+00  1.00E+00  5.82E-01  2.96E+00  3.27E+00  1.45E+00  1.99E+00
 


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
 #CPUT: Total CPU Time in Seconds,      128.836
Stop Time:
Sat Oct 23 22:23:35 CDT 2021
