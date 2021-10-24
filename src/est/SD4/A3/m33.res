Sun Oct 24 02:28:05 CDT 2021
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
$DATA ../../../../data/SD4/A3/dat33.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m33.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   723.483410760837        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4054E+02  8.9948E+01  1.0487E+02  1.5776E+01  1.7330E+02  4.4915E+01 -5.8088E+01 -5.0664E+01 -1.3537E+02 -1.2781E+02
            -4.2847E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1167.82689191144        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0065E+00  1.0090E+00  9.6647E-01  1.1238E+00  9.6933E-01  8.3056E-01  9.9497E-01  9.9673E-01  1.0865E+00  9.5082E-01
             4.6759E+00
 PARAMETER:  1.0652E-01  1.0898E-01  6.5900E-02  2.1669E-01  6.8845E-02 -8.5654E-02  9.4959E-02  9.6725E-02  1.8300E-01  4.9568E-02
             1.6424E+00
 GRADIENT:   9.8099E+01  8.2509E+00 -1.0699E+01  2.5828E+01 -4.7248E+00 -3.8888E+00  7.0733E+00  5.4152E+00  8.3372E+00  1.7370E+01
             7.6868E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1176.42723182049        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.8145E-01  1.0063E+00  6.5006E-01  1.0729E+00  7.7346E-01  8.8982E-01  5.9338E-01  4.6317E-01  1.3256E+00  4.4734E-01
             4.6429E+00
 PARAMETER:  8.1281E-02  1.0629E-01 -3.3069E-01  1.7035E-01 -1.5688E-01 -1.6741E-02 -4.2192E-01 -6.6965E-01  3.8183E-01 -7.0444E-01
             1.6353E+00
 GRADIENT:   2.9887E+01  2.0990E+00 -1.0636E+01  2.1115E+01  1.8377E+00  1.6490E+01  2.5361E+00  2.1189E+00  1.6280E+01  5.3415E+00
             8.8050E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1184.65961576361        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.5411E-01  6.8366E-01  8.3586E-01  1.2127E+00  7.4571E-01  8.4189E-01  6.1085E-01  1.8015E-01  1.0823E+00  3.1136E-01
             4.1301E+00
 PARAMETER:  5.3023E-02 -2.8030E-01 -7.9291E-02  2.9282E-01 -1.9341E-01 -7.2110E-02 -3.9291E-01 -1.6140E+00  1.7912E-01 -1.0668E+00
             1.5183E+00
 GRADIENT:   1.0082E+01  1.2665E+00  7.0469E+00 -6.6809E+00 -1.0298E+01  2.1947E+00  2.1355E-01  3.3031E-01 -1.2836E-01  1.6920E+00
            -7.0638E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1185.76526839421        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.5073E-01  5.6836E-01  7.6990E-01  1.2789E+00  6.8050E-01  8.3741E-01  6.2655E-01  7.7723E-02  1.0431E+00  1.2282E-01
             4.1669E+00
 PARAMETER:  4.9475E-02 -4.6500E-01 -1.6150E-01  3.4601E-01 -2.8492E-01 -7.7442E-02 -3.6753E-01 -2.4546E+00  1.4222E-01 -1.9970E+00
             1.5272E+00
 GRADIENT:  -7.2538E-01  1.7145E+00  1.9021E-01  3.8122E+00 -6.0667E-01  9.8325E-02 -1.1844E-01  7.7593E-02 -4.4572E-01  2.4821E-01
            -2.5318E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1186.39290021906        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      457
 NPARAMETR:  9.5979E-01  4.7956E-01  8.3380E-01  1.3535E+00  6.9340E-01  8.3934E-01  6.5451E-01  3.4595E-02  9.8815E-01  5.5858E-02
             4.2641E+00
 PARAMETER:  5.8961E-02 -6.3490E-01 -8.1767E-02  4.0267E-01 -2.6615E-01 -7.5139E-02 -3.2387E-01 -3.2640E+00  8.8080E-02 -2.7849E+00
             1.5502E+00
 GRADIENT:   2.8877E+00 -2.2845E-01  1.7244E-01 -2.6248E+00 -4.7927E-01 -1.5398E-02 -3.5714E-02  1.4277E-02 -7.8691E-02  4.8763E-02
             3.7610E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1186.47545468289        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      635
 NPARAMETR:  9.5588E-01  3.2886E-01  7.7810E-01  1.4292E+00  6.2161E-01  8.3706E-01  1.0704E+00  1.0000E-02  9.4341E-01  1.0000E-02
             4.2274E+00
 PARAMETER:  5.4882E-02 -1.0121E+00 -1.5090E-01  4.5713E-01 -3.7545E-01 -7.7856E-02  1.6799E-01 -5.3234E+00  4.1743E-02 -5.2230E+00
             1.5416E+00
 GRADIENT:  -3.1557E-01  6.1886E-04  3.8604E-01 -3.7560E-01 -2.5443E-01 -5.9645E-01 -2.7029E-02  0.0000E+00  3.1309E-01  0.0000E+00
            -6.1477E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1186.51056931457        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  9.5473E-01  2.1693E-01  6.6892E-01  1.4665E+00  5.3264E-01  8.3732E-01  2.2268E+00  1.0000E-02  9.2095E-01  1.0000E-02
             4.2065E+00
 PARAMETER:  5.3674E-02 -1.4282E+00 -3.0209E-01  4.8289E-01 -5.2992E-01 -7.7550E-02  9.0056E-01 -8.0222E+00  1.7649E-02 -8.6060E+00
             1.5366E+00
 GRADIENT:   7.7653E-01  5.2615E-01  7.0913E-01  2.1551E+00 -1.3372E+00 -6.0117E-01  1.4533E-01  0.0000E+00  1.5254E-01  0.0000E+00
            -4.7507E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1186.54728136298        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      991
 NPARAMETR:  9.5173E-01  1.4260E-01  6.1855E-01  1.4882E+00  4.8954E-01  8.4153E-01  2.8252E+00  1.0000E-02  9.1506E-01  1.0000E-02
             4.1968E+00
 PARAMETER:  5.0528E-02 -1.8477E+00 -3.8038E-01  4.9758E-01 -6.1428E-01 -7.2538E-02  1.1386E+00 -1.1099E+01  1.1236E-02 -1.2122E+01
             1.5343E+00
 GRADIENT:  -2.3054E+00  3.9118E-01  1.3658E+00  3.4282E+00 -2.9566E+00  6.2309E-01 -9.2137E-02  0.0000E+00 -7.8308E-01  0.0000E+00
            -1.4994E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1186.58200454560        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1168
 NPARAMETR:  9.4962E-01  6.9486E-02  5.9811E-01  1.5134E+00  4.6682E-01  8.3878E-01  4.5939E+00  1.0000E-02  9.1057E-01  1.0000E-02
             4.1893E+00
 PARAMETER:  4.8306E-02 -2.5666E+00 -4.1397E-01  5.1434E-01 -6.6182E-01 -7.5807E-02  1.6247E+00 -1.6722E+01  6.3141E-03 -1.8580E+01
             1.5325E+00
 GRADIENT:   6.1954E-02  5.2460E-02 -4.2700E-02  5.1117E-01  1.2314E-01 -5.7087E-02  3.7281E-02  0.0000E+00  5.3582E-02  0.0000E+00
             3.9937E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1186.58856733361        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1352             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4849E-01  4.1365E-02  5.7473E-01  1.5161E+00  4.4910E-01  8.3947E-01  6.0843E+00  1.0000E-02  9.1434E-01  1.0000E-02
             4.1831E+00
 PARAMETER:  4.7120E-02 -3.0853E+00 -4.5385E-01  5.1616E-01 -7.0051E-01 -7.4989E-02  1.9057E+00 -2.0922E+01  1.0448E-02 -2.3363E+01
             1.5311E+00
 GRADIENT:   1.9544E+01  3.6735E-01  2.0041E+00  3.6642E+01  7.7423E+00  1.3477E+00  1.5518E-01  0.0000E+00  7.6458E-01  0.0000E+00
             1.2265E+01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1186.58856733361        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1417
 NPARAMETR:  9.4841E-01  4.1365E-02  5.7485E-01  1.5163E+00  4.4899E-01  8.3925E-01  5.9694E+00  1.0000E-02  9.1376E-01  1.0000E-02
             4.1823E+00
 PARAMETER:  4.7120E-02 -3.0853E+00 -4.5385E-01  5.1616E-01 -7.0051E-01 -7.4989E-02  1.9057E+00 -2.0922E+01  1.0448E-02 -2.3363E+01
             1.5311E+00
 GRADIENT:   1.1700E-01  2.3478E-07 -8.3612E-02 -9.1662E-02  1.8661E-01  3.9016E-02  5.7897E-03  0.0000E+00  5.9424E-02  0.0000E+00
             8.4861E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1417
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.2244E-04 -1.2735E-03  1.2602E-04 -1.7137E-02 -2.1692E-05
 SE:             2.7610E-02  3.2072E-03  1.6692E-04  2.4877E-02  2.3771E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8490E-01  6.9132E-01  4.5027E-01  4.9091E-01  9.2729E-01

 ETASHRINKSD(%)  7.5026E+00  8.9255E+01  9.9441E+01  1.6660E+01  9.9204E+01
 ETASHRINKVR(%)  1.4442E+01  9.8846E+01  9.9997E+01  3.0544E+01  9.9994E+01
 EBVSHRINKSD(%)  7.0570E+00  8.9254E+01  9.9381E+01  1.6439E+01  9.9165E+01
 EBVSHRINKVR(%)  1.3616E+01  9.8845E+01  9.9996E+01  3.0175E+01  9.9993E+01
 RELATIVEINF(%)  6.3888E+01  4.1485E-02  1.5304E-04  4.8692E+00  2.0010E-04
 EPSSHRINKSD(%)  1.9657E+01
 EPSSHRINKVR(%)  3.5449E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1186.5885673336138     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -451.43774076987563     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1186.589       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.48E-01  4.14E-02  5.75E-01  1.52E+00  4.49E-01  8.39E-01  6.08E+00  1.00E-02  9.14E-01  1.00E-02  4.18E+00
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       41.668
Stop Time:
Sun Oct 24 02:28:15 CDT 2021
