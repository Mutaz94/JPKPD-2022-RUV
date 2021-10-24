Sat Oct 23 22:21:45 CDT 2021
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
$DATA ../../../../data/SD3/A2/dat85.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1547.73201754641        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5893E+02 -1.3050E+01  1.5555E+01 -1.6716E+01  1.6562E+01  6.6410E+01 -4.8206E+01  1.0617E+00 -3.0330E+01 -9.4922E+00
            -1.0063E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1790.10817407794        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.6368E-01  1.1056E+00  1.1244E+00  1.0348E+00  1.0806E+00  8.6069E-01  1.1618E+00  9.1295E-01  1.0567E+00  8.7509E-01
             2.4959E+00
 PARAMETER:  6.3000E-02  2.0036E-01  2.1724E-01  1.3423E-01  1.7751E-01 -5.0016E-02  2.4996E-01  8.9289E-03  1.5513E-01 -3.3425E-02
             1.0146E+00
 GRADIENT:  -6.4840E-01  1.7552E+01 -2.4929E+00  1.6066E+01 -3.6435E+01 -1.2684E+01  5.5300E+00  5.7011E+00  1.3707E+01  2.1914E+01
             1.8738E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1806.19332775719        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.5322E-01  8.0411E-01  8.3496E-01  1.2231E+00  8.0533E-01  9.2788E-01  1.8320E+00  3.8594E-01  8.0971E-01  4.9614E-01
             2.3006E+00
 PARAMETER:  5.2095E-02 -1.1802E-01 -8.0377E-02  3.0136E-01 -1.1651E-01  2.5144E-02  7.0539E-01 -8.5207E-01 -1.1108E-01 -6.0089E-01
             9.3318E-01
 GRADIENT:  -4.0907E+00  3.5680E+01 -1.1724E+01  1.0970E+02 -1.0324E+01  1.3611E+01  1.9290E+01  2.3113E+00  1.7516E+00  8.3930E+00
             1.3061E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1816.73381467881        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      252
 NPARAMETR:  9.4750E-01  9.9166E-01  7.6398E-01  1.0456E+00  8.5668E-01  9.0118E-01  1.3074E+00  3.6731E-01  9.3188E-01  3.8749E-01
             1.9773E+00
 PARAMETER:  4.6070E-02  9.1623E-02 -1.6921E-01  1.4457E-01 -5.4686E-02 -4.0499E-03  3.6803E-01 -9.0155E-01  2.9452E-02 -8.4805E-01
             7.8174E-01
 GRADIENT:  -9.0250E+01 -1.1415E+01  2.1696E+00 -1.5071E+01 -8.1803E+00 -8.8645E+00 -1.1037E+01  8.5411E-01  4.5895E-01  1.1422E+00
            -9.6567E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1819.32683410350        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  9.7993E-01  9.2495E-01  7.6508E-01  1.0957E+00  8.3055E-01  9.1339E-01  1.4995E+00  1.9765E-01  8.8212E-01  3.1227E-01
             2.0136E+00
 PARAMETER:  7.9727E-02  2.1980E-02 -1.6778E-01  1.9139E-01 -8.5665E-02  9.4048E-03  5.0513E-01 -1.5212E+00 -2.5429E-02 -1.0639E+00
             7.9995E-01
 GRADIENT:  -4.3029E-01  9.1163E-02 -2.4734E-01  1.1209E+00 -3.3483E-02 -2.4543E-01 -8.6896E-03  8.2629E-02 -2.0786E-01  2.0607E-02
            -2.1040E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1819.38309446465        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  9.8066E-01  1.0132E+00  7.3555E-01  1.0421E+00  8.5316E-01  9.1422E-01  1.3974E+00  6.1385E-02  9.1259E-01  2.8161E-01
             2.0185E+00
 PARAMETER:  8.0473E-02  1.1310E-01 -2.0714E-01  1.4120E-01 -5.8803E-02  1.0314E-02  4.3463E-01 -2.6906E+00  8.5352E-03 -1.1672E+00
             8.0236E-01
 GRADIENT:  -7.8765E-01  1.8247E-01 -1.7745E-01  2.8054E-01  3.0905E-01 -2.4006E-01  3.4940E-02  7.1688E-03 -2.0165E-01  1.5795E-02
             3.1992E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1819.38672296292        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  9.8123E-01  1.0165E+00  7.3243E-01  1.0397E+00  8.5275E-01  9.1502E-01  1.3933E+00  1.8650E-02  9.1504E-01  2.7824E-01
             2.0177E+00
 PARAMETER:  8.1047E-02  1.1641E-01 -2.1139E-01  1.3891E-01 -5.9290E-02  1.1188E-02  4.3171E-01 -3.8819E+00  1.1218E-02 -1.1793E+00
             8.0197E-01
 GRADIENT:   6.0457E-01  1.0786E-01 -5.9555E-02  1.8229E-01  1.0300E-01  5.8057E-02 -8.7673E-03  6.6336E-04 -1.1658E-02 -5.7040E-03
            -6.3904E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1819.38709921136        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      957
 NPARAMETR:  9.8104E-01  1.0143E+00  7.3353E-01  1.0410E+00  8.5244E-01  9.1489E-01  1.3957E+00  1.0000E-02  9.1441E-01  2.8046E-01
             2.0177E+00
 PARAMETER:  8.0854E-02  1.1420E-01 -2.0989E-01  1.4015E-01 -5.9655E-02  1.1049E-02  4.3340E-01 -5.7324E+00  1.0526E-02 -1.1713E+00
             8.0197E-01
 GRADIENT:   1.7042E-01  7.9320E-03 -6.4384E-02  4.7617E-02  1.6046E-01  1.8518E-02  8.1263E-03  0.0000E+00  9.1488E-03 -3.6106E-04
             1.8390E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1819.38711698331        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1133
 NPARAMETR:  9.8087E-01  1.0121E+00  7.3443E-01  1.0422E+00  8.5199E-01  9.1476E-01  1.3978E+00  1.0000E-02  9.1365E-01  2.8196E-01
             2.0176E+00
 PARAMETER:  8.0681E-02  1.1207E-01 -2.0865E-01  1.4136E-01 -6.0176E-02  1.0910E-02  4.3493E-01 -7.2248E+00  9.6961E-03 -1.1660E+00
             8.0190E-01
 GRADIENT:  -2.1767E-01 -3.5350E-02 -1.4828E-02 -4.3494E-02  1.2475E-01 -2.3645E-02 -1.0499E-02  0.0000E+00  4.6883E-04 -7.1881E-04
             1.2026E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1819.38712002200        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1308
 NPARAMETR:  9.8079E-01  1.0109E+00  7.3471E-01  1.0429E+00  8.5156E-01  9.1471E-01  1.3990E+00  1.0000E-02  9.1322E-01  2.8280E-01
             2.0174E+00
 PARAMETER:  8.0608E-02  1.1086E-01 -2.0828E-01  1.4203E-01 -6.0690E-02  1.0854E-02  4.3577E-01 -8.0468E+00  9.2190E-03 -1.1630E+00
             8.0181E-01
 GRADIENT:  -3.7786E-01 -2.9083E-02  4.0111E-02 -1.0256E-01  3.2258E-02 -4.0914E-02 -2.4378E-02  0.0000E+00 -4.1749E-03 -2.7250E-04
            -3.1987E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1819.38713847697        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1486
 NPARAMETR:  9.8071E-01  1.0085E+00  7.3456E-01  1.0442E+00  8.5027E-01  9.1466E-01  1.4014E+00  1.0000E-02  9.1228E-01  2.8356E-01
             2.0171E+00
 PARAMETER:  8.0522E-02  1.0849E-01 -2.0849E-01  1.4327E-01 -6.2207E-02  1.0794E-02  4.3751E-01 -9.5767E+00  8.1904E-03 -1.1603E+00
             8.0164E-01
 GRADIENT:  -5.5980E-01 -2.8500E-02  7.4040E-02 -1.6760E-01 -9.4699E-02 -6.0408E-02 -3.9493E-02  0.0000E+00 -9.7718E-03 -6.4967E-05
            -7.0927E-02

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1819.38720346422        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1580
 NPARAMETR:  9.8080E-01  1.0090E+00  7.3370E-01  1.0439E+00  8.4997E-01  9.1473E-01  1.4012E+00  1.0000E-02  9.1235E-01  2.8226E-01
             2.0172E+00
 PARAMETER:  8.0617E-02  1.0892E-01 -2.0966E-01  1.4296E-01 -6.2559E-02  1.0872E-02  4.3730E-01 -9.2186E+00  8.2636E-03 -1.1649E+00
             8.0169E-01
 GRADIENT:  -3.4553E-01 -4.6017E-02 -3.2814E-02 -8.6727E-02  8.7647E-03 -3.7236E-02 -1.9237E-02  0.0000E+00 -3.5483E-03 -9.0182E-04
            -2.4287E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1580
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5109E-03  1.0833E-02 -2.7874E-04 -1.3530E-02 -3.9711E-04
 SE:             2.9493E-02  2.3635E-02  1.6751E-04  2.3283E-02  8.0740E-03
 N:                     100         100         100         100         100

 P VAL.:         9.5914E-01  6.4672E-01  9.6124E-02  5.6116E-01  9.6077E-01

 ETASHRINKSD(%)  1.1936E+00  2.0818E+01  9.9439E+01  2.2000E+01  7.2951E+01
 ETASHRINKVR(%)  2.3730E+00  3.7303E+01  9.9997E+01  3.9159E+01  9.2684E+01
 EBVSHRINKSD(%)  1.5052E+00  2.0457E+01  9.9441E+01  2.1919E+01  7.3361E+01
 EBVSHRINKVR(%)  2.9877E+00  3.6730E+01  9.9997E+01  3.9034E+01  9.2904E+01
 RELATIVEINF(%)  9.6583E+01  6.2717E+00  2.0007E-04  7.3860E+00  3.8853E-01
 EPSSHRINKSD(%)  2.5331E+01
 EPSSHRINKVR(%)  4.4245E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1819.3872034642202     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -900.44867025954750     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1819.387       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.01E+00  7.34E-01  1.04E+00  8.50E-01  9.15E-01  1.40E+00  1.00E-02  9.12E-01  2.82E-01  2.02E+00
 


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
 #CPUT: Total CPU Time in Seconds,      126.177
Stop Time:
Sat Oct 23 22:22:05 CDT 2021
