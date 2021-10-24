Sat Oct 23 17:33:54 CDT 2021
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
$DATA ../../../../data/SD2/A1/dat93.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      800
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

 TOT. NO. OF OBS RECS:      700
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
 RAW OUTPUT FILE (FILE): m93.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2712.35770815895        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4609E+02 -1.0148E+01  6.3922E+01  4.2492E+01  8.3428E+01  4.3730E+01  2.1237E+01 -9.0278E+01  4.7333E+01 -1.9276E+01
            -4.1007E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2779.94867514495        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.4463E-01  1.0536E+00  8.1444E-01  1.0389E+00  9.1127E-01  9.6728E-01  8.6657E-01  1.4014E+00  8.1637E-01  9.4004E-01
             1.2670E+00
 PARAMETER:  4.3034E-02  1.5222E-01 -1.0525E-01  1.3813E-01  7.0858E-03  6.6735E-02 -4.3210E-02  4.3749E-01 -1.0288E-01  3.8170E-02
             3.3661E-01
 GRADIENT:   1.5315E+02  8.6108E+01 -2.1290E+01  1.5812E+02  1.4006E+01  1.8656E+01  1.5911E+00 -1.3790E+01 -1.3975E-01  9.4983E-01
             3.4916E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2785.38219237557        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:      215             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4577E-01  1.0630E+00  9.0933E-01  1.0251E+00  9.4996E-01  9.8189E-01  9.3241E-01  1.6621E+00  8.3099E-01  9.2257E-01
             1.2606E+00
 PARAMETER:  4.4244E-02  1.6111E-01  4.9520E-03  1.2477E-01  4.8666E-02  8.1728E-02  3.0021E-02  6.0806E-01 -8.5133E-02  1.9410E-02
             3.3160E-01
 GRADIENT:   1.6360E+02  7.5516E+01 -8.4776E+00  1.1215E+02 -6.7297E+00  2.5424E+01  6.4941E+00  2.0682E+00  9.7200E+00  2.3146E+00
             3.4985E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2787.13721363309        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.4577E-01  1.0630E+00  1.0431E+00  1.0251E+00  1.0228E+00  9.8612E-01  9.4007E-01  1.6621E+00  7.8671E-01  9.6163E-01
             1.2606E+00
 PARAMETER:  4.4249E-02  1.6111E-01  1.4222E-01  1.2476E-01  1.2251E-01  8.6021E-02  3.8201E-02  6.0809E-01 -1.3990E-01  6.0878E-02
             3.3158E-01
 GRADIENT:  -8.6933E+01  4.4634E+00 -3.7312E+00  4.6410E+01  8.1900E+00  2.9558E-01  5.1938E-01 -1.9977E+01  7.8667E-01 -1.0919E+00
             3.1030E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2790.57250286516        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  9.8315E-01  1.0569E+00  1.0510E+00  1.0058E+00  1.0194E+00  9.8551E-01  9.3677E-01  1.8694E+00  7.8359E-01  9.6757E-01
             1.2333E+00
 PARAMETER:  8.3007E-02  1.5531E-01  1.4979E-01  1.0581E-01  1.1922E-01  8.5401E-02  3.4686E-02  7.2562E-01 -1.4387E-01  6.7028E-02
             3.0967E-01
 GRADIENT:   2.0963E+00 -2.3229E+01 -6.3430E+00 -9.2814E+00 -1.3750E+01  3.6668E+00 -1.0941E+00  2.0915E+00  9.0117E-01  2.0236E+00
            -2.4594E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2792.62286346838        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  9.7826E-01  1.1307E+00  1.1959E+00  9.7454E-01  1.1207E+00  9.6838E-01  9.4569E-01  2.0507E+00  7.7042E-01  1.0161E+00
             1.2338E+00
 PARAMETER:  7.8025E-02  2.2288E-01  2.7890E-01  7.4210E-02  2.1393E-01  6.7866E-02  4.4159E-02  8.1820E-01 -1.6082E-01  1.1598E-01
             3.1013E-01
 GRADIENT:  -9.2786E+00 -7.2106E+00 -8.8248E-02  3.9119E+00  2.4950E-01 -3.2603E+00  9.1533E-01 -6.1307E-01 -1.5446E+00 -1.4930E+00
            -4.1076E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2792.66734305734        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      914            RESET HESSIAN, TYPE II
 NPARAMETR:  9.7858E-01  1.1314E+00  1.1951E+00  9.7484E-01  1.1199E+00  9.7646E-01  9.2762E-01  2.0508E+00  7.8593E-01  1.0170E+00
             1.2347E+00
 PARAMETER:  7.8351E-02  2.2344E-01  2.7820E-01  7.4522E-02  2.1326E-01  7.6177E-02  2.4872E-02  8.1825E-01 -1.4089E-01  1.1690E-01
             3.1080E-01
 GRADIENT:   2.5815E+02  7.5020E+01  6.9615E+00  5.3067E+01  2.9915E+01  2.5787E+01  3.0461E+00  1.3681E+01  4.1470E+00  4.7625E-01
             5.3948E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2792.68924282745        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1098
 NPARAMETR:  9.8242E-01  1.1320E+00  1.1942E+00  9.7533E-01  1.1187E+00  9.7646E-01  9.2752E-01  2.0510E+00  7.8561E-01  1.0208E+00
             1.2348E+00
 PARAMETER:  8.2262E-02  2.2400E-01  2.7751E-01  7.5022E-02  2.1220E-01  7.6177E-02  2.4760E-02  8.1832E-01 -1.4129E-01  1.2060E-01
             3.1087E-01
 GRADIENT:   6.6227E-01 -6.1850E+00  3.5662E-01  7.7822E+00 -1.4301E+00  7.3819E-02  4.7417E-02 -5.7431E-02  3.0638E-02 -5.3792E-01
             1.1650E+00

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2792.68924282745        NO. OF FUNC. EVALS.:  25
 CUMULATIVE NO. OF FUNC. EVALS.:     1123
 NPARAMETR:  9.8242E-01  1.1320E+00  1.1942E+00  9.7533E-01  1.1187E+00  9.7646E-01  9.2752E-01  2.0510E+00  7.8561E-01  1.0208E+00
             1.2348E+00
 PARAMETER:  8.2262E-02  2.2400E-01  2.7751E-01  7.5022E-02  2.1220E-01  7.6177E-02  2.4760E-02  8.1832E-01 -1.4129E-01  1.2060E-01
             3.1087E-01
 GRADIENT:   2.9907E+04 -6.6858E+03  5.3956E+03 -2.9913E+04  1.4093E+04  1.1439E-02  4.4136E-02 -1.3127E+02  2.1168E+04  2.4798E+04
            -9.6135E+03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1123
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4748E-04 -1.4304E-02 -3.8349E-02  6.1164E-03 -4.0872E-02
 SE:             2.9860E-02  2.1038E-02  2.0497E-02  2.3619E-02  2.1654E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7736E-01  4.9655E-01  6.1352E-02  7.9567E-01  5.9097E-02

 ETASHRINKSD(%)  1.0000E-10  2.9520E+01  3.1332E+01  2.0872E+01  2.7455E+01
 ETASHRINKVR(%)  1.0000E-10  5.0326E+01  5.2847E+01  3.7388E+01  4.7373E+01
 EBVSHRINKSD(%)  4.2371E-01  2.9350E+01  3.3177E+01  2.2530E+01  2.5200E+01
 EBVSHRINKVR(%)  8.4563E-01  5.0086E+01  5.5347E+01  3.9984E+01  4.4050E+01
 RELATIVEINF(%)  9.9143E+01  9.6728E+00  2.4265E+01  1.2797E+01  2.0383E+01
 EPSSHRINKSD(%)  2.5392E+01
 EPSSHRINKVR(%)  4.4337E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2792.6892428274496     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1506.1752963409078     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2792.689       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  1.13E+00  1.19E+00  9.75E-01  1.12E+00  9.76E-01  9.28E-01  2.05E+00  7.86E-01  1.02E+00  1.23E+00
 


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
 #CPUT: Total CPU Time in Seconds,       71.986
Stop Time:
Sat Oct 23 17:34:07 CDT 2021
