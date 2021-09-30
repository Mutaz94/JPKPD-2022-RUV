Wed Sep 29 22:30:42 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat50.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1838.60397788108        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3467E+02  1.7297E+01  1.2394E+01  1.7086E+01  4.7480E+01  3.9550E+01  1.1802E+00  1.5474E+00  1.1273E+01 -9.6058E+00
            -5.2671E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1923.03830816570        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.6492E-01  9.4775E-01  1.0761E+00  1.1159E+00  9.4280E-01  9.8463E-01  8.1608E-01  7.2819E-01  8.1444E-01  6.5718E-01
             1.9606E+00
 PARAMETER:  6.4288E-02  4.6332E-02  1.7330E-01  2.0967E-01  4.1099E-02  8.4516E-02 -1.0324E-01 -2.1720E-01 -1.0525E-01 -3.1980E-01
             7.7327E-01
 GRADIENT:   2.2572E+01  7.0118E+01  2.5249E+01  9.6845E+01 -5.0181E+01  5.0834E+00 -7.0286E+00  2.9472E+00 -3.2661E+01  2.1969E+00
             1.3807E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1932.41112527235        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.5551E-01  7.7527E-01  6.1310E-01  1.1881E+00  6.5783E-01  9.6940E-01  9.8775E-01  2.6840E-01  8.2322E-01  6.0836E-01
             1.8262E+00
 PARAMETER:  5.4492E-02 -1.5455E-01 -3.8922E-01  2.7239E-01 -3.1882E-01  6.8920E-02  8.7671E-02 -1.2153E+00 -9.4533E-02 -3.9699E-01
             7.0224E-01
 GRADIENT:   5.8525E+00  2.9088E+01 -5.6659E+01  1.9509E+02  7.3807E+01 -3.4084E+00 -9.4892E+00  1.4061E+00 -1.3557E+01  1.5958E+00
             1.3235E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1937.34802538542        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      249
 NPARAMETR:  9.4309E-01  7.2636E-01  4.9241E-01  1.1686E+00  5.6568E-01  9.5787E-01  1.1603E+00  2.1484E-01  7.8033E-01  6.7340E-01
             1.6449E+00
 PARAMETER:  4.1409E-02 -2.1971E-01 -6.0844E-01  2.5584E-01 -4.6972E-01  5.6959E-02  2.4866E-01 -1.4379E+00 -1.4804E-01 -2.9542E-01
             5.9765E-01
 GRADIENT:  -1.4510E+02  6.0547E+00 -9.1108E+01  1.0106E+02  7.6266E+01 -2.6431E+01 -3.4287E+00  1.3826E+00 -2.1127E+01  1.6901E+01
             7.4078E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1956.65572253448        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      425
 NPARAMETR:  9.9998E-01  5.0509E-01  6.5692E-01  1.2938E+00  5.8703E-01  1.0218E+00  1.8182E+00  7.5674E-02  7.7040E-01  6.7012E-01
             1.5014E+00
 PARAMETER:  9.9982E-02 -5.8301E-01 -3.2019E-01  3.5755E-01 -4.3268E-01  1.2154E-01  6.9783E-01 -2.4813E+00 -1.6085E-01 -3.0030E-01
             5.0641E-01
 GRADIENT:   1.1954E+01  1.6896E+01 -5.1493E-01  2.4777E+01 -1.0841E+01  1.0436E+01  5.7485E+00  1.0175E-01  1.6871E+00 -2.2549E+00
            -2.7251E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1959.21446013903        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      600
 NPARAMETR:  9.8833E-01  3.0423E-01  7.9088E-01  1.4148E+00  6.2283E-01  9.8268E-01  2.4858E+00  1.0000E-02  7.1577E-01  7.7234E-01
             1.4999E+00
 PARAMETER:  8.8264E-02 -1.0900E+00 -1.3461E-01  4.4702E-01 -3.7348E-01  8.2528E-02  1.0106E+00 -5.6098E+00 -2.3439E-01 -1.5833E-01
             5.0538E-01
 GRADIENT:   2.5255E-03  3.3192E+00 -4.7331E+00  1.2438E+01  6.2500E+00 -1.9571E+00  3.3077E-01  0.0000E+00 -1.8638E+00 -8.2511E-01
            -6.2575E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1959.86680488292        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      776
 NPARAMETR:  9.8280E-01  1.1394E-01  9.1499E-01  1.5407E+00  6.4640E-01  9.7959E-01  4.0068E+00  1.0000E-02  6.7878E-01  8.7625E-01
             1.5092E+00
 PARAMETER:  8.2648E-02 -2.0720E+00  1.1160E-02  5.3223E-01 -3.3633E-01  7.9378E-02  1.4880E+00 -1.1659E+01 -2.8746E-01 -3.2104E-02
             5.1159E-01
 GRADIENT:  -4.3630E-01  1.8765E+00 -2.2038E-01  2.2466E+01 -1.4563E+00 -1.5347E+00  2.7545E-01  0.0000E+00 -2.7627E+00  1.5089E+00
            -4.2520E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1960.69907898205        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      952
 NPARAMETR:  9.8014E-01  3.2703E-02  9.2747E-01  1.5801E+00  6.3112E-01  9.8544E-01  6.6523E+00  1.0000E-02  6.6487E-01  8.8579E-01
             1.5220E+00
 PARAMETER:  7.9938E-02 -3.3203E+00  2.4709E-02  5.5751E-01 -3.6026E-01  8.5333E-02  1.9950E+00 -1.9224E+01 -3.0817E-01 -2.1278E-02
             5.2005E-01
 GRADIENT:  -2.6248E+00 -3.9910E+00  1.8618E+01  4.3007E+01 -3.2310E+01  1.5561E+00 -9.9029E+00  0.0000E+00  5.1779E+00  6.8943E+00
             5.6283E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1960.85605164272        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1128
 NPARAMETR:  9.8126E-01  4.1446E-02  9.2711E-01  1.5666E+00  6.3971E-01  9.8173E-01  6.0297E+00  1.0000E-02  6.5963E-01  8.8632E-01
             1.5194E+00
 PARAMETER:  8.1086E-02 -3.0834E+00  2.4313E-02  5.4889E-01 -3.4675E-01  8.1560E-02  1.8967E+00 -1.7870E+01 -3.1607E-01 -2.0673E-02
             5.1832E-01
 GRADIENT:  -5.9273E-02  1.5208E+00 -1.1899E+00 -9.1170E+00  1.6362E+00  2.1918E-01  3.7486E+00  0.0000E+00 -1.9482E-02 -3.1513E+00
             1.1681E+00

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1960.86564534012        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1185
 NPARAMETR:  9.8138E-01  4.1314E-02  9.3594E-01  1.5677E+00  6.4531E-01  9.8079E-01  6.0415E+00  1.0000E-02  6.5683E-01  8.9320E-01
             1.5173E+00
 PARAMETER:  8.1199E-02 -3.0866E+00  3.3796E-02  5.4963E-01 -3.3802E-01  8.0607E-02  1.8986E+00 -1.7938E+01 -3.2033E-01 -1.2949E-02
             5.1695E-01
 GRADIENT:   2.7761E-01 -6.5891E-01  1.8914E-01 -4.3039E-01  1.8303E-01 -2.1088E-01 -1.9542E+00  0.0000E+00 -2.3117E-01  7.2522E-01
             1.2181E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1185
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.7319E-04  2.5666E-02 -1.5974E-04 -1.7807E-02 -1.0294E-02
 SE:             2.9722E-02  1.2161E-02  1.9171E-04  2.7787E-02  2.3267E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7388E-01  3.4805E-02  4.0472E-01  5.2164E-01  6.5820E-01

 ETASHRINKSD(%)  4.2756E-01  5.9261E+01  9.9358E+01  6.9087E+00  2.2051E+01
 ETASHRINKVR(%)  8.5329E-01  8.3403E+01  9.9996E+01  1.3340E+01  3.9240E+01
 EBVSHRINKSD(%)  7.6103E-01  7.4849E+01  9.9328E+01  5.7451E+00  1.8643E+01
 EBVSHRINKVR(%)  1.5163E+00  9.3674E+01  9.9995E+01  1.1160E+01  3.3810E+01
 RELATIVEINF(%)  9.8302E+01  3.8221E+00  4.0010E-04  4.1190E+01  6.1605E+00
 EPSSHRINKSD(%)  3.0170E+01
 EPSSHRINKVR(%)  5.1238E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1960.8656453401190     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1041.9271121354464     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.65
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1960.866       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  4.13E-02  9.36E-01  1.57E+00  6.45E-01  9.81E-01  6.04E+00  1.00E-02  6.57E-01  8.93E-01  1.52E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.18E+03
 
 TH 2
+       -1.28E+01  1.85E+05
 
 TH 3
+       -8.69E+00 -3.48E+03  8.36E+02
 
 TH 4
+       -8.06E+00 -2.78E+04  3.39E+02  3.47E+03
 
 TH 5
+        1.86E+01  5.57E+04 -1.76E+03 -7.99E+03  4.73E+03
 
 TH 6
+        1.34E+00  7.36E+01 -4.11E+00 -2.08E+01  1.13E+01  2.00E+02
 
 TH 7
+       -8.11E-02  2.01E+03 -2.41E+01 -1.87E+02  5.21E+01  6.80E-01  2.70E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.59E+00 -1.43E+03  1.26E+02  1.02E+01 -2.94E+02 -1.20E+00  7.81E-01  0.00E+00  4.42E+02
 
 TH10
+       -3.00E-01 -4.19E+03  3.05E+02  5.67E+02 -9.80E+02 -3.31E+00 -2.54E+01  0.00E+00  1.63E+02  5.10E+02
 
 TH11
+       -1.12E+01 -7.29E+02  4.40E+01  8.87E+01 -1.83E+02  1.80E-01 -3.37E+00  0.00E+00  4.48E+01  1.09E+02  2.08E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       26.165
Stop Time:
Wed Sep 29 22:31:10 CDT 2021
