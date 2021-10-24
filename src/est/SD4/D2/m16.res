Sun Oct 24 04:31:15 CDT 2021
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
$DATA ../../../../data/SD4/D2/dat16.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m16.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1517.48259259434        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2129E+02 -9.3076E+01 -6.4697E+01 -2.1739E+01  1.3152E+02  2.6866E+01 -9.6063E+01 -5.2573E+00 -8.6408E+01 -4.1061E+01
            -7.9947E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1574.72047421473        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0512E+00  1.2053E+00  1.1496E+00  8.7770E-01  1.1238E+00  1.0089E+00  2.0669E+00  9.7941E-01  1.3028E+00  1.0412E+00
             1.2083E+00
 PARAMETER:  1.4994E-01  2.8672E-01  2.3943E-01 -3.0449E-02  2.1669E-01  1.0883E-01  8.2603E-01  7.9191E-02  3.6455E-01  1.4034E-01
             2.8924E-01
 GRADIENT:   4.5069E+02  7.8067E+01 -5.3730E+00  1.4101E+01  3.7579E+01  1.9564E+01  1.2664E+02  1.4192E+00  4.1366E+01 -1.2494E+01
             1.4772E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1577.70577574938        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      266
 NPARAMETR:  1.0698E+00  1.2715E+00  1.4933E+00  9.5236E-01  1.2867E+00  1.1025E+00  1.8166E+00  9.1263E-01  1.3061E+00  1.3562E+00
             1.2503E+00
 PARAMETER:  1.6746E-01  3.4018E-01  5.0096E-01  5.1189E-02  3.5211E-01  1.9756E-01  6.9697E-01  8.5714E-03  3.6705E-01  4.0467E-01
             3.2339E-01
 GRADIENT:   3.5080E+01  7.0835E+00 -1.6597E+01  2.6310E+01  4.2763E+01  2.0324E+01  1.7138E+01  1.4810E+00  1.5017E+01 -3.3455E-01
             2.0985E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1582.39822015961        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.0451E+00  1.3030E+00  1.3288E+00  9.0405E-01  1.1862E+00  1.0455E+00  1.6670E+00  6.4308E-01  1.1778E+00  1.3045E+00
             1.2031E+00
 PARAMETER:  1.4410E-01  3.6469E-01  3.8429E-01 -8.7360E-04  2.7071E-01  1.4454E-01  6.1105E-01 -3.4148E-01  2.6362E-01  3.6583E-01
             2.8493E-01
 GRADIENT:  -8.5281E+00  5.6601E+00  1.4339E+00  6.9769E+00 -4.6251E+00  1.1694E+00 -2.6233E-01  2.0810E-01  2.2347E+00  2.6531E+00
             3.7593E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1584.32970950471        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      625
 NPARAMETR:  1.0554E+00  1.6465E+00  6.6992E-01  6.5317E-01  1.1093E+00  1.0404E+00  1.4794E+00  4.9046E-01  1.0490E+00  1.1441E+00
             1.1615E+00
 PARAMETER:  1.5395E-01  5.9864E-01 -3.0060E-01 -3.2592E-01  2.0369E-01  1.3965E-01  4.9162E-01 -6.1241E-01  1.4780E-01  2.3458E-01
             2.4972E-01
 GRADIENT:   4.6961E+00 -6.3588E-01 -4.4410E-01  2.2899E+00 -7.1915E-01 -9.5862E-01 -7.0363E-01  1.5646E-01 -3.8364E-01 -2.8942E-01
            -2.0418E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1584.36532311294        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  1.0541E+00  1.7652E+00  5.9414E-01  5.8417E-01  1.1393E+00  1.0444E+00  1.4091E+00  2.4649E-01  1.0812E+00  1.1758E+00
             1.1703E+00
 PARAMETER:  1.5265E-01  6.6828E-01 -4.2064E-01 -4.3756E-01  2.3045E-01  1.4343E-01  4.4292E-01 -1.3004E+00  1.7804E-01  2.6197E-01
             2.5730E-01
 GRADIENT:   1.1129E+00  3.9634E+00  1.3685E-01  3.2662E+00 -8.2276E-01  6.7289E-01 -5.7143E-01  2.2575E-02 -9.8120E-01  3.9343E-01
             3.1436E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1584.41654513756        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      983
 NPARAMETR:  1.0541E+00  1.7643E+00  5.8806E-01  5.7661E-01  1.1416E+00  1.0431E+00  1.4055E+00  1.4795E-01  1.1112E+00  1.1732E+00
             1.1697E+00
 PARAMETER:  1.5272E-01  6.6774E-01 -4.3093E-01 -4.5059E-01  2.3245E-01  1.4221E-01  4.4037E-01 -1.8109E+00  2.0547E-01  2.5974E-01
             2.5673E-01
 GRADIENT:   1.3761E+00 -2.9395E+00  2.1299E-01 -5.6601E-01 -5.4006E-01  1.7317E-01  5.3221E-01  1.0284E-02 -1.1629E-01  3.0541E-01
             3.3723E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1584.42399288861        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1165
 NPARAMETR:  1.0549E+00  1.7626E+00  5.8715E-01  5.7674E-01  1.1413E+00  1.0431E+00  1.4039E+00  2.0106E-02  1.1184E+00  1.1712E+00
             1.1689E+00
 PARAMETER:  1.5343E-01  6.6680E-01 -4.3247E-01 -4.5036E-01  2.3217E-01  1.4224E-01  4.3922E-01 -3.8067E+00  2.1189E-01  2.5802E-01
             2.5604E-01
 GRADIENT:   2.8500E+00 -3.8636E+00  5.8869E-02 -5.9992E-01  1.2218E-01  1.6822E-01  3.2288E-01  2.0678E-04  3.9292E-02  1.6321E-01
             6.9276E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1584.42455685034        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1261
 NPARAMETR:  1.0549E+00  1.7613E+00  5.8669E-01  5.7730E-01  1.1409E+00  1.0432E+00  1.4038E+00  1.0000E-02  1.1170E+00  1.1688E+00
             1.1685E+00
 PARAMETER:  1.5343E-01  6.6608E-01 -4.3326E-01 -4.4940E-01  2.3186E-01  1.4225E-01  4.3915E-01 -5.5573E+00  2.1064E-01  2.5601E-01
             2.5570E-01
 GRADIENT:   2.8460E+00 -4.1207E+00 -1.0653E-01 -3.9933E-01  8.3201E-01  1.7041E-01  1.6993E-01  0.0000E+00 -1.8461E-02 -1.0152E-01
            -1.0388E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1261
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.6959E-04 -7.9313E-03 -2.6718E-04  7.1020E-03 -2.8020E-02
 SE:             2.9812E-02  2.7481E-02  8.9648E-05  1.7045E-02  2.2475E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8476E-01  7.7288E-01  2.8795E-03  6.7692E-01  2.1250E-01

 ETASHRINKSD(%)  1.2739E-01  7.9352E+00  9.9700E+01  4.2899E+01  2.4706E+01
 ETASHRINKVR(%)  2.5463E-01  1.5241E+01  9.9999E+01  6.7394E+01  4.3308E+01
 EBVSHRINKSD(%)  5.1689E-01  7.7840E+00  9.9753E+01  4.6898E+01  2.1665E+01
 EBVSHRINKVR(%)  1.0311E+00  1.4962E+01  9.9999E+01  7.1801E+01  3.8637E+01
 RELATIVEINF(%)  9.8756E+01  9.8280E+00  6.1903E-05  1.8588E+00  1.3956E+01
 EPSSHRINKSD(%)  4.2956E+01
 EPSSHRINKVR(%)  6.7460E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1584.4245568503411     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -849.27373028660293     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1584.425       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.76E+00  5.87E-01  5.77E-01  1.14E+00  1.04E+00  1.40E+00  1.00E-02  1.12E+00  1.17E+00  1.17E+00
 


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
 #CPUT: Total CPU Time in Seconds,       37.493
Stop Time:
Sun Oct 24 04:31:24 CDT 2021
