Sun Oct 24 03:21:07 CDT 2021
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
$DATA ../../../../data/SD4/SL2/dat71.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1589.21603449537        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7485E+02 -3.5756E+01  2.9668E+01 -7.3091E+01 -3.1484E+01  3.5837E+01 -1.3709E+01 -1.9233E+00 -1.0203E+01  2.0702E+00
            -5.1150E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1601.04190439114        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.5629E-01  1.0492E+00  9.9068E-01  1.0548E+00  1.0271E+00  9.9408E-01  1.0577E+00  1.0011E+00  1.0142E+00  9.9875E-01
             1.1330E+00
 PARAMETER:  5.5302E-02  1.4801E-01  9.0637E-02  1.5337E-01  1.2672E-01  9.4063E-02  1.5613E-01  1.0111E-01  1.1406E-01  9.8750E-02
             2.2484E-01
 GRADIENT:   1.5526E+00 -1.5300E+00  4.7276E+00 -1.6651E+00 -1.4311E+01  3.4254E-01 -9.5128E+00  3.8464E-01  1.9295E+00  3.1289E+00
             5.6964E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1602.17810687514        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  9.5741E-01  1.0889E+00  9.4121E-01  1.0309E+00  1.0232E+00  1.0063E+00  1.2487E+00  9.6027E-01  9.3323E-01  9.3193E-01
             1.1325E+00
 PARAMETER:  5.6475E-02  1.8518E-01  3.9410E-02  1.3045E-01  1.2297E-01  1.0629E-01  3.2208E-01  5.9463E-02  3.0897E-02  2.9498E-02
             2.2439E-01
 GRADIENT:   2.9643E+00  8.2063E+00  7.5033E+00  8.0219E-01 -6.6494E+00  4.9074E+00  7.6505E-01  7.5079E-02 -3.5363E+00 -7.2456E-01
             4.9294E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1602.94851737796        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  9.5864E-01  1.1706E+00  6.7536E-01  9.5653E-01  9.0373E-01  9.9279E-01  1.1853E+00  6.1651E-01  9.6781E-01  7.9859E-01
             1.1086E+00
 PARAMETER:  5.7757E-02  2.5750E-01 -2.9252E-01  5.5555E-02 -1.2243E-03  9.2764E-02  2.6998E-01 -3.8368E-01  6.7278E-02 -1.2491E-01
             2.0309E-01
 GRADIENT:   9.4321E-01  4.5849E-01 -1.7634E+00  4.6270E+00 -3.5942E-01 -1.3062E+00  4.9918E-01  6.3392E-01  1.4637E+00  1.7636E+00
            -1.5984E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1603.25751209142        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.5896E-01  1.3561E+00  4.7740E-01  8.1989E-01  8.5911E-01  9.9632E-01  1.0460E+00  3.3694E-01  1.0267E+00  6.6662E-01
             1.1150E+00
 PARAMETER:  5.8093E-02  4.0458E-01 -6.3940E-01 -9.8579E-02 -5.1862E-02  9.6310E-02  1.4496E-01 -9.8786E-01  1.2638E-01 -3.0554E-01
             2.0881E-01
 GRADIENT:   6.6377E-02  5.7444E+00  1.6540E+00  1.7195E+00 -4.3665E+00 -1.1271E-01 -6.8479E-01  2.6618E-01 -2.0436E-02 -2.5589E-01
             7.2595E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1603.29635090905        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  9.5883E-01  1.3931E+00  4.6053E-01  7.9201E-01  8.7406E-01  9.9652E-01  1.0254E+00  2.5246E-01  1.0488E+00  6.8157E-01
             1.1146E+00
 PARAMETER:  5.7954E-02  4.3154E-01 -6.7537E-01 -1.3318E-01 -3.4601E-02  9.6516E-02  1.2511E-01 -1.2765E+00  1.4762E-01 -2.8336E-01
             2.0847E-01
 GRADIENT:  -7.8683E-02 -1.7710E+00 -9.7634E-01 -7.6974E-01  1.3791E+00 -6.7860E-03  7.0738E-02  1.5853E-01  1.7382E-01  2.0961E-01
             1.3739E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1603.37572294656        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  9.5877E-01  1.3489E+00  4.5288E-01  8.1641E-01  8.3920E-01  9.9642E-01  1.0498E+00  9.5444E-02  1.0210E+00  6.4920E-01
             1.1138E+00
 PARAMETER:  5.7894E-02  3.9929E-01 -6.9212E-01 -1.0283E-01 -7.5310E-02  9.6409E-02  1.4859E-01 -2.2492E+00  1.2083E-01 -3.3201E-01
             2.0780E-01
 GRADIENT:   1.9281E-02  5.2563E-01  8.8860E-02  2.5581E-01 -4.5006E-01 -2.7190E-02  9.7947E-02  1.6126E-02 -2.2676E-02  7.1521E-02
             3.3359E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1603.38166777661        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  9.5874E-01  1.3577E+00  4.4753E-01  8.1014E-01  8.4046E-01  9.9644E-01  1.0437E+00  1.9746E-02  1.0257E+00  6.4776E-01
             1.1139E+00
 PARAMETER:  5.7864E-02  4.0583E-01 -7.0401E-01 -1.1055E-01 -7.3809E-02  9.6438E-02  1.4279E-01 -3.8248E+00  1.2534E-01 -3.3423E-01
             2.0784E-01
 GRADIENT:   2.8165E-02  2.1421E-01  3.7957E-02  4.9710E-02 -1.6504E-01 -3.6767E-03  5.5294E-02  7.5660E-04 -2.5099E-03  3.0816E-02
             3.4771E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1603.38201705969        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1399
 NPARAMETR:  9.5872E-01  1.3575E+00  4.4721E-01  8.1006E-01  8.4019E-01  9.9646E-01  1.0434E+00  1.0000E-02  1.0257E+00  6.4722E-01
             1.1138E+00
 PARAMETER:  5.7845E-02  4.0563E-01 -7.0472E-01 -1.1065E-01 -7.4132E-02  9.6449E-02  1.4247E-01 -4.5150E+00  1.2538E-01 -3.3507E-01
             2.0774E-01
 GRADIENT:  -1.0221E+00  3.1227E-01  9.6461E-03  2.4326E-03  1.8288E-02 -9.8561E-02 -2.1254E-02  3.4934E-05  3.8246E-03  1.3581E-02
            -3.3859E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1399
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.3824E-04 -9.6248E-03 -3.5731E-04  7.4240E-03 -1.8648E-02
 SE:             2.9827E-02  2.5475E-02  1.5506E-04  2.4586E-02  1.8199E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8025E-01  7.0557E-01  2.1205E-02  7.6268E-01  3.0552E-01

 ETASHRINKSD(%)  7.4280E-02  1.4656E+01  9.9481E+01  1.7635E+01  3.9030E+01
 ETASHRINKVR(%)  1.4851E-01  2.7165E+01  9.9997E+01  3.2160E+01  6.2827E+01
 EBVSHRINKSD(%)  5.1541E-01  1.4326E+01  9.9530E+01  1.7907E+01  3.9772E+01
 EBVSHRINKVR(%)  1.0282E+00  2.6599E+01  9.9998E+01  3.2607E+01  6.3725E+01
 RELATIVEINF(%)  9.8883E+01  5.8751E+00  1.7845E-04  5.2008E+00  3.1339E+00
 EPSSHRINKSD(%)  4.3758E+01
 EPSSHRINKVR(%)  6.8368E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1603.3820170596853     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -868.23119049594709     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1603.382       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.59E-01  1.36E+00  4.47E-01  8.10E-01  8.40E-01  9.96E-01  1.04E+00  1.00E-02  1.03E+00  6.47E-01  1.11E+00
 


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
 #CPUT: Total CPU Time in Seconds,       38.444
Stop Time:
Sun Oct 24 03:21:16 CDT 2021
