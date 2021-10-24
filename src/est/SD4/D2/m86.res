Sun Oct 24 04:44:32 CDT 2021
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
$DATA ../../../../data/SD4/D2/dat86.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1542.90543575895        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9452E+02  1.2893E+00 -3.9615E+00  2.1917E+01  4.3337E+01 -4.8727E+00 -6.8609E+01 -1.1838E+01 -9.4413E+01 -1.6099E+01
             6.3082E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1575.04316505218        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0091E+00  1.0508E+00  1.0169E+00  8.6693E-01  1.0728E+00  1.1443E+00  1.7812E+00  1.1087E+00  1.4220E+00  1.0305E+00
             9.2016E-01
 PARAMETER:  1.0909E-01  1.4954E-01  1.1677E-01 -4.2801E-02  1.7023E-01  2.3480E-01  6.7731E-01  2.0318E-01  4.5208E-01  1.3004E-01
             1.6791E-02
 GRADIENT:  -2.0036E+01 -5.8972E+01  1.6286E+00 -7.0799E+01  1.5946E+01 -1.0922E+01  9.4088E+00 -1.3058E+00  1.1762E+01  5.4743E+00
            -2.1531E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1577.29573154518        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0254E+00  1.3486E+00  9.0331E-01  7.3370E-01  1.1051E+00  1.1753E+00  1.5775E+00  1.6172E+00  1.5097E+00  9.4918E-01
             9.1400E-01
 PARAMETER:  1.2504E-01  3.9903E-01 -1.6893E-03 -2.0966E-01  1.9996E-01  2.6152E-01  5.5585E-01  5.8067E-01  5.1192E-01  4.7839E-02
             1.0076E-02
 GRADIENT:   4.8747E+00 -3.2594E+01  7.7058E-01 -4.3077E+01 -1.5296E+01 -6.1642E-02  2.4428E+01  9.4383E+00  5.6772E+00 -3.5008E+00
            -2.6222E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1577.65790717161        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      556
 NPARAMETR:  1.0257E+00  1.3765E+00  9.0103E-01  7.2788E-01  1.1099E+00  1.1702E+00  1.5624E+00  1.6540E+00  1.5183E+00  9.5550E-01
             9.1547E-01
 PARAMETER:  1.2534E-01  4.1954E-01 -4.2222E-03 -2.1762E-01  2.0429E-01  2.5716E-01  5.4625E-01  6.0317E-01  5.1757E-01  5.4475E-02
             1.1683E-02
 GRADIENT:   5.2252E+00 -2.5697E+01  1.4434E+00 -3.7810E+01 -1.8657E+01 -1.8768E+00  2.5248E+01  9.3972E+00  5.7341E+00 -3.1893E+00
            -2.5588E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1579.55854938610        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  1.0255E+00  1.3759E+00  9.0093E-01  7.2789E-01  1.1187E+00  1.1651E+00  1.4598E+00  1.6530E+00  1.5183E+00  9.6827E-01
             9.4017E-01
 PARAMETER:  1.2521E-01  4.1912E-01 -4.3228E-03 -2.1761E-01  2.1217E-01  2.5282E-01  4.7829E-01  6.0256E-01  5.1759E-01  6.7751E-02
             3.8307E-02
 GRADIENT:   4.5398E+00 -3.5733E+01 -6.6102E-01 -3.4806E+01 -1.1549E+01 -3.5500E+00  1.2190E+01  8.8796E+00  4.5891E+00 -2.4151E+00
            -1.3700E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1580.37649595977        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  1.0247E+00  1.3799E+00  9.1272E-01  7.3382E-01  1.1198E+00  1.1905E+00  1.3206E+00  1.6518E+00  1.5182E+00  9.9961E-01
             9.9156E-01
 PARAMETER:  1.2444E-01  4.2203E-01  8.6741E-03 -2.0949E-01  2.1311E-01  2.7438E-01  3.7809E-01  6.0186E-01  5.1752E-01  9.9614E-02
             9.1523E-02
 GRADIENT:   2.4266E+00 -4.0956E+01 -7.8419E-02 -3.1454E+01 -1.6555E+01  5.5394E+00 -5.9656E+00  7.8052E+00  2.0811E+00  1.4856E+00
             8.6772E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1581.15786566501        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1092
 NPARAMETR:  1.0255E+00  1.3803E+00  9.0898E-01  7.5049E-01  1.1199E+00  1.1774E+00  1.3674E+00  1.6512E+00  1.5187E+00  9.8925E-01
             9.7053E-01
 PARAMETER:  1.2521E-01  4.2229E-01  4.5639E-03 -1.8702E-01  2.1325E-01  2.6330E-01  4.1290E-01  6.0149E-01  5.1784E-01  8.9196E-02
             7.0090E-02
 GRADIENT:   3.8607E+00 -2.8759E+01 -3.7691E+00 -1.8368E+01 -8.8870E+00  9.1757E-01  1.3522E-01  8.6388E+00  5.6232E+00 -2.3663E-01
            -2.0169E-01

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1581.15786566501        NO. OF FUNC. EVALS.:  26
 CUMULATIVE NO. OF FUNC. EVALS.:     1118
 NPARAMETR:  1.0255E+00  1.3803E+00  9.0898E-01  7.5049E-01  1.1199E+00  1.1774E+00  1.3674E+00  1.6512E+00  1.5187E+00  9.8925E-01
             9.7053E-01
 PARAMETER:  1.2521E-01  4.2229E-01  4.5639E-03 -1.8702E-01  2.1325E-01  2.6330E-01  4.1290E-01  6.0149E-01  5.1784E-01  8.9196E-02
             7.0090E-02
 GRADIENT:   1.9206E+00 -3.2846E+01  3.4518E+01 -3.8931E+01 -4.8649E+00  6.7217E-01  4.4901E-01 -1.6524E+02 -3.0592E+01  6.2134E-01
            -1.2693E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1118
 NO. OF SIG. DIGITS IN FINAL EST.:  4.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.1312E-05 -2.0973E-03 -3.9430E-02  2.8365E-02 -4.6097E-02
 SE:             2.9950E-02  2.3375E-02  1.2963E-02  2.2070E-02  1.9636E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9970E-01  9.2851E-01  2.3519E-03  1.9871E-01  1.8893E-02

 ETASHRINKSD(%)  1.0000E-10  2.1690E+01  5.6573E+01  2.6061E+01  3.4219E+01
 ETASHRINKVR(%)  1.0000E-10  3.8676E+01  8.1141E+01  4.5331E+01  5.6728E+01
 EBVSHRINKSD(%)  3.0593E-01  2.1515E+01  5.4565E+01  2.4885E+01  3.1655E+01
 EBVSHRINKVR(%)  6.1091E-01  3.8401E+01  7.9357E+01  4.3578E+01  5.3290E+01
 RELATIVEINF(%)  9.9214E+01  4.7151E+00  2.5059E+00  4.3280E+00  1.3896E+01
 EPSSHRINKSD(%)  4.6590E+01
 EPSSHRINKVR(%)  7.1473E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1581.1578656650099     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -846.00703910127174     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.61
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1581.158       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.38E+00  9.09E-01  7.50E-01  1.12E+00  1.18E+00  1.37E+00  1.65E+00  1.52E+00  9.89E-01  9.71E-01
 


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
 #CPUT: Total CPU Time in Seconds,       42.153
Stop Time:
Sun Oct 24 04:44:42 CDT 2021
