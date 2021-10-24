Sun Oct 24 00:41:00 CDT 2021
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
$DATA ../../../../data/SD3/TD2/dat18.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m18.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2050.77070721721        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2849E+02 -4.6216E+01 -2.1851E+01 -7.2657E+00  6.7888E+01  4.3223E+01 -2.5437E+01  3.2744E+00 -2.4136E+01 -2.2446E+01
            -5.8527E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2060.76017696285        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9952E-01  1.0567E+00  1.0001E+00  1.0275E+00  9.9898E-01  9.8304E-01  1.1079E+00  9.8796E-01  1.0945E+00  1.0691E+00
             1.0718E+00
 PARAMETER:  9.9523E-02  1.5517E-01  1.0014E-01  1.2713E-01  9.8981E-02  8.2897E-02  2.0242E-01  8.7888E-02  1.9027E-01  1.6682E-01
             1.6938E-01
 GRADIENT:   7.4609E+00 -1.8513E+01 -1.5672E+01  9.4108E+00  2.3268E+01 -1.0770E+00 -1.1653E+01  4.4281E+00  1.9308E+00 -6.3557E+00
             2.4495E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2062.69122379005        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  1.0014E+00  1.1522E+00  9.9346E-01  9.6522E-01  1.0410E+00  9.9268E-01  1.2899E+00  9.1655E-01  1.0528E+00  1.1608E+00
             1.0545E+00
 PARAMETER:  1.0137E-01  2.4166E-01  9.3437E-02  6.4601E-02  1.4020E-01  9.2653E-02  3.5457E-01  1.2865E-02  1.5144E-01  2.4908E-01
             1.5310E-01
 GRADIENT:   1.1349E+01 -2.0065E+00 -1.3096E+00 -2.1726E+00  4.7203E+00  2.5034E+00  2.6234E+00  1.4403E+00  3.1506E+00  2.4804E+00
            -1.0013E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2063.19498956178        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.9549E-01  1.0464E+00  9.4468E-01  1.0252E+00  9.6115E-01  9.8531E-01  1.3715E+00  7.3774E-01  9.9209E-01  1.0733E+00
             1.0698E+00
 PARAMETER:  9.5483E-02  1.4537E-01  4.3096E-02  1.2485E-01  6.0376E-02  8.5206E-02  4.1588E-01 -2.0416E-01  9.2054E-02  1.7071E-01
             1.6743E-01
 GRADIENT:  -1.5026E+00 -1.8819E+00 -2.6961E-01 -2.1091E+00  1.9565E+00 -1.1307E-01 -1.7370E-01 -2.7462E-01 -6.8589E-02 -3.5905E-01
             1.5033E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2063.29418882519        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      711
 NPARAMETR:  9.9319E-01  8.6871E-01  1.1390E+00  1.1537E+00  9.7010E-01  9.8215E-01  1.6419E+00  9.6358E-01  8.7789E-01  1.1140E+00
             1.0735E+00
 PARAMETER:  9.3163E-02 -4.0748E-02  2.3012E-01  2.4298E-01  6.9641E-02  8.1984E-02  5.9584E-01  6.2901E-02 -3.0229E-02  2.0800E-01
             1.7093E-01
 GRADIENT:  -1.5332E+00  8.1478E+00  2.4836E+00  1.0519E+01 -5.0378E+00 -4.5500E-01 -4.7922E-01  2.7608E-02 -5.4211E-01  3.4298E-02
             1.2223E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2063.73102783204        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      889
 NPARAMETR:  9.8923E-01  5.7859E-01  1.4274E+00  1.3390E+00  9.8634E-01  9.7901E-01  2.2573E+00  1.2321E+00  7.5567E-01  1.1671E+00
             1.0739E+00
 PARAMETER:  8.9173E-02 -4.4716E-01  4.5582E-01  3.9192E-01  8.6247E-02  7.8789E-02  9.1418E-01  3.0871E-01 -1.8015E-01  2.5453E-01
             1.7127E-01
 GRADIENT:  -1.6037E+00  6.2963E+00  4.0766E-01  1.5476E+01 -1.8890E+00 -3.5949E-01 -2.0562E-01  4.6022E-01 -2.2724E+00  1.5768E-01
             7.6193E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2063.99376031680        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1066
 NPARAMETR:  9.8846E-01  4.0962E-01  1.5233E+00  1.4372E+00  9.6489E-01  9.7909E-01  2.7815E+00  1.3035E+00  7.3559E-01  1.1627E+00
             1.0734E+00
 PARAMETER:  8.8391E-02 -7.9253E-01  5.2090E-01  4.6271E-01  6.4262E-02  7.8873E-02  1.1230E+00  3.6504E-01 -2.0708E-01  2.5076E-01
             1.7086E-01
 GRADIENT:   2.1496E+00  2.7742E+00  7.7988E-01  6.0970E+00 -4.3786E+00  5.0191E-01  8.6976E-01  2.0868E-01 -6.5133E-01  4.6367E-01
            -9.8121E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2064.00619987759        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1241
 NPARAMETR:  9.8691E-01  3.5275E-01  1.5942E+00  1.4742E+00  9.7293E-01  9.7746E-01  2.9995E+00  1.3614E+00  7.3043E-01  1.1742E+00
             1.0746E+00
 PARAMETER:  8.6824E-02 -9.4200E-01  5.6638E-01  4.8810E-01  7.2557E-02  7.7205E-02  1.1984E+00  4.0853E-01 -2.1412E-01  2.6063E-01
             1.7193E-01
 GRADIENT:   4.3317E-01  1.8777E+00  5.6512E-01  5.9092E+00 -2.4246E+00  1.1011E-01  4.8724E-01 -7.7602E-02 -5.3563E-01  3.1619E-01
            -1.8729E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2064.04609904634        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1423
 NPARAMETR:  9.8710E-01  3.4715E-01  1.5974E+00  1.4714E+00  9.7468E-01  9.7729E-01  3.0421E+00  1.3647E+00  7.3001E-01  1.1744E+00
             1.0746E+00
 PARAMETER:  8.7014E-02 -9.5801E-01  5.6840E-01  4.8621E-01  7.4356E-02  7.7031E-02  1.2125E+00  4.1090E-01 -2.1469E-01  2.6072E-01
             1.7199E-01
 GRADIENT:   1.1080E+00  4.8441E-01  3.6565E-01 -6.5548E+00 -5.2517E-01  9.0440E-02  1.7202E+00  6.2547E-03  1.7802E-01 -2.7306E-02
             4.3996E-03

0ITERATION NO.:   41    OBJECTIVE VALUE:  -2064.04609904634        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     1447
 NPARAMETR:  9.8710E-01  3.4715E-01  1.5974E+00  1.4714E+00  9.7468E-01  9.7729E-01  3.0421E+00  1.3647E+00  7.3001E-01  1.1744E+00
             1.0746E+00
 PARAMETER:  8.7014E-02 -9.5801E-01  5.6840E-01  4.8621E-01  7.4356E-02  7.7031E-02  1.2125E+00  4.1090E-01 -2.1469E-01  2.6072E-01
             1.7199E-01
 GRADIENT:  -2.0399E-01  1.5982E-01  4.0885E-01  1.6412E+00 -2.8478E-01  1.8225E-02 -1.2630E-01 -2.1730E-03  1.2827E-01 -3.9404E-02
             5.5635E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1447
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6700E-03  3.5387E-02 -4.2140E-02 -3.4633E-02 -3.0635E-02
 SE:             2.9845E-02  1.9003E-02  1.7440E-02  2.3009E-02  2.0939E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5538E-01  6.2575E-02  1.5678E-02  1.3227E-01  1.4344E-01

 ETASHRINKSD(%)  1.5738E-02  3.6337E+01  4.1575E+01  2.2918E+01  2.9853E+01
 ETASHRINKVR(%)  3.1473E-02  5.9470E+01  6.5865E+01  4.0584E+01  5.0794E+01
 EBVSHRINKSD(%)  4.6123E-01  4.1712E+01  4.6685E+01  1.8154E+01  2.5179E+01
 EBVSHRINKVR(%)  9.2033E-01  6.6025E+01  7.1575E+01  3.3012E+01  4.4018E+01
 RELATIVEINF(%)  9.8623E+01  6.6204E+00  8.6849E+00  1.3088E+01  1.8615E+01
 EPSSHRINKSD(%)  3.4289E+01
 EPSSHRINKVR(%)  5.6821E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2064.0460990463412     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1145.1075658416685     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2064.046       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.87E-01  3.47E-01  1.60E+00  1.47E+00  9.75E-01  9.77E-01  3.04E+00  1.36E+00  7.30E-01  1.17E+00  1.07E+00
 


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
 #CPUT: Total CPU Time in Seconds,       48.849
Stop Time:
Sun Oct 24 00:41:10 CDT 2021
