Sat Oct 23 22:20:45 CDT 2021
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
$DATA ../../../../data/SD3/A2/dat81.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1001.85042324423        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7174E+02 -1.2757E+01  1.6482E+02 -3.2063E+01  1.8617E+02  2.6361E+01 -2.1131E+01 -3.2640E+02 -7.1788E+01 -7.6421E+01
            -1.6930E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1637.58524323892        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0290E+00  1.1220E+00  9.0810E-01  1.0554E+00  9.6220E-01  1.1034E+00  9.4973E-01  8.1856E-01  6.8653E-01  7.5633E-01
             2.2288E+00
 PARAMETER:  1.2855E-01  2.1509E-01  3.6027E-03  1.5391E-01  6.1466E-02  1.9839E-01  4.8422E-02 -1.0021E-01 -2.7611E-01 -1.7928E-01
             9.0148E-01
 GRADIENT:   9.7250E+01  9.4828E+01 -3.0803E+00  1.4722E+02  5.0373E+01  4.8047E+01 -2.3022E+01 -6.8246E+00 -4.5479E+01 -1.1711E+01
            -2.5897E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1653.13473980107        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0287E+00  1.1128E+00  2.9097E-01  1.0080E+00  4.9919E-01  1.0731E+00  1.4909E+00  3.9067E-01  7.0954E-01  3.3017E-01
             2.2814E+00
 PARAMETER:  1.2832E-01  2.0691E-01 -1.1345E+00  1.0798E-01 -5.9476E-01  1.7053E-01  4.9935E-01 -8.3989E-01 -2.4314E-01 -1.0082E+00
             9.2479E-01
 GRADIENT:   6.3188E+01  1.4713E+02 -4.0968E+01  1.6764E+02  3.2738E+01  2.1250E+01  4.9677E+01 -7.4261E+00 -5.2293E+01 -9.5550E-01
            -1.3004E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1690.52442931033        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0197E+00  7.1712E-01  2.5341E-01  1.0767E+00  3.7090E-01  9.6136E-01  1.2258E+00  5.0256E-01  8.6131E-01  1.6472E-01
             2.4343E+00
 PARAMETER:  1.1950E-01 -2.3251E-01 -1.2728E+00  1.7391E-01 -8.9182E-01  6.0594E-02  3.0356E-01 -5.8803E-01 -4.9304E-02 -1.7035E+00
             9.8965E-01
 GRADIENT:   1.5787E+01  5.0440E+01 -3.1135E+00  7.5836E+01  1.7536E+01 -3.2524E+01  3.5071E+00 -1.1593E+01 -4.2071E+00 -5.6127E-01
            -2.0332E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1699.47286735316        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      411
 NPARAMETR:  1.0658E+00  6.0075E-01  3.0747E-01  1.1359E+00  3.7159E-01  1.0257E+00  1.2164E+00  7.9392E-01  8.1000E-01  1.5036E-01
             2.4747E+00
 PARAMETER:  1.6369E-01 -4.0957E-01 -1.0794E+00  2.2738E-01 -8.8996E-01  1.2540E-01  2.9591E-01 -1.3077E-01 -1.1073E-01 -1.7947E+00
             1.0061E+00
 GRADIENT:   4.1789E+01  1.2514E+01  3.7049E+01  2.0654E+01 -4.9498E+01 -9.1035E+00 -7.7972E+00 -1.7403E+01 -1.0157E+01 -4.4196E-01
             2.9379E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1710.56779976381        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  1.0223E+00  6.7585E-01  2.8734E-01  1.1136E+00  3.6995E-01  1.0377E+00  8.1021E-01  1.5828E+00  9.3401E-01  1.1319E-01
             2.3174E+00
 PARAMETER:  1.2205E-01 -2.9178E-01 -1.1471E+00  2.0757E-01 -8.9438E-01  1.3698E-01 -1.1046E-01  5.5919E-01  3.1732E-02 -2.0786E+00
             9.4046E-01
 GRADIENT:  -4.1913E+01  5.8239E+01  3.7215E+01  4.9546E+01 -8.3226E+01 -6.7808E+00 -1.1183E+01  2.7753E+00 -1.8821E+01  1.0075E-01
             5.7520E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1719.26476686046        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      763
 NPARAMETR:  1.0390E+00  5.5127E-01  2.4825E-01  1.1161E+00  3.2072E-01  1.0454E+00  7.3855E-01  1.6261E+00  1.0524E+00  6.6047E-02
             2.0338E+00
 PARAMETER:  1.3829E-01 -4.9553E-01 -1.2933E+00  2.0986E-01 -1.0372E+00  1.4443E-01 -2.0307E-01  5.8619E-01  1.5108E-01 -2.6174E+00
             8.0989E-01
 GRADIENT:   9.9388E-01  1.8591E+01  2.4209E+01  2.2342E+01 -4.3608E+01 -4.0438E+00 -1.2005E+01 -6.5370E+00 -7.8202E+00 -2.4360E-01
             4.1718E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1722.31400799211        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      939
 NPARAMETR:  1.0389E+00  5.6104E-01  2.4053E-01  1.1012E+00  3.2372E-01  1.0576E+00  1.0360E+00  1.6236E+00  1.0487E+00  1.3332E-02
             1.9824E+00
 PARAMETER:  1.3818E-01 -4.7797E-01 -1.3249E+00  1.9642E-01 -1.0279E+00  1.5604E-01  1.3535E-01  5.8463E-01  1.4750E-01 -4.2176E+00
             7.8430E-01
 GRADIENT:   2.6519E+00  1.9103E+00  9.3630E-01  2.5734E+00 -1.0365E+00  1.1066E-01 -8.5629E-02 -7.9640E-01  1.1106E-01 -3.0565E-03
            -2.5184E+00

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1722.33623308663        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1031
 NPARAMETR:  1.0375E+00  5.5532E-01  2.3893E-01  1.0997E+00  3.2143E-01  1.0575E+00  1.0379E+00  1.6255E+00  1.0490E+00  1.2601E-02
             1.9876E+00
 PARAMETER:  1.3683E-01 -4.8821E-01 -1.3316E+00  1.9504E-01 -1.0350E+00  1.5595E-01  1.3720E-01  5.8579E-01  1.4787E-01 -4.2740E+00
             7.8694E-01
 GRADIENT:  -5.1950E-02  2.7733E-01  4.4344E-01 -3.1561E-01 -6.6361E-01  9.6795E-02  6.5335E-02 -3.4133E-03  7.0401E-04 -2.9405E-03
            -1.0703E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1031
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0661E-03 -1.3372E-03 -6.3682E-03 -2.9794E-03  3.7756E-05
 SE:             2.9621E-02  1.8049E-02  2.5473E-02  2.7465E-02  4.4847E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7129E-01  9.4094E-01  8.0259E-01  9.1361E-01  9.3291E-01

 ETASHRINKSD(%)  7.6672E-01  3.9533E+01  1.4663E+01  7.9888E+00  9.8498E+01
 ETASHRINKVR(%)  1.5276E+00  6.3437E+01  2.7176E+01  1.5339E+01  9.9977E+01
 EBVSHRINKSD(%)  1.1166E+00  3.9006E+01  1.4979E+01  8.4598E+00  9.8615E+01
 EBVSHRINKVR(%)  2.2207E+00  6.2797E+01  2.7714E+01  1.6204E+01  9.9981E+01
 RELATIVEINF(%)  9.7730E+01  8.2498E+00  2.5049E+01  7.0767E+01  2.5614E-03
 EPSSHRINKSD(%)  3.4259E+01
 EPSSHRINKVR(%)  5.6781E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1722.3362330866335     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -803.39769988196076     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.47
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1722.336       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  5.55E-01  2.39E-01  1.10E+00  3.21E-01  1.06E+00  1.04E+00  1.63E+00  1.05E+00  1.26E-02  1.99E+00
 


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
 #CPUT: Total CPU Time in Seconds,       81.848
Stop Time:
Sat Oct 23 22:20:58 CDT 2021
