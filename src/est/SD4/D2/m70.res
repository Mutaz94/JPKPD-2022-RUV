Sun Oct 24 04:41:15 CDT 2021
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
$DATA ../../../../data/SD4/D2/dat70.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m70.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1519.38761265459        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.9204E+02 -4.4006E+01 -2.7015E+01 -3.8683E+01  5.9842E+01 -6.8044E+01 -5.7833E+01 -1.0642E+01 -8.9417E+01 -3.2612E-02
            -6.0465E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1562.32142264059        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      191
 NPARAMETR:  1.0442E+00  1.1173E+00  1.3486E+00  1.0200E+00  1.1154E+00  1.4513E+00  1.9942E+00  1.2180E+00  1.7636E+00  8.6057E-01
             1.1224E+00
 PARAMETER:  1.4321E-01  2.1095E-01  3.9908E-01  1.1982E-01  2.0924E-01  4.7244E-01  7.9025E-01  2.9721E-01  6.6736E-01 -5.0156E-02
             2.1550E-01
 GRADIENT:  -2.9037E+01  9.5142E+00  1.9350E+01  3.7436E+00 -1.6620E+01  3.9633E+01  2.9860E+01 -1.4773E+01  6.0070E+01 -7.1714E+00
            -4.8311E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1567.49903014426        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.0337E+00  6.9630E-01  2.8313E+00  1.3901E+00  1.2317E+00  1.5757E+00  2.1505E+00  2.8442E+00  1.2512E+00  8.9943E-01
             1.1790E+00
 PARAMETER:  1.3315E-01 -2.6198E-01  1.1407E+00  4.2938E-01  3.0842E-01  5.5473E-01  8.6568E-01  1.1453E+00  3.2412E-01 -5.9943E-03
             2.6469E-01
 GRADIENT:  -3.1453E+01  2.4974E+01  5.6391E+00  4.1114E+01 -1.0623E+01  6.7423E+01  1.2250E+01  1.6295E+01  2.0025E+01 -1.2432E+01
             1.1663E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1585.02135135019        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0726E+00  7.6694E-01  1.4886E+00  1.2444E+00  1.0481E+00  1.3000E+00  2.2393E+00  1.5012E+00  9.6150E-01  8.7499E-01
             1.1054E+00
 PARAMETER:  1.7005E-01 -1.6534E-01  4.9783E-01  3.1862E-01  1.4696E-01  3.6235E-01  9.0614E-01  5.0624E-01  6.0738E-02 -3.3540E-02
             2.0024E-01
 GRADIENT:   2.3454E+00  4.5901E+00  9.2046E+00  1.2381E+01 -6.8570E+00 -1.1197E+00 -4.4074E+00 -3.7048E+00  4.9001E+00 -4.8757E+00
            -9.0928E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1586.38404623228        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  1.0721E+00  8.0129E-01  1.3722E+00  1.1926E+00  1.0524E+00  1.3023E+00  2.3405E+00  1.4440E+00  8.3562E-01  9.3741E-01
             1.1241E+00
 PARAMETER:  1.6963E-01 -1.2153E-01  4.1643E-01  2.7611E-01  1.5108E-01  3.6414E-01  9.5036E-01  4.6744E-01 -7.9580E-02  3.5364E-02
             2.1702E-01
 GRADIENT:   8.0746E-01 -9.9814E-01  1.0177E+00 -2.7893E+00 -5.2496E-01 -5.5132E-01  1.8761E-01 -4.4600E-01  1.5747E-01 -2.1262E-01
            -9.1066E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1586.42305951031        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      903             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0780E+00  8.4234E-01  1.3219E+00  1.1666E+00  1.0499E+00  1.3194E+00  2.2717E+00  1.4169E+00  8.3930E-01  9.3097E-01
             1.1263E+00
 PARAMETER:  1.7509E-01 -7.1572E-02  3.7908E-01  2.5408E-01  1.4871E-01  3.7714E-01  9.2054E-01  4.4848E-01 -7.5189E-02  2.8471E-02
             2.1890E-01
 GRADIENT:   6.3698E+02  2.1544E+01  4.6009E+00  2.0964E+02  7.0854E+00  2.2801E+02  1.1220E+02  1.2022E+00  5.6788E+00  3.9688E-01
             1.7628E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1586.43125168703        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  1.0753E+00  8.4559E-01  1.3207E+00  1.1670E+00  1.0501E+00  1.3133E+00  2.2715E+00  1.4176E+00  8.3894E-01  9.3133E-01
             1.1263E+00
 PARAMETER:  1.7261E-01 -6.7721E-02  3.7820E-01  2.5441E-01  1.4889E-01  3.7258E-01  9.2044E-01  4.4897E-01 -7.5612E-02  2.8853E-02
             2.1894E-01
 GRADIENT:   4.1024E+00  5.4131E-01 -3.7481E-02  2.1761E-01 -6.7126E-02  2.8550E+00  6.9650E-01  4.3264E-02  1.5614E-02 -1.6530E-02
             1.2715E-02

0ITERATION NO.:   33    OBJECTIVE VALUE:  -1586.43202343728        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1165
 NPARAMETR:  1.0753E+00  8.4418E-01  1.3208E+00  1.1664E+00  1.0501E+00  1.3134E+00  2.2703E+00  1.4173E+00  8.3891E-01  9.3141E-01
             1.1263E+00
 PARAMETER:  1.7262E-01 -6.9384E-02  3.7823E-01  2.5388E-01  1.4891E-01  3.7259E-01  9.1992E-01  4.4872E-01 -7.5650E-02  2.8940E-02
             2.1892E-01
 GRADIENT:  -2.0738E-04 -1.1291E-01  3.5084E-02  6.9821E-02  7.2823E-04  1.0150E-05  2.0836E-02  2.5150E-02  2.3049E-02 -1.0338E-02
             1.7339E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1165
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.6075E-04  1.9033E-02 -4.5018E-02 -2.8369E-02 -3.2075E-02
 SE:             2.9901E-02  2.3588E-02  1.6594E-02  1.9987E-02  1.8676E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7970E-01  4.1973E-01  6.6691E-03  1.5579E-01  8.5896E-02

 ETASHRINKSD(%)  1.0000E-10  2.0976E+01  4.4408E+01  3.3042E+01  3.7434E+01
 ETASHRINKVR(%)  1.0000E-10  3.7553E+01  6.9096E+01  5.5167E+01  6.0855E+01
 EBVSHRINKSD(%)  3.2303E-01  2.0294E+01  4.8903E+01  3.3974E+01  3.4073E+01
 EBVSHRINKVR(%)  6.4501E-01  3.6470E+01  7.3891E+01  5.6406E+01  5.6536E+01
 RELATIVEINF(%)  9.9095E+01  1.0035E+01  5.3424E+00  6.3641E+00  1.0487E+01
 EPSSHRINKSD(%)  4.4654E+01
 EPSSHRINKVR(%)  6.9368E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1586.4320234372829     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -851.28119687354467     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1586.432       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.08E+00  8.44E-01  1.32E+00  1.17E+00  1.05E+00  1.31E+00  2.27E+00  1.42E+00  8.39E-01  9.31E-01  1.13E+00
 


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
 #CPUT: Total CPU Time in Seconds,       41.828
Stop Time:
Sun Oct 24 04:41:24 CDT 2021
