Sun Oct 24 04:41:36 CDT 2021
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
$DATA ../../../../data/SD4/D2/dat72.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m72.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1659.79885282261        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2459E+02 -1.5126E+01  7.0334E+00  9.1162E+00 -1.6772E+01  6.6999E+00 -3.6626E+01 -5.2436E+00  6.6543E+00 -1.3321E+01
             1.7455E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1675.84837709795        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.1523E+00  1.2264E+00  1.0226E+00  9.6117E-01  1.1632E+00  1.3505E+00  1.5583E+00  1.0592E+00  9.3168E-01  1.1386E+00
             9.1619E-01
 PARAMETER:  2.4176E-01  3.0412E-01  1.2238E-01  6.0396E-02  2.5120E-01  4.0048E-01  5.4360E-01  1.5748E-01  2.9229E-02  2.2976E-01
             1.2473E-02
 GRADIENT:   1.0499E+02  3.6929E+01 -1.0711E+01  6.7356E+01  1.8933E+01  5.4530E+01  1.1730E+01 -1.3417E+00  2.0605E+01  1.1452E+00
            -1.3448E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1678.47200106120        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.1361E+00  9.4575E-01  1.4378E+00  1.1415E+00  1.2145E+00  1.3011E+00  1.8267E+00  1.4214E+00  7.1577E-01  1.1919E+00
             8.9060E-01
 PARAMETER:  2.2759E-01  4.4223E-02  4.6313E-01  2.3235E-01  2.9432E-01  3.6319E-01  7.0251E-01  4.5167E-01 -2.3440E-01  2.7551E-01
            -1.5860E-02
 GRADIENT:   9.7689E+01  3.4443E+01 -3.6064E+00  8.9461E+01  1.8620E+01  4.5629E+01 -1.1907E+00 -2.5797E+00 -2.7966E+00 -2.2768E+00
            -2.9738E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1690.81296092197        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  1.0648E+00  1.2771E+00  8.5096E-01  8.7047E-01  1.0948E+00  1.1286E+00  1.4575E+00  9.8924E-01  6.7594E-01  1.0572E+00
             9.4250E-01
 PARAMETER:  1.6282E-01  3.4456E-01 -6.1386E-02 -3.8726E-02  1.9060E-01  2.2098E-01  4.7673E-01  8.9179E-02 -2.9165E-01  1.5562E-01
             4.0781E-02
 GRADIENT:  -2.3614E+00  1.1843E+01  1.9267E+00  1.4033E+01 -5.5963E+00  6.9211E-01 -9.5649E-01  5.9333E-01 -2.2846E+00 -1.2912E-02
            -9.0941E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1691.33846935714        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  1.0693E+00  1.4746E+00  6.0233E-01  7.3121E-01  1.0631E+00  1.1262E+00  1.2734E+00  6.2511E-01  7.6134E-01  1.0116E+00
             9.4384E-01
 PARAMETER:  1.6700E-01  4.8838E-01 -4.0696E-01 -2.1305E-01  1.6118E-01  2.1882E-01  3.4172E-01 -3.6983E-01 -1.7268E-01  1.1156E-01
             4.2206E-02
 GRADIENT:   1.2979E+00  2.5197E+00  1.3849E-01  3.0481E+00 -1.7031E+00 -4.6497E-01 -4.7891E-01  1.0967E-01 -6.6967E-03 -2.6055E-01
            -1.1957E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1691.35954759954        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  1.0723E+00  1.4877E+00  5.9051E-01  7.1990E-01  1.0679E+00  1.1305E+00  1.2652E+00  5.9260E-01  7.6745E-01  1.0165E+00
             9.4420E-01
 PARAMETER:  1.6981E-01  4.9722E-01 -4.2678E-01 -2.2864E-01  1.6570E-01  2.2264E-01  3.3526E-01 -4.2324E-01 -1.6468E-01  1.1641E-01
             4.2581E-02
 GRADIENT:   6.2031E+00 -1.1502E+00 -3.0845E-01  7.3718E-01  2.8866E-01  1.0788E+00  2.5007E-01  1.3018E-02  8.0421E-02 -1.5206E-01
             2.3543E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1691.36055668594        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1085             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0721E+00  1.4865E+00  5.9081E-01  7.1956E-01  1.0680E+00  1.1299E+00  1.2650E+00  5.9088E-01  7.6640E-01  1.0178E+00
             9.4411E-01
 PARAMETER:  1.6962E-01  4.9640E-01 -4.2627E-01 -2.2912E-01  1.6576E-01  2.2212E-01  3.3511E-01 -4.2614E-01 -1.6605E-01  1.1765E-01
             4.2484E-02
 GRADIENT:   8.9154E+02  4.1696E+02  6.1520E+00  1.1932E+02  1.3163E+01  1.4576E+02  4.5232E+01  1.7543E-01  4.8501E+00  1.1340E+00
             7.7611E-01

0ITERATION NO.:   32    OBJECTIVE VALUE:  -1691.36055668594        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     1144
 NPARAMETR:  1.0721E+00  1.4865E+00  5.9081E-01  7.1956E-01  1.0680E+00  1.1299E+00  1.2650E+00  5.9088E-01  7.6640E-01  1.0178E+00
             9.4411E-01
 PARAMETER:  1.6962E-01  4.9640E-01 -4.2627E-01 -2.2912E-01  1.6576E-01  2.2212E-01  3.3511E-01 -4.2614E-01 -1.6605E-01  1.1765E-01
             4.2484E-02
 GRADIENT:   7.4791E-03  3.3863E-02  7.3849E-02 -2.4276E-01 -2.0576E-01  4.0502E-04 -4.3788E-02 -7.1736E-03 -2.4469E-02 -3.6325E-02
            -2.3464E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1144
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.7318E-04 -5.3496E-03 -2.1290E-02  3.1760E-03 -2.4145E-02
 SE:             2.9888E-02  2.7034E-02  7.8434E-03  1.7943E-02  2.2181E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9538E-01  8.4314E-01  6.6400E-03  8.5950E-01  2.7636E-01

 ETASHRINKSD(%)  1.0000E-10  9.4335E+00  7.3724E+01  3.9890E+01  2.5691E+01
 ETASHRINKVR(%)  1.0000E-10  1.7977E+01  9.3096E+01  6.3867E+01  4.4782E+01
 EBVSHRINKSD(%)  2.9425E-01  9.2666E+00  7.7185E+01  4.2260E+01  2.3034E+01
 EBVSHRINKVR(%)  5.8764E-01  1.7675E+01  9.4795E+01  6.6661E+01  4.0763E+01
 RELATIVEINF(%)  9.9266E+01  8.0706E+00  4.4303E-01  1.8868E+00  1.0188E+01
 EPSSHRINKSD(%)  4.4646E+01
 EPSSHRINKVR(%)  6.9359E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1691.3605566859417     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -956.20973012220350     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1691.361       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.07E+00  1.49E+00  5.91E-01  7.20E-01  1.07E+00  1.13E+00  1.27E+00  5.91E-01  7.66E-01  1.02E+00  9.44E-01
 


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
 #CPUT: Total CPU Time in Seconds,       39.778
Stop Time:
Sun Oct 24 04:41:45 CDT 2021
