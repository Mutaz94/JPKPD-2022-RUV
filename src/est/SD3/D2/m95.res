Sun Oct 24 01:24:48 CDT 2021
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
$DATA ../../../../data/SD3/D2/dat95.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m95.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1937.28416997984        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2118E+02 -8.4933E+00 -3.4917E+01 -2.0621E+01  6.5942E+01 -5.5530E+01 -6.5433E+01  3.9563E+00 -1.4809E+02 -2.1968E+00
             1.9821E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1986.86337264824        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      130
 NPARAMETR:  1.0067E+00  9.2996E-01  1.0603E+00  1.0540E+00  9.5068E-01  1.1359E+00  1.1878E+00  9.8728E-01  1.4552E+00  9.8879E-01
             9.5578E-01
 PARAMETER:  1.0671E-01  2.7386E-02  1.5856E-01  1.5259E-01  4.9421E-02  2.2739E-01  2.7207E-01  8.7201E-02  4.7513E-01  8.8732E-02
             5.4770E-02
 GRADIENT:   3.2002E+01 -1.1176E+01 -1.3395E+01 -3.2735E+01 -4.2675E+00 -4.3942E+01 -1.3429E+01  4.0651E-01 -3.1009E+00  5.1408E-01
            -8.7519E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1994.89515743688        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  1.0088E+00  6.3753E-01  1.3975E+00  1.3264E+00  9.7124E-01  1.2153E+00  2.2897E+00  1.0994E+00  1.1243E+00  1.0760E+00
             9.7569E-01
 PARAMETER:  1.0878E-01 -3.5015E-01  4.3471E-01  3.8243E-01  7.0820E-02  2.9502E-01  9.2841E-01  1.9476E-01  2.1716E-01  1.7324E-01
             7.5392E-02
 GRADIENT:   3.7817E+01  1.4208E+01  2.0958E+01  2.7917E+01 -1.8789E+01 -1.2378E+01 -7.7950E+00 -1.1485E-01 -1.1215E+01  1.0563E+01
             8.4577E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1999.88248763704        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      488
 NPARAMETR:  9.8303E-01  6.6330E-01  1.0993E+00  1.2536E+00  8.7591E-01  1.2552E+00  2.3647E+00  8.1329E-01  1.1600E+00  8.7439E-01
             9.6428E-01
 PARAMETER:  8.2885E-02 -3.1053E-01  1.9469E-01  3.2605E-01 -3.2495E-02  3.2731E-01  9.6064E-01 -1.0666E-01  2.4840E-01 -3.4226E-02
             6.3630E-02
 GRADIENT:  -2.5915E+00  5.0320E+00  7.8291E+00  2.9648E+00 -8.5525E+00  1.2226E+00 -1.9214E+00  3.8554E-02  1.6731E+00 -5.2808E-01
            -3.5089E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2000.29349889877        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      663
 NPARAMETR:  9.8238E-01  5.0958E-01  1.1126E+00  1.3356E+00  8.4154E-01  1.2492E+00  2.7830E+00  8.1348E-01  1.0995E+00  8.7834E-01
             9.6429E-01
 PARAMETER:  8.2224E-02 -5.7418E-01  2.0667E-01  3.8938E-01 -7.2526E-02  3.2253E-01  1.1235E+00 -1.0644E-01  1.9483E-01 -2.9724E-02
             6.3633E-02
 GRADIENT:  -4.3370E-01 -1.8452E-01 -2.6208E-01 -2.6984E-01  1.9038E-01 -2.5327E-01 -3.8512E-01  5.7209E-02 -4.1800E-01  6.8573E-02
             2.6283E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2000.31796583729        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      827             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8386E-01  5.1108E-01  1.1117E+00  1.3315E+00  8.4183E-01  1.2562E+00  2.8119E+00  8.1272E-01  1.1016E+00  8.7711E-01
             9.6389E-01
 PARAMETER:  8.3727E-02 -5.7123E-01  2.0592E-01  3.8627E-01 -7.2180E-02  3.2811E-01  1.1339E+00 -1.0737E-01  1.9675E-01 -3.1125E-02
             6.3227E-02
 GRADIENT:   4.2356E+02  7.9179E+01  6.4983E+00  5.3163E+02  8.9109E+00  2.0470E+02  1.3479E+02  2.8254E-01  2.9793E+01  7.5915E-01
             1.0819E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2000.31831960860        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  9.8370E-01  5.1138E-01  1.1114E+00  1.3310E+00  8.4183E-01  1.2562E+00  2.8089E+00  8.1208E-01  1.1017E+00  8.7696E-01
             9.6386E-01
 PARAMETER:  8.3567E-02 -5.7064E-01  2.0563E-01  3.8592E-01 -7.2174E-02  3.2813E-01  1.1328E+00 -1.0815E-01  1.9690E-01 -3.1289E-02
             6.3195E-02
 GRADIENT:   1.4671E+00  6.4030E-03 -2.7412E-02 -4.2047E+00  6.8950E-02  2.0177E+00  1.5819E+00  6.3181E-02  7.7149E-02  8.4751E-02
             4.7312E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2000.31872009562        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1196             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8372E-01  5.1234E-01  1.1111E+00  1.3305E+00  8.4180E-01  1.2563E+00  2.8061E+00  8.1000E-01  1.1020E+00  8.7629E-01
             9.6382E-01
 PARAMETER:  8.3583E-02 -5.6876E-01  2.0534E-01  3.8552E-01 -7.2209E-02  3.2814E-01  1.1318E+00 -1.1072E-01  1.9711E-01 -3.2055E-02
             6.3151E-02
 GRADIENT:   4.2340E+02  7.9084E+01  6.8900E+00  5.2981E+02  8.6425E+00  2.0483E+02  1.3435E+02  1.6461E-01  2.9820E+01  6.0839E-01
             9.4784E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2000.31896947211        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1383
 NPARAMETR:  9.8373E-01  5.1282E-01  1.1107E+00  1.3301E+00  8.4184E-01  1.2563E+00  2.8045E+00  8.0997E-01  1.1021E+00  8.7622E-01
             9.6384E-01
 PARAMETER:  8.3592E-02 -5.6783E-01  2.0500E-01  3.8528E-01 -7.2170E-02  3.2815E-01  1.1312E+00 -1.1075E-01  1.9726E-01 -3.2135E-02
             6.3173E-02
 GRADIENT:   1.4671E+00  6.3579E-02  2.3673E-01 -4.1656E+00 -1.2286E-01  2.0187E+00  1.5753E+00 -2.4462E-02  6.7816E-02 -2.7648E-02
            -2.5422E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -2000.31896947211        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     1407
 NPARAMETR:  9.8373E-01  5.1282E-01  1.1107E+00  1.3301E+00  8.4184E-01  1.2563E+00  2.8045E+00  8.0997E-01  1.1021E+00  8.7622E-01
             9.6384E-01
 PARAMETER:  8.3592E-02 -5.6783E-01  2.0500E-01  3.8528E-01 -7.2170E-02  3.2815E-01  1.1312E+00 -1.1075E-01  1.9726E-01 -3.2135E-02
             6.3173E-02
 GRADIENT:  -8.2163E-03 -6.4145E-02  2.4627E-01  2.1516E-01 -1.1110E-01 -1.6537E-03  6.1738E-02 -1.4343E-02 -2.9014E-02 -2.8425E-02
            -2.5073E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1407
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0820E-03  3.0165E-02 -3.6146E-02 -2.2961E-02 -1.0711E-02
 SE:             2.9957E-02  2.0128E-02  1.4156E-02  2.5287E-02  2.0265E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7119E-01  1.3396E-01  1.0668E-02  3.6388E-01  5.9714E-01

 ETASHRINKSD(%)  1.0000E-10  3.2569E+01  5.2575E+01  1.5284E+01  3.2108E+01
 ETASHRINKVR(%)  1.0000E-10  5.4531E+01  7.7508E+01  2.8232E+01  5.3907E+01
 EBVSHRINKSD(%)  2.1467E-01  3.5622E+01  5.4925E+01  1.3424E+01  2.8706E+01
 EBVSHRINKVR(%)  4.2888E-01  5.8554E+01  7.9683E+01  2.5046E+01  4.9172E+01
 RELATIVEINF(%)  9.9211E+01  9.6725E+00  5.0966E+00  2.0969E+01  1.1357E+01
 EPSSHRINKSD(%)  3.4461E+01
 EPSSHRINKVR(%)  5.7047E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2000.3189694721141     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1081.3804362674414     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2000.319       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  5.13E-01  1.11E+00  1.33E+00  8.42E-01  1.26E+00  2.80E+00  8.10E-01  1.10E+00  8.76E-01  9.64E-01
 


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
 #CPUT: Total CPU Time in Seconds,       51.768
Stop Time:
Sun Oct 24 01:24:58 CDT 2021
