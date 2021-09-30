Wed Sep 29 13:50:34 CDT 2021
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
$DATA ../../../../data/spa/A3/dat82.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m82.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -333.778993844099        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3489E+02  4.4877E+01  3.7477E+01  2.9685E+01  2.4602E+02  7.7300E+01 -7.3378E+01 -1.7421E+01 -1.1458E+02 -1.7849E+02
            -2.2410E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1273.41582855881        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8434E-01  9.6450E-01  9.9648E-01  1.0747E+00  9.2451E-01  7.0252E-01  1.0362E+00  9.6285E-01  1.1410E+00  9.6635E-01
             3.0628E+00
 PARAMETER:  8.4218E-02  6.3858E-02  9.6475E-02  1.7206E-01  2.1507E-02 -2.5308E-01  1.3555E-01  6.2146E-02  2.3189E-01  6.5773E-02
             1.2193E+00
 GRADIENT:  -3.2178E+01 -7.3386E+00 -2.3414E+01  1.5305E+01  4.9580E+01 -3.8275E+01 -2.6815E-01  4.1498E+00 -2.3259E+00  2.1461E+00
            -8.3507E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1289.90062758481        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.8679E-01  1.0163E+00  9.7694E-01  1.0338E+00  8.8931E-01  7.3680E-01  9.9277E-01  4.7469E-01  1.0787E+00  2.9904E-01
             3.6902E+00
 PARAMETER:  8.6699E-02  1.1622E-01  7.6670E-02  1.3325E-01 -1.7311E-02 -2.0544E-01  9.2741E-02 -6.4508E-01  1.7572E-01 -1.1072E+00
             1.4057E+00
 GRADIENT:  -4.3292E+01  3.6744E+00  2.0422E+00 -2.0492E+00 -4.0176E+00 -1.1743E+01  1.6766E+00  7.0716E-01  2.4006E+00  9.6728E-01
             9.5460E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1292.16068751071        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0007E+00  7.8412E-01  9.3001E-01  1.1801E+00  7.8080E-01  7.7953E-01  1.1444E+00  2.1511E-01  9.8706E-01  2.6115E-01
             3.6538E+00
 PARAMETER:  1.0072E-01 -1.4320E-01  2.7445E-02  2.6560E-01 -1.4744E-01 -1.4907E-01  2.3490E-01 -1.4366E+00  8.6973E-02 -1.2427E+00
             1.3958E+00
 GRADIENT:   9.0325E+00  7.0794E+00  6.6668E+00  9.3215E+00 -1.0304E+01  3.6517E+00  4.8139E-01  1.6424E-01 -4.4022E-01  4.4462E-01
             2.0894E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1292.82367782042        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      364
 NPARAMETR:  1.0104E+00  6.5027E-01  9.4731E-01  1.2907E+00  7.5320E-01  7.7493E-01  1.1348E+00  1.2588E-01  9.6046E-01  2.1876E-01
             3.7636E+00
 PARAMETER:  1.1036E-01 -3.3037E-01  4.5866E-02  3.5521E-01 -1.8342E-01 -1.5498E-01  2.2650E-01 -1.9724E+00  5.9656E-02 -1.4198E+00
             1.4254E+00
 GRADIENT:   2.8127E+00  1.1280E+01  2.9404E+00  2.2623E+01 -9.0093E+00 -8.5661E-01 -4.1922E-01  8.0676E-02  2.5038E-01  2.8823E-01
             7.0362E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1294.03728433913        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      546
 NPARAMETR:  1.0094E+00  3.4971E-01  8.2939E-01  1.4428E+00  6.1393E-01  7.8507E-01  1.7499E+00  1.0000E-02  8.7283E-01  1.1141E-01
             3.7567E+00
 PARAMETER:  1.0934E-01 -9.5065E-01 -8.7061E-02  4.6662E-01 -3.8788E-01 -1.4198E-01  6.5957E-01 -4.9154E+00 -3.6011E-02 -2.0946E+00
             1.4235E+00
 GRADIENT:   1.0656E+01  6.7726E+00  8.6562E+00  2.4632E+01 -1.7925E+01  2.3262E+00 -6.9988E-01  0.0000E+00 -1.6351E+00 -4.1151E-02
             5.2532E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1295.18722615440        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  9.9788E-01  1.4899E-01  8.3924E-01  1.5330E+00  5.8535E-01  7.7428E-01  3.2339E+00  1.0000E-02  8.2722E-01  4.1842E-02
             3.7237E+00
 PARAMETER:  9.7875E-02 -1.8039E+00 -7.5258E-02  5.2720E-01 -4.3554E-01 -1.5582E-01  1.2737E+00 -1.0232E+01 -8.9682E-02 -3.0738E+00
             1.4147E+00
 GRADIENT:  -2.4969E+00  7.9979E-01 -6.7891E-01  1.9133E+00  1.8466E+00 -8.9923E-01  3.4708E-01  0.0000E+00 -7.8978E-01 -6.9302E-03
             4.3281E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1295.66820586960        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      896
 NPARAMETR:  9.9643E-01  3.9855E-02  6.2876E-01  1.5347E+00  4.6164E-01  7.8175E-01  6.1409E+00  1.0000E-02  8.4381E-01  1.0000E-02
             3.7028E+00
 PARAMETER:  9.6425E-02 -3.1225E+00 -3.6400E-01  5.2835E-01 -6.7297E-01 -1.4622E-01  1.9150E+00 -1.9531E+01 -6.9833E-02 -4.9928E+00
             1.4091E+00
 GRADIENT:   2.9990E-01  2.9958E-01 -1.0078E-01  8.6142E+00 -2.5395E+00 -5.8616E-01 -7.7095E-02  0.0000E+00 -1.2651E+00  0.0000E+00
            -1.0542E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1295.81132761902        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1072
 NPARAMETR:  9.9471E-01  1.0016E-02  6.1897E-01  1.5370E+00  4.5304E-01  7.8290E-01  1.2706E+01  1.0000E-02  8.4682E-01  1.0000E-02
             3.7002E+00
 PARAMETER:  9.4692E-02 -4.5036E+00 -3.7969E-01  5.2986E-01 -6.9178E-01 -1.4475E-01  2.6421E+00 -2.9590E+01 -6.6269E-02 -6.8565E+00
             1.4084E+00
 GRADIENT:   9.6138E-01  2.1708E-01 -1.3497E+00 -3.8042E+00  2.1685E+00 -9.9231E-02  4.8750E-01  0.0000E+00  2.6427E-01  0.0000E+00
            -4.5607E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1295.81681170154        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1242
 NPARAMETR:  9.9449E-01  1.0000E-02  6.4181E-01  1.5472E+00  4.6489E-01  7.8201E-01  1.2440E+01  1.0000E-02  8.4150E-01  1.0000E-02
             3.7496E+00
 PARAMETER:  9.4684E-02 -4.5031E+00 -3.4239E-01  5.3649E-01 -6.6678E-01 -1.4554E-01  2.6474E+00 -2.9550E+01 -7.2889E-02 -6.8229E+00
             1.4096E+00
 GRADIENT:   3.3015E-01  6.4063E+00  4.4819E-01  2.8310E-02 -7.1205E-01  5.3287E-02  1.0947E+01  0.0000E+00 -3.2212E-02  0.0000E+00
            -3.8500E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1242
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2015E-03  5.2156E-04  1.1644E-04 -1.4573E-02 -1.6297E-05
 SE:             2.7862E-02  1.9840E-03  1.9138E-04  2.5363E-02  2.9078E-04
 N:                     100         100         100         100         100

 P VAL.:         9.6560E-01  7.9264E-01  5.4292E-01  5.6558E-01  9.5530E-01

 ETASHRINKSD(%)  6.6587E+00  9.3353E+01  9.9359E+01  1.5031E+01  9.9026E+01
 ETASHRINKVR(%)  1.2874E+01  9.9558E+01  9.9996E+01  2.7803E+01  9.9991E+01
 EBVSHRINKSD(%)  6.4217E+00  9.3801E+01  9.9353E+01  1.4760E+01  9.9092E+01
 EBVSHRINKVR(%)  1.2431E+01  9.9616E+01  9.9996E+01  2.7341E+01  9.9992E+01
 RELATIVEINF(%)  7.3275E+01  2.0842E-02  1.3857E-04  6.2135E+00  2.5673E-04
 EPSSHRINKSD(%)  2.0627E+01
 EPSSHRINKVR(%)  3.6999E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1295.8168117015434     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -560.66598513780525     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.23
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1295.817       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  1.00E-02  6.42E-01  1.55E+00  4.65E-01  7.82E-01  1.28E+01  1.00E-02  8.41E-01  1.00E-02  3.70E+00
 


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
+        1.61E+03
 
 TH 2
+       -9.41E+02  6.95E+05
 
 TH 3
+       -2.36E+01  2.61E+03  1.01E+03
 
 TH 4
+       -1.20E+02  3.53E+03 -3.30E+01  2.42E+03
 
 TH 5
+        1.44E+02 -1.29E+03 -1.93E+03  5.22E+03  3.99E+03
 
 TH 6
+       -1.12E+01  1.15E+03  2.08E+01 -1.65E+01 -6.79E+00  2.47E+02
 
 TH 7
+       -1.13E+00 -8.91E+02  3.34E+00  4.42E+00 -1.30E+00  1.51E+00  9.22E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.35E+01 -1.67E+03 -5.51E-01 -1.41E+01  6.29E+01 -3.37E-01 -2.22E+00  0.00E+00  1.43E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.84E+01  6.69E+03 -2.69E+01  2.90E+02  2.29E+01 -2.58E+00  1.22E+00  0.00E+00  2.70E+01  0.00E+00  7.70E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.848
Stop Time:
Wed Sep 29 13:50:57 CDT 2021
