Sat Sep 18 13:27:07 CDT 2021
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
$DATA ../../../../data/spa/S2/dat43.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1685.45065124899        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.3231E+01 -4.2593E+01 -1.6377E+01 -2.6234E+01  3.1212E+01  3.9439E+00  8.3118E+00  7.8984E+00  3.9239E+01 -5.1465E+00
             3.2097E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1692.19730855429        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9427E-01  1.0270E+00  9.9099E-01  9.9634E-01  9.8547E-01  9.8599E-01  9.5615E-01  9.5210E-01  8.1085E-01  1.0373E+00
             9.1722E-01
 PARAMETER:  9.4254E-02  1.2667E-01  9.0954E-02  9.6336E-02  8.5361E-02  8.5892E-02  5.5159E-02  5.0918E-02 -1.0967E-01  1.3664E-01
             1.3590E-02
 GRADIENT:   5.8443E+01 -1.1367E+01  4.2985E+00 -2.6033E+01 -4.7134E-01 -8.0907E-01 -2.4152E+00  4.6531E+00 -5.6795E-01 -7.2551E-01
            -4.0511E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1695.26920110788        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.8825E-01  7.9087E-01  8.2725E-01  1.1699E+00  8.0755E-01  9.5495E-01  1.2829E+00  3.6758E-01  6.8803E-01  9.4851E-01
             9.2253E-01
 PARAMETER:  8.8183E-02 -1.3463E-01 -8.9643E-02  2.5695E-01 -1.1375E-01  5.3904E-02  3.4914E-01 -9.0082E-01 -2.7393E-01  4.7138E-02
             1.9364E-02
 GRADIENT:   4.0275E+01  2.8871E+01 -2.5719E+01  1.0482E+02  3.3517E+01 -1.4816E+01  3.8241E-01  1.0686E+00 -4.3156E+00  2.9677E+00
            -1.6992E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1696.83827630597        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.7388E-01  8.1410E-01  6.9615E-01  1.1100E+00  7.2590E-01  9.8839E-01  1.2254E+00  2.3165E-01  7.0543E-01  8.2491E-01
             9.2803E-01
 PARAMETER:  7.3532E-02 -1.0567E-01 -2.6219E-01  2.0433E-01 -2.2034E-01  8.8322E-02  3.0326E-01 -1.3625E+00 -2.4895E-01 -9.2485E-02
             2.5311E-02
 GRADIENT:   2.9393E+00  6.0904E+00 -3.4431E+00  1.6256E+01  6.6381E+00 -1.5102E+00 -3.3486E-01  5.1131E-01  3.1772E-02  7.8492E-01
             1.9934E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1696.88006737580        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.7347E-01  7.8190E-01  6.3017E-01  1.1148E+00  6.6941E-01  9.9217E-01  1.2624E+00  1.7653E-01  6.9236E-01  7.5836E-01
             9.2625E-01
 PARAMETER:  7.3116E-02 -1.4603E-01 -3.6176E-01  2.0870E-01 -3.0135E-01  9.2136E-02  3.3305E-01 -1.6343E+00 -2.6764E-01 -1.7660E-01
             2.3394E-02
 GRADIENT:   3.3826E-01  5.3778E+00 -2.7815E+00  1.1157E+01  5.1163E+00 -3.0280E-01  8.1044E-02  3.4638E-01  3.5262E-01  7.6705E-01
             1.2773E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1696.88189685888        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      401
 NPARAMETR:  9.7364E-01  7.6818E-01  6.1378E-01  1.1191E+00  6.5317E-01  9.9251E-01  1.2777E+00  1.6461E-01  6.8734E-01  7.4048E-01
             9.2611E-01
 PARAMETER:  7.3290E-02 -1.6373E-01 -3.8813E-01  2.1250E-01 -3.2592E-01  9.2485E-02  3.4509E-01 -1.7042E+00 -2.7492E-01 -2.0045E-01
             2.3238E-02
 GRADIENT:  -4.3833E+01  2.0234E+00 -4.7093E+00 -8.3975E+00 -2.4577E-01 -5.3237E+00 -1.1690E+00  2.8689E-01 -7.5673E-01  5.3808E-01
             1.0305E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1699.69715002927        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      580
 NPARAMETR:  9.9095E-01  5.2660E-01  7.8375E-01  1.2931E+00  6.7685E-01  9.8875E-01  1.7117E+00  4.6862E-02  6.4233E-01  8.9075E-01
             9.2612E-01
 PARAMETER:  9.0909E-02 -5.4131E-01 -1.4366E-01  3.5704E-01 -2.9031E-01  8.8686E-02  6.3749E-01 -2.9605E+00 -3.4265E-01 -1.5686E-02
             2.3247E-02
 GRADIENT:   6.2674E+00  1.2240E+01  5.2202E+00  2.7594E+01 -1.0196E+01 -4.5132E+00 -1.0609E+00  7.8926E-03 -2.1046E+00  2.5905E+00
            -1.2800E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1701.67831178155        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      755
 NPARAMETR:  9.8414E-01  2.6251E-01  7.7356E-01  1.4199E+00  6.0638E-01  9.9843E-01  2.7015E+00  1.0000E-02  6.1404E-01  8.6623E-01
             9.2619E-01
 PARAMETER:  8.4016E-02 -1.2375E+00 -1.5675E-01  4.5061E-01 -4.0024E-01  9.8427E-02  1.0938E+00 -7.1236E+00 -3.8770E-01 -4.3610E-02
             2.3322E-02
 GRADIENT:   2.5190E+00  2.8210E+00  5.3158E+00  3.7735E+00 -6.7613E+00  1.3761E+00  1.3683E+00  0.0000E+00  9.8473E-01 -2.2725E+00
            -3.9452E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1701.92345304071        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      930
 NPARAMETR:  9.8105E-01  1.8902E-01  7.8590E-01  1.4571E+00  6.0029E-01  9.9289E-01  3.2335E+00  1.0000E-02  6.0416E-01  8.8256E-01
             9.3547E-01
 PARAMETER:  8.0869E-02 -1.5659E+00 -1.4092E-01  4.7644E-01 -4.1035E-01  9.2869E-02  1.2736E+00 -9.4586E+00 -4.0392E-01 -2.4932E-02
             3.3297E-02
 GRADIENT:   5.1636E-02  4.4310E-02  1.9194E-01 -3.0407E-02 -1.4456E-01  2.2076E-03  3.9930E-02  0.0000E+00 -7.4646E-02 -9.7915E-04
            -2.1478E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1701.92353538129        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1022
 NPARAMETR:  9.8103E-01  1.8839E-01  7.8508E-01  1.4572E+00  5.9974E-01  9.9288E-01  3.2383E+00  1.0000E-02  6.0425E-01  8.8194E-01
             9.3552E-01
 PARAMETER:  8.0844E-02 -1.5692E+00 -1.4197E-01  4.7654E-01 -4.1126E-01  9.2856E-02  1.2751E+00 -9.4832E+00 -4.0376E-01 -2.5637E-02
             3.3349E-02
 GRADIENT:   1.4327E-02  5.8080E-04 -1.0408E-02  5.6000E-02  1.5453E-02  7.5040E-04 -4.5316E-03  0.0000E+00 -2.1451E-03 -6.2675E-03
             9.4082E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1022
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.8859E-04  3.5200E-02 -3.9909E-04 -2.4181E-02  2.5315E-03
 SE:             2.9898E-02  1.6237E-02  2.5231E-04  2.5674E-02  2.5465E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7896E-01  3.0171E-02  1.1370E-01  3.4627E-01  9.2081E-01

 ETASHRINKSD(%)  1.0000E-10  4.5603E+01  9.9155E+01  1.3990E+01  1.4689E+01
 ETASHRINKVR(%)  1.0000E-10  7.0410E+01  9.9993E+01  2.6023E+01  2.7220E+01
 EBVSHRINKSD(%)  3.7082E-01  5.4508E+01  9.9203E+01  1.0367E+01  1.0916E+01
 EBVSHRINKVR(%)  7.4026E-01  7.9305E+01  9.9994E+01  1.9659E+01  2.0641E+01
 RELATIVEINF(%)  9.8495E+01  4.7006E+00  4.9336E-04  2.0046E+01  6.0734E+00
 EPSSHRINKSD(%)  4.3967E+01
 EPSSHRINKVR(%)  6.8603E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1701.9235353812876     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -966.77270881754941     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.37
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1701.924       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.88E-01  7.85E-01  1.46E+00  6.00E-01  9.93E-01  3.24E+00  1.00E-02  6.04E-01  8.82E-01  9.36E-01
 


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
+        1.16E+03
 
 TH 2
+       -3.00E+01  8.73E+02
 
 TH 3
+        1.51E+01  2.11E+02  1.20E+03
 
 TH 4
+       -8.95E+00  4.85E+02 -2.93E+02  1.20E+03
 
 TH 5
+       -2.96E+00 -5.37E+02 -1.81E+03  1.08E+02  3.20E+03
 
 TH 6
+       -5.65E+00 -3.25E+00  9.50E+00 -3.02E+00  1.89E+00  1.93E+02
 
 TH 7
+        7.04E-01  4.91E+01 -6.94E-01 -7.46E+00  4.24E+00  2.70E-01  5.45E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.55E+00 -1.32E+01 -1.20E+01 -2.07E+01  1.84E+01 -1.96E+00  3.94E+00  0.00E+00  3.55E+02
 
 TH10
+       -3.25E+00 -7.80E+00 -8.21E+01 -2.43E+01 -3.32E+01 -1.47E+00 -1.22E-01  0.00E+00  2.08E+00  1.40E+02
 
 TH11
+       -1.08E+01 -1.57E+01 -4.03E+01 -6.96E+00  2.36E+01 -7.94E-02 -1.44E+00  0.00E+00  1.33E+01  4.40E+01  2.41E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.063
Stop Time:
Sat Sep 18 13:27:27 CDT 2021
