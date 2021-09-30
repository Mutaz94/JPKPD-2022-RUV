Wed Sep 29 12:29:00 CDT 2021
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
$DATA ../../../../data/spa/A1/dat96.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1425.22189579382        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8548E+02  1.7370E+01 -2.8021E+01  7.0747E+01  4.3711E+01  2.5195E+01  1.0216E+01  1.2146E+01  1.5171E+01  1.0194E+01
            -4.4039E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1526.72683014471        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.6153E-01  1.0390E+00  1.1757E+00  9.7273E-01  1.0910E+00  1.0147E+00  8.7720E-01  8.8725E-01  8.6720E-01  7.7723E-01
             2.0409E+00
 PARAMETER:  6.0769E-02  1.3823E-01  2.6184E-01  7.2356E-02  1.8711E-01  1.1460E-01 -3.1025E-02 -1.9625E-02 -4.2484E-02 -1.5202E-01
             8.1341E-01
 GRADIENT:   7.6858E+01 -3.0487E+01 -2.0918E+01 -2.4491E+01  4.1840E+01  1.5846E+01  2.3740E+00  4.0458E+00 -4.7523E+00  6.6193E+00
             6.9739E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1529.34402248751        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.5273E-01  1.0077E+00  1.2823E+00  1.0156E+00  1.0261E+00  9.9873E-01  5.2617E-01  5.4009E-01  1.0872E+00  7.0486E-01
             1.9896E+00
 PARAMETER:  5.1577E-02  1.0764E-01  3.4864E-01  1.1544E-01  1.2572E-01  9.8725E-02 -5.4213E-01 -5.1602E-01  1.8361E-01 -2.4976E-01
             7.8795E-01
 GRADIENT:   6.4768E+01  3.3796E+01  2.2573E+01  3.6465E+01 -5.5279E+01  8.0356E+00  2.7712E+00  8.5022E-02  1.9995E+01 -1.8300E+00
             5.0098E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1531.53248499827        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      224
 NPARAMETR:  9.2601E-01  9.8792E-01  1.1057E+00  1.0062E+00  9.9219E-01  9.9329E-01  5.1853E-01  2.7335E-01  1.0478E+00  8.5035E-01
             1.8334E+00
 PARAMETER:  2.3126E-02  8.7850E-02  2.0051E-01  1.0614E-01  9.2162E-02  9.3268E-02 -5.5677E-01 -1.1970E+00  1.4668E-01 -6.2105E-02
             7.0620E-01
 GRADIENT:   1.9741E+01  6.0068E+00  6.7465E-01  1.1172E+01 -1.1710E+01  3.5281E+00  1.7219E+00  2.5614E-01  9.4040E+00  2.7779E+00
             2.9874E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1534.88905316723        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      402
 NPARAMETR:  9.6961E-01  9.8183E-01  1.1926E+00  1.0218E+00  1.0371E+00  1.0098E+00  5.3176E-01  2.7810E-01  1.0150E+00  8.9093E-01
             1.7568E+00
 PARAMETER:  6.9136E-02  8.1665E-02  2.7617E-01  1.2161E-01  1.3645E-01  1.0970E-01 -5.3157E-01 -1.1798E+00  1.1493E-01 -1.5487E-02
             6.6349E-01
 GRADIENT:  -8.0302E-02 -5.9677E-01 -4.2597E-01 -1.6510E-01  1.2329E+00 -1.1142E-01  1.0029E-01  5.2785E-02 -2.6305E-02 -1.8586E-01
            -3.3001E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1534.90546434941        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      579            RESET HESSIAN, TYPE II
 NPARAMETR:  9.6992E-01  9.9534E-01  1.2003E+00  1.0120E+00  1.0451E+00  1.0102E+00  4.6559E-01  1.5906E-01  1.0374E+00  9.0878E-01
             1.7568E+00
 PARAMETER:  6.9458E-02  9.5330E-02  2.8257E-01  1.1190E-01  1.4416E-01  1.1010E-01 -6.6445E-01 -1.7384E+00  1.3673E-01  4.3497E-03
             6.6350E-01
 GRADIENT:   1.3423E+02  1.1671E+01  1.3459E+00  2.2776E+01  2.5937E+00  1.5054E+01  1.5557E+00  1.6065E-02  3.4730E+00  1.4229E-01
             5.2694E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1534.91247566901        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  9.6937E-01  1.0006E+00  1.1787E+00  1.0099E+00  1.0395E+00  1.0098E+00  5.2292E-01  1.5212E-02  1.0254E+00  9.0232E-01
             1.7556E+00
 PARAMETER:  6.8895E-02  1.0064E-01  2.6443E-01  1.0980E-01  1.3874E-01  1.0978E-01 -5.4833E-01 -4.0857E+00  1.2506E-01 -2.7862E-03
             6.6279E-01
 GRADIENT:  -9.6433E-01  7.1635E-01  3.3274E-01  9.0501E-01 -6.9244E-02 -1.1532E-01 -5.8043E-03  8.1369E-05 -4.0160E-01  1.6684E-01
            -2.8960E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1534.91879399022        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      929
 NPARAMETR:  9.6947E-01  1.0101E+00  1.1340E+00  1.0032E+00  1.0254E+00  1.0103E+00  5.9827E-01  1.0000E-02  1.0158E+00  8.7909E-01
             1.7548E+00
 PARAMETER:  6.8992E-02  1.1009E-01  2.2571E-01  1.0318E-01  1.2508E-01  1.1021E-01 -4.1372E-01 -8.3804E+00  1.1569E-01 -2.8870E-02
             6.6233E-01
 GRADIENT:  -1.2174E+00 -4.7232E-01 -2.1400E-01 -1.8297E-01  1.3470E-01 -9.5568E-02 -1.2567E-02  0.0000E+00  1.2930E-02  3.5133E-02
            -5.1254E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1534.91879399022        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      951
 NPARAMETR:  9.6947E-01  1.0101E+00  1.1340E+00  1.0032E+00  1.0254E+00  1.0103E+00  5.9827E-01  1.0000E-02  1.0158E+00  8.7909E-01
             1.7548E+00
 PARAMETER:  6.8992E-02  1.1009E-01  2.2571E-01  1.0318E-01  1.2508E-01  1.1021E-01 -4.1372E-01 -8.3804E+00  1.1569E-01 -2.8870E-02
             6.6233E-01
 GRADIENT:  -1.2174E+00 -4.7232E-01 -2.1400E-01 -1.8297E-01  1.3470E-01 -9.5568E-02 -1.2567E-02  0.0000E+00  1.2930E-02  3.5133E-02
            -5.1254E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      951
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.9584E-04 -1.6723E-02 -3.1958E-05 -2.0173E-03 -2.7495E-02
 SE:             2.9579E-02  1.1740E-02  1.0893E-04  2.5636E-02  2.0840E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8123E-01  1.5432E-01  7.6923E-01  9.3728E-01  1.8706E-01

 ETASHRINKSD(%)  9.0818E-01  6.0669E+01  9.9635E+01  1.4117E+01  3.0184E+01
 ETASHRINKVR(%)  1.8081E+00  8.4530E+01  9.9999E+01  2.6240E+01  5.1257E+01
 EBVSHRINKSD(%)  1.1509E+00  6.1210E+01  9.9620E+01  1.3846E+01  3.0496E+01
 EBVSHRINKVR(%)  2.2885E+00  8.4953E+01  9.9999E+01  2.5775E+01  5.1692E+01
 RELATIVEINF(%)  9.7131E+01  4.8260E-01  2.2203E-04  3.2501E+00  4.7864E+00
 EPSSHRINKSD(%)  3.4956E+01
 EPSSHRINKVR(%)  5.7692E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1534.9187939902217     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -799.76796742648355     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.47
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1534.919       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.69E-01  1.01E+00  1.13E+00  1.00E+00  1.03E+00  1.01E+00  5.98E-01  1.00E-02  1.02E+00  8.79E-01  1.75E+00
 


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
+        1.13E+03
 
 TH 2
+       -2.58E+01  4.68E+02
 
 TH 3
+        9.65E+00  1.07E+02  1.25E+02
 
 TH 4
+       -3.45E+01  5.07E+02 -1.43E+01  7.98E+02
 
 TH 5
+       -3.99E+00 -2.95E+02 -2.54E+02 -3.26E+01  6.57E+02
 
 TH 6
+        1.85E+00 -3.92E+00  3.11E+00 -7.84E+00 -2.49E+00  1.86E+02
 
 TH 7
+        1.19E+00 -1.71E+01  9.24E+00 -1.39E+01 -5.92E+00  4.85E-02  1.23E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.45E+00 -4.02E+01  1.73E+00  2.25E+01  5.90E-03 -1.55E+00  2.66E+01  0.00E+00  1.08E+02
 
 TH10
+       -2.80E+00  8.88E+00  1.67E+00 -8.02E+00 -4.91E+01  1.81E-01  1.01E+01  0.00E+00  4.02E-01  6.27E+01
 
 TH11
+       -9.59E+00 -1.72E+01 -1.21E+01 -1.31E+01 -4.85E+00  3.49E+00  5.09E+00  0.00E+00  1.06E+01  2.76E+01  8.58E+01
 
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
 #CPUT: Total CPU Time in Seconds,       16.288
Stop Time:
Wed Sep 29 12:29:18 CDT 2021
