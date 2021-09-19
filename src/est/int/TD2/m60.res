Sat Sep 18 06:08:43 CDT 2021
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
$DATA ../../../../data/int/TD2/dat60.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3317.22958693291        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6280E+01 -5.2030E+01  3.6537E+01 -1.0148E+02  3.6860E+01  5.8750E-01 -2.8952E+00 -2.5195E+02 -5.9196E+01  1.0692E+01
            -8.4613E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3681.05777621393        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  9.2941E-01  1.0234E+00  1.3475E+00  1.1067E+00  1.0883E+00  1.0739E+00  9.5712E-01  9.9295E-01  1.0284E+00  9.5108E-01
             1.6809E+00
 PARAMETER:  2.6796E-02  1.2317E-01  3.9822E-01  2.0138E-01  1.8465E-01  1.7128E-01  5.6177E-02  9.2929E-02  1.2803E-01  4.9841E-02
             6.1931E-01
 GRADIENT:  -1.4074E+02 -1.6891E+01  2.0624E+01  1.1277E+02  2.0824E+01  1.2605E+01 -1.0084E+01 -2.1366E+00  2.0882E+01 -6.7918E-01
             7.0157E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3786.02110651840        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0032E+00  1.0400E+00  9.8887E-01  1.0708E+00  1.0578E+00  7.6641E-01  1.1309E+00  2.3526E-01  1.0470E+00  1.7098E+00
             1.2991E+00
 PARAMETER:  1.0322E-01  1.3920E-01  8.8808E-02  1.6841E-01  1.5616E-01 -1.6604E-01  2.2304E-01 -1.3471E+00  1.4596E-01  6.3637E-01
             3.6164E-01
 GRADIENT:   1.4730E+01  1.6493E+01 -4.1192E+01  9.6558E+01  5.7717E+01 -1.4271E+02  3.4245E+01 -2.6884E-01  1.1573E+01  9.8185E+01
             3.9983E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3876.13401680520        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  9.8431E-01  1.0549E+00  1.1714E+00  1.0330E+00  1.0980E+00  9.0109E-01  9.6578E-01  6.1021E-01  9.6413E-01  1.3247E+00
             1.0577E+00
 PARAMETER:  8.4181E-02  1.5345E-01  2.5816E-01  1.3250E-01  1.9345E-01 -4.1541E-03  6.5180E-02 -3.9395E-01  6.3469E-02  3.8117E-01
             1.5609E-01
 GRADIENT:  -1.4078E+01 -8.5376E+00  4.7751E+01  5.7941E+01  4.8930E+01 -5.0205E+01  3.3564E+00 -1.4606E+01 -1.5007E+01  5.6894E+01
             4.6432E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3892.89267706125        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:      385
 NPARAMETR:  9.8432E-01  1.0236E+00  1.0485E+00  1.0330E+00  9.5964E-01  1.0156E+00  1.0438E+00  6.1142E-01  9.6416E-01  9.6215E-01
             1.0576E+00
 PARAMETER:  8.4196E-02  1.2334E-01  1.4736E-01  1.3243E-01  5.8803E-02  1.1548E-01  1.4287E-01 -3.9197E-01  6.3502E-02  6.1413E-02
             1.5605E-01
 GRADIENT:  -2.0154E+00  3.9918E+01  4.0062E+01  2.7786E+00 -4.8087E+01  5.4446E+00  2.9924E+00 -1.2361E+01 -9.1345E+00 -2.6343E-01
             6.4474E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3894.22408864761        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      522
 NPARAMETR:  9.8432E-01  9.9981E-01  1.0064E+00  1.0330E+00  9.6331E-01  1.0282E+00  1.0244E+00  6.1142E-01  9.6416E-01  9.6216E-01
             1.0576E+00
 PARAMETER:  8.4196E-02  9.9806E-02  1.0639E-01  1.3243E-01  6.2616E-02  1.2778E-01  1.2407E-01 -3.9197E-01  6.3502E-02  6.1423E-02
             1.5605E-01
 GRADIENT:  -4.2445E+01  1.1411E+00 -5.3395E-01 -1.8101E+01 -8.8282E-01  1.8141E-02 -7.6196E-02 -9.7753E+00 -1.1196E+01  2.3144E-01
             6.8391E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3894.87930070104        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      701            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8556E-01  9.9929E-01  1.0075E+00  1.0335E+00  9.6404E-01  1.0281E+00  1.0259E+00  6.3895E-01  9.6624E-01  9.6067E-01
             1.0554E+00
 PARAMETER:  8.5456E-02  9.9294E-02  1.0751E-01  1.3299E-01  6.3378E-02  1.2774E-01  1.2557E-01 -3.4792E-01  6.5658E-02  5.9871E-02
             1.5390E-01
 GRADIENT:   1.8019E+00  6.9565E+00 -1.4096E+00 -3.7577E+00  4.4922E+00  1.1503E+01  8.9481E-01 -9.2742E+00 -9.4308E+00  1.3988E+00
             6.6870E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3895.24752548707        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      797
 NPARAMETR:  9.8571E-01  9.9913E-01  1.0077E+00  1.0337E+00  9.6399E-01  1.0273E+00  1.0258E+00  6.5662E-01  9.6743E-01  9.6014E-01
             1.0539E+00
 PARAMETER:  8.5612E-02  9.9132E-02  1.0772E-01  1.3313E-01  6.3329E-02  1.2689E-01  1.2543E-01 -3.2064E-01  6.6884E-02  5.9326E-02
             1.5246E-01
 GRADIENT:   2.1975E+00  6.6952E+00 -2.4881E+00 -3.4596E+00  4.1500E+00  1.1126E+01  8.6994E-01 -8.9158E+00 -9.0680E+00  1.9125E+00
             6.5580E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3895.26676822392        NO. OF FUNC. EVALS.: 211
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  9.8578E-01  9.9916E-01  1.0113E+00  1.0337E+00  9.6402E-01  1.0273E+00  1.0257E+00  6.5728E-01  9.6747E-01  9.6010E-01
             1.0538E+00
 PARAMETER:  8.5674E-02  9.9157E-02  1.1125E-01  1.3311E-01  6.3358E-02  1.2692E-01  1.2540E-01 -3.1965E-01  6.6924E-02  5.9282E-02
             1.5243E-01
 GRADIENT:  -3.9374E+01 -4.8205E-01 -4.4600E-02 -1.6813E+01 -3.0882E+00 -2.1738E-01 -2.8641E-02 -9.1798E+00 -9.9976E+00  1.5060E+00
             6.5138E+01

0ITERATION NO.:   41    OBJECTIVE VALUE:  -3895.26676822392        NO. OF FUNC. EVALS.:  31
 CUMULATIVE NO. OF FUNC. EVALS.:     1039
 NPARAMETR:  9.8578E-01  9.9916E-01  1.0113E+00  1.0337E+00  9.6402E-01  1.0273E+00  1.0257E+00  6.5728E-01  9.6747E-01  9.6010E-01
             1.0538E+00
 PARAMETER:  8.5674E-02  9.9157E-02  1.1125E-01  1.3311E-01  6.3358E-02  1.2692E-01  1.2540E-01 -3.1965E-01  6.6924E-02  5.9282E-02
             1.5243E-01
 GRADIENT:  -1.7029E+06 -1.7029E+06 -3.0613E+06  2.5585E+06 -1.7029E+06 -1.3417E+06  2.7160E+06  1.0653E+06  3.4057E+06  1.7029E+06
            -2.2347E+06

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1039
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8539E-02 -2.3318E-02 -1.5200E-02  1.9378E-02 -2.0904E-02
 SE:             2.9858E-02  2.3005E-02  1.5045E-02  2.8416E-02  2.4910E-02
 N:                     100         100         100         100         100

 P VAL.:         5.3466E-01  3.1077E-01  3.1235E-01  4.9529E-01  4.0138E-01

 ETASHRINKSD(%)  1.0000E-10  2.2931E+01  4.9597E+01  4.8013E+00  1.6549E+01
 ETASHRINKVR(%)  1.0000E-10  4.0604E+01  7.4596E+01  9.3720E+00  3.0359E+01
 EBVSHRINKSD(%)  2.6490E-01  2.2899E+01  5.5470E+01  8.2080E+00  1.6435E+01
 EBVSHRINKVR(%)  5.2910E-01  4.0555E+01  8.0171E+01  1.5742E+01  3.0170E+01
 RELATIVEINF(%)  9.9469E+01  2.7983E+01  1.0767E+01  5.7099E+01  2.4712E+01
 EPSSHRINKSD(%)  2.2926E+01
 EPSSHRINKVR(%)  4.0596E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3895.2667682239230     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2241.1774084555122     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.58
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3895.267       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  9.99E-01  1.01E+00  1.03E+00  9.64E-01  1.03E+00  1.03E+00  6.57E-01  9.67E-01  9.60E-01  1.05E+00
 


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
+        8.76E+09
 
 TH 2
+       -2.79E+04  8.53E+09
 
 TH 3
+       -3.39E+05  7.57E+09  6.73E+09
 
 TH 4
+        1.23E+04  2.03E+04  2.43E+05  4.50E+09
 
 TH 5
+        8.96E+09  8.84E+09 -3.47E+05  4.84E+04  9.16E+09
 
 TH 6
+       -2.78E+04  6.54E+09  5.80E+09  1.99E+04 -2.85E+04  5.01E+09
 
 TH 7
+       -6.72E+09 -6.63E+09  2.59E+05 -2.02E+04 -6.87E+09  2.13E+04  5.15E+09
 
 TH 8
+        8.04E+03  1.31E+04  1.59E+05 -8.92E+04  3.16E+04  1.31E+04 -1.32E+04  1.93E+09
 
 TH 9
+        1.75E+04  2.85E+04  3.45E+05 -1.94E+05  6.86E+04  2.84E+04 -2.88E+04 -1.34E+06  9.10E+09
 
 TH10
+       -9.00E+09  2.87E+04  3.48E+05 -4.01E+04 -9.20E+09  2.86E+04  6.89E+09 -2.63E+04 -5.71E+04  9.24E+09
 
 TH11
+       -1.05E+04 -1.72E+04 -2.08E+05  1.17E+05 -4.13E+04 -1.71E+04  1.73E+04  8.08E+05  2.85E+05  3.44E+04  3.30E+09
 
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
 #CPUT: Total CPU Time in Seconds,       39.251
Stop Time:
Sat Sep 18 06:09:24 CDT 2021
