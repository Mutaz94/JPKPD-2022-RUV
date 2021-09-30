Thu Sep 30 09:35:09 CDT 2021
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
$DATA ../../../../data/spa2/D/dat70.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   18914.5926863039        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.4912E+02  3.5099E+02 -2.4417E+01  1.2964E+02  3.4456E+02 -2.6164E+03 -1.1234E+03 -2.1592E+01 -1.8793E+03 -9.5029E+02
            -3.5835E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -655.244385546827        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.3169E+00  1.2339E+00  8.6770E-01  2.3805E+00  1.0884E+00  3.0339E+00  1.9748E+00  9.3696E-01  2.5524E+00  1.2158E+00
             1.2105E+01
 PARAMETER:  3.7525E-01  3.1020E-01 -4.1904E-02  9.6731E-01  1.8471E-01  1.2099E+00  7.8046E-01  3.4888E-02  1.0370E+00  2.9544E-01
             2.5936E+00
 GRADIENT:  -1.8367E+01  1.6169E+01 -6.0327E+01  1.1378E+02  3.7748E+01  7.7428E+01 -5.9841E+00  4.6404E+00 -3.2601E+01  1.5100E+01
             1.4490E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -750.582126528784        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.2832E+00  2.4065E+00  4.1248E+00  2.2508E+00  3.3718E+00  3.5887E+00  5.2121E+00  5.4987E-01  6.9587E+00  1.8399E+00
             1.0019E+01
 PARAMETER:  3.4935E-01  9.7817E-01  1.5170E+00  9.1128E-01  1.3154E+00  1.3778E+00  1.7510E+00 -4.9807E-01  2.0400E+00  7.0970E-01
             2.4045E+00
 GRADIENT:  -3.8149E+00  2.0276E+01 -1.7489E+01  5.8521E+01  2.4523E-01  1.1820E+02  8.3771E+01  3.8977E-02  7.2065E+01  1.2945E+01
             1.8827E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -786.594984504197        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.1819E+00  1.4413E+00  3.5258E+00  1.2163E+00  1.6189E+00  2.8723E+00  3.6944E+00  4.1655E-02  2.5919E+00  1.7276E+00
             1.0181E+01
 PARAMETER:  2.6715E-01  4.6558E-01  1.3601E+00  2.9581E-01  5.8177E-01  1.1551E+00  1.4068E+00 -3.0783E+00  1.0524E+00  6.4673E-01
             2.4205E+00
 GRADIENT:  -3.9969E+01 -2.3455E+01  3.4602E-01  1.1305E+01 -3.2281E+01  3.5762E+01 -4.7110E+01  7.0548E-04  3.1556E+01  3.4173E+01
             2.1297E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -845.109795631872        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  1.2797E+00  2.3592E+00  6.0615E+00  5.2275E-01  2.0996E+00  3.1282E+00  3.7460E+00  8.3214E-02  1.9065E+00  5.3527E-01
             8.1402E+00
 PARAMETER:  3.4661E-01  9.5832E-01  1.9020E+00 -5.4865E-01  8.4175E-01  1.2404E+00  1.4207E+00 -2.3863E+00  7.4529E-01 -5.2498E-01
             2.1968E+00
 GRADIENT:  -1.5425E+01 -1.0235E+01  7.1035E-01 -1.5056E+01 -9.3452E+00 -7.8291E+00 -1.0522E+01  4.1379E-05  8.3988E+00  3.4110E+00
            -5.6566E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -867.016480470450        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  1.3960E+00  8.5316E-01  2.1705E+01  1.2884E+00  2.2822E+00  3.2471E+00  6.5455E+00  1.0000E-02  9.9671E-01  3.5265E-01
             8.0840E+00
 PARAMETER:  4.3364E-01 -5.8807E-02  3.1775E+00  3.5339E-01  9.2516E-01  1.2777E+00  1.9788E+00 -5.0634E+00  9.6707E-02 -9.4228E-01
             2.1899E+00
 GRADIENT:   3.7610E+00 -2.8916E+00  4.8455E-01  1.2085E+01  2.6173E+00  6.1509E+00 -1.6605E+01  0.0000E+00  4.0979E+00  1.2637E+00
            -4.7932E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -871.191251879727        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      716             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3204E+00  7.7742E-01  1.9773E+01  1.2929E+00  2.1904E+00  3.0446E+00  7.3665E+00  1.0000E-02  7.0871E-01  1.1429E-01
             8.5688E+00
 PARAMETER:  3.7790E-01 -1.5178E-01  3.0843E+00  3.5688E-01  8.8408E-01  1.2134E+00  2.0969E+00 -6.2126E+00 -2.4431E-01 -2.0691E+00
             2.2481E+00
 GRADIENT:   2.5400E+01  3.3103E-01  9.2210E-01  2.0890E+01 -5.2043E+00  7.8434E+01  1.4929E+02  0.0000E+00 -9.8606E-01  1.3529E-01
             3.1138E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -871.889021133927        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      862
 NPARAMETR:  1.3198E+00  7.7760E-01  9.9346E+00  1.2930E+00  2.1897E+00  3.0320E+00  7.2742E+00  1.0000E-02  8.0115E-01  5.8636E-02
             8.5463E+00
 PARAMETER:  3.7744E-01 -1.5154E-01  2.3960E+00  3.5697E-01  8.8378E-01  1.2092E+00  2.0843E+00 -6.2126E+00 -1.2171E-01 -2.7364E+00
             2.2455E+00
 GRADIENT:  -5.2949E+00 -2.2632E+00  7.6882E-01 -3.8070E-01  2.9972E+00  8.5157E-01  1.2004E+00  0.0000E+00  2.5597E+00  3.4228E-02
             8.0636E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -872.811272821933        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1042             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3514E+00  8.3255E-01  6.7753E+00  1.2777E+00  2.0780E+00  3.0655E+00  7.1593E+00  1.0000E-02  7.2330E-01  2.4903E-02
             8.5349E+00
 PARAMETER:  4.0113E-01 -8.3260E-02  2.0133E+00  3.4508E-01  8.3139E-01  1.2202E+00  2.0684E+00 -6.2126E+00 -2.2393E-01 -3.5927E+00
             2.2442E+00
 GRADIENT:   3.2829E+01  7.3428E-01  1.2452E+00  2.0209E+01  6.7381E-01  8.0872E+01  1.4118E+02  0.0000E+00  2.2883E-01  6.8097E-03
             2.5488E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -873.930561200442        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1115
 NPARAMETR:  1.3885E+00  9.7719E-01  2.6253E+00  1.2075E+00  1.6770E+00  3.1072E+00  6.9687E+00  1.0000E-02  6.1972E-01  1.0000E-02
             8.5384E+00
 PARAMETER:  4.2825E-01  7.6924E-02  1.0652E+00  2.8854E-01  6.1698E-01  1.2337E+00  2.0414E+00 -6.2126E+00 -3.7848E-01 -5.8964E+00
             2.2446E+00
 GRADIENT:   4.1621E+01  4.1376E+00  5.9806E-01  1.3488E+01 -2.8585E+00  8.7454E+01  1.3727E+02  0.0000E+00  1.3840E+00  0.0000E+00
             2.2906E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -874.086432428559        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  1.3870E+00  9.5015E-01  2.5638E+00  1.1951E+00  1.6639E+00  3.1087E+00  6.9644E+00  1.0000E-02  5.8218E-01  1.0000E-02
             8.5579E+00
 PARAMETER:  4.2716E-01  4.8862E-02  1.0415E+00  2.7825E-01  6.0917E-01  1.2342E+00  2.0408E+00 -6.2126E+00 -4.4098E-01 -6.2916E+00
             2.2469E+00
 GRADIENT:   6.9079E+00  1.2315E+00  5.5456E-01  7.8878E-01 -3.9714E+00  1.0028E+01  7.5764E+00  0.0000E+00  9.5260E-01  0.0000E+00
             3.0985E-01

0ITERATION NO.:   53    OBJECTIVE VALUE:  -874.163010808364        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     1357
 NPARAMETR:  1.3709E+00  9.3007E-01  2.6088E+00  1.1913E+00  1.6926E+00  3.0886E+00  7.0236E+00  1.0000E-02  5.3984E-01  1.0000E-02
             8.5711E+00
 PARAMETER:  4.1555E-01  2.8510E-02  1.0593E+00  2.7498E-01  6.2569E-01  1.2280E+00  2.0489E+00 -6.2126E+00 -5.1410E-01 -6.6381E+00
             2.2485E+00
 GRADIENT:   1.2206E+02  3.6798E-01  1.3359E-01 -5.9768E-01 -1.2944E+00  1.9654E+00 -3.4890E+01  0.0000E+00  3.4732E-01  0.0000E+00
             2.9003E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1357
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.6383E-03  3.1012E-02 -3.7014E-06 -5.3138E-02  4.8804E-05
 SE:             2.8827E-02  2.3346E-02  1.9770E-05  9.8264E-03  6.1957E-05
 N:                     100         100         100         100         100

 P VAL.:         8.7217E-01  1.8406E-01  8.5149E-01  6.3990E-08  4.3087E-01

 ETASHRINKSD(%)  3.4244E+00  2.1789E+01  9.9934E+01  6.7080E+01  9.9792E+01
 ETASHRINKVR(%)  6.7315E+00  3.8830E+01  1.0000E+02  8.9163E+01  1.0000E+02
 EBVSHRINKSD(%)  1.9730E+00  1.1066E+01  9.9893E+01  7.4037E+01  9.9680E+01
 EBVSHRINKVR(%)  3.9071E+00  2.0908E+01  1.0000E+02  9.3259E+01  9.9999E+01
 RELATIVEINF(%)  9.5813E+01  4.8025E+01  2.3428E-05  4.1425E+00  2.1194E-04
 EPSSHRINKSD(%)  8.5522E+00
 EPSSHRINKVR(%)  1.6373E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -874.16301080836365     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       228.56322903724345     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.51
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -874.163       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.37E+00  9.31E-01  2.61E+00  1.19E+00  1.69E+00  3.09E+00  7.02E+00  1.00E-02  5.41E-01  1.00E-02  8.57E+00
 


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
+        7.54E+03
 
 TH 2
+        6.37E+00  2.11E+01
 
 TH 3
+       -4.63E-01  7.28E-01  2.13E+00
 
 TH 4
+       -1.32E+04  2.87E+01  1.46E-02  2.35E+04
 
 TH 5
+        9.68E+01 -5.36E+00 -7.68E+00 -1.15E+01  3.78E+01
 
 TH 6
+       -1.28E+03 -1.77E-03  1.41E-01  2.26E+03 -1.07E+00  2.35E+02
 
 TH 7
+       -2.96E+02  1.83E+00 -2.20E-01  5.10E+02 -2.90E+00  5.05E+01  1.38E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -3.48E-01 -1.71E+00 -2.12E+00 -4.94E+01  6.97E+00  2.69E-01  1.63E+00  0.00E+00  2.42E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.93E+02 -2.17E+00 -1.12E-01 -1.81E+01  1.44E+00  8.38E-01  8.13E+00  0.00E+00  4.89E+00  0.00E+00  9.27E+00
 
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
 #CPUT: Total CPU Time in Seconds,       40.026
Stop Time:
Thu Sep 30 09:35:50 CDT 2021
