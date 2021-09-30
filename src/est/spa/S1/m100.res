Wed Sep 29 14:46:06 CDT 2021
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
$DATA ../../../../data/spa/S1/dat100.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1666.66416882802        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0698E+02 -1.8143E+01 -2.9078E+01  5.3963E+01  4.5042E+01  3.3066E+01  1.1044E+00 -3.6651E-01  3.7558E+01  6.5130E+00
             1.7216E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1673.71039080608        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.6143E-01  1.0460E+00  9.6055E-01  9.9704E-01  9.5642E-01  1.0117E+00  1.1057E+00  1.0540E+00  7.8243E-01  9.4229E-01
             1.0232E+00
 PARAMETER:  6.0669E-02  1.4494E-01  5.9746E-02  9.7040E-02  5.5443E-02  1.1165E-01  2.0048E-01  1.5256E-01 -1.4535E-01  4.0562E-02
             1.2295E-01
 GRADIENT:   3.9670E+02  5.1717E+01 -2.8469E+00  8.8542E+01 -1.0020E+01  5.0301E+01  6.1114E-01  2.5247E+00  9.1008E+00  7.5751E+00
             2.5770E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1675.30349599338        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      247
 NPARAMETR:  9.6428E-01  1.1364E+00  1.0039E+00  9.3453E-01  1.0389E+00  1.0577E+00  1.2002E+00  1.0636E+00  7.7282E-01  9.7822E-01
             9.2845E-01
 PARAMETER:  6.3627E-02  2.2788E-01  1.0393E-01  3.2288E-02  1.3818E-01  1.5613E-01  2.8249E-01  1.6170E-01 -1.5770E-01  7.7979E-02
             2.5763E-02
 GRADIENT:   4.9245E+00  9.3361E-01  5.1701E-01  1.1535E+00  5.3108E+00  8.7096E+00  9.2567E+00 -3.5842E+00  5.8754E+00 -2.3114E+00
            -1.6227E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1676.67651636720        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      424
 NPARAMETR:  9.6415E-01  1.2035E+00  1.2026E+00  9.0193E-01  1.1407E+00  1.0338E+00  1.1129E+00  1.5127E+00  7.3896E-01  1.0661E+00
             9.6855E-01
 PARAMETER:  6.3487E-02  2.8519E-01  2.8447E-01 -3.2184E-03  2.3162E-01  1.3328E-01  2.0697E-01  5.1389E-01 -2.0251E-01  1.6403E-01
             6.8046E-02
 GRADIENT:   5.0826E+00  3.3410E+00  6.0071E-01  4.3055E+00 -3.7438E+00 -1.5532E-01  1.4538E+00  4.2735E-01  2.4002E-01  6.8835E-01
             1.1194E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1676.84644043770        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      602
 NPARAMETR:  9.6212E-01  1.3510E+00  1.1253E+00  8.0738E-01  1.1831E+00  1.0372E+00  1.0140E+00  1.6043E+00  7.7851E-01  1.0881E+00
             9.6476E-01
 PARAMETER:  6.1388E-02  4.0084E-01  2.1801E-01 -1.1396E-01  2.6813E-01  1.3656E-01  1.1388E-01  5.7266E-01 -1.5037E-01  1.8440E-01
             6.4125E-02
 GRADIENT:  -7.1872E-01  5.3705E+00  2.3889E+00  2.9754E+00 -3.5221E+00  7.5848E-01  3.1421E-01 -3.6333E-01 -6.9525E-01 -2.7259E-01
            -1.2700E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1676.96342907112        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      778
 NPARAMETR:  9.6201E-01  1.5164E+00  9.1514E-01  6.9829E-01  1.1898E+00  1.0376E+00  9.0804E-01  1.5449E+00  8.8147E-01  1.0764E+00
             9.6487E-01
 PARAMETER:  6.1274E-02  5.1634E-01  1.1321E-02 -2.5912E-01  2.7380E-01  1.3695E-01  3.5354E-03  5.3495E-01 -2.6161E-02  1.7358E-01
             6.4242E-02
 GRADIENT:  -3.1736E+00  6.1231E+00  2.5455E+00  3.0446E+00 -3.3131E+00  4.4677E-01 -1.2135E+00 -7.7667E-01 -3.1917E-01 -4.4286E-01
            -1.2008E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1677.02769711137        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      961             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6504E-01  1.6020E+00  7.7221E-01  6.3396E-01  1.1872E+00  1.0382E+00  8.7223E-01  1.4676E+00  9.4258E-01  1.0642E+00
             9.6750E-01
 PARAMETER:  6.4410E-02  5.7124E-01 -1.5849E-01 -3.5578E-01  2.7156E-01  1.3753E-01 -3.6697E-02  4.8362E-01  4.0869E-02  1.6220E-01
             6.6958E-02
 GRADIENT:   4.6053E+02  5.1127E+02  1.2515E+00  1.1727E+02  1.9047E+01  8.5192E+01  6.6847E+00  4.0967E-01  3.3663E+00  1.0879E+00
             1.2246E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1677.03492677998        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1144
 NPARAMETR:  9.6476E-01  1.6055E+00  7.7160E-01  6.3283E-01  1.1855E+00  1.0379E+00  8.7291E-01  1.4712E+00  9.2969E-01  1.0633E+00
             9.6694E-01
 PARAMETER:  6.4128E-02  5.7342E-01 -1.5929E-01 -3.5756E-01  2.7019E-01  1.3717E-01 -3.5926E-02  4.8604E-01  2.7096E-02  1.6139E-01
             6.6384E-02
 GRADIENT:   1.3849E+00 -4.0623E+00  1.3533E+00 -2.1815E+00 -9.3611E-01  3.1332E-01 -7.7500E-02 -5.7274E-01 -7.6416E-02 -3.6048E-01
            -5.3565E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1677.04033898882        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1334             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6496E-01  1.6047E+00  7.6677E-01  6.3310E-01  1.1851E+00  1.0381E+00  8.7355E-01  1.4699E+00  9.2864E-01  1.0649E+00
             9.6683E-01
 PARAMETER:  6.4326E-02  5.7295E-01 -1.6557E-01 -3.5713E-01  2.6986E-01  1.3739E-01 -3.5185E-02  4.8520E-01  2.5971E-02  1.6292E-01
             6.6262E-02
 GRADIENT:   4.6079E+02  5.1571E+02  1.3785E+00  1.1874E+02  1.7922E+01  8.5171E+01  6.3807E+00  4.6773E-01  2.5938E+00  1.3753E+00
             8.7254E-01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1677.04033898882        NO. OF FUNC. EVALS.:  62
 CUMULATIVE NO. OF FUNC. EVALS.:     1396
 NPARAMETR:  9.6498E-01  1.6051E+00  7.6711E-01  6.3339E-01  1.1860E+00  1.0381E+00  8.7361E-01  1.4680E+00  9.2772E-01  1.0654E+00
             9.6674E-01
 PARAMETER:  6.4326E-02  5.7295E-01 -1.6557E-01 -3.5713E-01  2.6986E-01  1.3739E-01 -3.5185E-02  4.8520E-01  2.5971E-02  1.6292E-01
             6.6262E-02
 GRADIENT:  -8.8859E-02 -9.5735E-01 -4.6062E+04 -7.9476E-01 -2.8274E+04 -1.8451E-02 -2.2933E-02  1.5597E+04  5.4510E-02 -4.6837E+04
             7.3035E-02
 NUMSIGDIG:         3.4         3.1         2.3         2.6         2.3         3.5         3.0         2.3         1.7         2.3
                    2.8

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1396
 NO. OF SIG. DIGITS IN FINAL EST.:  1.7
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9182E-04 -1.3971E-02 -3.6278E-02  1.2553E-02 -3.9522E-02
 SE:             2.9857E-02  2.4965E-02  1.2730E-02  1.8576E-02  2.1618E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9220E-01  5.7573E-01  4.3761E-03  4.9918E-01  6.7524E-02

 ETASHRINKSD(%)  1.0000E-10  1.6365E+01  5.7352E+01  3.7768E+01  2.7576E+01
 ETASHRINKVR(%)  1.0000E-10  3.0051E+01  8.1811E+01  6.1272E+01  4.7548E+01
 EBVSHRINKSD(%)  3.7755E-01  1.6196E+01  6.1791E+01  4.0684E+01  2.3667E+01
 EBVSHRINKVR(%)  7.5367E-01  2.9769E+01  8.5401E+01  6.4816E+01  4.1732E+01
 RELATIVEINF(%)  9.8931E+01  1.8877E+00  7.3388E-01  7.8362E-01  1.4898E+01
 EPSSHRINKSD(%)  4.4922E+01
 EPSSHRINKVR(%)  6.9664E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1677.0403389888197     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -941.88951242508153     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.22
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1677.040       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.65E-01  1.60E+00  7.67E-01  6.33E-01  1.19E+00  1.04E+00  8.74E-01  1.47E+00  9.29E-01  1.06E+00  9.67E-01
 


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
+        1.10E+03
 
 TH 2
+       -6.71E+00  2.26E+05
 
 TH 3
+       -2.20E+02  1.44E+03  1.18E+07
 
 TH 4
+       -1.40E+01 -9.18E+05  7.46E+03  3.73E+06
 
 TH 5
+       -9.14E+01  4.54E+02 -7.86E+02  2.24E+02  1.87E+06
 
 TH 6
+        2.03E-01 -1.73E+00 -2.18E+02 -3.80E+00 -8.82E+01  1.82E+02
 
 TH 7
+        1.65E+00  1.80E+01 -1.42E+03 -2.48E+01 -5.68E+02  2.36E-01  1.35E+02
 
 TH 8
+        4.00E+01 -2.56E+02 -3.50E+04  1.94E+03  8.30E+05  3.89E+01  2.55E+02  3.81E+05
 
 TH 9
+        5.99E-01 -1.00E+01  4.37E+03  1.99E+01  1.75E+03 -9.31E-02  3.44E+01 -7.60E+02  2.21E+07
 
 TH10
+       -1.65E+02  1.01E+03  1.29E+02  9.17E+00  3.44E+06 -1.60E+02 -1.04E+03 -1.94E+01  3.22E+03  6.34E+06
 
 TH11
+       -6.99E+00 -1.11E+01  5.03E+03  6.54E+00  2.00E+03  1.56E+00  7.93E+00 -8.87E+02 -2.12E+07  3.70E+03  2.04E+07
 
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
 #CPUT: Total CPU Time in Seconds,       25.764
Stop Time:
Wed Sep 29 14:46:33 CDT 2021
