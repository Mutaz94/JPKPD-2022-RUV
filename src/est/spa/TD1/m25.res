Sat Sep 18 13:56:49 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1700.02689517432        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.6029E+01 -7.7460E+01 -1.5777E+01 -1.1244E+02 -4.1972E+01  3.0232E+01 -3.9093E+00  1.8048E+01  2.9479E+00  2.9523E+01
             7.1503E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1710.83294953329        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0039E+00  1.1404E+00  1.2428E+00  1.0118E+00  1.1895E+00  8.9211E-01  1.0415E+00  8.5009E-01  9.5606E-01  8.4058E-01
             1.0284E+00
 PARAMETER:  1.0392E-01  2.3137E-01  3.1740E-01  1.1177E-01  2.7356E-01 -1.4162E-02  1.4071E-01 -6.2418E-02  5.5065E-02 -7.3668E-02
             1.2796E-01
 GRADIENT:   7.0368E+01  2.6198E+01  1.8402E+01  1.6653E+01  1.5563E+01 -1.3406E+01 -9.3917E-01 -4.2317E+00 -7.3942E+00 -2.0574E+01
             9.9944E-02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1714.22949724439        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.9612E-01  1.0503E+00  9.7939E-01  1.0436E+00  1.0447E+00  9.2734E-01  1.1245E+00  3.2389E-01  9.2545E-01  9.1439E-01
             9.6673E-01
 PARAMETER:  9.6117E-02  1.4906E-01  7.9170E-02  1.4272E-01  1.4374E-01  2.4567E-02  2.1730E-01 -1.0273E+00  2.2521E-02  1.0501E-02
             6.6162E-02
 GRADIENT:   4.9855E+01 -2.9378E+00 -3.0658E+00  1.1349E+01  1.6312E+01  2.5690E+00 -4.1427E-01 -1.9970E-01 -3.8423E+00  2.8928E+00
            -1.3121E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1714.94098808158        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      272
 NPARAMETR:  9.9042E-01  1.1010E+00  9.1548E-01  1.0082E+00  1.0196E+00  9.3110E-01  1.0987E+00  3.0971E-01  9.5809E-01  8.6052E-01
             9.9457E-01
 PARAMETER:  9.0371E-02  1.9620E-01  1.1689E-02  1.0820E-01  1.1942E-01  2.8607E-02  1.9416E-01 -1.0721E+00  5.7181E-02 -5.0216E-02
             9.4553E-02
 GRADIENT:  -9.6725E+00 -4.0598E+00 -1.7882E+00 -1.4790E+00 -1.1304E+00  3.7171E-01  7.0955E-02  5.3347E-01 -4.6850E-01  2.6975E+00
             7.0812E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1715.55795210488        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  9.9483E-01  1.3881E+00  8.1359E-01  8.2888E-01  1.1234E+00  9.3145E-01  9.2397E-01  2.1678E-01  1.1105E+00  8.8769E-01
             9.9528E-01
 PARAMETER:  9.4815E-02  4.2791E-01 -1.0630E-01 -8.7684E-02  2.1640E-01  2.8990E-02  2.0925E-02 -1.4289E+00  2.0481E-01 -1.9132E-02
             9.5265E-02
 GRADIENT:  -1.7681E+00  2.1450E+00 -1.8009E+00  3.5405E+00 -1.8001E-01  3.5060E-02 -2.8576E-01  3.5153E-01 -4.8035E-01  4.5870E-01
             5.0483E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1715.67974531922        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      626
 NPARAMETR:  9.9554E-01  1.5407E+00  7.6294E-01  7.2698E-01  1.1964E+00  9.3193E-01  8.5688E-01  1.3160E-01  1.2273E+00  9.2938E-01
             9.9706E-01
 PARAMETER:  9.5534E-02  5.3225E-01 -1.7057E-01 -2.1885E-01  2.7932E-01  2.9500E-02 -5.4458E-02 -1.9280E+00  3.0481E-01  2.6761E-02
             9.7053E-02
 GRADIENT:  -7.3147E-01 -6.8322E-01 -6.5957E-01 -6.1848E-01 -8.5533E-02  1.4471E-01  3.6150E-02  1.2497E-01  7.0784E-02  3.8824E-01
             1.1577E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1715.75402759772        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      804
 NPARAMETR:  9.9573E-01  1.4682E+00  7.8637E-01  7.7408E-01  1.1610E+00  9.3139E-01  8.8839E-01  4.4976E-02  1.1747E+00  9.0741E-01
             9.9603E-01
 PARAMETER:  9.5720E-02  4.8407E-01 -1.4033E-01 -1.5609E-01  2.4932E-01  2.8919E-02 -1.8342E-02 -3.0016E+00  2.6098E-01  2.8365E-03
             9.6020E-02
 GRADIENT:   1.2918E-01 -3.9582E-01 -1.8650E-01 -7.9768E-02 -1.2063E-01 -4.1978E-02  1.4500E-01  1.4779E-02  3.4673E-01  1.8874E-02
             1.0685E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1715.76143509995        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  9.9568E-01  1.4729E+00  7.8572E-01  7.7135E-01  1.1639E+00  9.3149E-01  8.8581E-01  1.0000E-02  1.1762E+00  9.0978E-01
             9.9593E-01
 PARAMETER:  9.5673E-02  4.8724E-01 -1.4115E-01 -1.5961E-01  2.5174E-01  2.9032E-02 -2.1249E-02 -4.6207E+00  2.6228E-01  5.4427E-03
             9.5917E-02
 GRADIENT:  -8.1915E-03  7.6849E-02  2.3195E-02  2.9608E-02 -3.5793E-02  2.8369E-03 -6.2368E-03  0.0000E+00 -1.8487E-03  1.5369E-03
            -3.4924E-03

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1715.76144848737        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  9.9568E-01  1.4722E+00  7.8588E-01  7.7179E-01  1.1635E+00  9.3148E-01  8.8617E-01  1.0000E-02  1.1757E+00  9.0952E-01
             9.9592E-01
 PARAMETER:  9.5675E-02  4.8673E-01 -1.4096E-01 -1.5905E-01  2.5143E-01  2.9022E-02 -2.0841E-02 -4.5969E+00  2.6183E-01  5.1567E-03
             9.5908E-02
 GRADIENT:   9.7192E-04 -2.0454E-03 -1.1338E-03 -1.9236E-04  1.6027E-03 -5.3597E-04 -4.8730E-04  0.0000E+00  4.8484E-04 -4.5348E-04
            -9.4483E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1071
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.1124E-04 -1.9961E-02 -3.3779E-04  1.2994E-02 -3.1793E-02
 SE:             2.9814E-02  2.2307E-02  1.3480E-04  2.3468E-02  2.1836E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9702E-01  3.7086E-01  1.2218E-02  5.7980E-01  1.4540E-01

 ETASHRINKSD(%)  1.1984E-01  2.5270E+01  9.9548E+01  2.1379E+01  2.6845E+01
 ETASHRINKVR(%)  2.3954E-01  4.4154E+01  9.9998E+01  3.8188E+01  4.6484E+01
 EBVSHRINKSD(%)  4.8587E-01  2.4401E+01  9.9578E+01  2.2268E+01  2.5860E+01
 EBVSHRINKVR(%)  9.6938E-01  4.2848E+01  9.9998E+01  3.9578E+01  4.5033E+01
 RELATIVEINF(%)  9.8845E+01  3.0571E+00  2.4665E-04  3.6983E+00  8.4107E+00
 EPSSHRINKSD(%)  4.2642E+01
 EPSSHRINKVR(%)  6.7101E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1715.7614484873675     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -980.61062192362931     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.55
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1715.761       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  1.47E+00  7.86E-01  7.72E-01  1.16E+00  9.31E-01  8.86E-01  1.00E-02  1.18E+00  9.10E-01  9.96E-01
 


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
+        1.28E+03
 
 TH 2
+       -6.41E+00  3.91E+02
 
 TH 3
+        9.05E+00  1.64E+02  3.62E+02
 
 TH 4
+       -1.31E+01  3.23E+02 -2.09E+02  8.24E+02
 
 TH 5
+       -3.53E+00 -2.11E+02 -3.29E+02  2.18E+02  5.46E+02
 
 TH 6
+        1.45E+00 -1.01E+00  2.13E+00 -2.87E+00 -5.63E-01  2.26E+02
 
 TH 7
+       -1.93E-01  1.95E+01 -2.06E+00 -1.31E+01 -9.67E+00 -6.11E-03  8.99E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.67E+00 -2.22E+01 -2.78E+01  4.42E+01 -1.33E+00 -9.08E-01  1.77E+01  0.00E+00  6.58E+01
 
 TH10
+       -1.25E+00 -1.03E+01 -3.18E+01 -1.31E+01 -6.25E+01 -2.02E+00  1.59E+01  0.00E+00  5.64E+00  8.36E+01
 
 TH11
+       -1.16E+01 -1.96E+01 -3.21E+01  2.75E+00  1.02E+00  3.75E+00  1.32E+01  0.00E+00  6.88E+00  2.54E+01  2.25E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.627
Stop Time:
Sat Sep 18 13:57:09 CDT 2021
