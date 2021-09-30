Thu Sep 30 01:27:44 CDT 2021
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
$DATA ../../../../data/spa1/TD1/dat51.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m51.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1604.52498531318        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1582E+02  2.1954E+01  3.8751E+01  7.6794E+01  9.5340E+01  6.1581E+01  7.9459E+00 -2.1206E+02 -3.0651E+01  1.1112E+01
            -7.0946E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1958.77816324183        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:      115
 NPARAMETR:  7.2602E-01  9.8028E-01  9.7399E-01  9.3373E-01  9.2682E-01  1.0156E+00  9.8810E-01  1.1894E+00  1.0239E+00  9.7263E-01
             1.5405E+00
 PARAMETER: -2.2018E-01  8.0082E-02  7.3644E-02  3.1436E-02  2.4008E-02  1.1551E-01  8.8026E-02  2.7347E-01  1.2362E-01  7.2251E-02
             5.3213E-01
 GRADIENT:  -5.8946E+02 -5.6169E+01 -9.1180E+00 -9.0119E+01 -1.4378E+01 -1.3462E+02  5.9248E+00  9.4870E+00 -2.9778E+00  3.0175E+01
             2.8073E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2020.42398446133        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  8.0832E-01  8.4569E-01  6.3758E-01  1.1121E+00  7.0304E-01  9.2541E-01  5.8638E-01  3.4503E-01  1.1425E+00  9.7346E-01
             1.3766E+00
 PARAMETER: -1.1280E-01 -6.7605E-02 -3.5007E-01  2.0627E-01 -2.5234E-01  2.2487E-02 -4.3379E-01 -9.6411E-01  2.3320E-01  7.3104E-02
             4.1960E-01
 GRADIENT:  -4.4135E+02  6.1268E+00 -6.1823E+01  1.1802E+02  2.3861E+01 -8.2702E+01 -7.3178E-01  5.2147E-01  1.8908E+01  3.4172E+01
             2.2889E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2083.28684576171        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      474
 NPARAMETR:  8.7429E-01  1.1521E+00  8.5884E-01  1.0332E+00  9.3512E-01  8.7310E-01  1.3260E+00  5.2692E-01  1.1270E+00  8.5128E-01
             1.0795E+00
 PARAMETER: -3.4345E-02  2.4158E-01 -5.2169E-02  1.3261E-01  3.2919E-02 -3.5704E-02  3.8220E-01 -5.4071E-01  2.1958E-01 -6.1010E-02
             1.7650E-01
 GRADIENT:  -2.6264E+02  1.1691E+02 -5.0204E-01  1.7080E+02  4.9532E+01 -5.6534E+01  2.4299E+01 -1.0944E+01  2.3579E+00 -2.7703E+00
             1.0698E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2085.20015109061        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      658
 NPARAMETR:  8.7475E-01  1.1429E+00  8.7310E-01  1.0396E+00  9.3666E-01  8.7302E-01  1.3211E+00  5.7996E-01  1.1173E+00  8.5461E-01
             1.0713E+00
 PARAMETER: -3.3813E-02  2.3358E-01 -3.5711E-02  1.3887E-01  3.4565E-02 -3.5800E-02  3.7847E-01 -4.4479E-01  2.1095E-01 -5.7110E-02
             1.6883E-01
 GRADIENT:  -2.6070E+02  1.1768E+02 -2.4439E+00  1.7270E+02  4.7429E+01 -5.6214E+01  2.3014E+01 -1.0902E+01  1.7088E+00 -2.1443E+00
             1.0354E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2086.39853846525        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      843
 NPARAMETR:  8.7521E-01  1.1370E+00  8.8322E-01  1.0434E+00  9.3822E-01  8.7295E-01  1.3188E+00  6.1427E-01  1.1115E+00  8.5672E-01
             1.0667E+00
 PARAMETER: -3.3287E-02  2.2836E-01 -2.4183E-02  1.4248E-01  3.6228E-02 -3.5873E-02  3.7673E-01 -3.8731E-01  2.0573E-01 -5.4647E-02
             1.6453E-01
 GRADIENT:  -2.5893E+02  1.1773E+02 -3.4768E+00  1.7323E+02  4.6229E+01 -5.5914E+01  2.2349E+01 -1.0795E+01  1.4506E+00 -1.7963E+00
             1.0154E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2101.37521491147        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1012
 NPARAMETR:  9.2723E-01  1.1364E+00  8.8401E-01  1.0436E+00  9.3834E-01  9.4482E-01  1.3186E+00  6.1688E-01  1.1111E+00  8.8989E-01
             1.0291E+00
 PARAMETER:  2.4449E-02  2.2789E-01 -2.3282E-02  1.4272E-01  3.6361E-02  4.3244E-02  3.7654E-01 -3.8308E-01  2.0532E-01 -1.6656E-02
             1.2864E-01
 GRADIENT:  -7.7569E+01  1.1697E+02 -2.1788E+00  1.7290E+02  4.0973E+01 -1.8145E+00  2.2619E+01 -1.0808E+01  2.0122E+00  2.3710E+00
             7.9695E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2106.40972717790        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1197             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5893E-01  1.1357E+00  8.8426E-01  1.0432E+00  9.3808E-01  9.4596E-01  1.3171E+00  6.1621E-01  1.1104E+00  8.8513E-01
             9.6825E-01
 PARAMETER:  5.8063E-02  2.2723E-01 -2.3000E-02  1.4230E-01  3.6079E-02  4.4443E-02  3.7544E-01 -3.8416E-01  2.0473E-01 -2.2016E-02
             6.7740E-02
 GRADIENT:   4.7101E+02  2.3077E+02  3.1487E+00  2.9442E+02  5.2167E+01  6.2985E+01  4.6211E+01 -1.1654E+01  2.4528E+01  5.5602E-01
             3.8976E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2106.79259885291        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1361            RESET HESSIAN, TYPE II
 NPARAMETR:  9.5657E-01  1.1345E+00  8.8426E-01  1.0423E+00  9.3792E-01  9.4311E-01  1.3153E+00  6.1859E-01  1.1104E+00  8.8549E-01
             9.6791E-01
 PARAMETER:  5.5598E-02  2.2619E-01 -2.3009E-02  1.4139E-01  3.5905E-02  4.1425E-02  3.7408E-01 -3.8032E-01  2.0468E-01 -2.1613E-02
             6.7383E-02
 GRADIENT:   4.6527E+02  2.2929E+02  3.0241E+00  2.9216E+02  5.1895E+01  6.2264E+01  4.6110E+01 -1.1569E+01  2.4771E+01  7.2052E-01
             3.8877E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2107.85050122191        NO. OF FUNC. EVALS.: 108
 CUMULATIVE NO. OF FUNC. EVALS.:     1469
 NPARAMETR:  9.5380E-01  1.1308E+00  8.8424E-01  1.0398E+00  9.3733E-01  9.4026E-01  1.3082E+00  6.2430E-01  1.1074E+00  8.8482E-01
             9.6698E-01
 PARAMETER:  5.2696E-02  2.2290E-01 -2.3028E-02  1.3904E-01  3.5277E-02  3.8404E-02  3.6867E-01 -3.7112E-01  2.0199E-01 -2.2367E-02
             6.6425E-02
 GRADIENT:   6.9989E+00  1.0860E+02 -2.7326E+05  9.8439E+04  5.6693E+01  1.5399E+01  2.2166E+01  7.3542E+04 -1.1663E+02  1.3666E+05
            -7.5245E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1469
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.8855E-03 -7.2881E-02 -2.3095E-02 -8.1429E-02 -3.8267E-02
 SE:             2.9994E-02  1.7965E-02  1.2334E-02  2.3067E-02  2.1971E-02
 N:                     100         100         100         100         100

 P VAL.:         8.9693E-01  4.9788E-05  6.1146E-02  4.1541E-04  8.1559E-02

 ETASHRINKSD(%)  1.0000E-10  3.9814E+01  5.8678E+01  2.2723E+01  2.6394E+01
 ETASHRINKVR(%)  1.0000E-10  6.3776E+01  8.2925E+01  4.0283E+01  4.5822E+01
 EBVSHRINKSD(%)  3.5241E-01  2.7719E+01  6.9865E+01  1.7392E+01  2.4692E+01
 EBVSHRINKVR(%)  7.0359E-01  4.7755E+01  9.0919E+01  3.1760E+01  4.3288E+01
 RELATIVEINF(%)  9.8884E+01  4.3541E+00  1.5725E+00  6.5913E+00  1.0944E+01
 EPSSHRINKSD(%)  3.7282E+01
 EPSSHRINKVR(%)  6.0664E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2107.8505012219098     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1188.9119680172371     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.92
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2107.851       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.54E-01  1.13E+00  8.84E-01  1.04E+00  9.37E-01  9.40E-01  1.31E+00  6.24E-01  1.11E+00  8.85E-01  9.67E-01
 


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
+        7.51E+07
 
 TH 2
+        2.84E+07  1.08E+07
 
 TH 3
+       -4.52E+03  3.07E+07  8.74E+07
 
 TH 4
+        2.77E+03 -1.88E+07  5.35E+07  5.94E+02
 
 TH 5
+        7.64E+07 -2.23E+02  8.25E+07 -1.31E+04  1.56E+08
 
 TH 6
+       -4.75E+03 -2.88E+07 -5.12E+03  3.13E+03  7.75E+07  1.55E+08
 
 TH 7
+        1.49E+07 -2.17E+01 -7.69E+01 -9.80E+06  3.41E+01  4.35E-01  2.94E+06
 
 TH 8
+        1.71E+03 -5.13E+02  3.34E+07  2.04E+07  1.64E+03  1.94E+03  2.34E+01  1.27E+07
 
 TH 9
+        6.41E+07 -5.74E+01  3.46E+07 -2.11E+07  3.26E+07  1.05E+01  6.34E+06 -1.14E+04  1.36E+07
 
 TH10
+        1.29E+00 -1.31E+03 -8.73E+07  5.34E+07 -1.12E+02 -1.96E+03 -1.60E+07 -7.75E+02  2.09E+01  8.73E+07
 
 TH11
+       -9.00E+00 -8.10E+04 -3.09E+04  4.88E+07  7.53E+07  1.58E+00  1.46E+07 -8.81E+04  3.16E+07 -8.00E+07  7.33E+07
 
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
 #CPUT: Total CPU Time in Seconds,       30.545
Stop Time:
Thu Sep 30 01:28:16 CDT 2021
