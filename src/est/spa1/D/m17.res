Thu Sep 30 02:44:55 CDT 2021
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
$DATA ../../../../data/spa1/D/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1234.02510772566        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.0735E+01 -1.1772E+02 -3.5557E+01 -3.0349E+02  1.9964E+02 -5.3022E+02 -2.7980E+02 -1.7732E+01 -5.8968E+02 -2.3439E+02
            -1.0292E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1716.14531296008        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      108
 NPARAMETR:  8.6560E-01  8.8611E-01  9.8863E-01  1.1920E+00  1.2030E+00  1.3554E+00  2.7820E+00  1.0363E+00  9.6787E-01  2.6991E+00
             9.9123E-01
 PARAMETER: -4.4336E-02 -2.0917E-02  8.8568E-02  2.7563E-01  2.8484E-01  4.0411E-01  1.1232E+00  1.3564E-01  6.7346E-02  1.0929E+00
             9.1191E-02
 GRADIENT:  -3.6489E+02 -5.0218E+01 -9.8545E+01  1.3308E+02  5.0481E+00 -3.7069E+02 -2.0584E+02  1.2059E+01 -4.8617E+01  9.7554E+01
            -5.0349E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1831.67278785530        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      289
 NPARAMETR:  1.0814E+00  1.2102E+00  1.2353E+00  1.0029E+00  1.2087E+00  1.2875E+00  3.5054E+00  9.7218E-01  1.1311E+00  1.3629E+00
             1.0818E+00
 PARAMETER:  1.7828E-01  2.9080E-01  3.1133E-01  1.0288E-01  2.8952E-01  3.5270E-01  1.3543E+00  7.1785E-02  2.2322E-01  4.0963E-01
             1.7865E-01
 GRADIENT:  -1.0774E+02 -7.9516E-01  2.5966E+01 -3.6630E+01 -3.4368E+01 -3.5060E+02 -6.3097E+01 -5.9128E+00 -1.7678E+01  4.5539E+00
             1.4754E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1899.11686523269        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      466
 NPARAMETR:  1.1570E+00  1.2052E+00  1.0635E+00  9.2673E-01  1.2316E+00  1.6952E+00  3.6986E+00  1.0607E+00  1.2028E+00  1.2430E+00
             1.0556E+00
 PARAMETER:  2.4587E-01  2.8665E-01  1.6160E-01  2.3902E-02  3.0830E-01  6.2777E-01  1.4080E+00  1.5893E-01  2.8466E-01  3.1756E-01
             1.5409E-01
 GRADIENT:  -4.3199E+00 -1.4941E+01 -1.0647E+01 -4.3760E+01  9.7084E+00 -6.1144E+01 -4.9114E+01  5.3486E+00 -4.0777E+00 -8.3158E+00
            -2.0716E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1916.88573429191        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:      637
 NPARAMETR:  1.2058E+00  8.7780E-01  1.4335E+00  1.1988E+00  1.2089E+00  2.0772E+00  4.3169E+00  1.2163E+00  1.1732E+00  1.2971E+00
             1.0600E+00
 PARAMETER:  2.8712E-01 -3.0339E-02  4.6014E-01  2.8133E-01  2.8969E-01  8.3101E-01  1.5625E+00  2.9579E-01  2.5976E-01  3.6014E-01
             1.5827E-01
 GRADIENT:   2.6394E+01 -2.2126E+00  7.0487E-01 -2.3551E+00  3.5385E+00  8.5215E+01 -1.0412E+02  8.5316E-01 -6.2846E+00 -4.6288E+00
             2.3643E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1917.36272559358        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  1.1853E+00  8.9781E-01  1.4284E+00  1.1977E+00  1.2058E+00  2.0894E+00  4.3306E+00  1.1912E+00  1.2028E+00  1.3390E+00
             1.0593E+00
 PARAMETER:  2.7004E-01 -7.7954E-03  4.5656E-01  2.8044E-01  2.8715E-01  8.3687E-01  1.5657E+00  2.7493E-01  2.8465E-01  3.9190E-01
             1.5761E-01
 GRADIENT:   1.0758E+03  4.3771E+01  6.9600E+00  2.6913E+02  1.1915E+01  2.0737E+03  1.7041E+03  7.3989E-01  3.4795E+01  5.4372E+00
             2.3854E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1917.67461419841        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  1.1603E+00  8.9978E-01  1.4218E+00  1.2000E+00  1.2068E+00  2.1835E+00  4.3289E+00  1.1774E+00  1.2135E+00  1.3346E+00
             1.0580E+00
 PARAMETER:  2.4870E-01 -5.6066E-03  4.5196E-01  2.8229E-01  2.8794E-01  8.8091E-01  1.5653E+00  2.6332E-01  2.9350E-01  3.8861E-01
             1.5641E-01
 GRADIENT:   9.8208E+02  4.4110E+01  5.8760E+00  2.7408E+02  1.4781E+01  2.2535E+03  1.7033E+03  3.5406E-01  3.7743E+01  4.6991E+00
             1.5051E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1917.70308930787        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1071             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1807E+00  8.9829E-01  1.4213E+00  1.2015E+00  1.2039E+00  2.1697E+00  4.3273E+00  1.1820E+00  1.2157E+00  1.3321E+00
             1.0581E+00
 PARAMETER:  2.6614E-01 -7.2640E-03  4.5156E-01  2.8359E-01  2.8554E-01  8.7459E-01  1.5650E+00  2.6721E-01  2.9528E-01  3.8678E-01
             1.5651E-01
 GRADIENT:   1.0599E+03  4.4012E+01  5.7643E+00  2.7605E+02  1.3700E+01  2.2007E+03  1.7017E+03  7.0847E-01  3.8398E+01  4.6838E+00
             1.5120E+00

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1917.70742279407        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  1.1806E+00  9.0414E-01  1.4210E+00  1.2015E+00  1.2041E+00  2.1684E+00  4.3261E+00  1.1821E+00  1.2241E+00  1.3326E+00
             1.0582E+00
 PARAMETER:  2.6613E-01  2.3115E-04  4.5153E-01  2.8358E-01  2.8574E-01  8.7394E-01  1.5650E+00  2.6723E-01  3.0231E-01  3.8723E-01
             1.5654E-01
 GRADIENT:   1.0110E+04  6.6858E-01  1.1940E+04 -1.3979E+01  1.0943E-01 -2.4330E-01  8.5438E+03 -2.0212E+04  1.7785E+04  1.3893E+04
            -5.1635E-02
 NUMSIGDIG:         2.3         0.8         2.3         5.7         3.3         3.2         2.6         2.3         2.3         2.3
                    3.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1202
 NO. OF SIG. DIGITS IN FINAL EST.:  0.8
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.9261E-03  2.7589E-02 -4.9258E-02 -3.9325E-02 -2.5688E-02
 SE:             2.9982E-02  2.3603E-02  1.3593E-02  1.9851E-02  2.0769E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4878E-01  2.4245E-01  2.9040E-04  4.7590E-02  2.1614E-01

 ETASHRINKSD(%)  1.0000E-10  2.0929E+01  5.4462E+01  3.3497E+01  3.0420E+01
 ETASHRINKVR(%)  1.0000E-10  3.7477E+01  7.9263E+01  5.5774E+01  5.1586E+01
 EBVSHRINKSD(%)  9.2436E-02  1.5092E+01  6.1653E+01  3.4673E+01  2.5664E+01
 EBVSHRINKVR(%)  1.8479E-01  2.7906E+01  8.5295E+01  5.7324E+01  4.4741E+01
 RELATIVEINF(%)  9.9727E+01  2.7596E+01  5.2849E+00  1.3047E+01  2.2675E+01
 EPSSHRINKSD(%)  3.3822E+01
 EPSSHRINKVR(%)  5.6204E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1917.7074227940684     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -998.76888958939571     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1917.707       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.18E+00  9.05E-01  1.42E+00  1.20E+00  1.20E+00  2.17E+00  4.33E+00  1.18E+00  1.22E+00  1.33E+00  1.06E+00
 


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
+        1.36E+06
 
 TH 2
+        1.78E+07  1.64E+07
 
 TH 3
+        3.37E+00  5.33E+00  3.28E+05
 
 TH 4
+        1.26E+06  7.68E+06  1.08E+06  1.16E+06
 
 TH 5
+       -5.47E+00  4.32E+06  1.97E+00 -2.02E+06  1.14E+06
 
 TH 6
+        2.26E+05 -6.58E-04  6.97E+00  2.94E-01 -8.17E-02  4.24E+01
 
 TH 7
+        2.49E-01  4.70E+00  3.54E+01 -5.09E+01  1.02E+05  6.23E-01  7.98E+03
 
 TH 8
+       -1.36E+06 -4.72E+06 -5.35E+03  2.20E+06  1.24E+06 -2.26E+05  1.10E+05  1.36E+06
 
 TH 9
+       -1.16E+06 -1.18E+01 -6.66E+02 -2.95E+06  1.48E+02  1.17E+01 -3.02E+02  1.15E+06  9.81E+05
 
 TH10
+       -8.29E+05 -2.88E+06  4.07E+05 -2.11E+06 -7.57E+05  8.14E+00 -1.16E+02  8.22E+05  1.35E+01  5.05E+05
 
 TH11
+        2.58E+06  8.98E+06 -1.11E+03 -6.58E+06 -5.07E+00  3.73E-01  2.11E+05 -2.58E+06 -1.90E+03 -1.36E+03  4.90E+06
 
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
 #CPUT: Total CPU Time in Seconds,       32.047
Stop Time:
Thu Sep 30 02:45:29 CDT 2021
