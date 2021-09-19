Sat Sep 18 02:39:26 CDT 2021
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
$DATA ../../../../data/int/S1/dat49.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3598.37504127227        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5214E+02  1.6911E+01  2.2655E+01 -9.0743E+00  7.3816E+01 -1.5060E+01 -6.7278E+01 -4.3915E+01 -1.9623E+01 -3.0859E+01
            -3.1022E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3634.47008215608        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  9.6260E-01  8.9510E-01  1.0153E+00  1.0528E+00  9.4541E-01  9.7291E-01  1.4201E+00  1.2217E+00  9.5747E-01  9.3426E-01
             1.1524E+00
 PARAMETER:  6.1877E-02 -1.0824E-02  1.1515E-01  1.5143E-01  4.3862E-02  7.2533E-02  4.5072E-01  3.0023E-01  5.6535E-02  3.2002E-02
             2.4184E-01
 GRADIENT:   5.6475E+01 -3.2372E+01  1.5290E+01  2.8581E+01  4.8390E+01 -2.2856E+01 -8.7910E-01 -9.8857E+00 -6.5286E+00 -1.1995E+01
             4.8130E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3634.54134133128        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      216
 NPARAMETR:  9.6269E-01  8.9572E-01  1.0147E+00  1.0523E+00  9.4583E-01  9.7441E-01  1.4203E+00  1.2265E+00  9.5693E-01  9.3216E-01
             1.1523E+00
 PARAMETER:  6.1979E-02 -1.0123E-02  1.1459E-01  1.5097E-01  4.4304E-02  7.4074E-02  4.5090E-01  3.0415E-01  5.5971E-02  2.9748E-02
             2.4179E-01
 GRADIENT:   2.3887E+01 -3.7379E+01  1.3602E+01  1.5196E+01  4.4852E+01 -2.7126E+01 -5.3772E+00 -9.6990E+00 -7.6906E+00 -1.2803E+01
             4.7760E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3635.69925423368        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      402             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6267E-01  9.0006E-01  1.0147E+00  1.0523E+00  9.4279E-01  1.0413E+00  1.4203E+00  1.2264E+00  9.5690E-01  9.3218E-01
             1.1524E+00
 PARAMETER:  6.1954E-02 -5.2973E-03  1.1462E-01  1.5097E-01  4.1084E-02  1.4047E-01  4.5090E-01  3.0407E-01  5.5946E-02  2.9771E-02
             2.4185E-01
 GRADIENT:   5.3835E+01 -2.6948E+01  1.6131E+01  2.8082E+01  4.1901E+01  7.5842E+00 -3.9413E-01 -9.2939E+00 -6.8128E+00 -1.1976E+01
             4.8433E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3635.72171147666        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  9.6261E-01  9.0021E-01  1.0146E+00  1.0523E+00  9.4268E-01  1.0396E+00  1.4203E+00  1.2268E+00  9.5703E-01  9.3252E-01
             1.1524E+00
 PARAMETER:  6.1893E-02 -5.1322E-03  1.1451E-01  1.5093E-01  4.0970E-02  1.3882E-01  4.5090E-01  3.0442E-01  5.6074E-02  3.0133E-02
             2.4185E-01
 GRADIENT:   5.3776E+01 -2.6786E+01  1.6101E+01  2.8013E+01  4.1644E+01  6.8557E+00 -3.6275E-01 -9.2332E+00 -6.7859E+00 -1.1904E+01
             4.8483E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3635.85356524879        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  9.6261E-01  9.0240E-01  1.0146E+00  1.0523E+00  9.4118E-01  1.0395E+00  1.4203E+00  1.2268E+00  9.5703E-01  9.3252E-01
             1.1524E+00
 PARAMETER:  6.1894E-02 -2.7011E-03  1.1451E-01  1.5093E-01  3.9380E-02  1.3871E-01  4.5090E-01  3.0442E-01  5.6075E-02  3.0132E-02
             2.4184E-01
 GRADIENT:   2.1027E+01 -2.9439E+01  1.6170E+01  1.5267E+01  3.4163E+01 -5.7891E-02 -4.5594E+00 -9.4415E+00 -7.8506E+00 -1.2139E+01
             4.7986E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3636.19691258730        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      918             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6261E-01  9.1506E-01  1.0146E+00  1.0523E+00  9.4118E-01  1.0395E+00  1.4203E+00  1.2280E+00  9.5703E-01  9.3252E-01
             1.1524E+00
 PARAMETER:  6.1894E-02  1.1239E-02  1.1451E-01  1.5093E-01  3.9380E-02  1.3876E-01  4.5090E-01  3.0542E-01  5.6075E-02  3.0132E-02
             2.4184E-01
 GRADIENT:   5.3748E+01 -1.3148E+01  1.7136E+01  3.2700E+01  3.1852E+01  6.8072E+00  9.5360E-01 -8.9667E+00 -7.0823E+00 -1.2054E+01
             4.7634E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3636.41334770665        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1055
 NPARAMETR:  9.6261E-01  9.3676E-01  1.0146E+00  1.0523E+00  9.4118E-01  1.0398E+00  1.4203E+00  1.2280E+00  9.5703E-01  9.3252E-01
             1.1524E+00
 PARAMETER:  6.1894E-02  3.4667E-02  1.1451E-01  1.5093E-01  3.9380E-02  1.3901E-01  4.5090E-01  3.0543E-01  5.6075E-02  3.0132E-02
             2.4184E-01
 GRADIENT:   2.0937E+01  1.5420E-03  1.6911E+01  2.7742E+01  1.7122E+01 -1.7503E-03 -2.1102E+00 -9.1179E+00 -8.5306E+00 -1.3029E+01
             4.5710E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3636.42122644017        NO. OF FUNC. EVALS.: 206
 CUMULATIVE NO. OF FUNC. EVALS.:     1261
 NPARAMETR:  9.6259E-01  9.3676E-01  1.0146E+00  1.0523E+00  9.4116E-01  1.0398E+00  1.4203E+00  1.2286E+00  9.5703E-01  9.3264E-01
             1.1524E+00
 PARAMETER:  6.1874E-02  3.4667E-02  1.1445E-01  1.5093E-01  3.9357E-02  1.3901E-01  4.5090E-01  3.0590E-01  5.6075E-02  3.0267E-02
             2.4184E-01
 GRADIENT:  -2.3574E+06  1.0135E-02  1.0299E+06 -1.5619E+06  1.1788E+06 -6.0866E-03 -2.6142E+05  7.7054E+05 -1.1788E+06  2.3574E+06
             4.8728E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1261
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.2494E-03 -1.8694E-02 -3.7549E-02  8.5878E-04 -4.1771E-02
 SE:             2.9855E-02  2.5041E-02  2.0744E-02  2.7902E-02  2.3392E-02
 N:                     100         100         100         100         100

 P VAL.:         7.5670E-01  4.5535E-01  7.0274E-02  9.7545E-01  7.4149E-02

 ETASHRINKSD(%)  1.0000E-10  1.6108E+01  3.0506E+01  6.5233E+00  2.1634E+01
 ETASHRINKVR(%)  1.0000E-10  2.9622E+01  5.1705E+01  1.2621E+01  3.8588E+01
 EBVSHRINKSD(%)  3.3406E-01  1.5693E+01  3.4888E+01  1.0349E+01  2.5413E+01
 EBVSHRINKVR(%)  6.6700E-01  2.8923E+01  5.7604E+01  1.9626E+01  4.4368E+01
 RELATIVEINF(%)  9.9331E+01  4.1904E+01  2.8133E+01  5.7622E+01  2.4310E+01
 EPSSHRINKSD(%)  2.3395E+01
 EPSSHRINKVR(%)  4.1317E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3636.4212264401699     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1982.3318666717591     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.89
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3636.421       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.63E-01  9.37E-01  1.01E+00  1.05E+00  9.41E-01  1.04E+00  1.42E+00  1.23E+00  9.57E-01  9.33E-01  1.15E+00
 


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
+        6.36E+09
 
 TH 2
+       -1.53E+04  4.49E+02
 
 TH 3
+        4.50E+03  1.27E+04  4.37E+09
 
 TH 4
+       -3.29E+03 -9.09E+03 -9.04E+04  2.34E+09
 
 TH 5
+        5.54E+03  1.54E+04  5.39E+09 -4.77E+04  6.65E+09
 
 TH 6
+       -8.38E+03  1.69E+00  6.93E+03 -5.07E+03  8.55E+03  1.84E+02
 
 TH 7
+       -8.15E+02 -2.28E+03 -2.24E+04  1.46E+04 -1.19E+04 -1.26E+03  1.44E+08
 
 TH 8
+        4.87E+05 -1.67E+09 -4.04E+05  2.95E+05 -4.98E+05  2.14E+03  7.32E+04  4.17E+08
 
 TH 9
+       -5.46E+03 -1.54E+04 -1.50E+05  3.88E+09 -6.54E+09 -8.41E+03  1.01E+04  4.90E+05  6.43E+09
 
 TH10
+        5.60E+03  1.57E+04  5.44E+09  4.70E+04  6.71E+09  8.63E+03  1.17E+04 -5.03E+05  7.80E+04  6.78E+09
 
 TH11
+        1.87E+03  5.24E+03  5.14E+04 -3.36E+04  2.73E+04  2.89E+03 -5.79E+07 -1.68E+05 -2.33E+04 -2.68E+04  7.58E+08
 
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
 #CPUT: Total CPU Time in Seconds,       55.138
Stop Time:
Sat Sep 18 02:40:23 CDT 2021
