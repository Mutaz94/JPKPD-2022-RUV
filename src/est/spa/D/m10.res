Wed Sep 29 19:40:02 CDT 2021
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
$DATA ../../../../data/spa/D/dat10.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m10.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1161.87373039605        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2385E+02 -8.5607E+01 -4.2759E+01 -6.7339E+01  8.7358E+01 -2.5633E+02 -2.7177E+02 -1.1246E+01 -2.8661E+02 -3.3212E+01
            -5.8662E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1406.61622829927        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  1.1075E+00  1.4073E+00  1.1460E+00  8.7899E-01  1.2001E+00  2.3288E+00  4.0990E+00  1.0363E+00  1.7719E+00  9.2012E-01
             1.0587E+00
 PARAMETER:  2.0207E-01  4.4164E-01  2.3624E-01 -2.8984E-02  2.8239E-01  9.4537E-01  1.5107E+00  1.3569E-01  6.7208E-01  1.6751E-02
             1.5703E-01
 GRADIENT:   3.0369E+01  2.0165E+01 -2.3748E+01 -1.7168E+01  5.5006E+01  1.2431E+02  8.8655E+01  2.1457E-01  2.7548E+01 -6.2096E+00
            -1.4350E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1414.95208737957        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  1.1068E+00  7.3896E-01  2.1215E+00  1.5073E+00  1.0923E+00  2.0663E+00  5.9509E+00  1.3659E+00  1.7125E+00  9.3349E-01
             1.0696E+00
 PARAMETER:  2.0151E-01 -2.0251E-01  8.5211E-01  5.1035E-01  1.8830E-01  8.2577E-01  1.8835E+00  4.1183E-01  6.3793E-01  3.1177E-02
             1.6733E-01
 GRADIENT:   3.9901E+01  3.2860E+01  2.8649E+01  2.8144E+01 -4.1456E+01  8.4823E+01  3.4145E+01 -1.4735E+01  5.0307E+01 -9.1525E+00
            -1.8746E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1438.81235362161        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  1.0484E+00  5.4160E-01  2.8155E+00  1.5192E+00  1.2332E+00  1.8834E+00  5.3809E+00  2.4056E+00  1.6431E+00  1.2054E+00
             1.0953E+00
 PARAMETER:  1.4723E-01 -5.1322E-01  1.1351E+00  5.1819E-01  3.0960E-01  7.3305E-01  1.7828E+00  9.7782E-01  5.9656E-01  2.8682E-01
             1.9102E-01
 GRADIENT:   9.6101E+00  1.7776E+01 -8.4626E-01  1.0391E+01 -8.9072E+00  5.8375E+01  2.3549E+01  1.0212E+00  4.1127E+01  1.7240E+01
             2.0631E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1453.50435732699        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:      678
 NPARAMETR:  1.0253E+00  4.1853E-01  2.7725E+00  1.4654E+00  1.2396E+00  1.6048E+00  4.6323E+00  2.4063E+00  1.2332E+00  1.1578E+00
             1.0428E+00
 PARAMETER:  1.2500E-01 -7.7101E-01  1.1197E+00  4.8214E-01  3.1477E-01  5.7303E-01  1.6331E+00  9.7811E-01  3.0963E-01  2.4652E-01
             1.4190E-01
 GRADIENT:   4.8339E+02  1.1828E+02  8.9532E+00  6.5541E+02  1.4388E+01  6.2288E+02  8.2815E+02  8.6909E+00  4.6664E+01  1.1951E+01
            -1.6495E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1454.24074288253        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      860
 NPARAMETR:  1.0338E+00  4.1101E-01  2.7512E+00  1.4769E+00  1.2304E+00  1.6411E+00  4.7787E+00  2.3860E+00  1.2245E+00  1.0910E+00
             1.0493E+00
 PARAMETER:  1.3326E-01 -7.8914E-01  1.1120E+00  4.8994E-01  3.0737E-01  5.9539E-01  1.6642E+00  9.6962E-01  3.0250E-01  1.8707E-01
             1.4815E-01
 GRADIENT:  -1.0846E-01  1.3271E+00  3.6055E-01  8.0710E-01  5.7138E+00  1.3667E+00 -1.9858E+01  2.2850E+00  1.3479E-01  3.7184E+00
            -1.9215E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1454.57958486312        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  1.0366E+00  3.8624E-01  2.7581E+00  1.4849E+00  1.2161E+00  1.7047E+00  4.8437E+00  2.3684E+00  1.2220E+00  1.0662E+00
             1.0585E+00
 PARAMETER:  1.3591E-01 -8.5129E-01  1.1145E+00  4.9533E-01  2.9561E-01  6.3338E-01  1.6777E+00  9.6222E-01  3.0050E-01  1.6410E-01
             1.5689E-01
 GRADIENT:   2.3201E+00  7.9244E-02  2.4177E+00 -3.6226E+00  4.4432E-01  1.8765E+01 -2.2804E+01  1.8099E+00 -4.3754E-01  2.6825E+00
             1.4160E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1454.75221795890        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1194
 NPARAMETR:  1.0354E+00  3.8353E-01  2.7202E+00  1.4761E+00  1.2136E+00  1.6957E+00  4.7652E+00  2.3621E+00  1.2273E+00  1.0476E+00
             1.0559E+00
 PARAMETER:  1.3484E-01 -8.5833E-01  1.1007E+00  4.8939E-01  2.9362E-01  6.2807E-01  1.6613E+00  9.5955E-01  3.0480E-01  1.4647E-01
             1.5443E-01
 GRADIENT:   5.1756E+02  1.1977E+02  9.6136E+00  6.5494E+02  1.5927E+01  6.8307E+02  8.5143E+02  9.0044E+00  4.4547E+01  1.5952E+00
             1.1730E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1454.78237465148        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1382
 NPARAMETR:  1.0367E+00  3.8809E-01  2.7178E+00  1.4743E+00  1.2111E+00  1.7009E+00  4.7551E+00  2.3706E+00  1.2257E+00  1.0393E+00
             1.0561E+00
 PARAMETER:  1.3606E-01 -8.4652E-01  1.0998E+00  4.8815E-01  2.9150E-01  6.3116E-01  1.6592E+00  9.6314E-01  3.0351E-01  1.3853E-01
             1.5462E-01
 GRADIENT:   2.4324E+00 -1.5702E+00  1.5022E+00 -7.7334E+00  1.7546E+00  1.7784E+01 -2.4200E+01  3.5763E+00  2.6643E-01  2.1370E-01
            -2.4724E-01

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1454.78237465148        NO. OF FUNC. EVALS.:  27
 CUMULATIVE NO. OF FUNC. EVALS.:     1409
 NPARAMETR:  1.0364E+00  3.8809E-01  2.7250E+00  1.4735E+00  1.2102E+00  1.7022E+00  4.7551E+00  2.3761E+00  1.2248E+00  1.0378E+00
             1.0564E+00
 PARAMETER:  1.3606E-01 -8.4652E-01  1.0998E+00  4.8815E-01  2.9150E-01  6.3116E-01  1.6592E+00  9.6314E-01  3.0351E-01  1.3853E-01
             1.5462E-01
 GRADIENT:   1.8857E+04 -1.4898E+00 -2.3223E+03  9.5898E-01  4.4024E+03 -6.1551E-01  2.6585E+00 -2.7031E+03  8.4519E+03  1.7125E-01
            -2.4632E-01
 NUMSIGDIG:         2.3         5.9         2.3         2.7         2.3         2.6         5.4         2.3         2.3         1.7
                    2.4

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1409
 NO. OF SIG. DIGITS IN FINAL EST.:  1.7
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2054E-03  4.9866E-02 -6.0521E-02 -5.2115E-02 -4.2642E-02
 SE:             3.0022E-02  2.0352E-02  1.7732E-02  2.1180E-02  1.8273E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6797E-01  1.4281E-02  6.4249E-04  1.3873E-02  1.9613E-02

 ETASHRINKSD(%)  1.0000E-10  3.1817E+01  4.0594E+01  2.9043E+01  3.8785E+01
 ETASHRINKVR(%)  1.0000E-10  5.3511E+01  6.4710E+01  4.9651E+01  6.2527E+01
 EBVSHRINKSD(%)  1.8374E-01  3.4088E+01  4.6276E+01  2.2330E+01  3.2544E+01
 EBVSHRINKVR(%)  3.6714E-01  5.6556E+01  7.1138E+01  3.9673E+01  5.4497E+01
 RELATIVEINF(%)  9.9543E+01  1.5111E+01  1.2567E+01  2.1077E+01  2.0921E+01
 EPSSHRINKSD(%)  4.5644E+01
 EPSSHRINKVR(%)  7.0454E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1454.7823746514841     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -719.63154808774595     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.42
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1454.782       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  3.88E-01  2.72E+00  1.47E+00  1.21E+00  1.70E+00  4.76E+00  2.37E+00  1.23E+00  1.04E+00  1.06E+00
 


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
+        3.22E+06
 
 TH 2
+       -1.95E+02  5.94E+05
 
 TH 3
+       -1.52E+05  6.56E+04  7.13E+03
 
 TH 4
+        6.32E+05  6.89E+01  9.53E+00  1.24E+05
 
 TH 5
+        2.97E+02 -1.21E+02  6.10E+04  8.17E+01  5.15E+05
 
 TH 6
+        1.06E+02 -3.21E-01 -4.79E+00  2.53E-01  4.18E+01  6.87E+01
 
 TH 7
+       -5.77E+04  1.49E+01  2.73E+03 -1.13E+04 -2.30E+04 -9.34E-02  2.06E+03
 
 TH 8
+        1.96E+05  8.42E+04 -8.04E+01  3.84E+04  7.83E+04 -6.58E+00  3.51E+03  1.27E+04
 
 TH 9
+        2.26E+01 -5.25E+05  4.39E+02 -1.80E+02 -4.88E+05  4.02E+01 -2.19E+04  2.95E+01  4.63E+05
 
 TH10
+        2.10E+02  1.90E+02 -8.90E+00 -6.19E+05 -1.26E+06  6.14E-02  5.65E+04 -1.21E+01  8.07E+01  5.50E+01
 
 TH11
+       -1.78E+03 -1.68E+00  8.26E+01  5.45E+05 -7.15E+02  1.07E+00  3.23E-01  1.10E+02 -6.68E+02  1.57E+01  2.40E+06
 
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
 #CPUT: Total CPU Time in Seconds,       32.638
Stop Time:
Wed Sep 29 19:40:37 CDT 2021
