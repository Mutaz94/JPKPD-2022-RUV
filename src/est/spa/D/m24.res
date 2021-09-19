Sat Sep 18 15:13:08 CDT 2021
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
$DATA ../../../../data/spa/D/dat24.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m24.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   4754.59907136209        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.1101E+02 -8.7087E+01  9.0023E+00 -2.0617E+02  2.0026E+02 -1.2310E+03 -4.7325E+02 -8.0068E+01 -6.5197E+02 -4.3665E+02
            -1.0051E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -704.955513428403        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.6336E+00  1.2740E+00  8.6168E-01  1.9496E+00  1.3983E+00  2.9269E+00  2.9588E+00  1.0948E+00  2.4348E+00  1.3970E+00
             1.1068E+01
 PARAMETER:  5.9077E-01  3.4217E-01 -4.8877E-02  7.6763E-01  4.3524E-01  1.1739E+00  1.1848E+00  1.9058E-01  9.8987E-01  4.3432E-01
             2.5041E+00
 GRADIENT:   4.8225E+01  9.9317E+00 -4.7787E+01  3.6546E+01 -5.3316E+00  7.7705E+01  1.9089E+01  3.2899E+00  3.9423E+01  6.0435E+00
             2.0945E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -739.632104582133        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.5412E+00  7.1944E-01  3.3416E+00  3.0792E+00  4.5795E+00  2.4731E+00  1.1875E+01  1.5847E+00  3.0565E+00  5.4090E+00
             8.8480E+00
 PARAMETER:  5.3258E-01 -2.2928E-01  1.3065E+00  1.2247E+00  1.6216E+00  1.0055E+00  2.5744E+00  5.6042E-01  1.2173E+00  1.7881E+00
             2.2802E+00
 GRADIENT:   7.2364E+01  1.5881E+01  4.8836E+00  6.5402E+01  2.5729E+00  1.1918E+01  1.8393E+01 -4.6736E+00  6.5872E+01  2.1431E-02
             1.8259E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -828.965166162151        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0633E+00  4.5785E-01  3.7577E+00  1.7337E+00  2.9065E+00  2.8621E+00  5.0212E+00  4.9044E+00  1.6466E+00  3.2030E+00
             7.3901E+00
 PARAMETER:  1.6138E-01 -6.8120E-01  1.4238E+00  6.5027E-01  1.1669E+00  1.1516E+00  1.7137E+00  1.6901E+00  5.9868E-01  1.2641E+00
             2.1001E+00
 GRADIENT:  -2.5653E+01 -1.4401E+00  7.4033E+00 -5.6218E+01 -1.5904E+01  8.8338E+01  1.0591E+01 -1.2692E+00  3.2292E-02 -1.6177E+00
             1.6289E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -887.714571773588        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.1141E+00  4.4598E-01  2.5370E+00  1.6622E+00  6.0087E+00  1.9656E+00  4.8992E+00  4.3472E+00  1.4181E+00  5.3187E+00
             5.7046E+00
 PARAMETER:  2.0801E-01 -7.0748E-01  1.0310E+00  6.0814E-01  1.8932E+00  7.7578E-01  1.6891E+00  1.5695E+00  4.4933E-01  1.7712E+00
             1.8413E+00
 GRADIENT:   1.4725E+01  9.0440E+00  9.7321E+00 -1.1287E+01 -3.4194E+00  4.3501E+00 -4.1751E+00  5.9256E+00 -7.1464E+00 -2.7023E+00
            -1.4038E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -891.266578920700        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      401
 NPARAMETR:  1.0940E+00  4.1165E-01  2.4075E+00  1.7032E+00  5.8837E+00  1.9137E+00  5.1487E+00  4.2570E+00  1.3947E+00  5.3648E+00
             5.7340E+00
 PARAMETER:  1.8986E-01 -7.8759E-01  9.7858E-01  6.3250E-01  1.8722E+00  7.4902E-01  1.7387E+00  1.5486E+00  4.3270E-01  1.7799E+00
             1.8464E+00
 GRADIENT:   3.3744E+00  9.8622E+00  8.1233E+00  4.8919E+00 -4.6624E+00 -7.4638E+00 -8.3663E+00  9.7315E-01 -4.1256E+00 -1.5729E+00
            -1.1649E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -893.581935414591        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      502
 NPARAMETR:  1.0770E+00  4.1198E-01  2.4075E+00  1.6941E+00  2.0846E+04  1.9191E+00  5.2262E+00  4.2602E+00  1.3932E+00  5.4969E+00
             5.7390E+00
 PARAMETER:  1.7420E-01 -7.8679E-01  9.7858E-01  6.2715E-01  1.0045E+01  7.5185E-01  1.7537E+00  1.5493E+00  4.3163E-01  1.8042E+00
             1.8473E+00
 GRADIENT:   4.5486E+00  1.1381E+01  1.4175E+01  1.2711E+01 -2.2351E-04 -1.6533E+00 -1.6347E+00 -5.6668E+00 -5.2181E+00 -4.6378E-08
            -6.1944E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -893.855439800599        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  1.0601E+00  4.1193E-01  2.4001E+00  1.6711E+00  2.1222E+12  1.9304E+00  5.4269E+00  4.2916E+00  1.3907E+00  5.7721E+00
             5.7847E+00
 PARAMETER:  1.5839E-01 -7.8689E-01  9.7550E-01  6.1346E-01  2.8483E+01  7.5775E-01  1.7914E+00  1.5567E+00  4.2983E-01  1.8530E+00
             1.8552E+00
 GRADIENT:  -4.9012E+00  9.5300E+00  1.2906E+01  5.3569E+00  0.0000E+00  1.0949E+00 -1.3015E+00 -3.2686E+00 -4.0489E-01  0.0000E+00
             2.6784E+00

0ITERATION NO.:   38    OBJECTIVE VALUE:  -893.856168873560        NO. OF FUNC. EVALS.:  88
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  1.0602E+00  4.1193E-01  2.4001E+00  1.6710E+00  2.2191E+12  1.9305E+00  5.4274E+00  4.2917E+00  1.3907E+00  5.7728E+00
             5.7848E+00
 PARAMETER:  1.5842E-01 -7.8689E-01  9.7550E-01  6.1342E-01  2.8528E+01  7.5777E-01  1.7915E+00  1.5567E+00  4.2983E-01  1.8531E+00
             1.8552E+00
 GRADIENT:  -1.4356E+04  2.8970E+03  2.3421E+03 -7.4095E+03  1.2661E-05  5.9907E+03  2.5329E+03 -1.4650E+03 -5.2892E+03  4.6602E-05
            -2.4512E+03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      667
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.1447E-03 -7.8017E-03 -5.1193E-02 -5.4273E-02 -2.0648E-15
 SE:             2.8563E-02  1.5808E-02  1.5464E-02  1.9413E-02  5.5834E-15
 N:                     100         100         100         100         100

 P VAL.:         8.5706E-01  6.2165E-01  9.3160E-04  5.1791E-03  7.1153E-01

 ETASHRINKSD(%)  4.3091E+00  4.7040E+01  4.8193E+01  3.4963E+01  1.0000E+02
 ETASHRINKVR(%)  8.4325E+00  7.1952E+01  7.3161E+01  5.7702E+01  1.0000E+02
 EBVSHRINKSD(%)  3.2840E+00  5.1836E+01  6.6330E+01  2.7333E+01  1.0000E+02
 EBVSHRINKVR(%)  6.4601E+00  7.6802E+01  8.8663E+01  4.7196E+01  1.0000E+02
 RELATIVEINF(%)  8.4913E+01  9.4406E+00  6.2627E+00  1.4868E+01  0.0000E+00
 EPSSHRINKSD(%)  2.1310E+01
 EPSSHRINKVR(%)  3.8080E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -893.85616887355957     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -158.70534230982139     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.90
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -893.856       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  4.12E-01  2.40E+00  1.67E+00  2.22E+12  1.93E+00  5.43E+00  4.29E+00  1.39E+00  5.77E+00  5.78E+00
 


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
+        4.03E+07
 
 TH 2
+       -2.60E+03  1.08E+07
 
 TH 3
+       -3.61E+02 -8.24E+02  2.07E+05
 
 TH 4
+        8.06E+02 -6.49E+02 -9.86E+01  1.08E+06
 
 TH 5
+       -7.09E-15 -7.74E-15 -2.41E-16 -5.97E-16  7.07E-30
 
 TH 6
+       -5.71E+02 -4.76E+03 -6.59E+02  1.50E+03 -2.66E-16  5.30E+05
 
 TH 7
+       -8.64E+01 -1.90E+02 -2.79E+01 -2.17E+01  6.97E-17 -1.59E+02  1.20E+04
 
 TH 8
+        1.25E+02  5.25E+05  3.70E+01  3.14E+01 -1.26E-17  2.32E+02  2.33E+01  2.54E+04
 
 TH 9
+        1.39E+03 -4.47E+02 -5.99E+01  1.15E+02 -1.31E-15  2.57E+03 -1.04E+01  1.98E+01  3.18E+06
 
 TH10
+       -5.01E-02 -1.80E-02  5.32E-04 -8.51E-03  2.45E-17  4.50E-03  3.82E-04 -1.39E-03  1.09E-02  5.44E-04
 
 TH11
+        7.22E+01 -5.41E+02 -7.37E+01  1.65E+02 -6.08E-17  1.45E+02 -1.77E+01  2.55E+01  1.75E+01 -7.24E-04  9.89E+03
 
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
 #CPUT: Total CPU Time in Seconds,       20.119
Stop Time:
Sat Sep 18 15:13:29 CDT 2021
