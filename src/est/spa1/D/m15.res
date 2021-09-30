Thu Sep 30 02:41:31 CDT 2021
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
$DATA ../../../../data/spa1/D/dat15.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1347.94263739407        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2133E+02 -1.4872E+02 -6.1746E+01 -2.5811E+02  2.4403E+02 -3.6714E+02 -2.8060E+02 -6.5904E+01 -5.5582E+02 -1.7849E+02
            -1.3930E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1803.92256279331        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  1.0692E+00  1.1033E+00  1.2330E+00  1.1473E+00  8.1510E-01  1.6351E+00  1.9573E+00  1.3285E+00  2.2281E+00  1.5586E+00
             1.0734E+00
 PARAMETER:  1.6691E-01  1.9826E-01  3.0944E-01  2.3744E-01 -1.0444E-01  5.9173E-01  7.7159E-01  3.8408E-01  9.0115E-01  5.4376E-01
             1.7086E-01
 GRADIENT:  -5.5618E+01  1.6498E+01  3.2825E+01 -1.4505E+01 -1.3583E+02 -3.9826E+01 -4.0644E+01 -3.0240E+01  1.7759E+01 -7.8063E-01
            -2.8443E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1873.64020319783        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      336
 NPARAMETR:  1.0624E+00  1.2977E+00  1.8227E+00  1.1061E+00  1.1786E+00  1.6209E+00  3.4095E+00  2.2199E+00  1.7096E+00  1.3898E+00
             1.0891E+00
 PARAMETER:  1.6057E-01  3.6058E-01  7.0030E-01  2.0082E-01  2.6433E-01  5.8301E-01  1.3266E+00  8.9748E-01  6.3625E-01  4.2916E-01
             1.8534E-01
 GRADIENT:  -6.2203E+01  1.6612E+01  1.2798E+01  1.6772E+01 -4.2539E+01 -4.1624E+01  7.0810E-01 -1.1493E+01  2.1694E+01  2.6329E+01
            -9.2873E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1880.10446465769        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      505             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0691E+00  1.2704E+00  1.8238E+00  1.0989E+00  1.1873E+00  1.8611E+00  3.3784E+00  2.2362E+00  1.6700E+00  1.3471E+00
             1.0896E+00
 PARAMETER:  1.6679E-01  3.3933E-01  7.0092E-01  1.9431E-01  2.7165E-01  7.2115E-01  1.3174E+00  9.0476E-01  6.1281E-01  3.9796E-01
             1.8582E-01
 GRADIENT:   5.7520E+02  1.7362E+02  1.8936E+01  1.6260E+02 -2.2937E+01  8.9961E+02  6.5457E+02 -4.4503E+00  9.8139E+01  2.5470E+01
            -7.4908E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1880.86171031754        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:      678            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0702E+00  1.2502E+00  1.8240E+00  1.0975E+00  1.1876E+00  1.8520E+00  3.2590E+00  2.2318E+00  1.6446E+00  1.3395E+00
             1.1019E+00
 PARAMETER:  1.6780E-01  3.2329E-01  7.0104E-01  1.9304E-01  2.7193E-01  7.1626E-01  1.2814E+00  9.0283E-01  5.9752E-01  3.9230E-01
             1.9701E-01
 GRADIENT:   5.6534E+02  1.5845E+02  1.9434E+01  1.5742E+02 -2.2937E+01  8.7050E+02  6.1689E+02 -5.2897E+00  9.1156E+01  2.4926E+01
             1.4380E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1884.04019454394        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      844
 NPARAMETR:  1.1176E+00  1.2494E+00  1.8263E+00  1.0971E+00  1.2331E+00  1.8189E+00  3.2665E+00  2.2283E+00  1.4984E+00  1.3404E+00
             1.1028E+00
 PARAMETER:  2.1120E-01  3.2270E-01  7.0230E-01  1.9269E-01  3.0951E-01  6.9822E-01  1.2837E+00  9.0122E-01  5.0440E-01  3.9300E-01
             1.9783E-01
 GRADIENT:  -1.2465E+01  9.7322E+00  1.0581E+01  2.1017E+01 -7.1470E+00  2.0512E+01 -1.3360E+01 -1.4412E+01  6.8941E+00  1.8625E+01
            -8.6516E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1885.79182590623        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1025             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1362E+00  1.2296E+00  1.7946E+00  1.0837E+00  1.2431E+00  1.7914E+00  3.2846E+00  2.2963E+00  1.4475E+00  1.2956E+00
             1.1045E+00
 PARAMETER:  2.2770E-01  3.0668E-01  6.8479E-01  1.8034E-01  3.1763E-01  6.8303E-01  1.2892E+00  9.3128E-01  4.6984E-01  3.5896E-01
             1.9937E-01
 GRADIENT:   8.1911E+02  1.4396E+02  1.1497E+01  1.4970E+02  1.7346E+01  7.6295E+02  6.4532E+02 -3.5227E+00  5.2637E+01  1.6084E+01
             1.2568E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1886.14920249957        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1212
 NPARAMETR:  1.1359E+00  1.2281E+00  1.7898E+00  1.0821E+00  1.2425E+00  1.8124E+00  3.2979E+00  2.3018E+00  1.4025E+00  1.2890E+00
             1.1051E+00
 PARAMETER:  2.2741E-01  3.0550E-01  6.8213E-01  1.7889E-01  3.1714E-01  6.9463E-01  1.2933E+00  9.3370E-01  4.3826E-01  3.5390E-01
             1.9995E-01
 GRADIENT:  -8.1475E-01  6.0747E+00  5.8377E+00  1.7142E+01  3.4131E+00  1.7972E+01 -1.3631E+01 -9.8344E+00  2.2999E+00  1.2578E+01
            -6.9060E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1886.37587527352        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1385
 NPARAMETR:  1.1352E+00  1.2259E+00  1.7819E+00  1.0801E+00  1.2412E+00  1.8376E+00  3.2965E+00  2.3083E+00  1.4020E+00  1.2801E+00
             1.1056E+00
 PARAMETER:  2.2677E-01  3.0370E-01  6.7769E-01  1.7702E-01  3.1609E-01  7.0844E-01  1.2929E+00  9.3650E-01  4.3789E-01  3.4691E-01
             2.0041E-01
 GRADIENT:   1.1638E+04 -8.6907E+03  1.9583E+03 -1.4905E+04  4.1809E+03 -7.5066E-01 -1.0246E+03  2.7695E+03 -6.0318E+03  3.8184E+03
            -4.5408E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1385
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.2836E-03  5.7989E-03 -6.4050E-02 -3.3135E-02 -4.5364E-02
 SE:             2.9992E-02  2.5085E-02  1.7253E-02  1.9185E-02  1.8827E-02
 N:                     100         100         100         100         100

 P VAL.:         8.0812E-01  8.1718E-01  2.0532E-04  8.4138E-02  1.5972E-02

 ETASHRINKSD(%)  1.0000E-10  1.5963E+01  4.2201E+01  3.5729E+01  3.6928E+01
 ETASHRINKVR(%)  1.0000E-10  2.9378E+01  6.6593E+01  5.8692E+01  6.0219E+01
 EBVSHRINKSD(%)  1.3680E-01  1.4699E+01  5.1434E+01  3.8660E+01  2.6451E+01
 EBVSHRINKVR(%)  2.7341E-01  2.7237E+01  7.6414E+01  6.2374E+01  4.5906E+01
 RELATIVEINF(%)  9.9640E+01  2.2396E+01  9.8564E+00  9.7723E+00  2.7150E+01
 EPSSHRINKSD(%)  3.5099E+01
 EPSSHRINKVR(%)  5.7879E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1886.3758752735205     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -967.43734206884778     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.22
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1886.376       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.14E+00  1.23E+00  1.78E+00  1.08E+00  1.24E+00  1.84E+00  3.30E+00  2.31E+00  1.40E+00  1.28E+00  1.11E+00
 


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
+        9.97E+05
 
 TH 2
+       -5.30E+01  4.76E+05
 
 TH 3
+        1.03E+03 -6.80E+00  4.55E+04
 
 TH 4
+       -5.15E+02  1.23E+02 -2.87E+05  1.81E+06
 
 TH 5
+        5.14E+01 -4.97E+01 -1.95E+01 -8.03E+00  4.29E+05
 
 TH 6
+        2.11E+01 -1.19E+01  3.82E+00 -2.31E+01  1.11E+01  5.91E+01
 
 TH 7
+       -5.25E+01  4.16E+04  1.28E+04  1.70E+01  3.49E-01 -1.22E+00  7.28E+03
 
 TH 8
+       -3.51E+03 -8.28E+00  2.50E+04  4.73E+03  7.67E+04  2.07E+00  2.12E+02  1.37E+04
 
 TH 9
+       -3.76E+02  2.53E+02 -8.48E+01  5.00E+02 -2.46E+02 -8.27E+04 -2.52E+04 -4.10E+01  1.75E+05
 
 TH10
+        6.29E+01 -4.28E+01  1.46E+01 -8.90E+01  6.33E-01  1.00E+01 -4.02E+00  8.29E+00 -2.42E+05  3.35E+05
 
 TH11
+       -1.06E+03  7.29E+02 -2.26E+02  1.42E+03 -6.98E+02  3.83E-01  6.43E+01 -1.24E+02  4.86E+05 -6.71E+05  3.26E+02
 
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
 #CPUT: Total CPU Time in Seconds,       37.352
Stop Time:
Thu Sep 30 02:42:10 CDT 2021
