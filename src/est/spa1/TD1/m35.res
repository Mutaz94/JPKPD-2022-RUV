Thu Sep 30 01:20:59 CDT 2021
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
$DATA ../../../../data/spa1/TD1/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2147.95370415589        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4800E+02  1.1773E+01 -4.2700E+01  8.7558E+01  7.2736E+01  6.3610E+01  5.1835E-01  1.3772E+01  6.6412E+00 -1.2510E+01
            -2.0779E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2156.70606318058        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0580E+00  9.9171E-01  1.2168E+00  1.0422E+00  1.0532E+00  9.3849E-01  1.0216E+00  7.7828E-01  1.0084E+00  1.1722E+00
             9.9325E-01
 PARAMETER:  1.5636E-01  9.1678E-02  2.9620E-01  1.4136E-01  1.5184E-01  3.6516E-02  1.2136E-01 -1.5067E-01  1.0840E-01  2.5887E-01
             9.3232E-02
 GRADIENT:   4.5973E+00  2.1684E+01  9.8471E-01  2.7153E+01 -5.1913E-01 -1.4917E+01  1.5496E+00  5.2088E-01  2.0225E+00 -3.6145E+00
            -9.8734E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2157.50467814652        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0576E+00  8.8491E-01  1.2653E+00  1.0972E+00  1.0228E+00  9.5743E-01  9.8936E-01  7.0329E-01  9.8179E-01  1.1980E+00
             1.0024E+00
 PARAMETER:  1.5596E-01 -2.2269E-02  3.3528E-01  1.9277E-01  1.2256E-01  5.6496E-02  8.9300E-02 -2.5199E-01  8.1625E-02  2.8066E-01
             1.0237E-01
 GRADIENT:   6.1530E+00  1.1783E+01  8.3069E+00  9.8971E+00 -1.1577E+01 -6.0622E+00 -1.8391E-01 -2.4652E+00  1.3742E+00  7.8705E-01
            -2.7168E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2158.06333822825        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0476E+00  7.1626E-01  1.5232E+00  1.2121E+00  1.0534E+00  9.6401E-01  1.0673E+00  9.7407E-01  9.1861E-01  1.2425E+00
             1.0084E+00
 PARAMETER:  1.4653E-01 -2.3371E-01  5.2084E-01  2.9235E-01  1.5206E-01  6.3347E-02  1.6518E-01  7.3723E-02  1.5107E-02  3.1709E-01
             1.0835E-01
 GRADIENT:  -1.1510E+01  9.1213E+00  4.8103E+00  1.0267E+01 -8.2294E+00 -2.3617E+00  9.3754E-01  3.0663E-01  9.4302E-01 -2.2316E-01
             2.8065E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2158.59226105833        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0516E+00  4.8417E-01  1.6705E+00  1.3618E+00  1.0299E+00  9.7239E-01  1.0511E+00  1.0490E+00  8.5080E-01  1.2460E+00
             1.0037E+00
 PARAMETER:  1.5036E-01 -6.2531E-01  6.1313E-01  4.0880E-01  1.2949E-01  7.2001E-02  1.4988E-01  1.4788E-01 -6.1579E-02  3.1996E-01
             1.0369E-01
 GRADIENT:   4.9319E+00  5.4920E+00  1.0751E+00  1.4683E+01 -1.1560E+00  2.2898E+00 -2.0961E-02 -8.5514E-01 -1.0211E+00 -9.4523E-01
            -1.7861E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2158.67812626180        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  1.0489E+00  3.5608E-01  1.7875E+00  1.4453E+00  1.0233E+00  9.6601E-01  8.8855E-01  1.1479E+00  8.2243E-01  1.2611E+00
             1.0049E+00
 PARAMETER:  1.4778E-01 -9.3259E-01  6.8084E-01  4.6834E-01  1.2299E-01  6.5419E-02 -1.8160E-02  2.3795E-01 -9.5488E-02  3.3202E-01
             1.0485E-01
 GRADIENT:   2.7660E+00  4.3490E+00  4.4476E-01  1.5195E+01 -2.8578E+00  3.5120E-01 -1.1010E-02 -3.9391E-01  2.4510E-03  6.6203E-01
            -8.5324E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2158.76867507897        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1081             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0479E+00  2.5762E-01  1.8746E+00  1.4914E+00  1.0245E+00  9.6427E-01  6.8219E-01  1.2239E+00  7.9864E-01  1.2614E+00
             1.0057E+00
 PARAMETER:  1.4677E-01 -1.2563E+00  7.2838E-01  4.9969E-01  1.2417E-01  6.3620E-02 -2.8244E-01  3.0207E-01 -1.2485E-01  3.3225E-01
             1.0566E-01
 GRADIENT:   6.9231E+02  3.7227E+01  7.8695E+00  9.3539E+02  1.2643E+01  5.0364E+01  4.0547E-01  6.5468E-01  1.3622E+01  3.9468E+00
             1.3427E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2158.82283736854        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1267
 NPARAMETR:  1.0473E+00  2.5950E-01  1.8766E+00  1.4969E+00  1.0207E+00  9.6421E-01  2.5493E-01  1.2257E+00  7.9774E-01  1.2609E+00
             1.0057E+00
 PARAMETER:  1.4617E-01 -1.2490E+00  7.2949E-01  5.0340E-01  1.2054E-01  6.3550E-02 -1.2668E+00  3.0354E-01 -1.2597E-01  3.3185E-01
             1.0565E-01
 GRADIENT:   2.0822E+00  6.8090E-01  1.6480E-01 -1.0434E+01 -1.0578E-01  1.2425E-01  4.1644E-03 -5.3205E-02 -6.1418E-01  1.6568E-01
            -7.9494E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2158.82712480761        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1451
 NPARAMETR:  1.0480E+00  2.5692E-01  1.8763E+00  1.4972E+00  1.0197E+00  9.6420E-01  4.3851E-02  1.2264E+00  7.9887E-01  1.2598E+00
             1.0057E+00
 PARAMETER:  1.4686E-01 -1.2590E+00  7.2932E-01  5.0361E-01  1.1949E-01  6.3547E-02 -3.0270E+00  3.0408E-01 -1.2455E-01  3.3093E-01
             1.0569E-01
 GRADIENT:   3.8589E+00  3.2186E-01  2.9736E-01 -1.3096E+01 -1.9811E-01  1.3928E-01  5.6036E-04 -2.9028E-02  4.1868E-02  1.3083E-01
            -3.0376E-02

0ITERATION NO.:   44    OBJECTIVE VALUE:  -2158.82746770341        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1586
 NPARAMETR:  1.0479E+00  2.5645E-01  1.8756E+00  1.4982E+00  1.0186E+00  9.6417E-01  1.0000E-02  1.2251E+00  7.9836E-01  1.2588E+00
             1.0057E+00
 PARAMETER:  1.4686E-01 -1.2624E+00  7.2844E-01  5.0382E-01  1.1958E-01  6.3541E-02 -6.0851E+00  3.0449E-01 -1.2471E-01  3.3037E-01
             1.0573E-01
 GRADIENT:   3.8259E-02 -3.7846E-02 -7.5673E-02 -7.7265E-01  5.5460E-01  5.0601E-03  0.0000E+00  4.9819E-02  9.0183E-02  2.0976E-02
             1.7707E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1586
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.7333E-04 -1.0682E-04 -3.1948E-02 -4.8062E-03 -3.7538E-02
 SE:             2.9896E-02  4.8128E-05  1.4936E-02  2.9444E-02  2.2890E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8737E-01  2.6446E-02  3.2443E-02  8.7034E-01  1.0102E-01

 ETASHRINKSD(%)  1.0000E-10  9.9839E+01  4.9961E+01  1.3578E+00  2.3317E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  7.4961E+01  2.6971E+00  4.1197E+01
 EBVSHRINKSD(%)  3.3778E-01  9.9854E+01  5.5105E+01  1.6345E+00  1.8565E+01
 EBVSHRINKVR(%)  6.7441E-01  1.0000E+02  7.9844E+01  3.2422E+00  3.3683E+01
 RELATIVEINF(%)  9.8164E+01  1.5977E-05  6.3519E+00  7.8624E+00  1.8825E+01
 EPSSHRINKSD(%)  3.3639E+01
 EPSSHRINKVR(%)  5.5962E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2158.8274677034060     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1239.8889344987333     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.24
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.19
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2158.827       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  2.56E-01  1.87E+00  1.50E+00  1.02E+00  9.64E-01  1.00E-02  1.23E+00  7.99E-01  1.26E+00  1.01E+00
 


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
+        1.08E+03
 
 TH 2
+       -2.42E+01  3.81E+02
 
 TH 3
+       -3.13E+00  1.84E+01  4.64E+01
 
 TH 4
+       -7.20E+00  4.75E+02 -2.46E+01  7.52E+02
 
 TH 5
+        3.03E+00 -1.38E+02 -1.05E+02 -2.68E+01  4.53E+02
 
 TH 6
+        7.94E-02 -3.05E+00 -2.98E-01 -1.38E+00  6.41E-01  2.12E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -1.34E-01 -4.95E-01 -1.67E+01 -3.06E+00  4.34E-01 -6.22E-02  0.00E+00  2.27E+01
 
 TH 9
+        1.55E+00 -1.06E+02  3.11E+00  1.26E+00 -4.07E-01 -8.95E-01  0.00E+00 -2.16E-01  2.94E+02
 
 TH10
+        2.18E+00  3.21E+00 -6.18E+00  1.22E+00 -5.36E+01  2.91E-01  0.00E+00  7.56E+00 -1.53E-01  5.98E+01
 
 TH11
+       -8.27E+00 -1.23E+01 -5.03E+00 -1.23E+01  3.01E+00  2.15E+00  0.00E+00  5.24E+00  7.85E+00  1.01E+01  4.07E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.497
Stop Time:
Thu Sep 30 01:21:33 CDT 2021
