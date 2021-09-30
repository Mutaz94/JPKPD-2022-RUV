Wed Sep 29 22:49:01 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat79.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1548.92739942235        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0344E+02  2.4996E+01 -7.8022E+00  7.1777E+01  1.4624E+02  5.9165E+01 -4.5746E+01 -1.2052E+01 -3.5270E+01 -1.5572E+01
            -1.0606E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1819.67899454987        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1061E+00  1.0101E+00  1.0371E+00  1.0502E+00  9.9257E-01  8.8359E-01  1.0725E+00  9.2591E-01  1.0471E+00  7.8881E-01
             2.4531E+00
 PARAMETER:  2.0087E-01  1.1003E-01  1.3646E-01  1.4903E-01  9.2540E-02 -2.3762E-02  1.6997E-01  2.3021E-02  1.4599E-01 -1.3723E-01
             9.9736E-01
 GRADIENT:   3.2098E+02 -1.2359E+01 -3.1550E+01  3.3777E+01  3.9320E+01 -2.2037E+01  6.2050E+00  6.5229E+00  1.6838E+01  2.0742E+01
             1.2587E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1828.85790915221        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      177
 NPARAMETR:  1.1018E+00  7.6341E-01  4.0975E-01  1.1571E+00  4.6886E-01  9.4432E-01  1.6821E+00  5.2381E-01  9.0271E-01  2.3426E-01
             2.2078E+00
 PARAMETER:  1.9698E-01 -1.6996E-01 -7.9220E-01  2.4595E-01 -6.5744E-01  4.2711E-02  6.2006E-01 -5.4662E-01 -2.3503E-03 -1.3513E+00
             8.9199E-01
 GRADIENT:   1.4186E+02  9.0007E+01 -3.3747E+01  1.2664E+02 -1.2473E+01 -1.8038E+00  1.7860E+01  8.4132E+00 -1.2190E+01 -2.0997E+00
             9.6983E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1851.40542041335        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      355
 NPARAMETR:  1.0352E+00  6.8676E-01  5.0380E-01  1.1486E+00  5.3669E-01  8.8067E-01  1.8281E+00  3.5249E-01  8.1366E-01  3.5669E-01
             1.9365E+00
 PARAMETER:  1.3463E-01 -2.7577E-01 -5.8558E-01  2.3858E-01 -5.2234E-01 -2.7074E-02  7.0329E-01 -9.4272E-01 -1.0621E-01 -9.3088E-01
             7.6088E-01
 GRADIENT:   8.8402E+00  2.5919E+01 -5.6430E+01  5.1045E+01  8.8687E+01 -2.0011E+01  1.6966E+01  2.6557E+00 -3.4616E+01 -1.0769E+01
             8.3026E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1863.39011293389        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  1.0343E+00  4.9581E-01  6.0626E-01  1.2687E+00  5.3672E-01  8.9794E-01  1.9802E+00  2.2572E-01  8.8994E-01  5.6540E-01
             1.8922E+00
 PARAMETER:  1.3376E-01 -6.0156E-01 -4.0044E-01  3.3801E-01 -5.2228E-01 -7.6468E-03  7.8322E-01 -1.3884E+00 -1.6597E-02 -4.7022E-01
             7.3775E-01
 GRADIENT:   1.9724E+01  1.3756E+01 -3.5761E+00  4.3376E+01  6.4532E+00 -1.0629E+01 -2.4290E+00  5.9096E-01 -9.5020E+00  8.6707E-01
             3.9538E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1865.73289159172        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  1.0253E+00  3.9565E-01  5.5157E-01  1.2730E+00  4.8335E-01  9.2288E-01  2.2566E+00  1.4658E-01  9.0550E-01  5.3384E-01
             1.8846E+00
 PARAMETER:  1.2495E-01 -8.2723E-01 -4.9499E-01  3.4138E-01 -6.2702E-01  1.9748E-02  9.1385E-01 -1.8202E+00  7.2730E-04 -5.2767E-01
             7.3370E-01
 GRADIENT:   1.7386E+00  1.3326E+00 -1.6516E+00  2.1653E+00  5.4799E-01  1.0499E+00 -8.1656E-01  1.3567E-01  1.9480E-01  1.1914E+00
             2.1352E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1865.83670887942        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  1.0236E+00  3.5390E-01  5.6422E-01  1.2932E+00  4.8222E-01  9.1869E-01  2.4598E+00  5.8478E-02  8.9505E-01  5.4545E-01
             1.8833E+00
 PARAMETER:  1.2328E-01 -9.3874E-01 -4.7231E-01  3.5710E-01 -6.2935E-01  1.5191E-02  1.0001E+00 -2.7391E+00 -1.0881E-02 -5.0615E-01
             7.3305E-01
 GRADIENT:   1.1305E+00  3.1174E-01  2.5606E+00 -2.4806E+00 -2.9169E+00  7.8986E-02  5.2462E-01  4.3524E-03  1.1604E-01 -3.5946E-01
            -3.0565E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1865.84284818214        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  1.0232E+00  3.5577E-01  5.6379E-01  1.2934E+00  4.8278E-01  9.1866E-01  2.4395E+00  1.8533E-02  8.9568E-01  5.4891E-01
             1.8830E+00
 PARAMETER:  1.2295E-01 -9.3348E-01 -4.7308E-01  3.5726E-01 -6.2820E-01  1.5161E-02  9.9180E-01 -3.8882E+00 -1.0174E-02 -4.9983E-01
             7.3286E-01
 GRADIENT:   6.5944E-02 -9.6207E-02  3.3189E-02 -1.4719E-01 -5.6637E-02  7.0617E-03 -7.4141E-02  6.9226E-04  3.7606E-02  2.4297E-02
            -2.0724E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1865.84392149098        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  1.0245E+00  3.5840E-01  5.6381E-01  1.2884E+00  4.8294E-01  9.1895E-01  2.4507E+00  1.0000E-02  8.9639E-01  5.4878E-01
             1.8825E+00
 PARAMETER:  1.2297E-01 -9.2881E-01 -4.7315E-01  3.5676E-01 -6.2751E-01  1.5175E-02  9.8978E-01 -4.6457E+00 -1.0151E-02 -5.0090E-01
             7.3319E-01
 GRADIENT:  -5.9501E-01 -7.5525E-02 -3.5274E-02  1.4394E+00  1.9050E-01 -2.1297E-02 -1.5567E-01  0.0000E+00 -4.0360E-02 -1.2629E-02
             1.0198E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1163
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1123E-03  3.8387E-02 -3.7839E-04 -2.3960E-02  1.3266E-02
 SE:             2.9501E-02  1.9300E-02  2.2417E-04  2.6586E-02  1.7390E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4292E-01  4.6708E-02  9.1420E-02  3.6746E-01  4.4556E-01

 ETASHRINKSD(%)  1.1679E+00  3.5342E+01  9.9249E+01  1.0935E+01  4.1740E+01
 ETASHRINKVR(%)  2.3223E+00  5.8194E+01  9.9994E+01  2.0673E+01  6.6058E+01
 EBVSHRINKSD(%)  1.4018E+00  4.0841E+01  9.9190E+01  9.7512E+00  3.7958E+01
 EBVSHRINKVR(%)  2.7840E+00  6.5003E+01  9.9993E+01  1.8552E+01  6.1508E+01
 RELATIVEINF(%)  9.6486E+01  9.5507E+00  4.6812E-04  3.3103E+01  2.5974E+00
 EPSSHRINKSD(%)  2.8917E+01
 EPSSHRINKVR(%)  4.9472E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1865.8439214909786     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -946.90538828630588     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.20
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.35
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1865.844       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  3.57E-01  5.64E-01  1.29E+00  4.83E-01  9.19E-01  2.43E+00  1.00E-02  8.96E-01  5.48E-01  1.88E+00
 


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
+        1.22E+03
 
 TH 2
+       -3.65E+01  5.79E+02
 
 TH 3
+        1.30E+01  5.47E+02  2.47E+03
 
 TH 4
+       -1.12E+01  2.33E+02 -3.57E+02  6.83E+02
 
 TH 5
+        1.24E+01 -1.13E+03 -3.69E+03  2.29E+02  6.24E+03
 
 TH 6
+        3.09E+00 -6.60E+00  4.06E+00 -6.59E+00 -5.34E-01  2.22E+02
 
 TH 7
+        1.93E+00  4.40E+01 -5.42E+00 -4.61E+00 -4.23E+00  8.36E-01  1.06E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.34E+00 -2.37E+01 -1.16E+00 -3.81E+00  3.61E+01 -1.46E+00  4.81E+00  0.00E+00  1.72E+02
 
 TH10
+       -1.01E+00  1.48E+01 -1.44E+02 -2.11E+01  6.41E+01  5.28E-01  3.36E+00  0.00E+00 -8.31E+00  1.33E+02
 
 TH11
+       -1.48E+01 -2.12E+00 -4.11E+01 -1.40E+01  4.81E+01  2.94E+00  1.28E+00  0.00E+00  8.85E+00  2.60E+01  1.26E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       27.630
Stop Time:
Wed Sep 29 22:49:30 CDT 2021
