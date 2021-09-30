Wed Sep 29 12:45:16 CDT 2021
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
$DATA ../../../../data/spa/A2/dat34.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m34.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1300.19355738295        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7328E+02 -7.9944E+01 -8.0396E+00 -6.8732E+01  1.0072E+02 -1.3594E+01 -8.8389E+00 -5.6158E+00  1.4257E+01 -5.1774E+01
            -5.5682E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1456.21460208074        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9518E-01  1.1242E+00  1.0121E+00  1.0356E+00  9.4725E-01  1.1347E+00  9.6441E-01  9.7623E-01  7.6699E-01  1.0744E+00
             2.1103E+00
 PARAMETER:  9.5167E-02  2.1708E-01  1.1206E-01  1.3499E-01  4.5811E-02  2.2638E-01  6.3761E-02  7.5944E-02 -1.6528E-01  1.7177E-01
             8.4683E-01
 GRADIENT:   1.3601E+02  3.0652E+01  1.3253E+01  1.7645E+01 -3.5896E+01  3.6518E+01 -3.9275E+00  1.7396E+00 -1.2147E+00  6.1981E+00
             4.7442E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1459.79407499169        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.4800E-01  1.0541E+00  1.0764E+00  1.0669E+00  9.7154E-01  1.1021E+00  1.1757E+00  7.0448E-01  7.2029E-01  1.0755E+00
             2.0101E+00
 PARAMETER:  4.6595E-02  1.5265E-01  1.7361E-01  1.6474E-01  7.1130E-02  1.9721E-01  2.6185E-01 -2.5030E-01 -2.2809E-01  1.7274E-01
             7.9819E-01
 GRADIENT:   6.6421E+01  1.4662E+01  1.3021E+01  5.2664E+00 -1.7012E+01  2.9155E+01  7.0791E+00 -2.4229E-01  1.2593E+00 -1.3155E+00
             2.4475E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1461.01504831983        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      275
 NPARAMETR:  9.5062E-01  1.0238E+00  8.5534E-01  1.0785E+00  8.6404E-01  1.0778E+00  1.1489E+00  4.7615E-01  7.4259E-01  9.9684E-01
             1.9559E+00
 PARAMETER:  4.9364E-02  1.2352E-01 -5.6262E-02  1.7554E-01 -4.6140E-02  1.7494E-01  2.3882E-01 -6.4202E-01 -1.9761E-01  9.6836E-02
             7.7086E-01
 GRADIENT:  -2.3237E+01  1.7302E+00 -2.1762E-02 -4.0528E+00 -1.8957E+00  1.8899E+00  1.3709E+00  3.1893E-01  1.2333E-01  4.8283E-01
             1.1025E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1461.80471137894        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  9.6009E-01  7.7694E-01  9.8344E-01  1.2510E+00  8.2918E-01  1.0787E+00  1.3608E+00  1.8091E-01  7.1255E-01  1.0570E+00
             1.8770E+00
 PARAMETER:  5.9274E-02 -1.5240E-01  8.3297E-02  3.2398E-01 -8.7323E-02  1.7579E-01  4.0807E-01 -1.6097E+00 -2.3890E-01  1.5542E-01
             7.2967E-01
 GRADIENT:   1.0115E+00  1.4142E+01  4.6825E+00  2.5516E+01 -8.3466E+00  2.5738E+00 -1.3466E+00  5.2481E-02 -1.2126E+00 -1.0510E+00
            -1.1230E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1462.75376155512        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      627
 NPARAMETR:  9.6070E-01  5.4215E-01  9.9285E-01  1.3728E+00  7.5823E-01  1.0555E+00  1.7984E+00  1.7131E-02  6.6835E-01  1.0332E+00
             1.9362E+00
 PARAMETER:  5.9908E-02 -5.1222E-01  9.2828E-02  4.1686E-01 -1.7676E-01  1.5405E-01  6.8690E-01 -3.9668E+00 -3.0294E-01  1.3266E-01
             7.6075E-01
 GRADIENT:   6.1905E+00  2.8860E+00  3.8820E+00 -4.4029E+00 -7.1627E+00 -3.8366E+00  9.6884E-01  1.1431E-03  5.3607E-01  5.9622E-01
             6.5169E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1463.56648445089        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      802
 NPARAMETR:  9.5041E-01  2.9536E-01  1.1361E+00  1.5280E+00  7.7006E-01  1.0729E+00  2.5747E+00  1.0000E-02  6.4390E-01  1.1053E+00
             1.8847E+00
 PARAMETER:  4.9133E-02 -1.1196E+00  2.2763E-01  5.2398E-01 -1.6128E-01  1.7039E-01  1.0457E+00 -8.8418E+00 -3.4022E-01  2.0015E-01
             7.3375E-01
 GRADIENT:  -4.0767E+00  2.5335E-01 -3.5360E+00 -1.7990E+00  7.4471E+00  3.1468E+00 -1.4835E-01  0.0000E+00 -6.6984E-01 -2.0370E+00
            -6.6028E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1463.92981598743        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      978
 NPARAMETR:  9.5069E-01  1.5965E-01  1.2951E+00  1.6277E+00  7.9326E-01  1.0513E+00  3.5176E+00  1.0000E-02  6.3321E-01  1.1844E+00
             1.9343E+00
 PARAMETER:  4.9435E-02 -1.7348E+00  3.5858E-01  5.8715E-01 -1.3161E-01  1.5004E-01  1.3578E+00 -1.4262E+01 -3.5695E-01  2.6922E-01
             7.5976E-01
 GRADIENT:   9.3958E-03  1.9681E+00  3.3271E+00  1.5587E+01 -6.8780E+00 -3.2104E+00 -4.5799E-01  0.0000E+00  6.9055E-01  2.0036E+00
             5.0846E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1464.92812032183        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1156
 NPARAMETR:  9.4840E-01  6.2049E-02  1.2359E+00  1.6750E+00  7.5223E-01  1.0648E+00  5.1586E+00  1.0000E-02  6.3064E-01  1.1739E+00
             1.8893E+00
 PARAMETER:  4.7023E-02 -2.6798E+00  3.1176E-01  6.1579E-01 -1.8471E-01  1.6275E-01  1.7407E+00 -2.3030E+01 -3.6101E-01  2.6036E-01
             7.3622E-01
 GRADIENT:  -4.4202E-01 -1.9978E+00  1.0703E+00  2.4061E+01 -4.6502E+00  1.4900E+00 -8.2871E+00  0.0000E+00  4.8169E+00  2.7955E+00
            -3.4043E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1466.64026706240        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1329
 NPARAMETR:  9.4640E-01  2.0049E-02  1.3515E+00  1.7109E+00  7.9151E-01  1.0597E+00  8.8507E+00  1.0000E-02  6.0656E-01  1.2011E+00
             1.9035E+00
 PARAMETER:  4.3909E-02 -3.8061E+00  4.0043E-01  6.3658E-01 -1.3389E-01  1.5796E-01  2.2824E+00 -3.3690E+01 -4.0033E-01  2.8546E-01
             7.4121E-01
 GRADIENT:  -5.1891E+00  3.8155E+01 -8.8824E-01 -1.9627E+02 -4.0396E-01 -1.4535E-01  6.0457E+01  0.0000E+00 -9.6466E+01  1.3802E+00
            -3.6690E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1329
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1169E-03  1.2104E-02 -1.8477E-05 -2.7933E-02 -3.0412E-02
 SE:             2.9427E-02  7.7891E-03  1.1307E-04  2.6423E-02  2.2586E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1565E-01  1.2021E-01  8.7020E-01  2.9045E-01  1.7814E-01

 ETASHRINKSD(%)  1.4147E+00  7.3905E+01  9.9621E+01  1.1479E+01  2.4334E+01
 ETASHRINKVR(%)  2.8094E+00  9.3191E+01  9.9999E+01  2.1640E+01  4.2746E+01
 EBVSHRINKSD(%)  1.2734E+00  8.2416E+01  9.9592E+01  1.0882E+01  2.1880E+01
 EBVSHRINKVR(%)  2.5306E+00  9.6908E+01  9.9998E+01  2.0580E+01  3.8973E+01
 RELATIVEINF(%)  9.7318E+01  1.9782E+00  1.1909E-04  3.2825E+01  4.4272E+00
 EPSSHRINKSD(%)  3.4826E+01
 EPSSHRINKVR(%)  5.7523E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1466.6402670623982     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -731.48944049865997     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.49
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.85
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1466.640       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.45E-01  2.01E-02  1.35E+00  1.71E+00  7.91E-01  1.06E+00  8.87E+00  1.00E-02  6.06E-01  1.20E+00  1.90E+00
 


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
+       -1.08E+02  9.69E+05
 
 TH 3
+        7.38E-01  1.42E+02  1.14E+02
 
 TH 4
+       -2.61E+01  2.35E+03 -5.39E+01  3.91E+03
 
 TH 5
+        7.92E+00 -1.35E+03 -2.86E+02  1.30E+01  8.40E+02
 
 TH 6
+        8.81E+00 -8.99E+01  2.34E+00 -5.57E+00 -2.34E+00  1.65E+02
 
 TH 7
+       -2.54E-01  2.45E+03  4.66E-01 -8.16E+01 -4.53E+00 -3.34E-01  1.42E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.81E+00  8.99E+02  1.46E+04 -1.43E+02  1.49E+02  2.45E+00  3.38E+00  0.00E+00  3.27E+04
 
 TH10
+       -1.37E+00  8.20E+02  1.76E+01 -6.74E+01 -1.03E+02 -1.29E+00  3.13E+00  0.00E+00 -2.61E+01  8.57E+01
 
 TH11
+       -1.24E+01 -2.25E+01 -5.03E+00 -1.85E+01 -9.26E+00  3.44E+00 -3.96E-02  0.00E+00  5.61E+03  2.18E+01  7.67E+01
 
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
 #CPUT: Total CPU Time in Seconds,       23.412
Stop Time:
Wed Sep 29 12:45:41 CDT 2021
