Wed Sep 29 20:18:57 CDT 2021
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
$DATA ../../../../data/spa/D/dat75.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m75.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   14874.0521841245        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9483E+02  1.9820E+02 -9.1849E+01 -3.3315E+01  3.0175E+02 -2.0453E+03 -8.9059E+02 -5.1778E+01 -1.7144E+03 -5.4536E+02
            -2.7544E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -590.121881473224        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3062E+00  1.0647E+00  9.7090E-01  1.6527E+00  1.2876E+00  1.9310E+00  1.2629E+00  9.7431E-01  1.4185E+00  1.0704E+00
             1.4256E+01
 PARAMETER:  3.6711E-01  1.6273E-01  7.0473E-02  6.0239E-01  3.5281E-01  7.5802E-01  3.3339E-01  7.3969E-02  4.4959E-01  1.6800E-01
             2.7571E+00
 GRADIENT:  -1.6218E+01  1.3544E+01 -6.0827E+00  1.7639E+01 -7.7220E+00  3.8008E+01 -5.3900E-01  4.9150E+00  5.5901E+00  2.2064E+00
             1.0681E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -601.991493767753        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2953E+00  8.9900E-01  1.1668E+00  1.7683E+00  2.3099E+00  1.7722E+00  2.6803E+00  4.1518E-01  1.3954E+00  2.5628E+00
             1.3160E+01
 PARAMETER:  3.5874E-01 -6.4739E-03  2.5430E-01  6.7002E-01  9.3722E-01  6.7222E-01  1.0859E+00 -7.7905E-01  4.3318E-01  1.0411E+00
             2.6772E+00
 GRADIENT:  -2.9615E+00  1.1347E+01 -6.0506E+00  2.9216E+01 -3.9747E+00  5.3577E+00  3.5435E+00  5.0249E-01  1.2606E+01  3.8810E+00
             6.5493E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -610.157582221259        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.1578E+00  4.5998E-01  1.4708E+00  1.6576E+00  3.8128E+00  1.5538E+00  2.7464E+00  1.4522E+00  8.3212E-01  7.6561E+00
             1.2069E+01
 PARAMETER:  2.4652E-01 -6.7657E-01  4.8580E-01  6.0539E-01  1.4384E+00  5.4068E-01  1.1103E+00  4.7305E-01 -8.3784E-02  2.1355E+00
             2.5906E+00
 GRADIENT:  -9.0222E+00  4.9611E+00  1.3577E+01 -3.9979E+00 -6.7472E+00 -2.1326E+00  1.8843E+00 -7.6648E-01  7.1837E-01  8.7128E+00
             3.9207E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -638.958642421585        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0454E+00  1.1994E-01  4.0423E-01  1.4632E+00  1.0785E+01  1.2013E+00  1.2333E+00  1.1390E-01  3.0772E-01  6.4639E+00
             1.2235E+01
 PARAMETER:  1.4442E-01 -2.0208E+00 -8.0577E-01  4.8065E-01  2.4782E+00  2.8339E-01  3.0971E-01 -2.0724E+00 -1.0786E+00  1.9662E+00
             2.6043E+00
 GRADIENT:  -1.1852E+01  2.7493E+00  8.3669E+00  6.7449E+01  1.7544E+01 -8.9477E+01  5.9119E-02 -5.4713E-03  3.8038E+00 -2.5964E+01
            -5.1762E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -672.106882557619        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      425
 NPARAMETR:  1.1363E+00  1.1767E-01  4.0212E-01  1.4453E+00  1.4148E+01  1.5603E+00  4.1931E-01  6.7921E+00  8.2357E-01  6.4531E+00
             1.2042E+01
 PARAMETER:  2.2778E-01 -2.0399E+00 -8.1100E-01  4.6832E-01  2.7495E+00  5.4490E-01 -7.6914E-01  2.0158E+00 -9.4108E-02  1.9646E+00
             2.5884E+00
 GRADIENT:   5.2577E+01  4.1086E+01 -1.0793E+01 -3.1831E+01 -3.7990E+00 -6.5364E+00  8.2921E-01  8.2687E+01  1.4098E+01  2.4112E-01
             1.0450E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -705.660691715357        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.0462E+00  1.0755E-01  4.0061E-01  1.4271E+00  1.3214E+01  1.7939E+00  1.0000E-02  3.5565E+00  8.0229E-01  6.4747E+00
             1.0229E+01
 PARAMETER:  1.4520E-01 -2.1298E+00 -8.1477E-01  4.5566E-01  2.6813E+00  6.8439E-01 -1.6221E+01  1.3688E+00 -1.2028E-01  1.9679E+00
             2.4252E+00
 GRADIENT:   2.7220E+01  1.9586E+01  1.0328E+01  5.5043E+01  6.6256E-01  2.2723E+01  0.0000E+00 -3.3346E+00  4.5790E+00 -3.1357E+00
            -5.5469E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -709.884217437042        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      570
 NPARAMETR:  9.1927E-01  8.0311E-02  3.9540E-01  1.3442E+00  1.3085E+01  1.6736E+00  1.0000E-02  3.4331E+00  6.0612E-01  6.5432E+00
             9.8730E+00
 PARAMETER:  1.5822E-02 -2.4218E+00 -8.2786E-01  3.9580E-01  2.6715E+00  6.1498E-01 -2.0365E+01  1.3335E+00 -4.0068E-01  1.9784E+00
             2.3898E+00
 GRADIENT:  -4.0338E+01  8.1164E+00  2.5222E+01  7.0402E+01  1.2906E+00 -1.7319E+00  0.0000E+00 -1.3220E+01  9.4613E-01 -3.0931E+00
            -2.5311E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -714.845283891221        NO. OF FUNC. EVALS.: 112
 CUMULATIVE NO. OF FUNC. EVALS.:      682
 NPARAMETR:  9.1041E-01  5.0005E-02  3.8347E-01  1.2176E+00  1.1830E+01  1.6011E+00  1.0000E-02  3.7135E+00  4.0228E-01  6.5901E+00
             9.9532E+00
 PARAMETER:  6.8513E-03 -2.8800E+00 -8.5747E-01  2.9407E-01  2.5542E+00  5.7427E-01 -3.0956E+01  1.4129E+00 -8.1663E-01  1.9975E+00
             2.3741E+00
 GRADIENT:   2.3038E+02  6.0209E+00  4.5288E+00 -5.3317E+01 -7.9094E+00  3.5662E+01  0.0000E+00  1.5829E+00 -2.9074E+01  9.2774E+00
            -1.3406E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      682
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3285E-02 -2.9010E-05 -7.7386E-02 -2.0929E-02 -6.1302E-03
 SE:             2.9201E-02  5.8375E-06  2.2234E-02  7.8962E-03  3.7789E-03
 N:                     100         100         100         100         100

 P VAL.:         4.2522E-01  6.7210E-07  5.0048E-04  8.0359E-03  1.0476E-01

 ETASHRINKSD(%)  2.1731E+00  9.9980E+01  2.5514E+01  7.3547E+01  8.7340E+01
 ETASHRINKVR(%)  4.2990E+00  1.0000E+02  4.4518E+01  9.3002E+01  9.8397E+01
 EBVSHRINKSD(%)  5.5138E+00  9.9973E+01  2.3633E+01  7.3493E+01  9.1788E+01
 EBVSHRINKVR(%)  1.0724E+01  1.0000E+02  4.1681E+01  9.2974E+01  9.9326E+01
 RELATIVEINF(%)  4.2968E+01  4.8574E-06  1.4338E+01  1.1523E+00  4.1624E-01
 EPSSHRINKSD(%)  9.7537E+00
 EPSSHRINKVR(%)  1.8556E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -714.84528389122124     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       20.305542672516935     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.64
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -714.845       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.11E-01  5.08E-02  3.84E-01  1.21E+00  1.16E+01  1.61E+00  1.00E-02  3.72E+00  4.00E-01  6.67E+00  9.72E+00
 


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
+        7.48E+04
 
 TH 2
+       -1.78E+02  2.79E+04
 
 TH 3
+       -2.75E+02 -1.98E+02  3.31E+02
 
 TH 4
+       -1.06E+02  2.50E+02 -1.67E+02  3.62E+02
 
 TH 5
+        9.03E-01 -1.43E+02  9.49E-01  5.64E+01  6.82E-01
 
 TH 6
+        9.97E+01  1.02E+01 -2.14E+01 -4.12E+01  9.70E-02  7.88E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.37E+00 -5.98E+00 -4.03E+01  7.43E+00 -2.03E-01  3.32E+00  0.00E+00  2.99E+01
 
 TH 9
+        7.38E+01  3.84E+01  6.72E+01 -4.95E+01 -2.88E-01  6.68E+00  0.00E+00  2.97E+00  5.80E+03
 
 TH10
+       -9.63E-01 -1.91E+01 -1.41E+00  8.06E+00  8.64E-02  4.86E+01  0.00E+00  8.62E+00  7.74E-01  3.21E+00
 
 TH11
+       -8.59E+00  1.08E+01  7.56E+00 -1.80E+01  7.99E-01  1.84E+00  0.00E+00  3.08E-01  2.10E+00  1.89E-01  5.78E+00
 
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
 #CPUT: Total CPU Time in Seconds,       17.400
Stop Time:
Wed Sep 29 20:19:16 CDT 2021
