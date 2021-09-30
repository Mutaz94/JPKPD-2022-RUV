Wed Sep 29 11:52:13 CDT 2021
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
$DATA ../../../../data/spa/A1/dat3.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1126.71285051310        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0476E+02  1.7500E+01 -1.3686E+01  6.7208E+01  1.4958E+02  6.1952E+01 -2.6434E+01 -6.2699E+00 -3.3611E+01 -6.1346E+01
            -9.6799E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1416.67467292143        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.1995E+00  1.0031E+00  9.9271E-01  1.0434E+00  8.9701E-01  1.1760E+00  1.0378E+00  9.9012E-01  1.0695E+00  1.0648E+00
             2.3446E+00
 PARAMETER:  2.8194E-01  1.0306E-01  9.2686E-02  1.4245E-01 -8.6920E-03  2.6216E-01  1.3707E-01  9.0069E-02  1.6714E-01  1.6283E-01
             9.5212E-01
 GRADIENT:   4.6892E+02  1.7472E+01  6.0686E+00  2.1299E+01 -9.7183E+00  4.7421E+01  1.0483E+00  4.2034E+00  4.2075E+00  1.2994E+00
             1.6030E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1438.04533369374        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0754E+00  7.4911E-01  5.8146E-01  1.1726E+00  6.1837E-01  1.0208E+00  1.1242E+00  6.8508E-02  9.5964E-01  7.5275E-01
             2.1269E+00
 PARAMETER:  1.7271E-01 -1.8887E-01 -4.4221E-01  2.5926E-01 -3.8067E-01  1.2057E-01  2.1707E-01 -2.5808E+00  5.8802E-02 -1.8403E-01
             8.5466E-01
 GRADIENT:   2.3967E+02  3.7124E+00 -5.5223E+01  1.2673E+02  1.1362E+02  2.4732E+01 -3.5998E+00  4.4288E-02 -7.1267E+00 -1.2764E+01
            -3.9116E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1446.13459350362        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.0288E+00  6.9179E-01  6.8645E-01  1.2027E+00  6.2960E-01  9.5642E-01  1.2461E+00  4.5024E-02  9.7364E-01  8.2446E-01
             2.3038E+00
 PARAMETER:  1.2843E-01 -2.6847E-01 -2.7622E-01  2.8460E-01 -3.6267E-01  5.5442E-02  3.2003E-01 -3.0006E+00  7.3286E-02 -9.3033E-02
             9.3456E-01
 GRADIENT:  -1.1853E+01  2.0188E+01  8.5707E-01  2.5668E+01 -4.3834E+00 -4.5974E+00  1.2224E+00  1.6353E-02  2.9367E-01 -2.9895E+00
             2.2035E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1450.64495850980        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      471
 NPARAMETR:  1.0351E+00  2.9646E-01  6.8392E-01  1.4059E+00  5.2082E-01  9.7237E-01  1.6370E+00  1.8737E-02  8.9826E-01  8.1093E-01
             2.2883E+00
 PARAMETER:  1.3449E-01 -1.1159E+00 -2.7992E-01  4.4069E-01 -5.5235E-01  7.1978E-02  5.9287E-01 -3.8773E+00 -7.2966E-03 -1.0958E-01
             9.2780E-01
 GRADIENT:   1.5946E+01  7.7620E+00  1.9504E+01  2.4152E+01 -3.0878E+01  3.3927E+00 -1.1188E+00  3.0918E-03  1.4063E+00 -6.4836E-01
            -2.8250E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1452.20469527550        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      646
 NPARAMETR:  1.0188E+00  1.1436E-01  7.1543E-01  1.4968E+00  5.1456E-01  9.6147E-01  3.1287E+00  1.0000E-02  8.4515E-01  8.3488E-01
             2.2825E+00
 PARAMETER:  1.1858E-01 -2.0684E+00 -2.3488E-01  5.0333E-01 -5.6444E-01  6.0712E-02  1.2406E+00 -6.6815E+00 -6.8243E-02 -8.0464E-02
             9.2527E-01
 GRADIENT:  -5.6548E+00  1.0752E+00 -1.5591E+00  1.2407E+01 -1.6305E-03  1.4885E+00 -1.4552E-02  0.0000E+00 -1.1491E+00  9.0557E-01
            -3.7708E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1452.44134857531        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      821
 NPARAMETR:  1.0193E+00  4.2069E-02  7.2412E-01  1.5312E+00  5.0506E-01  9.5385E-01  5.8661E+00  1.0000E-02  8.3184E-01  8.2471E-01
             2.3096E+00
 PARAMETER:  1.1910E-01 -3.0684E+00 -2.2280E-01  5.2604E-01 -5.8307E-01  5.2755E-02  1.8692E+00 -9.9021E+00 -8.4119E-02 -9.2719E-02
             9.3709E-01
 GRADIENT:   9.0363E-01  5.3658E-01  3.4374E+00  6.3104E+00 -6.0619E+00 -5.1274E-01  3.5745E-01  0.0000E+00 -2.1794E-01 -3.0244E-01
             6.6826E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1452.52600974051        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      998
 NPARAMETR:  1.0176E+00  1.0377E-02  7.3213E-01  1.5456E+00  5.0537E-01  9.5446E-01  1.2592E+01  1.0000E-02  8.3053E-01  8.3391E-01
             2.3083E+00
 PARAMETER:  1.1745E-01 -4.4682E+00 -2.1180E-01  5.3540E-01 -5.8247E-01  5.3389E-02  2.6330E+00 -1.4605E+01 -8.5694E-02 -8.1634E-02
             9.3650E-01
 GRADIENT:   1.1743E-01  2.8416E-01  5.9872E-01 -1.5088E-01 -1.1718E+00  9.7435E-02  3.4289E-01  0.0000E+00  1.2285E+00  3.5597E-01
             7.7187E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1452.53608806829        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1134
 NPARAMETR:  1.0178E+00  1.0217E-02  7.3186E-01  1.5449E+00  5.0557E-01  9.5441E-01  1.2477E+01  1.0000E-02  8.2670E-01  8.3173E-01
             2.3064E+00
 PARAMETER:  1.1762E-01 -4.4837E+00 -2.1216E-01  5.3493E-01 -5.8207E-01  5.3339E-02  2.6239E+00 -1.4605E+01 -9.0310E-02 -8.4242E-02
             9.3568E-01
 GRADIENT:   6.0945E-01  2.8487E-01 -1.8839E-01 -1.3347E+00  2.5954E-01  4.6604E-02  2.1448E-01  0.0000E+00 -6.0414E-02 -2.9663E-02
            -3.0246E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1134
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.2002E-04  1.4678E-03  8.5434E-06 -1.0919E-02 -1.5500E-02
 SE:             2.9247E-02  1.9699E-03  1.8494E-04  2.7317E-02  2.1767E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8854E-01  4.5621E-01  9.6315E-01  6.8935E-01  4.7641E-01

 ETASHRINKSD(%)  2.0203E+00  9.3400E+01  9.9380E+01  8.4854E+00  2.7079E+01
 ETASHRINKVR(%)  3.9998E+00  9.9564E+01  9.9996E+01  1.6251E+01  4.6826E+01
 EBVSHRINKSD(%)  2.0926E+00  9.3615E+01  9.9362E+01  7.8418E+00  2.6103E+01
 EBVSHRINKVR(%)  4.1414E+00  9.9592E+01  9.9996E+01  1.5069E+01  4.5393E+01
 RELATIVEINF(%)  8.6591E+01  1.9218E-02  1.9905E-04  6.0100E+00  2.2790E+00
 EPSSHRINKSD(%)  3.3904E+01
 EPSSHRINKVR(%)  5.6313E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1452.5360880682947     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -717.38526150455652     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.27
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1452.536       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.02E-02  7.32E-01  1.54E+00  5.06E-01  9.54E-01  1.25E+01  1.00E-02  8.27E-01  8.32E-01  2.31E+00
 


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
+        1.13E+03
 
 TH 2
+       -2.84E+01  1.12E+06
 
 TH 3
+        1.62E+00  5.10E+01  8.55E+02
 
 TH 4
+       -3.83E+01 -6.14E+04 -9.95E+01  5.81E+02
 
 TH 5
+        3.79E+01  1.73E+05 -1.60E+03 -1.27E+02  3.42E+03
 
 TH 6
+       -2.75E+00 -1.23E+01  8.54E+00 -8.33E+00 -7.64E-01  1.99E+02
 
 TH 7
+       -3.10E-01  1.58E+03  2.69E+00 -8.66E+01  2.44E+02  1.21E-01  2.22E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.87E+00 -1.32E+02  2.50E+01 -9.06E+00 -7.79E+00  5.49E+00  5.63E-01  0.00E+00  2.07E+02
 
 TH10
+       -4.84E+00 -1.01E+02 -1.39E+01 -7.57E-01 -6.36E+01 -7.08E-01  1.38E+00  0.00E+00 -1.57E+00  9.81E+01
 
 TH11
+       -1.29E+01 -4.49E+01 -1.10E+01 -9.41E+00  5.06E+00  2.16E+00 -3.35E+01  0.00E+00  1.06E+01  2.02E+01  5.32E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.944
Stop Time:
Wed Sep 29 11:52:35 CDT 2021
