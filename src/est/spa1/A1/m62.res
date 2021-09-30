Wed Sep 29 22:41:16 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat62.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1665.15077051203        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0543E+02  5.8932E+01  1.7243E+01  1.0821E+02  1.2325E+02  7.7348E+01 -3.0132E+01 -1.2569E+01 -7.4719E+00 -4.8156E+01
            -8.1919E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1860.90574855501        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0108E+00  9.4850E-01  9.1545E-01  1.0191E+00  8.9153E-01  8.8340E-01  1.0773E+00  9.7572E-01  1.0462E+00  1.0505E+00
             1.5640E+00
 PARAMETER:  1.1073E-01  4.7128E-02  1.1657E-02  1.1889E-01 -1.4815E-02 -2.3974E-02  1.7450E-01  7.5416E-02  1.4513E-01  1.4925E-01
             5.4724E-01
 GRADIENT:   1.5694E+02  1.8888E+01  1.9876E+00  4.1295E+01  2.1895E+01  1.4554E+00 -5.8150E+00  5.9926E+00  9.8800E+00  7.4448E-01
            -9.2847E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1869.53823408484        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      224
 NPARAMETR:  1.0095E+00  7.2020E-01  6.3489E-01  1.1528E+00  6.5171E-01  8.9486E-01  1.4461E+00  5.4913E-01  8.9492E-01  8.2844E-01
             1.5262E+00
 PARAMETER:  1.0942E-01 -2.2823E-01 -3.5431E-01  2.4220E-01 -3.2816E-01 -1.1092E-02  4.6885E-01 -4.9942E-01 -1.1023E-02 -8.8213E-02
             5.2280E-01
 GRADIENT:  -5.5324E+01  1.6294E+01 -3.0606E+01  7.3442E+01  5.0071E+01 -1.4520E+01 -6.5449E+00  4.9855E+00 -8.0971E+00 -2.3700E+00
            -1.0410E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1881.72517773745        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      400
 NPARAMETR:  1.0351E+00  5.5173E-01  6.1926E-01  1.2163E+00  5.7237E-01  9.2277E-01  1.8440E+00  1.7425E-01  8.5328E-01  7.9348E-01
             1.7107E+00
 PARAMETER:  1.3454E-01 -4.9470E-01 -3.7923E-01  2.9584E-01 -4.5797E-01  1.9630E-02  7.1194E-01 -1.6473E+00 -5.8666E-02 -1.3132E-01
             6.3691E-01
 GRADIENT:   1.4615E+01  1.5785E+01  1.0812E+01  2.5787E+01 -1.4233E+01  8.7123E-01  1.9524E+00  2.6721E-01 -3.5919E+00 -2.4335E+00
             3.4460E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1882.76870652605        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      575
 NPARAMETR:  1.0293E+00  4.3701E-01  6.1302E-01  1.2540E+00  5.4156E-01  9.1921E-01  1.9819E+00  7.9771E-02  8.5908E-01  8.3419E-01
             1.6978E+00
 PARAMETER:  1.2887E-01 -7.2780E-01 -3.8936E-01  3.2634E-01 -5.1331E-01  1.5756E-02  7.8406E-01 -2.4286E+00 -5.1898E-02 -8.1290E-02
             6.2932E-01
 GRADIENT:   6.5744E+00  1.1787E+00  4.2488E+00 -2.6124E+00 -7.2950E+00  4.3665E-01 -9.6980E-01  6.5484E-02  1.3504E+00  2.9898E-01
            -2.8821E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1883.10443099374        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      752
 NPARAMETR:  1.0179E+00  2.6721E-01  7.0208E-01  1.3699E+00  5.5599E-01  9.0920E-01  2.6189E+00  1.4689E-02  8.2196E-01  9.2377E-01
             1.7157E+00
 PARAMETER:  1.1771E-01 -1.2197E+00 -2.5370E-01  4.1477E-01 -4.8701E-01  4.8056E-03  1.0627E+00 -4.1207E+00 -9.6065E-02  2.0706E-02
             6.3982E-01
 GRADIENT:  -7.9707E+00  3.7390E+00  3.5870E-01  1.3934E+01 -4.9317E-01 -1.4530E+00  1.5585E+00  4.5905E-04 -2.7536E+00  2.1224E-01
             1.3448E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1884.13820330562        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      930
 NPARAMETR:  1.0166E+00  1.1095E-01  6.8631E-01  1.4420E+00  5.1418E-01  9.0901E-01  4.0424E+00  1.0000E-02  8.0258E-01  9.1293E-01
             1.7138E+00
 PARAMETER:  1.1648E-01 -2.0987E+00 -2.7643E-01  4.6601E-01 -5.6518E-01  4.5991E-03  1.4968E+00 -7.3943E+00 -1.1993E-01  8.9063E-03
             6.3873E-01
 GRADIENT:   2.8503E+00  1.4063E+00  1.2527E+01  2.3195E+01 -2.3965E+01 -2.1329E-01 -2.3594E+00  0.0000E+00 -6.1107E-01 -1.3035E+00
            -2.4036E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1885.38511157594        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1107
 NPARAMETR:  1.0125E+00  3.8879E-02  7.1014E-01  1.4686E+00  5.2139E-01  9.0947E-01  6.4979E+00  1.0000E-02  7.8909E-01  9.2251E-01
             1.7272E+00
 PARAMETER:  1.1243E-01 -3.1473E+00 -2.4229E-01  4.8430E-01 -5.5126E-01  5.1073E-03  1.9715E+00 -1.1537E+01 -1.3688E-01  1.9348E-02
             6.4649E-01
 GRADIENT:  -5.1675E-01 -5.9029E-01  4.0852E+00 -2.1961E+00 -3.3239E+00  8.3203E-01 -1.5411E+00  0.0000E+00  1.2683E+00 -1.9344E+00
             4.0123E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1885.46148042658        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1172
 NPARAMETR:  1.0121E+00  3.6824E-02  7.0427E-01  1.4722E+00  5.1952E-01  9.0662E-01  6.7012E+00  1.0000E-02  7.8558E-01  9.4254E-01
             1.7136E+00
 PARAMETER:  1.1195E-01 -3.2008E+00 -2.5073E-01  4.8665E-01 -5.5472E-01  1.8904E-03  2.0027E+00 -1.1797E+01 -1.4140E-01  4.1824E-02
             6.3854E-01
 GRADIENT:  -1.0958E+00  1.9324E+01 -1.5336E+00 -1.2622E+02  1.2386E+02 -2.8623E-01  2.7959E+01  0.0000E+00 -2.1847E-01  2.9967E+00
            -1.0845E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1172
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4227E-03  1.5533E-02 -1.0590E-04 -1.3378E-02 -1.0461E-02
 SE:             2.9645E-02  8.8102E-03  2.0822E-04  2.8311E-02  2.4317E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6172E-01  7.7892E-02  6.1103E-01  6.3655E-01  6.6706E-01

 ETASHRINKSD(%)  6.8596E-01  7.0485E+01  9.9302E+01  5.1550E+00  1.8535E+01
 ETASHRINKVR(%)  1.3672E+00  9.1289E+01  9.9995E+01  1.0044E+01  3.3634E+01
 EBVSHRINKSD(%)  1.0684E+00  8.0443E+01  9.9314E+01  4.6542E+00  1.5818E+01
 EBVSHRINKVR(%)  2.1253E+00  9.6175E+01  9.9995E+01  9.0917E+00  2.9134E+01
 RELATIVEINF(%)  9.7276E+01  2.0726E+00  5.0287E-04  4.1638E+01  7.8609E+00
 EPSSHRINKSD(%)  3.0103E+01
 EPSSHRINKVR(%)  5.1143E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1885.4614804265802     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -966.52294722190754     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1885.461       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  3.69E-02  7.04E-01  1.47E+00  5.20E-01  9.07E-01  6.70E+00  1.00E-02  7.86E-01  9.43E-01  1.71E+00
 


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
+        1.29E+03
 
 TH 2
+       -6.26E+01  2.60E+05
 
 TH 3
+       -2.24E+01 -4.11E+02  1.15E+03
 
 TH 4
+       -1.22E+01  2.94E+01 -1.42E+02  4.12E+03
 
 TH 5
+        5.24E+04  2.60E+02 -3.56E+04 -8.42E+03  2.42E+04
 
 TH 6
+        2.02E+00  5.28E+00  7.91E+00 -8.21E+00 -6.54E+04  2.33E+02
 
 TH 7
+        1.49E-01  2.45E+03 -5.17E+00 -4.45E+00  8.52E+00  1.02E-01  3.18E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.09E+00 -3.64E+01  9.63E+01 -4.56E+01 -5.37E+04  5.61E+00 -6.66E-01  0.00E+00  3.17E+02
 
 TH10
+       -5.92E+00 -2.31E+02  4.89E+01  1.04E+01 -6.33E+04  7.84E+00 -3.09E+00  0.00E+00  5.15E+01  1.76E+02
 
 TH11
+       -1.53E+01 -2.69E+01  3.70E+00 -1.76E+01 -5.45E+03  3.48E+00 -2.80E-01  0.00E+00  2.25E+01  3.51E+01  1.52E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.017
Stop Time:
Wed Sep 29 22:41:44 CDT 2021
