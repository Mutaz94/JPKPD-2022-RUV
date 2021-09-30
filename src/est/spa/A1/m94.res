Wed Sep 29 12:28:11 CDT 2021
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
$DATA ../../../../data/spa/A1/dat94.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1482.34916881866        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2037E+02 -1.7311E+01 -2.9690E+01  1.2622E+01  1.1649E+02  4.2133E+01 -1.2960E+01  4.7489E+00 -3.4359E+01 -1.5727E+01
            -2.9120E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1511.84020167772        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0072E+00  9.2721E-01  1.0187E+00  1.0850E+00  8.4313E-01  9.1881E-01  1.0102E+00  8.8604E-01  1.1908E+00  7.7666E-01
             2.3760E+00
 PARAMETER:  1.0718E-01  2.4423E-02  1.1857E-01  1.8162E-01 -7.0640E-02  1.5322E-02  1.1016E-01 -2.0998E-02  2.7465E-01 -1.5275E-01
             9.6541E-01
 GRADIENT:   6.2731E+01  1.9388E+01  2.4415E+01  1.3225E+01 -6.4619E+01 -1.3518E+01  9.3254E+00  5.1677E+00  3.0182E+01  1.6057E+01
             1.5123E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1521.32529654952        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.9583E-01  7.7260E-01  5.0280E-01  1.1409E+00  5.7168E-01  9.9470E-01  1.2876E+00  4.6062E-01  9.7518E-01  3.9796E-01
             2.2517E+00
 PARAMETER:  9.5816E-02 -1.5799E-01 -5.8757E-01  2.3181E-01 -4.5917E-01  9.4689E-02  3.5278E-01 -6.7517E-01  7.4869E-02 -8.2141E-01
             9.1168E-01
 GRADIENT:   2.5955E+01 -1.2691E+01 -6.3234E+01  8.8354E+01  9.1814E+01  1.0865E+01  6.8679E+00  4.4139E+00  1.6036E+01  6.2235E+00
             1.4316E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1540.86298540528        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      277
 NPARAMETR:  9.8487E-01  8.9106E-01  4.4330E-01  1.0368E+00  5.6911E-01  9.8408E-01  1.1278E+00  2.4468E-01  1.0413E+00  4.6110E-01
             1.8692E+00
 PARAMETER:  8.4750E-02 -1.5349E-02 -7.1350E-01  1.3611E-01 -4.6369E-01  8.3954E-02  2.2031E-01 -1.3078E+00  1.4046E-01 -6.7415E-01
             7.2553E-01
 GRADIENT:  -6.6091E+01 -1.4048E+01 -3.6736E+01  1.9002E+01  3.5331E+01 -4.5837E+00  1.4659E+00  1.0107E+00  1.4579E+01  6.1685E+00
             8.3818E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1555.12985363413        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  1.0246E+00  5.9946E-01  4.9592E-01  1.1998E+00  4.8985E-01  9.9767E-01  1.5867E+00  1.9667E-01  9.1229E-01  5.7609E-01
             1.4927E+00
 PARAMETER:  1.2433E-01 -4.1173E-01 -6.0134E-01  2.8215E-01 -6.1365E-01  9.7670E-02  5.6165E-01 -1.5262E+00  8.2038E-03 -4.5149E-01
             5.0058E-01
 GRADIENT:   2.9370E+01  1.9251E+01  1.5606E+01 -5.0623E+00 -3.9778E+01 -3.1510E-01  1.4202E+00  5.1496E-01  3.8842E+00  4.2468E+00
             4.8641E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1557.73548241877        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      627
 NPARAMETR:  1.0086E+00  4.0318E-01  5.5323E-01  1.3109E+00  4.8258E-01  9.9475E-01  2.1186E+00  5.2358E-02  8.5616E-01  6.3732E-01
             1.4806E+00
 PARAMETER:  1.0859E-01 -8.0838E-01 -4.9198E-01  3.7070E-01 -6.2860E-01  9.4737E-02  8.5076E-01 -2.8496E+00 -5.5304E-02 -3.5048E-01
             4.9246E-01
 GRADIENT:   2.0969E+00  2.1499E+00  5.6197E+00 -3.4597E+00 -9.5343E+00  2.4343E-01 -6.9539E-01  3.0560E-03 -8.1747E-01 -9.8848E-01
             1.4424E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1557.85453642599        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      786
 NPARAMETR:  1.0064E+00  3.5524E-01  5.7428E-01  1.3413E+00  4.8674E-01  9.9265E-01  2.3231E+00  3.6287E-02  8.5015E-01  6.6867E-01
             1.4741E+00
 PARAMETER:  1.0640E-01 -9.3495E-01 -4.5463E-01  3.9361E-01 -6.2003E-01  9.2623E-02  9.4291E-01 -3.2163E+00 -6.2344E-02 -3.0247E-01
             4.8808E-01
 GRADIENT:   1.8650E+02  2.3744E+01  1.7085E+01  2.2791E+02  5.4963E+01  1.6276E+01  1.2505E+01  1.4031E-02  5.4794E+00  1.1891E+00
             3.8703E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1557.85578468814        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      970             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0065E+00  3.5532E-01  5.7395E-01  1.3406E+00  4.8670E-01  9.9307E-01  2.3302E+00  1.0000E-02  8.4949E-01  6.6827E-01
             1.4742E+00
 PARAMETER:  1.0646E-01 -9.3472E-01 -4.5522E-01  3.9311E-01 -6.2010E-01  9.3043E-02  9.4597E-01 -5.0133E+00 -6.3121E-02 -3.0306E-01
             4.8813E-01
 GRADIENT:   1.8672E+02  2.3796E+01  1.6883E+01  2.2668E+02  5.5438E+01  1.6433E+01  1.2960E+01  0.0000E+00  5.3569E+00  1.1525E+00
             3.9033E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1557.85578468814        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  1.0065E+00  3.5532E-01  5.7395E-01  1.3406E+00  4.8670E-01  9.9307E-01  2.3302E+00  1.0000E-02  8.4949E-01  6.6827E-01
             1.4742E+00
 PARAMETER:  1.0646E-01 -9.3472E-01 -4.5522E-01  3.9311E-01 -6.2010E-01  9.3043E-02  9.4597E-01 -5.0133E+00 -6.3121E-02 -3.0306E-01
             4.8813E-01
 GRADIENT:   6.6444E-01  1.4845E-02  9.0967E-02 -2.5944E+00  3.8125E-01  4.1764E-02  1.3668E-01  0.0000E+00  2.8886E-02  1.4212E-02
            -1.6928E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1028
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4327E-04  3.5315E-02 -4.4991E-04 -2.1752E-02  1.1900E-02
 SE:             2.9648E-02  1.8494E-02  2.4418E-04  2.6406E-02  2.1099E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7731E-01  5.6185E-02  6.5393E-02  4.1007E-01  5.7277E-01

 ETASHRINKSD(%)  6.7618E-01  3.8044E+01  9.9182E+01  1.1538E+01  2.9315E+01
 ETASHRINKVR(%)  1.3478E+00  6.1615E+01  9.9993E+01  2.1745E+01  5.0036E+01
 EBVSHRINKSD(%)  9.3892E-01  4.4348E+01  9.9167E+01  9.9060E+00  2.5238E+01
 EBVSHRINKVR(%)  1.8690E+00  6.9028E+01  9.9993E+01  1.8831E+01  4.4107E+01
 RELATIVEINF(%)  9.7390E+01  6.6578E+00  3.9134E-04  2.8281E+01  2.7387E+00
 EPSSHRINKSD(%)  4.1092E+01
 EPSSHRINKVR(%)  6.5299E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1557.8557846881381     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -822.70495812439992     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.65
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.32
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1557.856       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  3.55E-01  5.74E-01  1.34E+00  4.87E-01  9.93E-01  2.33E+00  1.00E-02  8.49E-01  6.68E-01  1.47E+00
 


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
+        1.09E+03
 
 TH 2
+       -2.73E+01  5.63E+02
 
 TH 3
+        2.48E+01  5.45E+02  2.49E+03
 
 TH 4
+       -1.32E+01  2.63E+02 -3.58E+02  7.05E+02
 
 TH 5
+        7.72E+00 -1.12E+03 -3.65E+03  2.55E+02  5.88E+03
 
 TH 6
+        5.31E-01 -4.95E+00  7.40E+00 -4.69E+00 -2.83E+00  1.94E+02
 
 TH 7
+        1.30E+00  3.95E+01  3.32E-02 -3.36E+00 -1.39E+01 -4.39E-02  9.17E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.46E+00 -1.80E+01 -2.26E+01 -6.19E+00  4.04E+01 -1.38E+00  6.89E+00  0.00E+00  1.88E+02
 
 TH10
+       -2.09E+00  2.68E+01 -1.35E+02 -2.90E+01  7.63E+01  5.18E-01  5.04E+00  0.00E+00 -4.46E+00  1.47E+02
 
 TH11
+       -1.06E+01 -2.82E+00 -5.66E+01 -9.49E+00  3.64E+01  3.45E+00  2.08E+00  0.00E+00  8.01E+00  2.85E+01  1.04E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.041
Stop Time:
Wed Sep 29 12:28:35 CDT 2021
