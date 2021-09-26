Sat Sep 25 08:14:54 CDT 2021
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
$DATA ../../../../data/spa/A1/dat77.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1530.02015743383        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.9429E+00  1.7710E+01  4.1950E+01 -1.5556E+01 -3.6067E+01 -6.5006E+00 -8.9610E+00 -9.7611E+00  1.2936E+01 -8.3030E+00
            -2.5839E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1576.64852610037        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0225E+00  1.0113E+00  9.7207E-01  1.0091E+00  1.0132E+00  1.0159E+00  1.0451E+00  9.9965E-01  8.7508E-01  9.7013E-01
             1.3882E+00
 PARAMETER:  1.2229E-01  1.1120E-01  7.1675E-02  1.0909E-01  1.1312E-01  1.1582E-01  1.4413E-01  9.9653E-02 -3.3440E-02  6.9672E-02
             4.2804E-01
 GRADIENT:   2.0756E+01  1.3434E+01  1.2791E+01  2.9653E+00 -7.1407E+00  1.5444E+00 -7.6776E+00 -2.0367E+00 -5.1633E+00 -8.7979E-01
            -3.4708E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1580.16028341332        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0281E+00  7.5250E-01  5.4916E-01  1.1637E+00  6.1988E-01  1.0234E+00  1.7163E+00  4.9688E-01  6.9194E-01  4.8370E-01
             1.4058E+00
 PARAMETER:  1.2773E-01 -1.8436E-01 -4.9936E-01  2.5163E-01 -3.7823E-01  1.2309E-01  6.4015E-01 -5.9941E-01 -2.6826E-01 -6.2629E-01
             4.4063E-01
 GRADIENT:   2.0846E+01  6.1449E+01 -1.8498E+01  1.4673E+02  3.2559E+01  2.3755E+00  1.5452E+01  1.4313E+00 -1.2083E+01 -8.3785E+00
            -2.3637E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1589.18583254566        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0154E+00  6.8089E-01  2.8728E-01  1.0408E+00  4.1586E-01  1.0262E+00  1.3480E+00  1.7524E-01  7.1144E-01  3.0655E-01
             1.5024E+00
 PARAMETER:  1.1525E-01 -2.8435E-01 -1.1473E+00  1.3997E-01 -7.7741E-01  1.2587E-01  3.9864E-01 -1.6416E+00 -2.4047E-01 -1.0824E+00
             5.0706E-01
 GRADIENT:  -9.3238E+00  3.0386E+01 -1.5026E+01  5.1946E+01  5.3671E+00  2.6781E+00  1.3922E+00  4.4520E-01  1.2311E+00  5.6661E-01
             2.1632E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1591.50971968134        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      325
 NPARAMETR:  1.0189E+00  4.8438E-01  2.6906E-01  1.0847E+00  3.4758E-01  1.0264E+00  1.6276E+00  8.8012E-02  6.8490E-01  3.5798E-01
             1.4019E+00
 PARAMETER:  1.1874E-01 -6.2490E-01 -1.2128E+00  1.8135E-01 -9.5675E-01  1.2605E-01  5.8712E-01 -2.3303E+00 -2.7848E-01 -9.2728E-01
             4.3781E-01
 GRADIENT:  -2.2331E+01  5.0037E+00 -3.0196E+00 -6.7093E+00 -2.0935E+01 -7.0488E-01 -1.0010E+00  1.5438E-03 -7.5216E+00 -2.3737E+00
            -2.2025E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1594.09186531093        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      502
 NPARAMETR:  1.0319E+00  4.3234E-01  3.6364E-01  1.1797E+00  4.0217E-01  1.0207E+00  1.9364E+00  9.4063E-02  6.9499E-01  4.7548E-01
             1.4072E+00
 PARAMETER:  1.3136E-01 -7.3854E-01 -9.1159E-01  2.6525E-01 -8.1089E-01  1.2050E-01  7.6084E-01 -2.2638E+00 -2.6386E-01 -6.4343E-01
             4.4161E-01
 GRADIENT:   3.3586E+00  3.0189E+00  7.8685E+00  3.4163E+00 -1.0194E+01 -1.2621E+00 -1.0216E+00  4.8038E-02  3.1264E+00 -1.4680E+00
            -1.2269E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1594.19930295653        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      677
 NPARAMETR:  1.0299E+00  3.9911E-01  3.6177E-01  1.1903E+00  3.9426E-01  1.0235E+00  2.0504E+00  7.6759E-02  6.8543E-01  4.9370E-01
             1.4064E+00
 PARAMETER:  1.2942E-01 -8.1852E-01 -9.1674E-01  2.7417E-01 -8.3075E-01  1.2326E-01  8.1802E-01 -2.4671E+00 -2.7772E-01 -6.0583E-01
             4.4103E-01
 GRADIENT:   7.9114E-02 -1.0210E-01 -1.8065E-01 -6.0616E-02  1.0456E-01 -1.1871E-02  4.0436E-02  5.5887E-02  3.3057E-02  7.7182E-02
             3.9618E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1594.22421896373        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      857
 NPARAMETR:  1.0299E+00  4.0053E-01  3.6150E-01  1.1896E+00  3.9442E-01  1.0235E+00  2.0439E+00  2.0283E-02  6.8569E-01  4.9622E-01
             1.4065E+00
 PARAMETER:  1.2947E-01 -8.1497E-01 -9.1748E-01  2.7358E-01 -8.3033E-01  1.2327E-01  8.1484E-01 -3.7980E+00 -2.7733E-01 -6.0073E-01
             4.4110E-01
 GRADIENT:   1.1584E-01 -6.3321E-02 -3.6622E-01 -1.2597E-01  3.9125E-01 -1.7921E-02  4.8075E-02  3.7321E-03  4.6406E-02  1.3989E-01
             4.3175E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1594.22572921539        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      984
 NPARAMETR:  1.0299E+00  4.0128E-01  3.6178E-01  1.1895E+00  3.9477E-01  1.0236E+00  2.0416E+00  1.0000E-02  6.8553E-01  4.9560E-01
             1.4065E+00
 PARAMETER:  1.2942E-01 -8.1310E-01 -9.1671E-01  2.7350E-01 -8.2944E-01  1.2333E-01  8.1375E-01 -4.7286E+00 -2.7757E-01 -6.0199E-01
             4.4113E-01
 GRADIENT:   7.6106E-04  4.6089E-03  1.2053E-02  2.3972E-02 -6.8100E-03  3.0424E-03 -5.4435E-03  0.0000E+00 -2.1553E-02 -1.6056E-03
            -5.5466E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      984
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.5115E-04  4.0575E-02 -6.8673E-04 -2.2724E-02  2.8483E-02
 SE:             2.9768E-02  2.1213E-02  2.7654E-04  2.6355E-02  1.8199E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7451E-01  5.5777E-02  1.3017E-02  3.8856E-01  1.1755E-01

 ETASHRINKSD(%)  2.7328E-01  2.8935E+01  9.9074E+01  1.1706E+01  3.9033E+01
 ETASHRINKVR(%)  5.4581E-01  4.9498E+01  9.9991E+01  2.2042E+01  6.2830E+01
 EBVSHRINKSD(%)  7.3726E-01  2.9818E+01  9.9054E+01  1.1100E+01  3.6558E+01
 EBVSHRINKVR(%)  1.4691E+00  5.0745E+01  9.9991E+01  2.0968E+01  5.9751E+01
 RELATIVEINF(%)  9.8371E+01  1.0747E+01  3.3059E-04  2.6023E+01  1.3149E+00
 EPSSHRINKSD(%)  4.1313E+01
 EPSSHRINKVR(%)  6.5558E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1594.2257292153906     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -859.07490265165245     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1594.226       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  4.01E-01  3.62E-01  1.19E+00  3.95E-01  1.02E+00  2.04E+00  1.00E-02  6.86E-01  4.96E-01  1.41E+00
 


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
+        9.85E+02
 
 TH 2
+       -2.33E+01  8.55E+02
 
 TH 3
+        1.43E+01  1.28E+03  8.30E+03
 
 TH 4
+       -2.17E+01  2.93E+02 -1.40E+03  1.34E+03
 
 TH 5
+        2.41E+01 -2.18E+03 -9.48E+03  8.49E+02  1.20E+04
 
 TH 6
+        3.12E-01 -4.25E+00  7.21E+00 -4.80E+00  2.08E+00  1.83E+02
 
 TH 7
+        1.96E+00  5.75E+01 -4.53E+01 -8.67E+00  1.18E+01  5.56E-02  1.75E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.04E+00 -1.11E+01 -3.33E+01 -9.28E+00  7.20E+01 -3.00E-01  6.56E+00  0.00E+00  2.87E+02
 
 TH10
+       -3.86E+00  1.09E+01 -4.08E+02 -4.74E+01  3.71E+02  9.30E-01  1.19E+01  0.00E+00 -6.69E+00  1.85E+02
 
 TH11
+       -7.58E+00 -8.50E+00 -9.10E+01 -1.49E+01  6.82E+01  2.06E+00  3.05E+00  0.00E+00  1.49E+01  2.99E+01  1.15E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.666
Stop Time:
Sat Sep 25 08:15:13 CDT 2021
