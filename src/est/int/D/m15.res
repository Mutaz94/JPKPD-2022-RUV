Sat Sep 25 05:26:17 CDT 2021
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
$DATA ../../../../data/int/D/dat15.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3199.99219807202        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.5892E+02 -1.6434E+02 -4.4350E+01 -3.1028E+02  2.1118E+02 -4.8976E+02 -2.5594E+02 -8.3215E+01 -5.9292E+02 -1.7031E+02
            -1.7266E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3624.19688713309        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.2377E+00  1.3510E+00  1.2507E+00  1.2010E+00  1.0464E+00  1.5612E+00  1.8686E+00  1.9064E+00  2.0239E+00  1.5222E+00
             1.1231E+00
 PARAMETER:  3.1325E-01  4.0086E-01  3.2367E-01  2.8317E-01  1.4532E-01  5.4548E-01  7.2518E-01  7.4521E-01  8.0504E-01  5.2015E-01
             2.1613E-01
 GRADIENT:   1.7948E+02  7.3245E+01 -7.9506E+00  8.0443E+01 -4.1719E+01 -2.7089E+01 -1.7245E-01 -7.0717E-01  1.5717E+01  5.0728E+01
             1.1865E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3654.66733573490        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.2245E+00  1.3678E+00  1.8625E+00  1.0552E+00  1.2288E+00  1.6604E+00  2.4359E+00  2.2250E+00  1.9237E+00  1.1742E+00
             1.1082E+00
 PARAMETER:  3.0251E-01  4.1319E-01  7.2192E-01  1.5375E-01  3.0604E-01  6.0707E-01  9.9032E-01  8.9975E-01  7.5424E-01  2.6060E-01
             2.0278E-01
 GRADIENT:   1.6282E+02  4.2760E+01  4.4648E+01  3.0831E+01 -2.7410E+01  1.5902E+01  4.2151E+01 -1.5423E+01  4.9130E+01  2.2570E+01
             6.7290E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3669.80909234121        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0969E+00  1.4199E+00  1.3178E+00  9.5014E-01  1.1990E+00  1.8363E+00  2.3489E+00  2.1164E+00  1.6890E+00  1.0257E+00
             1.0816E+00
 PARAMETER:  1.9247E-01  4.5062E-01  3.7598E-01  4.8853E-02  2.8146E-01  7.0775E-01  9.5393E-01  8.4971E-01  6.2416E-01  1.2537E-01
             1.7847E-01
 GRADIENT:   4.3931E+01  4.6129E+01 -4.3434E+00  1.0332E+01 -4.5698E+01  8.3597E+01  3.3822E+01  1.1518E+01  3.4607E+01 -1.3073E+01
             3.3291E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3675.20855267395        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      391
 NPARAMETR:  1.1112E+00  1.4226E+00  1.4116E+00  9.3585E-01  1.2580E+00  1.8433E+00  2.4842E+00  2.1144E+00  1.5704E+00  1.0984E+00
             1.0802E+00
 PARAMETER:  2.0548E-01  4.5250E-01  4.4472E-01  3.3697E-02  3.2956E-01  7.1158E-01  1.0100E+00  8.4878E-01  5.5131E-01  1.9387E-01
             1.7717E-01
 GRADIENT:  -2.1069E+01 -1.9281E+00  4.9193E+00  5.8038E+00 -2.1541E+01  1.9646E+00 -1.3498E+01  2.2935E+00  1.3519E+01 -5.5406E+00
            -3.9518E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3677.24646431361        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      573
 NPARAMETR:  1.1319E+00  1.4813E+00  1.4221E+00  8.9283E-01  1.3152E+00  1.8361E+00  2.5456E+00  2.0246E+00  1.4982E+00  1.1606E+00
             1.0902E+00
 PARAMETER:  2.2392E-01  4.9293E-01  4.5213E-01 -1.3360E-02  3.7402E-01  7.0764E-01  1.0344E+00  8.0537E-01  5.0428E-01  2.4897E-01
             1.8637E-01
 GRADIENT:  -9.2232E+00  3.6691E-01  1.2840E+01 -3.0582E+00  4.4103E+00  9.1753E-01 -5.9486E+00 -9.8505E+00  7.4248E+00 -3.7062E+00
             1.0083E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3678.07836010247        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      758             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1435E+00  1.4799E+00  1.3778E+00  8.9514E-01  1.3103E+00  1.8328E+00  2.5920E+00  2.0999E+00  1.4503E+00  1.1853E+00
             1.0843E+00
 PARAMETER:  2.3406E-01  4.9198E-01  4.2049E-01 -1.0776E-02  3.7023E-01  7.0585E-01  1.0524E+00  8.4187E-01  4.7178E-01  2.7003E-01
             1.8097E-01
 GRADIENT:   8.6861E+01  4.9352E+01  4.2250E+00  6.6914E+00  1.4892E+01  8.0886E+01  7.7006E+01 -9.4445E-01  6.9762E+00  7.4882E-01
             1.8490E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3678.08669309378        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      924
 NPARAMETR:  1.1463E+00  1.4797E+00  1.3773E+00  8.9477E-01  1.3103E+00  1.8323E+00  2.5917E+00  2.1009E+00  1.4501E+00  1.1878E+00
             1.0838E+00
 PARAMETER:  2.3658E-01  4.9181E-01  4.2013E-01 -1.1189E-02  3.7022E-01  7.0555E-01  1.0523E+00  8.4237E-01  4.7160E-01  2.7214E-01
             1.8045E-01
 GRADIENT:  -8.7173E-01  2.9195E-01  2.6062E+00  1.3140E-01  1.4070E+00  2.0223E-01 -1.0079E+00 -2.2942E+00  1.5847E+00 -2.2399E-01
             5.6300E-01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -3678.08669309378        NO. OF FUNC. EVALS.:  25
 CUMULATIVE NO. OF FUNC. EVALS.:      949
 NPARAMETR:  1.1463E+00  1.4797E+00  1.3773E+00  8.9477E-01  1.3103E+00  1.8323E+00  2.5917E+00  2.1009E+00  1.4501E+00  1.1878E+00
             1.0838E+00
 PARAMETER:  2.3658E-01  4.9181E-01  4.2013E-01 -1.1189E-02  3.7022E-01  7.0555E-01  1.0523E+00  8.4237E-01  4.7160E-01  2.7214E-01
             1.8045E-01
 GRADIENT:   1.3035E+05  6.2329E+04 -7.3846E+04  9.1742E+02  8.2778E+04 -4.9453E-02  2.9303E+04  3.6794E+04 -6.5788E+04  6.7118E+02
            -1.0383E+03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      949
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.3100E-04 -1.8450E-02 -3.3460E-02  2.8710E-02 -5.0622E-02
 SE:             2.9964E-02  2.8183E-02  1.9722E-02  2.5209E-02  2.2995E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7521E-01  5.1270E-01  8.9774E-02  2.5475E-01  2.7702E-02

 ETASHRINKSD(%)  1.0000E-10  5.5842E+00  3.3930E+01  1.5547E+01  2.2965E+01
 ETASHRINKVR(%)  1.0000E-10  1.0857E+01  5.6347E+01  2.8677E+01  4.0656E+01
 EBVSHRINKSD(%)  9.2639E-02  6.3241E+00  3.5985E+01  1.9704E+01  2.1652E+01
 EBVSHRINKVR(%)  1.8519E-01  1.2248E+01  5.9021E+01  3.5525E+01  3.8615E+01
 RELATIVEINF(%)  9.9814E+01  7.4161E+01  3.5592E+01  4.9360E+01  4.5104E+01
 EPSSHRINKSD(%)  2.2762E+01
 EPSSHRINKVR(%)  4.0342E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3678.0866930937764     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2023.9973333253656     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.02
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    18.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3678.087       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.15E+00  1.48E+00  1.38E+00  8.95E-01  1.31E+00  1.83E+00  2.59E+00  2.10E+00  1.45E+00  1.19E+00  1.08E+00
 


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
+        1.05E+08
 
 TH 2
+       -1.45E+02  1.46E+07
 
 TH 3
+       -1.02E+02  6.44E+03  2.32E+07
 
 TH 4
+        3.18E+08 -7.15E+05 -2.16E+03  9.69E+08
 
 TH 5
+       -9.79E+00 -3.85E+02  1.37E+04 -8.34E+03  3.30E+07
 
 TH 6
+        1.69E+02  5.51E+01 -8.09E+01  9.42E+01  9.14E+01  5.94E+01
 
 TH 7
+        5.98E+04 -5.56E+02 -4.88E+06  2.67E+03  5.83E+06  2.27E+01  1.04E+06
 
 TH 8
+        2.90E+01 -1.52E+04  1.69E+04  3.39E+02 -2.49E+04  2.63E+01 -4.00E+03  2.47E+06
 
 TH 9
+       -1.50E+02 -1.51E+02  1.95E+07 -6.99E+03  3.49E+03 -7.10E+01 -4.13E+06  1.56E+04  1.66E+07
 
 TH10
+       -4.58E+03 -6.18E+03  1.51E+03 -2.65E+08 -4.45E+03  3.92E+01  5.16E+04 -5.65E+02  1.76E+02  4.42E+05
 
 TH11
+        1.45E+08  2.18E+02  1.39E+04  4.41E+08  4.02E+03 -1.29E+01 -1.52E+03 -3.92E+03  1.66E+04 -1.22E+08  2.03E+08
 
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
 #CPUT: Total CPU Time in Seconds,       47.260
Stop Time:
Sat Sep 25 05:27:06 CDT 2021
