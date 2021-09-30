Wed Sep 29 21:57:18 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat26.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1812.79386262058        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1457E+02 -4.9003E-01  4.2875E+01  4.6328E+01  4.9685E+01  5.4504E+01 -1.1289E+01 -1.4023E+02  2.6758E+01 -2.0219E+01
            -4.2602E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1917.49980844688        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.7458E-01  1.1387E+00  8.8055E-01  9.2338E-01  9.8920E-01  8.3888E-01  1.0401E+00  1.5079E+00  8.0696E-01  8.9053E-01
             1.6693E+00
 PARAMETER:  7.4247E-02  2.2989E-01 -2.7207E-02  2.0286E-02  8.9145E-02 -7.5690E-02  1.3935E-01  5.1074E-01 -1.1448E-01 -1.5939E-02
             6.1240E-01
 GRADIENT:   6.3409E+01 -8.5481E+00 -3.9580E+00  3.6212E+00  4.9298E+00 -2.7616E+01 -7.0676E+00 -3.4154E+00 -3.9499E+00  7.7871E+00
             7.9945E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1930.81426126631        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      222
 NPARAMETR:  9.6476E-01  1.8128E+00  8.7097E-01  6.8555E-01  1.2007E+00  9.3127E-01  1.0775E+00  2.2786E+00  9.2985E-01  1.0201E+00
             1.5751E+00
 PARAMETER:  6.4128E-02  6.9487E-01 -3.8152E-02 -2.7753E-01  2.8288E-01  2.8791E-02  1.7468E-01  9.2357E-01  2.7273E-02  1.1992E-01
             5.5434E-01
 GRADIENT:  -9.7241E+01  2.1536E+02  1.1324E+01  1.2814E+02 -5.8028E+01 -9.9573E-01  1.2975E+01 -2.8726E+01 -4.8290E+00  9.7359E-01
             9.3875E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1961.45553526618        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      402
 NPARAMETR:  9.8804E-01  1.5280E+00  2.0709E+00  8.1276E-01  1.5125E+00  9.6332E-01  1.0261E+00  3.6878E+00  8.7463E-01  1.4244E+00
             1.4378E+00
 PARAMETER:  8.7970E-02  5.2397E-01  8.2799E-01 -1.0732E-01  5.1376E-01  6.2632E-02  1.2577E-01  1.4050E+00 -3.3952E-02  4.5373E-01
             4.6315E-01
 GRADIENT:  -2.3522E+01  1.0640E+02 -9.2654E+00  1.1201E+02 -1.3103E+01  1.7242E+01  4.2178E-01  7.3227E+00 -4.9746E+00  3.8905E+00
             3.5904E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1962.72741128665        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      566
 NPARAMETR:  9.9683E-01  1.5305E+00  2.0780E+00  8.0956E-01  1.5536E+00  9.0143E-01  1.0163E+00  3.7071E+00  9.6689E-01  1.4287E+00
             1.4317E+00
 PARAMETER:  9.6820E-02  5.2559E-01  8.3142E-01 -1.1127E-01  5.4060E-01 -3.7774E-03  1.1619E-01  1.4102E+00  6.6327E-02  4.5674E-01
             4.5887E-01
 GRADIENT:  -2.6453E+00  9.3376E+01 -1.1192E+01  1.0334E+02  2.0206E+00 -7.4190E+00  6.2308E+00  7.9859E+00  1.4228E+00  3.5775E-01
             3.1695E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1963.04836387594        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      744
 NPARAMETR:  9.9872E-01  1.5280E+00  2.0786E+00  8.0951E-01  1.5531E+00  9.0502E-01  9.4572E-01  3.7047E+00  9.8349E-01  1.4284E+00
             1.4302E+00
 PARAMETER:  9.8715E-02  5.2395E-01  8.3171E-01 -1.1133E-01  5.4023E-01  1.9820E-04  4.4189E-02  1.4096E+00  8.3349E-02  4.5657E-01
             4.5781E-01
 GRADIENT:   2.2270E+00  9.4741E+01 -1.1514E+01  1.0822E+02  1.1177E+00 -5.6994E+00 -3.5398E+00  8.6162E+00 -1.9130E+00  6.6824E-01
             2.9326E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1967.49457974325        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      923
 NPARAMETR:  9.9635E-01  1.4183E+00  2.1079E+00  8.0729E-01  1.5237E+00  9.1823E-01  1.0100E+00  3.5924E+00  8.8780E-01  1.4179E+00
             1.3634E+00
 PARAMETER:  9.6345E-02  4.4944E-01  8.4570E-01 -1.1408E-01  5.2114E-01  1.4689E-02  1.0996E-01  1.3788E+00 -1.9011E-02  4.4918E-01
             4.1000E-01
 GRADIENT:   1.9697E-01 -2.6656E+00 -8.0328E+00  2.9036E+01  2.8547E-02 -4.7790E-02  1.0610E-01  5.6048E+00 -1.2953E-01 -2.3186E-01
            -2.8802E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1967.49497649849        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1098
 NPARAMETR:  9.9643E-01  1.4181E+00  2.1083E+00  8.0727E-01  1.5236E+00  9.1850E-01  1.0056E+00  3.5910E+00  8.9562E-01  1.4179E+00
             1.3633E+00
 PARAMETER:  9.6425E-02  4.4934E-01  8.4589E-01 -1.1409E-01  5.2108E-01  1.4984E-02  1.0555E-01  1.3784E+00 -1.0240E-02  4.4916E-01
             4.0993E-01
 GRADIENT:   4.2088E-01 -3.3711E+00 -8.0654E+00  2.8687E+01 -3.0158E-02  7.4380E-02 -3.0906E-02  5.6837E+00  1.8798E-01 -1.8913E-01
            -2.9821E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1967.58413374195        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1280             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9642E-01  1.4183E+00  2.1245E+00  8.0629E-01  1.5234E+00  9.1872E-01  1.0014E+00  3.5882E+00  9.0329E-01  1.4178E+00
             1.3635E+00
 PARAMETER:  9.6417E-02  4.4944E-01  8.5355E-01 -1.1531E-01  5.2095E-01  1.5230E-02  1.0135E-01  1.3777E+00 -1.7171E-03  4.4911E-01
             4.1006E-01
 GRADIENT:   2.2213E+02  1.5394E+02 -4.2532E+00  5.4944E+01  1.0850E+01  1.6702E+01  3.1784E+00  8.8519E+01  2.2538E+00  2.7422E+00
             1.0686E+00

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1967.71925834527        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     1424
 NPARAMETR:  9.9641E-01  1.4186E+00  2.1231E+00  8.0205E-01  1.5233E+00  9.1867E-01  1.0014E+00  3.5921E+00  9.0122E-01  1.4178E+00
             1.3641E+00
 PARAMETER:  9.6432E-02  4.4982E-01  8.5314E-01 -1.2062E-01  5.2105E-01  1.5195E-02  1.0138E-01  1.3783E+00 -3.0050E-03  4.4927E-01
             4.1059E-01
 GRADIENT:   2.0892E+04  4.6283E+03  2.4473E+03 -1.7305E+04  4.0090E+03  1.3365E-01  1.3332E-02 -7.7937E+02  5.7274E-01  4.6495E+03
             5.0811E+03
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.3         2.5         3.1         2.3         0.8         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1424
 NO. OF SIG. DIGITS IN FINAL EST.:  0.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3822E-03 -7.0548E-03 -5.5574E-02 -4.2579E-03 -7.0663E-02
 SE:             2.9744E-02  2.3032E-02  1.8045E-02  1.8762E-02  1.9700E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6293E-01  7.5938E-01  2.0721E-03  8.2047E-01  3.3458E-04

 ETASHRINKSD(%)  3.5403E-01  2.2839E+01  3.9546E+01  3.7146E+01  3.4004E+01
 ETASHRINKVR(%)  7.0680E-01  4.0461E+01  6.3453E+01  6.0494E+01  5.6445E+01
 EBVSHRINKSD(%)  6.9633E-01  2.2149E+01  5.2446E+01  3.9928E+01  2.6311E+01
 EBVSHRINKVR(%)  1.3878E+00  3.9392E+01  7.7387E+01  6.3913E+01  4.5699E+01
 RELATIVEINF(%)  9.8359E+01  2.6282E+00  7.7522E+00  1.5049E+00  3.1606E+01
 EPSSHRINKSD(%)  3.3559E+01
 EPSSHRINKVR(%)  5.5856E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1967.7192583452679     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1048.7807251405952     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.24
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.00
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1967.719       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  1.42E+00  2.12E+00  8.02E-01  1.52E+00  9.19E-01  1.00E+00  3.59E+00  9.02E-01  1.42E+00  1.36E+00
 


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
+        5.26E+06
 
 TH 2
+        1.80E+01  1.28E+05
 
 TH 3
+        1.02E+01 -1.28E+02  1.60E+04
 
 TH 4
+       -2.05E+02 -8.45E+05 -2.55E+02  5.59E+06
 
 TH 5
+        6.60E+05  1.03E+05 -7.19E+00 -6.80E+05  8.30E+04
 
 TH 6
+        2.02E+02  2.87E+01  1.09E+01 -2.07E+02  2.45E+01  2.28E+02
 
 TH 7
+        5.16E+06  1.59E+02  4.87E+01 -9.00E+02  1.08E+02 -1.33E+00  8.16E+01
 
 TH 8
+       -4.11E+00  4.44E+01 -3.12E+01  9.39E+01 -4.65E+00 -3.83E+00 -1.86E+01  2.25E+03
 
 TH 9
+       -5.81E+06  1.48E+01  8.16E+00 -1.97E+02  2.06E+01 -2.34E-01  3.74E+01 -1.42E+00  4.37E+01
 
 TH10
+        3.29E+01  1.28E+05 -8.78E+00 -8.48E+05  1.03E+05  3.20E+01  1.33E+02  3.83E+00  2.80E+01  1.29E+05
 
 TH11
+        2.33E+01 -4.32E+02 -9.82E+01 -7.81E+02  9.19E+00  3.88E+01  1.60E+02  3.66E+01  3.27E+01 -1.40E+01  1.66E+05
 
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
 #CPUT: Total CPU Time in Seconds,       31.312
Stop Time:
Wed Sep 29 21:57:51 CDT 2021
