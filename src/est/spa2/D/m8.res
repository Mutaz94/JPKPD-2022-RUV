Thu Sep 30 08:32:43 CDT 2021
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
$DATA ../../../../data/spa2/D/dat8.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m8.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2074.57909496186        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1535E+02 -3.1333E+01  3.7058E+01  4.8439E+01  1.0401E+02 -2.2183E+02 -1.8861E+02 -2.4029E+01 -1.8426E+02 -2.7509E+01
            -4.6474E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2200.49246414837        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      170
 NPARAMETR:  1.1062E+00  1.2337E+00  8.2646E-01  8.9619E-01  1.0526E+00  1.1043E+00  2.4637E+00  1.1199E+00  1.4459E+00  9.3974E-01
             9.6910E-01
 PARAMETER:  2.0092E-01  3.1006E-01 -9.0605E-02 -9.6047E-03  1.5123E-01  1.9926E-01  1.0017E+00  2.1322E-01  4.6874E-01  3.7845E-02
             6.8612E-02
 GRADIENT:   1.7502E+02  7.6874E+00 -3.5032E+01  3.8075E+01  3.3894E+01 -2.5326E+02  5.0037E+01 -3.8449E+00  2.9142E+01  9.3477E+00
            -6.8175E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2209.50126928062        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      347
 NPARAMETR:  1.1076E+00  1.1641E+00  9.8404E-01  9.4685E-01  1.0373E+00  1.1376E+00  2.6883E+00  1.5604E+00  1.2983E+00  7.6238E-01
             9.7823E-01
 PARAMETER:  2.0217E-01  2.5194E-01  8.3913E-02  4.5386E-02  1.3658E-01  2.2896E-01  1.0889E+00  5.4496E-01  3.6105E-01 -1.7131E-01
             7.7991E-02
 GRADIENT:   1.6826E+02  1.6227E+01 -1.0609E+01  2.3364E+01  2.1906E+01 -2.2723E+02  5.6915E+01  3.3121E+00  1.7203E+01 -7.2094E+00
            -4.8944E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2256.94988117059        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  1.0214E+00  1.1796E+00  9.3385E-01  8.9132E-01  1.0504E+00  1.4347E+00  2.0690E+00  1.3039E+00  1.2299E+00  8.1052E-01
             1.0275E+00
 PARAMETER:  1.2122E-01  2.6518E-01  3.1565E-02 -1.5047E-02  1.4917E-01  4.6095E-01  8.2706E-01  3.6539E-01  3.0694E-01 -1.1008E-01
             1.2715E-01
 GRADIENT:   1.6252E+01 -1.5912E+00  2.3050E+00  1.2187E+00  1.9273E+00 -5.2185E+01 -2.8953E+00 -2.8516E+00 -3.8895E+00 -3.2789E+00
            -4.7771E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2260.95485289014        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  1.0165E+00  1.2171E+00  9.4446E-01  8.7430E-01  1.0711E+00  1.6726E+00  2.0663E+00  1.3555E+00  1.2635E+00  8.5000E-01
             1.0282E+00
 PARAMETER:  1.1639E-01  2.9645E-01  4.2857E-02 -3.4333E-02  1.6871E-01  6.1437E-01  8.2574E-01  4.0418E-01  3.3386E-01 -6.2522E-02
             1.2785E-01
 GRADIENT:   8.5975E+00  2.0344E+00  3.6842E+00 -8.1363E-01 -1.4053E+00  2.3584E+01 -7.3344E-03 -9.7766E-01  1.3077E+00 -1.4155E-01
            -2.1964E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2260.97592415981        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  1.0202E+00  1.2277E+00  9.3646E-01  8.7259E-01  1.0739E+00  1.6592E+00  2.0664E+00  1.3617E+00  1.2616E+00  8.5335E-01
             1.0297E+00
 PARAMETER:  1.1995E-01  3.0510E-01  3.4357E-02 -3.6293E-02  1.7134E-01  6.0631E-01  8.2580E-01  4.0876E-01  3.3241E-01 -5.8583E-02
             1.2929E-01
 GRADIENT:   1.1660E+01  3.8450E+00  7.6518E-01  2.4780E+00  1.0404E+00  1.9784E+01  7.5790E-01 -1.4193E-01  1.0240E+00 -1.6899E-01
            -7.0306E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2261.03705590924        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1072             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0172E+00  1.2112E+00  9.3465E-01  8.7691E-01  1.0708E+00  1.6562E+00  2.0762E+00  1.3576E+00  1.2487E+00  8.5366E-01
             1.0302E+00
 PARAMETER:  1.1706E-01  2.9162E-01  3.2412E-02 -3.1348E-02  1.6838E-01  6.0453E-01  8.3052E-01  4.0572E-01  3.2212E-01 -5.8227E-02
             1.2972E-01
 GRADIENT:   4.8935E+02  2.0506E+02  2.5209E+00  6.7330E+01  3.9689E+01  6.5592E+02  2.6054E+02  2.2988E+00  2.8350E+01  1.6123E+00
             1.6807E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2261.09343681987        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1232
 NPARAMETR:  1.0098E+00  1.2113E+00  9.3648E-01  8.7539E-01  1.0683E+00  1.6371E+00  2.0811E+00  1.3599E+00  1.2476E+00  8.5207E-01
             1.0301E+00
 PARAMETER:  1.0972E-01  2.9168E-01  3.4376E-02 -3.3085E-02  1.6603E-01  5.9295E-01  8.3290E-01  4.0743E-01  3.2119E-01 -6.0088E-02
             1.2965E-01
 GRADIENT:   3.1011E+00  7.4212E-01  6.0721E-01 -9.2259E-01 -6.0221E-01  1.4463E+01  1.3410E+00 -1.0051E-01 -3.6652E-01  2.3882E-01
             1.5704E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -2261.09578127444        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1333
 NPARAMETR:  1.0084E+00  1.2088E+00  9.3595E-01  8.7603E-01  1.0677E+00  1.6390E+00  2.0898E+00  1.3607E+00  1.2490E+00  8.5006E-01
             1.0299E+00
 PARAMETER:  1.0852E-01  2.9003E-01  3.2805E-02 -3.2374E-02  1.6622E-01  5.9370E-01  8.3499E-01  4.0835E-01  3.2314E-01 -6.2101E-02
             1.2957E-01
 GRADIENT:   2.2462E-01  3.5044E-01 -3.3664E-01 -2.5438E-02  9.4258E-01 -3.0367E-01 -7.8832E+03  4.5656E-02  2.0350E+04  3.6124E-02
             1.6698E-01
 NUMSIGDIG:         2.6         2.5         1.7         3.5         2.1         2.9         2.3         2.7         2.3         2.2
                    2.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1333
 NO. OF SIG. DIGITS IN FINAL EST.:  1.7
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5893E-04 -5.7040E-03 -3.3380E-02  8.5533E-03 -3.3886E-02
 SE:             2.9978E-02  2.6500E-02  1.6795E-02  2.4478E-02  1.9359E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9577E-01  8.2958E-01  4.6870E-02  7.2677E-01  8.0043E-02

 ETASHRINKSD(%)  1.0000E-10  1.1223E+01  4.3733E+01  1.7995E+01  3.5146E+01
 ETASHRINKVR(%)  1.0000E-10  2.1186E+01  6.8341E+01  3.2751E+01  5.7939E+01
 EBVSHRINKSD(%)  1.2872E-01  9.7786E+00  4.6402E+01  2.0708E+01  3.3572E+01
 EBVSHRINKVR(%)  2.5728E-01  1.8601E+01  7.1272E+01  3.7127E+01  5.5873E+01
 RELATIVEINF(%)  9.9740E+01  4.3054E+01  1.9363E+01  2.8314E+01  1.7459E+01
 EPSSHRINKSD(%)  3.1523E+01
 EPSSHRINKVR(%)  5.3110E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2261.0957812744387     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1158.3695414288316     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.08
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2261.096       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.21E+00  9.35E-01  8.76E-01  1.07E+00  1.64E+00  2.09E+00  1.36E+00  1.25E+00  8.50E-01  1.03E+00
 


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
+        4.06E+02
 
 TH 2
+       -4.89E-01  1.34E+06
 
 TH 3
+        1.10E+00  3.61E+01  1.91E+02
 
 TH 4
+       -1.28E+00 -5.35E+06 -1.40E+02  5.65E+02
 
 TH 5
+       -1.37E+00 -1.03E+02 -1.79E+02  3.02E+02  6.04E+02
 
 TH 6
+       -5.75E-02 -6.03E-02  1.54E-01 -4.28E-01 -3.72E-01  7.42E+01
 
 TH 7
+       -5.00E+00  3.21E+01 -7.62E+02  3.35E+02  8.16E+01 -7.95E+00  5.43E+04
 
 TH 8
+        1.83E-02 -1.68E+00 -2.58E+01 -1.95E+00 -9.59E+00 -9.74E-03  1.98E+03  1.62E+01
 
 TH 9
+       -3.72E+06 -1.16E+06  3.25E+03  4.65E+06 -3.03E+02  3.37E+01  4.12E+02  7.28E+05  1.01E+06
 
 TH10
+        1.38E-01 -7.10E+00  3.47E+00 -3.78E+00 -7.39E+01 -6.58E-02 -9.84E+00  5.02E+00  4.86E+01  7.21E+01
 
 TH11
+       -2.12E+00 -8.77E+00 -9.94E+00  1.10E+01 -3.45E+01  8.63E-01  4.93E+02  1.67E+01 -3.05E+06  1.96E+01  4.79E+02
 
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
 #CPUT: Total CPU Time in Seconds,       40.709
Stop Time:
Thu Sep 30 08:33:25 CDT 2021
