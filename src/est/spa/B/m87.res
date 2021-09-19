Sat Sep 18 08:47:08 CDT 2021
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
$DATA ../../../../data/spa/B/dat87.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1675.02599035543        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0304E+01 -1.0119E+02 -7.8740E+01 -4.2152E+01  1.1944E+02 -1.1124E+01 -4.5169E+00  8.7329E+00  5.1454E+00  2.9634E+00
            -1.9703E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1685.37616602268        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  1.0006E+00  1.0680E+00  1.1369E+00  9.8195E-01  1.0001E+00  1.0408E+00  9.9106E-01  9.5722E-01  9.8115E-01  9.2953E-01
             1.0787E+00
 PARAMETER:  1.0057E-01  1.6580E-01  2.2834E-01  8.1781E-02  1.0010E-01  1.4003E-01  9.1019E-02  5.6274E-02  8.0973E-02  2.6928E-02
             1.7577E-01
 GRADIENT:   8.1987E+00 -3.3591E+01 -8.1925E+00 -3.8823E+01  2.2744E+01  7.3798E+00  1.6254E-01  4.4794E-01 -3.7547E+00 -5.3375E+00
             9.2491E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1686.13639123991        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  1.0162E+00  1.0243E+00  1.2298E+00  1.0221E+00  1.0132E+00  1.0423E+00  8.9716E-01  9.9619E-01  1.0078E+00  9.5872E-01
             1.0683E+00
 PARAMETER:  1.1609E-01  1.2401E-01  3.0686E-01  1.2189E-01  1.1314E-01  1.4143E-01 -8.5259E-03  9.6187E-02  1.0772E-01  5.7846E-02
             1.6606E-01
 GRADIENT:   4.8045E+01 -2.6187E+01 -7.0805E+00 -2.3724E+01  2.3346E+01  8.4014E+00 -2.9536E+00 -1.2146E+00 -5.5365E-02 -5.2065E+00
             4.6359E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1686.69160466596        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  1.0016E+00  1.0622E+00  1.2916E+00  1.0128E+00  1.0372E+00  1.0295E+00  9.0970E-01  1.1401E+00  1.0123E+00  9.8929E-01
             1.0554E+00
 PARAMETER:  1.0159E-01  1.6030E-01  3.5587E-01  1.1269E-01  1.3648E-01  1.2907E-01  5.3608E-03  2.3109E-01  1.1227E-01  8.9232E-02
             1.5387E-01
 GRADIENT:   1.3328E+01 -7.6051E+00 -3.3460E+00 -6.9699E+00  8.6861E+00  2.9060E+00 -1.2623E+00 -1.8234E-02 -6.6415E-02 -1.6300E+00
             8.5230E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1686.88754238521        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      395
 NPARAMETR:  1.0040E+00  1.0638E+00  1.3253E+00  1.0170E+00  1.0455E+00  1.0314E+00  9.1818E-01  1.1842E+00  1.0105E+00  9.9982E-01
             1.0540E+00
 PARAMETER:  1.0404E-01  1.6185E-01  3.8162E-01  1.1684E-01  1.4453E-01  1.3095E-01  1.4633E-02  2.6910E-01  1.1047E-01  9.9820E-02
             1.5256E-01
 GRADIENT:  -2.3388E+01 -8.9274E+00 -2.6167E+00 -9.2160E+00  4.8298E+00 -4.0695E+00 -9.3147E-01  1.6468E-03 -5.7751E-01 -8.8454E-01
             3.7162E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1687.03398487821        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  1.0086E+00  1.0683E+00  1.3501E+00  1.0167E+00  1.0528E+00  1.0353E+00  9.2025E-01  1.2170E+00  1.0122E+00  1.0060E+00
             1.0534E+00
 PARAMETER:  1.0854E-01  1.6604E-01  4.0021E-01  1.1657E-01  1.5146E-01  1.3470E-01  1.6889E-02  2.9635E-01  1.1210E-01  1.0602E-01
             1.5205E-01
 GRADIENT:  -1.3875E+01 -6.5326E+00 -1.6104E+00 -6.9803E+00  3.0418E+00 -2.3939E+00 -5.5998E-01 -7.2784E-02 -2.8370E-01 -6.7692E-01
             1.9144E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1687.05166528094        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  1.0094E+00  1.0683E+00  1.3516E+00  1.0167E+00  1.0528E+00  1.0367E+00  9.3197E-01  1.2197E+00  1.0122E+00  1.0060E+00
             1.0535E+00
 PARAMETER:  1.0932E-01  1.6604E-01  4.0126E-01  1.1657E-01  1.5146E-01  1.3603E-01  2.9544E-02  2.9863E-01  1.1210E-01  1.0602E-01
             1.5208E-01
 GRADIENT:   3.3294E+01  2.7848E-01 -9.4246E-01  9.2785E-01  3.3075E+00  6.3480E+00  3.8238E-01  2.7383E-02  1.2109E+00 -4.9452E-01
             4.2056E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1687.05380115956        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  1.0094E+00  1.0683E+00  1.3528E+00  1.0167E+00  1.0528E+00  1.0372E+00  9.3040E-01  1.2197E+00  1.0122E+00  1.0060E+00
             1.0533E+00
 PARAMETER:  1.0932E-01  1.6604E-01  4.0216E-01  1.1657E-01  1.5146E-01  1.3657E-01  2.7863E-02  2.9863E-01  1.1210E-01  1.0602E-01
             1.5195E-01
 GRADIENT:   3.3329E+01  4.2566E-01 -7.0302E-01  8.3394E-01  2.9307E+00  6.5976E+00  3.0529E-01 -3.8294E-02  1.1115E+00 -5.2765E-01
             3.3484E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1687.05380115956        NO. OF FUNC. EVALS.:  60
 CUMULATIVE NO. OF FUNC. EVALS.:      960
 NPARAMETR:  1.0094E+00  1.0683E+00  1.3528E+00  1.0167E+00  1.0528E+00  1.0372E+00  9.3040E-01  1.2197E+00  1.0122E+00  1.0060E+00
             1.0533E+00
 PARAMETER:  1.0932E-01  1.6604E-01  4.0216E-01  1.1657E-01  1.5146E-01  1.3657E-01  2.7863E-02  2.9863E-01  1.1210E-01  1.0602E-01
             1.5195E-01
 GRADIENT:  -9.4123E+00 -1.5989E+01  5.4189E+05 -2.6075E+01  1.9025E+01 -7.9780E+05  1.8795E-03 -3.6490E+05 -6.4471E+00 -5.1930E+00
            -1.4345E+06

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      960
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8351E-03 -1.8128E-02 -2.8259E-02  1.0297E-02 -3.6476E-02
 SE:             2.9947E-02  1.6318E-02  1.3248E-02  2.5213E-02  2.2328E-02
 N:                     100         100         100         100         100

 P VAL.:         8.4551E-01  2.6660E-01  3.2924E-02  6.8298E-01  1.0233E-01

 ETASHRINKSD(%)  1.0000E-10  4.5332E+01  5.5616E+01  1.5534E+01  2.5199E+01
 ETASHRINKVR(%)  1.0000E-10  7.0114E+01  8.0300E+01  2.8655E+01  4.4048E+01
 EBVSHRINKSD(%)  4.4613E-01  4.4802E+01  5.9235E+01  1.6396E+01  2.2852E+01
 EBVSHRINKVR(%)  8.9027E-01  6.9531E+01  8.3382E+01  3.0104E+01  4.0481E+01
 RELATIVEINF(%)  9.8378E+01  6.9877E-01  2.1713E+00  1.8265E+00  1.0940E+01
 EPSSHRINKSD(%)  4.4095E+01
 EPSSHRINKVR(%)  6.8746E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1687.0538011595556     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -951.90297459581745     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.34
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1687.054       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.07E+00  1.35E+00  1.02E+00  1.05E+00  1.04E+00  9.30E-01  1.22E+00  1.01E+00  1.01E+00  1.05E+00
 


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
+        8.95E+09
 
 TH 2
+       -1.01E+01  3.46E+09
 
 TH 3
+        9.08E+08  5.67E+01  1.84E+08
 
 TH 4
+       -7.54E+00 -2.59E+09 -4.52E+01  7.27E+02
 
 TH 5
+       -3.10E+09 -1.83E+02 -1.55E+02  4.93E+01  4.29E+09
 
 TH 6
+        1.02E+01 -3.01E+00  2.06E+03  1.71E+00 -1.68E+00  2.71E+09
 
 TH 7
+        3.45E+00  5.57E+00  4.03E+02  6.29E-01 -7.98E+00 -4.13E+09  2.78E+01
 
 TH 8
+       -1.36E+09 -1.55E+01 -2.75E+08  3.07E+01  9.38E+08 -3.07E+03 -5.96E+02  8.21E+08
 
 TH 9
+        4.13E+00 -1.45E+01 -7.46E+00 -4.05E+09 -3.01E+09  3.16E+00  2.87E+01  3.63E+01  8.46E+09
 
 TH10
+       -4.63E+09  5.22E+00 -6.63E+00 -2.75E+00  3.20E+09  5.72E-01  3.13E+00  1.40E+09 -4.50E+09  9.58E+09
 
 TH11
+       -3.09E+09 -7.82E+00 -5.86E+04 -1.32E+01 -9.33E-01 -6.99E+03 -1.35E+03  9.35E+08 -7.26E+00  1.40E+01  2.13E+09
 
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
 #CPUT: Total CPU Time in Seconds,       18.996
Stop Time:
Sat Sep 18 08:47:29 CDT 2021
