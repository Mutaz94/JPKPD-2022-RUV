Wed Sep 29 19:51:10 CDT 2021
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
$DATA ../../../../data/spa/D/dat23.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   2874.36046474777        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5695E+02  3.9024E+01 -9.7689E+01 -2.2505E+01  2.7046E+02 -8.2666E+02 -3.4118E+02 -1.9401E+01 -6.9262E+02 -3.8907E+02
            -6.6875E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -879.056128338987        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2698E+00  7.9665E-01  1.5140E+00  1.6522E+00  1.5157E+00  2.6621E+00  2.1766E+00  9.6461E-01  4.0653E+00  1.8012E+00
             6.7936E+00
 PARAMETER:  3.3885E-01 -1.2735E-01  5.1473E-01  6.0213E-01  5.1588E-01  1.0791E+00  8.7777E-01  6.3973E-02  1.5025E+00  6.8844E-01
             2.0160E+00
 GRADIENT:   3.9984E+01 -1.7793E+01 -3.0835E+01  3.0131E+01 -2.2992E+01  1.1966E+02  1.2904E+01  3.2166E+00  1.2321E+02  5.8238E+00
             1.8784E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -935.422323343019        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.4396E+00  2.6196E-01  2.3562E+00  2.0993E+00  1.0818E+01  2.3602E+00  9.5641E+00  9.1334E-01  2.4086E+00  1.3095E+01
             5.5533E+00
 PARAMETER:  4.6434E-01 -1.2396E+00  9.5706E-01  8.4161E-01  2.4812E+00  9.5874E-01  2.3580E+00  9.3545E-03  9.7903E-01  2.6722E+00
             1.8144E+00
 GRADIENT:   1.5652E+02  9.5663E+00  6.4313E+00  1.0031E+02 -1.6112E+00  8.1381E+01  1.1349E+01 -1.6963E+00  7.6380E+01  1.6941E+01
             9.0926E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1026.49878443472        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0634E+00  1.2971E-01  6.1830E-01  1.3788E+00  5.6251E+01  1.5484E+00  8.1256E+00  2.7537E+00  1.0259E+00  1.1700E+01
             5.8372E+00
 PARAMETER:  1.6150E-01 -1.9424E+00 -3.8078E-01  4.2120E-01  4.1298E+00  5.3720E-01  2.1950E+00  1.1130E+00  1.2560E-01  2.5596E+00
             1.8642E+00
 GRADIENT:   7.9174E+01  1.1262E+01  1.7941E+01  6.3733E+01  2.6026E-02 -4.1440E+00  1.2977E+01 -6.5670E+00  1.4748E+01 -2.4968E-02
             1.1132E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1062.65845860437        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  7.3441E-01  4.7777E-02  1.0232E-01  7.4824E-01  1.9702E+02  1.0674E+00  4.5418E+00  1.8846E+00  3.5380E-01  1.4909E+01
             4.8180E+00
 PARAMETER: -2.0869E-01 -2.9412E+00 -2.1796E+00 -1.9003E-01  5.3833E+00  1.6522E-01  1.6133E+00  7.3373E-01 -9.3902E-01  2.8020E+00
             1.6724E+00
 GRADIENT:   6.3615E+01  2.0977E+01 -4.6170E+01  2.2300E+02 -2.5049E-02 -1.4520E+02  1.3427E+01 -1.1262E+01 -3.4868E+00  1.6963E-03
            -5.1889E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1099.90365754715        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      478
 NPARAMETR:  7.3814E-01  2.6867E-02  1.2957E-01  7.5784E-01  1.6227E+02  1.3805E+00  3.7205E+00  2.2627E+00  2.6040E-01  1.1779E+01
             4.6033E+00
 PARAMETER: -2.0362E-01 -3.5169E+00 -1.9435E+00 -1.7728E-01  5.1893E+00  4.2247E-01  1.4139E+00  9.1656E-01 -1.2455E+00  2.5663E+00
             1.6268E+00
 GRADIENT:   1.4065E+01  1.3374E-01 -1.6459E+00  7.6932E+00  1.3545E-02 -2.2564E+01  2.5933E-02  1.9621E+01  9.0064E-01 -2.8065E-03
            -1.7645E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1106.67051937116        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      653
 NPARAMETR:  6.2829E-01  1.1094E-02  7.2984E-02  5.6009E-01  4.3473E+02  1.4664E+00  2.4130E+00  1.6172E+00  1.2181E-01  1.0069E+01
             4.6927E+00
 PARAMETER: -3.6476E-01 -4.4013E+00 -2.5175E+00 -4.7966E-01  6.1747E+00  4.8280E-01  9.8089E-01  5.8070E-01 -2.0053E+00  2.4095E+00
             1.6460E+00
 GRADIENT:   3.4693E+00 -4.9415E-02 -2.9884E+00  3.7813E+00  4.2372E-03  1.6675E+00 -2.1814E-03 -5.0222E-01 -4.8211E-02 -1.1182E-04
            -1.4504E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1106.68963654372        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  6.2662E-01  1.3397E-02  7.3614E-02  5.6147E-01  4.2961E+02  1.4583E+00  2.4440E+00  1.6089E+00  1.7624E-01  3.6505E+01
             4.6953E+00
 PARAMETER: -3.6741E-01 -4.2127E+00 -2.5089E+00 -4.7720E-01  6.1629E+00  4.7729E-01  9.9362E-01  5.7555E-01 -1.6359E+00  3.6974E+00
             1.6466E+00
 GRADIENT:  -1.5730E+00 -5.1644E-02  4.6660E-01  6.5740E-01  4.2429E-03 -5.6941E-02 -2.5195E-03 -2.9173E-01 -4.9007E-03 -1.0325E-03
            -6.7228E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1106.72293851706        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1012
 NPARAMETR:  6.2868E-01  1.8421E-02  7.3458E-02  5.5991E-01  4.1860E+02  1.4577E+00  2.4602E+00  1.6128E+00  1.7526E-01  3.7014E+01
             4.6980E+00
 PARAMETER: -3.6414E-01 -3.8943E+00 -2.5110E+00 -4.7997E-01  6.1369E+00  4.7687E-01  1.0002E+00  5.7796E-01 -1.6415E+00  3.7113E+00
             1.6471E+00
 GRADIENT:   6.3873E-01 -4.5361E-02  9.1695E-01 -1.4392E+00  3.3699E-03  1.2552E-01  1.4566E-03  9.4837E-02  6.1937E-03  9.9127E-04
             2.9904E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1106.72785628378        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1194
 NPARAMETR:  6.2892E-01  1.8987E-02  7.3342E-02  5.5969E-01  3.9797E+02  1.4582E+00  2.4499E+00  1.6127E+00  1.7845E-01  3.6444E+01
             4.6972E+00
 PARAMETER: -3.6375E-01 -3.8640E+00 -2.5126E+00 -4.8037E-01  6.0864E+00  4.7718E-01  9.9606E-01  5.7788E-01 -1.6234E+00  3.6958E+00
             1.6470E+00
 GRADIENT:   1.0028E+00 -1.5064E-02  2.1413E-01 -5.7636E-01  2.4590E-03  2.7609E-01  3.8440E-03  2.5392E-01  1.2743E-02  3.2819E-03
             1.2258E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1106.73171968119        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1359
 NPARAMETR:  6.2738E-01  1.9057E-02  7.2997E-02  5.5911E-01  3.7137E+02  1.4556E+00  2.3719E+00  1.6106E+00  1.7477E-01  3.5739E+01
             4.6949E+00
 PARAMETER: -3.6621E-01 -3.8603E+00 -2.5173E+00 -4.8141E-01  6.0172E+00  4.7545E-01  9.6370E-01  5.7660E-01 -1.6443E+00  3.6762E+00
             1.6465E+00
 GRADIENT:   8.6605E+01  1.8816E-01  7.9878E+01  5.2064E+01 -2.5605E-02  2.3329E+01  1.3822E-02  4.7849E+00  4.9651E-01  7.5448E-02
             1.2717E+01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1106.73171968119        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1424
 NPARAMETR:  6.2789E-01  1.9026E-02  7.3405E-02  5.5844E-01  3.7700E+02  1.4570E+00  2.3492E+00  1.6129E+00  1.7549E-01  3.5412E+01
             4.7131E+00
 PARAMETER: -3.6621E-01 -3.8603E+00 -2.5173E+00 -4.8141E-01  6.0172E+00  4.7545E-01  9.6370E-01  5.7660E-01 -1.6443E+00  3.6762E+00
             1.6465E+00
 GRADIENT:  -1.6460E+00  4.1407E-03 -6.3222E+01  3.2036E+02 -2.6421E+01 -7.0792E-01  4.0814E-03 -2.7574E+02 -4.8139E+01  2.1983E+01
            -9.5752E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1424
 NO. OF SIG. DIGITS UNREPORTABLE

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3941E-03 -1.5387E-03 -7.4393E-04 -8.4866E-03 -6.4794E-04
 SE:             2.9511E-02  5.5857E-04  2.6202E-02  5.3082E-03  3.1037E-04
 N:                     100         100         100         100         100

 P VAL.:         8.8163E-01  5.8750E-03  9.7735E-01  1.0987E-01  3.6828E-02

 ETASHRINKSD(%)  1.1335E+00  9.8129E+01  1.2221E+01  8.2217E+01  9.8960E+01
 ETASHRINKVR(%)  2.2542E+00  9.9965E+01  2.2949E+01  9.6838E+01  9.9989E+01
 EBVSHRINKSD(%)  1.6785E+00  9.8092E+01  1.1086E+01  8.2816E+01  9.9436E+01
 EBVSHRINKVR(%)  3.3289E+00  9.9964E+01  2.0942E+01  9.7047E+01  9.9997E+01
 RELATIVEINF(%)  8.8248E+00  6.0289E-03  1.9461E+00  4.4780E-02  6.1268E-04
 EPSSHRINKSD(%)  2.1319E+01
 EPSSHRINKVR(%)  3.8093E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1106.7317196811900     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -371.58089311745186     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.65
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1106.732       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.27E-01  1.91E-02  7.30E-02  5.59E-01  3.71E+02  1.46E+00  2.37E+00  1.61E+00  1.75E-01  3.57E+01  4.69E+00
 


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
+        7.65E+04
 
 TH 2
+        2.35E+05  7.36E+05
 
 TH 3
+       -7.18E+01 -8.36E+02  1.79E+05
 
 TH 4
+       -9.39E+02  2.37E+02 -1.30E+04  3.18E+04
 
 TH 5
+        2.47E-02  3.98E-02 -9.42E-02  5.44E-01  7.94E-04
 
 TH 6
+        5.08E+00  1.55E+01  4.73E+01 -7.45E+00 -1.34E-03  8.79E+01
 
 TH 7
+       -1.06E-02  8.71E-01 -3.97E-01  1.98E-01 -2.91E-05  1.10E-02 -6.33E-04
 
 TH 8
+        6.29E+01  5.81E+01 -2.28E+02 -3.83E+01 -9.37E-03 -2.23E+00 -5.98E-02  4.64E+03
 
 TH 9
+       -5.98E+04 -1.87E+05 -9.23E+02  9.81E+02 -8.61E-02 -9.30E+00 -2.22E-01  1.48E+04  4.76E+04
 
 TH10
+       -4.31E-01 -6.74E-01  1.48E+00 -9.67E+00  5.41E-05  2.39E-02  4.85E-04  1.51E-01  1.51E+00  2.37E-01
 
 TH11
+       -4.25E+00  6.62E+00  6.36E+01  1.72E+02 -6.22E-03  3.17E-01 -9.13E-03  1.08E+00 -2.67E+01  1.10E-01  8.31E+01
 
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
 #CPUT: Total CPU Time in Seconds,       26.576
Stop Time:
Wed Sep 29 19:51:38 CDT 2021
