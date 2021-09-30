Wed Sep 29 13:29:49 CDT 2021
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
$DATA ../../../../data/spa/A3/dat38.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1160.40316178538        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5944E+02  1.5106E+02  5.0307E+01  1.6564E+02  2.2335E+02  3.1032E+01 -8.7598E+01 -3.8960E+01 -1.6550E+02 -1.9133E+02
            -5.0912E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -791.258940527781        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  1.2281E+00  1.0072E+00  9.3658E-01  1.1537E+00  7.1996E-01  6.8102E-01  1.1047E+00  1.0322E+00  1.1730E+00  1.2529E+00
             1.5728E+01
 PARAMETER:  3.0550E-01  1.0713E-01  3.4475E-02  2.4294E-01 -2.2856E-01 -2.8416E-01  1.9961E-01  1.3165E-01  2.5953E-01  3.2542E-01
             2.8554E+00
 GRADIENT:  -1.0367E+02 -3.9085E+01 -6.9179E+00 -5.7519E+01 -2.4349E-01  1.2278E+00  9.3614E+00  2.8559E+00  2.2573E+01  2.0384E+01
             4.9085E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -885.214317035054        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0957E+00  1.7188E-01  8.7956E-02  2.0762E+00  1.0199E-01  6.6164E-01  4.3086E+00  7.5610E-01  6.8676E-01  7.2633E-01
             1.0165E+01
 PARAMETER:  1.9140E-01 -1.6610E+00 -2.3309E+00  8.3054E-01 -2.1829E+00 -3.1304E-01  1.5606E+00 -1.7958E-01 -2.7577E-01 -2.1976E-01
             2.4190E+00
 GRADIENT:   1.1274E+02  4.2379E+01  1.3558E+02  5.8640E+01 -2.4927E+02 -1.1135E+01  2.5386E+01  5.6630E-01 -4.6525E+00  5.5871E+00
             1.1507E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1082.45796310198        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  8.9092E-01  1.9717E-01  3.1197E-02  1.7382E+00  1.0926E-01  1.0551E+00  3.7010E-01  1.0000E-02  1.0242E+00  1.2142E-01
             5.7482E+00
 PARAMETER: -1.5504E-02 -1.5237E+00 -3.3674E+00  6.5288E-01 -2.1140E+00  1.5364E-01 -8.9397E-01 -6.9242E+00  1.2391E-01 -2.0085E+00
             1.8489E+00
 GRADIENT:  -1.3219E+01 -3.3055E+01 -1.4878E+01  9.4868E+00  9.6685E+01 -1.3782E+01  1.1156E+00  0.0000E+00 -1.0362E+01 -2.8678E+00
            -5.4216E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1093.33966684417        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      413
 NPARAMETR:  9.2812E-01  2.5205E-01  3.7629E-02  1.6446E+00  1.2291E-01  1.1054E+00  2.9924E-01  1.0000E-02  1.2909E+00  1.4179E-01
             6.0680E+00
 PARAMETER:  2.5404E-02 -1.2781E+00 -3.1800E+00  5.9749E-01 -1.9963E+00  2.0024E-01 -1.1065E+00 -7.2041E+00  3.5534E-01 -1.8534E+00
             1.9030E+00
 GRADIENT:   2.2085E+01 -1.9396E+01 -1.3339E+00  7.0548E+00  3.2764E+01  4.3947E+00  1.0025E+00  0.0000E+00 -7.3619E+00 -2.0390E+00
            -3.4401E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1117.68899077490        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      594
 NPARAMETR:  9.1581E-01  1.7285E-01  2.3276E-02  1.8785E+00  8.7187E-02  8.3228E-01  1.5573E-01  1.0000E-02  7.4203E+00  8.0622E-01
             5.1664E+00
 PARAMETER:  1.2053E-02 -1.6553E+00 -3.6603E+00  7.3047E-01 -2.3397E+00 -8.3584E-02 -1.7596E+00 -1.4290E+01  2.1042E+00 -1.1540E-01
             1.7422E+00
 GRADIENT:  -2.2970E+01  1.9413E+01  3.9550E+01 -2.3573E+01 -1.9095E+02 -8.7301E+00  2.9739E-01  0.0000E+00  1.0677E+01  4.5924E+00
             4.3625E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1136.01307655577        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      769
 NPARAMETR:  9.3416E-01  1.8748E-01  2.2721E-02  2.0360E+00  9.8112E-02  9.0185E-01  7.1776E-02  1.0000E-02  6.0189E+00  7.1843E-01
             4.8858E+00
 PARAMETER:  3.1895E-02 -1.5741E+00 -3.6845E+00  8.1099E-01 -2.2216E+00 -3.3017E-03 -2.5342E+00 -1.6427E+01  1.8949E+00 -2.3068E-01
             1.6863E+00
 GRADIENT:   2.3767E+01 -1.0527E+01  1.0399E+01 -1.3195E+01  1.7760E+01 -5.2684E+00  8.6866E-02  0.0000E+00  9.3592E+00  5.4404E-01
            -3.5678E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1138.50710865984        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      945
 NPARAMETR:  9.2029E-01  1.8883E-01  2.1991E-02  2.6971E+00  9.5571E-02  9.3572E-01  3.7299E-02  1.0000E-02  5.6783E+00  7.1782E-01
             4.8960E+00
 PARAMETER:  1.6935E-02 -1.5669E+00 -3.7171E+00  1.0922E+00 -2.2479E+00  3.3561E-02 -3.1888E+00 -1.9317E+01  1.8367E+00 -2.3154E-01
             1.6884E+00
 GRADIENT:   2.6451E+00  9.6135E+00 -1.4054E+01 -2.0243E+00 -2.3633E+01  6.3634E+00  2.4294E-02  0.0000E+00 -8.2493E+00 -3.3182E+00
             9.3701E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1145.05759211166        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1125
 NPARAMETR:  8.3455E-01  2.0212E-01  3.5966E-02  8.0281E+00  9.9338E-02  9.1140E-01  1.0000E-02  1.0000E-02  5.8381E+00  6.2233E-01
             4.8652E+00
 PARAMETER: -8.0867E-02 -1.4989E+00 -3.2252E+00  2.1830E+00 -2.2092E+00  7.2279E-03 -5.7658E+00 -2.8643E+01  1.8644E+00 -3.7429E-01
             1.6821E+00
 GRADIENT:   6.6351E+01  2.0518E+01  3.4401E+00  6.3146E+00 -3.3637E+01 -2.5801E+00  0.0000E+00  0.0000E+00 -7.2946E+00 -2.9571E+00
            -4.0162E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1147.73474313232        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1300
 NPARAMETR:  7.6573E-01  1.9651E-01  4.1086E-02  9.1613E+00  1.0094E-01  9.1373E-01  1.0000E-02  1.0000E-02  6.0070E+00  6.1070E-01
             4.8246E+00
 PARAMETER: -1.6693E-01 -1.5271E+00 -3.0921E+00  2.3150E+00 -2.1932E+00  9.7820E-03 -5.9915E+00 -3.0203E+01  1.8929E+00 -3.9315E-01
             1.6737E+00
 GRADIENT:  -1.1223E+00 -1.5576E+00  1.3533E+00 -9.4443E-01  1.9193E+00 -5.8309E-01  0.0000E+00  0.0000E+00  1.7443E+00 -1.1449E-01
            -5.6566E-01

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1147.75975392389        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:     1426
 NPARAMETR:  7.6458E-01  1.9735E-01  4.1120E-02  9.3483E+00  1.0113E-01  9.1884E-01  1.0000E-02  1.0000E-02  5.9592E+00  6.1160E-01
             4.8306E+00
 PARAMETER: -1.6843E-01 -1.5228E+00 -3.0913E+00  2.3352E+00 -2.1913E+00  1.5360E-02 -6.0624E+00 -3.0476E+01  1.8849E+00 -3.9168E-01
             1.6750E+00
 GRADIENT:  -1.0940E+00 -1.2662E+00 -9.1094E-01  8.0828E-02  6.7100E+00  9.1312E-02  0.0000E+00  0.0000E+00  6.6949E-01  2.5006E-01
             7.0755E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1426
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.2485E-03  6.9511E-05  1.9047E-04 -3.6954E-02  3.6422E-02
 SE:             2.2985E-02  8.8632E-05  1.2837E-04  2.1207E-02  2.1704E-02
 N:                     100         100         100         100         100

 P VAL.:         6.8742E-01  4.3289E-01  1.3788E-01  8.1409E-02  9.3324E-02

 ETASHRINKSD(%)  2.2996E+01  9.9703E+01  9.9570E+01  2.8954E+01  2.7289E+01
 ETASHRINKVR(%)  4.0704E+01  9.9999E+01  9.9998E+01  4.9525E+01  4.7132E+01
 EBVSHRINKSD(%)  2.4009E+01  9.9580E+01  9.9556E+01  2.3910E+01  2.9404E+01
 EBVSHRINKVR(%)  4.2253E+01  9.9998E+01  9.9998E+01  4.2103E+01  5.0163E+01
 RELATIVEINF(%)  1.4614E+01  8.2476E-04  4.6328E-04  2.2579E+01  1.6938E+01
 EPSSHRINKSD(%)  2.1095E+01
 EPSSHRINKVR(%)  3.7740E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1147.7597539238870     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -412.60892736014887     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.38
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
 





 #OBJV:********************************************    -1147.760       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         7.65E-01  1.97E-01  4.11E-02  9.35E+00  1.01E-01  9.19E-01  1.00E-02  1.00E-02  5.96E+00  6.12E-01  4.83E+00
 


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
+        1.30E+03
 
 TH 2
+        1.44E+02  5.25E+03
 
 TH 3
+        8.49E+03  1.58E+03  1.27E+05
 
 TH 4
+        5.96E+00  1.29E+00 -5.29E+01  2.20E-01
 
 TH 5
+       -6.09E+04  8.49E+03 -1.45E+04  4.10E+02  1.97E+05
 
 TH 6
+        7.97E+00 -4.37E+00 -2.47E+01 -1.40E-01  7.68E+02  8.14E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.84E+00  8.55E+00  5.80E+01 -2.51E-01  6.72E+02  1.95E+00  0.00E+00  0.00E+00  1.72E+00
 
 TH10
+        3.34E+01 -4.38E+01  8.06E+02 -2.58E-01  3.24E+04  1.94E+01  0.00E+00  0.00E+00  1.53E+00  1.53E+02
 
 TH11
+       -2.32E+01 -1.82E+01 -2.05E+02 -5.41E-02  1.50E+02  4.60E+00  0.00E+00  0.00E+00  6.92E-01  7.20E+00  1.79E+01
 
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
 #CPUT: Total CPU Time in Seconds,       31.105
Stop Time:
Wed Sep 29 13:30:22 CDT 2021
