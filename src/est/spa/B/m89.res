Sat Sep 18 08:47:47 CDT 2021
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
$DATA ../../../../data/spa/B/dat89.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1642.54745696230        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9796E+01 -7.0873E+01 -1.9005E+01 -6.0693E+01  6.6034E+01  3.7462E+01 -5.3488E+00 -2.5888E+00  9.7898E+00 -1.4927E+01
            -2.1300E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1649.05559149980        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9332E-01  1.0439E+00  9.6976E-01  1.0289E+00  9.6028E-01  8.9251E-01  1.0179E+00  1.0211E+00  9.3619E-01  1.0339E+00
             1.0497E+00
 PARAMETER:  9.3300E-02  1.4297E-01  6.9294E-02  1.2853E-01  5.9470E-02 -1.3719E-02  1.1778E-01  1.2083E-01  3.4067E-02  1.3338E-01
             1.4848E-01
 GRADIENT:   2.6964E+01  4.0786E+00 -5.5794E+00  1.6219E+01  1.1276E+01 -3.7185E+00 -2.9213E+00  1.7648E-01  8.4856E-01  6.7642E-01
            -1.9705E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1649.53839307160        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.9176E-01  9.3710E-01  8.8972E-01  1.0958E+00  8.6104E-01  9.0627E-01  1.2325E+00  9.0164E-01  8.2177E-01  9.0873E-01
             1.0647E+00
 PARAMETER:  9.1731E-02  3.5033E-02 -1.6847E-02  1.9152E-01 -4.9609E-02  1.5861E-03  3.0907E-01 -3.5354E-03 -9.6301E-02  4.2887E-03
             1.6270E-01
 GRADIENT:   2.0377E+01  1.7395E+01 -2.1131E+00  3.2336E+01  2.5922E+00  2.0678E+00 -7.4541E-01  1.7018E+00 -6.4347E+00 -7.4575E-01
             4.6997E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1650.02465279806        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.8559E-01  9.2087E-01  6.3758E-01  1.0642E+00  7.1562E-01  9.0480E-01  1.2691E+00  5.2779E-01  8.2397E-01  7.6640E-01
             1.0480E+00
 PARAMETER:  8.5484E-02  1.7561E-02 -3.5008E-01  1.6225E-01 -2.3461E-01 -3.8341E-05  3.3827E-01 -5.3905E-01 -9.3619E-02 -1.6605E-01
             1.4686E-01
 GRADIENT:  -2.3443E+00  4.6997E+00 -1.1057E+01  1.5157E+01  1.0646E+01  1.5480E-01  1.9045E+00  2.6335E+00 -2.3995E+00  3.2500E+00
            -1.3099E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1650.18469177105        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.8439E-01  8.2097E-01  5.1743E-01  1.0886E+00  5.9780E-01  9.0421E-01  1.3761E+00  3.1773E-01  7.9169E-01  6.1971E-01
             1.0507E+00
 PARAMETER:  8.4271E-02 -9.7272E-02 -5.5889E-01  1.8492E-01 -4.1450E-01 -6.9468E-04  4.1925E-01 -1.0466E+00 -1.3359E-01 -3.7851E-01
             1.4950E-01
 GRADIENT:  -8.8797E+00  3.6165E+00 -7.9504E+00  9.6654E+00  6.4363E+00 -7.0560E-01  1.2566E+00  1.6132E+00 -7.1933E-01  2.9489E+00
            -2.7004E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1650.18834118322        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:      407
 NPARAMETR:  9.8479E-01  7.9521E-01  5.0013E-01  1.0969E+00  5.7700E-01  9.0441E-01  1.4081E+00  2.8808E-01  7.8328E-01  5.9235E-01
             1.0531E+00
 PARAMETER:  8.4669E-02 -1.2915E-01 -5.9289E-01  1.9250E-01 -4.4991E-01 -4.7521E-04  4.4226E-01 -1.1445E+00 -1.4427E-01 -4.2366E-01
             1.5174E-01
 GRADIENT:  -4.0878E+01  1.8487E+00 -1.0044E+01 -2.1317E+00 -7.5688E-01 -3.0132E+00 -4.9440E-01  1.3704E+00 -1.2746E+00  2.5285E+00
            -2.6916E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1651.65037081220        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.0003E+00  6.8865E-01  6.0703E-01  1.1824E+00  6.0974E-01  9.1045E-01  1.6093E+00  2.8387E-01  7.6192E-01  6.8992E-01
             1.0618E+00
 PARAMETER:  1.0034E-01 -2.7302E-01 -3.9917E-01  2.6758E-01 -3.9472E-01  6.1795E-03  5.7580E-01 -1.1592E+00 -1.7191E-01 -2.7119E-01
             1.5996E-01
 GRADIENT:   4.5572E+00  2.5520E+00  1.7897E+00  2.4614E+00 -1.6556E+00  4.0166E-01  1.9562E-01 -1.3657E-01 -7.9702E-01 -7.0715E-01
             7.6533E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1651.67780949051        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  9.9877E-01  6.6804E-01  6.1109E-01  1.1921E+00  6.0704E-01  9.0989E-01  1.6434E+00  2.8985E-01  7.5886E-01  6.9807E-01
             1.0605E+00
 PARAMETER:  9.8767E-02 -3.0341E-01 -3.9252E-01  2.7572E-01 -3.9916E-01  5.5691E-03  5.9678E-01 -1.1384E+00 -1.7594E-01 -2.5943E-01
             1.5876E-01
 GRADIENT:   9.0976E-01 -4.2292E-01 -1.4596E+00 -7.7404E-01  2.3703E+00  1.8587E-01  3.4743E-01 -9.2807E-02 -1.1739E-01  2.9479E-01
             5.5106E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1651.81682518318        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      940
 NPARAMETR:  9.9546E-01  5.6847E-01  7.0641E-01  1.2682E+00  6.2851E-01  9.0758E-01  1.8538E+00  5.1943E-01  7.3467E-01  7.3944E-01
             1.0549E+00
 PARAMETER:  9.5453E-02 -4.6481E-01 -2.4756E-01  3.3757E-01 -3.6440E-01  3.0276E-03  7.1725E-01 -5.5502E-01 -2.0833E-01 -2.0186E-01
             1.5340E-01
 GRADIENT:  -1.5409E+00  3.1735E+00  3.6382E+00  2.9214E+00 -8.3663E+00 -2.7070E-01 -1.8215E-03  1.3089E-01  4.6947E-01  3.6127E-01
            -7.1138E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1652.15015483247        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1118
 NPARAMETR:  9.9079E-01  4.2549E-01  9.8900E-01  1.3861E+00  7.3762E-01  9.0508E-01  2.2589E+00  8.7126E-01  6.9839E-01  8.5844E-01
             1.0527E+00
 PARAMETER:  9.0749E-02 -7.5453E-01  8.8943E-02  4.2650E-01 -2.0433E-01  2.6836E-04  9.1487E-01 -3.7820E-02 -2.5897E-01 -5.2633E-02
             1.5138E-01
 GRADIENT:   6.7053E-01  1.9280E+00  2.7984E+00  1.6057E+00 -1.9221E+00  1.6835E-01  8.1138E-02 -5.3158E-01 -5.1782E-01 -6.9609E-01
             2.3077E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1652.31459616344        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1293
 NPARAMETR:  9.8741E-01  2.9973E-01  1.1200E+00  1.4679E+00  7.6343E-01  9.0281E-01  2.7705E+00  1.0287E+00  6.7860E-01  8.9334E-01
             1.0500E+00
 PARAMETER:  8.7334E-02 -1.1049E+00  2.1334E-01  4.8381E-01 -1.6994E-01 -2.2442E-03  1.1190E+00  1.2833E-01 -2.8773E-01 -1.2788E-02
             1.4878E-01
 GRADIENT:  -4.1048E-01 -1.6953E-01  9.0252E-01 -1.6922E+00 -3.8076E-01 -1.0293E-01 -2.9042E-01 -1.6742E-02  7.4618E-02 -6.3277E-01
            -4.0531E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1652.32325994086        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1468
 NPARAMETR:  9.8781E-01  3.0078E-01  1.0759E+00  1.4635E+00  7.4572E-01  9.0329E-01  2.7750E+00  9.7972E-01  6.7910E-01  8.8404E-01
             1.0510E+00
 PARAMETER:  8.7731E-02 -1.1014E+00  1.7320E-01  4.8080E-01 -1.9341E-01 -1.7135E-03  1.1207E+00  7.9508E-02 -2.8698E-01 -2.3249E-02
             1.4977E-01
 GRADIENT:   5.5447E-03 -4.6388E-03 -2.7398E-02 -3.3765E-02  3.9516E-02  9.9770E-04  6.7080E-03  3.8879E-03  6.4344E-03  4.3014E-03
             2.3782E-03

0ITERATION NO.:   56    OBJECTIVE VALUE:  -1652.32325994086        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1490
 NPARAMETR:  9.8781E-01  3.0078E-01  1.0759E+00  1.4635E+00  7.4572E-01  9.0329E-01  2.7750E+00  9.7972E-01  6.7910E-01  8.8404E-01
             1.0510E+00
 PARAMETER:  8.7731E-02 -1.1014E+00  1.7320E-01  4.8080E-01 -1.9341E-01 -1.7135E-03  1.1207E+00  7.9508E-02 -2.8698E-01 -2.3249E-02
             1.4977E-01
 GRADIENT:   5.5447E-03 -4.6388E-03 -2.7398E-02 -3.3765E-02  3.9516E-02  9.9770E-04  6.7080E-03  3.8879E-03  6.4344E-03  4.3014E-03
             2.3782E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1490
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0333E-03  3.1593E-02 -3.0865E-02 -2.6589E-02 -1.8674E-02
 SE:             2.9813E-02  1.7389E-02  1.6679E-02  2.4583E-02  2.0383E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7235E-01  6.9237E-02  6.4247E-02  2.7943E-01  3.5958E-01

 ETASHRINKSD(%)  1.2337E-01  4.1745E+01  4.4122E+01  1.7643E+01  3.1714E+01
 ETASHRINKVR(%)  2.4658E-01  6.6064E+01  6.8777E+01  3.2174E+01  5.3370E+01
 EBVSHRINKSD(%)  5.6991E-01  4.7937E+01  4.5920E+01  1.4807E+01  2.7984E+01
 EBVSHRINKVR(%)  1.1366E+00  7.2895E+01  7.0754E+01  2.7421E+01  4.8137E+01
 RELATIVEINF(%)  9.7994E+01  4.3208E+00  4.2021E+00  1.3279E+01  6.8327E+00
 EPSSHRINKSD(%)  4.4781E+01
 EPSSHRINKVR(%)  6.9509E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1652.3232599408639     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -917.17243337712569     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.78
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1652.323       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  3.01E-01  1.08E+00  1.46E+00  7.46E-01  9.03E-01  2.78E+00  9.80E-01  6.79E-01  8.84E-01  1.05E+00
 


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
+        1.38E+03
 
 TH 2
+       -2.59E+01  5.12E+02
 
 TH 3
+        6.31E+00  8.29E+01  3.10E+02
 
 TH 4
+       -8.69E+00  4.23E+02 -9.06E+01  8.45E+02
 
 TH 5
+        4.55E+00 -2.94E+02 -5.60E+02  3.46E+01  1.35E+03
 
 TH 6
+        3.57E+00 -4.86E+00  2.86E+00 -2.34E+00  2.56E+00  2.38E+02
 
 TH 7
+        5.54E-01  3.53E+01 -1.11E+00 -6.50E+00  2.36E+00 -3.61E-01  7.63E+00
 
 TH 8
+       -8.86E-01 -2.17E+00 -5.71E+01 -8.28E+00  5.95E+00  9.66E+00  4.86E-01  3.79E+01
 
 TH 9
+        3.44E+00 -3.39E+01  3.45E+00 -1.05E+01  7.64E+00 -2.24E+00  1.72E+00 -1.66E+00  2.73E+02
 
 TH10
+        1.72E+00  1.45E+01 -1.40E+01 -2.53E+01 -7.72E+01 -1.18E+00  2.00E+00  2.29E+01 -2.52E+00  7.14E+01
 
 TH11
+       -1.18E+01 -3.97E+00 -6.36E+00 -8.54E+00 -7.10E+00  5.45E+00  2.59E-01  9.05E+00  1.22E+01  1.34E+01  1.93E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.267
Stop Time:
Sat Sep 18 08:48:12 CDT 2021
