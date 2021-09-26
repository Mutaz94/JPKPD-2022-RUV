Sat Sep 25 02:09:55 CDT 2021
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
$DATA ../../../../data/int/SL3/dat32.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      981
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

 TOT. NO. OF OBS RECS:      881
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
 RAW OUTPUT FILE (FILE): m32.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1190.56758775099        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.8313E+01  7.9495E+01 -5.4182E+01 -6.8697E+00  1.3602E+02  3.6691E+01 -1.1487E+02 -1.6178E+02 -1.3048E+02 -4.9907E+01
            -4.8415E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2798.61950617341        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0532E+00  1.0129E+00  1.0501E+00  1.0502E+00  9.6431E-01  8.6354E-01  1.0146E+00  1.0715E+00  1.0055E+00  1.0777E+00
             2.1655E+00
 PARAMETER:  1.5188E-01  1.1286E-01  1.4892E-01  1.4893E-01  6.3653E-02 -4.6711E-02  1.1449E-01  1.6911E-01  1.0544E-01  1.7487E-01
             8.7266E-01
 GRADIENT:   2.4893E+01  1.1340E+01 -7.6349E+00  1.9187E+01  2.7151E+00 -1.4580E+01 -2.9786E+00 -1.6767E+01 -6.4809E+00 -6.1236E+00
            -3.0784E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2809.86784598109        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0527E+00  1.1519E+00  1.9652E+00  1.0265E+00  1.2936E+00  8.9929E-01  8.5720E-01  2.8677E+00  6.4868E-01  1.3869E+00
             2.1879E+00
 PARAMETER:  1.5137E-01  2.4143E-01  7.7560E-01  1.2618E-01  3.5741E-01 -6.1460E-03 -5.4079E-02  1.1535E+00 -3.3282E-01  4.2707E-01
             8.8293E-01
 GRADIENT:   1.8133E+01  5.5035E+01 -1.0226E+01  5.7308E+01  1.0959E+01  1.6678E+00 -3.1584E+01 -1.3386E+01 -5.0300E+01  9.1466E+00
            -2.3787E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2832.12270492780        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0500E+00  9.8908E-01  2.2824E+00  1.1176E+00  1.2062E+00  8.8468E-01  9.0370E-01  3.1276E+00  7.7402E-01  1.1158E+00
             2.4378E+00
 PARAMETER:  1.4879E-01  8.9021E-02  9.2524E-01  2.1118E-01  2.8751E-01 -2.2526E-02 -1.2627E-03  1.2403E+00 -1.5616E-01  2.0956E-01
             9.9109E-01
 GRADIENT:   2.2616E+00  1.4484E+01 -1.5337E+00  7.0367E+00  1.6094E+00 -3.4450E+00 -6.2783E+00  3.0307E+00 -1.7147E+01 -8.9561E+00
             3.5737E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2832.20422819137        NO. OF FUNC. EVALS.: 110
 CUMULATIVE NO. OF FUNC. EVALS.:      341
 NPARAMETR:  1.0499E+00  9.8868E-01  2.2816E+00  1.1177E+00  1.2057E+00  8.8489E-01  9.0105E-01  3.1247E+00  7.7720E-01  1.1169E+00
             2.4368E+00
 PARAMETER:  1.4874E-01  8.8612E-02  9.2490E-01  2.1124E-01  2.8709E-01 -2.2291E-02 -4.1999E-03  1.2393E+00 -1.5206E-01  2.1057E-01
             9.9069E-01
 GRADIENT:  -9.3898E+00  1.3262E+01 -2.4332E+00  2.8137E+00 -3.6583E-01 -4.1703E+00 -6.1203E+00  1.9967E+00 -1.6908E+01 -8.9439E+00
             3.3319E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2832.73310234608        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.0499E+00  9.8868E-01  2.2817E+00  1.1177E+00  1.2057E+00  8.8574E-01  1.0064E+00  3.1239E+00  7.7722E-01  1.1169E+00
             2.4298E+00
 PARAMETER:  1.4872E-01  8.8613E-02  9.2490E-01  2.1125E-01  2.8706E-01 -2.1332E-02  1.0637E-01  1.2391E+00 -1.5203E-01  2.1057E-01
             9.8782E-01
 GRADIENT:   2.7910E+00  1.4244E+01 -1.6317E+00  5.5940E+00  5.6218E-01 -3.1419E+00 -1.2684E+00  3.1190E+00 -8.4534E+00 -7.7563E+00
             3.1345E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2832.98569215466        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  1.0497E+00  9.8872E-01  2.2808E+00  1.1178E+00  1.2053E+00  8.9351E-01  1.0340E+00  3.1207E+00  7.7738E-01  1.1168E+00
             2.3931E+00
 PARAMETER:  1.4852E-01  8.8657E-02  9.2454E-01  2.1138E-01  2.8673E-01 -1.2596E-02  1.3341E-01  1.2381E+00 -1.5183E-01  2.1049E-01
             9.7260E-01
 GRADIENT:   4.1247E+00  1.6194E+01 -1.6882E+00  6.4949E+00  7.9377E-01 -4.9596E-02  2.2535E-02  2.2277E+00 -7.0372E+00 -8.0703E+00
            -1.6096E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2832.99938558652        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      597
 NPARAMETR:  1.0490E+00  9.8819E-01  2.2917E+00  1.1168E+00  1.2043E+00  8.9007E-01  1.0523E+00  3.0797E+00  7.7789E-01  1.1182E+00
             2.3925E+00
 PARAMETER:  1.4782E-01  8.8119E-02  9.2932E-01  2.1047E-01  2.8588E-01 -1.6453E-02  1.5100E-01  1.2248E+00 -1.5117E-01  2.1171E-01
             9.7232E-01
 GRADIENT:   2.2540E+00  1.5735E+01 -2.3450E-01  3.4184E+00 -5.9982E-01 -1.5099E+00  1.1616E+00  3.2295E-01 -5.8846E+00 -7.7837E+00
            -9.5291E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2833.01061031816        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  1.0484E+00  9.8773E-01  2.3012E+00  1.1159E+00  1.2035E+00  8.9354E-01  1.0328E+00  3.0455E+00  7.7831E-01  1.1194E+00
             2.3944E+00
 PARAMETER:  1.4725E-01  8.7651E-02  9.3342E-01  2.0969E-01  2.8520E-01 -1.2562E-02  1.3231E-01  1.2137E+00 -1.5063E-01  2.1276E-01
             9.7315E-01
 GRADIENT:   7.0294E-01  1.4760E+01  1.1107E+00  1.2511E+00 -1.6143E+00 -6.4999E-03  7.0795E-03 -1.3529E+00 -7.1958E+00 -7.7838E+00
             8.1913E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2833.48510070904        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      827            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0528E+00  9.7097E-01  2.3080E+00  1.1173E+00  1.2076E+00  8.9559E-01  1.0349E+00  3.0949E+00  8.1553E-01  1.1833E+00
             2.3959E+00
 PARAMETER:  1.5143E-01  7.0539E-02  9.3640E-01  2.1092E-01  2.8862E-01 -1.0273E-02  1.3432E-01  1.2298E+00 -1.0391E-01  2.6829E-01
             9.7374E-01
 GRADIENT:   1.3071E+01 -2.4428E+00 -9.1655E-01 -9.9458E+00  3.6139E+00  1.1133E+00  3.4023E+00  1.5236E+00  9.2607E-01  1.5228E+00
             6.8467E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2833.74422228379        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      897
 NPARAMETR:  1.0479E+00  9.5694E-01  2.2884E+00  1.1289E+00  1.1864E+00  8.9301E-01  8.5460E-01  3.0620E+00  8.6655E-01  1.1649E+00
             2.3897E+00
 PARAMETER:  1.4678E-01  5.5990E-02  9.2784E-01  2.2122E-01  2.7094E-01 -1.3152E-02 -5.7120E-02  1.2191E+00 -4.3233E-02  2.5262E-01
             9.7118E-01
 GRADIENT:   2.8861E-01  4.3542E-02  1.8919E-01  1.1071E+00 -3.6383E-01  3.1199E-02  2.8345E-02  7.2574E-02  6.2105E-02  3.4254E-02
             2.6868E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2833.78272620285        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1062
 NPARAMETR:  1.0526E+00  9.5844E-01  2.2974E+00  1.1293E+00  1.1891E+00  8.9490E-01  8.6293E-01  3.0705E+00  8.6449E-01  1.1652E+00
             2.3913E+00
 PARAMETER:  1.5123E-01  5.7551E-02  9.3179E-01  2.2164E-01  2.7324E-01 -1.1042E-02 -4.7420E-02  1.2218E+00 -4.5621E-02  2.5287E-01
             9.7184E-01
 GRADIENT:   1.2559E+01  9.5857E-01  2.6543E-01  2.8155E+00  2.9035E-01  8.7912E-01  1.2053E-01  2.7839E-01  2.1099E-01 -3.1588E-02
             1.1646E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2833.78282758502        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1221
 NPARAMETR:  1.0525E+00  9.5844E-01  2.2974E+00  1.1293E+00  1.1891E+00  8.9480E-01  8.6065E-01  3.0705E+00  8.6468E-01  1.1652E+00
             2.3913E+00
 PARAMETER:  1.5115E-01  5.7551E-02  9.3179E-01  2.2164E-01  2.7324E-01 -1.1152E-02 -5.0068E-02  1.2218E+00 -4.5395E-02  2.5287E-01
             9.7184E-01
 GRADIENT:   6.0336E-02 -4.3768E-02 -7.6823E-01 -1.7202E+00 -1.5657E+00 -2.3776E-04  6.8920E-03 -7.3089E-01 -3.8955E-02 -2.8526E-01
            -5.8165E-01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -2833.78283160152        NO. OF FUNC. EVALS.:  60
 CUMULATIVE NO. OF FUNC. EVALS.:     1281
 NPARAMETR:  1.0525E+00  9.5844E-01  2.2974E+00  1.1293E+00  1.1891E+00  8.9480E-01  8.6051E-01  3.0705E+00  8.6473E-01  1.1652E+00
             2.3913E+00
 PARAMETER:  1.5115E-01  5.7551E-02  9.3179E-01  2.2164E-01  2.7324E-01 -1.1153E-02 -5.0233E-02  1.2218E+00 -4.5334E-02  2.5287E-01
             9.7184E-01
 GRADIENT:   5.0147E-02 -5.7643E+03  1.2338E+03 -2.6032E+03  4.2132E+03 -3.0896E-04 -1.1529E+04  9.4292E+02  1.1530E+04  4.5599E+03
            -1.1919E+03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1281
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7128E-03 -2.7236E-02 -2.8923E-02  9.9770E-03 -3.5320E-02
 SE:             2.9494E-02  1.5733E-02  2.1740E-02  2.5504E-02  2.1640E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5369E-01  8.3430E-02  1.8338E-01  6.9566E-01  1.0266E-01

 ETASHRINKSD(%)  1.1926E+00  4.7292E+01  2.7168E+01  1.4557E+01  2.7502E+01
 ETASHRINKVR(%)  2.3710E+00  7.2219E+01  4.6955E+01  2.6995E+01  4.7440E+01
 EBVSHRINKSD(%)  1.5327E+00  4.7707E+01  2.8022E+01  1.5768E+01  2.5883E+01
 EBVSHRINKVR(%)  3.0419E+00  7.2655E+01  4.8192E+01  2.9050E+01  4.5067E+01
 RELATIVEINF(%)  9.6881E+01  5.8571E+00  3.2987E+01  1.8155E+01  2.5660E+01
 EPSSHRINKSD(%)  1.8381E+01
 EPSSHRINKVR(%)  3.3383E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          881
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1619.1696955066332     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2833.7828316015225     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1214.6131360948893     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.30
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.00
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2833.783       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  9.58E-01  2.30E+00  1.13E+00  1.19E+00  8.95E-01  8.61E-01  3.07E+00  8.65E-01  1.17E+00  2.39E+00
 


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
+        1.14E+07
 
 TH 2
+       -2.13E+03  3.14E+07
 
 TH 3
+        9.49E+01  4.69E+02  6.26E+04
 
 TH 4
+       -8.29E+02  1.20E+07 -2.05E+02  4.60E+06
 
 TH 5
+        6.20E+02 -9.24E+06 -1.05E+03 -1.25E+03  2.72E+06
 
 TH 6
+        3.81E+00 -1.82E+03  8.17E+01 -7.05E+02  5.38E+02  2.37E+02
 
 TH 7
+       -2.35E+03 -6.21E+03  2.77E+02 -2.38E+03  1.83E+03 -2.04E+03  3.89E+07
 
 TH 8
+        5.39E+01  2.58E+02 -6.34E+00 -1.11E+02 -5.86E+02  4.64E+01  1.60E+02  2.05E+04
 
 TH 9
+        2.10E+07  2.73E+03 -1.22E+02  1.08E+03 -8.08E+02 -3.73E+07  3.12E+03 -6.98E+01  3.85E+07
 
 TH10
+        6.88E+02 -1.57E+03  6.85E+01 -5.81E+02  4.33E+02  5.92E+02  2.02E+03  4.01E+01 -8.95E+02  3.32E+06
 
 TH11
+       -1.03E+02 -4.38E+02 -4.78E+02  1.63E+02  9.36E+02 -7.21E+01 -2.51E+02 -2.68E+02  1.21E+02 -5.81E+01  5.40E+04
 
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
 #CPUT: Total CPU Time in Seconds,       44.408
Stop Time:
Sat Sep 25 02:10:41 CDT 2021
