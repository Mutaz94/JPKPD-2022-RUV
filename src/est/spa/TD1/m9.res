Wed Sep 29 17:54:10 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat9.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m9.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1632.99472768407        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0953E+02 -2.9032E+01 -2.6146E+01  8.9397E+00  7.1418E+01  2.1301E+01  2.8803E+00  3.0482E+00  9.5927E+00 -5.7934E+00
             5.2544E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1640.70834140596        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  9.2707E-01  9.9310E-01  9.9211E-01  9.9852E-01  9.8166E-01  9.3566E-01  9.9241E-01  9.9656E-01  9.8965E-01  9.9820E-01
             9.2437E-01
 PARAMETER:  2.4274E-02  9.3073E-02  9.2077E-02  9.8518E-02  8.1495E-02  3.3499E-02  9.2379E-02  9.6550E-02  8.9595E-02  9.8201E-02
             2.1352E-02
 GRADIENT:  -2.1057E+01 -5.6281E+01 -2.1305E+01 -5.5594E+01  4.9943E+01 -2.3274E+01 -8.8198E-01  2.6562E+00  1.2294E+00 -5.6304E+00
             2.6731E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1644.73735638297        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      291
 NPARAMETR:  9.3282E-01  9.3210E-01  9.8276E-01  1.0809E+00  9.1122E-01  1.0159E+00  1.0080E+00  8.8457E-01  9.5167E-01  9.9837E-01
             8.6274E-01
 PARAMETER:  3.0462E-02  2.9688E-02  8.2608E-02  1.7775E-01  7.0325E-03  1.1581E-01  1.0796E-01 -2.2657E-02  5.0464E-02  9.8367E-02
            -4.7642E-02
 GRADIENT:  -3.6453E+00  4.3050E+00  1.7233E-01  6.9195E+00  7.4518E-01  9.8834E+00 -2.0036E+00 -3.5175E-01  3.1727E-01 -5.3177E-01
             7.9459E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1645.17473949244        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      468
 NPARAMETR:  9.3512E-01  8.4530E-01  8.9200E-01  1.1276E+00  8.2848E-01  9.9039E-01  1.2415E+00  7.2470E-01  8.7094E-01  9.2021E-01
             8.6500E-01
 PARAMETER:  3.2917E-02 -6.8058E-02 -1.4287E-02  2.2007E-01 -8.8168E-02  9.0342E-02  3.1632E-01 -2.2200E-01 -3.8182E-02  1.6843E-02
            -4.5026E-02
 GRADIENT:   7.2856E-01  5.6184E+00  2.0877E+00  7.4905E+00 -2.6899E+00  2.6203E-02 -2.5129E-01 -6.9094E-01 -4.4711E-01 -7.5333E-01
             2.5444E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1645.27632046196        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      644
 NPARAMETR:  9.3340E-01  6.9765E-01  1.0343E+00  1.2275E+00  8.3964E-01  9.8982E-01  1.3621E+00  8.6753E-01  8.3784E-01  9.6831E-01
             8.5269E-01
 PARAMETER:  3.1078E-02 -2.6004E-01  1.3377E-01  3.0496E-01 -7.4778E-02  8.9765E-02  4.0903E-01 -4.2104E-02 -7.6929E-02  6.7799E-02
            -5.9364E-02
 GRADIENT:   1.9832E+00  6.0946E+00  4.2488E+00  6.0432E+00 -9.0022E+00  6.2115E-01  1.6343E-01  6.2213E-02 -6.0575E-01  7.8336E-01
            -3.2601E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1645.35791191784        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      821
 NPARAMETR:  9.3309E-01  5.7571E-01  1.1546E+00  1.3126E+00  8.5357E-01  9.9033E-01  1.5088E+00  9.7446E-01  8.1070E-01  1.0038E+00
             8.5973E-01
 PARAMETER:  3.0747E-02 -4.5215E-01  2.4373E-01  3.7200E-01 -5.8328E-02  9.0283E-02  5.1134E-01  7.4132E-02 -1.0985E-01  1.0377E-01
            -5.1141E-02
 GRADIENT:   5.8235E+00  7.0646E+00  2.6426E+00  1.2225E+01 -6.7788E+00  1.5956E+00  1.0706E+00  6.6470E-01 -4.1434E-01  1.3349E+00
            -7.9925E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1645.54018851019        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      996
 NPARAMETR:  9.3014E-01  4.2402E-01  1.2569E+00  1.4118E+00  8.5289E-01  9.8603E-01  1.7497E+00  1.0482E+00  7.7914E-01  1.0223E+00
             8.6677E-01
 PARAMETER:  2.7577E-02 -7.5797E-01  3.2865E-01  4.4486E-01 -5.9124E-02  8.5936E-02  6.5944E-01  1.4710E-01 -1.4957E-01  1.2202E-01
            -4.2986E-02
 GRADIENT:   4.1103E+00  5.7186E+00  9.4433E-01  1.5259E+01 -1.8097E+00  7.6720E-01  1.0987E+00 -5.7063E-02 -4.6379E-01  3.1030E-01
             2.8603E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1645.69079114678        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1171
 NPARAMETR:  9.2695E-01  2.9334E-01  1.3447E+00  1.4969E+00  8.5108E-01  9.8211E-01  2.0451E+00  1.1309E+00  7.5429E-01  1.0330E+00
             8.6941E-01
 PARAMETER:  2.4143E-02 -1.1264E+00  3.9620E-01  5.0340E-01 -6.1255E-02  8.1949E-02  8.1546E-01  2.2302E-01 -1.8198E-01  1.3250E-01
            -3.9942E-02
 GRADIENT:   1.2691E+00  4.3030E+00  8.6241E-02  1.7771E+01  3.5247E-01 -6.9944E-02  6.2673E-01 -4.7751E-01 -6.1661E-01 -2.5955E-01
             3.8894E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1645.73869140456        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1346
 NPARAMETR:  9.2456E-01  2.1101E-01  1.4163E+00  1.5515E+00  8.5211E-01  9.7989E-01  2.3109E+00  1.2125E+00  7.4179E-01  1.0423E+00
             8.6561E-01
 PARAMETER:  2.1563E-02 -1.4558E+00  4.4806E-01  5.3920E-01 -6.0037E-02  7.9689E-02  9.3762E-01  2.9270E-01 -1.9869E-01  1.4139E-01
            -4.4324E-02
 GRADIENT:  -1.2882E+00  3.4766E+00  1.3750E+00  1.9270E+01 -2.8100E+00 -5.3445E-01  3.4935E-01 -2.0779E-01  2.8015E-01  3.3982E-01
             2.1986E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1645.79455819335        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1521
 NPARAMETR:  9.2313E-01  1.4336E-01  1.4675E+00  1.5936E+00  8.5031E-01  9.7910E-01  2.6216E+00  1.2762E+00  7.2987E-01  1.0441E+00
             8.6147E-01
 PARAMETER:  2.0014E-02 -1.8424E+00  4.8358E-01  5.6597E-01 -6.2152E-02  7.8883E-02  1.0638E+00  3.4390E-01 -2.1489E-01  1.4315E-01
            -4.9116E-02
 GRADIENT:  -2.0558E+00  2.3507E+00  1.8129E+00  1.5200E+01 -4.3707E+00 -5.0411E-01  1.8770E-01  6.8950E-02  2.0657E-01  5.3029E-01
             3.2373E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1645.86236170925        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1696
 NPARAMETR:  9.2222E-01  8.2810E-02  1.5174E+00  1.6303E+00  8.5120E-01  9.7861E-01  3.0392E+00  1.3368E+00  7.1869E-01  1.0452E+00
             8.5960E-01
 PARAMETER:  1.9033E-02 -2.3912E+00  5.1697E-01  5.8878E-01 -6.1106E-02  7.8374E-02  1.2116E+00  3.9031E-01 -2.3032E-01  1.4422E-01
            -5.1292E-02
 GRADIENT:  -1.8080E+00  1.2235E+00  1.1795E+00  8.1451E+00 -3.5327E+00 -3.8769E-01  9.9910E-02  4.2369E-01  9.1305E-02  4.5066E-01
            -5.0105E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1646.03308623318        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1878
 NPARAMETR:  9.2286E-01  6.6277E-02  1.5184E+00  1.6319E+00  8.5074E-01  9.7933E-01  9.2930E-01  1.3323E+00  7.1943E-01  1.0466E+00
             8.6050E-01
 PARAMETER:  1.9727E-02 -2.6139E+00  5.1766E-01  5.8977E-01 -6.1647E-02  7.9116E-02  2.6676E-02  3.8688E-01 -2.2929E-01  1.4558E-01
            -5.0243E-02
 GRADIENT:   5.2367E-01  2.5335E-01  1.2830E-03 -1.1722E+01  1.0847E+00  3.6453E-02  2.2149E-02 -1.2815E-01  6.8753E-02 -1.4826E-01
            -2.1140E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1646.04920132754        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     2048             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2343E-01  6.4780E-02  1.5179E+00  1.6306E+00  8.4873E-01  9.7940E-01  3.1306E-01  1.3336E+00  7.1910E-01  1.0474E+00
             8.6040E-01
 PARAMETER:  2.0338E-02 -2.6368E+00  5.1730E-01  5.8893E-01 -6.4013E-02  7.9182E-02 -1.0614E+00  3.8786E-01 -2.2975E-01  1.4632E-01
            -5.0364E-02
 GRADIENT:   4.8251E+02  4.8607E+00  1.2926E+01  1.1463E+03  8.4553E+00  3.9722E+01  1.9604E-02  1.9105E+00  2.3138E+01  1.8525E+00
             8.0281E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1646.05105491736        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2227
 NPARAMETR:  9.2319E-01  6.3923E-02  1.5149E+00  1.6309E+00  8.4908E-01  9.7937E-01  1.2946E-01  1.3325E+00  7.1927E-01  1.0450E+00
             8.6038E-01
 PARAMETER:  2.0085E-02 -2.6501E+00  5.1535E-01  5.8913E-01 -6.3605E-02  7.9155E-02 -1.9444E+00  3.8704E-01 -2.2951E-01  1.4401E-01
            -5.0383E-02
 GRADIENT:   1.4892E+00  7.6535E-02 -1.4101E-01 -1.7273E+01  1.3227E+00  7.8227E-02  5.6165E-04  8.1064E-02  1.1463E-01 -1.4540E-01
            -6.4901E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1646.05381510383        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2417             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2326E-01  6.5306E-02  1.5126E+00  1.6300E+00  8.4752E-01  9.7938E-01  6.2638E-02  1.3298E+00  7.1940E-01  1.0453E+00
             8.6037E-01
 PARAMETER:  2.0150E-02 -2.6287E+00  5.1386E-01  5.8855E-01 -6.5436E-02  7.9167E-02 -2.6704E+00  3.8501E-01 -2.2933E-01  1.4431E-01
            -5.0397E-02
 GRADIENT:   4.8201E+02  4.8789E+00  1.2590E+01  1.1446E+03  8.9613E+00  3.9731E+01  1.8746E-03  1.9031E+00  2.3108E+01  1.6499E+00
             7.9574E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1646.05542851228        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2603
 NPARAMETR:  9.2327E-01  6.6052E-02  1.5104E+00  1.6293E+00  8.4719E-01  9.7939E-01  2.9690E-02  1.3279E+00  7.1966E-01  1.0448E+00
             8.6035E-01
 PARAMETER:  2.0165E-02 -2.6173E+00  5.1237E-01  5.8818E-01 -6.5836E-02  7.9173E-02 -3.4169E+00  3.8359E-01 -2.2897E-01  1.4385E-01
            -5.0418E-02
 GRADIENT:   1.5981E+00  1.1461E-01  5.2867E-01 -1.6922E+01 -1.8578E-01  8.0535E-02  4.2922E-05  8.7275E-02 -4.0096E-02  5.9654E-02
            -8.6616E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1646.05764096294        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2789             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2329E-01  6.7043E-02  1.5069E+00  1.6282E+00  8.4686E-01  9.7940E-01  1.8491E-02  1.3250E+00  7.2015E-01  1.0443E+00
             8.6035E-01
 PARAMETER:  2.0192E-02 -2.6024E+00  5.1008E-01  5.8749E-01 -6.6220E-02  7.9187E-02 -3.8905E+00  3.8143E-01 -2.2830E-01  1.4333E-01
            -5.0414E-02
 GRADIENT:   4.8197E+02  5.0412E+00  1.2020E+01  1.1395E+03  9.9304E+00  3.9725E+01  2.4908E-04  1.9030E+00  2.3054E+01  1.5175E+00
             8.0898E-01

0ITERATION NO.:   84    OBJECTIVE VALUE:  -1646.05874361378        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     2925
 NPARAMETR:  9.2330E-01  6.8425E-02  1.5057E+00  1.6277E+00  8.4609E-01  9.7941E-01  1.3305E-02  1.3224E+00  7.2038E-01  1.0446E+00
             8.6031E-01
 PARAMETER:  2.0200E-02 -2.5820E+00  5.0928E-01  5.8717E-01 -6.7125E-02  7.9192E-02 -4.2196E+00  3.7947E-01 -2.2797E-01  1.4364E-01
            -5.0458E-02
 GRADIENT:  -5.6447E-02  1.6730E-02  7.5691E-01  1.7982E+00 -5.7618E-01 -5.1195E-03  5.2399E-06  1.7145E-02 -1.5281E-01  9.6516E-02
            -2.2411E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2925
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7755E-04 -3.0262E-05 -3.4664E-02 -7.0505E-03 -3.8531E-02
 SE:             2.9863E-02  1.6152E-05  1.8484E-02  2.9375E-02  2.1237E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9526E-01  6.0990E-02  6.0750E-02  8.1032E-01  6.9624E-02

 ETASHRINKSD(%)  1.0000E-10  9.9946E+01  3.8075E+01  1.5897E+00  2.8854E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.1653E+01  3.1542E+00  4.9383E+01
 EBVSHRINKSD(%)  3.2354E-01  9.9948E+01  4.1645E+01  2.0813E+00  2.4967E+01
 EBVSHRINKVR(%)  6.4603E-01  1.0000E+02  6.5947E+01  4.1193E+00  4.3701E+01
 RELATIVEINF(%)  9.7211E+01  1.4702E-06  7.6202E+00  6.2617E+00  8.2444E+00
 EPSSHRINKSD(%)  4.6344E+01
 EPSSHRINKVR(%)  7.1210E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1646.0587436137805     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -910.90791705004233     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    38.20
 Elapsed covariance  time in seconds:     5.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1646.059       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.23E-01  6.84E-02  1.51E+00  1.63E+00  8.46E-01  9.79E-01  1.33E-02  1.32E+00  7.20E-01  1.04E+00  8.60E-01
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         2.70E-02  3.16E-01  3.91E-01  2.05E-01  1.18E-01  8.17E-02  1.26E-01  4.63E-01  1.17E-01  1.40E-01  5.92E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        7.30E-04
 
 TH 2
+        3.60E-04  9.97E-02
 
 TH 3
+        1.72E-03 -5.30E-02  1.53E-01
 
 TH 4
+       -1.27E-05 -6.34E-02  4.13E-02  4.21E-02
 
 TH 5
+        6.76E-04  1.07E-02  3.19E-02 -4.39E-03  1.39E-02
 
 TH 6
+       -5.89E-04 -3.19E-03 -1.21E-03  1.98E-03 -1.50E-03  6.68E-03
 
 TH 7
+        2.97E-04 -2.52E-02  1.03E-02  1.58E-02 -4.77E-03  6.83E-04  1.59E-02
 
 TH 8
+        1.11E-03 -7.57E-02  1.58E-01  5.54E-02  2.65E-02 -1.46E-03  3.15E-02  2.14E-01
 
 TH 9
+        9.66E-05  3.40E-02 -2.03E-02 -2.15E-02  2.96E-03 -8.08E-04 -7.90E-03 -2.79E-02  1.36E-02
 
 TH10
+        1.42E-04  1.17E-02  1.15E-02 -6.22E-03  8.48E-03 -2.12E-03 -1.44E-02 -3.63E-03  2.81E-03  1.97E-02
 
 TH11
+        1.81E-04  1.05E-03 -2.02E-03 -8.17E-04 -5.92E-04 -1.89E-04 -8.11E-04 -3.77E-03  1.83E-04 -3.63E-05  3.50E-03
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.70E-02
 
 TH 2
+        4.22E-02  3.16E-01
 
 TH 3
+        1.63E-01 -4.30E-01  3.91E-01
 
 TH 4
+       -2.29E-03 -9.79E-01  5.15E-01  2.05E-01
 
 TH 5
+        2.12E-01  2.87E-01  6.92E-01 -1.81E-01  1.18E-01
 
 TH 6
+       -2.67E-01 -1.23E-01 -3.80E-02  1.18E-01 -1.56E-01  8.17E-02
 
 TH 7
+        8.72E-02 -6.34E-01  2.10E-01  6.14E-01 -3.21E-01  6.64E-02  1.26E-01
 
 TH 8
+        8.91E-02 -5.18E-01  8.73E-01  5.84E-01  4.85E-01 -3.87E-02  5.40E-01  4.63E-01
 
 TH 9
+        3.06E-02  9.23E-01 -4.45E-01 -8.99E-01  2.14E-01 -8.47E-02 -5.38E-01 -5.17E-01  1.17E-01
 
 TH10
+        3.74E-02  2.65E-01  2.09E-01 -2.16E-01  5.12E-01 -1.84E-01 -8.15E-01 -5.60E-02  1.72E-01  1.40E-01
 
 TH11
+        1.13E-01  5.62E-02 -8.75E-02 -6.74E-02 -8.47E-02 -3.91E-02 -1.09E-01 -1.38E-01  2.65E-02 -4.37E-03  5.92E-02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.91E+05
 
 TH 2
+       -9.18E+04  2.21E+04
 
 TH 3
+       -9.45E+04  2.23E+04  2.30E+04
 
 TH 4
+        5.79E+03 -8.59E+02 -1.46E+03  9.01E+02
 
 TH 5
+       -5.10E+03  9.39E+02  9.75E+02 -1.26E+02  9.34E+02
 
 TH 6
+       -4.39E+04  1.04E+04  1.07E+04 -7.04E+02  5.73E+02  5.16E+03
 
 TH 7
+       -6.70E+05  1.58E+05  1.63E+05 -1.02E+04  8.58E+03  7.58E+04  1.15E+06
 
 TH 8
+        1.32E+05 -3.10E+04 -3.20E+04  2.00E+03 -1.70E+03 -1.49E+04 -2.27E+05  4.46E+04
 
 TH 9
+        5.18E+04 -1.24E+04 -1.25E+04  6.82E+02 -6.86E+02 -5.86E+03 -8.91E+04  1.75E+04  7.42E+03
 
 TH10
+       -3.67E+05  8.65E+04  8.91E+04 -5.62E+03  4.59E+03  4.16E+04  6.32E+05 -1.24E+05 -4.88E+04  3.46E+05
 
 TH11
+       -6.90E+04  1.62E+04  1.67E+04 -1.04E+03  9.50E+02  7.80E+03  1.19E+05 -2.33E+04 -9.13E+03  6.50E+04  1.25E+04
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       43.899
Stop Time:
Wed Sep 29 17:54:56 CDT 2021
