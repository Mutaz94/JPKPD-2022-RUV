Thu Sep 30 00:32:59 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat77.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -758.456214534339        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9300E+02  4.2893E+01  2.2525E+02 -1.4475E+01  8.2618E+01  4.6531E+01 -1.9941E+00 -1.9293E+02  1.2927E+01 -5.1018E+01
            -2.4158E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1653.11033789780        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0463E+00  1.0221E+00  8.9634E-01  1.0515E+00  1.0042E+00  9.0089E-01  9.1218E-01  1.0175E+00  8.3281E-01  8.9991E-01
             2.3306E+00
 PARAMETER:  1.4522E-01  1.2189E-01 -9.4395E-03  1.5020E-01  1.0417E-01 -4.3728E-03  8.0855E-03  1.1736E-01 -8.2954E-02 -5.4581E-03
             9.4612E-01
 GRADIENT:   1.8398E+02  8.1824E+00  1.6686E+00  2.2710E+01  3.1240E+01 -2.1291E+01  1.5228E+00 -1.4343E+00 -4.9588E+00  4.6478E+00
            -1.6846E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1668.69295681332        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      180
 NPARAMETR:  1.0425E+00  8.8123E-01  3.7196E-01  1.0807E+00  5.0032E-01  9.2246E-01  7.5036E-01  7.6493E-01  1.0017E+00  2.1128E-01
             2.2878E+00
 PARAMETER:  1.4165E-01 -2.6442E-02 -8.8898E-01  1.7765E-01 -5.9251E-01  1.9287E-02 -1.8720E-01 -1.6797E-01  1.0167E-01 -1.4546E+00
             9.2760E-01
 GRADIENT:   4.0849E+01  1.6261E+02  1.3626E+02  4.1532E+01 -2.2555E+02 -2.4390E+01 -3.2926E+01 -1.5958E+01  1.8344E+01 -4.4959E+00
            -1.3350E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1690.86715825572        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      356
 NPARAMETR:  1.0208E+00  8.2195E-01  3.0632E-01  1.0899E+00  4.7269E-01  9.6036E-01  9.2844E-01  4.9075E-01  9.2807E-01  1.4335E-01
             2.6672E+00
 PARAMETER:  1.2063E-01 -9.6070E-02 -1.0831E+00  1.8608E-01 -6.4933E-01  5.9558E-02  2.5749E-02 -6.1182E-01  2.5354E-02 -1.8424E+00
             1.0810E+00
 GRADIENT:  -1.8655E+01  2.0902E+01 -1.9158E+01  6.3582E+01  1.2555E+01 -5.6863E+00 -1.0060E+01 -3.3839E-02  8.0822E+00 -1.0334E-01
             3.3740E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1699.50514272995        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  1.0250E+00  5.3038E-01  2.6780E-01  1.1338E+00  3.5083E-01  9.7876E-01  1.3196E+00  4.4300E-01  8.8185E-01  1.1944E-01
             2.4488E+00
 PARAMETER:  1.2467E-01 -5.3416E-01 -1.2175E+00  2.2554E-01 -9.4744E-01  7.8526E-02  3.7736E-01 -7.1418E-01 -2.5733E-02 -2.0249E+00
             9.9558E-01
 GRADIENT:   9.5684E-01  3.2323E-01  1.2318E+01 -2.3587E+00 -8.2760E+00 -9.5390E-01 -8.9439E+00 -4.0789E+00  5.3255E+00 -8.0023E-01
            -1.3077E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1701.30524119157        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  1.0245E+00  4.5362E-01  2.4538E-01  1.1376E+00  3.1456E-01  9.8438E-01  1.4604E+00  6.8255E-01  8.7719E-01  8.2858E-02
             2.4119E+00
 PARAMETER:  1.2421E-01 -6.9050E-01 -1.3049E+00  2.2893E-01 -1.0566E+00  8.4258E-02  4.7870E-01 -2.8192E-01 -3.1029E-02 -2.3906E+00
             9.8041E-01
 GRADIENT:   2.7339E-01 -6.4736E-01  5.4770E-01 -1.8691E+00  2.0381E+00  2.4869E-01 -5.4335E-01 -1.1926E-01  4.6015E-01 -1.9421E-01
            -1.5961E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1702.22773683037        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      827
 NPARAMETR:  1.0251E+00  4.4496E-01  2.4268E-01  1.1368E+00  3.0841E-01  9.8163E-01  1.4104E+00  6.0691E-01  8.7449E-01  3.5208E-01
             2.4097E+00
 PARAMETER:  1.2481E-01 -7.0977E-01 -1.3160E+00  2.2825E-01 -1.0763E+00  8.1460E-02  4.4387E-01 -3.9937E-01 -3.4110E-02 -9.4391E-01
             9.7951E-01
 GRADIENT:   7.1225E+01  1.6400E+01  3.1493E+01  2.6127E+01  7.3975E+01  4.7374E+00  5.4619E+00  3.0286E+00  1.2329E+00  1.5795E+00
             1.3788E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1703.18221203189        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      897
 NPARAMETR:  1.0235E+00  4.1673E-01  2.3462E-01  1.1396E+00  2.9631E-01  9.8000E-01  1.1164E+00  2.9035E-01  8.8105E-01  5.7965E-01
             2.3941E+00
 PARAMETER:  1.2326E-01 -7.7531E-01 -1.3498E+00  2.3071E-01 -1.1163E+00  7.9797E-02  2.1013E-01 -1.1367E+00 -2.6641E-02 -4.4533E-01
             9.7301E-01
 GRADIENT:   6.7649E+01  1.0515E+01  3.1535E+01  3.7633E+01  9.1430E+01  4.4674E+00  2.3403E+00  9.8045E-01 -1.6537E+00  7.1313E+00
             5.2105E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1703.92575292701        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      968
 NPARAMETR:  1.0218E+00  3.3877E-01  2.0163E-01  1.1279E+00  2.5950E-01  9.8107E-01  5.7967E-01  3.3213E-02  9.2018E-01  7.0133E-01
             2.3899E+00
 PARAMETER:  1.2155E-01 -9.8243E-01 -1.5013E+00  2.2034E-01 -1.2490E+00  8.0893E-02 -4.4529E-01 -3.3048E+00  1.6814E-02 -2.5478E-01
             9.7124E-01
 GRADIENT:   6.0679E+01 -8.2805E+00  1.2739E+01  6.8408E+01  1.4900E+02  4.4881E+00  1.6194E+00  2.3343E-02 -1.6643E+00  1.9139E+01
             7.2560E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1705.99434412715        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1147
 NPARAMETR:  1.0217E+00  3.2175E-01  1.9074E-01  1.0972E+00  2.4475E-01  9.8153E-01  3.9528E-01  1.2592E-02  9.3907E-01  6.1405E-01
             2.4250E+00
 PARAMETER:  1.2143E-01 -1.0340E+00 -1.5568E+00  1.9276E-01 -1.3075E+00  8.1357E-02 -8.2817E-01 -4.2747E+00  3.7133E-02 -3.8767E-01
             9.8581E-01
 GRADIENT:  -5.9770E+00 -8.4445E+00  1.2275E+01  5.5139E+00 -5.1340E+00 -3.6570E-01  3.8737E-01  1.4073E-04 -1.1578E+00 -4.8156E+00
             4.1087E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1707.26198129592        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1322
 NPARAMETR:  1.0243E+00  3.3364E-01  1.6962E-01  1.0562E+00  2.3199E-01  9.8654E-01  2.5661E-01  1.0000E-02  9.8510E-01  6.3696E-01
             2.3872E+00
 PARAMETER:  1.2396E-01 -9.9769E-01 -1.6742E+00  1.5468E-01 -1.3610E+00  8.6448E-02 -1.2602E+00 -5.3707E+00  8.4992E-02 -3.5104E-01
             9.7012E-01
 GRADIENT:  -3.1376E-01  1.3604E-01 -5.6586E-01 -7.0446E-01 -1.0898E+00 -1.9463E-02  3.4159E-01  0.0000E+00  1.6683E-01 -2.6943E-01
             3.8649E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1707.43086696152        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1461
 NPARAMETR:  1.0237E+00  3.3244E-01  1.6944E-01  1.0563E+00  2.3170E-01  9.8632E-01  2.6333E-02  1.0000E-02  9.8438E-01  6.4142E-01
             2.3857E+00
 PARAMETER:  1.2345E-01 -1.0013E+00 -1.6752E+00  1.5477E-01 -1.3623E+00  8.6228E-02 -3.5369E+00 -5.3707E+00  8.4255E-02 -3.4408E-01
             9.6947E-01
 GRADIENT:  -1.3965E+00 -5.5486E-01 -3.2235E-01  4.0415E-01 -1.0312E+00 -6.7548E-02  3.7381E-03  0.0000E+00 -1.6474E-01  9.5424E-02
            -4.7695E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1707.43408571740        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1635
 NPARAMETR:  1.0240E+00  3.3340E-01  1.6935E-01  1.0558E+00  2.3180E-01  9.8640E-01  1.0000E-02  1.0000E-02  9.8512E-01  6.4106E-01
             2.3867E+00
 PARAMETER:  1.2375E-01 -9.9843E-01 -1.6758E+00  1.5427E-01 -1.3619E+00  8.6306E-02 -4.6610E+00 -5.3707E+00  8.5004E-02 -3.4463E-01
             9.6993E-01
 GRADIENT:  -7.8418E-01 -9.1739E-02 -3.9640E-01 -1.6897E-01 -1.1416E+00 -5.9009E-02  0.0000E+00  0.0000E+00 -2.4069E-02 -2.5830E-02
            -7.0149E-02

0ITERATION NO.:   61    OBJECTIVE VALUE:  -1707.43408571740        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1657
 NPARAMETR:  1.0240E+00  3.3340E-01  1.6935E-01  1.0558E+00  2.3180E-01  9.8640E-01  1.0000E-02  1.0000E-02  9.8512E-01  6.4106E-01
             2.3867E+00
 PARAMETER:  1.2375E-01 -9.9843E-01 -1.6758E+00  1.5427E-01 -1.3619E+00  8.6306E-02 -4.6610E+00 -5.3707E+00  8.5004E-02 -3.4463E-01
             9.6993E-01
 GRADIENT:  -7.8418E-01 -9.1739E-02 -3.9640E-01 -1.6897E-01 -1.1416E+00 -5.9009E-02  0.0000E+00  0.0000E+00 -2.4069E-02 -2.5830E-02
            -7.0149E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1657
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.0398E-04 -7.5268E-05  1.8007E-04 -9.4669E-03  1.0179E-03
 SE:             2.9361E-02  1.2439E-04  2.3240E-04  2.6707E-02  2.5090E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7544E-01  5.4512E-01  4.3843E-01  7.2298E-01  9.6764E-01

 ETASHRINKSD(%)  1.6368E+00  9.9583E+01  9.9221E+01  1.0530E+01  1.5945E+01
 ETASHRINKVR(%)  3.2468E+00  9.9998E+01  9.9994E+01  1.9951E+01  2.9347E+01
 EBVSHRINKSD(%)  1.7077E+00  9.9559E+01  9.9229E+01  9.3244E+00  1.6262E+01
 EBVSHRINKVR(%)  3.3862E+00  9.9998E+01  9.9994E+01  1.7779E+01  2.9880E+01
 RELATIVEINF(%)  9.6493E+01  4.1659E-04  4.9481E-04  3.5414E+01  3.4772E+00
 EPSSHRINKSD(%)  2.7658E+01
 EPSSHRINKVR(%)  4.7666E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1707.4340857174041     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -788.49555251273136     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.15
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1707.434       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  3.33E-01  1.69E-01  1.06E+00  2.32E-01  9.86E-01  1.00E-02  1.00E-02  9.85E-01  6.41E-01  2.39E+00
 


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
+        1.05E+03
 
 TH 2
+       -4.33E+01  1.87E+03
 
 TH 3
+       -5.63E+01  3.58E+03  2.25E+04
 
 TH 4
+       -1.81E+01  1.44E+02 -1.23E+03  8.47E+02
 
 TH 5
+        1.56E+02 -6.59E+03 -2.43E+04 -4.64E+02  3.53E+04
 
 TH 6
+        2.39E+00 -1.41E+01  3.92E+01 -8.26E+00  6.83E+00  1.91E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.22E+00 -3.04E+01  2.70E+02 -2.13E+01  5.76E+01 -1.56E+00  0.00E+00  0.00E+00  1.34E+02
 
 TH10
+       -2.37E+00 -3.74E+01 -8.09E+01  2.01E+01  4.86E+01  1.71E+00  0.00E+00  0.00E+00  6.59E+00  2.35E+02
 
 TH11
+       -1.43E+01 -8.35E+00 -4.77E+01 -9.89E+00  4.43E+01  2.57E+00  0.00E+00  0.00E+00  1.00E+01  2.36E+01  7.89E+01
 
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
 #CPUT: Total CPU Time in Seconds,       29.898
Stop Time:
Thu Sep 30 00:33:30 CDT 2021
