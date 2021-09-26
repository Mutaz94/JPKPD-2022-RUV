Sat Sep 25 08:06:28 CDT 2021
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
$DATA ../../../../data/spa/A1/dat53.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m53.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1488.53814537336        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.9229E+01 -6.8302E+01 -2.3663E+00 -7.4641E+01  3.9254E+01 -2.2549E+01  1.0076E+01 -3.4257E+00  3.7148E+01  7.8366E+00
            -3.3020E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1563.24039768740        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.4664E-01  1.1005E+00  1.0266E+00  1.0231E+00  9.9258E-01  1.0598E+00  8.5828E-01  9.8045E-01  6.9285E-01  8.1833E-01
             1.7794E+00
 PARAMETER:  4.5160E-02  1.9576E-01  1.2627E-01  1.2284E-01  9.2553E-02  1.5808E-01 -5.2820E-02  8.0260E-02 -2.6694E-01 -1.0049E-01
             6.7629E-01
 GRADIENT:  -6.2246E+01  3.6230E+01  1.5140E+01  2.5988E+01 -2.7674E+01  3.1072E+00 -3.8312E+00 -1.6556E+00 -7.7466E+00  5.6462E+00
             3.3655E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1567.50900137352        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.5254E-01  1.1196E+00  5.7321E-01  9.7669E-01  7.5603E-01  1.0666E+00  9.4425E-01  7.4961E-01  7.5642E-01  4.0863E-01
             1.7097E+00
 PARAMETER:  5.1377E-02  2.1293E-01 -4.5651E-01  7.6419E-02 -1.7967E-01  1.6449E-01  4.2631E-02 -1.8820E-01 -1.7916E-01 -7.9494E-01
             6.3631E-01
 GRADIENT:  -5.5627E+01  2.8811E+01 -2.5808E+00  4.8465E+01  4.5160E+00  3.1238E+00  1.8099E+00  8.0853E-01  8.7580E+00 -5.8352E-02
             2.6308E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1574.25151552975        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.8978E-01  1.0176E+00  2.6771E-01  9.1692E-01  5.0270E-01  1.0848E+00  9.5616E-01  2.6622E-01  6.9065E-01  2.7848E-01
             1.5081E+00
 PARAMETER:  8.9731E-02  1.1743E-01 -1.2178E+00  1.3267E-02 -5.8776E-01  1.8138E-01  5.5167E-02 -1.2234E+00 -2.7013E-01 -1.1784E+00
             5.1087E-01
 GRADIENT:   3.5811E+01  2.1988E+01 -9.7845E+00  4.2208E+01  5.0505E+00  1.0824E+01  7.5344E+00  3.4394E-01  3.4396E+00  3.5738E+00
            -8.8991E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1575.27125707498        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.7108E-01  9.5915E-01  2.5332E-01  9.1662E-01  4.7337E-01  1.0569E+00  9.6501E-01  2.4793E-01  6.7106E-01  1.7170E-01
             1.5386E+00
 PARAMETER:  7.0650E-02  5.8294E-02 -1.2731E+00  1.2936E-02 -6.4788E-01  1.5538E-01  6.4382E-02 -1.2946E+00 -2.9890E-01 -1.6620E+00
             5.3086E-01
 GRADIENT:   1.1959E+00 -1.9398E+00  1.3071E+00 -4.7599E+00 -3.5389E+00  2.6292E-01  9.4771E-01  7.8234E-02  6.1956E-01  1.0708E+00
             9.3916E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1576.19159706986        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.6893E-01  1.0787E+00  2.3983E-01  8.5252E-01  5.1340E-01  1.0542E+00  8.9343E-01  3.9110E-01  6.9558E-01  4.0577E-02
             1.5464E+00
 PARAMETER:  6.8433E-02  1.7576E-01 -1.3278E+00 -5.9558E-02 -5.6671E-01  1.5280E-01 -1.2692E-02 -8.3879E-01 -2.6300E-01 -3.1046E+00
             5.3591E-01
 GRADIENT:   1.1422E+00 -2.4063E+00 -1.0018E+00  1.7215E-01  3.5089E+00  5.2037E-01  3.5381E-01  5.1108E-02 -7.9760E-02  5.5383E-02
            -2.7626E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1576.73899873674        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      552
 NPARAMETR:  9.6826E-01  1.3352E+00  2.1567E-01  7.2425E-01  6.0741E-01  1.0567E+00  7.7606E-01  6.9410E-01  7.5087E-01  1.0000E-02
             1.5578E+00
 PARAMETER:  6.7743E-02  3.8906E-01 -1.4340E+00 -2.2261E-01 -3.9855E-01  1.5519E-01 -1.5353E-01 -2.6514E-01 -1.8652E-01 -5.3471E+00
             5.4331E-01
 GRADIENT:  -1.0767E+01  2.5652E+01  5.9426E+00  5.6575E+00 -1.9655E+01  2.8437E-01 -2.2700E+00 -6.0397E-01 -2.4063E+00  0.0000E+00
            -3.3455E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1577.71199012819        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      735
 NPARAMETR:  9.7402E-01  1.5843E+00  1.7896E-01  6.0037E-01  7.0085E-01  1.0475E+00  6.9742E-01  1.2567E+00  8.4702E-01  1.0000E-02
             1.5742E+00
 PARAMETER:  7.3678E-02  5.6016E-01 -1.6206E+00 -4.1021E-01 -2.5546E-01  1.4640E-01 -2.6037E-01  3.2851E-01 -6.6030E-02 -8.6565E+00
             5.5376E-01
 GRADIENT:   3.7459E+00  7.8163E+01  5.9839E+00  3.7252E+01 -3.5520E+01 -1.4857E+00 -4.5165E+00 -4.1684E+00  1.9252E+00  0.0000E+00
             1.3218E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1578.55275831747        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      917
 NPARAMETR:  9.7408E-01  1.5868E+00  1.7865E-01  5.9220E-01  7.0186E-01  1.0465E+00  6.9667E-01  1.4969E+00  8.3747E-01  1.0000E-02
             1.5741E+00
 PARAMETER:  7.3738E-02  5.6174E-01 -1.6223E+00 -4.2391E-01 -2.5402E-01  1.4543E-01 -2.6144E-01  5.0338E-01 -7.7369E-02 -8.6907E+00
             5.5367E-01
 GRADIENT:   2.0075E+01  8.7240E+01  9.1206E+00  2.7204E+01 -5.5015E+01  1.9402E+00 -3.2973E+00  3.5898E-01  2.4275E+00  0.0000E+00
             2.3061E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1578.73303550527        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:     1031
 NPARAMETR:  9.7408E-01  1.5868E+00  1.7865E-01  5.8382E-01  7.0186E-01  1.0478E+00  6.9667E-01  1.5024E+00  8.2238E-01  1.0000E-02
             1.5741E+00
 PARAMETER:  7.3738E-02  5.6174E-01 -1.6223E+00 -4.3817E-01 -2.5402E-01  1.4670E-01 -2.6144E-01  5.0709E-01 -9.5558E-02 -8.6907E+00
             5.5367E-01
 GRADIENT:   6.0137E+00  6.1782E+01  1.4202E+01 -7.3184E-01 -7.1101E+01 -6.3132E-01 -3.5844E+00 -1.5365E-01 -8.3019E-03  0.0000E+00
             2.2194E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1578.77695934982        NO. OF FUNC. EVALS.: 204
 CUMULATIVE NO. OF FUNC. EVALS.:     1235             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7405E-01  1.5868E+00  1.7856E-01  5.8399E-01  7.0202E-01  1.0490E+00  6.9800E-01  1.5053E+00  8.2207E-01  1.0000E-02
             1.5732E+00
 PARAMETER:  7.3711E-02  5.6172E-01 -1.6228E+00 -4.3787E-01 -2.5379E-01  1.4781E-01 -2.5954E-01  5.0899E-01 -9.5931E-02 -8.6907E+00
             5.5313E-01
 GRADIENT:   2.0961E+01  7.8394E+01  1.4941E+01  4.8266E+00 -6.8109E+01  3.1415E+00 -2.7340E+00  2.2174E-02  7.1363E-02  0.0000E+00
             2.2544E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1578.78216676580        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:     1389
 NPARAMETR:  9.7357E-01  1.5868E+00  1.7855E-01  5.8400E-01  7.0203E-01  1.0488E+00  6.9813E-01  1.5053E+00  8.2207E-01  1.0000E-02
             1.5733E+00
 PARAMETER:  7.3213E-02  5.6170E-01 -1.6229E+00 -4.3786E-01 -2.5378E-01  1.4764E-01 -2.5935E-01  5.0899E-01 -9.5928E-02 -8.6907E+00
             5.5315E-01
 GRADIENT:   4.9798E+00  6.1057E+01  1.3628E+01  2.5551E-01 -6.9536E+01 -2.7412E-01 -3.0517E+00 -6.2800E-02  9.8441E-03  0.0000E+00
             2.2142E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1578.91658779122        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1552
 NPARAMETR:  9.7396E-01  1.5832E+00  1.7739E-01  5.8502E-01  7.0185E-01  1.0498E+00  6.9831E-01  1.5022E+00  8.2240E-01  1.0000E-02
             1.5767E+00
 PARAMETER:  7.3613E-02  5.5945E-01 -1.6294E+00 -4.3611E-01 -2.5403E-01  1.4860E-01 -2.5909E-01  5.0695E-01 -9.5528E-02 -8.6907E+00
             5.5536E-01
 GRADIENT:   5.8760E+00  5.2330E+01  1.0455E+01  4.2699E+00 -5.9661E+01  9.2094E-02 -2.6726E+00  2.4469E-01  3.0550E-01  0.0000E+00
             2.3081E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1579.49924699569        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1713             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7292E-01  1.5788E+00  1.7635E-01  5.8449E-01  7.0417E-01  1.0492E+00  6.9964E-01  1.5021E+00  8.2235E-01  1.0000E-02
             1.5593E+00
 PARAMETER:  7.2549E-02  5.5664E-01 -1.6353E+00 -4.3701E-01 -2.5074E-01  1.4802E-01 -2.5719E-01  5.0687E-01 -9.5589E-02 -8.6907E+00
             5.4427E-01
 GRADIENT:   1.9469E+01  5.1147E+01  7.4529E+00  1.2209E+01 -3.9063E+01  3.2549E+00 -1.7716E+00 -2.1455E-01 -2.0986E-01  0.0000E+00
             1.9187E+01

0ITERATION NO.:   69    OBJECTIVE VALUE:  -1579.51780724311        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:     1867
 NPARAMETR:  9.7302E-01  1.5787E+00  1.7606E-01  5.8475E-01  7.0399E-01  1.0491E+00  6.9982E-01  1.5013E+00  8.2243E-01  1.0000E-02
             1.5602E+00
 PARAMETER:  7.2550E-02  5.5607E-01 -1.6353E+00 -4.3700E-01 -2.5074E-01  1.4804E-01 -2.5719E-01  5.0686E-01 -9.5588E-02 -8.6907E+00
             5.4427E-01
 GRADIENT:  -3.7895E+05 -6.8125E+04  1.1601E+04 -8.6721E+04  1.5109E+05  2.5598E+05 -7.3674E+04  7.4704E+04 -1.8948E+05  0.0000E+00
            -6.9621E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1867
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.7423E-03 -1.7622E-02 -9.7345E-03  2.7452E-03 -2.4072E-04
 SE:             2.9716E-02  2.6786E-02  1.5960E-02  2.1045E-02  3.4011E-04
 N:                     100         100         100         100         100

 P VAL.:         9.5325E-01  5.1060E-01  5.4190E-01  8.9622E-01  4.7909E-01

 ETASHRINKSD(%)  4.4679E-01  1.0265E+01  4.6533E+01  2.9495E+01  9.8861E+01
 ETASHRINKVR(%)  8.9159E-01  1.9476E+01  7.1413E+01  5.0291E+01  9.9987E+01
 EBVSHRINKSD(%)  8.1132E-01  1.0992E+01  4.6790E+01  2.9731E+01  9.8861E+01
 EBVSHRINKVR(%)  1.6161E+00  2.0775E+01  7.1687E+01  5.0623E+01  9.9987E+01
 RELATIVEINF(%)  9.7225E+01  5.7922E+00  6.7450E+00  3.7128E+00  9.9284E-04
 EPSSHRINKSD(%)  4.0979E+01
 EPSSHRINKVR(%)  6.5165E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1579.5178072431140     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -844.36698067937584     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.57
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1579.518       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.73E-01  1.58E+00  1.76E-01  5.84E-01  7.04E-01  1.05E+00  7.00E-01  1.50E+00  8.22E-01  1.00E-02  1.56E+00
 


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
+        1.00E+09
 
 TH 2
+        2.52E+04  1.97E+07
 
 TH 3
+        2.20E+03 -7.98E+03  1.14E+08
 
 TH 4
+       -2.61E+03  9.96E+03 -1.97E+05  1.45E+08
 
 TH 5
+        3.66E+03 -1.49E+04  2.80E+05  2.10E+08  3.04E+08
 
 TH 6
+       -9.13E+03 -1.02E+03  3.08E+03 -3.49E+03  5.02E+03  7.85E+08
 
 TH 7
+       -1.95E+04 -2.17E+03  6.49E+03 -7.45E+03  1.08E+04 -4.93E+03  2.93E+08
 
 TH 8
+        8.61E+02 -3.23E+03  6.53E+04 -1.29E+04 -6.07E+03  1.17E+03  2.49E+03  1.63E+07
 
 TH 9
+       -1.18E+09 -1.31E+08 -7.11E+03  8.08E+03 -1.16E+04 -1.08E+04 -6.40E+08 -2.68E+03  1.40E+09
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        4.78E+04  2.88E+03 -1.61E+04  1.82E+04 -2.64E+04 -1.04E+03 -2.21E+03 -6.10E+03 -1.36E+08  0.00E+00  1.32E+07
 
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
 #CPUT: Total CPU Time in Seconds,       28.623
Stop Time:
Sat Sep 25 08:06:58 CDT 2021
