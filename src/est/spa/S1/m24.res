Wed Sep 29 14:07:30 CDT 2021
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
$DATA ../../../../data/spa/S1/dat24.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m24.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1648.05252518744        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4498E+02 -3.9770E+01  3.5023E+01 -9.9044E+01 -3.8463E+01  3.2089E+01 -2.9668E+00  4.7672E+00 -1.0910E+01 -1.3026E+00
             1.9691E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1658.17023535414        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.7458E-01  1.0307E+00  9.8259E-01  1.0833E+00  1.0307E+00  1.0240E+00  1.0172E+00  9.7213E-01  1.0101E+00  1.0259E+00
             9.3851E-01
 PARAMETER:  7.4256E-02  1.3023E-01  8.2433E-02  1.7998E-01  1.3023E-01  1.2371E-01  1.1707E-01  7.1730E-02  1.1005E-01  1.2558E-01
             3.6540E-02
 GRADIENT:   7.9096E-01  8.4133E-01  5.7206E-01 -7.7029E-01 -1.0171E+00  2.3142E+00 -2.9384E+00  5.2763E+00  2.6015E-01 -4.3678E+00
            -5.5551E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1660.09198886253        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  9.7390E-01  8.6874E-01  9.0079E-01  1.2020E+00  9.1997E-01  1.0675E+00  1.4058E+00  5.0574E-01  8.7378E-01  1.0034E+00
             9.7315E-01
 PARAMETER:  7.3550E-02 -4.0711E-02 -4.4847E-03  2.8399E-01  1.6587E-02  1.6534E-01  4.4058E-01 -5.8172E-01 -3.4924E-02  1.0344E-01
             7.2785E-02
 GRADIENT:  -3.3923E-01  2.0366E+01 -7.3555E+00  4.2256E+01  1.4587E+01  1.8062E+01  1.8958E+00  4.2208E-01 -1.9168E+00  2.9743E-01
             7.9481E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1662.08856508803        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      538
 NPARAMETR:  9.7545E-01  7.3478E-01  7.9009E-01  1.2673E+00  7.7348E-01  1.0168E+00  1.5505E+00  3.0761E-01  8.3101E-01  8.7225E-01
             9.4990E-01
 PARAMETER:  7.5142E-02 -2.0819E-01 -1.3561E-01  3.3688E-01 -1.5685E-01  1.1662E-01  5.3857E-01 -1.0789E+00 -8.5109E-02 -3.6685E-02
             4.8601E-02
 GRADIENT:   2.6977E+00  2.3268E+01  7.6350E+00  4.3201E+01 -1.1098E+01 -5.3330E-01 -1.9397E+00  3.2508E-02 -2.5491E+00 -3.5211E+00
            -1.0487E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1664.36532502139        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      713
 NPARAMETR:  9.6998E-01  4.2492E-01  7.5896E-01  1.4138E+00  6.6039E-01  1.0197E+00  2.3169E+00  6.9567E-02  7.6265E-01  8.9297E-01
             9.4764E-01
 PARAMETER:  6.9523E-02 -7.5585E-01 -1.7580E-01  4.4627E-01 -3.1493E-01  1.1952E-01  9.4021E-01 -2.5655E+00 -1.7096E-01 -1.3204E-02
             4.6217E-02
 GRADIENT:  -5.0593E-01  7.1793E+00  5.5087E+00  1.4422E+01 -9.8346E+00  2.3782E+00  1.5641E+00 -2.1967E-02  3.1550E-02  8.5960E-01
             7.5594E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1664.61258912726        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      889
 NPARAMETR:  9.6789E-01  3.1530E-01  7.5692E-01  1.4636E+00  6.3567E-01  1.0106E+00  2.7318E+00  2.8394E-02  7.4414E-01  9.0175E-01
             9.4534E-01
 PARAMETER:  6.7363E-02 -1.0542E+00 -1.7849E-01  4.8093E-01 -3.5308E-01  1.1055E-01  1.1050E+00 -3.4616E+00 -1.9553E-01 -3.4229E-03
             4.3791E-02
 GRADIENT:   1.5152E-02  1.2709E+00  1.4563E+00  3.5928E+00 -2.5003E+00 -1.3341E-01  1.2501E-01 -6.5426E-03 -4.2625E-01 -2.0404E-01
            -7.8067E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1664.64309309181        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1072
 NPARAMETR:  9.6842E-01  3.1486E-01  7.5659E-01  1.4586E+00  6.3589E-01  1.0112E+00  2.7571E+00  3.2978E-02  7.4447E-01  9.0155E-01
             9.4492E-01
 PARAMETER:  6.7913E-02 -1.0556E+00 -1.7894E-01  4.7750E-01 -3.5273E-01  1.1112E-01  1.1142E+00 -3.3119E+00 -1.9508E-01 -3.6391E-03
             4.3349E-02
 GRADIENT:   1.2996E+00  5.8809E-01  1.7613E+00 -8.0280E+00 -1.8108E+00  1.1455E-01  1.3512E+00 -8.2830E-03  4.9184E-02  3.8538E-02
            -6.7436E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1664.64786379505        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1260             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6848E-01  3.1562E-01  7.5535E-01  1.4574E+00  6.3626E-01  1.0112E+00  2.7432E+00  4.1391E-02  7.4477E-01  9.0066E-01
             9.4501E-01
 PARAMETER:  6.7971E-02 -1.0532E+00 -1.8058E-01  4.7663E-01 -3.5215E-01  1.1118E-01  1.1091E+00 -3.0847E+00 -1.9467E-01 -4.6312E-03
             4.3439E-02
 GRADIENT:   4.0160E+02  4.6962E+01  1.0272E+01  6.5651E+02  4.7577E+01  4.7912E+01  4.4105E+01  1.7108E-02  1.7042E+01  1.2398E+00
             9.4725E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1664.64884369398        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1443
 NPARAMETR:  9.6850E-01  3.1632E-01  7.5467E-01  1.4569E+00  6.3636E-01  1.0113E+00  2.7378E+00  4.6769E-02  7.4496E-01  9.0001E-01
             9.4500E-01
 PARAMETER:  6.7996E-02 -1.0510E+00 -1.8147E-01  4.7629E-01 -3.5199E-01  1.1121E-01  1.1072E+00 -2.9625E+00 -1.9442E-01 -5.3503E-03
             4.3425E-02
 GRADIENT:   1.3144E+00 -2.3068E-01 -2.4118E+00 -8.5925E+00  2.7946E+00  1.2150E-01  8.8349E-01 -1.3662E-02  1.1343E-01  2.1237E-01
             2.2936E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1664.65458680679        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1631             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6853E-01  3.1826E-01  7.5531E-01  1.4563E+00  6.3586E-01  1.0113E+00  2.7351E+00  5.4108E-02  7.4504E-01  8.9891E-01
             9.4463E-01
 PARAMETER:  6.8026E-02 -1.0449E+00 -1.8063E-01  4.7589E-01 -3.5278E-01  1.1124E-01  1.1062E+00 -2.8168E+00 -1.9432E-01 -6.5748E-03
             4.3038E-02
 GRADIENT:   4.0199E+02  4.7618E+01  1.2610E+01  6.5512E+02  4.4784E+01  4.7988E+01  4.4141E+01  2.4479E-02  1.6924E+01  1.0761E+00
             7.2677E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1664.65721993118        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1811
 NPARAMETR:  9.6855E-01  3.1913E-01  7.5527E-01  1.4559E+00  6.3582E-01  1.0113E+00  2.7316E+00  6.0493E-02  7.4516E-01  8.9825E-01
             9.4453E-01
 PARAMETER:  6.8046E-02 -1.0422E+00 -1.8068E-01  4.7562E-01 -3.5284E-01  1.1126E-01  1.1049E+00 -2.7052E+00 -1.9416E-01 -7.3078E-03
             4.2930E-02
 GRADIENT:   1.3215E+00  5.4725E-01  1.6381E+00 -8.6513E+00 -2.0759E+00  1.1983E-01  1.1987E+00 -2.4930E-02 -6.9057E-03 -2.6067E-02
            -1.1153E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1664.66209411535        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1997             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6857E-01  3.1940E-01  7.5451E-01  1.4554E+00  6.3624E-01  1.0113E+00  2.7254E+00  6.9106E-02  7.4541E-01  8.9782E-01
             9.4463E-01
 PARAMETER:  6.8068E-02 -1.0413E+00 -1.8169E-01  4.7531E-01 -3.5219E-01  1.1128E-01  1.1026E+00 -2.5721E+00 -1.9382E-01 -7.7907E-03
             4.3042E-02
 GRADIENT:   4.0196E+02  4.7448E+01  1.0446E+01  6.5350E+02  4.7018E+01  4.8011E+01  4.3869E+01  3.8415E-02  1.6927E+01  1.2040E+00
             9.0366E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1664.66556252403        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2183
 NPARAMETR:  9.6859E-01  3.1981E-01  7.5444E-01  1.4551E+00  6.3638E-01  1.0114E+00  2.7213E+00  7.6516E-02  7.4558E-01  8.9722E-01
             9.4460E-01
 PARAMETER:  6.8085E-02 -1.0400E+00 -1.8177E-01  4.7507E-01 -3.5196E-01  1.1130E-01  1.1011E+00 -2.4703E+00 -1.9360E-01 -8.4524E-03
             4.3008E-02
 GRADIENT:   1.3288E+00  2.8770E-02 -1.4365E+00 -8.7328E+00  1.2001E+00  1.2273E-01  9.2695E-01 -3.2990E-02  9.2214E-02  1.3793E-01
             1.3541E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1664.66934643196        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2369             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6861E-01  3.2158E-01  7.5568E-01  1.4547E+00  6.3593E-01  1.0114E+00  2.7203E+00  8.8479E-02  7.4557E-01  8.9566E-01
             9.4416E-01
 PARAMETER:  6.8104E-02 -1.0345E+00 -1.8014E-01  4.7479E-01 -3.5267E-01  1.1132E-01  1.1007E+00 -2.3250E+00 -1.9360E-01 -1.0193E-02
             4.2536E-02
 GRADIENT:   4.0254E+02  4.8191E+01  1.4362E+01  6.5250E+02  4.2283E+01  4.8102E+01  4.3900E+01  4.8993E-02  1.6801E+01  8.9324E-01
             5.6820E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1664.67327659567        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2548
 NPARAMETR:  9.6810E-01  3.2038E-01  7.5466E-01  1.4564E+00  6.3643E-01  1.0113E+00  2.7093E+00  9.6990E-02  7.4581E-01  8.9589E-01
             9.4432E-01
 PARAMETER:  6.7584E-02 -1.0383E+00 -1.8149E-01  4.7595E-01 -3.5189E-01  1.1120E-01  1.0967E+00 -2.2331E+00 -1.9329E-01 -9.9405E-03
             4.2712E-02
 GRADIENT:   2.0284E-01  1.9130E-01 -1.9711E+00 -5.2863E+00  1.1318E+00  6.9218E-02  4.5828E-01 -4.7468E-02  2.3897E-02  1.0376E-01
             8.1713E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1664.68107099519        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2731
 NPARAMETR:  9.6810E-01  3.2024E-01  7.5578E-01  1.4558E+00  6.3651E-01  1.0111E+00  2.7090E+00  1.0692E-01  7.4558E-01  8.9500E-01
             9.4414E-01
 PARAMETER:  6.7584E-02 -1.0387E+00 -1.8000E-01  4.7557E-01 -3.5175E-01  1.1102E-01  1.0966E+00 -2.1357E+00 -1.9359E-01 -1.0929E-02
             4.2518E-02
 GRADIENT:   2.5292E-01  1.6694E-01  1.5822E-02 -7.3913E+00 -8.8881E-01  7.2537E-03  4.1385E-01 -6.0627E-02 -6.6941E-02 -7.3215E-02
            -6.7265E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1664.69698433500        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2914
 NPARAMETR:  9.6844E-01  3.2044E-01  7.5658E-01  1.4554E+00  6.3723E-01  1.0113E+00  2.7047E+00  1.3349E-01  7.4505E-01  8.9209E-01
             9.4345E-01
 PARAMETER:  6.7930E-02 -1.0381E+00 -1.7894E-01  4.7528E-01 -3.5062E-01  1.1125E-01  1.0950E+00 -1.9137E+00 -1.9430E-01 -1.4190E-02
             4.1789E-02
 GRADIENT:   1.0197E+00 -8.4417E-02 -9.0017E-01 -8.2104E+00  4.1058E-02  9.8340E-02  2.6441E-01 -7.8011E-02 -2.8794E-01 -2.6221E-01
            -2.5467E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1664.72260122013        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3097
 NPARAMETR:  9.6748E-01  3.1989E-01  7.5811E-01  1.4570E+00  6.3708E-01  1.0107E+00  2.6957E+00  2.0468E-01  7.4597E-01  8.8535E-01
             9.4254E-01
 PARAMETER:  6.6943E-02 -1.0398E+00 -1.7693E-01  4.7639E-01 -3.5086E-01  1.1065E-01  1.0917E+00 -1.4863E+00 -1.9307E-01 -2.1777E-02
             4.0825E-02
 GRADIENT:  -1.0355E+00  6.8563E-02 -2.3370E+00 -6.3106E+00 -8.9176E-01 -1.3651E-01 -2.4466E-01 -1.7228E-02 -1.1705E-01  2.1727E-01
             7.4839E-03

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1664.72520712452        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     3272
 NPARAMETR:  9.6562E-01  3.1478E-01  7.6265E-01  1.4597E+00  6.3830E-01  1.0095E+00  2.7103E+00  2.3492E-01  7.4539E-01  8.8234E-01
             9.4217E-01
 PARAMETER:  6.5017E-02 -1.0559E+00 -1.7096E-01  4.7824E-01 -3.4895E-01  1.0941E-01  1.0971E+00 -1.3485E+00 -1.9385E-01 -2.5182E-02
             4.0430E-02
 GRADIENT:  -4.8068E+00 -3.0557E-01 -1.2787E+00 -8.0845E+00 -1.8898E+00 -5.6958E-01 -7.2662E-01 -4.6004E-03 -1.4211E-01 -2.4772E-01
            -1.1771E-01

0ITERATION NO.:   94    OBJECTIVE VALUE:  -1664.75165598619        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:     3404
 NPARAMETR:  9.6807E-01  3.1471E-01  7.6459E-01  1.4605E+00  6.3920E-01  1.0109E+00  2.7525E+00  2.4211E-01  7.4535E-01  8.8369E-01
             9.4225E-01
 PARAMETER:  6.7546E-02 -1.0561E+00 -1.6841E-01  4.7875E-01 -3.4754E-01  1.1084E-01  1.1125E+00 -1.3184E+00 -1.9390E-01 -2.3652E-02
             4.0513E-02
 GRADIENT:  -5.4181E-01  4.2602E-01 -7.9481E-01  1.0993E+00 -2.6260E+00 -9.4914E-02  1.1285E-01  5.5005E-03  1.3834E-01  3.6531E-02
            -1.5979E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3404
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2342E-03  3.5046E-02 -1.3221E-02 -2.5283E-02  5.1594E-03
 SE:             2.9870E-02  1.7974E-02  6.1825E-03  2.5727E-02  2.3805E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6704E-01  5.1198E-02  3.2484E-02  3.2573E-01  8.2841E-01

 ETASHRINKSD(%)  1.0000E-10  3.9785E+01  7.9288E+01  1.3811E+01  2.0251E+01
 ETASHRINKVR(%)  1.0000E-10  6.3741E+01  9.5710E+01  2.5715E+01  3.6400E+01
 EBVSHRINKSD(%)  4.0082E-01  4.7042E+01  8.0376E+01  1.0837E+01  1.5683E+01
 EBVSHRINKVR(%)  8.0004E-01  7.1955E+01  9.6149E+01  2.0499E+01  2.8906E+01
 RELATIVEINF(%)  9.8408E+01  6.1422E+00  4.0034E-01  2.2609E+01  6.3635E+00
 EPSSHRINKSD(%)  4.4455E+01
 EPSSHRINKVR(%)  6.9148E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1664.7516559861926     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -929.60082942245447     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    50.61
 Elapsed covariance  time in seconds:     6.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1664.752       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  3.15E-01  7.65E-01  1.46E+00  6.39E-01  1.01E+00  2.75E+00  2.42E-01  7.45E-01  8.84E-01  9.42E-01
 


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
 
         2.99E-02  1.37E-01  1.75E-01  8.13E-02  1.17E-01  6.06E-02  7.75E-01  4.99E-01  6.78E-02  1.32E-01  6.57E-02
 


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
+        8.92E-04
 
 TH 2
+        2.50E-04  1.87E-02
 
 TH 3
+       -8.28E-04 -1.27E-03  3.05E-02
 
 TH 4
+       -3.35E-04 -9.08E-03  5.87E-03  6.61E-03
 
 TH 5
+       -4.79E-04  3.38E-03  1.94E-02  1.64E-03  1.36E-02
 
 TH 6
+        3.35E-04  4.87E-04 -8.54E-04 -2.18E-04 -5.45E-04  3.67E-03
 
 TH 7
+       -2.91E-03 -9.60E-02  2.67E-02  5.14E-02 -4.02E-03 -4.57E-03  6.01E-01
 
 TH 8
+       -2.33E-03 -8.61E-03  7.63E-02  1.72E-02  4.69E-02 -2.82E-03  9.58E-02  2.49E-01
 
 TH 9
+        3.05E-04  4.78E-03 -1.49E-03 -2.12E-03  2.45E-05  2.30E-04 -2.71E-02 -2.37E-03  4.59E-03
 
 TH10
+       -2.30E-04 -1.57E-03  1.05E-02  2.38E-03  6.59E-03 -6.78E-04  1.63E-02  1.34E-02 -1.63E-03  1.73E-02
 
 TH11
+        4.05E-04  1.09E-03  6.60E-04 -1.56E-04  8.70E-04 -2.52E-04 -5.77E-03  1.73E-03  7.79E-04 -7.44E-04  4.32E-03
 
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
+        2.99E-02
 
 TH 2
+        6.12E-02  1.37E-01
 
 TH 3
+       -1.59E-01 -5.29E-02  1.75E-01
 
 TH 4
+       -1.38E-01 -8.16E-01  4.13E-01  8.13E-02
 
 TH 5
+       -1.37E-01  2.11E-01  9.50E-01  1.72E-01  1.17E-01
 
 TH 6
+        1.85E-01  5.87E-02 -8.07E-02 -4.43E-02 -7.70E-02  6.06E-02
 
 TH 7
+       -1.26E-01 -9.05E-01  1.97E-01  8.16E-01 -4.44E-02 -9.73E-02  7.75E-01
 
 TH 8
+       -1.56E-01 -1.26E-01  8.74E-01  4.23E-01  8.03E-01 -9.31E-02  2.48E-01  4.99E-01
 
 TH 9
+        1.51E-01  5.16E-01 -1.26E-01 -3.84E-01  3.09E-03  5.59E-02 -5.17E-01 -6.99E-02  6.78E-02
 
 TH10
+       -5.86E-02 -8.70E-02  4.58E-01  2.23E-01  4.29E-01 -8.51E-02  1.60E-01  2.04E-01 -1.83E-01  1.32E-01
 
 TH11
+        2.06E-01  1.21E-01  5.75E-02 -2.91E-02  1.13E-01 -6.34E-02 -1.13E-01  5.26E-02  1.75E-01 -8.61E-02  6.57E-02
 
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
+        1.30E+03
 
 TH 2
+        1.06E+02  5.94E+02
 
 TH 3
+       -4.45E+01  2.23E+02  1.44E+03
 
 TH 4
+        1.33E+02  3.13E+02 -3.34E+02  8.72E+02
 
 TH 5
+        7.08E+01 -5.20E+02 -1.78E+03  2.58E+02  2.68E+03
 
 TH 6
+       -1.20E+02 -1.89E+01 -3.65E+01 -4.90E+01  4.31E+01  2.95E+02
 
 TH 7
+        8.64E+00  4.83E+01  4.34E+00 -1.52E+01 -8.55E+00  3.98E+00  1.10E+01
 
 TH 8
+        3.60E+00  9.12E+00 -7.51E+01  9.99E+00  1.14E+01  5.48E+00 -9.82E-01  2.36E+01
 
 TH 9
+       -7.54E+01 -1.12E+02  4.32E+01 -1.12E+02  3.55E+01  4.17E+00  7.37E+00 -1.62E+01  3.42E+02
 
 TH10
+       -2.00E+01  1.08E+01 -7.55E+01  6.40E+00 -2.21E+01  1.46E+01 -1.54E+00  2.24E+01  1.33E+00  9.78E+01
 
 TH11
+       -1.38E+02  4.11E+00  8.64E+01 -6.68E+01 -1.58E+02  3.32E+01  1.30E+00  3.27E+00 -2.73E+01  2.18E+01  2.71E+02
 
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
 #CPUT: Total CPU Time in Seconds,       57.146
Stop Time:
Wed Sep 29 14:08:28 CDT 2021
