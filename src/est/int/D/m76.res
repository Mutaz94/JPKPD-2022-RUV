Sat Sep 25 06:18:51 CDT 2021
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
$DATA ../../../../data/int/D/dat76.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m76.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   34313.6048990571        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.9892E+02  4.2253E+02 -4.9847E+01  2.3226E+02  1.7780E+02 -3.4655E+03 -1.5351E+03 -5.9321E+01 -2.2101E+03 -8.1356E+02
            -6.8061E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -830.039626180732        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  2.1779E+00  1.3813E+00  7.2116E-01  2.2599E+00  1.1493E+00  7.1600E+00  4.2440E+00  1.0451E+00  3.3302E+00  1.6238E+00
             1.1088E+01
 PARAMETER:  8.7836E-01  4.2306E-01 -2.2690E-01  9.1534E-01  2.3914E-01  2.0685E+00  1.5455E+00  1.4414E-01  1.3030E+00  5.8477E-01
             2.5058E+00
 GRADIENT:   2.3846E+01 -2.9353E+01 -5.9509E+01  1.1547E+02  3.2961E+01  1.5168E+02  6.0963E+00  3.9223E+00  2.7793E+01  3.2267E+01
             2.4459E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -890.836491236167        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  2.2134E+00  1.0316E+00  1.5006E+01  3.0886E+00  2.8115E+00  5.3683E+00  5.8217E+00  2.2112E+00  3.3144E+00  3.8753E+00
             1.0813E+01
 PARAMETER:  8.9455E-01  1.3114E-01  2.8085E+00  1.2277E+00  1.1337E+00  1.7805E+00  1.8616E+00  8.9353E-01  1.2983E+00  1.4546E+00
             2.4808E+00
 GRADIENT:   3.9547E+01  1.4292E+01 -2.0359E+01  8.0540E+01  1.9113E+01  1.0554E+02  3.0251E+00  1.8008E+00  7.7851E+00  8.6528E+01
             2.0637E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1050.28014468929        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1417E+00  1.2752E+00  3.1699E+01  1.5251E+00  2.4886E+00  2.8718E+00  6.1880E+00  1.3034E+00  2.2051E+00  1.1950E+00
             1.1155E+01
 PARAMETER:  2.3248E-01  3.4310E-01  3.5563E+00  5.2204E-01  1.0117E+00  1.1549E+00  1.9226E+00  3.6496E-01  8.9076E-01  2.7812E-01
             2.5119E+00
 GRADIENT:  -2.4594E+01  6.7573E+00 -5.6751E-01 -8.0272E+00  8.0975E+00  6.9130E+00  2.0787E+01  2.1359E-02  3.2633E+01  1.9684E+01
             2.8682E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1093.91667688620        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.1131E+00  6.7519E-01  6.6209E+01  1.4195E+00  2.4617E+00  2.8875E+00  6.5053E+00  2.3701E+00  1.3350E+00  5.9934E-01
             9.4068E+00
 PARAMETER:  2.0715E-01 -2.9277E-01  4.2928E+00  4.5030E-01  1.0008E+00  1.1604E+00  1.9726E+00  9.6291E-01  3.8897E-01 -4.1192E-01
             2.3414E+00
 GRADIENT:  -1.2780E+01 -9.5378E+00  3.4338E-01 -7.2315E+00 -5.7663E-02  8.7633E+00 -2.8726E+00  3.6092E-03  7.6158E+00  3.7324E+00
             2.9732E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1096.69114235695        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.1952E+00  1.3357E+00  2.6397E+01  1.1056E+00  2.4301E+00  2.7760E+00  5.5061E+00  2.7943E+00  9.2677E-01  3.4066E-01
             9.2596E+00
 PARAMETER:  2.7829E-01  3.8949E-01  3.3733E+00  2.0037E-01  9.8793E-01  1.1210E+00  1.8059E+00  1.1276E+00  2.3955E-02 -9.7688E-01
             2.3257E+00
 GRADIENT:   7.0711E+00  2.0469E+00  8.1872E-01 -8.1866E+00 -2.6844E+00 -5.0756E+00  3.4666E+00 -1.9119E-02  3.2849E+00  9.1687E-01
            -1.3999E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1097.91045689196        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      489
 NPARAMETR:  1.1749E+00  1.1869E+00  3.0292E+01  1.1561E+00  2.4605E+00  2.8618E+00  5.8626E+00  3.1059E+00  7.9009E-01  3.3583E-01
             9.3511E+00
 PARAMETER:  2.6120E-01  2.7134E-01  3.5109E+00  2.4504E-01  1.0004E+00  1.1515E+00  1.8686E+00  1.2333E+00 -1.3561E-01 -9.9115E-01
             2.3355E+00
 GRADIENT:   1.9915E-01  5.9148E-02  6.9487E-01 -1.8406E+00 -1.1971E+00  6.7787E-01 -5.7848E+00 -1.9470E-02 -7.4424E-01  8.1913E-01
            -5.0556E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1099.00752124845        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  1.1660E+00  1.0595E+00  1.2448E+01  1.2093E+00  2.3609E+00  2.8766E+00  6.1787E+00  8.4550E-01  8.6642E-01  7.4180E-02
             9.4288E+00
 PARAMETER:  2.5360E-01  1.5784E-01  2.6216E+00  2.9001E-01  9.5902E-01  1.1566E+00  1.9211E+00 -6.7831E-02 -4.3384E-02 -2.5013E+00
             2.3438E+00
 GRADIENT:  -2.3837E+00 -1.4104E+00 -6.7129E-01 -3.6796E+00  9.4352E+00  2.7327E+00 -7.2192E-01  2.6978E-02  8.2809E-01  4.1885E-02
             1.2812E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1099.22761530290        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      842
 NPARAMETR:  1.1743E+00  1.0720E+00  9.9877E+00  1.2055E+00  2.2503E+00  2.8504E+00  6.2015E+00  6.9632E-01  8.5410E-01  5.7178E-02
             9.3802E+00
 PARAMETER:  2.6069E-01  1.6953E-01  2.4014E+00  2.8690E-01  9.1107E-01  1.1475E+00  1.9248E+00 -2.6195E-01 -5.7708E-02 -2.7616E+00
             2.3386E+00
 GRADIENT:  -9.9046E-02 -5.8886E-02 -7.9094E-02 -4.0871E-02  1.8696E-01 -5.0221E-02  6.0995E-02  2.0368E-02  2.4946E-02  2.1950E-02
             9.1382E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1099.24762984299        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1022
 NPARAMETR:  1.1751E+00  1.0762E+00  1.0009E+01  1.2038E+00  2.2522E+00  2.8506E+00  6.1924E+00  1.1102E-01  8.5051E-01  2.0687E-02
             9.3815E+00
 PARAMETER:  2.6131E-01  1.7344E-01  2.4035E+00  2.8550E-01  9.1190E-01  1.1475E+00  1.9233E+00 -2.0980E+00 -6.1916E-02 -3.7782E+00
             2.3387E+00
 GRADIENT:   6.1068E-02  3.7054E-03  1.0942E-02  5.9396E-02 -9.3951E-02 -1.9792E-02 -3.6432E-02  5.8080E-04 -2.1337E-02  2.7832E-03
            -5.1749E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1099.24901587174        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  1.1748E+00  1.0762E+00  1.0022E+01  1.2040E+00  2.2528E+00  2.8513E+00  6.1957E+00  2.2221E-02  8.5117E-01  1.0000E-02
             9.3813E+00
 PARAMETER:  2.6111E-01  1.7344E-01  2.4048E+00  2.8561E-01  9.1217E-01  1.1478E+00  1.9239E+00 -3.7067E+00 -6.1139E-02 -4.6727E+00
             2.3387E+00
 GRADIENT:   1.3113E-02  2.3189E-02  4.7003E-03 -4.2150E-02 -3.3428E-02  5.7183E-02  9.0590E-02  2.3225E-05  4.5092E-03  0.0000E+00
            -1.7451E-03

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1099.24902846359        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1294
 NPARAMETR:  1.1748E+00  1.0755E+00  1.0020E+01  1.2042E+00  2.2529E+00  2.8511E+00  6.1959E+00  2.2192E-02  8.5111E-01  1.0000E-02
             9.3814E+00
 PARAMETER:  2.6110E-01  1.7275E-01  2.4046E+00  2.8578E-01  9.1222E-01  1.1477E+00  1.9239E+00 -3.7080E+00 -6.1213E-02 -4.6727E+00
             2.3387E+00
 GRADIENT:   7.8788E-03  8.1721E-04 -2.5079E-03 -7.7321E-03  1.0039E-02  2.8650E-02  5.0722E-02  2.3404E-05 -5.2001E-03  0.0000E+00
            -1.3771E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1294
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2216E-02  2.9718E-02 -2.6209E-05 -6.5443E-02  3.6598E-05
 SE:             2.8753E-02  2.4966E-02  3.7069E-05  1.2432E-02  1.3406E-04
 N:                     100         100         100         100         100

 P VAL.:         6.7095E-01  2.3392E-01  4.7954E-01  1.4134E-07  7.8486E-01

 ETASHRINKSD(%)  3.6725E+00  1.6361E+01  9.9876E+01  5.8350E+01  9.9551E+01
 ETASHRINKVR(%)  7.2101E+00  3.0044E+01  1.0000E+02  8.2653E+01  9.9998E+01
 EBVSHRINKSD(%)  3.4271E+00  1.0198E+01  9.9862E+01  6.5243E+01  9.9527E+01
 EBVSHRINKVR(%)  6.7368E+00  1.9356E+01  1.0000E+02  8.7920E+01  9.9998E+01
 RELATIVEINF(%)  9.3192E+01  4.3991E+01  4.0280E-05  6.6869E+00  4.7324E-04
 EPSSHRINKSD(%)  5.7413E+00
 EPSSHRINKVR(%)  1.1153E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1099.2490284635871     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       554.84033130482362     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.65
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1099.249       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.17E+00  1.08E+00  1.00E+01  1.20E+00  2.25E+00  2.85E+00  6.20E+00  2.22E-02  8.51E-01  1.00E-02  9.38E+00
 


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
+        9.27E+01
 
 TH 2
+       -5.56E-01  2.71E+01
 
 TH 3
+        8.07E-02  1.20E-01  5.73E-02
 
 TH 4
+       -3.70E+00  2.57E+01 -1.62E-01  1.78E+02
 
 TH 5
+       -1.42E+00 -5.04E+00 -1.38E+00 -5.63E+00  4.44E+01
 
 TH 6
+        2.91E+00 -4.58E-01  1.73E-02 -1.21E-02 -1.29E+00  2.26E+01
 
 TH 7
+       -3.02E-02  2.45E+00 -4.46E-02 -1.09E+01  1.11E+00 -1.83E-01  3.11E+00
 
 TH 8
+        3.72E+00  9.11E+00  1.45E-02  8.77E+00  2.31E-01 -5.01E-01  3.20E-01  1.02E+01
 
 TH 9
+       -5.16E+00 -2.91E+00 -2.43E-01 -4.06E+01  7.37E+00 -2.66E-01  2.29E+00 -2.41E+01  1.71E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.31E+00 -2.02E+00 -1.66E-02 -1.30E+01  9.41E-01  3.46E-01  6.45E-01 -1.54E-01  4.24E+00  0.00E+00  1.16E+01
 
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
 #CPUT: Total CPU Time in Seconds,       52.940
Stop Time:
Sat Sep 25 06:19:46 CDT 2021
