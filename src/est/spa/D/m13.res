Wed Sep 29 19:43:58 CDT 2021
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
$DATA ../../../../data/spa/D/dat13.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m13.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1161.85713048040        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.2182E+02 -1.8329E+02 -4.5601E+01 -1.7830E+02  1.1575E+02 -2.0330E+02 -3.0401E+02 -3.6825E+01 -3.1542E+02 -5.1173E+01
            -6.0855E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1424.46382049977        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  8.9334E-01  1.2340E+00  1.0318E+00  9.5338E-01  1.0843E+00  1.1829E+00  2.0107E+00  1.1615E+00  1.0851E+00  8.7230E-01
             1.0432E+00
 PARAMETER: -1.2786E-02  3.1026E-01  1.3135E-01  5.2256E-02  1.8098E-01  2.6795E-01  7.9848E-01  2.4967E-01  1.8166E-01 -3.6618E-02
             1.4233E-01
 GRADIENT:   6.6957E+01  1.0718E+02  2.0643E+01  7.6253E+01 -2.4058E+01  1.7162E+02  9.0945E+01 -1.7983E+01 -3.3706E+00  6.4901E+00
            -8.0685E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1428.21184824838        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      181
 NPARAMETR:  8.9208E-01  1.2481E+00  9.1879E-01  9.2774E-01  1.0464E+00  1.1901E+00  1.9867E+00  1.4310E+00  1.0987E+00  7.4037E-01
             1.0489E+00
 PARAMETER: -1.4200E-02  3.2162E-01  1.5301E-02  2.4993E-02  1.4538E-01  2.7407E-01  7.8647E-01  4.5835E-01  1.9409E-01 -2.0060E-01
             1.4775E-01
 GRADIENT:  -3.1408E+02 -4.6685E+01  3.3501E-01  2.2344E+01 -2.2965E+01 -2.2468E+02 -1.0580E+02 -1.0382E+00 -8.4786E+00  3.0348E+00
            -3.9540E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1473.23679136085        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      361
 NPARAMETR:  1.0989E+00  1.1398E+00  1.2512E+00  9.5155E-01  1.1708E+00  1.2547E+00  2.3215E+00  2.1115E+00  1.0024E+00  8.9301E-01
             1.2025E+00
 PARAMETER:  1.9434E-01  2.3085E-01  3.2409E-01  5.0339E-02  2.5769E-01  3.2690E-01  9.4221E-01  8.4738E-01  1.0240E-01 -1.3154E-02
             2.8441E-01
 GRADIENT:   1.3896E+01 -5.4327E+01 -9.6680E+00 -2.5107E+00 -2.4546E+01 -1.2604E+02 -6.7092E+01  1.4590E+01 -3.3584E+00  9.1048E+00
             4.6466E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1501.86187062362        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      538
 NPARAMETR:  1.0707E+00  1.5571E+00  1.2013E+00  7.9549E-01  1.2847E+00  1.5655E+00  2.5131E+00  1.9649E+00  1.1091E+00  9.5187E-01
             1.0637E+00
 PARAMETER:  1.6832E-01  5.4283E-01  2.8344E-01 -1.2879E-01  3.5050E-01  5.4819E-01  1.0215E+00  7.7545E-01  2.0356E-01  5.0668E-02
             1.6175E-01
 GRADIENT:  -1.3632E+01  1.2495E-01  1.9366E+00 -6.7197E+00  1.3247E+00  1.1052E+00  9.7437E+00 -8.0878E-01  1.5210E+00 -3.1069E+00
             4.3069E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1502.22466915610        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  1.0880E+00  1.6152E+00  1.1165E+00  7.6057E-01  1.2898E+00  1.5713E+00  2.3635E+00  1.9636E+00  1.0809E+00  9.8306E-01
             1.0632E+00
 PARAMETER:  1.8433E-01  5.7946E-01  2.1017E-01 -1.7369E-01  3.5446E-01  5.5190E-01  9.6014E-01  7.7480E-01  1.7778E-01  8.2918E-02
             1.6132E-01
 GRADIENT:   1.3005E+00 -1.6669E+00  6.5894E-01 -1.3365E+00 -7.0549E-01  2.6676E+00 -9.9536E-01 -2.3950E-01 -4.6763E-01  2.5807E-02
             8.7546E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1502.42648482292        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      858
 NPARAMETR:  1.0921E+00  1.6174E+00  1.1128E+00  7.6299E-01  1.2895E+00  1.6272E+00  2.3714E+00  1.9597E+00  1.0995E+00  9.8239E-01
             1.0594E+00
 PARAMETER:  1.8812E-01  5.8082E-01  2.0687E-01 -1.7050E-01  3.5424E-01  5.8683E-01  9.6349E-01  7.7280E-01  1.9485E-01  8.2236E-02
             1.5769E-01
 GRADIENT:   4.9509E+00 -7.6104E-01 -4.7061E-01  5.3968E-01  1.1360E+00  1.8758E+01  2.3417E-01  2.4124E-01  2.7591E-01  6.1586E-03
            -1.6141E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1502.44021110791        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1042
 NPARAMETR:  1.0929E+00  1.6082E+00  1.1138E+00  7.6711E-01  1.2884E+00  1.6163E+00  2.3817E+00  1.9588E+00  1.0904E+00  9.8231E-01
             1.0600E+00
 PARAMETER:  1.8881E-01  5.7512E-01  2.0782E-01 -1.6512E-01  3.5341E-01  5.8013E-01  9.6781E-01  7.7233E-01  1.8658E-01  8.2154E-02
             1.5830E-01
 GRADIENT:   5.5793E+00 -1.1730E+00 -1.0885E+00  1.0078E+00  2.2502E+00  1.5639E+01  3.9588E-01  4.2025E-01  2.0697E-01  6.1780E-02
             8.4403E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1502.45179411815        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1231             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0924E+00  1.5999E+00  1.1245E+00  7.6808E-01  1.2859E+00  1.6200E+00  2.3872E+00  1.9468E+00  1.0837E+00  9.8102E-01
             1.0599E+00
 PARAMETER:  1.8837E-01  5.6995E-01  2.1733E-01 -1.6387E-01  3.5146E-01  5.8243E-01  9.7012E-01  7.6618E-01  1.8042E-01  8.0834E-02
             1.5815E-01
 GRADIENT:   7.1346E+02  3.6114E+02  2.1451E+00  6.1159E+01  1.3692E+01  6.5921E+02  2.6860E+02  2.1639E+00  5.9832E+00  4.6038E-01
             1.1149E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1502.45716693572        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1418
 NPARAMETR:  1.0923E+00  1.5962E+00  1.1271E+00  7.7019E-01  1.2851E+00  1.6204E+00  2.3913E+00  1.9453E+00  1.0825E+00  9.8020E-01
             1.0598E+00
 PARAMETER:  1.8830E-01  5.6765E-01  2.1961E-01 -1.6112E-01  3.5085E-01  5.8268E-01  9.7183E-01  7.6539E-01  1.7930E-01  8.0001E-02
             1.5813E-01
 GRADIENT:   5.1958E+00 -1.9194E+00  8.0507E-01 -2.2015E+00 -1.0395E+00  1.6846E+01  2.2079E-01 -4.0007E-01 -2.1957E-01 -7.2780E-03
            -1.6115E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1502.46453069707        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1605             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0923E+00  1.5912E+00  1.1241E+00  7.7512E-01  1.2850E+00  1.6205E+00  2.3975E+00  1.9542E+00  1.0871E+00  9.7966E-01
             1.0600E+00
 PARAMETER:  1.8825E-01  5.6447E-01  2.1696E-01 -1.5474E-01  3.5074E-01  5.8275E-01  9.7442E-01  7.7000E-01  1.8348E-01  7.9454E-02
             1.5831E-01
 GRADIENT:   7.1237E+02  3.5418E+02  2.4672E-01  6.1329E+01  1.6612E+01  6.5910E+02  2.6938E+02  3.0534E+00  6.6846E+00  5.3230E-01
             1.4278E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1502.46933149879        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1789
 NPARAMETR:  1.0922E+00  1.5878E+00  1.1270E+00  7.7694E-01  1.2843E+00  1.6206E+00  2.4013E+00  1.9531E+00  1.0857E+00  9.7919E-01
             1.0600E+00
 PARAMETER:  1.8824E-01  5.6238E-01  2.1955E-01 -1.5239E-01  3.5025E-01  5.8277E-01  9.7599E-01  7.6941E-01  1.8226E-01  7.8967E-02
             1.5830E-01
 GRADIENT:   5.1312E+00 -1.7598E+00 -1.0357E+00  5.8080E-01  1.7411E+00  1.6835E+01  4.0681E-01  4.2556E-01  2.6145E-01  8.2678E-02
             1.4997E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1502.47654508368        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1976             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0922E+00  1.5830E+00  1.1366E+00  7.7807E-01  1.2828E+00  1.6206E+00  2.4064E+00  1.9446E+00  1.0791E+00  9.7839E-01
             1.0598E+00
 PARAMETER:  1.8822E-01  5.5930E-01  2.2806E-01 -1.5093E-01  3.4906E-01  5.8280E-01  9.7815E-01  7.6503E-01  1.7616E-01  7.8157E-02
             1.5812E-01
 GRADIENT:   7.1267E+02  3.4840E+02  1.9085E+00  5.8181E+01  1.3956E+01  6.5955E+02  2.7101E+02  2.3333E+00  6.0021E+00  4.3959E-01
             1.1568E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1502.48105826012        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2165
 NPARAMETR:  1.0922E+00  1.5800E+00  1.1384E+00  7.7998E-01  1.2824E+00  1.6206E+00  2.4099E+00  1.9460E+00  1.0789E+00  9.7815E-01
             1.0599E+00
 PARAMETER:  1.8820E-01  5.5740E-01  2.2965E-01 -1.4848E-01  3.4875E-01  5.8280E-01  9.7960E-01  7.6577E-01  1.7595E-01  7.7911E-02
             1.5813E-01
 GRADIENT:   5.1614E+00 -1.8679E+00  2.6700E-01 -1.3149E+00 -4.5480E-01  1.6862E+01  3.4560E-01 -1.4483E-01 -7.6267E-02  2.5974E-02
            -6.5802E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1502.48605588471        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2354             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0922E+00  1.5747E+00  1.1390E+00  7.8432E-01  1.2822E+00  1.6206E+00  2.4163E+00  1.9503E+00  1.0803E+00  9.7761E-01
             1.0599E+00
 PARAMETER:  1.8817E-01  5.5407E-01  2.3014E-01 -1.4294E-01  3.4858E-01  5.8281E-01  9.8225E-01  7.6796E-01  1.7725E-01  7.7359E-02
             1.5822E-01
 GRADIENT:   7.1204E+02  3.4201E+02  5.9857E-01  5.7797E+01  1.6061E+01  6.5911E+02  2.7192E+02  2.9450E+00  6.4361E+00  4.8086E-01
             1.3614E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1502.49003982363        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2540
 NPARAMETR:  1.0922E+00  1.5720E+00  1.1425E+00  7.8561E-01  1.2817E+00  1.6206E+00  2.4195E+00  1.9486E+00  1.0792E+00  9.7737E-01
             1.0599E+00
 PARAMETER:  1.8816E-01  5.5233E-01  2.3318E-01 -1.4130E-01  3.4820E-01  5.8282E-01  9.8354E-01  7.6710E-01  1.7619E-01  7.7109E-02
             1.5821E-01
 GRADIENT:   5.1313E+00 -1.7495E+00 -5.9618E-01  7.5112E-02  9.7703E-01  1.6837E+01  4.3719E-01  2.0027E-01  1.3313E-01  3.4102E-02
             6.7644E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1502.49470573476        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2727             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0921E+00  1.5668E+00  1.1512E+00  7.8739E-01  1.2806E+00  1.6207E+00  2.4252E+00  1.9450E+00  1.0745E+00  9.7692E-01
             1.0598E+00
 PARAMETER:  1.8814E-01  5.4904E-01  2.4077E-01 -1.3903E-01  3.4729E-01  5.8284E-01  9.8591E-01  7.6525E-01  1.7190E-01  7.6647E-02
             1.5811E-01
 GRADIENT:   7.1220E+02  3.3647E+02  1.9246E+00  5.4946E+01  1.3847E+01  6.5941E+02  2.7349E+02  2.4256E+00  5.9635E+00  4.4499E-01
             1.1828E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1502.49789222890        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2916
 NPARAMETR:  1.0921E+00  1.5644E+00  1.1525E+00  7.8904E-01  1.2804E+00  1.6207E+00  2.4280E+00  1.9462E+00  1.0747E+00  9.7676E-01
             1.0598E+00
 PARAMETER:  1.8812E-01  5.4751E-01  2.4192E-01 -1.3694E-01  3.4715E-01  5.8285E-01  9.8708E-01  7.6590E-01  1.7202E-01  7.6488E-02
             1.5812E-01
 GRADIENT:   5.1494E+00 -1.8337E+00  1.5758E-01 -1.0259E+00 -3.7722E-01  1.6851E+01  4.3183E-01 -1.0357E-01 -4.4557E-02  1.9205E-02
            -4.3560E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1502.50189253914        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3101             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0921E+00  1.5600E+00  1.1546E+00  7.9246E-01  1.2803E+00  1.6207E+00  2.4334E+00  1.9494E+00  1.0754E+00  9.7650E-01
             1.0599E+00
 PARAMETER:  1.8809E-01  5.4469E-01  2.4374E-01 -1.3261E-01  3.4711E-01  5.8285E-01  9.8929E-01  7.6751E-01  1.7268E-01  7.6216E-02
             1.5820E-01
 GRADIENT:   7.1171E+02  3.3130E+02  1.0179E+00  5.4393E+01  1.5410E+01  6.5908E+02  2.7426E+02  2.8436E+00  6.2759E+00  4.5670E-01
             1.3254E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1502.50401576514        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3285
 NPARAMETR:  1.0921E+00  1.5569E+00  1.1571E+00  7.9441E-01  1.2801E+00  1.6207E+00  2.4371E+00  1.9506E+00  1.0746E+00  9.7628E-01
             1.0599E+00
 PARAMETER:  1.8807E-01  5.4267E-01  2.4588E-01 -1.3015E-01  3.4694E-01  5.8286E-01  9.9083E-01  7.6814E-01  1.7194E-01  7.5990E-02
             1.5821E-01
 GRADIENT:   5.1190E+00 -1.7199E+00 -6.9267E-01  3.4574E-01  1.0935E+00  1.6828E+01  5.0539E-01  2.4190E-01  1.4047E-01  9.9830E-03
             7.7851E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1502.50707407242        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3474             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0921E+00  1.5526E+00  1.1651E+00  7.9530E-01  1.2790E+00  1.6207E+00  2.4417E+00  1.9469E+00  1.0708E+00  9.7608E-01
             1.0598E+00
 PARAMETER:  1.8806E-01  5.3996E-01  2.5285E-01 -1.2904E-01  3.4608E-01  5.8288E-01  9.9271E-01  7.6622E-01  1.6839E-01  7.5787E-02
             1.5808E-01
 GRADIENT:   7.1185E+02  3.2624E+02  2.0303E+00  5.1912E+01  1.3651E+01  6.5935E+02  2.7583E+02  2.4855E+00  5.9269E+00  4.4975E-01
             1.1853E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1502.50958099880        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     3665
 NPARAMETR:  1.0920E+00  1.5509E+00  1.1662E+00  7.9689E-01  1.2790E+00  1.6207E+00  2.4439E+00  1.9479E+00  1.0708E+00  9.7597E-01
             1.0598E+00
 PARAMETER:  1.8805E-01  5.3886E-01  2.5377E-01 -1.2704E-01  3.4604E-01  5.8289E-01  9.9359E-01  7.6677E-01  1.6844E-01  7.5676E-02
             1.5811E-01
 GRADIENT:   5.1400E+00 -1.7949E+00  1.5691E-01 -8.8163E-01 -4.3774E-01  1.6844E+01  4.8414E-01 -1.0028E-01 -4.4623E-02  1.3027E-02
            -4.3375E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1502.51227543322        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3852             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0920E+00  1.5474E+00  1.1681E+00  7.9964E-01  1.2791E+00  1.6208E+00  2.4482E+00  1.9512E+00  1.0715E+00  9.7583E-01
             1.0599E+00
 PARAMETER:  1.8802E-01  5.3659E-01  2.5541E-01 -1.2359E-01  3.4612E-01  5.8289E-01  9.9536E-01  7.6846E-01  1.6906E-01  7.5529E-02
             1.5818E-01
 GRADIENT:   7.1140E+02  3.2224E+02  1.2139E+00  5.1543E+01  1.5128E+01  6.5904E+02  2.7628E+02  2.8509E+00  6.1798E+00  4.4332E-01
             1.3093E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1502.51334692160        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     4036
 NPARAMETR:  1.0920E+00  1.5448E+00  1.1700E+00  8.0144E-01  1.2790E+00  1.6208E+00  2.4514E+00  1.9527E+00  1.0711E+00  9.7571E-01
             1.0599E+00
 PARAMETER:  1.8800E-01  5.3493E-01  2.5703E-01 -1.2135E-01  3.4609E-01  5.8289E-01  9.9664E-01  7.6922E-01  1.6867E-01  7.5414E-02
             1.5820E-01
 GRADIENT:   5.1142E+00 -1.6748E+00 -6.5325E-01  4.3725E-01  1.0011E+00  1.6824E+01  5.3799E-01  2.2912E-01  1.2525E-01 -2.6919E-03
             7.0664E-02

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1502.51555974235        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     4225             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0920E+00  1.5414E+00  1.1771E+00  8.0196E-01  1.2781E+00  1.6208E+00  2.4551E+00  1.9495E+00  1.0680E+00  9.7567E-01
             1.0598E+00
 PARAMETER:  1.8799E-01  5.3269E-01  2.6302E-01 -1.2069E-01  3.4538E-01  5.8291E-01  9.9817E-01  7.6758E-01  1.6577E-01  7.5373E-02
             1.5809E-01
 GRADIENT:   7.1152E+02  3.1815E+02  2.0355E+00  4.9488E+01  1.3710E+01  6.5927E+02  2.7758E+02  2.5609E+00  5.9117E+00  4.4704E-01
             1.1988E+00

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1502.51725362028        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     4416
 NPARAMETR:  1.0920E+00  1.5400E+00  1.1780E+00  8.0327E-01  1.2781E+00  1.6208E+00  2.4569E+00  1.9505E+00  1.0681E+00  9.7561E-01
             1.0598E+00
 PARAMETER:  1.8798E-01  5.3177E-01  2.6382E-01 -1.1907E-01  3.4538E-01  5.8292E-01  9.9890E-01  7.6810E-01  1.6586E-01  7.5304E-02
             1.5811E-01
 GRADIENT:   5.1330E+00 -1.7842E+00  1.3259E-01 -7.4816E-01 -4.2140E-01  1.6839E+01  5.3921E-01 -8.2883E-02 -3.2982E-02  9.9573E-03
            -3.5941E-02

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1502.51908145905        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     4603             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0919E+00  1.5372E+00  1.1799E+00  8.0541E-01  1.2782E+00  1.6208E+00  2.4603E+00  1.9533E+00  1.0684E+00  9.7554E-01
             1.0599E+00
 PARAMETER:  1.8796E-01  5.2999E-01  2.6539E-01 -1.1640E-01  3.4548E-01  5.8292E-01  1.0003E+00  7.6954E-01  1.6619E-01  7.5233E-02
             1.5817E-01
 GRADIENT:   7.1116E+02  3.1500E+02  1.4183E+00  4.9129E+01  1.4848E+01  6.5903E+02  2.7796E+02  2.8426E+00  6.0952E+00  4.3632E-01
             1.2919E+00

0ITERATION NO.:  132    OBJECTIVE VALUE:  -1502.51908145905        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     4662
 NPARAMETR:  1.0919E+00  1.5372E+00  1.1799E+00  8.0541E-01  1.2782E+00  1.6208E+00  2.4603E+00  1.9533E+00  1.0684E+00  9.7554E-01
             1.0599E+00
 PARAMETER:  1.8796E-01  5.2999E-01  2.6539E-01 -1.1640E-01  3.4548E-01  5.8292E-01  1.0003E+00  7.6954E-01  1.6619E-01  7.5233E-02
             1.5817E-01
 GRADIENT:   2.0932E-03  3.1740E-01 -1.7883E-01  4.4878E-03  3.8831E-01 -3.8031E-03 -1.8572E-01  6.0537E-02  2.8137E-02  2.0828E-04
             2.7619E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     4662
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.5162E-04  4.3940E-03 -5.1514E-02 -1.2003E-02 -4.3751E-02
 SE:             2.9924E-02  2.7264E-02  1.5480E-02  1.6012E-02  1.8805E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8796E-01  8.7197E-01  8.7585E-04  4.5349E-01  1.9989E-02

 ETASHRINKSD(%)  1.0000E-10  8.6613E+00  4.8139E+01  4.6357E+01  3.7001E+01
 ETASHRINKVR(%)  1.0000E-10  1.6572E+01  7.3104E+01  7.1224E+01  6.0312E+01
 EBVSHRINKSD(%)  1.8768E-01  7.0279E+00  5.4227E+01  5.0971E+01  3.3188E+01
 EBVSHRINKVR(%)  3.7502E-01  1.3562E+01  7.9048E+01  7.5962E+01  5.5361E+01
 RELATIVEINF(%)  9.9544E+01  2.1232E+01  4.0387E+00  3.7024E+00  1.6609E+01
 EPSSHRINKSD(%)  4.5209E+01
 EPSSHRINKVR(%)  6.9979E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1502.5190814590485     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -767.36825489531032     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    76.05
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1502.519       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.09E+00  1.54E+00  1.18E+00  8.05E-01  1.28E+00  1.62E+00  2.46E+00  1.95E+00  1.07E+00  9.76E-01  1.06E+00
 


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
+        3.53E+02
 
 TH 2
+       -7.36E-01  6.85E+01
 
 TH 3
+        1.34E+00  7.91E+00  5.30E+01
 
 TH 4
+       -8.50E-01  8.25E+01 -8.96E+01  3.93E+02
 
 TH 5
+       -1.49E+00 -1.96E+01 -8.33E+01  1.55E+02  3.37E+02
 
 TH 6
+       -1.79E-02  1.98E-03  2.56E-01 -6.55E-01 -4.11E-01  7.56E+01
 
 TH 7
+        2.85E-01  5.82E+00 -3.09E+00 -2.85E+01  4.44E+00 -2.85E-01  2.38E+01
 
 TH 8
+       -1.18E-01 -2.13E+00 -1.21E+01  1.51E+01 -3.58E+00  8.79E-03 -1.84E+04  8.39E+00
 
 TH 9
+       -4.20E-02 -5.31E+00 -8.93E+00  1.84E+00  1.93E+00 -2.80E-01  6.26E+00  4.07E+00  1.84E+01
 
 TH10
+        2.50E-01 -5.07E-01 -2.71E+00  1.20E-01 -6.69E+01 -2.06E-01 -7.05E-01  3.95E+00  4.97E+00  5.75E+01
 
 TH11
+       -2.13E+00 -2.27E+00 -6.28E+00  6.28E+00 -1.68E+01  1.40E+00  9.39E-01  2.77E+00  5.91E+00  1.66E+01  1.82E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       83.068
Stop Time:
Wed Sep 29 19:45:23 CDT 2021
