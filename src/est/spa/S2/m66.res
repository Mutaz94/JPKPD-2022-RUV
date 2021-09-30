Wed Sep 29 17:35:54 CDT 2021
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
$DATA ../../../../data/spa/S2/dat66.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1627.18027950955        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8318E+02 -4.5980E+01  1.3950E+00 -3.2171E+01  1.2983E+01  4.0587E+01 -7.0128E+00  2.7882E+00  2.5965E+01 -9.5487E+00
            -1.3793E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1635.02987585433        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.4864E-01  1.1813E+00  9.8992E-01  9.6762E-01  1.0741E+00  9.1665E-01  1.1163E+00  9.6884E-01  7.2950E-01  1.1400E+00
             1.1347E+00
 PARAMETER:  4.7276E-02  2.6665E-01  8.9873E-02  6.7084E-02  1.7150E-01  1.2967E-02  2.1003E-01  6.8349E-02 -2.1540E-01  2.3100E-01
             2.2636E-01
 GRADIENT:  -1.3929E+01  2.4243E+01  4.4155E+00  1.9269E+01 -1.0530E+01 -2.2960E+01 -7.2555E+00  5.4630E-01 -8.9546E+00  4.1291E+00
             3.0139E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1635.93794993487        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.6107E-01  9.3072E-01  1.1285E+00  1.1248E+00  1.0305E+00  9.5760E-01  1.5666E+00  8.9445E-01  5.5203E-01  1.1048E+00
             1.0617E+00
 PARAMETER:  6.0296E-02  2.8208E-02  2.2091E-01  2.1762E-01  1.3002E-01  5.6677E-02  5.4889E-01 -1.1542E-02 -4.9416E-01  1.9962E-01
             1.5990E-01
 GRADIENT:   2.5629E+01  3.0948E+01  1.4982E+01  2.0156E+01 -4.1409E+00 -4.4268E+00  4.5131E+00 -5.6553E+00 -1.3818E+01 -5.2136E+00
            -1.3754E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1639.52623942080        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.4894E-01  9.0420E-01  8.6598E-01  1.1200E+00  8.8102E-01  9.5717E-01  1.4133E+00  6.3547E-01  7.3840E-01  9.4702E-01
             1.0480E+00
 PARAMETER:  4.7593E-02 -7.0863E-04 -4.3889E-02  2.1336E-01 -2.6675E-02  5.6229E-02  4.4590E-01 -3.5340E-01 -2.0327E-01  4.5569E-02
             1.4687E-01
 GRADIENT:  -9.6614E+00  7.5985E+00 -1.0309E+01  1.9408E+01  9.5233E+00 -5.3374E+00  4.2772E-01  2.0911E+00  4.1276E+00  5.0039E-01
             3.9736E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1641.23650208367        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.5359E-01  6.6872E-01  7.7642E-01  1.2300E+00  7.2719E-01  9.7787E-01  1.8031E+00  2.7639E-01  6.5952E-01  8.4706E-01
             1.0394E+00
 PARAMETER:  5.2483E-02 -3.0240E-01 -1.5306E-01  3.0701E-01 -2.1857E-01  7.7624E-02  6.8949E-01 -1.1859E+00 -3.1624E-01 -6.5988E-02
             1.3865E-01
 GRADIENT:   4.4032E+00  4.4039E+00  1.3674E+01 -1.3292E+01 -1.8109E+01  3.7640E+00 -9.3443E-01  5.5130E-02 -1.4006E+00 -1.2724E+00
            -3.1416E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1641.69334551345        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      884             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5056E-01  5.2793E-01  7.8535E-01  1.3065E+00  6.9684E-01  9.6621E-01  2.1666E+00  1.2761E-01  6.3703E-01  8.6709E-01
             1.0495E+00
 PARAMETER:  4.9295E-02 -5.3880E-01 -1.4163E-01  3.6738E-01 -2.6119E-01  6.5622E-02  8.7318E-01 -1.9587E+00 -3.5095E-01 -4.2614E-02
             1.4827E-01
 GRADIENT:   3.2177E+02  4.2631E+01  6.5738E+00  3.3220E+02  2.4864E+01  3.2477E+01  3.8180E+01  1.0690E-01  1.2811E+01  8.0186E-01
             1.1502E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1641.69446456597        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1064             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5056E-01  5.2894E-01  7.8497E-01  1.3063E+00  6.9677E-01  9.6653E-01  2.1667E+00  1.3037E-01  6.3634E-01  8.6604E-01
             1.0495E+00
 PARAMETER:  4.9294E-02 -5.3687E-01 -1.4210E-01  3.6717E-01 -2.6130E-01  6.5952E-02  8.7321E-01 -1.9373E+00 -3.5202E-01 -4.3818E-02
             1.4829E-01
 GRADIENT:   3.2170E+02  4.2891E+01  6.7018E+00  3.3248E+02  2.4716E+01  3.2599E+01  3.8330E+01  1.0937E-01  1.2644E+01  6.8788E-01
             1.1224E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1641.69512697965        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1246             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5056E-01  5.2958E-01  7.8455E-01  1.3062E+00  6.9678E-01  9.6652E-01  2.1659E+00  1.3076E-01  6.3640E-01  8.6561E-01
             1.0494E+00
 PARAMETER:  4.9296E-02 -5.3567E-01 -1.4264E-01  3.6709E-01 -2.6128E-01  6.5950E-02  8.7283E-01 -1.9344E+00 -3.5193E-01 -4.4322E-02
             1.4825E-01
 GRADIENT:   3.2167E+02  4.2999E+01  6.3269E+00  3.3313E+02  2.5163E+01  3.2593E+01  3.8365E+01  1.1070E-01  1.2643E+01  6.6041E-01
             1.1155E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1641.69543743650        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1427             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5058E-01  5.3014E-01  7.8427E-01  1.3058E+00  6.9667E-01  9.6654E-01  2.1642E+00  1.2797E-01  6.3654E-01  8.6539E-01
             1.0494E+00
 PARAMETER:  4.9312E-02 -5.3461E-01 -1.4300E-01  3.6680E-01 -2.6144E-01  6.5966E-02  8.7205E-01 -1.9559E+00 -3.5171E-01 -4.4574E-02
             1.4824E-01
 GRADIENT:   3.2170E+02  4.3016E+01  6.6303E+00  3.3244E+02  2.4863E+01  3.2598E+01  3.8320E+01  1.0621E-01  1.2624E+01  6.2975E-01
             1.0909E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1641.69594263950        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1609             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5058E-01  5.3050E-01  7.8391E-01  1.3056E+00  6.9668E-01  9.6654E-01  2.1628E+00  1.2881E-01  6.3661E-01  8.6523E-01
             1.0494E+00
 PARAMETER:  4.9319E-02 -5.3394E-01 -1.4346E-01  3.6664E-01 -2.6142E-01  6.5972E-02  8.7142E-01 -1.9494E+00 -3.5159E-01 -4.4755E-02
             1.4825E-01
 GRADIENT:   3.2167E+02  4.2959E+01  6.2800E+00  3.3235E+02  2.5278E+01  3.2594E+01  3.8271E+01  1.0890E-01  1.2630E+01  6.6454E-01
             1.1282E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1641.69629711608        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1791
 NPARAMETR:  9.5059E-01  5.3090E-01  7.8359E-01  1.3053E+00  6.9668E-01  9.6655E-01  2.1615E+00  1.2735E-01  6.3667E-01  8.6504E-01
             1.0495E+00
 PARAMETER:  4.9325E-02 -5.3318E-01 -1.4387E-01  3.6647E-01 -2.6143E-01  6.5978E-02  8.7082E-01 -1.9608E+00 -3.5151E-01 -4.4975E-02
             1.4827E-01
 GRADIENT:   1.0706E+00  1.4610E-02  4.6109E-01 -2.7788E+00 -1.7618E-01  7.8609E-02  4.1699E-01  6.4714E-04  3.0047E-02 -5.3485E-02
            -5.0267E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1641.69663998926        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1975             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5060E-01  5.3150E-01  7.8306E-01  1.3049E+00  6.9666E-01  9.6656E-01  2.1595E+00  1.3415E-01  6.3688E-01  8.6531E-01
             1.0495E+00
 PARAMETER:  4.9339E-02 -5.3206E-01 -1.4454E-01  3.6616E-01 -2.6146E-01  6.5991E-02  8.6986E-01 -1.9088E+00 -3.5117E-01 -4.4666E-02
             1.4836E-01
 GRADIENT:   3.2156E+02  4.2833E+01  5.5253E+00  3.3168E+02  2.5974E+01  3.2583E+01  3.8173E+01  1.2251E-01  1.2676E+01  8.9964E-01
             1.3099E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1641.69694567329        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2157            RESET HESSIAN, TYPE II
 NPARAMETR:  9.5061E-01  5.3192E-01  7.8263E-01  1.3047E+00  6.9667E-01  9.6657E-01  2.1579E+00  1.2951E-01  6.3706E-01  8.6518E-01
             1.0496E+00
 PARAMETER:  4.9348E-02 -5.3127E-01 -1.4510E-01  3.6596E-01 -2.6144E-01  6.5999E-02  8.6915E-01 -1.9440E+00 -3.5090E-01 -4.4815E-02
             1.4842E-01
 GRADIENT:   3.2149E+02  4.2773E+01  5.2723E+00  3.3147E+02  2.6389E+01  3.2577E+01  3.8123E+01  1.1592E-01  1.2680E+01  8.6685E-01
             1.3247E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1641.69738267662        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2340             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5062E-01  5.3229E-01  7.8249E-01  1.3044E+00  6.9656E-01  9.6658E-01  2.1569E+00  1.2623E-01  6.3711E-01  8.6499E-01
             1.0496E+00
 PARAMETER:  4.9360E-02 -5.3057E-01 -1.4528E-01  3.6577E-01 -2.6161E-01  6.6011E-02  8.6866E-01 -1.9697E+00 -3.5081E-01 -4.5036E-02
             1.4837E-01
 GRADIENT:   3.2155E+02  4.2825E+01  5.6773E+00  3.3103E+02  2.5983E+01  3.2584E+01  3.8094E+01  1.0973E-01  1.2649E+01  8.0856E-01
             1.2613E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1641.69766814432        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2524             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5063E-01  5.3265E-01  7.8229E-01  1.3042E+00  6.9649E-01  9.6659E-01  2.1558E+00  1.2373E-01  6.3715E-01  8.6480E-01
             1.0495E+00
 PARAMETER:  4.9368E-02 -5.2990E-01 -1.4553E-01  3.6561E-01 -2.6171E-01  6.6018E-02  8.6818E-01 -1.9897E+00 -3.5075E-01 -4.5254E-02
             1.4833E-01
 GRADIENT:   3.2158E+02  4.2852E+01  5.8436E+00  3.3077E+02  2.5847E+01  3.2588E+01  3.8063E+01  1.0555E-01  1.2625E+01  7.6762E-01
             1.2225E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1641.69782038103        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2710
 NPARAMETR:  9.5064E-01  5.3306E-01  7.8190E-01  1.3039E+00  6.9647E-01  9.6660E-01  2.1544E+00  1.2642E-01  6.3729E-01  8.6480E-01
             1.0496E+00
 PARAMETER:  4.9377E-02 -5.2913E-01 -1.4603E-01  3.6540E-01 -2.6173E-01  6.6027E-02  8.6752E-01 -1.9682E+00 -3.5053E-01 -4.5261E-02
             1.4840E-01
 GRADIENT:   1.0553E+00 -9.7761E-02 -2.8378E-01 -2.6996E+00  4.5645E-01  7.4834E-02  4.0612E-01  6.1903E-03  1.1690E-01  1.6761E-01
             1.2897E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1641.69810141224        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2892
 NPARAMETR:  9.5065E-01  5.3339E-01  7.8172E-01  1.3037E+00  6.9639E-01  9.6661E-01  2.1534E+00  1.2481E-01  6.3734E-01  8.6468E-01
             1.0495E+00
 PARAMETER:  4.9387E-02 -5.2850E-01 -1.4626E-01  3.6524E-01 -2.6184E-01  6.6036E-02  8.6706E-01 -1.9809E+00 -3.5045E-01 -4.5394E-02
             1.4836E-01
 GRADIENT:   1.0628E+00 -5.9436E-02 -1.1133E-01 -2.7725E+00  2.4901E-01  7.6205E-02  4.0552E-01  5.4841E-03  1.0416E-01  1.5172E-01
             9.6467E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1641.69836375279        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3076
 NPARAMETR:  9.5065E-01  5.3377E-01  7.8136E-01  1.3035E+00  6.9637E-01  9.6661E-01  2.1521E+00  1.2141E-01  6.3743E-01  8.6450E-01
             1.0496E+00
 PARAMETER:  4.9394E-02 -5.2780E-01 -1.4672E-01  3.6507E-01 -2.6187E-01  6.6043E-02  8.6647E-01 -2.0086E+00 -3.5030E-01 -4.5600E-02
             1.4838E-01
 GRADIENT:   1.0546E+00 -7.9942E-02 -2.7635E-01 -2.6066E+00  4.9839E-01  7.4702E-02  3.9588E-01  4.9462E-03  1.0241E-01  1.2037E-01
             9.3394E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1641.69868155226        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     3257
 NPARAMETR:  9.5066E-01  5.3407E-01  7.8121E-01  1.3033E+00  6.9628E-01  9.6662E-01  2.1513E+00  1.1609E-01  6.3750E-01  8.6439E-01
             1.0495E+00
 PARAMETER:  4.9404E-02 -5.2722E-01 -1.4691E-01  3.6491E-01 -2.6201E-01  6.6052E-02  8.6606E-01 -2.0534E+00 -3.5020E-01 -4.5729E-02
             1.4836E-01
 GRADIENT:   1.0669E+00 -3.3028E-02  7.0337E-02 -2.7957E+00  1.6039E-01  7.7675E-02  3.9829E-01  2.8273E-03  8.0182E-02  5.3086E-02
             3.6244E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1641.69898481399        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3445             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5068E-01  5.3464E-01  7.8094E-01  1.3030E+00  6.9615E-01  9.6664E-01  2.1497E+00  1.1082E-01  6.3752E-01  8.6397E-01
             1.0495E+00
 PARAMETER:  4.9417E-02 -5.2616E-01 -1.4726E-01  3.6466E-01 -2.6219E-01  6.6066E-02  8.6533E-01 -2.0999E+00 -3.5017E-01 -4.6218E-02
             1.4827E-01
 GRADIENT:   3.2162E+02  4.2889E+01  6.2789E+00  3.2921E+02  2.5603E+01  3.2596E+01  3.7884E+01  8.7296E-02  1.2542E+01  6.4727E-01
             1.1237E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1641.69917156903        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3632
 NPARAMETR:  9.5068E-01  5.3493E-01  7.8073E-01  1.3028E+00  6.9611E-01  9.6664E-01  2.1488E+00  1.0974E-01  6.3759E-01  8.6386E-01
             1.0495E+00
 PARAMETER:  4.9424E-02 -5.2562E-01 -1.4753E-01  3.6452E-01 -2.6225E-01  6.6073E-02  8.6491E-01 -2.1096E+00 -3.5007E-01 -4.6345E-02
             1.4827E-01
 GRADIENT:   1.0762E+00  5.4396E-02  4.5491E-01 -2.8762E+00 -2.1126E-01  7.9480E-02  3.9276E-01  7.5397E-04  2.7672E-02 -5.6289E-02
            -6.0750E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1641.69927758940        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3817
 NPARAMETR:  9.5069E-01  5.3530E-01  7.8053E-01  1.3025E+00  6.9598E-01  9.6665E-01  2.1477E+00  1.0511E-01  6.3769E-01  8.6376E-01
             1.0494E+00
 PARAMETER:  4.9436E-02 -5.2493E-01 -1.4779E-01  3.6432E-01 -2.6243E-01  6.6085E-02  8.6442E-01 -2.1528E+00 -3.4990E-01 -4.6466E-02
             1.4825E-01
 GRADIENT:   1.0901E+00  1.0070E-01  7.8057E-01 -3.0912E+00 -5.6563E-01  8.2508E-02  3.9854E-01 -2.6957E-04  1.7988E-02 -9.5406E-02
            -9.9682E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1641.69963959665        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     4002             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5070E-01  5.3558E-01  7.8008E-01  1.3024E+00  6.9601E-01  9.6666E-01  2.1466E+00  1.0988E-01  6.3782E-01  8.6388E-01
             1.0495E+00
 PARAMETER:  4.9441E-02 -5.2440E-01 -1.4835E-01  3.6418E-01 -2.6239E-01  6.6088E-02  8.6387E-01 -2.1083E+00 -3.4971E-01 -4.6320E-02
             1.4835E-01
 GRADIENT:   3.2154E+02  4.2808E+01  5.9319E+00  3.2842E+02  2.5980E+01  3.2589E+01  3.7793E+01  8.8395E-02  1.2563E+01  7.6307E-01
             1.2215E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1641.69979765833        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     4189
 NPARAMETR:  9.5070E-01  5.3584E-01  7.7989E-01  1.3022E+00  6.9597E-01  9.6666E-01  2.1458E+00  1.0885E-01  6.3788E-01  8.6375E-01
             1.0495E+00
 PARAMETER:  4.9447E-02 -5.2391E-01 -1.4861E-01  3.6405E-01 -2.6245E-01  6.6094E-02  8.6349E-01 -2.1178E+00 -3.4961E-01 -4.6474E-02
             1.4835E-01
 GRADIENT:   1.0672E+00 -4.3750E-03  7.1330E-02 -2.8170E+00  1.2096E-01  7.7365E-02  3.8900E-01  2.7682E-03  7.2704E-02  5.3326E-02
             3.2139E-02

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1641.70002250781        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     4377             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5072E-01  5.3631E-01  7.7964E-01  1.3019E+00  6.9585E-01  9.6667E-01  2.1444E+00  1.0228E-01  6.3791E-01  8.6338E-01
             1.0495E+00
 PARAMETER:  4.9459E-02 -5.2304E-01 -1.4892E-01  3.6384E-01 -2.6263E-01  6.6106E-02  8.6288E-01 -2.1800E+00 -3.4956E-01 -4.6902E-02
             1.4827E-01
 GRADIENT:   3.2160E+02  4.2863E+01  6.3357E+00  3.2785E+02  2.5680E+01  3.2597E+01  3.7731E+01  7.7404E-02  1.2509E+01  6.4247E-01
             1.1223E+00

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1641.70017288332        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     4564
 NPARAMETR:  9.5072E-01  5.3655E-01  7.7944E-01  1.3018E+00  6.9580E-01  9.6668E-01  2.1437E+00  1.0120E-01  6.3798E-01  8.6332E-01
             1.0495E+00
 PARAMETER:  4.9465E-02 -5.2260E-01 -1.4918E-01  3.6372E-01 -2.6269E-01  6.6112E-02  8.6252E-01 -2.1906E+00 -3.4945E-01 -4.6966E-02
             1.4828E-01
 GRADIENT:   1.0761E+00  6.5025E-02  4.1651E-01 -2.9020E+00 -1.9339E-01  7.9534E-02  3.8408E-01  8.9183E-04  3.1674E-02 -5.0861E-02
            -5.4440E-02

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1641.70033283911        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     4752             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5073E-01  5.3685E-01  7.7899E-01  1.3015E+00  6.9579E-01  9.6669E-01  2.1425E+00  1.0557E-01  6.3817E-01  8.6350E-01
             1.0496E+00
 PARAMETER:  4.9473E-02 -5.2203E-01 -1.4976E-01  3.6354E-01 -2.6271E-01  6.6118E-02  8.6196E-01 -2.1483E+00 -3.4916E-01 -4.6758E-02
             1.4840E-01
 GRADIENT:   3.2148E+02  4.2744E+01  5.7415E+00  3.2737E+02  2.6264E+01  3.2585E+01  3.7676E+01  8.4363E-02  1.2562E+01  8.1953E-01
             1.2768E+00

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1641.70050270356        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     4937             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5073E-01  5.3707E-01  7.7885E-01  1.3014E+00  6.9571E-01  9.6669E-01  2.1418E+00  1.0301E-01  6.3821E-01  8.6339E-01
             1.0496E+00
 PARAMETER:  4.9480E-02 -5.2163E-01 -1.4994E-01  3.6343E-01 -2.6283E-01  6.6125E-02  8.6166E-01 -2.1730E+00 -3.4909E-01 -4.6886E-02
             1.4837E-01
 GRADIENT:   3.2151E+02  4.2769E+01  5.9130E+00  3.2718E+02  2.6107E+01  3.2588E+01  3.7657E+01  8.0665E-02  1.2547E+01  7.9171E-01
             1.2445E+00

0ITERATION NO.:  138    OBJECTIVE VALUE:  -1641.70050776134        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     5031
 NPARAMETR:  9.5069E-01  5.3707E-01  7.7885E-01  1.3015E+00  6.9570E-01  9.6667E-01  2.1413E+00  9.8743E-02  6.3817E-01  8.6334E-01
             1.0495E+00
 PARAMETER:  4.9438E-02 -5.2163E-01 -1.4994E-01  3.6352E-01 -2.6284E-01  6.6107E-02  8.6142E-01 -2.2152E+00 -3.4914E-01 -4.6941E-02
             1.4836E-01
 GRADIENT:  -1.1285E-01 -6.6977E-02  5.7725E-02  5.8598E-01  2.3966E-01 -9.5239E-03 -1.3783E-02  1.7693E-04 -1.3651E-02  1.4068E-02
             9.7413E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     5031
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1762E-04  2.7796E-02 -4.9240E-03 -2.9381E-02  3.7225E-03
 SE:             2.9824E-02  2.1245E-02  2.2602E-03  2.3016E-02  2.3419E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8080E-01  1.9076E-01  2.9363E-02  2.0176E-01  8.7371E-01

 ETASHRINKSD(%)  8.7406E-02  2.8826E+01  9.2428E+01  2.2892E+01  2.1542E+01
 ETASHRINKVR(%)  1.7474E-01  4.9343E+01  9.9427E+01  4.0544E+01  3.8444E+01
 EBVSHRINKSD(%)  5.0828E-01  2.9901E+01  9.2948E+01  2.1501E+01  1.9217E+01
 EBVSHRINKVR(%)  1.0140E+00  5.0861E+01  9.9503E+01  3.8380E+01  3.4741E+01
 RELATIVEINF(%)  9.8353E+01  8.1105E+00  4.3508E-02  1.1075E+01  5.1945E+00
 EPSSHRINKSD(%)  4.2823E+01
 EPSSHRINKVR(%)  6.7308E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1641.7005077613417     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -906.54968119760349     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    68.99
 Elapsed covariance  time in seconds:     5.77
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1641.701       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.51E-01  5.37E-01  7.79E-01  1.30E+00  6.96E-01  9.67E-01  2.14E+00  9.87E-02  6.38E-01  8.63E-01  1.05E+00
 


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
 
         2.81E-02  1.59E-01  2.62E-01  7.19E-02  1.83E-01  6.35E-02  4.88E-01  2.46E+00  9.14E-02  2.11E-01  7.38E-02
 


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
+        7.91E-04
 
 TH 2
+       -1.53E-05  2.54E-02
 
 TH 3
+       -1.24E-03  1.47E-02  6.85E-02
 
 TH 4
+       -1.72E-04 -9.49E-03  1.07E-03  5.17E-03
 
 TH 5
+       -7.73E-04  1.60E-02  4.64E-02 -2.05E-03  3.36E-02
 
 TH 6
+        1.80E-04 -2.02E-03 -1.19E-03  7.56E-04 -1.21E-03  4.03E-03
 
 TH 7
+        1.72E-04 -6.95E-02 -3.52E-02  2.67E-02 -4.13E-02  9.55E-03  2.39E-01
 
 TH 8
+       -4.81E-03  1.01E-01  5.39E-01  4.65E-03  3.57E-01  2.51E-03 -2.08E-01  6.06E+00
 
 TH 9
+        2.37E-04  6.51E-03 -5.84E-03 -3.19E-03 -1.79E-03 -1.10E-03 -2.39E-02 -8.05E-02  8.36E-03
 
 TH10
+       -1.15E-03  1.09E-02  3.98E-02  5.44E-04  2.79E-02 -1.56E-03 -2.31E-02  1.85E-01 -2.25E-03  4.47E-02
 
 TH11
+        3.77E-04 -8.03E-04 -3.22E-03  4.13E-04 -2.03E-03 -7.64E-04 -9.24E-04 -4.15E-02  1.63E-03 -1.64E-03  5.45E-03
 
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
+        2.81E-02
 
 TH 2
+       -3.41E-03  1.59E-01
 
 TH 3
+       -1.68E-01  3.52E-01  2.62E-01
 
 TH 4
+       -8.49E-02 -8.29E-01  5.71E-02  7.19E-02
 
 TH 5
+       -1.50E-01  5.49E-01  9.67E-01 -1.56E-01  1.83E-01
 
 TH 6
+        1.01E-01 -1.99E-01 -7.15E-02  1.66E-01 -1.04E-01  6.35E-02
 
 TH 7
+        1.25E-02 -8.92E-01 -2.75E-01  7.61E-01 -4.61E-01  3.08E-01  4.88E-01
 
 TH 8
+       -6.95E-02  2.58E-01  8.37E-01  2.63E-02  7.91E-01  1.60E-02 -1.73E-01  2.46E+00
 
 TH 9
+        9.20E-02  4.47E-01 -2.44E-01 -4.86E-01 -1.07E-01 -1.90E-01 -5.36E-01 -3.58E-01  9.14E-02
 
 TH10
+       -1.93E-01  3.23E-01  7.19E-01  3.58E-02  7.20E-01 -1.16E-01 -2.24E-01  3.57E-01 -1.16E-01  2.11E-01
 
 TH11
+        1.82E-01 -6.82E-02 -1.67E-01  7.79E-02 -1.50E-01 -1.63E-01 -2.56E-02 -2.28E-01  2.42E-01 -1.05E-01  7.38E-02
 
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
+        1.41E+03
 
 TH 2
+       -1.08E+01  5.20E+02
 
 TH 3
+        1.58E+01  2.08E+02  1.11E+03
 
 TH 4
+        8.91E+01  4.54E+02 -4.09E+02  1.31E+03
 
 TH 5
+        4.27E+01 -3.77E+02 -1.42E+03  3.89E+02  2.17E+03
 
 TH 6
+       -7.45E+01 -3.59E+01  7.12E+01 -5.21E+01 -1.08E+02  3.02E+02
 
 TH 7
+       -8.01E+00  5.54E+01 -3.53E+00 -8.34E+00  2.38E+01 -2.40E+01  2.68E+01
 
 TH 8
+       -4.14E+00 -2.63E+00 -1.61E+01  4.51E+00  5.59E+00 -4.54E-01 -4.87E-01  1.17E+00
 
 TH 9
+       -2.69E+01 -5.22E+01 -3.63E+01 -3.34E+01  3.72E+01 -1.66E+00  2.20E+01  5.11E+00  2.57E+02
 
 TH10
+        1.81E+00 -4.58E+01 -7.93E+01 -3.36E+01 -2.28E+01  1.32E+01 -9.03E+00  6.52E+00  2.31E+01  8.91E+01
 
 TH11
+       -1.15E+02  1.15E+01  6.17E+01 -1.10E+02 -1.03E+02  4.49E+01  4.30E+00  1.28E+00 -3.85E+01  1.02E+01  2.31E+02
 
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
 #CPUT: Total CPU Time in Seconds,       74.854
Stop Time:
Wed Sep 29 17:37:11 CDT 2021
