Sat Sep 25 14:12:19 CDT 2021
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
$DATA ../../../../data/spa/D/dat28.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   9795.25259274023        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5677E+01 -1.3914E+02 -3.0207E+01 -1.9777E+02  1.6350E+02 -1.0491E+03 -6.1607E+02 -4.7546E+01 -7.9385E+02 -2.0215E+02
            -2.0165E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -669.942417794743        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4881E+00  1.3711E+00  9.7008E-01  1.8770E+00  1.1055E+00  2.1277E+00  1.3450E+00  9.7879E-01  1.2942E+00  1.0227E+00
             1.3863E+01
 PARAMETER:  4.9753E-01  4.1563E-01  6.9621E-02  7.2965E-01  2.0028E-01  8.5503E-01  3.9642E-01  7.8564E-02  3.5791E-01  1.2247E-01
             2.7292E+00
 GRADIENT:   1.9846E+01  3.6253E+01  1.3328E+00  5.4276E+01 -1.4860E+01  4.8623E+01 -2.8553E+00  4.2983E+00  3.6073E+00  3.4374E+00
             1.5367E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -690.447286371399        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.4454E+00  8.3015E-01  1.5612E+00  2.4422E+00  1.4587E+00  1.9515E+00  3.9911E+00  3.8628E-01  1.5197E+00  1.3940E+00
             1.2409E+01
 PARAMETER:  4.6837E-01 -8.6151E-02  5.4546E-01  9.9290E-01  4.7754E-01  7.6858E-01  1.4841E+00 -8.5120E-01  5.1854E-01  4.3214E-01
             2.6184E+00
 GRADIENT:   1.9301E+01  2.5551E+01  3.7456E+00  8.1288E+01 -1.0487E+01  7.1826E+00  3.7849E+00  5.5536E-01  8.3550E+00  1.2584E+00
             1.3724E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -746.387615524907        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.2060E+00  5.2654E-01  1.1407E+00  1.7277E+00  2.2000E+00  1.7308E+00  2.5195E+00  2.3488E-02  1.0720E+00  1.2696E+00
             9.1454E+00
 PARAMETER:  2.8727E-01 -5.4142E-01  2.3162E-01  6.4681E-01  8.8844E-01  6.4861E-01  1.0241E+00 -3.6513E+00  1.6955E-01  3.3869E-01
             2.3132E+00
 GRADIENT:   2.6557E+01  7.7551E+00 -1.0193E+00 -1.6767E+01 -1.1197E+01  4.5333E+00  1.8369E+00  1.5704E-03 -8.6869E-02  1.1364E+00
             6.7972E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -773.858308512699        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0377E+00  3.2160E-01  3.7325E-01  1.5082E+00  6.5627E+00  1.5542E+00  1.1381E+00  1.0000E-02  6.9623E-01  4.3918E+00
             9.0820E+00
 PARAMETER:  1.3701E-01 -1.0344E+00 -8.8551E-01  5.1092E-01  1.9814E+00  5.4099E-01  2.2937E-01 -1.0500E+01 -2.6207E-01  1.5797E+00
             2.3063E+00
 GRADIENT:   1.0301E+01  5.2823E+01 -5.4406E+01  7.2873E+01 -3.4645E+01 -2.1072E+01  5.2777E+00  0.0000E+00 -4.5268E+00  1.8998E+01
            -1.8518E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -863.563288919988        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  5.5312E-01  2.2987E-02  4.0544E-02  4.5616E-01  7.2818E+00  1.4144E+00  7.5937E-02  1.0000E-02  1.7653E-01  3.4178E-01
             8.7842E+00
 PARAMETER: -4.9217E-01 -3.6728E+00 -3.1054E+00 -6.8492E-01  2.0854E+00  4.4673E-01 -2.4779E+00 -2.5155E+01 -1.6343E+00 -9.7359E-01
             2.2730E+00
 GRADIENT:   2.7589E+01 -7.4050E-02 -3.0426E+01  4.4219E+01 -1.1158E+00  3.4092E+00  6.2853E-03  0.0000E+00  6.4866E-01  1.4160E-01
             1.0408E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -865.182681331819        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      468
 NPARAMETR:  4.3706E-01  1.1545E-02  2.3763E-02  3.0400E-01  8.2819E+00  1.3140E+00  3.3904E-02  1.0000E-02  1.1409E-01  2.3348E-01
             8.4864E+00
 PARAMETER: -7.2768E-01 -4.3615E+00 -3.6396E+00 -1.0907E+00  2.2141E+00  3.7307E-01 -3.2842E+00 -2.8985E+01 -2.0708E+00 -1.3547E+00
             2.2385E+00
 GRADIENT:   3.7891E+00  2.7553E-01 -4.1987E+01  5.1398E+01  5.8045E-01 -7.4286E+00  6.3116E-05  0.0000E+00  1.0706E-01  8.6294E-03
            -1.8003E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -865.984886502376        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      644
 NPARAMETR:  4.5399E-01  1.2519E-02  2.6424E-02  3.2312E-01  9.0093E+00  1.3499E+00  3.7990E-02  1.0000E-02  1.1841E-01  2.9773E-01
             8.6238E+00
 PARAMETER: -6.8968E-01 -4.2805E+00 -3.5335E+00 -1.0297E+00  2.2983E+00  4.0000E-01 -3.1704E+00 -2.8620E+01 -2.0336E+00 -1.1116E+00
             2.2545E+00
 GRADIENT:  -6.0822E-01  5.9044E-01  4.2580E-01 -1.5656E+00 -3.5099E-01 -1.7929E-01  5.3255E-05  0.0000E+00  1.2078E-01  3.6280E-03
             1.5065E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -866.090167791057        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      820
 NPARAMETR:  4.5584E-01  1.0000E-02  2.7063E-02  3.2845E-01  1.4753E+01  1.3543E+00  3.8518E-02  1.0000E-02  9.1248E-02  7.6964E-01
             8.6194E+00
 PARAMETER: -6.8560E-01 -4.5893E+00 -3.5096E+00 -1.0134E+00  2.7914E+00  4.0330E-01 -3.1566E+00 -2.9962E+01 -2.2942E+00 -1.6184E-01
             2.2540E+00
 GRADIENT:  -7.4275E-01  0.0000E+00 -5.6834E-02  3.9867E-01 -1.4493E-03  2.5607E-01  6.4280E-07  0.0000E+00  5.9222E-02 -3.4951E-05
            -2.4984E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -866.116852244251        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      998
 NPARAMETR:  4.5405E-01  1.0000E-02  2.6969E-02  3.2735E-01  1.1590E+02  1.3510E+00  3.4696E-02  1.0000E-02  2.6978E-02  3.7463E+01
             8.6381E+00
 PARAMETER: -6.8954E-01 -6.0051E+00 -3.5131E+00 -1.0167E+00  4.8527E+00  4.0084E-01 -3.2611E+00 -3.6281E+01 -3.5127E+00  3.7234E+00
             2.2562E+00
 GRADIENT:  -2.0910E-01  0.0000E+00  8.4888E-02 -2.6186E-02  1.2073E-03 -8.2297E-02  2.3981E-07  0.0000E+00  6.0457E-03  2.0506E-04
             1.3322E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -866.119492544111        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1157
 NPARAMETR:  4.5407E-01  1.0000E-02  2.6982E-02  3.2752E-01  1.0489E+02  1.3516E+00  3.4707E-02  1.0000E-02  1.0000E-02  3.5646E+01
             8.6352E+00
 PARAMETER: -6.8949E-01 -6.0195E+00 -3.5126E+00 -1.0162E+00  4.7529E+00  4.0131E-01 -3.2608E+00 -3.6344E+01 -5.5474E+00  3.6737E+00
             2.2558E+00
 GRADIENT:  -2.7116E-01  0.0000E+00 -2.1178E-01  3.9879E-01  1.3907E-03  5.5725E-02  2.4667E-07  0.0000E+00  0.0000E+00  2.2510E-04
            -3.2184E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -866.124091067078        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1340
 NPARAMETR:  4.5420E-01  1.0000E-02  2.6862E-02  3.2664E-01  1.9709E+01  1.3512E+00  3.4707E-02  1.0000E-02  1.0000E-02  1.1955E+01
             8.6377E+00
 PARAMETER: -6.8921E-01 -6.0195E+00 -3.5171E+00 -1.0189E+00  3.0811E+00  4.0100E-01 -3.2608E+00 -3.6344E+01 -1.2329E+01  2.5812E+00
             2.2561E+00
 GRADIENT:  -2.2076E-01  0.0000E+00 -3.0665E-02  6.9750E-03  4.8783E-03 -4.1093E-02  4.5035E-07  0.0000E+00  0.0000E+00 -3.7956E-03
             6.3380E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -866.124373692655        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1515
 NPARAMETR:  4.5441E-01  1.0000E-02  2.6873E-02  3.2675E-01  1.8761E+01  1.3513E+00  3.4707E-02  1.0000E-02  1.0000E-02  1.2127E+01
             8.6386E+00
 PARAMETER: -6.8876E-01 -6.0195E+00 -3.5166E+00 -1.0186E+00  3.0318E+00  4.0105E-01 -3.2608E+00 -3.6344E+01 -1.1291E+01  2.5955E+00
             2.2562E+00
 GRADIENT:  -3.3512E-02  0.0000E+00 -1.0231E-02 -1.1793E-01 -1.2677E-03 -3.8024E-02  5.0324E-07  0.0000E+00  0.0000E+00  1.1784E-03
             6.6674E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -866.124507634987        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1697
 NPARAMETR:  4.5508E-01  1.0000E-02  2.6978E-02  3.2770E-01  1.8577E+01  1.3509E+00  3.4705E-02  1.0000E-02  1.0000E-02  1.2018E+01
             8.6403E+00
 PARAMETER: -6.8729E-01 -6.0195E+00 -3.5127E+00 -1.0156E+00  3.0219E+00  4.0077E-01 -3.2609E+00 -3.6344E+01 -4.5179E+00  2.5864E+00
             2.2564E+00
 GRADIENT:  -8.0622E-02  0.0000E+00  5.5557E-02 -1.0432E-01 -2.0400E-03 -2.0900E-01  4.8583E-07  0.0000E+00  0.0000E+00  1.6341E-03
             2.3687E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -866.124620606349        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1867
 NPARAMETR:  4.5510E-01  1.0000E-02  2.6978E-02  3.2772E-01  1.8215E+01  1.3518E+00  3.4705E-02  1.0000E-02  1.0000E-02  1.1664E+01
             8.6380E+00
 PARAMETER: -6.8723E-01 -6.0195E+00 -3.5128E+00 -1.0156E+00  3.0023E+00  4.0146E-01 -3.2609E+00 -3.6344E+01 -4.5093E+00  2.5565E+00
             2.2562E+00
 GRADIENT:  -2.5041E-02  0.0000E+00 -3.9579E-02 -2.2825E-02 -1.2389E-04 -1.4549E-03  4.9545E-07  0.0000E+00  3.6177E-05  9.3029E-05
            -1.2030E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1867
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9150E-03  6.4341E-06  1.0321E-04 -2.1820E-04 -3.9710E-04
 SE:             2.8892E-02  4.4280E-06  3.3070E-04  4.1914E-04  7.9515E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4715E-01  1.4621E-01  7.5496E-01  6.0265E-01  6.1750E-01

 ETASHRINKSD(%)  3.2079E+00  9.9985E+01  9.8892E+01  9.8596E+01  9.7336E+01
 ETASHRINKVR(%)  6.3129E+00  1.0000E+02  9.9988E+01  9.9980E+01  9.9929E+01
 EBVSHRINKSD(%)  3.5414E+00  9.9984E+01  9.8896E+01  9.8604E+01  9.7506E+01
 EBVSHRINKVR(%)  6.9574E+00  1.0000E+02  9.9988E+01  9.9981E+01  9.9938E+01
 RELATIVEINF(%)  4.2145E-01  9.9842E-08  2.6887E-05  2.8173E-05  5.6076E-04
 EPSSHRINKSD(%)  8.3833E+00
 EPSSHRINKVR(%)  1.6064E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -866.12462060634925     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -130.97379404261108     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -866.125       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.55E-01  1.00E-02  2.70E-02  3.28E-01  1.82E+01  1.35E+00  3.47E-02  1.00E-02  1.00E-02  1.17E+01  8.64E+00
 


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
+        2.68E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.33E+04  0.00E+00  1.18E+06
 
 TH 4
+       -3.89E+02  0.00E+00 -1.20E+05  1.36E+04
 
 TH 5
+        1.55E-01  0.00E+00 -2.73E+00  2.19E-01  7.05E-04
 
 TH 6
+       -1.52E+01  0.00E+00  3.38E+02 -5.86E+01  1.76E-03  8.49E+01
 
 TH 7
+        6.07E-02  0.00E+00  1.20E-01  7.95E-04  1.38E-04 -5.91E-02  6.65E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.91E-01  0.00E+00 -1.09E+00 -1.45E-01  1.08E-03 -2.84E-03 -1.64E+00  0.00E+00  2.58E+02
 
 TH10
+        1.93E-03  0.00E+00 -1.28E-01 -1.31E-02 -7.90E-04  9.62E-05  1.99E-03  0.00E+00 -2.70E-03  9.56E-04
 
 TH11
+       -2.62E+01  0.00E+00  4.27E+02 -3.01E+01 -1.90E-03  3.51E-01 -2.01E-03  0.00E+00  1.84E-02 -2.14E-04  6.46E+00
 
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
 #CPUT: Total CPU Time in Seconds,       30.733
Stop Time:
Sat Sep 25 14:12:55 CDT 2021
