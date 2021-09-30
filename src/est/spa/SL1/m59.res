Wed Sep 29 15:10:25 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat59.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m59.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1696.86354877397        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4623E+02 -7.9876E+00 -2.7375E+01  2.7518E+01  8.0481E+01  3.5406E+01  5.5527E+00  5.7432E+00 -1.5296E+01  4.7237E-01
             1.6358E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1708.71780047002        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0659E+00  9.0619E-01  9.7003E-01  1.1233E+00  8.8668E-01  1.0673E+00  9.1678E-01  9.2866E-01  1.1405E+00  9.0746E-01
             9.4159E-01
 PARAMETER:  1.6386E-01  1.4959E-03  6.9574E-02  2.1623E-01 -2.0267E-02  1.6515E-01  1.3111E-02  2.5992E-02  2.3146E-01  2.8909E-03
             3.9818E-02
 GRADIENT:   6.1388E-01  1.4600E+01 -3.6651E+00  3.7819E+01  7.5098E+00 -1.0032E+01  2.6416E+00  4.4143E+00  1.9253E+01 -4.9070E-01
            -6.6459E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1710.48034001910        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0747E+00  7.6528E-01  8.5111E-01  1.2032E+00  7.7952E-01  1.0753E+00  9.1489E-01  4.9363E-01  1.0439E+00  9.2721E-01
             9.3107E-01
 PARAMETER:  1.7205E-01 -1.6752E-01 -6.1211E-02  2.8498E-01 -1.4907E-01  1.7261E-01  1.1045E-02 -6.0596E-01  1.4293E-01  2.4426E-02
             2.8576E-02
 GRADIENT:   1.5963E+01  1.5678E+01 -1.0360E+01  5.0626E+01  1.0067E+01 -7.1613E+00 -2.0368E+00 -5.9013E-01  1.0857E+01  4.0850E+00
            -1.1443E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1712.60294519025        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      553
 NPARAMETR:  1.0646E+00  6.5353E-01  8.1955E-01  1.2357E+00  7.1521E-01  1.0968E+00  1.1846E+00  4.5899E-01  9.3408E-01  8.6446E-01
             9.5564E-01
 PARAMETER:  1.6258E-01 -3.2536E-01 -9.9005E-02  3.1167E-01 -2.3517E-01  1.9243E-01  2.6937E-01 -6.7874E-01  3.1811E-02 -4.5645E-02
             5.4625E-02
 GRADIENT:  -2.0014E+00  7.9954E+00  4.7356E+00  3.6781E+00 -1.0449E+01  1.4685E+00 -3.3917E-01 -2.6577E-01 -2.1574E+00  1.7253E+00
             1.3386E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1713.38704374372        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  1.0604E+00  3.7683E-01  8.8655E-01  1.4057E+00  6.7246E-01  1.0840E+00  1.6157E+00  5.5207E-01  8.6013E-01  8.7249E-01
             9.4654E-01
 PARAMETER:  1.5861E-01 -8.7595E-01 -2.0419E-02  4.4051E-01 -2.9681E-01  1.8067E-01  5.7978E-01 -4.9408E-01 -5.0666E-02 -3.6400E-02
             4.5062E-02
 GRADIENT:   7.6474E-02  6.3868E+00  8.1190E+00  1.2726E+01 -1.3411E+01 -1.2365E+00  1.8869E-01 -5.2563E-01  5.0682E-01 -6.5819E-02
            -1.4418E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1713.44840639461        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      905
 NPARAMETR:  1.0575E+00  2.7010E-01  9.3551E-01  1.4731E+00  6.7200E-01  1.0849E+00  1.8545E+00  6.5072E-01  8.2721E-01  8.8044E-01
             9.4842E-01
 PARAMETER:  1.5590E-01 -1.2089E+00  3.3342E-02  4.8739E-01 -2.9750E-01  1.8150E-01  7.1760E-01 -3.2967E-01 -8.9691E-02 -2.7332E-02
             4.7046E-02
 GRADIENT:   1.8460E-01  5.0310E+00  6.2660E+00  1.5416E+01 -1.1482E+01  7.7892E-02 -4.5452E-04  2.3639E-01 -2.0858E+00  6.3095E-01
            -3.9006E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1713.49367145112        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1080
 NPARAMETR:  1.0558E+00  2.0966E-01  9.5573E-01  1.5077E+00  6.6943E-01  1.0840E+00  2.0404E+00  6.8558E-01  8.1352E-01  8.8447E-01
             9.4844E-01
 PARAMETER:  1.5431E-01 -1.4623E+00  5.4719E-02  5.1057E-01 -3.0133E-01  1.8065E-01  8.1316E-01 -2.7749E-01 -1.0638E-01 -2.2763E-02
             4.7068E-02
 GRADIENT:   1.2962E-01  3.4771E+00  4.2575E+00  1.1507E+01 -7.8642E+00  2.4532E-01 -6.4406E-02  2.9151E-01 -1.9042E+00  5.1908E-01
             1.8343E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1713.69021281755        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1264
 NPARAMETR:  1.0581E+00  1.9452E-01  9.5410E-01  1.5005E+00  6.6836E-01  1.0847E+00  2.1256E+00  6.7727E-01  8.1299E-01  8.8837E-01
             9.4761E-01
 PARAMETER:  1.5646E-01 -1.5372E+00  5.3015E-02  5.0580E-01 -3.0294E-01  1.8133E-01  8.5405E-01 -2.8969E-01 -1.0704E-01 -1.8368E-02
             4.6184E-02
 GRADIENT:   5.1919E+00  4.5812E-01  1.8245E+00 -2.0014E+01  3.5050E-01  6.7665E-01  3.1813E-02  1.9833E-01 -4.1076E-01  4.8503E-01
             1.0092E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1713.69707514814        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1446
 NPARAMETR:  1.0578E+00  1.9371E-01  9.5097E-01  1.5017E+00  6.6672E-01  1.0846E+00  2.1420E+00  6.7446E-01  8.1317E-01  8.8900E-01
             9.4767E-01
 PARAMETER:  1.5615E-01 -1.5414E+00  4.9732E-02  5.0658E-01 -3.0538E-01  1.8126E-01  8.6176E-01 -2.9384E-01 -1.0682E-01 -1.7654E-02
             4.6253E-02
 GRADIENT:   4.5754E+00  6.4281E-01  1.3360E+00 -1.7603E+01  9.3165E-02  6.4634E-01  5.0449E-02  3.0273E-01 -2.2812E-01  8.0715E-01
             2.2690E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1713.70508118878        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1627
 NPARAMETR:  1.0577E+00  1.9043E-01  9.4715E-01  1.5025E+00  6.6543E-01  1.0846E+00  2.1596E+00  6.7191E-01  8.1292E-01  8.8718E-01
             9.4761E-01
 PARAMETER:  1.5611E-01 -1.5584E+00  4.5697E-02  5.0712E-01 -3.0732E-01  1.8122E-01  8.6992E-01 -2.9763E-01 -1.0713E-01 -1.9705E-02
             4.6192E-02
 GRADIENT:   4.6164E+00  3.3431E-01 -1.1933E+00 -1.8619E+01  3.6046E+00  6.4915E-01  3.8307E-02  3.6198E-01 -7.8643E-02  5.1117E-01
             2.7548E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1713.71194087146        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1810             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0577E+00  1.8753E-01  9.4487E-01  1.5027E+00  6.6375E-01  1.0846E+00  2.1819E+00  6.6876E-01  8.1283E-01  8.8691E-01
             9.4762E-01
 PARAMETER:  1.5611E-01 -1.5738E+00  4.3294E-02  5.0723E-01 -3.0986E-01  1.8123E-01  8.8021E-01 -3.0233E-01 -1.0724E-01 -2.0013E-02
             4.6198E-02
 GRADIENT:   7.9906E+02  3.2917E+01  4.1156E+00  1.0789E+03  4.3358E+01  1.2395E+02  4.9610E+00  9.9167E-01  1.6968E+01  1.5406E+00
             1.1218E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1713.71574764673        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1992
 NPARAMETR:  1.0572E+00  1.8652E-01  9.4514E-01  1.5049E+00  6.6258E-01  1.0843E+00  2.1803E+00  6.6233E-01  8.1249E-01  8.8494E-01
             9.4725E-01
 PARAMETER:  1.5564E-01 -1.5792E+00  4.3574E-02  5.0873E-01 -3.1161E-01  1.8091E-01  8.7947E-01 -3.1200E-01 -1.0766E-01 -2.2237E-02
             4.5808E-02
 GRADIENT:   3.8485E+00  5.7016E-01  1.5079E+00 -1.7297E+01  5.2961E-01  5.4719E-01  1.3644E-02 -1.1605E-01  1.2654E-01 -1.6264E-01
            -1.3088E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1713.71949677943        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2177
 NPARAMETR:  1.0576E+00  1.8465E-01  9.4402E-01  1.5045E+00  6.6093E-01  1.0845E+00  2.2072E+00  6.5928E-01  8.1188E-01  8.8450E-01
             9.4711E-01
 PARAMETER:  1.5602E-01 -1.5893E+00  4.2396E-02  5.0846E-01 -3.1411E-01  1.8115E-01  8.9172E-01 -3.1661E-01 -1.0841E-01 -2.2737E-02
             4.5661E-02
 GRADIENT:   4.6602E+00  5.2689E-01  3.2670E+00 -1.9811E+01 -1.5714E+00  6.6096E-01  3.3739E-02 -2.2490E-01  1.1596E-01 -1.8167E-01
            -2.2073E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1713.72611822874        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2360             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0576E+00  1.8051E-01  9.4189E-01  1.5057E+00  6.5975E-01  1.0845E+00  2.2406E+00  6.6048E-01  8.1127E-01  8.8492E-01
             9.4725E-01
 PARAMETER:  1.5596E-01 -1.6120E+00  4.0131E-02  5.0926E-01 -3.1590E-01  1.8112E-01  9.0673E-01 -3.1479E-01 -1.0916E-01 -2.2256E-02
             4.5805E-02
 GRADIENT:   7.9808E+02  3.1801E+01  6.7613E+00  1.0894E+03  4.1224E+01  1.2374E+02  4.9820E+00  6.5565E-01  1.7421E+01  1.1204E+00
             8.2102E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1713.72850763894        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2543
 NPARAMETR:  1.0575E+00  1.8111E-01  9.4026E-01  1.5058E+00  6.5859E-01  1.0845E+00  2.2401E+00  6.5578E-01  8.1099E-01  8.8349E-01
             9.4703E-01
 PARAMETER:  1.5595E-01 -1.6086E+00  3.8396E-02  5.0930E-01 -3.1765E-01  1.8110E-01  9.0650E-01 -3.2194E-01 -1.0950E-01 -2.3878E-02
             4.5575E-02
 GRADIENT:   4.6557E+00  4.7640E-01  2.6632E+00 -2.0089E+01 -1.0869E+00  6.5892E-01  3.7604E-02 -1.9737E-01  8.6480E-02 -2.0735E-01
            -1.9664E-01

0ITERATION NO.:   72    OBJECTIVE VALUE:  -1713.72869445447        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     2602
 NPARAMETR:  1.0570E+00  1.8032E-01  9.3932E-01  1.5077E+00  6.5866E-01  1.0841E+00  2.2356E+00  6.5746E-01  8.1095E-01  8.8390E-01
             9.4712E-01
 PARAMETER:  1.5546E-01 -1.6130E+00  3.7396E-02  5.1056E-01 -3.1754E-01  1.8077E-01  9.0451E-01 -3.1937E-01 -1.0955E-01 -2.3417E-02
             4.5670E-02
 GRADIENT:  -9.1023E-01  1.3398E-01  4.3889E-01  3.5097E+00  9.4611E-01 -1.3014E-01 -9.5588E-03 -3.9114E-02  4.4098E-02 -4.9488E-02
            -6.6664E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2602
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.5778E-04 -2.0631E-03 -1.9100E-02 -5.0081E-03 -2.0276E-02
 SE:             2.9893E-02  6.5349E-03  1.4000E-02  2.8851E-02  2.3826E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8511E-01  7.5223E-01  1.7249E-01  8.6219E-01  3.9477E-01

 ETASHRINKSD(%)  1.0000E-10  7.8107E+01  5.3097E+01  3.3471E+00  2.0179E+01
 ETASHRINKVR(%)  1.0000E-10  9.5207E+01  7.8001E+01  6.5822E+00  3.6286E+01
 EBVSHRINKSD(%)  3.2458E-01  7.8637E+01  5.5478E+01  3.5276E+00  1.8259E+01
 EBVSHRINKVR(%)  6.4810E-01  9.5436E+01  8.0178E+01  6.9307E+00  3.3184E+01
 RELATIVEINF(%)  9.5114E+01  2.0987E-01  2.5555E+00  6.4514E+00  4.3776E+00
 EPSSHRINKSD(%)  4.5391E+01
 EPSSHRINKVR(%)  7.0179E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1713.7286944544746     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -978.57786789073646     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.39
 Elapsed covariance  time in seconds:     5.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1713.729       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.80E-01  9.39E-01  1.51E+00  6.59E-01  1.08E+00  2.24E+00  6.57E-01  8.11E-01  8.84E-01  9.47E-01
 


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
 
         3.52E-02  2.59E-01  1.83E-01  1.33E-01  1.14E-01  6.93E-02  1.96E+00  2.97E-01  8.23E-02  1.08E-01  6.63E-02
 


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
+        1.24E-03
 
 TH 2
+        1.99E-03  6.73E-02
 
 TH 3
+        9.16E-04  7.63E-03  3.35E-02
 
 TH 4
+       -8.92E-04 -3.26E-02  9.37E-04  1.78E-02
 
 TH 5
+        8.93E-04  1.81E-02  1.78E-02 -6.50E-03  1.30E-02
 
 TH 6
+        1.88E-04  1.78E-03 -4.98E-04 -1.11E-03  2.70E-04  4.80E-03
 
 TH 7
+       -1.59E-02 -4.92E-01 -5.58E-02  2.34E-01 -1.32E-01 -1.42E-02  3.86E+00
 
 TH 8
+        1.08E-04  2.32E-03  4.08E-02  4.47E-03  1.95E-02 -9.97E-04 -5.63E-03  8.84E-02
 
 TH 9
+        3.11E-04  1.25E-02 -2.02E-03 -6.44E-03  1.93E-03  7.95E-04 -8.15E-02 -5.20E-03  6.78E-03
 
 TH10
+        7.41E-04  7.82E-04  6.06E-03  6.73E-05  3.94E-03 -5.32E-04 -7.55E-03 -2.07E-03 -2.56E-04  1.17E-02
 
 TH11
+        1.79E-04  2.23E-03 -2.62E-04 -1.22E-03  4.66E-04  1.23E-04 -1.78E-02 -3.98E-03  5.08E-04  7.80E-04  4.40E-03
 
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
+        3.52E-02
 
 TH 2
+        2.18E-01  2.59E-01
 
 TH 3
+        1.42E-01  1.61E-01  1.83E-01
 
 TH 4
+       -1.90E-01 -9.41E-01  3.84E-02  1.33E-01
 
 TH 5
+        2.23E-01  6.12E-01  8.55E-01 -4.28E-01  1.14E-01
 
 TH 6
+        7.70E-02  9.89E-02 -3.92E-02 -1.20E-01  3.42E-02  6.93E-02
 
 TH 7
+       -2.31E-01 -9.66E-01 -1.55E-01  8.93E-01 -5.90E-01 -1.04E-01  1.96E+00
 
 TH 8
+        1.04E-02  3.00E-02  7.50E-01  1.13E-01  5.76E-01 -4.84E-02 -9.65E-03  2.97E-01
 
 TH 9
+        1.07E-01  5.85E-01 -1.34E-01 -5.86E-01  2.06E-01  1.39E-01 -5.04E-01 -2.12E-01  8.23E-02
 
 TH10
+        1.95E-01  2.79E-02  3.06E-01  4.67E-03  3.20E-01 -7.10E-02 -3.56E-02 -6.46E-02 -2.88E-02  1.08E-01
 
 TH11
+        7.68E-02  1.30E-01 -2.16E-02 -1.38E-01  6.17E-02  2.67E-02 -1.37E-01 -2.02E-01  9.31E-02  1.09E-01  6.63E-02
 
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
+        8.97E+02
 
 TH 2
+        7.65E+00  7.62E+02
 
 TH 3
+       -4.99E+01  1.63E+02  6.50E+02
 
 TH 4
+        1.49E+01  5.53E+02 -9.17E+01  8.10E+02
 
 TH 5
+        6.04E+01 -5.56E+02 -1.13E+03 -1.02E+02  2.46E+03
 
 TH 6
+       -3.16E+01  5.65E+01  2.87E+01  4.18E+01 -8.21E+01  2.20E+02
 
 TH 7
+        4.65E+00  4.44E+01 -1.45E+00  1.59E+01 -2.83E-01  2.51E+00  4.69E+00
 
 TH 8
+        5.28E+00 -3.48E+00 -4.47E+01  9.54E+00 -1.52E+01  3.13E+00 -1.62E+00  3.75E+01
 
 TH 9
+       -1.07E+01 -1.48E+02  7.79E+01 -5.66E+01 -1.24E+02 -2.36E+01 -1.25E+01  1.46E+01  2.90E+02
 
 TH10
+       -4.90E+01  7.43E+01  3.18E+01  5.07E+01 -2.14E+02  2.24E+01  2.81E-02  3.27E+01  1.17E+01  1.46E+02
 
 TH11
+       -1.12E+01  1.52E+01 -9.25E+00  1.90E+01 -3.75E+01  3.03E-01  6.28E-01  2.30E+01  5.05E+00  4.90E+00  2.51E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.04
 #CPUT: Total CPU Time in Seconds,       40.168
Stop Time:
Wed Sep 29 15:11:08 CDT 2021
