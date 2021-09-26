Sat Sep 25 09:48:09 CDT 2021
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
$DATA ../../../../data/spa/S1/dat28.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1648.12027112526        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.5656E+01 -1.1192E+02 -5.2622E+01 -1.1259E+02  8.6163E+01  1.3470E+01  3.4828E-01  4.8232E+00 -1.6287E+01  2.0923E+01
            -1.5952E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1658.99890283075        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0067E+00  1.0118E+00  1.1089E+00  1.0899E+00  9.3055E-01  9.7496E-01  8.8074E-01  1.0140E+00  1.0957E+00  7.2858E-01
             1.0340E+00
 PARAMETER:  1.0669E-01  1.1177E-01  2.0334E-01  1.8605E-01  2.8023E-02  7.4644E-02 -2.6989E-02  1.1386E-01  1.9144E-01 -2.1666E-01
             1.3342E-01
 GRADIENT:   9.2753E+01  3.1072E+01  3.4717E+01  1.6973E+01 -4.1886E+01  3.0419E+00  2.1695E+00 -8.0226E+00  9.6527E+00 -1.3089E+01
             3.5125E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1661.35572124131        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.8569E-01  1.1006E+00  1.2751E+00  1.0267E+00  1.0428E+00  9.5244E-01  5.5573E-01  1.3183E+00  1.2188E+00  9.2094E-01
             1.0107E+00
 PARAMETER:  8.5584E-02  1.9590E-01  3.4300E-01  1.2633E-01  1.4189E-01  5.1276E-02 -4.8747E-01  3.7633E-01  2.9785E-01  1.7635E-02
             1.1061E-01
 GRADIENT:   4.5008E+01  1.8322E+01  1.2367E+01  1.0439E+01 -1.9804E+01 -4.3614E+00  1.2957E+00 -2.6344E+00  1.0039E+01  2.1101E+00
            -1.0639E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1661.77671724779        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.7311E-01  1.0412E+00  1.2202E+00  1.0571E+00  1.0112E+00  9.5977E-01  6.2675E-01  1.2690E+00  1.1385E+00  8.8028E-01
             1.0079E+00
 PARAMETER:  7.2744E-02  1.4036E-01  2.9902E-01  1.5554E-01  1.1112E-01  5.8939E-02 -3.6721E-01  3.3825E-01  2.2974E-01 -2.7516E-02
             1.0787E-01
 GRADIENT:   1.4528E+01  3.5391E+00  2.6407E+00  7.0892E-01 -4.4691E+00 -1.3178E+00  1.0213E+00  4.9592E-01  2.5713E+00  3.2299E-01
            -1.9983E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1661.77717752514        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  9.7220E-01  1.0357E+00  1.2193E+00  1.0603E+00  1.0093E+00  9.6030E-01  6.2151E-01  1.2630E+00  1.1352E+00  8.7968E-01
             1.0086E+00
 PARAMETER:  7.1808E-02  1.3509E-01  2.9827E-01  1.5854E-01  1.0926E-01  5.9487E-02 -3.7561E-01  3.3351E-01  2.2679E-01 -2.8195E-02
             1.0858E-01
 GRADIENT:   1.2298E+01  3.0379E+00  2.2558E+00  5.9553E-01 -3.8165E+00 -1.1234E+00  8.8304E-01  4.1230E-01  2.1930E+00  2.9420E-01
            -1.7035E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1661.77722051129        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      395
 NPARAMETR:  9.7197E-01  1.0342E+00  1.2194E+00  1.0612E+00  1.0089E+00  9.6043E-01  6.1960E-01  1.2619E+00  1.1344E+00  8.7959E-01
             1.0088E+00
 PARAMETER:  7.1573E-02  1.3361E-01  2.9836E-01  1.5938E-01  1.0885E-01  5.9628E-02 -3.7868E-01  3.3259E-01  2.2606E-01 -2.8297E-02
             1.0876E-01
 GRADIENT:  -2.5264E+01 -1.4238E+00  1.7183E+00 -8.1195E+00 -4.3706E+00 -4.2681E+00  4.3245E-01  3.0456E-01  6.5359E-01  2.2141E-01
            -1.7146E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1662.45802451881        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      573
 NPARAMETR:  9.7956E-01  7.7720E-01  1.3389E+00  1.2380E+00  9.5930E-01  9.6470E-01  5.2226E-01  1.2391E+00  1.0069E+00  8.5432E-01
             1.0114E+00
 PARAMETER:  7.9346E-02 -1.5205E-01  3.9184E-01  3.1346E-01  5.8444E-02  6.4059E-02 -5.4960E-01  3.1442E-01  1.0688E-01 -5.7450E-02
             1.1137E-01
 GRADIENT:  -2.1085E+00  5.6241E+00  1.4873E+00  7.2148E+00 -2.5084E+00 -1.5892E+00 -2.4729E-01  8.0217E-02 -1.0909E+00 -1.5698E+00
            -1.2676E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1662.79819770675        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      748
 NPARAMETR:  9.7821E-01  5.2860E-01  1.3464E+00  1.3944E+00  8.8773E-01  9.6598E-01  4.3340E-01  1.1619E+00  8.9891E-01  8.2653E-01
             1.0128E+00
 PARAMETER:  7.7971E-02 -5.3752E-01  3.9742E-01  4.3249E-01 -1.9088E-02  6.5384E-02 -7.3610E-01  2.5003E-01 -6.5683E-03 -9.0520E-02
             1.1271E-01
 GRADIENT:   3.9157E-01  3.3235E+00 -7.6566E-01  9.8681E+00  6.8074E-01 -4.2857E-01 -1.5783E-01 -4.4771E-01 -2.1314E+00 -5.1204E-01
            -6.9696E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1662.90135436308        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      925
 NPARAMETR:  9.7666E-01  4.1349E-01  1.3628E+00  1.4623E+00  8.5965E-01  9.6628E-01  4.0258E-01  1.1731E+00  8.5761E-01  8.1000E-01
             1.0139E+00
 PARAMETER:  7.6384E-02 -7.8312E-01  4.0956E-01  4.7999E-01 -5.1231E-02  6.5697E-02 -8.0985E-01  2.5962E-01 -5.3601E-02 -1.1072E-01
             1.1383E-01
 GRADIENT:   2.9167E-01  8.9708E-01  3.2313E-01  2.9599E+00 -8.2204E-01  5.6324E-02 -5.5414E-02 -3.1916E-02 -4.5819E-01 -1.2700E-01
            -9.6900E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1662.92090769087        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1101
 NPARAMETR:  9.7557E-01  3.4104E-01  1.3602E+00  1.5046E+00  8.3948E-01  9.6570E-01  4.0373E-01  1.1737E+00  8.3157E-01  7.9657E-01
             1.0142E+00
 PARAMETER:  7.5263E-02 -9.7577E-01  4.0761E-01  5.0855E-01 -7.4970E-02  6.5100E-02 -8.0701E-01  2.6014E-01 -8.4445E-02 -1.2744E-01
             1.1411E-01
 GRADIENT:   1.2554E-01 -1.1181E-02 -1.5763E-02 -4.0677E-01  3.8596E-03  3.9806E-02 -3.0631E-02  1.9251E-02 -5.0772E-03 -5.5967E-02
             1.3050E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1662.92257774873        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1277
 NPARAMETR:  9.7496E-01  3.0261E-01  1.3613E+00  1.5283E+00  8.2976E-01  9.6531E-01  4.1850E-01  1.1778E+00  8.1739E-01  7.8963E-01
             1.0142E+00
 PARAMETER:  7.4639E-02 -1.0953E+00  4.0846E-01  5.2418E-01 -8.6618E-02  6.4694E-02 -7.7107E-01  2.6362E-01 -1.0164E-01 -1.3619E-01
             1.1415E-01
 GRADIENT:   4.5915E-02  4.0251E-02 -1.3380E-02  2.5635E-01  3.1202E-02 -3.4466E-03 -2.3759E-02 -3.0165E-03 -8.9416E-02 -2.8692E-02
            -1.8815E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1662.92269065820        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1452
 NPARAMETR:  9.7477E-01  2.9074E-01  1.3600E+00  1.5353E+00  8.2624E-01  9.6523E-01  4.3278E-01  1.1779E+00  8.1323E-01  7.8704E-01
             1.0143E+00
 PARAMETER:  7.4449E-02 -1.1353E+00  4.0748E-01  5.2870E-01 -9.0875E-02  6.4613E-02 -7.3753E-01  2.6375E-01 -1.0674E-01 -1.3947E-01
             1.1419E-01
 GRADIENT:   3.6812E-02 -2.7004E-02 -5.1707E-02 -1.2873E-01  9.3649E-02 -4.2837E-04 -2.2049E-02  3.2050E-03 -3.6060E-02 -2.2356E-02
            -6.8491E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1662.92281750650        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1627
 NPARAMETR:  9.7451E-01  2.7392E-01  1.3581E+00  1.5452E+00  8.2121E-01  9.6512E-01  4.6939E-01  1.1783E+00  8.0719E-01  7.8319E-01
             1.0144E+00
 PARAMETER:  7.4182E-02 -1.1949E+00  4.0606E-01  5.3516E-01 -9.6975E-02  6.4493E-02 -6.5633E-01  2.6404E-01 -1.1419E-01 -1.4438E-01
             1.1425E-01
 GRADIENT:   3.9649E-02 -6.6659E-02 -6.6059E-02 -3.6767E-01  1.0304E-01  2.6536E-03 -2.0699E-02  7.9856E-03  4.9482E-04 -2.0413E-02
             5.1715E-04

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1662.92332188658        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1804
 NPARAMETR:  9.7380E-01  2.3003E-01  1.3518E+00  1.5716E+00  8.0769E-01  9.6479E-01  6.9197E-01  1.1787E+00  7.9084E-01  7.7231E-01
             1.0145E+00
 PARAMETER:  7.3448E-02 -1.3695E+00  4.0144E-01  5.5207E-01 -1.1357E-01  6.4158E-02 -2.6822E-01  2.6444E-01 -1.3466E-01 -1.5836E-01
             1.1438E-01
 GRADIENT:   6.2462E-03 -1.7719E-02  2.4179E-02  6.0464E-02 -1.3808E-01 -2.7450E-03 -1.7914E-02  1.3686E-02 -1.1215E-02 -2.3245E-03
             1.7290E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1662.92358860837        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1981
 NPARAMETR:  9.7356E-01  2.1492E-01  1.3475E+00  1.5803E+00  8.0249E-01  9.6471E-01  8.7229E-01  1.1772E+00  7.8477E-01  7.6773E-01
             1.0145E+00
 PARAMETER:  7.3200E-02 -1.4375E+00  3.9828E-01  5.5759E-01 -1.2004E-01  6.4074E-02 -3.6638E-02  2.6310E-01 -1.4236E-01 -1.6432E-01
             1.1443E-01
 GRADIENT:   3.1685E-02 -8.2291E-02 -7.5764E-02 -4.3435E-01  1.5720E-01  4.3139E-03 -9.8799E-03  8.6405E-03  3.6489E-02 -8.6068E-03
             6.0453E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1662.93307012517        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2166
 NPARAMETR:  9.7423E-01  2.5771E-01  1.3383E+00  1.5552E+00  8.0959E-01  9.6511E-01  9.0345E-01  1.1634E+00  7.9676E-01  7.7217E-01
             1.0144E+00
 PARAMETER:  7.3897E-02 -1.2559E+00  3.9143E-01  5.4161E-01 -1.1122E-01  6.4491E-02 -1.5375E-03  2.5134E-01 -1.2720E-01 -1.5855E-01
             1.1431E-01
 GRADIENT:  -8.8611E-02  1.9124E-01  2.5362E-01  1.6138E+00 -6.6998E-01  1.5080E-02 -1.7354E-02  3.1700E-03 -2.7199E-01 -3.4355E-02
            -2.5074E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1662.94381167584        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2344
 NPARAMETR:  9.7520E-01  3.1681E-01  1.3289E+00  1.5188E+00  8.2169E-01  9.6557E-01  8.7824E-01  1.1475E+00  8.1549E-01  7.8070E-01
             1.0142E+00
 PARAMETER:  7.4884E-02 -1.0495E+00  3.8433E-01  5.1790E-01 -9.6397E-02  6.4965E-02 -2.9839E-02  2.3762E-01 -1.0397E-01 -1.4756E-01
             1.1413E-01
 GRADIENT:  -3.9541E-02  5.4942E-02  1.4601E-02  2.9205E-01 -1.3071E-01 -2.0516E-04 -1.0730E-02 -3.0661E-02 -1.9214E-01 -6.8410E-03
            -1.2137E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1662.94459109966        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2521
 NPARAMETR:  9.7537E-01  3.2777E-01  1.3302E+00  1.5122E+00  8.2501E-01  9.6567E-01  8.7278E-01  1.1481E+00  8.1945E-01  7.8276E-01
             1.0142E+00
 PARAMETER:  7.5064E-02 -1.0154E+00  3.8534E-01  5.1353E-01 -9.2361E-02  6.5064E-02 -3.6074E-02  2.3811E-01 -9.9118E-02 -1.4493E-01
             1.1409E-01
 GRADIENT:  -5.1551E-03  4.1350E-03  2.5376E-02 -9.5980E-03 -7.1181E-02  2.7937E-03  5.1603E-04  4.6972E-03  2.0468E-02 -2.0385E-03
             2.6049E-03

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1662.94459164040        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     2686
 NPARAMETR:  9.7537E-01  3.2770E-01  1.3305E+00  1.5122E+00  8.2510E-01  9.6566E-01  8.7163E-01  1.1483E+00  8.1943E-01  7.8285E-01
             1.0142E+00
 PARAMETER:  7.5065E-02 -1.0157E+00  3.8554E-01  5.1357E-01 -9.2250E-02  6.5059E-02 -3.7388E-02  2.3829E-01 -9.9150E-02 -1.4482E-01
             1.1409E-01
 GRADIENT:  -1.0887E-02  2.1225E-04  1.4783E-02  6.3569E-02 -4.1069E-02  8.0409E-04  6.6752E-04  2.7663E-03  1.2116E-02 -1.8097E-03
             1.8978E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2686
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0033E-04 -7.4494E-03 -2.7034E-02 -2.9946E-03 -3.1989E-02
 SE:             2.9845E-02  5.0242E-03  1.8486E-02  2.8931E-02  2.0052E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9197E-01  1.3815E-01  1.4363E-01  9.1756E-01  1.1065E-01

 ETASHRINKSD(%)  1.7128E-02  8.3168E+01  3.8069E+01  3.0783E+00  3.2822E+01
 ETASHRINKVR(%)  3.4253E-02  9.7167E+01  6.1645E+01  6.0618E+00  5.4871E+01
 EBVSHRINKSD(%)  4.4135E-01  8.3953E+01  4.0060E+01  3.4022E+00  3.1336E+01
 EBVSHRINKVR(%)  8.8075E-01  9.7425E+01  6.4072E+01  6.6887E+00  5.2852E+01
 RELATIVEINF(%)  9.7222E+01  1.1829E-01  6.6545E+00  5.8834E+00  4.5567E+00
 EPSSHRINKSD(%)  4.4731E+01
 EPSSHRINKVR(%)  6.9453E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1662.9445916404006     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -927.79376507666245     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.20
 Elapsed covariance  time in seconds:     5.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1662.945       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  3.28E-01  1.33E+00  1.51E+00  8.25E-01  9.66E-01  8.72E-01  1.15E+00  8.19E-01  7.83E-01  1.01E+00
 


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
 
         2.91E-02  5.53E-01  3.41E-01  3.20E-01  2.46E-01  5.69E-02  2.17E-01  2.74E-01  1.62E-01  2.32E-01  6.02E-02
 


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
+        8.47E-04
 
 TH 2
+        3.70E-03  3.06E-01
 
 TH 3
+        1.61E-03  1.30E-01  1.16E-01
 
 TH 4
+       -1.94E-03 -1.75E-01 -7.20E-02  1.02E-01
 
 TH 5
+        1.46E-03  1.26E-01  7.52E-02 -7.15E-02  6.06E-02
 
 TH 6
+       -3.00E-04  1.05E-03 -7.16E-05 -7.51E-04  6.17E-06  3.24E-03
 
 TH 7
+        2.15E-03  1.00E-01  4.91E-02 -5.69E-02  4.41E-02  1.60E-03  4.70E-02
 
 TH 8
+        7.78E-04  9.02E-02  7.71E-02 -5.00E-02  5.11E-02 -8.68E-04  2.68E-02  7.52E-02
 
 TH 9
+        1.06E-03  8.45E-02  3.50E-02 -4.85E-02  3.48E-02 -2.35E-04  2.86E-02  2.42E-02  2.63E-02
 
 TH10
+        1.83E-03  1.06E-01  6.02E-02 -6.03E-02  4.97E-02  5.47E-04  4.64E-02  3.61E-02  2.91E-02  5.41E-02
 
 TH11
+       -1.22E-04 -2.47E-03 -6.47E-04  1.33E-03 -6.17E-04 -2.16E-04 -1.53E-03  5.80E-04 -2.81E-04 -1.40E-03  3.63E-03
 
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
+        2.91E-02
 
 TH 2
+        2.30E-01  5.53E-01
 
 TH 3
+        1.63E-01  6.91E-01  3.41E-01
 
 TH 4
+       -2.08E-01 -9.92E-01 -6.61E-01  3.20E-01
 
 TH 5
+        2.03E-01  9.29E-01  8.97E-01 -9.08E-01  2.46E-01
 
 TH 6
+       -1.81E-01  3.34E-02 -3.69E-03 -4.13E-02  4.40E-04  5.69E-02
 
 TH 7
+        3.41E-01  8.35E-01  6.65E-01 -8.21E-01  8.26E-01  1.30E-01  2.17E-01
 
 TH 8
+        9.75E-02  5.95E-01  8.25E-01 -5.70E-01  7.57E-01 -5.56E-02  4.51E-01  2.74E-01
 
 TH 9
+        2.24E-01  9.41E-01  6.34E-01 -9.36E-01  8.72E-01 -2.54E-02  8.12E-01  5.44E-01  1.62E-01
 
 TH10
+        2.71E-01  8.26E-01  7.60E-01 -8.11E-01  8.69E-01  4.13E-02  9.21E-01  5.67E-01  7.72E-01  2.32E-01
 
 TH11
+       -6.93E-02 -7.42E-02 -3.15E-02  6.92E-02 -4.16E-02 -6.28E-02 -1.17E-01  3.51E-02 -2.87E-02 -1.00E-01  6.02E-02
 
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
+        1.60E+03
 
 TH 2
+       -1.85E+02  4.10E+02
 
 TH 3
+       -8.71E+01  1.23E+02  2.05E+02
 
 TH 4
+       -1.74E+02  4.57E+02  1.43E+01  7.25E+02
 
 TH 5
+        3.19E+02 -4.27E+02 -4.83E+02 -1.31E+02  1.47E+03
 
 TH 6
+        2.56E+02 -4.73E+01 -4.61E+01  9.45E+00  1.43E+02  4.00E+02
 
 TH 7
+       -2.04E+02  7.17E+00  1.22E+01 -1.38E+01 -5.91E+01 -1.14E+02  2.38E+02
 
 TH 8
+       -1.50E+01 -3.89E+00 -2.61E+01 -4.67E+00 -2.20E+01 -3.56E+00  3.00E+01  4.97E+01
 
 TH 9
+        6.85E+01 -6.93E+01  2.62E+01  2.63E+01 -6.27E+01  9.63E+01 -9.55E+01 -6.03E+00  3.93E+02
 
 TH10
+        6.45E+01  2.24E+00 -1.23E+01  2.10E+01 -4.03E+01  5.83E+01 -1.53E+02 -3.59E+00  6.39E+01  1.85E+02
 
 TH11
+       -7.86E+00  5.06E+01  3.39E+01  2.54E+01 -1.01E+02 -4.56E+00  1.75E+01 -7.22E+00 -3.91E+01  2.81E+00  2.96E+02
 
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
 #CPUT: Total CPU Time in Seconds,       34.564
Stop Time:
Sat Sep 25 09:48:45 CDT 2021
