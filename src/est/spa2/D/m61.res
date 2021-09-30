Thu Sep 30 09:24:28 CDT 2021
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
$DATA ../../../../data/spa2/D/dat61.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   18295.2935472466        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.8340E+02  4.2196E+02  5.2278E+01  1.6835E+02  3.4188E+02 -2.4194E+03 -8.3140E+02 -8.0900E+01 -1.4046E+03 -7.9814E+02
            -3.5609E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -668.142385239331        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.3112E+00  1.1027E+00  8.8032E-01  1.6891E+00  1.0430E+00  2.5529E+00  1.7165E+00  9.6030E-01  1.7646E+00  1.2652E+00
             1.3133E+01
 PARAMETER:  3.7096E-01  1.9776E-01 -2.7467E-02  6.2418E-01  1.4212E-01  1.0372E+00  6.4031E-01  5.9489E-02  6.6790E-01  3.3524E-01
             2.6751E+00
 GRADIENT:  -2.5689E+01 -5.0275E+00 -2.1757E+01  1.6389E+01  3.2222E+01  3.3038E+01 -7.4679E+00  3.5390E+00 -4.0820E+01  1.9723E+01
             2.5592E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -759.928010068516        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.2689E+00  7.6313E-01  8.4940E+00  2.4953E+00  2.8332E+00  3.6836E+00  5.6225E+00  3.2624E-01  4.4579E+00  2.8049E+00
             1.0368E+01
 PARAMETER:  3.3814E-01 -1.7033E-01  2.2394E+00  1.0144E+00  1.1414E+00  1.4039E+00  1.8268E+00 -1.0201E+00  1.5947E+00  1.1314E+00
             2.4387E+00
 GRADIENT:  -2.9959E+00 -1.5920E+00 -4.7318E+00  4.1799E+01 -3.2017E+01  1.5054E+02  2.0963E+01  9.1412E-02  1.2345E+02  4.5438E+01
             2.1346E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -797.590323700762        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.2779E+00  8.8081E-01  8.8727E+00  1.8725E+00  3.4945E+00  2.6951E+00  4.4836E+00  8.9079E-01  2.8030E+00  2.2518E+00
             9.9351E+00
 PARAMETER:  3.4522E-01 -2.6917E-02  2.2830E+00  7.2728E-01  1.3512E+00  1.0914E+00  1.6004E+00 -1.5642E-02  1.1307E+00  9.1175E-01
             2.3961E+00
 GRADIENT:  -1.0143E+01  6.6542E+00  4.9106E+00 -4.5970E+00 -1.7896E+01  2.5104E+01  1.8923E+01 -3.1379E-01 -1.1964E+00  1.2712E+01
             1.8920E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -836.430893709352        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      392
 NPARAMETR:  1.3005E+00  7.5180E-01  1.8436E+01  1.9212E+00  6.5867E+00  2.6940E+00  4.4970E+00  8.2211E+00  2.7315E+00  2.8502E+00
             8.5433E+00
 PARAMETER:  3.6273E-01 -1.8528E-01  3.0143E+00  7.5293E-01  1.9851E+00  1.0910E+00  1.6034E+00  2.2067E+00  1.1048E+00  1.1474E+00
             2.2452E+00
 GRADIENT:  -2.4815E-01  1.3505E+01  1.0197E+00  5.4930E+00 -9.8939E+00 -9.1323E+00 -5.5476E+00  1.3404E+01  1.6340E+01 -3.7134E+00
             1.3733E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -843.743626193222        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      568
 NPARAMETR:  1.2942E+00  3.4978E-01  1.7930E+01  2.1747E+00  6.3896E+00  2.6884E+00  4.6591E+00  5.9725E+00  2.3372E+00  2.9713E+00
             8.4006E+00
 PARAMETER:  3.5787E-01 -9.5046E-01  2.9865E+00  8.7691E-01  1.9547E+00  1.0889E+00  1.6388E+00  1.8872E+00  9.4894E-01  1.1890E+00
             2.2283E+00
 GRADIENT:   2.0707E+00  1.1211E+01  1.2071E+01  1.3893E+01 -1.5211E+01 -5.5798E+00 -1.9064E+00 -2.4992E+01  8.0310E+00 -5.1804E+00
            -1.2161E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -938.550385195560        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      747
 NPARAMETR:  1.1459E+00  1.0000E-02  2.0314E-01  1.3887E+00  1.9448E+01  2.3075E+00  3.4177E+00  1.0456E+00  1.5437E+00  3.6792E+00
             8.7744E+00
 PARAMETER:  2.3618E-01 -5.0631E+00 -1.4939E+00  4.2833E-01  3.0677E+00  9.3617E-01  1.3290E+00  1.4460E-01  5.3419E-01  1.4027E+00
             2.2718E+00
 GRADIENT:   9.0850E+01  0.0000E+00 -6.3958E+01  1.4949E+02  7.9858E-02 -1.6881E+01  6.3847E-03 -8.0618E-01  4.0350E+01  3.4811E-03
            -7.3184E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1020.02901537961        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      923
 NPARAMETR:  5.6834E-01  1.0000E-02  4.5215E-02  4.6197E-01  3.7939E+01  2.1913E+00  1.4418E+00  4.2929E-01  8.5372E-01  2.3262E+00
             8.3861E+00
 PARAMETER: -4.6504E-01 -5.5319E+00 -2.9963E+00 -6.7225E-01  3.7360E+00  8.8449E-01  4.6591E-01 -7.4562E-01 -5.8148E-02  9.4424E-01
             2.2266E+00
 GRADIENT:  -8.8750E+00  0.0000E+00 -5.6678E+00  1.6384E+01 -2.5392E-01  2.7626E+01  1.8123E-01 -1.2161E+00  6.9224E+00  2.7692E-03
             1.3396E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1022.68587099623        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1108
 NPARAMETR:  5.5091E-01  1.0000E-02  3.9496E-02  4.1396E-01  4.1649E+01  2.0043E+00  1.2420E+00  6.5471E-01  8.2205E-01  2.2363E+00
             8.3351E+00
 PARAMETER: -4.9618E-01 -5.4075E+00 -3.1315E+00 -7.8199E-01  3.8293E+00  7.9531E-01  3.1671E-01 -3.2357E-01 -9.5949E-02  9.0481E-01
             2.2205E+00
 GRADIENT:  -1.3270E+00  0.0000E+00 -4.7863E+00 -1.4947E+00 -5.4742E-01  4.6485E+00  3.5478E-01  1.4647E+00  1.5630E+01  4.7111E-03
             8.3748E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1025.68847027377        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1287
 NPARAMETR:  5.5893E-01  1.0000E-02  4.1033E-02  4.2263E-01  3.3117E+02  1.9466E+00  3.3091E-01  1.0463E+00  4.9825E-01  8.2318E-01
             8.2789E+00
 PARAMETER: -4.8173E-01 -5.4075E+00 -3.0934E+00 -7.6125E-01  5.9026E+00  7.6608E-01 -1.0059E+00  1.4522E-01 -5.9665E-01 -9.4585E-02
             2.2137E+00
 GRADIENT:   1.4741E+00  0.0000E+00 -1.0055E+01  5.3604E+00 -2.2330E-02 -2.7493E+00  1.9413E-02  3.2097E+00  2.3652E+00  2.4401E-06
            -3.4266E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1025.96760040923        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1462
 NPARAMETR:  5.6393E-01  1.0000E-02  4.2705E-02  4.3171E-01  2.4846E+03  1.9596E+00  8.0088E-02  1.1589E+00  3.1037E-01  4.6902E-01
             8.3074E+00
 PARAMETER: -4.7283E-01 -5.4075E+00 -3.0534E+00 -7.4000E-01  7.9179E+00  7.7275E-01 -2.4246E+00  2.4747E-01 -1.0700E+00 -6.5711E-01
             2.2171E+00
 GRADIENT:  -3.7086E-01  0.0000E+00 -4.4789E+00 -1.0655E+00 -2.8689E-03 -6.0706E-01  1.5430E-03  1.1638E+00  2.4373E-02  1.7474E-08
            -9.5526E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1025.98169753294        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1643
 NPARAMETR:  5.6460E-01  1.0000E-02  4.2877E-02  4.3267E-01  2.2357E+04  1.9635E+00  2.3548E-02  1.1500E+00  3.1225E-01  4.6322E-01
             8.3135E+00
 PARAMETER: -4.7164E-01 -5.4075E+00 -3.0494E+00 -7.3778E-01  1.0115E+01  7.7474E-01 -3.6487E+00  2.3980E-01 -1.0640E+00 -6.6955E-01
             2.2179E+00
 GRADIENT:  -4.3624E-01  0.0000E+00 -3.1574E+00 -2.0564E+00 -2.9946E-04  4.5791E-02  1.2721E-04  4.6309E-02 -9.7980E-02  1.8678E-10
            -2.9737E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1025.98265138550        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1818
 NPARAMETR:  5.6391E-01  1.0000E-02  4.2938E-02  4.3324E-01  3.9378E+04  1.9591E+00  1.8113E-02  1.1436E+00  3.2815E-01  4.6327E-01
             8.3094E+00
 PARAMETER: -4.7287E-01 -5.4075E+00 -3.0480E+00 -7.3647E-01  1.0681E+01  7.7248E-01 -3.9111E+00  2.3422E-01 -1.0143E+00 -6.6945E-01
             2.2174E+00
 GRADIENT:  -1.3976E+00  0.0000E+00 -3.1333E+00 -1.2900E+00 -1.6429E-04 -8.1087E-01  7.1464E-05  1.4020E-01 -3.2481E-02  5.0947E-11
            -5.9840E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1025.99071231794        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2014             RESET HESSIAN, TYPE I
 NPARAMETR:  5.6595E-01  1.0000E-02  4.2942E-02  4.3365E-01  2.1405E+06  1.9675E+00  1.0000E-02  1.1400E+00  3.3929E-01  4.7252E-01
             8.3146E+00
 PARAMETER: -4.6925E-01 -5.4075E+00 -3.0479E+00 -7.3551E-01  1.4677E+01  7.7674E-01 -2.7378E+01  2.3101E-01 -9.8091E-01 -6.4968E-01
             2.2180E+00
 GRADIENT:   5.2840E+01  0.0000E+00  9.7033E+01  2.7026E+01 -2.6786E-06  3.3858E+01  0.0000E+00  1.3640E+00  5.3341E-01  0.0000E+00
             2.6680E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1025.99158677588        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     2181
 NPARAMETR:  5.6609E-01  1.0000E-02  4.2962E-02  4.3390E-01  2.7026E+06  1.9674E+00  1.0000E-02  1.1383E+00  3.3715E-01  4.7559E-01
             8.3133E+00
 PARAMETER: -4.6901E-01 -5.4075E+00 -3.0474E+00 -7.3493E-01  1.4910E+01  7.7674E-01 -2.7378E+01  2.2954E-01 -9.8723E-01 -6.4320E-01
             2.2179E+00
 GRADIENT:   4.2557E-01  0.0000E+00 -4.5153E+00 -1.1388E-02 -2.2654E-06  6.8992E-01  0.0000E+00  7.8849E-02 -1.1660E-02  0.0000E+00
            -3.7031E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2181
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.5724E-03 -3.4514E-05  7.8616E-03 -1.2485E-02  5.4208E-10
 SE:             2.9172E-02  8.8919E-06  2.3073E-02  1.0209E-02  3.6472E-10
 N:                     100         100         100         100         100

 P VAL.:         8.4851E-01  1.0386E-04  7.3331E-01  2.2137E-01  1.3721E-01

 ETASHRINKSD(%)  2.2688E+00  9.9970E+01  2.2702E+01  6.5798E+01  1.0000E+02
 ETASHRINKVR(%)  4.4862E+00  1.0000E+02  4.0250E+01  8.8302E+01  1.0000E+02
 EBVSHRINKSD(%)  2.1315E+00  9.9959E+01  2.2529E+01  6.6671E+01  1.0000E+02
 EBVSHRINKVR(%)  4.2176E+00  1.0000E+02  3.9982E+01  8.8892E+01  1.0000E+02
 RELATIVEINF(%)  2.9673E+01  1.3554E-05  1.7347E+00  3.1303E-01  0.0000E+00
 EPSSHRINKSD(%)  9.8180E+00
 EPSSHRINKVR(%)  1.8672E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1025.9915867758807     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       76.734653069726392     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    51.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1025.992       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.66E-01  1.00E-02  4.30E-02  4.34E-01  2.70E+06  1.97E+00  1.00E-02  1.14E+00  3.37E-01  4.76E-01  8.31E+00
 


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
+        8.49E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.67E+03  0.00E+00  2.67E+05
 
 TH 4
+       -2.98E+02  0.00E+00 -3.57E+04  5.42E+03
 
 TH 5
+       -7.83E-12  0.00E+00  5.15E-11  6.35E-12 -9.43E-20
 
 TH 6
+        4.68E+00  0.00E+00  1.17E+02 -3.30E+01 -5.47E-12  4.56E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.22E+00  0.00E+00 -2.62E+02 -3.61E+01  2.18E-11  1.73E+00  0.00E+00  5.27E+01
 
 TH 9
+        1.26E-01  0.00E+00  1.83E+02 -4.91E+01 -1.61E-11  1.11E+00  0.00E+00  2.89E+01  2.20E+01
 
 TH10
+        9.47E-03  0.00E+00  1.57E-02 -1.61E-03  5.12E-11 -8.78E-04  0.00E+00  3.09E-03  1.10E-03 -1.42E-03
 
 TH11
+       -1.12E+01  0.00E+00  1.94E+02 -2.53E+01 -7.47E-13  1.07E+00  0.00E+00  3.03E+00  2.48E+00  8.50E-05  9.17E+00
 
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
 #CPUT: Total CPU Time in Seconds,       64.033
Stop Time:
Thu Sep 30 09:25:34 CDT 2021
