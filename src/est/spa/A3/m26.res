Sat Sep 25 09:09:40 CDT 2021
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
$DATA ../../../../data/spa/A3/dat26.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -55.8037417547771        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.5487E+01 -7.5855E+00  1.4336E+02 -1.5677E+02  5.0986E+01  4.6300E+01 -2.9622E+01 -1.3230E+02 -1.2927E+02 -1.2036E+02
            -2.7805E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1219.34407633748        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0790E+00  1.0179E+00  9.1510E-01  1.1429E+00  1.0290E+00  7.4809E-01  9.5092E-01  1.0516E+00  9.6498E-01  1.0436E+00
             5.3759E+00
 PARAMETER:  1.7607E-01  1.1773E-01  1.1274E-02  2.3355E-01  1.2858E-01 -1.9024E-01  4.9676E-02  1.5035E-01  6.4352E-02  1.4271E-01
             1.7819E+00
 GRADIENT:   1.4451E+02 -3.5018E+01 -2.1406E+01 -3.1092E+01  4.2641E+00 -9.1830E+00  1.1581E+01  7.1571E+00  2.5236E+01  1.8065E+01
             1.9439E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1246.51533681226        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0121E+00  1.4168E+00  1.4064E+00  9.6177E-01  1.5702E+00  7.5846E-01  1.3568E-01  2.2080E+00  9.5959E-01  1.8061E+00
             4.6857E+00
 PARAMETER:  1.1202E-01  4.4842E-01  4.4104E-01  6.1017E-02  5.5123E-01 -1.7647E-01 -1.8974E+00  8.9211E-01  5.8751E-02  6.9114E-01
             1.6445E+00
 GRADIENT:  -4.3921E+01  5.9466E+01 -9.6943E+00  6.1507E+01 -3.1434E+01 -1.0817E+01 -1.9906E-01  5.6149E+00 -1.9343E+00  3.0955E+01
             8.8961E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1270.33830544746        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.0035E+00  1.1350E+00  1.3653E+00  1.0477E+00  1.2969E+00  7.9405E-01  3.1961E-01  8.6454E-01  9.5312E-01  6.8562E-01
             3.9909E+00
 PARAMETER:  1.0353E-01  2.2668E-01  4.1135E-01  1.4660E-01  3.6000E-01 -1.3061E-01 -1.0407E+00 -4.5563E-02  5.1985E-02 -2.7743E-01
             1.4840E+00
 GRADIENT:   4.4321E+00  2.3811E+00 -1.9375E+00  2.5597E+00 -4.1288E+00  3.6723E+00 -3.4360E-01  1.9518E+00  2.4035E-01  2.9678E+00
            -5.0640E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1272.72615430705        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.9881E-01  8.7700E-01  1.3577E+00  1.2086E+00  1.1709E+00  7.8048E-01  7.9149E-01  2.1744E-01  7.3083E-01  2.5071E-01
             4.0584E+00
 PARAMETER:  9.8808E-02 -3.1248E-02  4.0576E-01  2.8944E-01  2.5777E-01 -1.4785E-01 -1.3384E-01 -1.4258E+00 -2.1358E-01 -1.2835E+00
             1.5008E+00
 GRADIENT:  -3.4849E+00 -2.4440E-02 -3.3917E-01 -1.4939E+00 -3.3844E-01 -2.5849E-01  1.6052E-01  1.6776E-01 -1.0541E+00  3.5384E-01
             5.1624E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1273.09966868724        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      367
 NPARAMETR:  9.9799E-01  7.0876E-01  1.3444E+00  1.3115E+00  1.0965E+00  7.8174E-01  8.1475E-01  4.4269E-02  7.2124E-01  7.8138E-02
             4.0319E+00
 PARAMETER:  9.7984E-02 -2.4424E-01  3.9596E-01  3.7116E-01  1.9210E-01 -1.4624E-01 -1.0488E-01 -3.0175E+00 -2.2679E-01 -2.4493E+00
             1.4942E+00
 GRADIENT:   2.3665E-01 -1.6868E-01 -2.2549E-01 -4.9194E-01  1.5507E-01  1.2247E-01 -1.3002E-01  8.7925E-03 -2.5336E-01  2.7605E-02
             2.4485E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1273.18627335520        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      486
 NPARAMETR:  9.9682E-01  5.4163E-01  1.3643E+00  1.4228E+00  1.0392E+00  7.8132E-01  1.0782E+00  1.0000E-02  6.6847E-01  1.5414E-02
             4.0358E+00
 PARAMETER:  9.6814E-02 -5.1318E-01  4.1065E-01  4.5260E-01  1.3842E-01 -1.4677E-01  1.7527E-01 -5.1702E+00 -3.0277E-01 -4.0725E+00
             1.4952E+00
 GRADIENT:  -2.6810E+00  3.3480E+00  7.4332E-01  1.1008E+01 -2.2804E+00 -4.0650E-01  3.8219E-02  0.0000E+00 -4.5600E-01  1.0708E-03
            -6.5203E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1273.25393209310        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      665
 NPARAMETR:  9.9585E-01  4.5901E-01  1.3648E+00  1.4680E+00  1.0121E+00  7.8140E-01  1.1315E+00  1.0000E-02  6.5931E-01  1.0000E-02
             4.0336E+00
 PARAMETER:  9.5840E-02 -6.7868E-01  4.1099E-01  4.8388E-01  1.1205E-01 -1.4667E-01  2.2358E-01 -6.7122E+00 -3.1656E-01 -5.2409E+00
             1.4947E+00
 GRADIENT:  -5.3528E-01  1.0518E+00  2.5432E-01  3.9141E+00 -7.9634E-01 -8.4350E-02 -2.9216E-02  0.0000E+00 -1.0415E-01  0.0000E+00
            -2.3149E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1273.29450358355        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      843
 NPARAMETR:  9.9419E-01  3.5900E-01  1.3746E+00  1.5257E+00  9.8386E-01  7.8090E-01  1.2644E+00  1.0000E-02  6.4220E-01  1.0000E-02
             4.0322E+00
 PARAMETER:  9.4176E-02 -9.2442E-01  4.1819E-01  5.2244E-01  8.3729E-02 -1.4730E-01  3.3457E-01 -9.0786E+00 -3.4286E-01 -7.0441E+00
             1.4943E+00
 GRADIENT:  -1.1809E-01  1.3988E-01  9.2996E-02  1.4764E-01 -1.6369E-01 -1.3373E-02 -6.1740E-02  0.0000E+00 -1.0029E-01  0.0000E+00
            -6.1795E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1273.32751030995        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1020
 NPARAMETR:  9.9245E-01  2.4903E-01  1.3582E+00  1.5906E+00  9.4209E-01  7.8068E-01  1.6955E+00  1.0000E-02  6.2233E-01  1.0000E-02
             4.0299E+00
 PARAMETER:  9.2418E-02 -1.2902E+00  4.0618E-01  5.6409E-01  4.0350E-02 -1.4759E-01  6.2800E-01 -1.2804E+01 -3.7428E-01 -9.8873E+00
             1.4937E+00
 GRADIENT:  -3.6568E-01  3.3755E-01  1.2873E-01  2.2289E+00 -4.8495E-01 -2.4560E-02  4.3134E-02  0.0000E+00 -1.8739E-02  0.0000E+00
            -1.0248E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1273.34563646019        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1196
 NPARAMETR:  9.9109E-01  1.7540E-01  1.3688E+00  1.6347E+00  9.2469E-01  7.8030E-01  1.9665E+00  1.0000E-02  6.1259E-01  1.0000E-02
             4.0292E+00
 PARAMETER:  9.1046E-02 -1.6407E+00  4.1397E-01  5.9145E-01  2.1704E-02 -1.4808E-01  7.7625E-01 -1.6505E+01 -3.9006E-01 -1.2723E+01
             1.4936E+00
 GRADIENT:  -6.5531E-01  3.0645E-01  2.9271E-01  3.4662E+00 -9.2123E-01 -4.8869E-02 -1.2626E-02  0.0000E+00 -1.2418E-01  0.0000E+00
            -2.2577E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1273.35851244745        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1373
 NPARAMETR:  9.8980E-01  1.0290E-01  1.3923E+00  1.6788E+00  9.1469E-01  7.7984E-01  2.3231E+00  1.0000E-02  6.0402E-01  1.0000E-02
             4.0285E+00
 PARAMETER:  8.9744E-02 -2.1740E+00  4.3092E-01  6.1805E-01  1.0831E-02 -1.4867E-01  9.4290E-01 -2.2308E+01 -4.0415E-01 -1.7177E+01
             1.4934E+00
 GRADIENT:  -5.1550E-01  1.6238E-01  2.1848E-01  3.7840E+00 -7.6270E-01 -5.7383E-02 -3.5686E-02  0.0000E+00 -1.4309E-02  0.0000E+00
            -2.2410E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1273.36939554843        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1552
 NPARAMETR:  9.8890E-01  5.9964E-02  1.3992E+00  1.7033E+00  9.0612E-01  7.7946E-01  3.3517E+00  1.0000E-02  5.9590E-01  1.0000E-02
             4.0278E+00
 PARAMETER:  8.8842E-02 -2.7140E+00  4.3590E-01  6.3258E-01  1.4172E-03 -1.4915E-01  1.3095E+00 -2.8339E+01 -4.1769E-01 -2.1815E+01
             1.4932E+00
 GRADIENT:  -4.1481E-01  5.6790E-02  8.9760E-02  2.0771E+00 -3.7926E-01 -7.1419E-02 -8.4460E-03  0.0000E+00 -1.5180E-01  0.0000E+00
            -3.0591E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1273.37214296404        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1728
 NPARAMETR:  9.8851E-01  3.8364E-02  1.3995E+00  1.7154E+00  9.0038E-01  7.7941E-01  4.7561E+00  1.0000E-02  5.9301E-01  1.0000E-02
             4.0282E+00
 PARAMETER:  8.8442E-02 -3.1606E+00  4.3610E-01  6.3966E-01 -4.9371E-03 -1.4922E-01  1.6594E+00 -3.3395E+01 -4.2254E-01 -2.5706E+01
             1.4933E+00
 GRADIENT:  -2.1135E-01  3.6561E-02  8.9813E-02  1.1713E+00 -2.5428E-01 -1.5951E-02  2.3035E-02  0.0000E+00 -1.8559E-02  0.0000E+00
            -1.7676E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1273.37574763791        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1906
 NPARAMETR:  9.8805E-01  1.8525E-02  1.3994E+00  1.7268E+00  8.9493E-01  7.7928E-01  6.2132E+00  1.0000E-02  5.9135E-01  1.0000E-02
             4.0276E+00
 PARAMETER:  8.7983E-02 -3.8886E+00  4.3602E-01  6.4624E-01 -1.1006E-02 -1.4938E-01  1.9267E+00 -4.1671E+01 -4.2534E-01 -3.2066E+01
             1.4932E+00
 GRADIENT:  -3.4563E-01  1.1702E-02  1.2222E-01  1.1979E+00 -3.5668E-01 -2.8158E-02 -6.8592E-04  0.0000E+00 -4.6634E-02  0.0000E+00
            -1.0333E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1273.37665812852        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2082
 NPARAMETR:  9.8794E-01  1.0219E-02  1.4041E+00  1.7316E+00  8.9499E-01  7.7924E-01  7.1646E+00  1.0000E-02  5.9027E-01  1.0000E-02
             4.0279E+00
 PARAMETER:  8.7870E-02 -4.4835E+00  4.3937E-01  6.4905E-01 -1.0939E-02 -1.4944E-01  2.0691E+00 -4.8457E+01 -4.2718E-01 -3.7282E+01
             1.4932E+00
 GRADIENT:  -7.5564E-02  1.3314E-03  3.3316E-02  3.6696E-01 -1.1194E-01 -3.5120E-03 -4.0164E-03  0.0000E+00 -3.5339E-02  0.0000E+00
            -1.9528E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1273.37708940101        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2260
 NPARAMETR:  9.8795E-01  1.0000E-02  1.4060E+00  1.7317E+00  8.9588E-01  7.7921E-01  8.5209E+00  1.0000E-02  5.9005E-01  1.0000E-02
             4.0278E+00
 PARAMETER:  8.7882E-02 -4.7309E+00  4.4072E-01  6.4912E-01 -9.9542E-03 -1.4947E-01  2.2425E+00 -5.1301E+01 -4.2755E-01 -3.9472E+01
             1.4932E+00
 GRADIENT:   3.4144E-02  0.0000E+00 -1.4114E-02 -8.9419E-02  4.0586E-02  7.6436E-04  2.2899E-04  0.0000E+00  9.0259E-03  0.0000E+00
             1.5088E-03

0ITERATION NO.:   82    OBJECTIVE VALUE:  -1273.37709457267        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     2317
 NPARAMETR:  9.8795E-01  1.0000E-02  1.4050E+00  1.7317E+00  8.9541E-01  7.7922E-01  8.4717E+00  1.0000E-02  5.9009E-01  1.0000E-02
             4.0278E+00
 PARAMETER:  8.7873E-02 -4.7248E+00  4.4002E-01  6.4909E-01 -1.0470E-02 -1.4946E-01  2.2367E+00 -5.1232E+01 -4.2749E-01 -3.9418E+01
             1.4932E+00
 GRADIENT:  -3.1034E-03  0.0000E+00  1.0531E-03  1.5788E-02 -3.5510E-03 -5.9186E-04 -8.2979E-05  0.0000E+00 -8.7357E-04  0.0000E+00
            -5.8005E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2317
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.6903E-04 -9.8313E-04  8.0377E-05 -1.5000E-02 -4.0717E-06
 SE:             2.7824E-02  1.2293E-03  7.4996E-05  2.2329E-02  1.3184E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7508E-01  4.2386E-01  2.8383E-01  5.0172E-01  9.7536E-01

 ETASHRINKSD(%)  6.7844E+00  9.5882E+01  9.9749E+01  2.5195E+01  9.9558E+01
 ETASHRINKVR(%)  1.3109E+01  9.9830E+01  9.9999E+01  4.4042E+01  9.9998E+01
 EBVSHRINKSD(%)  6.6511E+00  9.5872E+01  9.9706E+01  2.5283E+01  9.9537E+01
 EBVSHRINKVR(%)  1.2860E+01  9.9830E+01  9.9999E+01  4.4173E+01  9.9998E+01
 RELATIVEINF(%)  4.3805E+01  4.3678E-04  2.8634E-05  1.7660E-01  5.2010E-05
 EPSSHRINKSD(%)  1.7806E+01
 EPSSHRINKVR(%)  3.2441E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1273.3770945726656     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -538.22626800892738     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.40
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1273.377       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  1.00E-02  1.40E+00  1.73E+00  8.95E-01  7.79E-01  8.47E+00  1.00E-02  5.90E-01  1.00E-02  4.03E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.85E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -8.12E+00  0.00E+00  4.57E+01
 
 TH 4
+       -1.59E+02  0.00E+00  2.14E+00  5.64E+02
 
 TH 5
+        3.72E+01  0.00E+00 -1.09E+02 -1.31E+02  2.90E+02
 
 TH 6
+        2.16E+02  0.00E+00  6.05E+00 -2.69E+00 -1.93E+01  2.93E+02
 
 TH 7
+        8.61E-02  0.00E+00  1.56E-02 -1.85E-03 -3.86E-02  9.17E-02  3.31E-05
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.85E+00  0.00E+00 -5.55E+00 -5.14E+01  2.47E+01  3.88E+00 -5.37E-04  0.00E+00  5.54E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.32E+01  0.00E+00 -1.40E+00 -1.42E+01  6.53E+00  9.85E+00  2.31E-03  0.00E+00  1.84E+00  0.00E+00  1.11E+00
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.65E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -3.44E+00  0.00E+00  4.40E+01
 
 TH 4
+       -1.76E+02  0.00E+00  4.49E+00  5.86E+02
 
 TH 5
+        3.47E+01  0.00E+00 -1.05E+02 -1.38E+02  2.80E+02
 
 TH 6
+       -2.82E+01  0.00E+00  2.79E+00 -3.10E+01 -3.38E+00  2.60E+02
 
 TH 7
+        1.31E-02  0.00E+00  5.07E-03 -1.66E-03 -3.83E-02  7.61E-02  6.12E-04
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.37E+00  0.00E+00  1.88E+00 -3.03E+01  2.10E+01  2.00E+00  6.52E-02  0.00E+00  1.69E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.37E+01  0.00E+00  2.20E-01 -1.29E+01  4.22E+00  5.79E+00  4.12E-03  0.00E+00  1.84E+01  0.00E+00  2.85E+01
 
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
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.63E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        1.67E+00  0.00E+00  3.60E+01
 
 TH 4
+       -1.65E+02  0.00E+00  4.41E+00  5.97E+02
 
 TH 5
+        7.42E-01  0.00E+00 -9.21E+01 -1.39E+02  2.71E+02
 
 TH 6
+       -2.77E+02  0.00E+00  5.44E+00 -4.59E+01  5.29E+00  2.54E+02
 
 TH 7
+       -9.85E-03  0.00E+00 -2.86E-03 -1.78E-02  1.54E-02  5.32E-03  6.34E-06
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.36E+00  0.00E+00  4.78E-01 -6.76E+01  2.39E+01  9.51E+00  2.57E-02  0.00E+00  1.29E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.02E+02  0.00E+00 -7.61E+00 -5.32E+01  2.64E+01  1.43E+01  3.51E-03  0.00E+00  2.00E+01  0.00E+00  1.02E+02
 
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
 #CPUT: Total CPU Time in Seconds,       33.596
Stop Time:
Sat Sep 25 09:10:15 CDT 2021
