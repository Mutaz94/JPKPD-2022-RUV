Thu Sep 30 03:49:02 CDT 2021
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
$DATA ../../../../data/spa1/D/dat96.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   32040.5964079775        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0061E+03  8.3058E+02  4.4206E+01  7.9861E+02  2.1876E+01 -3.1648E+03 -1.4730E+03 -1.3641E+02 -2.0100E+03 -7.6550E+02
            -5.9904E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -357.591287396578        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0889E+00  8.7248E-01  8.8015E-01  1.2924E+00  1.5180E+00  2.2076E+00  1.2102E+00  9.6084E-01  1.2673E+00  9.3486E-01
             1.4173E+01
 PARAMETER:  1.8518E-01 -3.6419E-02 -2.7665E-02  3.5648E-01  5.1739E-01  8.9189E-01  2.9080E-01  6.0048E-02  3.3691E-01  3.2643E-02
             2.7513E+00
 GRADIENT:  -1.9479E+01  2.1104E+01 -3.0567E+00  5.7280E+00 -7.3658E+00  5.4181E+01 -5.6800E+00  4.5020E+00 -2.0010E+01  6.9997E-01
            -3.9562E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -406.968044173410        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2136E+00  7.7371E-01  2.2486E+00  2.1506E+00  3.1799E+01  1.6651E+00  1.4536E+00  2.2492E-01  4.4956E+00  3.7245E-01
             1.3734E+01
 PARAMETER:  2.9361E-01 -1.5655E-01  9.1032E-01  8.6575E-01  3.5594E+00  6.0990E-01  4.7406E-01 -1.3920E+00  1.6031E+00 -8.8765E-01
             2.7199E+00
 GRADIENT:   1.0005E+02 -1.6817E+01 -5.2090E+00  1.9821E+01 -8.3525E-02 -5.1033E+01  7.0106E+00  3.5967E-02  2.0935E+00  2.1236E-05
            -2.2582E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -418.336478337589        NO. OF FUNC. EVALS.:  86
 CUMULATIVE NO. OF FUNC. EVALS.:      241
 NPARAMETR:  1.1956E+00  9.9957E-01  2.5243E+00  1.8090E+00  2.6438E+01  1.6543E+00  1.1401E+00  1.0968E-01  4.6135E+00  2.5077E-01
             1.4149E+01
 PARAMETER:  2.7861E-01  9.9567E-02  1.0260E+00  6.9276E-01  3.3748E+00  6.0336E-01  2.3112E-01 -2.1102E+00  1.6290E+00 -1.2832E+00
             2.7496E+00
 GRADIENT:   8.1026E+01 -1.9451E+01 -8.7952E+00  8.7700E+00 -3.6363E-01 -4.4383E+01  7.6754E+00  1.5908E-02 -1.4177E+01  4.8078E-05
             1.3197E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -423.600157853175        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      367
 NPARAMETR:  1.1927E+00  1.0001E+00  2.5127E+00  1.8153E+00  2.5991E+01  1.6513E+00  1.1851E-01  1.0000E-02  4.6621E+00  4.7963E-02
             1.4025E+01
 PARAMETER:  2.7623E-01  1.0012E-01  1.0213E+00  6.9625E-01  3.3577E+00  6.0154E-01 -2.0328E+00 -5.1483E+00  1.6395E+00 -2.9373E+00
             2.7408E+00
 GRADIENT:   8.4162E+01 -2.7224E+01 -6.9038E+00  1.0267E+01 -2.9549E-01 -4.5514E+01  1.9710E-01  0.0000E+00 -1.6528E+01  1.9471E-05
             8.3472E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -423.633184869540        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      437
 NPARAMETR:  1.1908E+00  1.0002E+00  2.5144E+00  1.8152E+00  2.5986E+01  1.6543E+00  5.8385E-02  1.0000E-02  4.6748E+00  1.5485E-02
             1.3807E+01
 PARAMETER:  2.7464E-01  1.0020E-01  1.0220E+00  6.9622E-01  3.3576E+00  6.0341E-01 -2.7407E+00 -6.5985E+00  1.6422E+00 -4.0679E+00
             2.7252E+00
 GRADIENT:   8.7410E+01 -2.6718E+01 -6.8038E+00  1.1519E+01 -2.7565E-01 -4.6817E+01  5.0956E-02  0.0000E+00 -1.5602E+01  2.2230E-06
            -4.8878E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -425.256408656556        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      514
 NPARAMETR:  1.1509E+00  1.0018E+00  2.5561E+00  1.8118E+00  2.6020E+01  1.7236E+00  1.0000E-02  1.0000E-02  4.9309E+00  1.0000E-02
             1.2671E+01
 PARAMETER:  2.4056E-01  1.0178E-01  1.0385E+00  6.9431E-01  3.3589E+00  6.4444E-01 -1.3557E+01 -2.7229E+01  1.6955E+00 -2.0730E+01
             2.6393E+00
 GRADIENT:   8.4717E+01 -2.6862E+01 -6.6058E+00  1.7966E+01 -2.4292E-01 -4.0082E+01  0.0000E+00  0.0000E+00  4.5505E+00  0.0000E+00
            -7.1491E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -437.228132434352        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      652
 NPARAMETR:  1.0446E+00  1.0062E+00  2.6834E+00  1.7991E+00  2.6271E+01  1.9387E+00  1.0000E-02  1.0000E-02  5.6713E+00  1.0000E-02
             1.3501E+01
 PARAMETER:  1.4364E-01  1.0615E-01  1.0871E+00  6.8731E-01  3.3684E+00  7.6204E-01 -4.1004E+01 -7.4391E+01  1.8354E+00 -6.0885E+01
             2.7027E+00
 GRADIENT:   5.4561E+00 -4.8702E+01 -4.6346E+00  7.0445E+00 -4.1413E-01 -2.9290E+00  0.0000E+00  0.0000E+00  1.1009E+01  0.0000E+00
             1.6414E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -452.847673992374        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      747
 NPARAMETR:  1.0294E+00  1.4612E+00  3.8226E+00  1.2045E+00  8.3949E+01  1.8681E+00  1.0000E-02  1.0000E-02  5.6400E+00  1.0000E-02
             1.3506E+01
 PARAMETER:  1.2899E-01  4.7927E-01  1.4409E+00  2.8603E-01  4.5302E+00  7.2493E-01 -4.1004E+01 -7.4391E+01  1.8299E+00 -6.0885E+01
             2.7031E+00
 GRADIENT:   8.6321E-01  6.8922E+00 -2.6996E+00  5.3975E+00 -7.5519E-02  4.4746E+00  0.0000E+00  0.0000E+00 -3.3075E+00  0.0000E+00
            -1.7659E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -465.038887080806        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      927
 NPARAMETR:  1.0191E+00  1.8126E+00  5.2893E+00  7.7771E-01  3.4629E+02  1.7105E+00  1.0000E-02  1.0000E-02  7.7605E+00  1.0000E-02
             1.4545E+01
 PARAMETER:  1.1894E-01  6.9477E-01  1.7657E+00 -1.5140E-01  5.9473E+00  6.3677E-01 -4.1004E+01 -7.4391E+01  2.1490E+00 -6.0885E+01
             2.7772E+00
 GRADIENT:  -1.0800E+00  1.5788E+00 -6.6439E-04  1.7460E+00 -1.4511E-02  3.9502E+00  0.0000E+00  0.0000E+00  2.6198E+00  0.0000E+00
            -5.4851E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -465.121941041122        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:     1073
 NPARAMETR:  1.0182E+00  1.8220E+00  5.2630E+00  7.6097E-01  5.8966E+03  1.6999E+00  1.0000E-02  1.0000E-02  7.8038E+00  1.0000E-02
             1.4613E+01
 PARAMETER:  1.1808E-01  6.9996E-01  1.7607E+00 -1.7317E-01  8.7821E+00  6.3059E-01 -4.1004E+01 -7.4391E+01  2.1546E+00 -6.0885E+01
             2.7819E+00
 GRADIENT:   2.6120E+00  1.0502E+01  3.7826E-02  2.1277E+00 -8.6235E-04  1.0479E+01  0.0000E+00  0.0000E+00  8.1030E+01  0.0000E+00
             3.1216E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -465.154814994863        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:     1222
 NPARAMETR:  1.0182E+00  1.8214E+00  5.2646E+00  7.5374E-01  6.8042E+03  1.6930E+00  1.0000E-02  1.0000E-02  7.8091E+00  1.0000E-02
             1.4599E+01
 PARAMETER:  1.1808E-01  6.9963E-01  1.7610E+00 -1.8271E-01  8.9253E+00  6.2650E-01 -4.1004E+01 -7.4391E+01  2.1553E+00 -6.0885E+01
             2.7810E+00
 GRADIENT:  -1.2944E+00  2.3013E-01 -6.6206E-02  1.2459E+00 -7.8857E-04  2.8556E+00  0.0000E+00  0.0000E+00  3.1207E+00  0.0000E+00
             2.6155E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -465.182696546452        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     1419             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0180E+00  1.8248E+00  5.2429E+00  7.5270E-01  2.1866E+05  1.6944E+00  1.0000E-02  1.0000E-02  7.8510E+00  1.0000E-02
             1.4571E+01
 PARAMETER:  1.1785E-01  7.0144E-01  1.7569E+00 -1.8408E-01  1.2395E+01  6.2735E-01 -4.1004E+01 -7.4391E+01  2.1606E+00 -6.0885E+01
             2.7790E+00
 GRADIENT:   3.5232E+00  1.1369E+01 -1.1912E-02  2.2158E+00 -2.2838E-05  9.7250E+00  0.0000E+00  0.0000E+00  8.2737E+01  0.0000E+00
             2.7859E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -465.209463432973        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1586
 NPARAMETR:  1.0175E+00  1.8264E+00  5.2336E+00  7.4760E-01  5.3690E+06  1.6890E+00  1.0000E-02  1.0000E-02  7.8715E+00  1.0000E-02
             1.4553E+01
 PARAMETER:  1.1738E-01  7.0237E-01  1.7551E+00 -1.9089E-01  1.5596E+01  6.2416E-01 -4.1004E+01 -7.4391E+01  2.1632E+00 -6.0885E+01
             2.7778E+00
 GRADIENT:  -6.7313E-01  1.8760E+00 -1.6189E-01  1.4889E+00 -9.9880E-07  2.2990E+00  0.0000E+00  0.0000E+00  4.3489E+00  0.0000E+00
            -1.4955E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -465.226218277801        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1765
 NPARAMETR:  1.0174E+00  1.8302E+00  5.2155E+00  7.4311E-01  1.0544E+08  1.6918E+00  1.0000E-02  1.0000E-02  7.9112E+00  1.0000E-02
             1.4606E+01
 PARAMETER:  1.1723E-01  7.0442E-01  1.7516E+00 -1.9691E-01  1.8574E+01  6.2578E-01 -4.1004E+01 -7.4391E+01  2.1683E+00 -6.0885E+01
             2.7814E+00
 GRADIENT:   3.0179E+00  1.0412E+01 -9.2169E-02  2.1714E+00 -4.8064E-08  9.7753E+00  0.0000E+00  0.0000E+00  8.4781E+01  0.0000E+00
             3.0066E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -465.241438538629        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1936
 NPARAMETR:  1.0169E+00  1.8371E+00  5.1792E+00  7.3757E-01  8.8193E+07  1.6883E+00  1.0000E-02  1.0000E-02  8.0029E+00  1.0000E-02
             1.4728E+01
 PARAMETER:  1.1739E-01  7.0464E-01  1.7534E+00 -2.0543E-01  1.8581E+01  6.2597E-01 -4.1004E+01 -7.4391E+01  2.1689E+00 -6.0885E+01
             2.7807E+00
 GRADIENT:   6.4097E+03 -3.7053E+03  4.2924E+02 -3.6641E+03  3.8968E-06  8.5734E+02  0.0000E+00  0.0000E+00 -5.9609E+02  0.0000E+00
            -1.8537E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1936
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3139E-02 -1.0148E-03 -4.2427E-06  2.7231E-02 -2.6830E-13
 SE:             2.6233E-02  2.3111E-04  1.4560E-05  2.4650E-02  2.4707E-13
 N:                     100         100         100         100         100

 P VAL.:         2.0650E-01  1.1295E-05  7.7075E-01  2.6928E-01  2.7751E-01

 ETASHRINKSD(%)  1.2116E+01  9.9226E+01  9.9951E+01  1.7421E+01  1.0000E+02
 ETASHRINKVR(%)  2.2764E+01  9.9994E+01  1.0000E+02  3.1807E+01  1.0000E+02
 EBVSHRINKSD(%)  1.4086E+01  9.9601E+01  9.9904E+01  1.3645E+01  1.0000E+02
 EBVSHRINKVR(%)  2.6189E+01  9.9998E+01  1.0000E+02  2.5428E+01  1.0000E+02
 RELATIVEINF(%)  7.2324E+01  9.7145E-04  8.6465E-05  4.8292E+01  0.0000E+00
 EPSSHRINKSD(%)  2.5785E+00
 EPSSHRINKVR(%)  5.0906E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -465.24143853862898     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       453.69709466604371     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.00
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.85
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -465.241       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.83E+00  5.22E+00  7.37E-01  1.06E+08  1.69E+00  1.00E-02  1.00E-02  7.92E+00  1.00E-02  1.46E+01
 


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
+        4.57E+06
 
 TH 2
+        3.05E+02  3.94E+04
 
 TH 3
+       -1.21E+01 -8.57E+00  2.24E+02
 
 TH 4
+       -2.56E+06 -1.34E+02 -1.93E+03  2.85E+06
 
 TH 5
+        3.91E-12  1.80E-13 -3.80E-15  6.04E-12 -1.70E-22
 
 TH 6
+        4.01E+05 -3.07E+03 -3.39E+00 -1.01E+05  2.71E-13  2.07E+04
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.58E+01 -4.12E-01 -6.22E-01 -1.03E+01  6.63E-15  8.34E+02  0.00E+00  0.00E+00  3.78E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -9.67E+03 -2.92E+02 -1.26E+02  7.62E+03  3.74E-15 -1.18E+03  0.00E+00  0.00E+00 -2.08E+01  0.00E+00  4.50E+01
 
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
 #CPUT: Total CPU Time in Seconds,       50.925
Stop Time:
Thu Sep 30 03:49:54 CDT 2021
