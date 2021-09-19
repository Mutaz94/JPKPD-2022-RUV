Sat Sep 18 14:39:43 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat41.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1698.66669465770        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.7409E+01 -9.5835E+01 -2.9116E+01 -9.0294E+01  6.7767E+01  1.2610E+01  5.4459E+00  2.0457E+00  2.5874E+01 -4.9306E+00
             6.9663E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1707.29369880358        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      133
 NPARAMETR:  1.0388E+00  1.0689E+00  9.9173E-01  1.0286E+00  9.6875E-01  9.2845E-01  9.2000E-01  1.0014E+00  8.0240E-01  9.9707E-01
             1.0072E+00
 PARAMETER:  1.3811E-01  1.6662E-01  9.1698E-02  1.2819E-01  6.8254E-02  2.5761E-02  1.6619E-02  1.0142E-01 -1.2015E-01  9.7069E-02
             1.0716E-01
 GRADIENT:   1.3814E+01 -3.7566E+00  1.7982E+00 -1.2794E+01  5.7828E+00 -1.9203E+01 -5.1064E+00 -1.1913E+00 -1.5328E+01 -1.0800E+00
            -1.7129E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1707.75405385011        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      311
 NPARAMETR:  1.0376E+00  1.0944E+00  8.4078E-01  1.0069E+00  8.9323E-01  9.5542E-01  9.6573E-01  9.2542E-01  7.9713E-01  8.6909E-01
             1.0118E+00
 PARAMETER:  1.3689E-01  1.9022E-01 -7.3431E-02  1.0691E-01 -1.2909E-02  5.4398E-02  6.5131E-02  2.2491E-02 -1.2674E-01 -4.0312E-02
             1.1172E-01
 GRADIENT:   6.6822E+00  3.6807E+00  2.9982E+00 -4.8307E+00 -8.0706E-01 -7.7383E+00 -2.7032E+00  8.8796E-01 -1.2202E+01 -5.1469E+00
             6.4598E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1708.59884245766        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      489
 NPARAMETR:  1.0332E+00  9.9772E-01  9.1150E-01  1.0722E+00  8.9109E-01  9.7374E-01  1.0009E+00  8.4857E-01  8.2965E-01  9.3959E-01
             1.0046E+00
 PARAMETER:  1.3270E-01  9.7716E-02  7.3369E-03  1.6974E-01 -1.5313E-02  7.3388E-02  1.0090E-01 -6.4201E-02 -8.6750E-02  3.7688E-02
             1.0462E-01
 GRADIENT:  -9.4041E-01  2.6183E+00  2.3081E+00  4.2509E-01 -2.3784E+00  1.9499E-01 -3.5692E-03 -6.1228E-01 -4.1545E-01 -2.0872E-01
             1.5100E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1708.68949058071        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  1.0307E+00  8.2023E-01  1.0857E+00  1.1864E+00  9.0541E-01  9.6944E-01  1.0951E+00  9.7977E-01  7.8475E-01  9.8752E-01
             1.0029E+00
 PARAMETER:  1.3028E-01 -9.8175E-02  1.8224E-01  2.7095E-01  6.3026E-04  6.8958E-02  1.9084E-01  7.9565E-02 -1.4239E-01  8.7444E-02
             1.0291E-01
 GRADIENT:  -6.6834E-01 -2.8561E+00  8.3528E-01 -8.6212E+00 -5.9007E-01 -5.4578E-01 -3.8880E-01 -2.4610E-01 -3.6407E-01  8.3461E-01
            -2.7072E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1708.81709026662        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      844
 NPARAMETR:  1.0290E+00  6.6028E-01  1.1975E+00  1.2991E+00  8.9549E-01  9.6928E-01  1.2697E+00  1.0448E+00  7.4025E-01  9.8942E-01
             1.0034E+00
 PARAMETER:  1.2859E-01 -3.1509E-01  2.8024E-01  3.6165E-01 -1.0385E-02  6.8799E-02  3.3875E-01  1.4386E-01 -2.0076E-01  8.9366E-02
             1.0343E-01
 GRADIENT:   4.3943E-01  5.1276E+00  1.0976E+00  1.3278E+01 -1.2465E+00  2.3524E-01 -1.2953E-01 -7.3429E-01 -2.1576E-01 -9.2438E-01
            -3.5713E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1709.07148535526        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1021
 NPARAMETR:  1.0253E+00  4.3539E-01  1.4223E+00  1.4392E+00  9.1295E-01  9.6658E-01  1.6176E+00  1.2464E+00  6.9949E-01  1.0309E+00
             1.0048E+00
 PARAMETER:  1.2500E-01 -7.3152E-01  4.5225E-01  4.6408E-01  8.9271E-03  6.6012E-02  5.8094E-01  3.2027E-01 -2.5740E-01  1.3047E-01
             1.0482E-01
 GRADIENT:   1.2142E+00  1.0304E-01  1.1923E-01 -3.2408E+00 -1.9357E-01  6.6930E-01  7.0788E-01  6.5902E-01  2.8265E+00  1.7643E-02
             3.8747E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1709.19665984272        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1196
 NPARAMETR:  1.0222E+00  3.0769E-01  1.5301E+00  1.5261E+00  9.1196E-01  9.6275E-01  1.8599E+00  1.3429E+00  6.7328E-01  1.0428E+00
             1.0048E+00
 PARAMETER:  1.2197E-01 -1.0786E+00  5.2534E-01  5.2272E-01  7.8423E-03  6.2038E-02  7.2052E-01  3.9486E-01 -2.9559E-01  1.4188E-01
             1.0479E-01
 GRADIENT:  -1.0470E+00  1.5548E+00  2.7397E-02  1.0430E+01 -1.8475E+00 -7.9391E-02 -1.4839E-01  3.2561E-01  3.2470E-01  5.6421E-02
             1.9529E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1709.27335214339        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1371
 NPARAMETR:  1.0208E+00  2.0197E-01  1.6803E+00  1.5936E+00  9.2937E-01  9.6124E-01  2.2020E+00  1.4763E+00  6.5528E-01  1.0586E+00
             1.0054E+00
 PARAMETER:  1.2061E-01 -1.4996E+00  6.1898E-01  5.6597E-01  2.6749E-02  6.0473E-02  8.8938E-01  4.8956E-01 -3.2270E-01  1.5696E-01
             1.0535E-01
 GRADIENT:   4.3754E-02  2.0183E-01  5.2112E-01  1.5828E+00 -2.2129E-01  1.9352E-02 -2.2427E-01 -4.2510E-01 -1.4519E-01 -3.5278E-01
            -1.5934E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1709.38031322571        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1546
 NPARAMETR:  1.0193E+00  1.0332E-01  1.8470E+00  1.6598E+00  9.4793E-01  9.5923E-01  3.3122E+00  1.6308E+00  6.3671E-01  1.0756E+00
             1.0064E+00
 PARAMETER:  1.1908E-01 -2.1700E+00  7.1358E-01  6.0672E-01  4.6522E-02  5.8375E-02  1.2976E+00  5.8908E-01 -3.5143E-01  1.7289E-01
             1.0637E-01
 GRADIENT:   1.4322E-01  1.4139E-01  5.2021E-01 -3.5745E-01 -8.1595E-01 -1.8258E-01  1.1041E-02 -1.2906E-01  1.2847E-01  1.9683E-01
             6.3318E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1709.42013040407        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1722
 NPARAMETR:  1.0182E+00  5.0962E-02  1.9735E+00  1.7006E+00  9.6655E-01  9.5902E-01  4.9205E+00  1.7443E+00  6.2995E-01  1.0891E+00
             1.0072E+00
 PARAMETER:  1.1805E-01 -2.8767E+00  7.7980E-01  6.3100E-01  6.5979E-02  5.8153E-02  1.6934E+00  6.5637E-01 -3.6211E-01  1.8533E-01
             1.0713E-01
 GRADIENT:  -6.6793E-01  3.7642E-01 -4.8137E-01  1.0524E+01 -2.3335E-02  6.3869E-04  2.4217E-01  2.3023E-01  1.1399E+00  2.7810E-01
            -3.1512E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1709.46942765764        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1898
 NPARAMETR:  1.0180E+00  2.0965E-02  1.9941E+00  1.7165E+00  9.6196E-01  9.5843E-01  7.1797E+00  1.7607E+00  6.2334E-01  1.0843E+00
             1.0072E+00
 PARAMETER:  1.1781E-01 -3.7649E+00  7.9019E-01  6.4029E-01  6.1218E-02  5.7538E-02  2.0713E+00  6.6569E-01 -3.7267E-01  1.8097E-01
             1.0713E-01
 GRADIENT:  -9.2444E-02  3.7752E-02  4.2002E-01  2.1004E+00 -1.3156E+00 -5.1962E-02 -2.0430E-02  2.7543E-03 -5.7629E-01  8.6472E-02
            -3.4063E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1709.48146129107        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2073
 NPARAMETR:  1.0178E+00  1.0000E-02  2.0190E+00  1.7234E+00  9.6615E-01  9.5838E-01  9.7297E+00  1.7807E+00  6.2277E-01  1.0868E+00
             1.0073E+00
 PARAMETER:  1.1769E-01 -4.5304E+00  8.0259E-01  6.4431E-01  6.5564E-02  5.7492E-02  2.3752E+00  6.7702E-01 -3.7357E-01  1.8324E-01
             1.0731E-01
 GRADIENT:   1.1338E-02  0.0000E+00  4.6891E-02 -5.1711E-01  7.5435E-02 -2.1069E-03 -1.1988E-02 -5.8744E-02 -1.6911E-02 -4.6319E-02
             1.9465E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1709.48234755128        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2250
 NPARAMETR:  1.0178E+00  1.0000E-02  2.0183E+00  1.7236E+00  9.6603E-01  9.5844E-01  1.0700E+01  1.7810E+00  6.2262E-01  1.0870E+00
             1.0073E+00
 PARAMETER:  1.1766E-01 -4.7649E+00  8.0225E-01  6.4440E-01  6.5436E-02  5.7547E-02  2.4702E+00  6.7717E-01 -3.7382E-01  1.8340E-01
             1.0728E-01
 GRADIENT:  -5.1837E-02  0.0000E+00 -2.2405E-02  3.9281E-02  1.7337E-02 -4.8392E-03 -6.6093E-04  4.1750E-03  6.2329E-03  5.1222E-03
             6.6994E-04

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1709.48234755128        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     2272
 NPARAMETR:  1.0178E+00  1.0000E-02  2.0183E+00  1.7236E+00  9.6603E-01  9.5844E-01  1.0700E+01  1.7810E+00  6.2262E-01  1.0870E+00
             1.0073E+00
 PARAMETER:  1.1766E-01 -4.7649E+00  8.0225E-01  6.4440E-01  6.5436E-02  5.7547E-02  2.4702E+00  6.7717E-01 -3.7382E-01  1.8340E-01
             1.0728E-01
 GRADIENT:  -5.1837E-02  0.0000E+00 -2.2405E-02  3.9281E-02  1.7337E-02 -4.8392E-03 -6.6093E-04  4.1750E-03  6.2329E-03  5.1222E-03
             6.6994E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2272
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.5632E-05  4.4944E-04 -3.8243E-02 -9.5614E-03 -4.8490E-02
 SE:             2.9834E-02  1.9679E-03  1.9068E-02  2.8941E-02  1.9980E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9878E-01  8.1935E-01  4.4901E-02  7.4111E-01  1.5229E-02

 ETASHRINKSD(%)  5.1917E-02  9.3407E+01  3.6119E+01  3.0453E+00  3.3063E+01
 ETASHRINKVR(%)  1.0381E-01  9.9565E+01  5.9192E+01  5.9978E+00  5.5195E+01
 EBVSHRINKSD(%)  4.4399E-01  9.3725E+01  3.9678E+01  3.5716E+00  2.8984E+01
 EBVSHRINKVR(%)  8.8601E-01  9.9606E+01  6.3613E+01  7.0156E+00  4.9567E+01
 RELATIVEINF(%)  9.5304E+01  1.0679E-02  1.0399E+01  2.8292E+00  9.0783E+00
 EPSSHRINKSD(%)  4.5197E+01
 EPSSHRINKVR(%)  6.9966E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1709.4823475512806     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -974.33152098754238     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.86
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1709.482       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.00E-02  2.02E+00  1.72E+00  9.66E-01  9.58E-01  1.07E+01  1.78E+00  6.23E-01  1.09E+00  1.01E+00
 


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
+        1.17E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -3.03E+00  0.00E+00  3.36E+01
 
 TH 4
+       -1.63E+01  0.00E+00 -1.85E+01  8.99E+02
 
 TH 5
+       -2.46E+01  0.00E+00 -9.51E+01 -7.38E+01  5.49E+02
 
 TH 6
+        9.82E+00  0.00E+00 -1.93E+00 -3.50E+00  8.43E+00  2.09E+02
 
 TH 7
+        2.89E-02  0.00E+00  3.41E-03  1.21E-02  7.17E-02  1.43E-01  2.30E-03
 
 TH 8
+       -1.51E+00  0.00E+00 -1.30E+01 -3.60E+00 -1.13E+01  1.40E+00  3.98E-03  1.91E+01
 
 TH 9
+       -3.16E+00  0.00E+00  6.25E+00 -3.47E-01 -3.04E+00  9.28E+00  1.26E-01 -8.13E-01  4.52E+02
 
 TH10
+        7.93E+00  0.00E+00  6.10E-01 -1.74E+00 -8.16E+01 -6.23E+00  1.35E-02  1.05E+01  2.13E+00  5.41E+01
 
 TH11
+       -2.43E+00  0.00E+00 -4.10E+00 -1.36E+01 -1.33E+01 -4.37E+00  1.98E-01  3.98E+00  1.09E+01  2.22E+01  2.13E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.332
Stop Time:
Sat Sep 18 14:40:20 CDT 2021
