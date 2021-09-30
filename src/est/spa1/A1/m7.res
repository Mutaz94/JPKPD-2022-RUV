Wed Sep 29 21:48:11 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat7.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1543.67407631294        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8115E+02  4.4781E+01  4.3902E+01  1.2281E+02  1.0992E+02  4.7969E+00  1.1153E+01 -1.0189E+02  2.8203E+01 -1.1896E+01
            -9.7005E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1795.30458613902        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8467E-01  1.0281E+00  9.8241E-01  1.0038E+00  9.9460E-01  1.1269E+00  8.9917E-01  9.0372E-01  8.1089E-01  8.6172E-01
             1.5908E+00
 PARAMETER:  8.4548E-02  1.2774E-01  8.2253E-02  1.0379E-01  9.4582E-02  2.1949E-01 -6.2824E-03 -1.2411E-03 -1.0963E-01 -4.8826E-02
             5.6424E-01
 GRADIENT:   1.7339E+02  4.3820E+01 -3.2106E+00  9.2054E+01  7.0090E+01  5.6932E+01 -6.5632E+00 -1.3763E+01 -1.5608E+01 -1.5071E+01
            -2.0165E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1811.99192580667        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      185
 NPARAMETR:  9.9021E-01  8.9347E-01  4.3809E-01  1.0334E+00  6.0111E-01  1.1578E+00  9.9789E-01  1.9124E-01  8.6277E-01  5.5179E-01
             1.7984E+00
 PARAMETER:  9.0159E-02 -1.2645E-02 -7.2534E-01  1.3281E-01 -4.0898E-01  2.4651E-01  9.7889E-02 -1.5542E+00 -4.7602E-02 -4.9459E-01
             6.8690E-01
 GRADIENT:  -2.4241E+00 -1.1182E+01 -5.8119E+01  1.2477E+02  7.9972E+01  1.8299E+01 -8.8403E+00 -8.6544E-01  5.6830E+00 -5.0595E+00
            -5.6193E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1819.29725880235        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      363
 NPARAMETR:  9.9052E-01  9.8997E-01  4.8542E-01  9.3952E-01  6.6282E-01  1.0819E+00  1.0184E+00  1.8759E-01  8.1780E-01  5.7438E-01
             1.8840E+00
 PARAMETER:  9.0474E-02  8.9922E-02 -6.2274E-01  3.7617E-02 -3.1125E-01  1.7876E-01  1.1823E-01 -1.5735E+00 -1.0114E-01 -4.5446E-01
             7.3339E-01
 GRADIENT:  -2.3627E+00 -9.0585E+00  5.6771E+00 -1.5766E+01 -3.0449E+00 -5.3323E+00 -7.8469E-01 -6.7817E-01 -3.0849E+00 -2.1278E+00
            -1.5951E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1820.00807284869        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      538
 NPARAMETR:  9.9336E-01  1.1987E+00  4.4879E-01  8.3159E-01  7.3634E-01  1.0993E+00  8.8621E-01  1.9102E-01  9.1058E-01  5.9603E-01
             1.9140E+00
 PARAMETER:  9.3336E-02  2.8121E-01 -7.0120E-01 -8.4414E-02 -2.0607E-01  1.9469E-01 -2.0807E-02 -1.5554E+00  6.3303E-03 -4.1746E-01
             7.4917E-01
 GRADIENT:   5.1712E-01  4.8433E+00  1.3622E+00  2.6357E+00 -2.5598E+00  1.0204E-01  5.3025E-01 -5.1204E-01 -6.4228E-01 -8.5672E-02
            -3.3444E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1820.28274533329        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  9.9398E-01  1.4653E+00  4.0015E-01  6.8119E-01  8.4690E-01  1.1063E+00  7.7367E-01  3.6691E-01  1.0636E+00  5.9908E-01
             1.9139E+00
 PARAMETER:  9.3960E-02  4.8207E-01 -8.1592E-01 -2.8391E-01 -6.6172E-02  2.0107E-01 -1.5661E-01 -9.0263E-01  1.6166E-01 -4.1236E-01
             7.4914E-01
 GRADIENT:   5.3962E-01  2.5332E+01  6.3838E+00  6.2235E+00 -1.0028E+01  1.8420E+00  4.5467E-01 -1.7701E+00 -1.3413E+00 -1.1738E+00
            -7.6740E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1867.88953809103        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      905             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9745E-01  1.8177E+00  3.2345E-01  4.8085E-01  1.0915E+00  1.1122E+00  7.3931E-01  2.2974E+00  1.3019E+00  7.1950E-01
             1.9059E+00
 PARAMETER:  9.7450E-02  6.9756E-01 -1.0287E+00 -6.3219E-01  1.8760E-01  2.0633E-01 -2.0204E-01  9.3177E-01  3.6380E-01 -2.2920E-01
             7.4498E-01
 GRADIENT:   1.2598E+02  2.2267E+02 -1.1589E+01  7.5906E+01  6.4240E+01  3.1211E+01  2.3553E+01 -2.0308E+01  2.8799E+00  2.0900E+00
             1.6222E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1868.75599512614        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      999
 NPARAMETR:  9.7618E-01  1.7625E+00  3.2371E-01  4.8060E-01  1.0534E+00  1.0877E+00  7.1687E-01  2.2958E+00  1.3015E+00  6.4401E-01
             1.9068E+00
 PARAMETER:  7.5891E-02  6.6672E-01 -1.0279E+00 -6.3273E-01  1.5204E-01  1.8407E-01 -2.3287E-01  9.3109E-01  3.6352E-01 -3.4005E-01
             7.4540E-01
 GRADIENT:  -3.2457E+01 -7.9245E+01 -1.3146E+01  1.1556E+01  5.0334E+01 -6.8066E+00  1.4592E+01 -4.0177E+01 -6.1513E-03  1.0428E+00
             1.5651E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1872.16754832758        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1179
 NPARAMETR:  1.0045E+00  1.8043E+00  3.2423E-01  4.7963E-01  1.0376E+00  1.0681E+00  6.7098E-01  2.3013E+00  1.3010E+00  6.2719E-01
             1.8924E+00
 PARAMETER:  1.0445E-01  6.9015E-01 -1.0263E+00 -6.3474E-01  1.3693E-01  1.6584E-01 -2.9902E-01  9.3346E-01  3.6312E-01 -3.6651E-01
             7.3786E-01
 GRADIENT:   1.9525E+01  2.8513E-01 -5.0292E+00  2.2310E+01  1.4930E+01 -1.4138E+01  5.0800E+00 -3.9595E+01 -5.3481E+00  7.2244E-01
             1.4839E+02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1880.14221316238        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1369             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9392E-01  1.8212E+00  3.2861E-01  4.7239E-01  1.0416E+00  1.1055E+00  6.5821E-01  2.3330E+00  1.2963E+00  7.0766E-01
             1.8004E+00
 PARAMETER:  9.3906E-02  6.9948E-01 -1.0129E+00 -6.4995E-01  1.4075E-01  2.0031E-01 -3.1823E-01  9.4717E-01  3.5954E-01 -2.4579E-01
             6.8802E-01
 GRADIENT:   1.3548E+02  2.9161E+02  3.7442E+00  5.7986E+01  4.2393E+00  3.1380E+01  5.2986E+00 -1.3124E+01 -6.4298E+00  1.6229E+00
             1.3538E+02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1881.87011617914        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1549
 NPARAMETR:  9.8343E-01  1.8058E+00  3.3001E-01  4.6909E-01  1.0466E+00  1.1080E+00  6.4477E-01  2.3489E+00  1.2951E+00  7.9794E-01
             1.7749E+00
 PARAMETER:  8.3288E-02  6.9102E-01 -1.0086E+00 -6.5696E-01  1.4554E-01  2.0254E-01 -3.3886E-01  9.5396E-01  3.5859E-01 -1.2572E-01
             6.7374E-01
 GRADIENT:  -1.7022E+01 -1.3287E+01 -2.9534E+00  1.1554E+01  3.4003E+00  4.1111E-01 -3.0483E+00 -3.8714E+01 -1.0042E+01  6.6379E+00
             1.3008E+02

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1884.65245396047        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1687
 NPARAMETR:  9.8407E-01  1.7810E+00  3.3236E-01  4.6370E-01  1.0513E+00  1.1073E+00  6.4214E-01  2.3752E+00  1.2930E+00  8.0400E-01
             1.7337E+00
 PARAMETER:  8.2939E-02  6.7737E-01 -1.0018E+00 -6.6834E-01  1.5153E-01  2.0197E-01 -3.4304E-01  9.6535E-01  3.5709E-01 -1.1819E-01
             6.5008E-01
 GRADIENT:  -1.6061E+01  2.7622E+03 -9.5957E+02  2.8619E+03  1.6673E+01  2.0907E-01 -2.8036E+03  1.9181E+03  5.3667E+03 -1.6247E+04
            -2.8409E+03
 NUMSIGDIG:         0.7         2.3         2.3         2.3         0.7         2.6         2.3         2.3         2.3         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1687
 NO. OF SIG. DIGITS IN FINAL EST.:  0.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.3594E-03 -9.1228E-03 -1.8790E-02  3.0703E-02 -4.5041E-02
 SE:             2.9717E-02  2.3499E-02  2.0228E-02  2.2367E-02  1.7204E-02
 N:                     100         100         100         100         100

 P VAL.:         7.5280E-01  6.9785E-01  3.5293E-01  1.6984E-01  8.8447E-03

 ETASHRINKSD(%)  4.4274E-01  2.1275E+01  3.2232E+01  2.5068E+01  4.2364E+01
 ETASHRINKVR(%)  8.8351E-01  3.8024E+01  5.4076E+01  4.3852E+01  6.6781E+01
 EBVSHRINKSD(%)  7.8340E-01  2.1813E+01  4.8954E+01  3.2424E+01  3.9766E+01
 EBVSHRINKVR(%)  1.5607E+00  3.8868E+01  7.3943E+01  5.4335E+01  6.3719E+01
 RELATIVEINF(%)  9.8389E+01  7.6905E+00  8.6784E+00  5.7377E+00  9.7082E+00
 EPSSHRINKSD(%)  4.1240E+01
 EPSSHRINKVR(%)  6.5473E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1884.6524539604679     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -965.71392075579524     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1884.652       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  1.78E+00  3.32E-01  4.64E-01  1.05E+00  1.11E+00  6.42E-01  2.38E+00  1.29E+00  8.04E-01  1.73E+00
 


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
+        4.97E+06
 
 TH 2
+       -3.35E-01  3.34E+04
 
 TH 3
+       -4.85E+01  1.19E+03  4.30E+05
 
 TH 4
+        3.67E+01 -1.83E+01  3.63E+03  4.99E+05
 
 TH 5
+        3.56E+00 -1.84E+02 -4.41E+02  3.84E+02  5.56E+02
 
 TH 6
+        1.61E+01  7.08E+00 -3.77E+01  3.65E+01 -1.35E+06  9.60E+05
 
 TH 7
+        2.04E+00  1.43E+02 -4.13E+02  5.29E+02 -1.37E+06  9.75E+05  9.89E+05
 
 TH 8
+        7.85E+00 -7.54E+01  5.07E+02 -2.81E+02  2.88E+00  6.05E+00  6.39E+01  8.85E+03
 
 TH 9
+        3.80E+01 -4.27E+01  1.29E+02 -1.17E+02  8.64E+01  2.85E+01  3.66E+02 -1.87E+01  2.25E+05
 
 TH10
+       -3.13E-01  1.17E+01 -4.02E+01  5.52E+01 -4.56E+02 -1.38E+02 -1.55E+03  7.60E+00  4.02E+01  5.32E+06
 
 TH11
+       -2.37E+01  1.29E+02 -1.11E+03  5.91E+02 -4.55E+01 -1.04E+01 -1.23E+02 -6.14E+01  5.37E+01  6.82E+00  3.80E+04
 
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
 #CPUT: Total CPU Time in Seconds,       34.113
Stop Time:
Wed Sep 29 21:48:47 CDT 2021
