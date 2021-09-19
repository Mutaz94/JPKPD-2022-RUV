Sat Sep 18 14:49:30 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat67.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1652.04321951191        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2967E+01 -9.9839E+01  5.7275E+00 -1.6612E+02 -1.6968E+01 -2.5913E+00 -1.7867E+01  3.7763E+00 -1.9175E+01  5.1058E+00
             5.6347E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1663.39192072003        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.8254E-01  1.1052E+00  1.0473E+00  1.0530E+00  1.0816E+00  1.0140E+00  1.1439E+00  9.7180E-01  1.0016E+00  9.8840E-01
             9.7399E-01
 PARAMETER:  8.2384E-02  2.0001E-01  1.4624E-01  1.5169E-01  1.7846E-01  1.1393E-01  2.3442E-01  7.1390E-02  1.0156E-01  8.8330E-02
             7.3646E-02
 GRADIENT:   1.5837E+01  7.5489E+00  3.0212E+00  1.1725E+01  9.8679E+00  3.3246E+00 -4.0509E+00 -1.9690E+00  2.0183E+00 -8.9786E+00
            -8.7166E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1664.00417932335        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  9.8681E-01  1.0470E+00  1.0715E+00  1.0905E+00  1.0587E+00  1.0301E+00  1.2686E+00  9.2044E-01  9.6658E-01  1.0368E+00
             9.8950E-01
 PARAMETER:  8.6717E-02  1.4598E-01  1.6902E-01  1.8666E-01  1.5707E-01  1.2965E-01  3.3794E-01  1.7101E-02  6.6011E-02  1.3610E-01
             8.9440E-02
 GRADIENT:   2.5297E+01  1.0439E+01  7.9714E+00  1.1207E+01 -7.2511E+00  1.0479E+01  2.9794E+00 -1.5654E+00  3.4162E+00 -5.9263E-01
            -1.2248E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.25172133082        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      323
 NPARAMETR:  9.9687E-01  9.8201E-01  1.1013E+00  1.1332E+00  1.0509E+00  1.0173E+00  1.3084E+00  9.7603E-01  9.2797E-01  1.0349E+00
             9.9030E-01
 PARAMETER:  9.6869E-02  8.1848E-02  1.9646E-01  2.2505E-01  1.4963E-01  1.1711E-01  3.6879E-01  7.5735E-02  2.5240E-02  1.3429E-01
             9.0255E-02
 GRADIENT:   6.5464E+00  1.6140E+00 -7.4290E-01  1.5646E+00  6.9794E-01 -8.7969E-01 -3.2245E-01  2.4788E-01 -9.1841E-01 -1.0452E-01
            -1.0140E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1664.43712358322        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      498
 NPARAMETR:  9.9055E-01  7.4786E-01  1.2524E+00  1.2857E+00  1.0221E+00  1.0178E+00  1.5314E+00  1.0465E+00  8.7703E-01  1.0427E+00
             9.9133E-01
 PARAMETER:  9.0506E-02 -1.9054E-01  3.2509E-01  3.5132E-01  1.2189E-01  1.1767E-01  5.2621E-01  1.4545E-01 -3.1218E-02  1.4185E-01
             9.1297E-02
 GRADIENT:  -9.3677E-01  2.5153E+00  1.5562E+00  3.5898E+00 -1.9788E+00  4.9117E-01  2.2987E-01 -3.1736E-01 -3.6586E-01 -3.6919E-01
             7.7414E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1664.79175311200        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      677
 NPARAMETR:  9.7874E-01  3.9516E-01  1.5561E+00  1.5153E+00  1.0300E+00  1.0147E+00  1.8648E+00  1.3239E+00  8.3754E-01  1.0979E+00
             9.9170E-01
 PARAMETER:  7.8514E-02 -8.2847E-01  5.4217E-01  5.1559E-01  1.2956E-01  1.1457E-01  7.2316E-01  3.8060E-01 -7.7289E-02  1.9336E-01
             9.1666E-02
 GRADIENT:  -1.4631E+01  1.4908E+00  9.6046E-02  3.7611E+00 -9.5241E-01  1.2233E+00 -1.3795E-01  9.0239E-01  2.3703E+00  6.7249E-01
             5.7697E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1665.15818441289        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      857
 NPARAMETR:  9.8454E-01  2.3883E-01  1.6219E+00  1.6104E+00  1.0077E+00  1.0069E+00  2.1330E+00  1.3824E+00  8.0216E-01  1.0928E+00
             9.8959E-01
 PARAMETER:  8.4418E-02 -1.3320E+00  5.8357E-01  5.7651E-01  1.0764E-01  1.0683E-01  8.5752E-01  4.2384E-01 -1.2045E-01  1.8871E-01
             8.9534E-02
 GRADIENT:   3.5462E+00  4.0795E-01 -3.2583E-01  7.1848E-01 -9.9471E-02 -7.8587E-01 -6.4238E-01  1.6632E-01 -1.3862E+00  2.6696E-01
            -5.2581E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1665.45096760890        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1032
 NPARAMETR:  9.8017E-01  8.0339E-02  1.7438E+00  1.7164E+00  1.0010E+00  1.0061E+00  3.1365E+00  1.4960E+00  7.7228E-01  1.0922E+00
             9.9063E-01
 PARAMETER:  7.9971E-02 -2.4215E+00  6.5604E-01  6.4021E-01  1.0105E-01  1.0613E-01  1.2431E+00  5.0276E-01 -1.5840E-01  1.8819E-01
             9.0588E-02
 GRADIENT:  -2.9717E-01  5.0403E-01  1.7300E+00  1.1986E+01 -3.0388E+00 -1.3809E-01 -2.2655E-01 -3.5970E-01 -9.1095E-01 -6.8112E-02
            -3.2964E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1665.59626838025        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1207
 NPARAMETR:  9.7836E-01  2.1941E-02  1.7593E+00  1.7470E+00  9.9330E-01  1.0063E+00  5.8959E+00  1.5158E+00  7.6110E-01  1.0928E+00
             9.9117E-01
 PARAMETER:  7.8125E-02 -3.7194E+00  6.6491E-01  6.5788E-01  9.3282E-02  1.0626E-01  1.8743E+00  5.1593E-01 -1.7299E-01  1.8871E-01
             9.1127E-02
 GRADIENT:  -1.9123E+00  3.0934E-02  6.6332E-02  1.1020E+00 -7.7777E-01  2.8132E-01 -5.6737E-02  7.2648E-02 -3.6964E-01  2.4290E-01
             1.5815E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1665.62601929286        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1383
 NPARAMETR:  9.7868E-01  1.0000E-02  1.7758E+00  1.7541E+00  9.9639E-01  1.0058E+00  1.1217E+01  1.5295E+00  7.5899E-01  1.0929E+00
             9.9127E-01
 PARAMETER:  7.8453E-02 -4.9242E+00  6.7427E-01  6.6194E-01  9.6385E-02  1.0583E-01  2.5174E+00  5.2492E-01 -1.7576E-01  1.8888E-01
             9.1229E-02
 GRADIENT:  -7.5187E-01  0.0000E+00 -7.0545E-02 -8.8007E-01  4.5580E-01  1.7379E-01 -1.4233E-02 -6.8728E-02  1.3165E-01 -9.2928E-02
             9.8680E-02

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1665.62705578951        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1510
 NPARAMETR:  9.7902E-01  1.0000E-02  1.7741E+00  1.7544E+00  9.9542E-01  1.0054E+00  1.2075E+01  1.5288E+00  7.5857E-01  1.0926E+00
             9.9098E-01
 PARAMETER:  7.8795E-02 -5.0610E+00  6.7331E-01  6.6213E-01  9.5406E-02  1.0534E-01  2.5911E+00  5.2450E-01 -1.7632E-01  1.8856E-01
             9.0937E-02
 GRADIENT:  -6.3402E-02  0.0000E+00  1.8329E-02  9.0598E-02 -1.0766E-02 -3.6155E-04  7.1486E-04 -7.9792E-03 -3.0675E-02 -9.4515E-03
            -7.4194E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1510
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0635E-04  4.2702E-04 -3.6285E-02 -7.4107E-03 -5.0409E-02
 SE:             2.9843E-02  1.8084E-03  1.9188E-02  2.9239E-02  1.9579E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9448E-01  8.1333E-01  5.8617E-02  7.9992E-01  1.0036E-02

 ETASHRINKSD(%)  2.1627E-02  9.3942E+01  3.5718E+01  2.0472E+00  3.4406E+01
 ETASHRINKVR(%)  4.3248E-02  9.9633E+01  5.8678E+01  4.0524E+00  5.6975E+01
 EBVSHRINKSD(%)  4.0532E-01  9.4136E+01  3.8990E+01  2.3565E+00  3.1036E+01
 EBVSHRINKVR(%)  8.0899E-01  9.9656E+01  6.2778E+01  4.6575E+00  5.2440E+01
 RELATIVEINF(%)  9.3551E+01  9.2420E-03  1.0601E+01  3.0022E+00  7.2481E+00
 EPSSHRINKSD(%)  4.5359E+01
 EPSSHRINKVR(%)  7.0144E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1665.6270557895091     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -930.47622922577091     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.30
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1665.627       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.00E-02  1.77E+00  1.75E+00  9.95E-01  1.01E+00  1.21E+01  1.53E+00  7.59E-01  1.09E+00  9.91E-01
 


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
+        1.15E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.34E+00  0.00E+00  5.97E+01
 
 TH 4
+       -1.10E+01  0.00E+00 -1.41E+01  5.99E+02
 
 TH 5
+       -1.04E+01  0.00E+00 -1.23E+02 -3.99E+01  4.42E+02
 
 TH 6
+        2.16E+01  0.00E+00  1.90E+00 -1.45E+00 -3.61E+01  1.78E+02
 
 TH 7
+        1.78E-01  0.00E+00  4.52E-03 -4.56E-03  2.62E-02 -2.60E-02  1.16E-03
 
 TH 8
+       -2.19E-03  0.00E+00 -1.73E+01 -2.73E+00 -9.82E+00  3.24E-02  1.00E-02  2.11E+01
 
 TH 9
+        1.42E+01  0.00E+00  5.09E+00 -7.48E-01  9.03E+00  1.32E+00  7.99E-02 -2.10E+00  3.18E+02
 
 TH10
+        1.05E+01  0.00E+00  4.12E-01 -3.19E-01 -6.51E+01 -1.47E+01  1.30E-02  1.09E+01 -2.25E-01  5.41E+01
 
 TH11
+       -6.07E+00  0.00E+00 -8.33E+00 -1.20E+01  4.06E+00 -5.17E+00  2.83E-02  1.70E+01 -8.57E+00  1.19E+01  1.97E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.329
Stop Time:
Sat Sep 18 14:49:58 CDT 2021
