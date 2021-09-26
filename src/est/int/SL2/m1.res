Sat Sep 25 00:46:41 CDT 2021
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
$DATA ../../../../data/int/SL2/dat1.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      998
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      898
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2189.49043776466        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.1143E+01  6.3138E+01  1.3732E+02  9.8404E+00  2.2430E+01  8.7957E+00 -8.4416E+01 -1.7042E+02 -6.8870E+01 -1.8260E+01
            -3.0096E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3052.67932823655        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0240E+00  1.1136E+00  9.6523E-01  9.4294E-01  1.0821E+00  9.7280E-01  1.0500E+00  1.0158E+00  9.8401E-01  1.0352E+00
             1.9642E+00
 PARAMETER:  1.2370E-01  2.0757E-01  6.4609E-02  4.1245E-02  1.7887E-01  7.2419E-02  1.4876E-01  1.1563E-01  8.3882E-02  1.3455E-01
             7.7507E-01
 GRADIENT:   3.6467E+01 -1.1100E+01 -2.9382E+00 -1.8476E+01 -9.8717E+00 -1.9944E+00  4.1662E+00 -7.0640E-01 -7.0904E+00 -6.0468E+00
            -8.2670E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3057.13994140665        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0305E+00  1.4873E+00  1.1103E+00  7.9467E-01  1.4541E+00  1.0094E+00  7.9997E-01  8.6609E-01  1.2062E+00  1.3264E+00
             2.0021E+00
 PARAMETER:  1.3004E-01  4.9697E-01  2.0467E-01 -1.2983E-01  4.7436E-01  1.0933E-01 -1.2318E-01 -4.3763E-02  2.8749E-01  3.8243E-01
             7.9421E-01
 GRADIENT:   4.6593E+01  6.8498E+01  1.0438E+00  8.0329E+01  2.5575E+01  1.1549E+01  1.9546E+00 -3.4761E+00  4.1981E+00  4.0750E+00
            -6.4233E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3061.54268580142        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0087E+00  1.4401E+00  1.0822E+00  7.6033E-01  1.3929E+00  9.7466E-01  7.7950E-01  1.6209E+00  1.1626E+00  1.2762E+00
             2.0405E+00
 PARAMETER:  1.0864E-01  4.6472E-01  1.7899E-01 -1.7401E-01  4.3141E-01  7.4336E-02 -1.4910E-01  5.8299E-01  2.5069E-01  3.4389E-01
             8.1322E-01
 GRADIENT:  -3.1663E+00 -1.1251E+01 -3.1683E+00 -1.3474E+00 -2.4525E+00 -8.6744E-01 -9.1771E-02  3.5752E+00  1.8739E+00 -2.0699E+00
             1.4230E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3064.74101981162        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0174E+00  1.6907E+00  1.0822E+00  6.0885E-01  1.6245E+00  9.8354E-01  7.1323E-01  1.9804E+00  1.3399E+00  1.4002E+00
             2.0330E+00
 PARAMETER:  1.1727E-01  6.2511E-01  1.7903E-01 -3.9619E-01  5.8519E-01  8.3402E-02 -2.3795E-01  7.8331E-01  3.9263E-01  4.3664E-01
             8.0950E-01
 GRADIENT:   4.2755E+00 -1.0289E+01 -3.9583E-01  2.4677E+00  1.3406E+00  1.3669E+00  8.0434E-01  6.4838E-01  6.2321E-01 -4.0116E-01
             4.5245E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3065.63839592087        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      544
 NPARAMETR:  1.0157E+00  1.9232E+00  9.1947E-01  4.5649E-01  1.7951E+00  9.8134E-01  6.6804E-01  2.1632E+00  1.6101E+00  1.4918E+00
             2.0324E+00
 PARAMETER:  1.1560E-01  7.5400E-01  1.6044E-02 -6.8419E-01  6.8508E-01  8.1164E-02 -3.0341E-01  8.7158E-01  5.7631E-01  4.9999E-01
             8.0920E-01
 GRADIENT:   3.0973E-01  2.3725E+00  1.8880E+00 -3.7097E-01 -1.7232E+00  3.9170E-01 -6.2999E-01  1.6737E-01  2.1517E-01  9.1700E-04
            -9.1308E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3066.65645409621        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  1.0133E+00  2.1767E+00  4.9437E-01  2.9284E-01  1.9724E+00  9.7674E-01  6.5220E-01  1.1305E+00  2.1034E+00  1.5792E+00
             2.0462E+00
 PARAMETER:  1.1324E-01  8.7782E-01 -6.0448E-01 -1.1281E+00  7.7926E-01  7.6461E-02 -3.2740E-01  2.2265E-01  8.4357E-01  5.5694E-01
             8.1600E-01
 GRADIENT:  -5.8089E+00  1.6381E+01 -3.5732E+00  8.9096E+00  1.2018E+01 -1.7682E+00  6.8020E-01  7.1435E-01 -3.0499E+00  9.0151E-02
             3.3249E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3067.95685447421        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      897
 NPARAMETR:  1.0158E+00  2.2841E+00  4.0259E-01  2.0563E-01  2.0411E+00  9.8145E-01  6.3086E-01  2.9431E-01  2.7029E+00  1.6146E+00
             2.0476E+00
 PARAMETER:  1.1571E-01  9.2596E-01 -8.0983E-01 -1.4817E+00  8.1349E-01  8.1271E-02 -3.6068E-01 -1.1231E+00  1.0943E+00  5.7909E-01
             8.1667E-01
 GRADIENT:   5.8161E-01  8.1479E-01  1.1603E-01 -6.8896E-02 -1.0657E+00  2.0087E-01  1.0805E-01  8.2893E-02  3.0063E-01  7.6818E-02
             1.2273E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3067.99400985271        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  1.0156E+00  2.2740E+00  4.1201E-01  2.1250E-01  2.0370E+00  9.8085E-01  6.3163E-01  1.2435E-01  2.6506E+00  1.6118E+00
             2.0470E+00
 PARAMETER:  1.1548E-01  9.2152E-01 -7.8671E-01 -1.4488E+00  8.1146E-01  8.0659E-02 -3.5946E-01 -1.9846E+00  1.0748E+00  5.7735E-01
             8.1639E-01
 GRADIENT:   7.4362E-02  3.2290E-02 -6.2082E-02  2.5726E-01  5.4831E-01 -3.2304E-02  3.5824E-02  1.6233E-02  1.7952E-02  1.6775E-01
             6.7173E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3068.00292756170        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1251
 NPARAMETR:  1.0156E+00  2.2770E+00  4.0987E-01  2.1018E-01  2.0383E+00  9.8093E-01  6.3106E-01  1.1801E-02  2.6670E+00  1.6120E+00
             2.0470E+00
 PARAMETER:  1.1545E-01  9.2285E-01 -7.9192E-01 -1.4598E+00  8.1210E-01  8.0742E-02 -3.6036E-01 -4.3396E+00  1.0810E+00  5.7746E-01
             8.1636E-01
 GRADIENT:  -7.5383E-03  5.6656E-02  1.7027E-02 -7.5803E-03 -2.3950E-02  5.3069E-06 -2.1274E-02  1.4587E-04 -3.8876E-03  8.3719E-03
            -3.1112E-02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -3068.00295382456        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     1314
 NPARAMETR:  1.0155E+00  2.2767E+00  4.0999E-01  2.1014E-01  2.0381E+00  9.8089E-01  6.3097E-01  1.0000E-02  2.6660E+00  1.6120E+00
             2.0469E+00
 PARAMETER:  1.1546E-01  9.2283E-01 -7.9242E-01 -1.4597E+00  8.1210E-01  8.0747E-02 -3.6020E-01 -4.8449E+00  1.0809E+00  5.7742E-01
             8.1638E-01
 GRADIENT:   1.6890E-02  6.4499E-02 -4.9174E-03  6.9929E-03  8.0100E-03  2.6512E-03  1.1177E-02  0.0000E+00  7.4083E-03 -2.2552E-03
             2.9278E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1314
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0659E-03 -3.7078E-02 -3.8549E-05  4.3969E-02 -2.4233E-02
 SE:             2.9638E-02  2.4702E-02  5.1275E-05  1.9696E-02  2.6715E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7131E-01  1.3336E-01  4.5217E-01  2.5594E-02  3.6437E-01

 ETASHRINKSD(%)  7.0918E-01  1.7245E+01  9.9828E+01  3.4015E+01  1.0501E+01
 ETASHRINKVR(%)  1.4133E+00  3.1517E+01  1.0000E+02  5.6460E+01  1.9899E+01
 EBVSHRINKSD(%)  1.0197E+00  1.5921E+01  9.9842E+01  3.9925E+01  7.8737E+00
 EBVSHRINKVR(%)  2.0290E+00  2.9308E+01  1.0000E+02  6.3910E+01  1.5127E+01
 RELATIVEINF(%)  9.7944E+01  1.8264E+01  1.5145E-04  9.0016E+00  5.2827E+01
 EPSSHRINKSD(%)  1.7386E+01
 EPSSHRINKVR(%)  3.1749E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          898
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1650.4136056355921     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3068.0029538245612     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1417.5893481889691     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    32.88
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.83
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3068.003       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.28E+00  4.10E-01  2.10E-01  2.04E+00  9.81E-01  6.31E-01  1.00E-02  2.67E+00  1.61E+00  2.05E+00
 


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
+        1.11E+03
 
 TH 2
+       -9.51E+00  4.00E+02
 
 TH 3
+        1.18E+00  4.53E+01  1.08E+02
 
 TH 4
+       -1.26E+01  4.17E+02 -1.38E+02  1.42E+03
 
 TH 5
+       -1.78E+00 -3.19E+01 -1.92E+01  8.46E+01  8.76E+01
 
 TH 6
+       -4.28E+00 -2.30E+00  1.25E+00 -5.86E+00 -6.01E-01  2.03E+02
 
 TH 7
+       -3.35E-01 -9.12E+00  8.00E+00 -2.81E+01 -3.30E+00 -3.71E-01  2.89E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.98E-01 -7.03E+00 -4.37E+00  5.91E+01 -2.51E-01 -1.48E-01  5.53E+00  0.00E+00  9.18E+00
 
 TH10
+        1.27E-01 -6.39E+00  4.50E+00  1.30E+01 -6.45E+00  1.22E+00  3.71E+00  0.00E+00  7.18E-01  5.44E+01
 
 TH11
+       -1.24E+01 -1.95E+01 -1.14E+01 -2.49E+00 -6.16E-01  2.36E+00  8.53E+00  0.00E+00  2.37E+00  3.40E+00  2.78E+02
 
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
 #CPUT: Total CPU Time in Seconds,       45.807
Stop Time:
Sat Sep 25 00:47:28 CDT 2021
