Sat Sep 18 10:08:00 CDT 2021
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
$DATA ../../../../data/spa/A2/dat90.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m90.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1158.87181963606        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1662E+02 -3.7195E+00  2.9854E+01 -4.1235E+01  4.9910E+01 -1.3735E+00 -2.7424E+01 -6.0977E+00 -4.6072E+01 -5.4899E+01
            -8.3077E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1414.68834354424        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9205E-01  9.8470E-01  1.0112E+00  1.0524E+00  9.6795E-01  9.6249E-01  1.0702E+00  9.4475E-01  1.0820E+00  1.0482E+00
             2.1338E+00
 PARAMETER:  9.2020E-02  8.4580E-02  1.1118E-01  1.5112E-01  6.7429E-02  6.1772E-02  1.6783E-01  4.3163E-02  1.7880E-01  1.4709E-01
             8.5789E-01
 GRADIENT:   4.5730E+01  4.9913E+00  2.7309E-01  1.6349E+00 -2.1633E+00 -7.2551E+00 -2.4847E+00  5.1494E+00  1.0653E+00 -2.2847E+00
            -8.7350E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1417.11386698679        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.8738E-01  7.5110E-01  9.0648E-01  1.2272E+00  8.2547E-01  9.8702E-01  1.3079E+00  4.3894E-01  9.4413E-01  1.0417E+00
             2.1041E+00
 PARAMETER:  8.7296E-02 -1.8621E-01  1.8112E-03  3.0474E-01 -9.1806E-02  8.6937E-02  3.6842E-01 -7.2340E-01  4.2508E-02  1.4086E-01
             8.4388E-01
 GRADIENT:   3.2069E+01  2.3579E+01 -8.1489E+00  5.9689E+01  4.4355E+00  1.4729E+00 -2.7150E+00  1.2418E+00 -9.5490E+00  5.5779E+00
            -1.2140E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1420.07943083945        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  9.7033E-01  5.7489E-01  6.4397E-01  1.2626E+00  5.9761E-01  9.8339E-01  1.6359E+00  1.4727E-01  9.1311E-01  7.3042E-01
             2.1472E+00
 PARAMETER:  6.9881E-02 -4.5358E-01 -3.4010E-01  3.3317E-01 -4.1481E-01  8.3252E-02  5.9219E-01 -1.8155E+00  9.1039E-03 -2.1414E-01
             8.6415E-01
 GRADIENT:  -1.0387E+01  1.1641E+01  1.1616E+00  2.8007E+01 -4.7632E+00 -1.2367E-01 -1.5839E-01  2.5270E-01 -1.9206E+00 -2.2312E+00
             7.2291E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1420.88992554175        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  9.7232E-01  4.0324E-01  6.2061E-01  1.3270E+00  5.4096E-01  9.8243E-01  1.9751E+00  4.6976E-02  8.9674E-01  7.6553E-01
             2.1278E+00
 PARAMETER:  7.1930E-02 -8.0823E-01 -3.7706E-01  3.8292E-01 -5.1441E-01  8.2270E-02  7.8061E-01 -2.9581E+00 -8.9871E-03 -1.6718E-01
             8.5508E-01
 GRADIENT:   2.9128E-01 -1.1904E+00 -6.1214E+00  9.3333E-01  5.8428E+00  3.5679E-01 -6.4619E-01  3.0936E-02  1.8049E+00  1.2586E+00
             4.3559E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1421.53812076740        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      405
 NPARAMETR:  9.6886E-01  3.5068E-01  7.5603E-01  1.3825E+00  6.1077E-01  9.7681E-01  2.0365E+00  4.1041E-02  8.8420E-01  8.8286E-01
             2.1389E+00
 PARAMETER:  6.8361E-02 -9.4788E-01 -1.7968E-01  4.2388E-01 -3.9303E-01  7.6540E-02  8.1122E-01 -3.0932E+00 -2.3067E-02 -2.4587E-02
             8.6029E-01
 GRADIENT:  -8.3025E+00 -1.8729E+00 -3.2380E+00 -1.4721E+01  3.4363E+00 -9.1478E-01 -6.9310E-01  1.2778E-02  8.4977E-01  1.0263E+00
            -3.6536E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1422.24472871086        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      580
 NPARAMETR:  9.6618E-01  1.6719E-01  8.1227E-01  1.4985E+00  5.9968E-01  9.7349E-01  2.7984E+00  1.0000E-02  8.3666E-01  9.2255E-01
             2.1532E+00
 PARAMETER:  6.5594E-02 -1.6886E+00 -1.0792E-01  5.0450E-01 -4.1135E-01  7.3132E-02  1.1291E+00 -5.5621E+00 -7.8335E-02  1.9382E-02
             8.6696E-01
 GRADIENT:  -3.5468E+00  4.4257E-01  3.6590E+00 -3.6183E-01 -5.1099E+00 -5.2924E-01 -8.7769E-01  0.0000E+00 -2.6015E+00 -1.1598E-01
            -7.4338E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1422.50652075233        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      755
 NPARAMETR:  9.6524E-01  7.3608E-02  8.2420E-01  1.5491E+00  5.9088E-01  9.7167E-01  4.3191E+00  1.0000E-02  8.2281E-01  9.3494E-01
             2.1581E+00
 PARAMETER:  6.4622E-02 -2.5090E+00 -9.3339E-02  5.3767E-01 -4.2614E-01  7.1259E-02  1.5630E+00 -8.7785E+00 -9.5034E-02  3.2725E-02
             8.6924E-01
 GRADIENT:   5.8311E-01  1.4261E-02 -1.3498E+00 -2.0780E+00  2.4128E+00 -4.6544E-01 -1.0284E-01  0.0000E+00 -1.3841E+00 -1.3877E-01
             3.3452E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1422.61583285712        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      931
 NPARAMETR:  9.6342E-01  1.6018E-02  8.2660E-01  1.5789E+00  5.8101E-01  9.7165E-01  8.9319E+00  1.0000E-02  8.1774E-01  9.4006E-01
             2.1569E+00
 PARAMETER:  6.2731E-02 -4.0341E+00 -9.0439E-02  5.5673E-01 -4.4299E-01  7.1238E-02  2.2896E+00 -1.5138E+01 -1.0121E-01  3.8188E-02
             8.6865E-01
 GRADIENT:   3.8634E-01 -3.8042E-02 -5.9962E-01 -2.4958E+00  1.4036E+00 -5.5469E-02 -9.1846E-02  0.0000E+00  2.8811E-01 -2.2368E-02
             1.4538E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1422.62978057171        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1106
 NPARAMETR:  9.6312E-01  1.0000E-02  8.2356E-01  1.5825E+00  5.7787E-01  9.7179E-01  1.1766E+01  1.0000E-02  8.1632E-01  9.3882E-01
             2.1564E+00
 PARAMETER:  6.2423E-02 -4.6091E+00 -9.4117E-02  5.5900E-01 -4.4840E-01  7.1382E-02  2.5653E+00 -1.7581E+01 -1.0295E-01  3.6863E-02
             8.6842E-01
 GRADIENT:   2.5972E-03  0.0000E+00 -4.8360E-02 -2.5944E-01  1.1297E-01  7.0489E-03  8.2641E-03  0.0000E+00  5.1930E-02 -1.0666E-02
            -7.7331E-03

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1422.62980535386        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  9.6312E-01  1.0000E-02  8.2361E-01  1.5826E+00  5.7786E-01  9.7177E-01  1.1730E+01  1.0000E-02  8.1617E-01  9.3886E-01
             2.1564E+00
 PARAMETER:  6.2427E-02 -4.6028E+00 -9.4059E-02  5.5909E-01 -4.4842E-01  7.1366E-02  2.5621E+00 -1.7554E+01 -1.0313E-01  3.6906E-02
             8.6846E-01
 GRADIENT:  -9.5583E-04  0.0000E+00  2.2256E-04 -1.0679E-03 -4.0293E-04  4.3208E-05 -2.4327E-03  0.0000E+00  1.7973E-04  1.2396E-03
             1.0789E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1163
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.5233E-07  9.7758E-04 -1.7846E-05 -1.0736E-02 -2.0745E-02
 SE:             2.9281E-02  1.9580E-03  1.7341E-04  2.7572E-02  2.2370E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9998E-01  6.1758E-01  9.1804E-01  6.9701E-01  3.5373E-01

 ETASHRINKSD(%)  1.9052E+00  9.3441E+01  9.9419E+01  7.6294E+00  2.5059E+01
 ETASHRINKVR(%)  3.7742E+00  9.9570E+01  9.9997E+01  1.4677E+01  4.3839E+01
 EBVSHRINKSD(%)  1.8496E+00  9.4216E+01  9.9403E+01  7.2201E+00  2.3999E+01
 EBVSHRINKVR(%)  3.6649E+00  9.9665E+01  9.9996E+01  1.3919E+01  4.2238E+01
 RELATIVEINF(%)  8.8260E+01  1.4406E-02  2.0775E-04  5.7029E+00  2.3654E+00
 EPSSHRINKSD(%)  3.4938E+01
 EPSSHRINKVR(%)  5.7670E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1422.6298053538553     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -687.47897879011714     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.02
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1422.630       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.63E-01  1.00E-02  8.24E-01  1.58E+00  5.78E-01  9.72E-01  1.17E+01  1.00E-02  8.16E-01  9.39E-01  2.16E+00
 


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
+        1.25E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -4.60E+01  0.00E+00  6.20E+02
 
 TH 4
+       -3.17E+01  0.00E+00 -7.65E+01  5.76E+02
 
 TH 5
+        1.75E+01  0.00E+00 -1.06E+03 -9.66E+01  2.17E+03
 
 TH 6
+       -9.55E+00  0.00E+00 -1.21E+01 -6.39E+00 -6.24E+00  2.12E+02
 
 TH 7
+        4.06E-03  0.00E+00  1.50E-02 -5.35E-02  4.90E-02 -1.13E-02  1.16E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.82E+00  0.00E+00  3.40E+01 -8.48E+00 -2.46E+00 -3.93E+01  2.91E-02  0.00E+00  2.52E+02
 
 TH10
+        2.67E+01  0.00E+00 -1.46E+01 -1.61E+00 -7.41E+01  5.22E+00 -7.64E-02  0.00E+00  1.64E+01  9.12E+01
 
 TH11
+       -1.24E+01  0.00E+00 -1.37E+01 -8.98E+00  2.66E+00  4.07E+00 -7.99E-03  0.00E+00  8.55E+00  2.05E+01  5.78E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.181
Stop Time:
Sat Sep 18 10:08:21 CDT 2021
