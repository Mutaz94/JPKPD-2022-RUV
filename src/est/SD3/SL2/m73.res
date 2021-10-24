Sun Oct 24 00:02:16 CDT 2021
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
$DATA ../../../../data/SD3/SL2/dat73.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m73.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2132.87256348018        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0070E+02  3.6710E+01  6.4194E+01  1.0621E+01 -1.0273E+02  8.4005E+01 -1.2885E+01 -9.2955E+00  5.1559E+00  3.9693E+00
            -4.4390E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2144.75861045526        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  1.0242E+00  1.0494E+00  9.6743E-01  1.0005E+00  1.0934E+00  8.2929E-01  1.0892E+00  1.0241E+00  9.9085E-01  9.9029E-01
             1.0086E+00
 PARAMETER:  1.2393E-01  1.4826E-01  6.6888E-02  1.0048E-01  1.8932E-01 -8.7191E-02  1.8543E-01  1.2381E-01  9.0813E-02  9.0243E-02
             1.0852E-01
 GRADIENT:   2.4454E+01 -1.1566E+01  3.5360E+00  1.9341E+00 -3.5770E+00 -1.5742E+01 -1.3806E+01 -3.7373E+00  2.7851E+00  1.6013E+00
             3.0833E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2146.82132998286        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0247E+00  1.1820E+00  8.9765E-01  9.2177E-01  1.1231E+00  8.3904E-01  1.1865E+00  1.0951E+00  9.6691E-01  9.4545E-01
             9.9888E-01
 PARAMETER:  1.2437E-01  2.6724E-01 -7.9762E-03  1.8540E-02  2.1611E-01 -7.5497E-02  2.7097E-01  1.9087E-01  6.6355E-02  4.3910E-02
             9.8882E-02
 GRADIENT:   2.1870E+01  6.8655E+00  4.1889E+00  4.9072E+00 -2.3074E+00 -1.1444E+01  1.1262E+00  6.3724E-02 -3.3864E+00 -1.1857E+00
            -2.0867E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2147.83243285665        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  1.0169E+00  1.3387E+00  6.7861E-01  8.1719E-01  1.0714E+00  8.7521E-01  1.0711E+00  7.9099E-01  1.0640E+00  9.0747E-01
             9.9814E-01
 PARAMETER:  1.1676E-01  3.9172E-01 -2.8771E-01 -1.0189E-01  1.6896E-01 -3.3292E-02  1.6873E-01 -1.3447E-01  1.6205E-01  2.9081E-03
             9.8140E-02
 GRADIENT:  -8.2308E+00  1.3418E+01  2.3242E+00  1.1241E+01 -1.2262E+01  4.5979E+00 -6.3052E-01  1.3192E+00 -6.8045E-01  1.3679E+00
             1.7937E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2149.17230043632        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  1.0206E+00  1.6876E+00  5.0120E-01  5.9870E-01  1.2216E+00  8.6406E-01  9.0679E-01  3.9666E-01  1.3110E+00  1.0195E+00
             9.9636E-01
 PARAMETER:  1.2040E-01  6.2333E-01 -5.9075E-01 -4.1299E-01  3.0019E-01 -4.6110E-02  2.1515E-03 -8.2468E-01  3.7082E-01  1.1932E-01
             9.6356E-02
 GRADIENT:  -1.4555E+00  2.3280E+01  9.4747E-01  1.5115E+01  3.1943E-01 -9.8621E-01 -2.9991E+00  4.2345E-01 -2.4612E+00 -9.6357E-01
            -3.3373E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2149.49850338218        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  1.0224E+00  1.8164E+00  4.4355E-01  5.0624E-01  1.3044E+00  8.6340E-01  8.6544E-01  2.4272E-01  1.4797E+00  1.0868E+00
             9.9555E-01
 PARAMETER:  1.2215E-01  6.9685E-01 -7.1294E-01 -5.8074E-01  3.6576E-01 -4.6874E-02 -4.4519E-02 -1.3159E+00  4.9183E-01  1.8327E-01
             9.5538E-02
 GRADIENT:   3.7480E+00  9.2547E+00 -4.4052E-01  5.7404E+00 -1.2020E+00 -1.3375E+00 -1.4035E+00  2.7628E-01  2.6633E-02 -3.3088E-01
            -6.6083E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2149.72776776499        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1069
 NPARAMETR:  1.0219E+00  1.8033E+00  4.4375E-01  5.0363E-01  1.3054E+00  8.6591E-01  8.7117E-01  3.3794E-02  1.4826E+00  1.0900E+00
             9.9560E-01
 PARAMETER:  1.2167E-01  6.8965E-01 -7.1250E-01 -5.8591E-01  3.6652E-01 -4.3976E-02 -3.7922E-02 -3.2875E+00  4.9377E-01  1.8614E-01
             9.5592E-02
 GRADIENT:   2.5024E+00 -8.3989E+00 -2.0845E-01 -9.6451E-01  9.1092E-01 -1.5188E-01 -4.0744E-01  5.2350E-03  3.0100E-01  9.4056E-02
             5.0270E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2149.73704404126        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1255
 NPARAMETR:  1.0218E+00  1.8012E+00  4.4441E-01  5.0574E-01  1.3026E+00  8.6615E-01  8.7157E-01  1.0000E-02  1.4775E+00  1.0864E+00
             9.9542E-01
 PARAMETER:  1.2155E-01  6.8847E-01 -7.1102E-01 -5.8173E-01  3.6436E-01 -4.3697E-02 -3.7457E-02 -5.2415E+00  4.9033E-01  1.8284E-01
             9.5407E-02
 GRADIENT:   2.1425E+00 -6.8432E+00 -6.8487E-02 -5.4495E-01  7.8517E-01 -6.4657E-02 -5.8060E-01  0.0000E+00  1.3172E-01 -1.1730E-01
            -1.5263E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2149.74454650968        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1447             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0219E+00  1.7977E+00  4.4500E-01  5.0743E-01  1.2999E+00  8.6648E-01  8.7583E-01  1.0000E-02  1.4749E+00  1.0858E+00
             9.9557E-01
 PARAMETER:  1.2168E-01  6.8650E-01 -7.0967E-01 -5.7840E-01  3.6232E-01 -4.3316E-02 -3.2589E-02 -5.2415E+00  4.8862E-01  1.8234E-01
             9.5562E-02
 GRADIENT:   5.4434E+02  8.1255E+02  4.9029E+00  1.2349E+02  1.9360E+01  3.3299E+01  9.2068E+00  0.0000E+00  2.1438E+01  1.4463E+00
             1.1610E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2149.74894572047        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1635
 NPARAMETR:  1.0219E+00  1.7958E+00  4.4553E-01  5.0871E-01  1.2985E+00  8.6648E-01  8.7609E-01  1.0000E-02  1.4724E+00  1.0843E+00
             9.9553E-01
 PARAMETER:  1.2168E-01  6.8544E-01 -7.0850E-01 -5.7588E-01  3.6124E-01 -4.3316E-02 -3.2292E-02 -5.2415E+00  4.8692E-01  1.8091E-01
             9.5524E-02
 GRADIENT:   2.5149E+00 -7.7875E+00 -4.5079E-01 -6.7211E-01  6.6824E-01  8.4069E-02  2.6074E-01  0.0000E+00  4.2064E-01  1.5890E-01
             1.2860E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2149.75346787612        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1822
 NPARAMETR:  1.0217E+00  1.7942E+00  4.4627E-01  5.1076E-01  1.2965E+00  8.6644E-01  8.7699E-01  1.0000E-02  1.4687E+00  1.0826E+00
             9.9553E-01
 PARAMETER:  1.2149E-01  6.8459E-01 -7.0684E-01 -5.7186E-01  3.5964E-01 -4.3367E-02 -3.1264E-02 -5.2415E+00  4.8441E-01  1.7933E-01
             9.5521E-02
 GRADIENT:   1.4179E+00 -5.8783E+00 -2.9741E-01 -7.2604E-02  3.8370E-01  1.6343E-01  2.3014E-01  0.0000E+00  3.5835E-01  1.0646E-01
             6.4257E-02

0ITERATION NO.:   53    OBJECTIVE VALUE:  -2149.75616640411        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1916
 NPARAMETR:  1.0219E+00  1.7921E+00  4.4692E-01  5.1052E-01  1.2956E+00  8.6649E-01  8.7598E-01  1.0000E-02  1.4664E+00  1.0816E+00
             9.9539E-01
 PARAMETER:  1.2169E-01  6.8339E-01 -7.0539E-01 -5.7233E-01  3.5898E-01 -4.3301E-02 -3.2413E-02 -5.2415E+00  4.8283E-01  1.7840E-01
             9.5376E-02
 GRADIENT:   2.5459E+00 -8.0167E+00  8.9138E-02 -1.4833E+00  7.2841E-02  9.4689E-02 -1.3270E-01  0.0000E+00  7.1821E-02  6.0091E-02
            -3.2925E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1916
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.0541E-04 -2.8501E-02 -2.8971E-04  3.3130E-02 -4.1600E-02
 SE:             2.9852E-02  2.5669E-02  1.0929E-04  2.2681E-02  2.1318E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8649E-01  2.6685E-01  8.0285E-03  1.4410E-01  5.1009E-02

 ETASHRINKSD(%)  1.0000E-10  1.4007E+01  9.9634E+01  2.4016E+01  2.8581E+01
 ETASHRINKVR(%)  1.0000E-10  2.6052E+01  9.9999E+01  4.2264E+01  4.8994E+01
 EBVSHRINKSD(%)  4.2609E-01  1.4185E+01  9.9660E+01  2.5533E+01  2.6065E+01
 EBVSHRINKVR(%)  8.5037E-01  2.6357E+01  9.9999E+01  4.4546E+01  4.5335E+01
 RELATIVEINF(%)  9.9105E+01  1.0814E+01  3.3074E-04  7.1017E+00  1.4441E+01
 EPSSHRINKSD(%)  3.4006E+01
 EPSSHRINKVR(%)  5.6447E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2149.7561664041109     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1230.8176331994382     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2149.756       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.79E+00  4.47E-01  5.11E-01  1.30E+00  8.66E-01  8.76E-01  1.00E-02  1.47E+00  1.08E+00  9.95E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       68.379
Stop Time:
Sun Oct 24 00:02:29 CDT 2021
