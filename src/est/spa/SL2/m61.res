Wed Sep 29 16:00:52 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat61.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1698.48534122653        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6220E+02 -2.1206E+01  1.8079E+00 -1.9079E+01 -1.4691E+01  5.7249E+01  1.0386E+01  1.1880E+01  2.1143E+01  7.0750E+00
            -1.5378E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1705.19253109701        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      195
 NPARAMETR:  1.0307E+00  1.0529E+00  1.0485E+00  1.0416E+00  1.0475E+00  9.8597E-01  9.4980E-01  9.2938E-01  9.0756E-01  9.9289E-01
             1.0589E+00
 PARAMETER:  1.3020E-01  1.5155E-01  1.4736E-01  1.4077E-01  1.4641E-01  8.5873E-02  4.8499E-02  2.6763E-02  3.0096E-03  9.2864E-02
             1.5727E-01
 GRADIENT:  -1.9962E+00  1.4272E+01  4.1125E+00  1.0986E+01 -8.3295E+00 -6.5092E-01  2.5512E+00  5.4008E+00 -1.6409E+00 -4.6945E+00
             2.6052E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1706.63748093477        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0365E+00  9.8298E-01  1.0076E+00  1.0778E+00  1.0137E+00  9.8420E-01  8.3812E-01  6.6143E-01  9.2488E-01  1.0508E+00
             1.0409E+00
 PARAMETER:  1.3584E-01  8.2830E-02  1.0758E-01  1.7497E-01  1.1361E-01  8.4077E-02 -7.6594E-02 -3.1335E-01  2.1908E-02  1.4957E-01
             1.4010E-01
 GRADIENT:   1.1603E+01 -1.0172E-01 -7.0660E+00  9.4389E+00  4.3989E+00 -1.1210E+00 -2.3129E+00  2.5843E+00  4.9433E-02 -1.4326E+00
            -3.3950E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1707.46626953740        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0290E+00  8.8175E-01  9.5678E-01  1.1367E+00  9.3422E-01  9.8697E-01  1.0831E+00  3.6491E-01  8.4185E-01  1.0072E+00
             1.0540E+00
 PARAMETER:  1.2861E-01 -2.5842E-02  5.5813E-02  2.2811E-01  3.1956E-02  8.6880E-02  1.7979E-01 -9.0810E-01 -7.2158E-02  1.0720E-01
             1.5259E-01
 GRADIENT:  -4.7833E+00  6.0302E+00  3.2062E+00  7.8828E+00 -6.2614E+00  2.8815E-01 -6.2601E-01  5.0976E-01 -8.7195E-01 -1.1221E-02
             1.7369E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1707.71162588702        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.0310E+00  8.2946E-01  9.1874E-01  1.1594E+00  8.9491E-01  9.8624E-01  1.1816E+00  1.3927E-01  8.1959E-01  9.8428E-01
             1.0496E+00
 PARAMETER:  1.3049E-01 -8.6986E-02  1.5245E-02  2.4786E-01 -1.1033E-02  8.6144E-02  2.6685E-01 -1.8714E+00 -9.8956E-02  8.4154E-02
             1.4844E-01
 GRADIENT:  -2.7927E-01 -1.9575E-01  1.4396E-01 -3.6527E-01 -7.2275E-01  7.7969E-02  5.2058E-02  7.3779E-02  3.2732E-01  1.2491E-01
             2.3364E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1707.73384358194        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  1.0315E+00  8.1827E-01  9.1868E-01  1.1658E+00  8.9066E-01  9.8593E-01  1.1979E+00  4.6915E-02  8.1370E-01  9.8373E-01
             1.0497E+00
 PARAMETER:  1.3104E-01 -1.0057E-01  1.5179E-02  2.5344E-01 -1.5797E-02  8.5828E-02  2.8057E-01 -2.9594E+00 -1.0616E-01  8.3598E-02
             1.4851E-01
 GRADIENT:   1.2129E+00 -1.2342E-01  7.2916E-01 -4.8744E-01 -5.3644E-01  2.6325E-02 -1.7346E-03  7.3579E-03 -9.2631E-02 -2.4625E-01
            -9.3149E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1707.73978402022        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  1.0305E+00  8.2357E-01  9.1634E-01  1.1626E+00  8.9131E-01  9.8618E-01  1.1919E+00  1.0000E-02  8.1572E-01  9.8386E-01
             1.0492E+00
 PARAMETER:  1.3009E-01 -9.4101E-02  1.2637E-02  2.5066E-01 -1.5065E-02  8.6087E-02  2.7554E-01 -4.8504E+00 -1.0369E-01  8.3730E-02
             1.4798E-01
 GRADIENT:  -1.1276E+00  7.4780E-02  7.9839E-01 -3.5247E-01 -8.8319E-01  8.5845E-02  2.9334E-02  0.0000E+00 -4.9226E-02 -1.2882E-01
            -2.6043E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1707.76778151094        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  1.0318E+00  8.5831E-01  9.0267E-01  1.1410E+00  8.9753E-01  9.8652E-01  1.1514E+00  1.0000E-02  8.2625E-01  9.8135E-01
             1.0492E+00
 PARAMETER:  1.3134E-01 -5.2792E-02 -2.4036E-03  2.3189E-01 -8.1134E-03  8.6429E-02  2.4094E-01 -1.4142E+01 -9.0862E-02  8.1172E-02
             1.4804E-01
 GRADIENT:   8.2397E-01  8.3107E-02  1.3901E-01 -4.6229E-01 -1.0987E+00 -1.0898E-02 -6.8978E-02  0.0000E+00 -4.7002E-01  2.7579E-01
             7.3695E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1707.80748636554        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1418
 NPARAMETR:  1.0319E+00  9.4410E-01  8.8165E-01  1.0880E+00  9.2460E-01  9.8742E-01  1.0584E+00  1.0000E-02  8.6365E-01  9.8775E-01
             1.0460E+00
 PARAMETER:  1.3140E-01  4.2479E-02 -2.5958E-02  1.8437E-01  2.1604E-02  8.7345E-02  1.5677E-01 -3.6633E+01 -4.6592E-02  8.7678E-02
             1.4500E-01
 GRADIENT:  -6.2700E-01 -1.0776E+00 -1.6816E+00 -4.8785E-01  9.2040E-01 -1.0283E-01  9.8912E-02  0.0000E+00  1.8645E-01  5.5282E-01
            -8.4764E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1707.83033616173        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1596
 NPARAMETR:  1.0326E+00  1.0062E+00  8.8095E-01  1.0514E+00  9.5343E-01  9.8816E-01  9.9217E-01  1.0000E-02  8.9298E-01  1.0016E+00
             1.0510E+00
 PARAMETER:  1.3205E-01  1.0614E-01 -2.6757E-02  1.5008E-01  5.2307E-02  8.8089E-02  9.2142E-02 -5.1723E+01 -1.3187E-02  1.0164E-01
             1.4973E-01
 GRADIENT:   1.5220E-01  2.1841E-01  4.6440E-01 -4.8340E-02 -2.2848E-01  5.4442E-02 -2.3947E-02  0.0000E+00 -5.1222E-02 -2.6530E-01
             4.8999E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1707.83174791457        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1764
 NPARAMETR:  1.0341E+00  1.0037E+00  8.8026E-01  1.0525E+00  9.5226E-01  9.8837E-01  9.9478E-01  1.0000E-02  8.9203E-01  1.0026E+00
             1.0496E+00
 PARAMETER:  1.3349E-01  1.0374E-01 -2.7533E-02  1.5119E-01  5.1085E-02  8.8307E-02  9.4766E-02 -5.1197E+01 -1.4252E-02  1.0260E-01
             1.4845E-01
 GRADIENT:   3.4881E+00 -3.4735E-01 -9.5728E-02 -3.6460E-01  2.2066E-01  1.3610E-01  2.4961E-02  0.0000E+00  3.2975E-02  6.7074E-02
             7.3736E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1764
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -7.2198E-04 -1.0855E-02 -3.9082E-04  1.0135E-03 -2.3247E-02
 SE:             2.9828E-02  1.8390E-02  1.7342E-04  2.4952E-02  2.3821E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8069E-01  5.5501E-01  2.4223E-02  9.6760E-01  3.2911E-01

 ETASHRINKSD(%)  7.2776E-02  3.8391E+01  9.9419E+01  1.6407E+01  2.0198E+01
 ETASHRINKVR(%)  1.4550E-01  6.2043E+01  9.9997E+01  3.0122E+01  3.6316E+01
 EBVSHRINKSD(%)  4.6805E-01  3.8089E+01  9.9456E+01  1.6723E+01  1.8172E+01
 EBVSHRINKVR(%)  9.3392E-01  6.1671E+01  9.9997E+01  3.0649E+01  3.3041E+01
 RELATIVEINF(%)  9.8414E+01  1.2455E+00  3.4977E-04  2.8572E+00  6.3661E+00
 EPSSHRINKSD(%)  4.2495E+01
 EPSSHRINKVR(%)  6.6932E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1707.8317479145680     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -972.68092135082986     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.33
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.76
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1707.832       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.00E+00  8.80E-01  1.05E+00  9.52E-01  9.88E-01  9.95E-01  1.00E-02  8.92E-01  1.00E+00  1.05E+00
 


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
+        1.05E+03
 
 TH 2
+       -1.01E+01  4.49E+02
 
 TH 3
+        1.62E+01  2.00E+02  4.52E+02
 
 TH 4
+       -1.31E+01  4.27E+02 -1.74E+02  8.76E+02
 
 TH 5
+       -3.09E+00 -3.10E+02 -5.12E+02  1.60E+02  8.14E+02
 
 TH 6
+       -1.15E+00 -2.11E+00  4.00E+00 -2.44E+00 -6.47E-01  2.00E+02
 
 TH 7
+        1.05E+00  1.92E+01  8.13E+00 -6.72E+00 -1.14E+01 -1.96E-01  3.41E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.52E+00 -2.61E+01 -2.46E+01  2.91E+01  3.91E+00 -6.84E-01  2.51E+01  0.00E+00  1.27E+02
 
 TH10
+       -6.80E-01 -3.15E+00 -4.66E+01 -1.43E+01 -5.19E+01  4.00E-01  1.30E+01  0.00E+00  5.78E+00  8.64E+01
 
 TH11
+       -7.64E+00 -1.94E+01 -4.18E+01 -6.26E+00  9.15E+00  2.28E+00  6.17E+00  0.00E+00  1.16E+01  2.09E+01  2.02E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.115
Stop Time:
Wed Sep 29 16:01:22 CDT 2021
