Wed Sep 29 14:19:28 CDT 2021
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
$DATA ../../../../data/spa/S1/dat50.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1696.61860074373        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2910E+02 -5.8880E+00 -3.6979E+01  1.9534E+01  6.3021E+00  4.1894E+01  1.7522E+01  1.9227E+01  2.6033E+00  4.1886E+01
             3.2711E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1707.46101521502        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.9419E-01  1.0161E+00  1.2204E+00  1.0233E+00  1.0558E+00  1.0009E+00  8.8642E-01  8.8064E-01  1.0326E+00  7.8085E-01
             1.0141E+00
 PARAMETER:  9.4171E-02  1.1601E-01  2.9918E-01  1.2307E-01  1.5428E-01  1.0091E-01 -2.0562E-02 -2.7105E-02  1.3207E-01 -1.4738E-01
             1.1396E-01
 GRADIENT:   1.5303E+00  1.8644E+01  2.9896E+01 -1.6042E+01 -3.1325E+01  8.2881E-01  1.1370E+01  1.7500E+00 -4.2821E+00 -4.4824E+00
            -4.6799E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1711.15366599057        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  1.0065E+00  9.5304E-01  1.0760E+00  1.0621E+00  9.9122E-01  9.9209E-01  5.8712E-01  5.7648E-01  1.0394E+00  7.8578E-01
             9.8691E-01
 PARAMETER:  1.0650E-01  5.1900E-02  1.7324E-01  1.6023E-01  9.1180E-02  9.2057E-02 -4.3252E-01 -4.5081E-01  1.3864E-01 -1.4108E-01
             8.6823E-02
 GRADIENT:   2.8618E+01  1.0495E+01  5.4617E+00  5.4055E+00 -2.1502E+00 -3.2349E+00 -2.0571E+00 -4.0617E-02 -8.5386E+00 -4.8420E+00
            -1.7193E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1712.16375248640        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.9341E-01  8.4986E-01  1.0670E+00  1.1280E+00  9.4749E-01  9.9804E-01  7.0468E-01  3.9281E-01  9.9972E-01  8.0188E-01
             1.0280E+00
 PARAMETER:  9.3384E-02 -6.2683E-02  1.6487E-01  2.2043E-01  4.6060E-02  9.8036E-02 -2.5001E-01 -8.3442E-01  9.9717E-02 -1.2080E-01
             1.2764E-01
 GRADIENT:  -1.0608E-01  4.8758E+00  2.3544E-01  6.5030E+00 -1.7554E+00 -9.5317E-02  5.0282E-01  5.5873E-01 -1.9439E-01  1.2490E+00
             1.1859E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1712.58131358276        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      700
 NPARAMETR:  9.9093E-01  6.0218E-01  1.0806E+00  1.2780E+00  8.6325E-01  9.9582E-01  7.6374E-01  1.1598E-01  8.9718E-01  8.0918E-01
             1.0270E+00
 PARAMETER:  9.0886E-02 -4.0720E-01  1.7755E-01  3.4533E-01 -4.7056E-02  9.5812E-02 -1.6953E-01 -2.0543E+00 -8.4958E-03 -1.1173E-01
             1.2661E-01
 GRADIENT:  -2.1744E-01  3.1703E+00  3.0243E+00  3.2003E+00 -5.9543E+00  1.5327E-01 -1.4644E-01  6.2859E-02 -2.6402E-01  1.4689E-01
            -2.8836E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1712.59832457208        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      875
 NPARAMETR:  9.9036E-01  5.1513E-01  1.0831E+00  1.3321E+00  8.3774E-01  9.9400E-01  7.9524E-01  5.5709E-02  8.6488E-01  8.1154E-01
             1.0248E+00
 PARAMETER:  9.0317E-02 -5.6334E-01  1.7986E-01  3.8676E-01 -7.7051E-02  9.3984E-02 -1.2911E-01 -2.7876E+00 -4.5162E-02 -1.0883E-01
             1.2450E-01
 GRADIENT:   1.0262E+00  3.0898E+00  1.5836E+00  6.4602E+00 -3.4617E+00 -9.8986E-02 -2.0814E-01  1.4954E-02 -6.2376E-01  2.0660E-01
            -5.9056E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1712.62585090307        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1056
 NPARAMETR:  9.9007E-01  4.9103E-01  1.0827E+00  1.3411E+00  8.3195E-01  9.9395E-01  8.4805E-01  2.6592E-02  8.5850E-01  8.1051E-01
             1.0252E+00
 PARAMETER:  9.0018E-02 -6.1125E-01  1.7945E-01  3.9351E-01 -8.3979E-02  9.3928E-02 -6.4816E-02 -3.5271E+00 -5.2568E-02 -1.1009E-01
             1.2489E-01
 GRADIENT:   1.2375E+00 -6.2253E-01 -7.9245E-01 -4.7757E+00  1.8779E+00  3.7667E-02  2.8017E-02  3.6199E-03  7.6501E-01  1.9295E-01
             2.6110E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1712.62870184797        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  9.8998E-01  4.9042E-01  1.0830E+00  1.3423E+00  8.3122E-01  9.9405E-01  8.5025E-01  1.0000E-02  8.5601E-01  8.0960E-01
             1.0254E+00
 PARAMETER:  8.9930E-02 -6.1248E-01  1.7976E-01  3.9436E-01 -8.4862E-02  9.4031E-02 -6.2230E-02 -4.7037E+00 -5.5477E-02 -1.1122E-01
             1.2512E-01
 GRADIENT:   1.0392E+00  2.0248E-01  5.3642E-01 -3.2867E+00 -7.1727E-02  8.2506E-02 -1.2355E-02  0.0000E+00 -5.4173E-02  2.5743E-02
            -5.8843E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1712.62912425896        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1420
 NPARAMETR:  9.9014E-01  4.9055E-01  1.0825E+00  1.3419E+00  8.3092E-01  9.9411E-01  8.5147E-01  1.0000E-02  8.5596E-01  8.0926E-01
             1.0254E+00
 PARAMETER:  9.0087E-02 -6.1222E-01  1.7932E-01  3.9410E-01 -8.5216E-02  9.4090E-02 -6.0787E-02 -4.7037E+00 -5.5535E-02 -1.1163E-01
             1.2513E-01
 GRADIENT:   1.3837E+00  1.4009E-01  6.9286E-01 -3.7840E+00 -3.1298E-01  1.0504E-01 -1.0805E-02  0.0000E+00 -7.2330E-02  1.8430E-02
            -5.0290E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1712.62946074607        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1609
 NPARAMETR:  9.9015E-01  4.9065E-01  1.0820E+00  1.3417E+00  8.3070E-01  9.9412E-01  8.5368E-01  1.0000E-02  8.5602E-01  8.0899E-01
             1.0254E+00
 PARAMETER:  9.0097E-02 -6.1203E-01  1.7881E-01  3.9395E-01 -8.5483E-02  9.4101E-02 -5.8202E-02 -4.7037E+00 -5.5461E-02 -1.1196E-01
             1.2511E-01
 GRADIENT:   1.3971E+00  6.9414E-02  6.2013E-01 -4.0276E+00 -2.4015E-01  1.0755E-01 -3.3562E-03  0.0000E+00 -2.8239E-02  2.4489E-02
            -3.6641E-02

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1712.62957153999        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1708
 NPARAMETR:  9.9015E-01  4.9062E-01  1.0814E+00  1.3416E+00  8.3078E-01  9.9412E-01  8.5491E-01  1.0000E-02  8.5613E-01  8.0890E-01
             1.0255E+00
 PARAMETER:  9.0097E-02 -6.1208E-01  1.7828E-01  3.9387E-01 -8.5388E-02  9.4099E-02 -5.6757E-02 -4.7037E+00 -5.5334E-02 -1.1208E-01
             1.2519E-01
 GRADIENT:   1.3873E+00 -1.7614E-01 -1.4578E-01 -4.2129E+00  8.5474E-01  1.0444E-01  1.4156E-03  0.0000E+00  2.8583E-02  2.6336E-02
             5.0028E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1708
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.6626E-04 -1.0158E-02 -2.1180E-04 -2.4913E-03 -2.0178E-02
 SE:             2.9822E-02  8.0297E-03  2.1148E-04  2.8417E-02  2.4839E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9555E-01  2.0583E-01  3.1660E-01  9.3014E-01  4.1661E-01

 ETASHRINKSD(%)  9.2440E-02  7.3100E+01  9.9292E+01  4.7985E+00  1.6785E+01
 ETASHRINKVR(%)  1.8479E-01  9.2764E+01  9.9995E+01  9.3667E+00  3.0753E+01
 EBVSHRINKSD(%)  4.4201E-01  7.3662E+01  9.9281E+01  4.7283E+00  1.5469E+01
 EBVSHRINKVR(%)  8.8206E-01  9.3063E+01  9.9995E+01  9.2331E+00  2.8545E+01
 RELATIVEINF(%)  9.7459E+01  2.2482E-01  6.5715E-04  4.7134E+00  3.4295E+00
 EPSSHRINKSD(%)  4.1541E+01
 EPSSHRINKVR(%)  6.5825E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1712.6295715399890     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -977.47874497625082     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1712.630       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.90E-01  4.91E-01  1.08E+00  1.34E+00  8.31E-01  9.94E-01  8.55E-01  1.00E-02  8.56E-01  8.09E-01  1.03E+00
 


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
+        1.14E+03
 
 TH 2
+       -2.27E+01  4.49E+02
 
 TH 3
+        1.03E+01  2.29E+02  5.00E+02
 
 TH 4
+       -9.95E+00  4.49E+02 -5.25E+01  7.65E+02
 
 TH 5
+       -8.92E-01 -5.38E+02 -8.60E+02 -9.47E+00  1.76E+03
 
 TH 6
+        1.28E-01 -3.81E+00  2.27E+00 -2.51E+00 -1.54E+00  1.98E+02
 
 TH 7
+        2.38E-01 -3.10E+00  2.93E+00 -3.64E+00 -1.22E+00  2.16E-02  2.80E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.55E+00 -6.30E+01  6.75E+00  9.08E+00 -8.73E+00 -4.82E-01  9.69E+00  0.00E+00  2.24E+02
 
 TH10
+       -9.19E-01  2.35E+01 -2.02E+01 -5.07E+00 -8.30E+01  6.83E-01  5.86E+00  0.00E+00  1.17E-01  1.52E+02
 
 TH11
+       -7.71E+00 -1.74E+01 -3.95E+01 -8.81E+00  1.63E+01  1.57E+00  1.53E+00  0.00E+00  9.36E+00  4.02E+01  2.17E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.897
Stop Time:
Wed Sep 29 14:19:57 CDT 2021
