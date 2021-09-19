Sat Sep 18 08:32:52 CDT 2021
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
$DATA ../../../../data/spa/B/dat50.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1690.47054049652        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7427E+01 -4.2117E+01 -2.5863E+01 -5.2560E+01 -5.1120E+00 -1.6528E-01  7.8099E+00  1.8297E+01  6.2061E+00  2.4209E+01
             6.7955E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1695.70594681011        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7878E-01  1.0306E+00  1.1281E+00  1.0203E+00  1.0634E+00  1.0036E+00  9.5290E-01  9.0217E-01  9.6755E-01  8.9243E-01
             9.8260E-01
 PARAMETER:  7.8556E-02  1.3017E-01  2.2052E-01  1.2008E-01  1.6147E-01  1.0363E-01  5.1756E-02 -2.9529E-03  6.7011E-02 -1.3811E-02
             8.2447E-02
 GRADIENT:   3.4307E+00  6.3198E+00  1.7751E+00 -5.1623E+00  6.6107E+00  1.1573E+00  3.7507E+00  7.4075E+00 -5.7075E+00 -7.5954E+00
            -9.9554E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1697.14404130951        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  9.8579E-01  9.8060E-01  1.1618E+00  1.0557E+00  1.0390E+00  1.0063E+00  8.7267E-01  7.4104E-01  9.5829E-01  9.7277E-01
             9.5047E-01
 PARAMETER:  8.5683E-02  8.0405E-02  2.5000E-01  1.5424E-01  1.3822E-01  1.0626E-01 -3.6194E-02 -1.9970E-01  5.7395E-02  7.2396E-02
             4.9200E-02
 GRADIENT:  -2.2644E+01  1.8300E+01  1.8294E+01 -1.0538E+01 -3.2264E+01 -1.8924E+00  4.4065E-02  3.6449E+00 -1.1216E+01 -1.9526E-01
            -2.3851E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.18600895783        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  9.9702E-01  8.5398E-01  1.0514E+00  1.1277E+00  9.6534E-01  1.0096E+00  1.0072E+00  4.0131E-01  9.3184E-01  9.2891E-01
             9.9667E-01
 PARAMETER:  9.7020E-02 -5.7842E-02  1.5013E-01  2.2020E-01  6.4724E-02  1.0952E-01  1.0715E-01 -8.1301E-01  2.9402E-02  2.6261E-02
             9.6662E-02
 GRADIENT:   1.5332E+00 -6.3674E+00 -1.2207E+01 -5.9802E+00  1.5895E+01 -1.9343E-01  5.8525E-01  1.3700E+00  9.3472E-01  1.0294E+00
             1.1883E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1701.64191239608        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      621
 NPARAMETR:  9.9210E-01  5.4231E-01  1.2231E+00  1.3327E+00  9.2318E-01  1.0044E+00  8.3520E-01  1.3454E-01  8.6240E-01  1.0372E+00
             9.9823E-01
 PARAMETER:  9.2070E-02 -5.1192E-01  3.0136E-01  3.8717E-01  2.0073E-02  1.0440E-01 -8.0085E-02 -1.9059E+00 -4.8041E-02  1.3652E-01
             9.8223E-02
 GRADIENT:   2.6184E-01  6.2178E+00  6.2759E+00  1.0448E+01 -1.1947E+01 -1.0496E-01  1.5350E-01  7.9311E-03  9.9738E-01  3.5635E+00
            -5.7771E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1702.19783208861        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      796
 NPARAMETR:  9.8846E-01  3.2257E-01  1.2382E+00  1.4601E+00  8.6775E-01  1.0007E+00  6.0205E-01  1.6863E-02  7.9982E-01  1.0190E+00
             9.9857E-01
 PARAMETER:  8.8391E-02 -1.0314E+00  3.1368E-01  4.7853E-01 -4.1852E-02  1.0075E-01 -4.0742E-01 -3.9826E+00 -1.2337E-01  1.1880E-01
             9.8572E-02
 GRADIENT:  -4.7763E-01  1.1844E+00  1.7054E+00  4.2492E+00 -2.2844E+00 -1.0559E-01 -4.0739E-02 -2.6928E-04  8.2072E-01 -3.2132E-01
            -7.0328E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1702.30621502621        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      971
 NPARAMETR:  9.8606E-01  1.9650E-01  1.2446E+00  1.5332E+00  8.3706E-01  9.9843E-01  4.0258E-01  1.0000E-02  7.6128E-01  1.0174E+00
             9.9865E-01
 PARAMETER:  8.5961E-02 -1.5271E+00  3.1879E-01  5.2733E-01 -7.7856E-02  9.8433E-02 -8.0987E-01 -6.3076E+00 -1.7275E-01  1.1723E-01
             9.8649E-02
 GRADIENT:  -8.1677E-01  4.7741E-02 -2.7846E-02  4.3236E-01  2.2922E-01 -9.3311E-02 -2.2164E-03  0.0000E+00  4.9206E-01  1.3359E-01
            -1.9451E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1702.31047082318        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1148
 NPARAMETR:  9.8587E-01  1.6641E-01  1.2417E+00  1.5504E+00  8.2751E-01  9.9809E-01  3.4849E-01  1.0000E-02  7.5132E-01  1.0135E+00
             9.9909E-01
 PARAMETER:  8.5773E-02 -1.6933E+00  3.1645E-01  5.3849E-01 -8.9332E-02  9.8089E-02 -9.5415E-01 -7.1306E+00 -1.8592E-01  1.1342E-01
             9.9094E-02
 GRADIENT:   2.0773E-02  1.2617E-02  2.7897E-02  6.1871E-02 -8.1404E-02  1.1390E-03 -9.1997E-04  0.0000E+00 -1.0698E-02  3.2763E-03
             3.8414E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1702.31050916968        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1335
 NPARAMETR:  9.8583E-01  1.6616E-01  1.2419E+00  1.5506E+00  8.2764E-01  9.9809E-01  3.7550E-01  1.0000E-02  7.5128E-01  1.0136E+00
             9.9910E-01
 PARAMETER:  8.5729E-02 -1.6948E+00  3.1661E-01  5.3864E-01 -8.9180E-02  9.8088E-02 -8.7950E-01 -7.1306E+00 -1.8597E-01  1.1354E-01
             9.9098E-02
 GRADIENT:  -6.4532E-02  6.1656E-03 -1.0208E-01  1.9684E-01  1.1953E-01  2.3554E-03 -8.0894E-04  0.0000E+00  3.9371E-02  7.7331E-03
             1.3722E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1702.31063116977        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1513
 NPARAMETR:  9.8585E-01  1.6539E-01  1.2420E+00  1.5510E+00  8.2741E-01  9.9807E-01  4.9601E-01  1.0000E-02  7.5061E-01  1.0133E+00
             9.9910E-01
 PARAMETER:  8.5754E-02 -1.6994E+00  3.1672E-01  5.3892E-01 -8.9452E-02  9.8064E-02 -6.0116E-01 -7.1306E+00 -1.8687E-01  1.1324E-01
             9.9095E-02
 GRADIENT:   2.6550E-02 -4.2643E-03  3.8732E-02 -2.8556E-02  4.4294E-03 -1.3227E-03  2.9810E-05  0.0000E+00  3.4618E-03 -7.8146E-03
            -9.1350E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1702.31067343453        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1683
 NPARAMETR:  9.8586E-01  1.6574E-01  1.2408E+00  1.5507E+00  8.2696E-01  9.9808E-01  5.0345E-01  1.0000E-02  7.5069E-01  1.0129E+00
             9.9906E-01
 PARAMETER:  8.5757E-02 -1.6974E+00  3.1579E-01  5.3870E-01 -8.9993E-02  9.8078E-02 -5.8628E-01 -7.1306E+00 -1.8676E-01  1.1281E-01
             9.9061E-02
 GRADIENT:   1.5431E-02 -7.8670E-03  1.3488E-02 -1.2678E-01 -4.0267E-02  1.7153E-03 -3.3239E-05  0.0000E+00  5.4364E-04  9.0519E-04
            -1.2269E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1683
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.7330E-05 -2.8796E-03 -2.0696E-04 -5.4851E-03 -2.2380E-02
 SE:             2.9829E-02  1.5604E-03  1.8118E-04  2.9177E-02  2.5316E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9900E-01  6.4970E-02  2.5332E-01  8.5088E-01  3.7667E-01

 ETASHRINKSD(%)  7.0249E-02  9.4773E+01  9.9393E+01  2.2540E+00  1.5189E+01
 ETASHRINKVR(%)  1.4045E-01  9.9727E+01  9.9996E+01  4.4572E+00  2.8071E+01
 EBVSHRINKSD(%)  4.0898E-01  9.5103E+01  9.9400E+01  2.3729E+00  1.2311E+01
 EBVSHRINKVR(%)  8.1629E-01  9.9760E+01  9.9996E+01  4.6894E+00  2.3106E+01
 RELATIVEINF(%)  9.6826E+01  1.0838E-02  4.0377E-04  5.7205E+00  5.4493E+00
 EPSSHRINKSD(%)  4.1700E+01
 EPSSHRINKVR(%)  6.6011E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1702.3106734345279     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -967.15984687078969     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.97
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1702.311       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.66E-01  1.24E+00  1.55E+00  8.27E-01  9.98E-01  5.03E-01  1.00E-02  7.51E-01  1.01E+00  9.99E-01
 


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
+        1.13E+03
 
 TH 2
+       -3.22E+01  3.95E+02
 
 TH 3
+        7.64E-01  1.08E+02  2.63E+02
 
 TH 4
+       -8.89E+00  4.67E+02 -5.09E+01  7.82E+02
 
 TH 5
+        1.93E+00 -3.40E+02 -5.21E+02 -3.77E+01  1.22E+03
 
 TH 6
+        2.07E+00 -3.74E+00  2.75E+00 -3.10E+00  7.13E+00  1.89E+02
 
 TH 7
+        1.53E+00 -3.40E-01  5.50E-01 -5.54E-01 -8.03E-01  9.99E-01 -2.46E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.34E+00 -9.62E+01  7.35E+00 -1.02E+00 -4.44E+00  4.10E+00  1.06E+00  0.00E+00  3.20E+02
 
 TH10
+        2.74E-01  1.39E+01 -1.12E+01 -2.12E+00 -7.17E+01  6.13E+00  2.82E+00  0.00E+00 -3.53E+00  1.03E+02
 
 TH11
+       -1.24E+01 -1.80E+01 -3.42E+01 -1.04E+01  2.70E+01 -7.15E-01 -2.08E+00  0.00E+00  8.75E+00  3.11E+01  2.23E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.272
Stop Time:
Sat Sep 18 08:33:19 CDT 2021
