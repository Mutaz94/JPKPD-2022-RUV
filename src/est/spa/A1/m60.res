Sat Sep 18 09:20:17 CDT 2021
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
$DATA ../../../../data/spa/A1/dat60.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -948.852162235675        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3276E+01 -5.3334E+01  2.2780E+01 -6.6170E+01  7.6495E+01 -7.3130E+00 -2.4863E+01 -3.7881E+01 -1.7609E+01 -1.1453E+01
            -1.3123E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1405.96458141914        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0246E+00  1.1773E+00  1.0753E+00  9.7719E-01  1.1038E+00  8.9223E-01  1.0261E+00  1.0650E+00  7.8639E-01  6.2705E-01
             2.7551E+00
 PARAMETER:  1.2434E-01  2.6319E-01  1.7259E-01  7.6922E-02  1.9880E-01 -1.4026E-02  1.2580E-01  1.6301E-01 -1.4030E-01 -3.6672E-01
             1.1135E+00
 GRADIENT:   2.6004E+01 -4.2113E+00 -1.2615E+01  1.4058E+01  2.4040E+01 -2.7118E+01 -6.8057E+00  8.3952E-01 -5.5982E+00  4.4613E+00
            -8.3835E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1409.47361659446        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0228E+00  1.1128E+00  9.0963E-01  1.0078E+00  9.5543E-01  1.0059E+00  1.1130E+00  9.2223E-01  8.9590E-01  3.2813E-01
             2.7476E+00
 PARAMETER:  1.2259E-01  2.0688E-01  5.2790E-03  1.0774E-01  5.4411E-02  1.0584E-01  2.0704E-01  1.9041E-02 -9.9228E-03 -1.0143E+00
             1.1107E+00
 GRADIENT:   1.5278E+01  4.5706E+00  4.2945E+00  7.8377E+00 -4.4312E-01  1.0925E+01 -1.0580E+00 -2.8942E+00  5.2354E+00  1.0252E+00
             1.3094E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1410.96257883488        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      237
 NPARAMETR:  1.0098E+00  1.2012E+00  9.5479E-01  9.5684E-01  1.0264E+00  9.4887E-01  1.0624E+00  1.5254E+00  8.9455E-01  2.9492E-01
             2.6701E+00
 PARAMETER:  1.0974E-01  2.8330E-01  5.3740E-02  5.5882E-02  1.2603E-01  4.7512E-02  1.6052E-01  5.2223E-01 -1.1434E-02 -1.1210E+00
             1.0821E+00
 GRADIENT:  -1.5292E+01  1.4072E+00 -8.6315E+00  1.6748E+01  1.6722E+01 -1.5346E+01 -1.0230E-01  7.6447E-02  2.1164E+00  1.2949E+00
             1.1250E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1412.53380738391        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  1.0131E+00  1.2350E+00  9.4201E-01  9.1868E-01  1.0093E+00  9.8832E-01  1.0591E+00  1.6843E+00  8.8652E-01  6.1099E-02
             2.5919E+00
 PARAMETER:  1.1303E-01  3.1108E-01  4.0258E-02  1.5182E-02  1.0921E-01  8.8250E-02  1.5737E-01  6.2137E-01 -2.0451E-02 -2.6953E+00
             1.0524E+00
 GRADIENT:  -2.5913E+00 -4.8047E-01  1.5692E+00 -6.7529E-01 -1.0385E+00 -1.1356E+00  4.7054E-02  1.8001E-01  2.6966E-01  5.6634E-02
            -1.9305E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1415.67578643494        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      384
 NPARAMETR:  1.0147E+00  1.2943E+00  2.8582E-01  7.9562E-01  6.3340E-01  9.8979E-01  9.6842E-01  7.2967E-01  8.6984E-01  1.0000E-02
             2.4781E+00
 PARAMETER:  1.1459E-01  3.5801E-01 -1.1524E+00 -1.2863E-01 -3.5666E-01  8.9734E-02  6.7909E-02 -2.1517E-01 -3.9444E-02 -4.8114E+00
             1.0075E+00
 GRADIENT:   4.7284E+00  2.6159E+01  7.0023E+00  2.5673E+01 -1.4787E+01 -2.2631E+00 -1.1953E-02 -1.9491E+00  5.1082E+00  0.0000E+00
            -4.2208E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1460.84453990688        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      458
 NPARAMETR:  9.7526E-01  1.8197E+00  1.1707E-01  4.2633E-01  8.4583E-01  1.0027E+00  7.8995E-01  1.9215E+00  1.2201E+00  2.9346E-02
             2.2583E+00
 PARAMETER:  7.4954E-02  6.9869E-01 -2.0450E+00 -7.5255E-01 -6.7442E-02  1.0273E-01 -1.3578E-01  7.5312E-01  2.9893E-01 -3.4286E+00
             9.1461E-01
 GRADIENT:  -4.8944E+01 -4.0535E+01 -3.9416E+00  2.2708E+01  7.4726E+01  8.3484E+00  1.3278E+00 -2.0866E+01 -2.4906E+00  1.5048E-02
             5.6633E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1469.07635271053        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  9.9632E-01  1.7843E+00  1.2767E-01  4.3587E-01  8.0298E-01  9.6838E-01  8.1361E-01  2.1963E+00  1.3721E+00  1.0000E-02
             1.9421E+00
 PARAMETER:  9.6314E-02  6.7905E-01 -1.9583E+00 -7.3042E-01 -1.1943E-01  6.7874E-02 -1.0628E-01  8.8678E-01  4.1635E-01 -5.5842E+00
             7.6377E-01
 GRADIENT:   1.1440E+01  2.3506E+00  3.1145E+00 -6.2381E+00  2.8697E+00 -2.1020E+00  2.4328E+00 -5.1279E+00  4.8773E+00  0.0000E+00
             8.4559E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1469.80301309627        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      599
 NPARAMETR:  9.8782E-01  1.6277E+00  1.6054E-01  5.2915E-01  7.3043E-01  9.7243E-01  8.4697E-01  2.3180E+00  1.1013E+00  1.0000E-02
             1.9011E+00
 PARAMETER:  8.7740E-02  5.8717E-01 -1.7292E+00 -5.3649E-01 -2.1413E-01  7.2040E-02 -6.6092E-02  9.4071E-01  1.9651E-01 -6.9323E+00
             7.4242E-01
 GRADIENT:  -6.8432E+00  1.4283E+00  4.6390E+00 -9.7846E+00 -6.9492E+00  1.6595E-01 -1.1788E+00  1.4594E+00 -7.4405E-01  0.0000E+00
             2.1187E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1469.80456625929        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      673
 NPARAMETR:  9.8810E-01  1.6208E+00  1.6046E-01  5.3316E-01  7.2660E-01  9.7244E-01  8.4905E-01  2.3042E+00  1.0970E+00  1.0000E-02
             1.8972E+00
 PARAMETER:  8.8028E-02  5.8291E-01 -1.7297E+00 -5.2894E-01 -2.1939E-01  7.2056E-02 -6.3641E-02  9.3475E-01  1.9255E-01 -6.9512E+00
             7.4039E-01
 GRADIENT:  -5.8636E+00  1.1744E+00  3.9881E+00 -8.4404E+00 -5.9756E+00  1.4607E-01 -1.0051E+00  1.2398E+00 -6.7574E-01  0.0000E+00
             1.8375E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1470.05206890325        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      852
 NPARAMETR:  9.9637E-01  1.6705E+00  1.5644E-01  5.1494E-01  7.5119E-01  9.7470E-01  8.4269E-01  2.3643E+00  1.1353E+00  1.0000E-02
             1.9123E+00
 PARAMETER:  9.6360E-02  6.1314E-01 -1.7551E+00 -5.6371E-01 -1.8610E-01  7.4371E-02 -7.1154E-02  9.6049E-01  2.2690E-01 -6.6036E+00
             7.4831E-01
 GRADIENT:   2.5787E-01 -2.4206E+00  1.7554E-01 -1.4108E+00  4.4311E-02  5.9131E-02  1.3830E-01 -8.0645E-02 -2.4059E-01  0.0000E+00
             3.9295E-01

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1470.05206890325        NO. OF FUNC. EVALS.:  32
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.9635E-01  1.6705E+00  1.5645E-01  5.1492E-01  7.5119E-01  9.7467E-01  8.4261E-01  2.3642E+00  1.1353E+00  1.0000E-02
             1.9124E+00
 PARAMETER:  9.6360E-02  6.1314E-01 -1.7551E+00 -5.6371E-01 -1.8610E-01  7.4371E-02 -7.1154E-02  9.6049E-01  2.2690E-01 -6.6036E+00
             7.4831E-01
 GRADIENT:   1.9761E-01  2.1277E+04 -7.4268E+03  2.3141E+04 -7.0121E+04  4.9782E-02  1.3203E-01  1.3564E+04 -2.8757E+04  0.0000E+00
            -1.7442E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      884
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5443E-04 -1.0756E-02 -8.2740E-03  1.3260E-02 -5.0514E-04
 SE:             2.9471E-02  2.6681E-02  1.7560E-02  2.0520E-02  2.7144E-04
 N:                     100         100         100         100         100

 P VAL.:         9.9311E-01  6.8685E-01  6.3751E-01  5.1816E-01  6.2752E-02

 ETASHRINKSD(%)  1.2675E+00  1.0615E+01  4.1171E+01  3.1256E+01  9.9091E+01
 ETASHRINKVR(%)  2.5189E+00  2.0102E+01  6.5391E+01  5.2742E+01  9.9992E+01
 EBVSHRINKSD(%)  1.4538E+00  1.0656E+01  4.1137E+01  3.1728E+01  9.9088E+01
 EBVSHRINKVR(%)  2.8864E+00  2.0177E+01  6.5351E+01  5.3389E+01  9.9992E+01
 RELATIVEINF(%)  9.5655E+01  1.3054E+01  1.3408E+01  7.7348E+00  1.3287E-03
 EPSSHRINKSD(%)  3.9268E+01
 EPSSHRINKVR(%)  6.3117E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1470.0520689032462     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -734.90124233950803     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.29
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1470.052       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  1.67E+00  1.56E-01  5.15E-01  7.51E-01  9.75E-01  8.43E-01  2.36E+00  1.14E+00  1.00E-02  1.91E+00
 


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
+        3.29E+08
 
 TH 2
+        3.27E+02  3.11E+06
 
 TH 3
+       -1.36E+03  2.57E+04  4.32E+07
 
 TH 4
+       -1.13E+08 -5.99E+03  2.18E+04  3.87E+07
 
 TH 5
+       -2.60E+03 -5.07E+02 -1.37E+03  4.63E+04  1.67E+08
 
 TH 6
+       -8.93E-01  2.61E+02 -9.90E+02 -1.15E+08 -1.95E+03  1.93E+02
 
 TH 7
+        4.52E+00  5.13E+02 -1.90E+03  1.33E+08 -3.66E+03  1.60E+00  1.79E+02
 
 TH 8
+        1.58E+02 -4.64E+03  1.14E+04 -2.80E+03 -2.22E+01  1.21E+02  2.24E+02  6.31E+05
 
 TH 9
+        1.27E+08  9.89E+00 -7.43E+01 -4.36E+07 -1.58E+02 -1.06E+03 -1.98E+03  1.15E+01  4.92E+07
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.62E+02 -1.10E+03 -1.82E+04  4.46E+03 -4.34E+01 -1.89E+02 -3.49E+02 -4.92E+02 -3.70E-01  0.00E+00  1.59E+06
 
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
 #CPUT: Total CPU Time in Seconds,       14.988
Stop Time:
Sat Sep 18 09:20:33 CDT 2021
