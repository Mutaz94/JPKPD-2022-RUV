Sat Sep 18 11:07:56 CDT 2021
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1696.61860074373        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2273E+01 -3.7560E+01 -3.8536E+01 -4.1565E+01  1.9181E-01  4.1386E+00  1.5087E+01  1.9092E+01 -2.2831E+00  4.1362E+01
             2.5538E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1707.41573486331        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0101E+00  9.5768E-01  1.2057E+00  1.0269E+00  1.0579E+00  9.8323E-01  8.0979E-01  8.3858E-01  1.0577E+00  7.4290E-01
             1.0586E+00
 PARAMETER:  1.1009E-01  5.6756E-02  2.8706E-01  1.2652E-01  1.5629E-01  8.3092E-02 -1.1098E-01 -7.6049E-02  1.5610E-01 -1.9720E-01
             1.5696E-01
 GRADIENT:   7.8159E+01 -3.2199E+01  5.5947E+00 -4.7791E+01  2.7690E+01 -2.9851E+00  7.3489E+00  2.6634E+00  3.4907E+00 -1.0187E+01
             1.2337E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1709.06128514256        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0100E+00  8.9609E-01  1.1509E+00  1.0886E+00  1.0429E+00  9.9684E-01  4.9882E-01  4.4650E-01  1.0305E+00  9.2014E-01
             9.5980E-01
 PARAMETER:  1.0992E-01 -9.7198E-03  2.4051E-01  1.8488E-01  1.4204E-01  9.6833E-02 -5.9551E-01 -7.0631E-01  1.3002E-01  1.6771E-02
             5.8966E-02
 GRADIENT:   8.6629E+01 -1.6974E+01 -2.2408E+01  2.6297E+00  5.2564E+01  3.1790E+00 -1.0394E+00  4.6076E-01 -7.1646E+00  5.0204E+00
            -3.0948E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1711.97939534896        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  9.8318E-01  8.3662E-01  1.0649E+00  1.1306E+00  9.4932E-01  9.8961E-01  7.0016E-01  3.5967E-01  1.0024E+00  8.2153E-01
             1.0191E+00
 PARAMETER:  8.3034E-02 -7.8389E-02  1.6290E-01  2.2278E-01  4.7990E-02  8.9559E-02 -2.5644E-01 -9.2257E-01  1.0237E-01 -9.6587E-02
             1.1890E-01
 GRADIENT:  -2.4432E+01 -5.0972E+00 -1.0020E+01 -7.2408E-01  1.3237E+01 -3.7873E+00  7.1291E-01  6.4501E-01  1.6836E+00  3.4930E+00
            -1.4709E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1712.54181644509        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      430
 NPARAMETR:  9.9065E-01  6.6569E-01  1.0818E+00  1.2411E+00  8.8615E-01  9.9591E-01  7.3262E-01  1.3689E-01  9.2186E-01  8.1692E-01
             1.0294E+00
 PARAMETER:  9.0607E-02 -3.0693E-01  1.7865E-01  3.1603E-01 -2.0864E-02  9.5901E-02 -2.1113E-01 -1.8886E+00  1.8634E-02 -1.0222E-01
             1.2896E-01
 GRADIENT:  -3.8559E+00  4.1766E+00  1.6350E+00  8.0413E+00 -4.7970E+00 -2.1802E-01 -8.3124E-02  9.5791E-02 -4.4424E-01  1.3061E+00
             5.5232E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1712.63172642669        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  9.9072E-01  5.4833E-01  1.0789E+00  1.3079E+00  8.4691E-01  9.9481E-01  8.1314E-01  3.2892E-02  8.7652E-01  8.0921E-01
             1.0263E+00
 PARAMETER:  9.0673E-02 -5.0087E-01  1.7590E-01  3.6839E-01 -6.6165E-02  9.4795E-02 -1.0685E-01 -3.3145E+00 -3.1791E-02 -1.1170E-01
             1.2598E-01
 GRADIENT:  -5.6387E-01  4.9964E-01  1.9289E-01  1.3272E+00 -6.6956E-01 -6.4192E-02 -2.8419E-02  5.2940E-03 -9.9801E-02  1.4675E-01
            -5.1645E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1712.63330109153        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      780
 NPARAMETR:  9.9062E-01  5.2468E-01  1.0779E+00  1.3212E+00  8.3915E-01  9.9461E-01  8.3942E-01  2.2612E-02  8.6775E-01  8.0752E-01
             1.0259E+00
 PARAMETER:  9.0580E-02 -5.4496E-01  1.7497E-01  3.7858E-01 -7.5370E-02  9.4598E-02 -7.5038E-02 -3.6893E+00 -4.1848E-02 -1.1379E-01
             1.2561E-01
 GRADIENT:  -5.8062E-02 -9.3317E-03 -3.6176E-02  3.6617E-02  2.1365E-02 -1.2231E-02  1.3909E-02  2.4898E-03  8.2046E-03  4.2013E-02
             1.5531E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1712.63332188447        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      955
 NPARAMETR:  9.9063E-01  5.2356E-01  1.0781E+00  1.3219E+00  8.3892E-01  9.9463E-01  8.3611E-01  2.2091E-02  8.6753E-01  8.0758E-01
             1.0259E+00
 PARAMETER:  9.0590E-02 -5.4710E-01  1.7520E-01  3.7905E-01 -7.5645E-02  9.4618E-02 -7.8994E-02 -3.7126E+00 -4.2104E-02 -1.1371E-01
             1.2561E-01
 GRADIENT:   3.5871E-03 -3.6197E-02 -1.4733E-02 -4.4356E-02  1.5920E-02  2.8604E-03 -4.4557E-03  2.3628E-03  1.0366E-03 -3.9472E-03
             6.9789E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1712.63477041673        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1135
 NPARAMETR:  9.9082E-01  5.3808E-01  1.0787E+00  1.3133E+00  8.4380E-01  9.9482E-01  8.2544E-01  1.0000E-02  8.7282E-01  8.0837E-01
             1.0263E+00
 PARAMETER:  9.0779E-02 -5.1975E-01  1.7572E-01  3.7253E-01 -6.9840E-02  9.4807E-02 -9.1840E-02 -4.6682E+00 -3.6028E-02 -1.1273E-01
             1.2601E-01
 GRADIENT:  -1.1106E-02 -3.3089E-02 -6.1521E-02 -6.4418E-02  8.1306E-02 -8.5434E-04  1.3091E-03  0.0000E+00  9.6013E-03  2.6909E-02
             3.0621E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1712.63477432836        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1319
 NPARAMETR:  9.9083E-01  5.3802E-01  1.0787E+00  1.3133E+00  8.4376E-01  9.9483E-01  8.2577E-01  1.0000E-02  8.7277E-01  8.0826E-01
             1.0263E+00
 PARAMETER:  9.0787E-02 -5.1986E-01  1.7573E-01  3.7257E-01 -6.9892E-02  9.4812E-02 -9.1445E-02 -4.6683E+00 -3.6083E-02 -1.1287E-01
             1.2597E-01
 GRADIENT:   9.1618E-03 -8.2073E-03  3.7666E-03 -2.9117E-02  7.8708E-03  9.5511E-04  1.1350E-03  0.0000E+00  4.8121E-03  2.3313E-03
             1.7189E-03

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1712.63477432836        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1348
 NPARAMETR:  9.9084E-01  5.3809E-01  1.0787E+00  1.3133E+00  8.4374E-01  9.9483E-01  8.2568E-01  1.0000E-02  8.7273E-01  8.0822E-01
             1.0263E+00
 PARAMETER:  9.0787E-02 -5.1986E-01  1.7573E-01  3.7257E-01 -6.9892E-02  9.4812E-02 -9.1445E-02 -4.6683E+00 -3.6083E-02 -1.1287E-01
             1.2597E-01
 GRADIENT:  -5.1424E-03 -8.7427E-03  4.2349E-03  6.0668E-03  8.2944E-03 -1.5376E-04  1.5288E-03  0.0000E+00  4.6568E-03  2.5618E-03
             1.6921E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1348
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.7380E-04 -1.0998E-02 -2.1660E-04 -2.1172E-03 -2.0108E-02
 SE:             2.9822E-02  8.5598E-03  2.0920E-04  2.8316E-02  2.4772E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9535E-01  1.9883E-01  3.0050E-01  9.4040E-01  4.1695E-01

 ETASHRINKSD(%)  9.3701E-02  7.1323E+01  9.9299E+01  5.1371E+00  1.7012E+01
 ETASHRINKVR(%)  1.8731E-01  9.1777E+01  9.9995E+01  1.0010E+01  3.1130E+01
 EBVSHRINKSD(%)  4.4301E-01  7.1887E+01  9.9288E+01  5.0451E+00  1.5739E+01
 EBVSHRINKVR(%)  8.8405E-01  9.2097E+01  9.9995E+01  9.8356E+00  2.9001E+01
 RELATIVEINF(%)  9.7676E+01  2.5675E-01  6.7555E-04  4.7003E+00  3.5595E+00
 EPSSHRINKSD(%)  4.1526E+01
 EPSSHRINKVR(%)  6.5807E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1712.6347743283604     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -977.48394776462226     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.86
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1712.635       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  5.38E-01  1.08E+00  1.31E+00  8.44E-01  9.95E-01  8.26E-01  1.00E-02  8.73E-01  8.08E-01  1.03E+00
 


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
+       -2.08E+01  4.58E+02
 
 TH 3
+        1.14E+01  2.33E+02  4.93E+02
 
 TH 4
+       -9.78E+00  4.48E+02 -5.21E+01  7.63E+02
 
 TH 5
+       -7.43E+00 -5.35E+02 -8.33E+02 -3.94E+00  1.71E+03
 
 TH 6
+        5.16E+00 -3.16E+00  2.37E+00 -2.48E+00  1.81E+00  2.07E+02
 
 TH 7
+       -1.07E+00 -4.22E+00  2.65E+00 -4.25E+00 -4.76E+00 -1.76E+00  9.60E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.66E-01 -6.05E+01  4.20E+00  1.07E+01 -5.05E+00  1.59E+00  1.00E+01  0.00E+00  2.22E+02
 
 TH10
+        2.98E+00  2.22E+01 -2.22E+01 -4.99E+00 -8.58E+01  3.77E-03  7.47E+00  0.00E+00 -1.21E+00  1.55E+02
 
 TH11
+       -6.74E+00 -1.76E+01 -3.93E+01 -8.44E+00  1.63E+01  1.19E-02  2.76E+00  0.00E+00  9.79E+00  3.67E+01  2.21E+02
 
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
 #CPUT: Total CPU Time in Seconds,       20.171
Stop Time:
Sat Sep 18 11:08:18 CDT 2021
