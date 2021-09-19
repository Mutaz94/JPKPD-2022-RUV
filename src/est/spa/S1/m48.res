Sat Sep 18 11:07:09 CDT 2021
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
$DATA ../../../../data/spa/S1/dat48.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m48.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1646.88760538911        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6775E+01 -4.0224E+01  1.0044E+00 -5.0204E+01  9.4297E-01 -2.8906E+01 -5.2707E+00 -2.0817E+00  7.9606E+00  1.0345E+01
            -3.3462E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1651.74101046142        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8839E-01  1.0825E+00  9.4038E-01  9.7505E-01  1.0073E+00  1.0664E+00  1.1189E+00  1.0687E+00  9.1203E-01  8.5737E-01
             1.0908E+00
 PARAMETER:  8.8322E-02  1.7929E-01  3.8526E-02  7.4733E-02  1.0731E-01  1.6431E-01  2.1236E-01  1.6645E-01  7.9162E-03 -5.3883E-02
             1.8688E-01
 GRADIENT:   4.6500E+00 -8.5540E+00 -4.9927E+00 -7.9588E+00  1.8722E+01  1.7742E-01  1.8873E+00 -6.1692E-02 -3.3427E+00 -1.9245E+00
             3.3331E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1652.29451589991        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0021E+00  1.0913E+00  8.7490E-01  9.7559E-01  9.6233E-01  1.0762E+00  1.1232E+00  9.9820E-01  9.2899E-01  8.0012E-01
             1.0803E+00
 PARAMETER:  1.0211E-01  1.8740E-01 -3.3646E-02  7.5283E-02  6.1597E-02  1.7339E-01  2.1623E-01  9.8196E-02  2.6340E-02 -1.2299E-01
             1.7728E-01
 GRADIENT:   3.1783E+01  6.6714E+00  6.7536E-02  5.8766E+00  4.7810E+00  4.1957E+00  2.9341E+00  1.1592E-01  6.3996E-02 -3.2669E+00
            -7.3047E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1652.50475340247        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      346
 NPARAMETR:  1.0082E+00  1.2164E+00  7.5484E-01  8.9190E-01  9.6124E-01  1.0863E+00  1.0227E+00  8.4954E-01  9.8594E-01  8.2375E-01
             1.0801E+00
 PARAMETER:  1.0818E-01  2.9592E-01 -1.8125E-01 -1.4406E-02  6.0467E-02  1.8281E-01  1.2242E-01 -6.3060E-02  8.5844E-02 -9.3888E-02
             1.7705E-01
 GRADIENT:   2.7814E+00  4.4242E-01  7.8196E-02  9.1089E-01 -8.4938E-01  8.2145E-01 -2.0871E-03  1.4665E-02 -2.5785E-01  4.6250E-01
            -2.0877E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1652.56801964949        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  1.0077E+00  1.3757E+00  6.2466E-01  7.8403E-01  9.7389E-01  1.0850E+00  9.4045E-01  6.9753E-01  1.0666E+00  8.1261E-01
             1.0813E+00
 PARAMETER:  1.0772E-01  4.1896E-01 -3.7054E-01 -1.4331E-01  7.3541E-02  1.8155E-01  3.8604E-02 -2.6021E-01  1.6447E-01 -1.0751E-01
             1.7818E-01
 GRADIENT:   1.8277E-01  1.0883E+00  9.5405E-02  7.8715E-01 -3.3112E-01 -5.0063E-02 -5.7137E-02  2.5578E-02 -1.5274E-01  4.6615E-02
            -5.4665E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1652.57389586986        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      698
 NPARAMETR:  1.0080E+00  1.4396E+00  5.6529E-01  7.3847E-01  9.7557E-01  1.0851E+00  9.1239E-01  5.7874E-01  1.1041E+00  8.0625E-01
             1.0815E+00
 PARAMETER:  1.0793E-01  4.6436E-01 -4.7041E-01 -2.0318E-01  7.5270E-02  1.8170E-01  8.3111E-03 -4.4690E-01  1.9903E-01 -1.1536E-01
             1.7838E-01
 GRADIENT:   1.0561E-01 -3.0552E-01 -9.2752E-02 -2.8744E-01 -7.0471E-02 -1.0664E-01 -2.3055E-02  6.4009E-02  4.5543E-03  4.3970E-02
            -1.4697E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1652.57844052650        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      874
 NPARAMETR:  1.0081E+00  1.4760E+00  5.1804E-01  7.1170E-01  9.6653E-01  1.0853E+00  8.9945E-01  3.8826E-01  1.1237E+00  7.9660E-01
             1.0817E+00
 PARAMETER:  1.0803E-01  4.8932E-01 -5.5771E-01 -2.4009E-01  6.5956E-02  1.8188E-01 -5.9709E-03 -8.4607E-01  2.1663E-01 -1.2740E-01
             1.7849E-01
 GRADIENT:   3.3647E-02  7.2312E-01  3.9132E-03  2.1347E-01 -7.5824E-01 -1.2693E-01  9.8090E-02  3.0434E-02  3.1906E-02  1.4830E-01
             4.0748E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1652.58237734049        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1050
 NPARAMETR:  1.0080E+00  1.5054E+00  4.8990E-01  6.8987E-01  9.6751E-01  1.0863E+00  8.8748E-01  2.1828E-01  1.1433E+00  7.9388E-01
             1.0819E+00
 PARAMETER:  1.0793E-01  5.0907E-01 -6.1356E-01 -2.7125E-01  6.6967E-02  1.8278E-01 -1.9375E-02 -1.4220E+00  2.3392E-01 -1.3083E-01
             1.7873E-01
 GRADIENT:  -2.1254E-01 -1.1991E+00 -5.2031E-01 -1.3103E-01  7.4907E-01  1.9426E-01  5.4766E-02  1.4787E-02  9.8873E-02  8.8195E-02
             1.2280E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1652.58597033856        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1226
 NPARAMETR:  1.0080E+00  1.5076E+00  4.8421E-01  6.8837E-01  9.6419E-01  1.0860E+00  8.8733E-01  8.7400E-02  1.1435E+00  7.9171E-01
             1.0818E+00
 PARAMETER:  1.0801E-01  5.1050E-01 -6.2523E-01 -2.7343E-01  6.3531E-02  1.8251E-01 -1.9539E-02 -2.3373E+00  2.3410E-01 -1.3356E-01
             1.7866E-01
 GRADIENT:  -6.8808E-02  3.2353E-01  2.0699E-02  1.3673E-01 -2.0260E-01  7.7775E-02  6.1187E-03  1.2681E-03  4.4899E-03  1.6218E-02
             2.8115E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1652.58656561594        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1401
 NPARAMETR:  1.0081E+00  1.5096E+00  4.8264E-01  6.8681E-01  9.6463E-01  1.0858E+00  8.8654E-01  1.5660E-02  1.1451E+00  7.9194E-01
             1.0818E+00
 PARAMETER:  1.0805E-01  5.1184E-01 -6.2848E-01 -2.7570E-01  6.3986E-02  1.8233E-01 -2.0431E-02 -4.0566E+00  2.3550E-01 -1.3327E-01
             1.7862E-01
 GRADIENT:   8.2900E-04  1.1417E-02 -2.0739E-03  5.3351E-03 -2.2938E-02  7.1133E-03  3.5869E-03  3.9150E-05  7.5121E-03  7.4928E-03
             6.0764E-03

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1652.58657696891        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1500
 NPARAMETR:  1.0081E+00  1.5095E+00  4.8268E-01  6.8674E-01  9.6464E-01  1.0859E+00  8.8660E-01  1.0000E-02  1.1451E+00  7.9192E-01
             1.0818E+00
 PARAMETER:  1.0805E-01  5.1186E-01 -6.2850E-01 -2.7573E-01  6.4019E-02  1.8231E-01 -2.0463E-02 -4.5410E+00  2.3549E-01 -1.3327E-01
             1.7861E-01
 GRADIENT:  -1.6815E-02  2.3349E-02 -4.9297E-03  1.1607E-02  6.1143E-03 -6.6285E-03 -2.8320E-03  0.0000E+00 -7.9261E-04  5.0740E-04
            -1.2689E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1500
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7944E-04 -1.8592E-02 -3.2052E-04  1.5396E-02 -2.7499E-02
 SE:             2.9838E-02  2.4470E-02  1.3485E-04  2.3607E-02  2.0574E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9520E-01  4.4738E-01  1.7457E-02  5.1427E-01  1.8135E-01

 ETASHRINKSD(%)  4.0060E-02  1.8022E+01  9.9548E+01  2.0914E+01  3.1076E+01
 ETASHRINKVR(%)  8.0105E-02  3.2796E+01  9.9998E+01  3.7454E+01  5.2495E+01
 EBVSHRINKSD(%)  4.3072E-01  1.7873E+01  9.9606E+01  2.1593E+01  3.0829E+01
 EBVSHRINKVR(%)  8.5958E-01  3.2552E+01  9.9998E+01  3.8523E+01  5.2154E+01
 RELATIVEINF(%)  9.9121E+01  4.3279E+00  1.2851E-04  3.7958E+00  5.9494E+00
 EPSSHRINKSD(%)  4.3882E+01
 EPSSHRINKVR(%)  6.8507E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1652.5865769689055     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -917.43575040516737     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.58
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.04
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1652.587       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.51E+00  4.83E-01  6.87E-01  9.65E-01  1.09E+00  8.87E-01  1.00E-02  1.15E+00  7.92E-01  1.08E+00
 


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
+        9.20E+02
 
 TH 2
+       -5.47E+00  4.32E+02
 
 TH 3
+        5.49E+00  2.56E+02  8.01E+02
 
 TH 4
+       -1.31E+01  2.92E+02 -5.35E+02  1.12E+03
 
 TH 5
+       -6.00E+00 -3.15E+02 -7.06E+02  4.49E+02  9.14E+02
 
 TH 6
+        1.53E-01 -1.10E+00  2.52E+00 -4.39E+00 -1.37E+00  1.66E+02
 
 TH 7
+       -9.25E-01  1.74E+01 -3.82E+01 -8.65E+00  2.54E+00 -1.23E+00  1.20E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.57E+00 -1.93E+01 -4.30E+01  5.65E+01 -8.18E+00  1.35E-01  1.61E+01  0.00E+00  6.59E+01
 
 TH10
+        3.70E-01 -1.50E+01 -4.50E+01 -1.62E+01 -7.34E+01  3.43E-01  1.90E+01  0.00E+00  1.27E+01  8.79E+01
 
 TH11
+       -7.26E+00 -1.53E+01 -3.00E+01 -1.16E+00 -7.58E+00  1.95E+00  9.68E+00  0.00E+00  8.56E+00  1.89E+01  1.81E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.682
Stop Time:
Sat Sep 18 11:07:35 CDT 2021
