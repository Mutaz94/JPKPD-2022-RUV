Sat Sep 18 12:52:26 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat46.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1607.06331667885        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0733E+02 -1.0425E+02 -7.1293E+01 -4.6565E+01  1.4634E+02  2.0764E+01  4.0391E+00  7.5834E+00  1.8546E+01 -2.0989E+01
            -9.0730E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1627.26995290611        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.5106E-01  1.0671E+00  1.0692E+00  1.0267E+00  9.2917E-01  9.2483E-01  9.4612E-01  9.6816E-01  9.1782E-01  1.0039E+00
             1.2256E+00
 PARAMETER:  4.9820E-02  1.6494E-01  1.6688E-01  1.2638E-01  2.6538E-02  2.1855E-02  4.4614E-02  6.7644E-02  1.4249E-02  1.0389E-01
             3.0340E-01
 GRADIENT:  -2.8619E+01  2.3016E+01  4.3903E+00  1.9816E+01 -1.8501E+01 -7.5363E+00  5.0865E+00  2.8484E+00  3.4946E+00  7.5010E-01
             8.9348E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1628.02081228752        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.5938E-01  1.0862E+00  1.0696E+00  1.0179E+00  9.4791E-01  9.2262E-01  7.5002E-01  8.2009E-01  9.7302E-01  1.0531E+00
             1.1904E+00
 PARAMETER:  5.8536E-02  1.8264E-01  1.6732E-01  1.1778E-01  4.6502E-02  1.9464E-02 -1.8766E-01 -9.8345E-02  7.2652E-02  1.5171E-01
             2.7428E-01
 GRADIENT:  -4.7270E+00  2.7271E+01  3.9092E+00  3.0667E+01 -7.8439E+00 -8.0709E+00 -1.1621E+00 -1.9605E+00 -4.4687E-01  7.8777E-02
            -4.3412E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1628.61176418555        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.6021E-01  1.0229E+00  1.1733E+00  1.0417E+00  9.6999E-01  9.4049E-01  7.5007E-01  9.5684E-01  9.5353E-01  1.0745E+00
             1.1989E+00
 PARAMETER:  5.9400E-02  1.2268E-01  2.5980E-01  1.4088E-01  6.9525E-02  3.8643E-02 -1.8759E-01  5.5876E-02  5.2411E-02  1.7182E-01
             2.8140E-01
 GRADIENT:   7.4468E-01 -1.4419E+00 -5.8722E-01 -1.1911E+00  2.1763E+00  2.6810E-01  1.7653E-01 -2.8153E-01 -9.9099E-02  2.1650E-01
            -2.3184E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1628.65133601818        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.5930E-01  9.9839E-01  1.3095E+00  1.0570E+00  1.0049E+00  9.3916E-01  5.2848E-01  1.1038E+00  1.0053E+00  1.1197E+00
             1.2002E+00
 PARAMETER:  5.8444E-02  9.8385E-02  3.6962E-01  1.5546E-01  1.0484E-01  3.7235E-02 -5.3776E-01  1.9877E-01  1.0530E-01  2.1308E-01
             2.8250E-01
 GRADIENT:   6.2723E-02 -4.2026E-01  3.0065E-01 -1.6069E-01 -2.8259E-01 -4.7218E-02  3.0195E-01  7.6245E-02  1.2225E+00  1.9165E-01
             1.2403E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1628.99302715656        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      482
 NPARAMETR:  9.7199E-01  8.7211E-01  1.4506E+00  1.1474E+00  1.0041E+00  9.4668E-01  4.7500E-01  1.1653E+00  9.4742E-01  1.1340E+00
             1.2010E+00
 PARAMETER:  7.1595E-02 -3.6835E-02  4.7199E-01  2.3752E-01  1.0407E-01  4.5210E-02 -6.4443E-01  2.5296E-01  4.5984E-02  2.2576E-01
             2.8320E-01
 GRADIENT:   7.6982E+00  3.3261E+00  6.4449E-01  3.8025E+00 -9.1219E-01  9.1613E-01  1.7761E-02  2.6248E-02 -3.8516E-02 -5.2879E-02
            -4.0814E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1629.20541492891        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      658
 NPARAMETR:  9.6762E-01  6.7158E-01  1.5859E+00  1.2749E+00  9.8020E-01  9.4252E-01  3.1060E-01  1.1944E+00  8.7337E-01  1.1380E+00
             1.2005E+00
 PARAMETER:  6.7082E-02 -2.9812E-01  5.6115E-01  3.4287E-01  8.0006E-02  4.0806E-02 -1.0692E+00  2.7764E-01 -3.5398E-02  2.2923E-01
             2.8275E-01
 GRADIENT:   1.1873E+00 -7.9926E-01  2.5900E-01 -2.0258E+00  3.2557E-01  7.3803E-02  7.9423E-02 -2.7119E-01  1.1163E+00  1.5956E-01
             1.5031E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1629.24232269790        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      833
 NPARAMETR:  9.6623E-01  5.6183E-01  1.6663E+00  1.3484E+00  9.6866E-01  9.4112E-01  2.0434E-01  1.2392E+00  8.2653E-01  1.1357E+00
             1.2003E+00
 PARAMETER:  6.5647E-02 -4.7655E-01  6.1059E-01  3.9894E-01  6.8158E-02  3.9318E-02 -1.4880E+00  3.1445E-01 -9.0524E-02  2.2725E-01
             2.8258E-01
 GRADIENT:   3.0258E-01  1.4613E-01  1.2925E-01  7.2170E-01 -1.2897E-01  4.8108E-03  2.5647E-02 -8.8096E-02  1.1057E-02  7.9249E-03
            -1.9615E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1629.25127915774        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  9.6614E-01  5.6352E-01  1.6647E+00  1.3455E+00  9.6959E-01  9.4085E-01  9.4117E-02  1.2388E+00  8.2868E-01  1.1354E+00
             1.2004E+00
 PARAMETER:  6.5549E-02 -4.7356E-01  6.0964E-01  3.9678E-01  6.9119E-02  3.9029E-02 -2.2632E+00  3.1417E-01 -8.7920E-02  2.2697E-01
             2.8263E-01
 GRADIENT:   5.9291E-02 -8.6532E-01 -2.4130E-01 -2.6851E+00  9.7220E-01 -1.0913E-01  5.5384E-03 -5.7141E-02 -2.4592E-01 -1.8626E-01
             2.7345E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1629.25702254928        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1188
 NPARAMETR:  9.6632E-01  5.8399E-01  1.6534E+00  1.3330E+00  9.7200E-01  9.4131E-01  1.7904E-02  1.2342E+00  8.3864E-01  1.1370E+00
             1.2003E+00
 PARAMETER:  6.5735E-02 -4.3788E-01  6.0281E-01  3.8745E-01  7.1605E-02  3.9513E-02 -3.9227E+00  3.1040E-01 -7.5978E-02  2.2843E-01
             2.8260E-01
 GRADIENT:  -2.3768E-02 -7.3817E-02 -2.9865E-02 -2.4412E-01  1.4307E-01 -3.0117E-02  2.0224E-04 -1.5663E-02 -4.2775E-02 -3.1490E-02
             5.5563E-05

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1629.25710764732        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1356
 NPARAMETR:  9.6634E-01  5.8324E-01  1.6539E+00  1.3336E+00  9.7176E-01  9.4140E-01  1.0000E-02  1.2347E+00  8.3849E-01  1.1371E+00
             1.2003E+00
 PARAMETER:  6.5735E-02 -4.3916E-01  6.0307E-01  3.8790E-01  7.1399E-02  3.9570E-02 -4.5179E+00  3.1078E-01 -7.6258E-02  2.2849E-01
             2.8259E-01
 GRADIENT:  -1.6311E-02 -1.1031E-04 -8.0918E-03  1.9625E-02  1.4550E-02 -5.1564E-03  0.0000E+00 -8.2248E-04 -1.0771E-02  1.6688E-04
             3.2013E-05

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1356
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.4070E-04 -2.9146E-04 -2.7972E-02 -6.3633E-03 -3.3957E-02
 SE:             2.9766E-02  1.1895E-04  1.4118E-02  2.8780E-02  2.2568E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8819E-01  1.4271E-02  4.7555E-02  8.2501E-01  1.3241E-01

 ETASHRINKSD(%)  2.7893E-01  9.9602E+01  5.2703E+01  3.5840E+00  2.4396E+01
 ETASHRINKVR(%)  5.5708E-01  9.9998E+01  7.7630E+01  7.0396E+00  4.2840E+01
 EBVSHRINKSD(%)  6.5489E-01  9.9638E+01  5.7434E+01  3.9901E+00  2.1058E+01
 EBVSHRINKVR(%)  1.3055E+00  9.9999E+01  8.1881E+01  7.8211E+00  3.7682E+01
 RELATIVEINF(%)  9.8088E+01  7.7918E-05  4.8031E+00  6.2934E+00  1.2968E+01
 EPSSHRINKSD(%)  4.2771E+01
 EPSSHRINKVR(%)  6.7249E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1629.2571076473155     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -894.10628108357730     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.26
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
 





 #OBJV:********************************************    -1629.257       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.66E-01  5.83E-01  1.65E+00  1.33E+00  9.72E-01  9.41E-01  1.00E-02  1.23E+00  8.38E-01  1.14E+00  1.20E+00
 


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
+        1.33E+03
 
 TH 2
+       -2.80E+01  4.51E+02
 
 TH 3
+        2.41E+00  3.56E+01  5.42E+01
 
 TH 4
+       -1.90E+01  5.41E+02 -1.21E+01  8.20E+02
 
 TH 5
+        6.18E-01 -1.90E+02 -1.33E+02 -5.20E+01  5.80E+02
 
 TH 6
+        1.30E+00 -3.85E+00  1.43E+00 -5.17E+00  7.97E-01  2.16E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        7.98E-01 -6.08E+00 -1.80E+01 -2.33E+00 -2.95E-01 -5.63E-01  0.00E+00  2.06E+01
 
 TH 9
+        6.34E+00 -9.86E+01  6.23E+00  2.83E+00  1.04E+00  2.27E-01  0.00E+00 -1.83E+00  2.39E+02
 
 TH10
+       -7.81E-01  4.20E+00 -4.31E+00 -1.63E+00 -6.15E+01  4.61E-01  0.00E+00  8.40E+00  2.89E-01  6.63E+01
 
 TH11
+       -8.44E+00 -2.04E+01 -7.63E+00 -1.23E+01  5.32E-02  2.66E+00  0.00E+00  2.49E+00  1.44E+01  1.44E+01  1.52E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.074
Stop Time:
Sat Sep 18 12:52:48 CDT 2021
