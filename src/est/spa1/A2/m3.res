Wed Sep 29 22:59:39 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat3.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1213.08128978522        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8778E+02  9.9914E+00  4.2194E+01  1.4093E+01  1.4496E+02  5.5906E+01 -2.7620E+01 -3.8016E+01  1.2199E+01 -6.8994E+01
            -1.6868E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1751.75375796990        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0887E+00  9.9957E-01  1.0467E+00  1.0744E+00  9.5596E-01  1.0737E+00  8.8715E-01  8.2508E-01  8.8580E-01  7.9140E-01
             2.0782E+00
 PARAMETER:  1.8501E-01  9.9569E-02  1.4561E-01  1.7176E-01  5.4958E-02  1.7112E-01 -1.9747E-02 -9.2274E-02 -2.1260E-02 -1.3396E-01
             8.3149E-01
 GRADIENT:   2.8999E+02  1.1775E+01  7.0861E+00  2.8780E+01  2.2330E+01  4.8394E+01 -2.6218E+00 -2.5022E+00 -1.0735E+01 -5.7651E+00
            -1.4218E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1757.65310006054        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      185
 NPARAMETR:  1.0693E+00  9.8320E-01  4.1812E-01  1.0200E+00  5.8094E-01  1.0418E+00  6.8857E-01  3.0470E-01  1.0346E+00  5.3444E-01
             2.0626E+00
 PARAMETER:  1.6700E-01  8.3058E-02 -7.7198E-01  1.1976E-01 -4.4311E-01  1.4097E-01 -2.7315E-01 -1.0884E+00  1.3400E-01 -5.2653E-01
             8.2399E-01
 GRADIENT:   6.6529E+01  9.7348E+00 -1.9723E+01  3.6830E+01  3.2063E+01  2.1842E+01 -1.9628E+01 -1.1057E+00  1.1742E+01 -1.5835E+01
            -1.1937E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1770.99851747579        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      360
 NPARAMETR:  1.0326E+00  9.6160E-01  5.1153E-01  1.0442E+00  6.3314E-01  9.5489E-01  8.9910E-01  2.1769E-01  8.9182E-01  5.2290E-01
             2.3491E+00
 PARAMETER:  1.3209E-01  6.0841E-02 -5.7035E-01  1.4321E-01 -3.5707E-01  5.3845E-02 -6.3606E-03 -1.4247E+00 -1.4490E-02 -5.4837E-01
             9.5402E-01
 GRADIENT:  -1.0577E+01  7.7863E+00 -4.5491E+00  1.3524E+01  7.6438E+00 -4.2237E+00 -7.2047E+00 -4.7041E-02 -3.0542E+00 -2.2615E+00
            -2.1685E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1774.63503919877        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      536
 NPARAMETR:  1.0191E+00  5.7832E-01  5.3768E-01  1.2604E+00  5.1154E-01  9.4013E-01  1.3296E+00  4.7651E-01  8.3749E-01  6.3257E-01
             2.2086E+00
 PARAMETER:  1.1894E-01 -4.4763E-01 -5.2049E-01  3.3141E-01 -5.7034E-01  3.8263E-02  3.8491E-01 -6.4127E-01 -7.7349E-02 -3.5797E-01
             8.9236E-01
 GRADIENT:  -3.0906E+01  1.7384E+01 -6.7573E+00  5.7522E+01  6.5349E+00 -9.6097E+00  1.1569E-01  3.1293E+00  2.4250E+00 -6.7923E-01
            -1.3137E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1778.45676675823        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      711
 NPARAMETR:  1.0336E+00  4.3949E-01  4.4366E-01  1.2493E+00  4.1669E-01  9.6594E-01  1.4434E+00  1.6919E-01  8.3080E-01  6.8393E-01
             2.1949E+00
 PARAMETER:  1.3305E-01 -7.2214E-01 -7.1270E-01  3.2257E-01 -7.7542E-01  6.5345E-02  4.6699E-01 -1.6767E+00 -8.5369E-02 -2.7990E-01
             8.8612E-01
 GRADIENT:   6.9961E+00  2.8359E+00  9.9366E+00 -1.6281E+01 -1.3046E+01  1.7336E+00  5.4375E-01  3.8059E-01  9.2584E-01  2.5086E+00
            -2.3811E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1779.18636957905        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      888
 NPARAMETR:  1.0232E+00  3.3043E-01  5.6308E-01  1.3651E+00  4.6302E-01  9.5127E-01  1.6437E+00  1.0000E-02  7.9614E-01  7.5070E-01
             2.2669E+00
 PARAMETER:  1.2297E-01 -1.0074E+00 -4.7433E-01  4.1120E-01 -6.6998E-01  5.0039E-02  5.9696E-01 -4.6174E+00 -1.2798E-01 -1.8674E-01
             9.1841E-01
 GRADIENT:  -6.5800E+00  9.6483E+00  1.8876E+01  1.9180E+01 -3.1087E+01 -9.6687E-01 -1.4116E+00  0.0000E+00 -7.9327E-01  3.6270E+00
             1.3957E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1781.92893432938        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  1.0162E+00  1.1633E-01  5.7938E-01  1.4619E+00  4.4520E-01  9.4604E-01  3.0712E+00  1.0000E-02  7.6170E-01  7.6114E-01
             2.2198E+00
 PARAMETER:  1.1611E-01 -2.0513E+00 -4.4579E-01  4.7973E-01 -7.0924E-01  4.4527E-02  1.2221E+00 -1.4724E+01 -1.7220E-01 -1.7294E-01
             8.9742E-01
 GRADIENT:  -2.6475E-01  1.1274E+00  1.0974E+00  8.8278E+00 -2.2932E+00 -3.0312E-01 -6.7290E-01  0.0000E+00 -2.2163E+00 -1.0106E+00
            -2.9210E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1782.50641415014        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  1.0122E+00  2.3286E-02  5.9067E-01  1.5088E+00  4.3717E-01  9.4398E-01  7.7817E+00  1.0000E-02  7.5283E-01  7.7972E-01
             2.2097E+00
 PARAMETER:  1.1211E-01 -3.6599E+00 -4.2650E-01  5.1131E-01 -7.2744E-01  4.2345E-02  2.1518E+00 -3.1966E+01 -1.8391E-01 -1.4882E-01
             8.9286E-01
 GRADIENT:   3.5727E-01  4.7134E-01  6.1460E+00  1.2493E+01 -1.0419E+01  2.3706E-02  2.6255E-01  0.0000E+00 -1.2200E+00 -4.9168E-01
            -5.8526E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1782.66638668495        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1420             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0115E+00  1.0000E-02  5.8496E-01  1.5058E+00  4.3416E-01  9.4338E-01  1.1231E+01  1.0000E-02  7.5459E-01  7.8000E-01
             2.2209E+00
 PARAMETER:  1.1142E-01 -4.5283E+00 -4.3622E-01  5.0935E-01 -7.3434E-01  4.1712E-02  2.5186E+00 -4.1673E+01 -1.8159E-01 -1.4846E-01
             8.9792E-01
 GRADIENT:   8.7470E+01  1.0972E-01  8.6014E+00  1.8746E+02  3.1184E+01  7.5093E+00  1.6289E-01  0.0000E+00  5.2645E+00  3.6240E-01
             8.9467E+00

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1782.66638668495        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1484
 NPARAMETR:  1.0114E+00  1.0000E-02  5.8386E-01  1.5059E+00  4.3443E-01  9.4334E-01  1.1274E+01  1.0000E-02  7.5461E-01  7.8094E-01
             2.2209E+00
 PARAMETER:  1.1142E-01 -4.5283E+00 -4.3622E-01  5.0935E-01 -7.3434E-01  4.1712E-02  2.5186E+00 -4.1673E+01 -1.8159E-01 -1.4846E-01
             8.9792E-01
 GRADIENT:   7.3069E-02  6.9994E-03  1.1258E+00 -1.4882E-02 -6.5570E-01  6.8251E-03 -2.2995E-03  0.0000E+00 -4.8339E-03 -8.7567E-02
             1.0151E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1484
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.9239E-04  1.2937E-03  4.0915E-06 -7.8127E-03 -1.0622E-02
 SE:             2.9425E-02  1.9319E-03  2.1601E-04  2.7882E-02  2.2793E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8394E-01  5.0309E-01  9.8489E-01  7.7932E-01  6.4120E-01

 ETASHRINKSD(%)  1.4214E+00  9.3528E+01  9.9276E+01  6.5917E+00  2.3641E+01
 ETASHRINKVR(%)  2.8226E+00  9.9581E+01  9.9995E+01  1.2749E+01  4.1693E+01
 EBVSHRINKSD(%)  1.5977E+00  9.4135E+01  9.9260E+01  6.3641E+00  2.3526E+01
 EBVSHRINKVR(%)  3.1699E+00  9.9656E+01  9.9995E+01  1.2323E+01  4.1518E+01
 RELATIVEINF(%)  8.4182E+01  1.9541E-02  2.9690E-04  6.6644E+00  3.1436E+00
 EPSSHRINKSD(%)  2.7185E+01
 EPSSHRINKVR(%)  4.6979E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1782.6663866849542     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -863.72785348028151     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.84
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1782.666       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.00E-02  5.85E-01  1.51E+00  4.34E-01  9.43E-01  1.12E+01  1.00E-02  7.55E-01  7.80E-01  2.22E+00
 


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
+        1.18E+03
 
 TH 2
+       -1.85E+01  2.99E+03
 
 TH 3
+       -1.89E+01  7.16E+01  1.75E+03
 
 TH 4
+       -2.51E+01  8.70E+01 -2.43E+02  7.56E+02
 
 TH 5
+        5.37E+01 -2.05E+02 -2.92E+03 -1.83E+02  5.68E+03
 
 TH 6
+        8.86E-01 -1.77E+00  4.71E+00 -6.95E+00  1.38E+00  2.10E+02
 
 TH 7
+        3.49E-03  7.78E-01  6.79E-03 -1.83E-02  1.84E-02  3.18E-03  5.19E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.54E+00 -5.22E+00  4.43E+01 -9.95E+00 -7.22E-01  5.92E-01  2.13E-02  0.00E+00  2.69E+02
 
 TH10
+       -6.09E-01  4.68E+00 -2.37E+01  9.91E+00 -3.67E+01  2.26E+00 -5.47E-03  0.00E+00 -1.57E+00  1.19E+02
 
 TH11
+       -1.44E+01 -1.55E+00 -1.44E+01 -1.43E+01  2.62E+00  1.66E+00  2.16E-03  0.00E+00  1.12E+01  2.34E+01  9.33E+01
 
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
 #CPUT: Total CPU Time in Seconds,       29.332
Stop Time:
Wed Sep 29 23:00:10 CDT 2021
