Wed Sep 29 13:26:00 CDT 2021
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
$DATA ../../../../data/spa/A3/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   311.324451590530        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8566E+02  7.2634E+01  3.4565E+01  3.4497E+01  2.1234E+02  2.1600E+01 -5.9759E+01 -1.3118E+01 -1.5825E+02 -1.2138E+02
            -3.5108E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -971.184097565185        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1622E+00  1.0481E+00  9.4769E-01  1.3834E+00  7.2693E-01  7.5462E-01  1.1049E+00  9.7520E-01  1.4769E+00  1.0928E+00
             1.0119E+01
 PARAMETER:  2.5032E-01  1.4699E-01  4.6274E-02  4.2452E-01 -2.1892E-01 -1.8154E-01  1.9980E-01  7.4889E-02  4.8998E-01  1.8875E-01
             2.4144E+00
 GRADIENT:  -2.9138E+01  2.4328E+01  6.7490E+00  5.0001E+01 -4.2194E+01 -1.6216E+01  7.9592E+00  4.1677E+00  3.4806E+01  1.8413E+01
             3.7280E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1038.92776954631        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0407E+00  4.9426E-01  1.1316E-01  1.4155E+00  2.0242E-01  8.6946E-01  1.3351E+00  2.4051E-01  1.6385E+00  3.6375E-01
             7.0139E+00
 PARAMETER:  1.3988E-01 -6.0469E-01 -2.0789E+00  4.4745E-01 -1.4974E+00 -3.9882E-02  3.8898E-01 -1.3250E+00  5.9378E-01 -9.1128E-01
             2.0479E+00
 GRADIENT:  -6.2327E+01  5.2889E+01 -2.9802E+01  1.5961E+02 -2.6540E+01 -4.2284E+01  2.3089E+01  3.6811E-01 -1.6056E+01  4.7777E+00
             2.3494E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1201.61703883505        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.1602E-01  7.6259E-01  1.8503E-01  9.3104E-01  3.5524E-01  1.0521E+00  6.2243E-01  2.9271E-02  1.3913E+00  3.4892E-01
             4.0411E+00
 PARAMETER:  1.2283E-02 -1.7104E-01 -1.5873E+00  2.8547E-02 -9.3495E-01  1.5079E-01 -3.7412E-01 -3.4312E+00  4.3026E-01 -9.5291E-01
             1.4965E+00
 GRADIENT:  -1.2784E+02 -2.3913E+01 -1.0303E+01 -2.8933E-02  6.4016E+01  1.5264E+01 -3.3288E+00  8.4665E-03  1.0909E+01  3.5080E+00
             7.7549E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1217.51477278995        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      387
 NPARAMETR:  1.0019E+00  6.0331E-01  1.8117E-01  9.9577E-01  2.9250E-01  9.8146E-01  1.0120E+00  5.1743E-02  1.4009E+00  2.1560E-01
             3.5251E+00
 PARAMETER:  1.0187E-01 -4.0533E-01 -1.6083E+00  9.5757E-02 -1.1293E+00  8.1290E-02  1.1190E-01 -2.8615E+00  4.3709E-01 -1.4343E+00
             1.3599E+00
 GRADIENT:   1.6113E+01  1.2660E+01 -2.5831E+00  5.1038E+00 -4.4790E+00 -6.1762E+00 -4.0571E+00  9.9495E-03  1.0510E+01 -3.0523E-01
            -8.1852E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1219.83343215088        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      562
 NPARAMETR:  1.0036E+00  4.9985E-01  2.1720E-01  1.0631E+00  2.8947E-01  9.8669E-01  1.2724E+00  8.0013E-02  1.2385E+00  1.7644E-01
             3.6049E+00
 PARAMETER:  1.0357E-01 -5.9346E-01 -1.4269E+00  1.6117E-01 -1.1397E+00  8.6602E-02  3.4088E-01 -2.4256E+00  3.1391E-01 -1.6348E+00
             1.3823E+00
 GRADIENT:  -1.7910E+00  8.3579E-01 -2.7307E+00 -1.8828E+00  4.8674E+00 -2.3050E+00  1.2362E-01  1.8085E-02 -1.8659E+00 -6.4048E-01
            -2.2749E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1222.43382923471        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      741
 NPARAMETR:  1.0089E+00  3.5143E-01  2.5551E-01  1.1749E+00  2.7414E-01  9.8829E-01  1.4307E+00  1.0712E-02  1.2100E+00  4.4369E-01
             3.5791E+00
 PARAMETER:  1.0886E-01 -9.4576E-01 -1.2645E+00  2.6116E-01 -1.1941E+00  8.8218E-02  4.5817E-01 -4.4364E+00  2.9059E-01 -7.1263E-01
             1.3751E+00
 GRADIENT:  -5.7856E+00  1.2296E+01  6.7781E+00  1.5745E+01 -2.0931E+01 -1.0733E+00  9.2163E-01  9.3130E-04  4.3265E-01 -1.3669E+00
            -4.8546E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1224.48834029702        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      916
 NPARAMETR:  1.0077E+00  2.0032E-01  2.4716E-01  1.1822E+00  2.4646E-01  9.9864E-01  1.3413E+00  1.0000E-02  1.2307E+00  5.1011E-01
             3.6068E+00
 PARAMETER:  1.0770E-01 -1.5078E+00 -1.2977E+00  2.6742E-01 -1.3006E+00  9.8642E-02  3.9362E-01 -7.9981E+00  3.0759E-01 -5.7314E-01
             1.3828E+00
 GRADIENT:  -2.2732E+00  2.5108E+00  1.3813E-01 -7.6032E+00 -4.0226E+00  4.0342E+00 -1.8353E+00  0.0000E+00  2.1372E+00 -1.1763E+00
             1.4988E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1228.21224614804        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1093
 NPARAMETR:  9.9124E-01  3.9643E-02  2.5935E-01  1.2455E+00  2.4151E-01  9.9363E-01  1.4920E+00  1.0000E-02  1.1626E+00  5.1196E-01
             3.5804E+00
 PARAMETER:  9.1206E-02 -3.1278E+00 -1.2496E+00  3.1957E-01 -1.3209E+00  9.3605E-02  5.0013E-01 -1.7081E+01  2.5064E-01 -5.6951E-01
             1.3755E+00
 GRADIENT:  -9.5546E+00  1.0820E+00  1.6797E+00 -3.8213E+00 -1.2009E+01  4.3853E+00 -6.9877E-02  0.0000E+00 -4.4514E+00 -6.6125E-01
            -2.7493E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1229.96325244413        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1269            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9107E-01  1.0000E-02  3.6171E-01  1.3831E+00  3.0457E-01  9.7214E-01  1.6478E+00  1.0000E-02  1.0689E+00  4.2651E-01
             3.7027E+00
 PARAMETER:  9.1033E-02 -4.6920E+00 -9.1692E-01  4.2436E-01 -1.0889E+00  7.1749E-02  5.9944E-01 -2.5942E+01  1.6658E-01 -7.5212E-01
             1.4091E+00
 GRADIENT:   2.4631E+01  0.0000E+00  6.6228E+00  4.0301E+01  3.1315E+01  2.1861E+00 -3.6267E-03  0.0000E+00  2.8319E+00  4.7090E-01
             1.2514E+01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1229.96516854191        NO. OF FUNC. EVALS.:  90
 CUMULATIVE NO. OF FUNC. EVALS.:     1359
 NPARAMETR:  9.9115E-01  1.0000E-02  3.6202E-01  1.3841E+00  3.0441E-01  9.7194E-01  1.6506E+00  1.0000E-02  1.0680E+00  4.2331E-01
             3.6995E+00
 PARAMETER:  9.1115E-02 -4.6920E+00 -9.1606E-01  4.2503E-01 -1.0894E+00  7.1543E-02  6.0117E-01 -2.5942E+01  1.6576E-01 -7.5964E-01
             1.4082E+00
 GRADIENT:   1.2779E-01  0.0000E+00  3.1021E-01  1.0689E-01 -8.1404E-01 -2.2726E-02 -4.1994E-03  0.0000E+00 -5.7402E-02 -9.9346E-02
             1.3737E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1359
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.3595E-04 -5.5319E-04  1.2879E-04 -1.4917E-02  4.3266E-04
 SE:             2.8366E-02  2.6743E-04  2.0532E-04  2.6281E-02  1.3057E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7649E-01  3.8587E-02  5.3051E-01  5.7030E-01  9.7357E-01

 ETASHRINKSD(%)  4.9688E+00  9.9104E+01  9.9312E+01  1.1955E+01  5.6257E+01
 ETASHRINKVR(%)  9.6907E+00  9.9992E+01  9.9995E+01  2.2481E+01  8.0865E+01
 EBVSHRINKSD(%)  4.6477E+00  9.9285E+01  9.9263E+01  1.1253E+01  5.6204E+01
 EBVSHRINKVR(%)  9.0794E+00  9.9995E+01  9.9995E+01  2.1239E+01  8.0819E+01
 RELATIVEINF(%)  7.4225E+01  4.8343E-04  1.6983E-04  1.9448E+01  4.0626E-01
 EPSSHRINKSD(%)  2.4722E+01
 EPSSHRINKVR(%)  4.3332E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1229.9651685419124     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -494.81434197817418     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.72
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1229.965       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  1.00E-02  3.62E-01  1.38E+00  3.04E-01  9.72E-01  1.65E+00  1.00E-02  1.07E+00  4.23E-01  3.70E+00
 


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
+        1.09E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -8.32E+01  0.00E+00  4.29E+03
 
 TH 4
+       -5.04E+01  0.00E+00 -1.72E+02  3.99E+02
 
 TH 5
+        2.70E+02  0.00E+00 -6.81E+03 -3.46E+02  1.22E+04
 
 TH 6
+       -1.99E+00  0.00E+00  2.67E+01 -1.32E+01 -3.15E+00  1.71E+02
 
 TH 7
+       -5.34E-03  0.00E+00 -1.24E-02 -2.73E-03  3.99E-02 -1.73E-02 -4.28E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.08E+01  0.00E+00  2.23E+01 -1.07E+01  1.02E+02  6.17E-02  1.03E-02  0.00E+00  1.03E+02
 
 TH10
+       -1.17E+01  0.00E+00 -1.09E+02 -9.35E+00  2.32E+02  2.99E+00  6.44E-03  0.00E+00  7.99E+00  4.32E+01
 
 TH11
+       -1.84E+01  0.00E+00 -4.46E+00 -5.69E+00  1.23E+00  3.94E+00  1.20E-03  0.00E+00  6.42E+00  1.56E+01  2.75E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.936
Stop Time:
Wed Sep 29 13:26:26 CDT 2021
