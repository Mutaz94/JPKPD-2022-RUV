Thu Sep 30 06:20:24 CDT 2021
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
$DATA ../../../../data/spa2/A3/dat30.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   475.157503033878        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5320E+02  1.3478E+02  3.1942E+02  2.2166E+01  2.2938E+02  6.0490E+01 -1.3044E+02 -3.3617E+02 -1.7437E+01 -1.1330E+02
            -5.2586E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1626.23514550503        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0221E+00  1.0152E+00  8.7191E-01  1.0942E+00  9.2664E-01  8.5039E-01  1.0092E+00  9.9897E-01  9.0133E-01  8.6500E-01
             5.2819E+00
 PARAMETER:  1.2187E-01  1.1509E-01 -3.7071E-02  1.9006E-01  2.3810E-02 -6.2060E-02  1.0918E-01  9.8970E-02 -3.8884E-03 -4.5030E-02
             1.7643E+00
 GRADIENT:  -6.8945E+01 -2.8576E+01 -2.9483E+01  1.1510E+01  3.1446E+01 -3.3143E+01  5.5782E+00  7.1179E+00  1.1549E+01  2.2823E+01
             4.4229E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1648.58845077070        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.9842E-01  6.4956E-01  3.8251E-01  1.2816E+00  4.5360E-01  9.3704E-01  1.3893E+00  4.9004E-01  1.0500E+00  2.8120E-01
             4.8179E+00
 PARAMETER:  9.8422E-02 -3.3147E-01 -8.6101E-01  3.4812E-01 -6.9054E-01  3.4967E-02  4.2882E-01 -6.1327E-01  1.4877E-01 -1.1687E+00
             1.6723E+00
 GRADIENT:  -1.0328E+02  4.8222E+01 -2.7581E+01  1.2463E+02 -5.8343E+00 -1.1296E+01  1.4228E+01  2.1797E+00  1.7475E+01  1.0619E+00
             3.7888E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1728.38831322412        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  1.0003E+00  5.1979E-01  3.0168E-01  1.2311E+00  3.5603E-01  9.0960E-01  1.3336E+00  7.1998E-02  1.0429E+00  5.1785E-01
             3.4200E+00
 PARAMETER:  1.0032E-01 -5.5434E-01 -1.0984E+00  3.0795E-01 -9.3275E-01  5.2494E-03  3.8786E-01 -2.5311E+00  1.4202E-01 -5.5807E-01
             1.3296E+00
 GRADIENT:  -5.5490E+01  7.9245E+01 -1.3144E+01  9.4354E+01 -8.5709E+01 -2.6109E+01  3.2189E+00 -1.1961E-02 -7.0295E+00 -9.4447E+00
             7.4634E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1744.61975323830        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      428
 NPARAMETR:  1.0181E+00  4.8419E-01  3.8610E-01  1.2144E+00  4.0539E-01  9.6234E-01  1.1659E+00  3.3397E-02  1.0090E+00  7.1916E-01
             3.1627E+00
 PARAMETER:  1.1791E-01 -6.2528E-01 -8.5166E-01  2.9429E-01 -8.0291E-01  6.1612E-02  2.5345E-01 -3.2993E+00  1.0899E-01 -2.2967E-01
             1.2514E+00
 GRADIENT:   2.9412E+00  1.0577E+00 -1.3320E+00  7.1875E+00  3.4888E+00 -4.3198E-01  8.1820E-01  1.0340E-02  1.2948E+00 -6.5579E-01
            -8.2806E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1746.28047703941        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  1.0184E+00  3.4701E-01  2.5165E-01  1.1882E+00  2.7935E-01  9.6396E-01  6.1743E-01  2.3927E-02  1.1654E+00  8.4978E-01
             3.0967E+00
 PARAMETER:  1.1826E-01 -9.5840E-01 -1.2797E+00  2.7248E-01 -1.1753E+00  6.3295E-02 -3.8219E-01 -3.6327E+00  2.5303E-01 -6.2774E-02
             1.2303E+00
 GRADIENT:   7.3409E-01  1.9334E+01  6.0141E+00  6.2400E+01 -3.3099E+01 -2.5008E+00  5.2006E-01  8.7579E-03 -1.0333E+01  1.0861E+01
             1.4972E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1749.54958342380        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      783
 NPARAMETR:  1.0163E+00  3.0705E-01  2.0816E-01  1.0880E+00  2.4689E-01  9.7035E-01  3.8060E-01  2.3304E-02  1.2794E+00  8.3403E-01
             3.0732E+00
 PARAMETER:  1.1614E-01 -1.0807E+00 -1.4694E+00  1.8431E-01 -1.2988E+00  6.9897E-02 -8.6600E-01 -3.6591E+00  3.4642E-01 -8.1484E-02
             1.2227E+00
 GRADIENT:   1.0611E-01  4.1127E+00  7.0024E-01  1.1791E+00 -3.6067E+00 -2.5260E-01  2.7876E-01  5.8713E-03  4.8062E-01  1.9517E+00
             3.6798E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1749.64751524770        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      966
 NPARAMETR:  1.0162E+00  3.0041E-01  2.0353E-01  1.0802E+00  2.4233E-01  9.7132E-01  2.4593E-01  1.0331E-02  1.2867E+00  8.2934E-01
             3.0627E+00
 PARAMETER:  1.1609E-01 -1.1026E+00 -1.4919E+00  1.7715E-01 -1.3174E+00  7.0897E-02 -1.3027E+00 -4.4726E+00  3.5210E-01 -8.7124E-02
             1.2193E+00
 GRADIENT:   2.0488E-01  3.0321E+00  2.2364E+00  4.9650E-01 -5.7684E+00  1.4879E-02  4.0039E-02  1.1207E-03 -6.6911E-02 -1.9370E+00
            -2.1542E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1749.69108776842        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1143
 NPARAMETR:  1.0163E+00  2.9852E-01  2.0320E-01  1.0788E+00  2.4220E-01  9.7103E-01  3.9491E-02  1.0000E-02  1.2872E+00  8.4103E-01
             3.0699E+00
 PARAMETER:  1.1621E-01 -1.1089E+00 -1.4935E+00  1.7586E-01 -1.3180E+00  7.0597E-02 -3.1317E+00 -7.1768E+00  3.5246E-01 -7.3132E-02
             1.2217E+00
 GRADIENT:   1.7045E-01  2.4918E-03  4.1351E-01 -1.2923E-02  3.0650E-01  3.8942E-02  1.3125E-03  0.0000E+00  1.0697E-01  2.1042E-01
             2.5462E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1749.69754858104        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1320
 NPARAMETR:  1.0163E+00  2.9675E-01  2.0119E-01  1.0760E+00  2.4058E-01  9.7109E-01  1.0000E-02  1.0000E-02  1.2904E+00  8.4056E-01
             3.0687E+00
 PARAMETER:  1.1616E-01 -1.1149E+00 -1.5035E+00  1.7321E-01 -1.3247E+00  7.0669E-02 -1.2663E+01 -1.9068E+01  3.5498E-01 -7.3681E-02
             1.2212E+00
 GRADIENT:   5.2869E-02 -1.6207E-01 -1.2633E-01 -1.9618E-04  2.4124E-01 -1.7521E-02  0.0000E+00  0.0000E+00 -4.9168E-02  4.0259E-02
             1.5549E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1749.69754858104        NO. OF FUNC. EVALS.:  31
 CUMULATIVE NO. OF FUNC. EVALS.:     1351
 NPARAMETR:  1.0165E+00  2.9780E-01  2.0169E-01  1.0755E+00  2.3991E-01  9.7154E-01  1.0000E-02  1.0000E-02  1.2941E+00  8.3972E-01
             3.0674E+00
 PARAMETER:  1.1616E-01 -1.1149E+00 -1.5035E+00  1.7321E-01 -1.3247E+00  7.0669E-02 -1.2663E+01 -1.9068E+01  3.5498E-01 -7.3681E-02
             1.2212E+00
 GRADIENT:  -6.0658E-02 -3.3457E-01 -2.8585E-01  7.1181E-02  1.0642E+00 -2.3962E-02  0.0000E+00  0.0000E+00 -1.0415E-01  3.9637E-02
             7.7109E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1351
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2590E-03 -5.0058E-05  1.8597E-04 -9.3249E-03  1.5359E-03
 SE:             2.8977E-02  1.6355E-04  1.7872E-04  2.6631E-02  2.7159E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6535E-01  7.5955E-01  2.9806E-01  7.2622E-01  9.5490E-01

 ETASHRINKSD(%)  2.9221E+00  9.9452E+01  9.9401E+01  1.0783E+01  9.0137E+00
 ETASHRINKVR(%)  5.7587E+00  9.9997E+01  9.9996E+01  2.0403E+01  1.7215E+01
 EBVSHRINKSD(%)  2.8548E+00  9.9454E+01  9.9373E+01  8.2981E+00  9.4815E+00
 EBVSHRINKVR(%)  5.6280E+00  9.9997E+01  9.9996E+01  1.5908E+01  1.8064E+01
 RELATIVEINF(%)  9.4211E+01  5.2640E-04  3.7752E-04  3.8821E+01  4.6405E+00
 EPSSHRINKSD(%)  2.2902E+01
 EPSSHRINKVR(%)  4.0559E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1749.6975485810424     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -646.97130873543529     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.91
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1749.698       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.97E-01  2.01E-01  1.08E+00  2.41E-01  9.71E-01  1.00E-02  1.00E-02  1.29E+00  8.41E-01  3.07E+00
 


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
+        1.08E+03
 
 TH 2
+        3.18E+01  3.66E+03
 
 TH 3
+       -5.60E+00  1.99E+03  9.64E+03
 
 TH 4
+       -1.94E+01  2.72E+01 -5.67E+02  4.84E+02
 
 TH 5
+        4.32E+01 -6.38E+03 -1.19E+04 -1.98E+02  2.23E+04
 
 TH 6
+        2.13E+00 -1.43E+01  5.11E+01 -9.73E+00 -3.25E+00  1.89E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.08E+01 -2.96E+01  1.60E+02 -6.46E+00  2.85E+01  1.31E+00  0.00E+00  0.00E+00  7.46E+01
 
 TH10
+       -6.94E-01 -2.63E+01 -9.39E+00  1.02E+01  1.76E+01  1.47E+00  0.00E+00  0.00E+00  4.12E+00  1.90E+02
 
 TH11
+       -1.95E+01 -1.21E+01 -4.18E+01 -7.12E+00  5.43E+01  1.71E+00  0.00E+00  0.00E+00  5.24E+00  8.39E+00  6.86E+01
 
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
 #CPUT: Total CPU Time in Seconds,       35.723
Stop Time:
Thu Sep 30 06:21:01 CDT 2021
