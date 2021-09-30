Wed Sep 29 23:29:44 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat59.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m59.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -830.059147987675        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3864E+02  5.7699E+01  1.9966E+02  1.2235E+01  1.6565E+02  3.0916E+01 -1.6337E+01 -3.1741E+02 -4.0059E+01 -9.8887E+01
            -2.1112E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1694.60321427564        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      152
 NPARAMETR:  1.0313E+00  8.6432E-01  9.2527E-01  1.1112E+00  9.1968E-01  8.9522E-01  9.2691E-01  9.3319E-01  1.0337E+00  9.6471E-01
             2.1693E+00
 PARAMETER:  1.3083E-01 -4.5809E-02  2.2326E-02  2.0543E-01  1.6274E-02 -1.0687E-02  2.4097E-02  3.0850E-02  1.3313E-01  6.4073E-02
             8.7443E-01
 GRADIENT:  -9.7549E+01 -1.1661E+01 -3.2046E+01  2.6407E+00  8.5568E+01 -7.6612E+01  1.1629E+00  1.0324E+01  9.4886E-01 -1.5292E+01
            -1.7105E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1728.36923033830        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      327
 NPARAMETR:  1.0351E+00  6.6017E-01  4.6627E-01  1.2138E+00  5.2768E-01  9.4731E-01  4.0688E-01  1.9404E-01  1.2001E+00  7.6497E-01
             2.1092E+00
 PARAMETER:  1.3453E-01 -3.1526E-01 -6.6299E-01  2.9374E-01 -5.3926E-01  4.5873E-02 -7.9923E-01 -1.5397E+00  2.8237E-01 -1.6792E-01
             8.4630E-01
 GRADIENT:  -8.6609E+01  1.3729E+01 -6.1514E+01  1.0635E+02  1.2228E+02 -5.3185E+01 -4.6766E+00  2.2517E-01  4.2650E+00 -7.0891E+00
            -1.0109E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1754.84667508361        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      507
 NPARAMETR:  1.0881E+00  5.1376E-01  3.5894E-01  1.1534E+00  3.9019E-01  1.0759E+00  1.0655E+00  3.8095E-02  1.0361E+00  6.6341E-01
             2.2855E+00
 PARAMETER:  1.8439E-01 -5.6600E-01 -9.2461E-01  2.4272E-01 -8.4113E-01  1.7316E-01  1.6344E-01 -3.1677E+00  1.3545E-01 -3.1037E-01
             9.2657E-01
 GRADIENT:   2.1622E+01  1.5870E+01 -7.4585E-01  2.6766E+00 -3.1057E+00  3.0571E+00  6.4417E+00  1.3169E-03 -6.4260E+00  4.8631E+00
            -2.6819E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1757.16557653286        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      682
 NPARAMETR:  1.0728E+00  4.1122E-01  3.4669E-01  1.1781E+00  3.5945E-01  1.0635E+00  7.4789E-01  2.3114E-02  1.0524E+00  6.5121E-01
             2.2971E+00
 PARAMETER:  1.7027E-01 -7.8863E-01 -9.5932E-01  2.6392E-01 -9.2317E-01  1.6154E-01 -1.9051E-01 -3.6673E+00  1.5105E-01 -3.2892E-01
             9.3165E-01
 GRADIENT:  -2.1074E+00  5.3412E-01  9.1696E-04 -7.7798E-01 -1.7822E+00 -6.1382E-01 -2.5982E-01  8.4232E-04 -1.0667E+00  5.0956E-02
             3.9195E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1757.98702133203        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      861
 NPARAMETR:  1.0572E+00  2.1126E-01  4.4305E-01  1.3141E+00  3.9356E-01  1.0525E+00  1.6169E+00  1.0000E-02  9.6022E-01  6.9821E-01
             2.3119E+00
 PARAMETER:  1.5560E-01 -1.4547E+00 -7.1408E-01  3.7313E-01 -8.3252E-01  1.5118E-01  5.8052E-01 -5.7224E+00  5.9408E-02 -2.5923E-01
             9.3808E-01
 GRADIENT:  -9.1293E+00 -1.2973E+00 -2.4920E+01 -6.4121E+00  3.3210E+01  2.7347E-01  1.1206E-01  0.0000E+00  9.1922E-01  2.4972E+00
            -3.5095E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1758.64577549275        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1038
 NPARAMETR:  1.0562E+00  1.0724E-01  4.7397E-01  1.3800E+00  3.9414E-01  1.0470E+00  2.7366E+00  1.0000E-02  9.1986E-01  6.9964E-01
             2.3360E+00
 PARAMETER:  1.5470E-01 -2.1327E+00 -6.4662E-01  4.2211E-01 -8.3105E-01  1.4590E-01  1.1067E+00 -8.2363E+00  1.6462E-02 -2.5720E-01
             9.4844E-01
 GRADIENT:  -3.1200E-01  7.4978E-01  2.2250E+00  4.9328E+00 -5.2793E+00 -1.0694E-02  2.9908E-01  0.0000E+00 -1.7410E-01  1.7822E-01
             7.8732E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1758.69473402809        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1213
 NPARAMETR:  1.0533E+00  6.1228E-02  4.8623E-01  1.4032E+00  3.9661E-01  1.0449E+00  3.4105E+00  1.0000E-02  9.0777E-01  7.0885E-01
             2.3353E+00
 PARAMETER:  1.5193E-01 -2.6931E+00 -6.2108E-01  4.3873E-01 -8.2479E-01  1.4396E-01  1.3268E+00 -1.0491E+01  3.2305E-03 -2.4411E-01
             9.4815E-01
 GRADIENT:  -4.2338E-01  1.6217E-01  3.5573E-01  7.9517E-01 -9.2057E-01 -8.0611E-02  1.4096E-01  0.0000E+00  6.2780E-02  9.6772E-02
            -1.5963E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1758.71372832531        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1401
 NPARAMETR:  1.0537E+00  5.8452E-02  4.8616E-01  1.4024E+00  3.9658E-01  1.0453E+00  2.3357E+00  1.0000E-02  9.0901E-01  7.0945E-01
             2.3357E+00
 PARAMETER:  1.5228E-01 -2.7396E+00 -6.2121E-01  4.3822E-01 -8.2488E-01  1.4432E-01  9.4830E-01 -1.0491E+01  4.6009E-03 -2.4327E-01
             9.4831E-01
 GRADIENT:   6.3998E-01 -2.6711E-02 -6.2771E-01 -2.2496E+00  1.2356E+00  4.6398E-02  1.3146E-03  0.0000E+00  1.0860E-01 -1.7278E-02
             1.0142E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1758.71457709281        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1581
 NPARAMETR:  1.0536E+00  6.0606E-02  4.8611E-01  1.4022E+00  3.9647E-01  1.0453E+00  2.3263E+00  1.0000E-02  9.0880E-01  7.0953E-01
             2.3353E+00
 PARAMETER:  1.5226E-01 -2.7034E+00 -6.2132E-01  4.3807E-01 -8.2516E-01  1.4431E-01  9.4428E-01 -1.0491E+01  4.3735E-03 -2.4315E-01
             9.4814E-01
 GRADIENT:   3.0253E-01  4.6760E-02  6.3220E-01 -8.5958E-01 -9.3257E-01  4.0184E-03  4.3089E-04  0.0000E+00 -1.5798E-01  4.2306E-02
            -6.2101E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1758.71457709281        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1609
 NPARAMETR:  1.0538E+00  5.9809E-02  4.8579E-01  1.4012E+00  3.9666E-01  1.0455E+00  2.3484E+00  1.0000E-02  9.0953E-01  7.0925E-01
             2.3354E+00
 PARAMETER:  1.5226E-01 -2.7034E+00 -6.2132E-01  4.3807E-01 -8.2516E-01  1.4431E-01  9.4428E-01 -1.0491E+01  4.3735E-03 -2.4315E-01
             9.4814E-01
 GRADIENT:  -2.3712E-01  1.8588E-02  6.3727E-01  1.2837E+00 -8.4992E-01 -3.6196E-02 -3.9250E-04  0.0000E+00 -9.0874E-02  4.0115E-02
            -2.5004E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1609
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.9440E-04 -1.7813E-03  2.8936E-05 -6.8988E-03 -7.3977E-03
 SE:             2.9466E-02  2.0637E-03  2.3161E-04  2.8120E-02  2.2126E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8661E-01  3.8806E-01  9.0057E-01  8.0619E-01  7.3812E-01

 ETASHRINKSD(%)  1.2848E+00  9.3086E+01  9.9224E+01  5.7959E+00  2.5874E+01
 ETASHRINKVR(%)  2.5532E+00  9.9522E+01  9.9994E+01  1.1256E+01  4.5054E+01
 EBVSHRINKSD(%)  1.4524E+00  9.3243E+01  9.9213E+01  5.2656E+00  2.5513E+01
 EBVSHRINKVR(%)  2.8838E+00  9.9543E+01  9.9994E+01  1.0254E+01  4.4517E+01
 RELATIVEINF(%)  8.5589E+01  4.1646E-02  4.0063E-04  1.4087E+01  2.8542E+00
 EPSSHRINKSD(%)  2.6900E+01
 EPSSHRINKVR(%)  4.6564E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1758.7145770928103     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -839.77604388813756     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.26
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1758.715       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  6.06E-02  4.86E-01  1.40E+00  3.96E-01  1.05E+00  2.33E+00  1.00E-02  9.09E-01  7.10E-01  2.34E+00
 


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
+        8.90E+02
 
 TH 2
+       -6.77E+01  2.67E+02
 
 TH 3
+       -1.81E+01  4.19E+02  2.84E+03
 
 TH 4
+       -1.74E+01  2.69E+02 -2.53E+02  6.16E+02
 
 TH 5
+        5.72E+01 -9.68E+02 -4.38E+03 -1.70E+02  7.89E+03
 
 TH 6
+        5.97E-01 -8.38E+00  5.26E+00 -5.05E+00  3.55E+00  1.72E+02
 
 TH 7
+        2.26E-02 -3.15E-02 -1.04E-02 -8.69E-02  1.25E-01  2.59E-02  1.04E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.45E+00 -4.63E+01  4.92E+01 -5.86E+00  1.62E+00  1.21E+00  2.94E-01  0.00E+00  1.91E+02
 
 TH10
+       -1.54E+00  1.53E+01 -5.32E+01  6.45E+00 -4.36E+01  3.04E-01  1.05E-01  0.00E+00 -2.67E+00  1.38E+02
 
 TH11
+       -1.16E+01 -1.12E+00 -1.13E+01 -1.07E+01  1.31E+01  1.57E+00  7.04E-03  0.00E+00  7.88E+00  2.35E+01  8.55E+01
 
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
 
 Elapsed finaloutput time in seconds:    93.06
 #CPUT: Total CPU Time in Seconds,       34.080
Stop Time:
Wed Sep 29 23:31:58 CDT 2021
