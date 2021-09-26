Sat Sep 25 14:36:24 CDT 2021
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
$DATA ../../../../data/spa/D/dat72.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m72.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13272.4860343244        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.2971E+02  2.1563E+02  2.9893E+01 -8.2068E+01  1.8178E+02 -2.4818E+03 -8.5653E+02 -9.7036E+01 -1.5854E+03 -6.9338E+02
            -2.4260E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -634.614754562223        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.7127E+00  9.3765E-01  8.2188E-01  1.7333E+00  1.3706E+00  2.3969E+00  1.1906E+00  9.8899E-01  1.3510E+00  1.1424E+00
             1.3815E+01
 PARAMETER:  6.3808E-01  3.5622E-02 -9.6162E-02  6.5003E-01  4.1524E-01  9.7416E-01  2.7445E-01  8.8931E-02  4.0085E-01  2.3316E-01
             2.7258E+00
 GRADIENT:   6.0083E+01  1.5206E+01 -1.0231E+00 -1.6205E+00 -7.3916E+00  3.1703E+01 -1.4120E+00  7.0283E+00  3.1723E+00  9.4982E-01
             6.2665E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -654.005389060022        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.6198E+00  9.0146E-01  1.2889E+00  1.8293E+00  3.6874E+00  2.2405E+00  2.2556E+00  2.7447E-01  1.8634E+00  5.7247E+00
             1.2179E+01
 PARAMETER:  5.8227E-01 -3.7391E-03  3.5382E-01  7.0392E-01  1.4049E+00  9.0670E-01  9.1344E-01 -1.1929E+00  7.2242E-01  1.8448E+00
             2.5997E+00
 GRADIENT:   5.3229E+01  2.1465E+01  1.3700E+01  1.3826E+01 -5.4477E+00 -1.9547E+01  3.7655E+00  7.1790E-03 -7.2402E-02 -3.6802E+00
             4.9634E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -655.333792970921        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:      265
 NPARAMETR:  1.6029E+00  9.2772E-01  1.2718E+00  1.7918E+00  3.8965E+00  2.2187E+00  2.2313E+00  2.5388E-01  1.8294E+00  6.9483E+00
             1.2136E+01
 PARAMETER:  5.7184E-01  2.4974E-02  3.4043E-01  6.8321E-01  1.4601E+00  8.9693E-01  9.0256E-01 -1.2709E+00  7.0401E-01  2.0385E+00
             2.5962E+00
 GRADIENT:   4.8455E+01  2.2138E+01  1.4507E+01  1.0343E+01 -3.8909E+00 -2.1711E+01  3.6874E+00  5.6094E-03 -1.0351E+00 -4.4304E+00
             4.4739E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -800.902483304015        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  8.8430E-01  9.5658E-02  1.1085E-01  8.0733E-01  2.8230E+02  2.4234E+00  1.2197E-01  1.2451E+00  8.9308E-01  6.9329E+00
             8.8040E+00
 PARAMETER: -2.2964E-02 -2.2470E+00 -2.0996E+00 -1.1402E-01  5.7430E+00  9.8517E-01 -2.0040E+00  3.1924E-01 -1.3081E-02  2.0363E+00
             2.2752E+00
 GRADIENT:   2.8398E+01  3.3826E+01 -2.2955E+01  1.1490E+01 -1.0148E-01  4.6722E+01  1.2638E-01  7.8224E+00  2.4908E+01  4.3175E-03
            -3.8407E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -831.177595300426        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  5.4656E-01  2.6988E-02  3.6242E-02  3.6574E-01  2.3752E+03  1.9521E+00  3.0244E-02  1.0085E+00  3.2045E-01  6.9260E+00
             8.7616E+00
 PARAMETER: -5.0410E-01 -3.5124E+00 -3.2175E+00 -9.0583E-01  7.8729E+00  7.6891E-01 -3.3984E+00  1.0844E-01 -1.0380E+00  2.0353E+00
             2.2704E+00
 GRADIENT:  -3.7431E+00  4.1847E+00  1.4325E+00  1.4727E+00 -1.7070E-03  3.2177E-02  2.2343E-03  1.7987E+00  2.1002E+00  1.5701E-05
            -9.1432E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -831.855813330312        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      618
 NPARAMETR:  5.8068E-01  2.8806E-02  4.1669E-02  4.0607E-01  1.4906E+03  1.9688E+00  3.5484E-02  1.0478E+00  3.1075E-01  6.9267E+00
             8.9265E+00
 PARAMETER: -4.4356E-01 -3.4472E+00 -3.0780E+00 -8.0123E-01  7.4069E+00  7.7744E-01 -3.2387E+00  1.4671E-01 -1.0688E+00  2.0354E+00
             2.2890E+00
 GRADIENT:  -4.1089E+00 -1.9486E-01 -4.5986E+00  6.2510E+00  2.9582E-05 -4.8496E-01  1.9617E-03  1.9099E+00  2.2786E+00  2.5376E-05
             9.1678E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -832.925636085625        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      793
 NPARAMETR:  6.0818E-01  3.3271E-02  4.6598E-02  4.3867E-01  6.4570E+02  1.9755E+00  3.4210E-02  1.1329E+00  4.9799E-02  6.9265E+00
             8.9494E+00
 PARAMETER: -3.9728E-01 -3.3031E+00 -2.9662E+00 -7.2400E-01  6.5703E+00  7.8080E-01 -3.2752E+00  2.2482E-01 -2.8998E+00  2.0354E+00
             2.2916E+00
 GRADIENT:  -7.4121E-01  1.7095E+00 -3.5971E+00  3.3712E+00 -2.2078E-03 -7.2737E-02  3.1502E-03 -2.5257E-01  5.0388E-02  2.3283E-04
            -9.5676E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -832.984652513037        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      968
 NPARAMETR:  6.0655E-01  3.1880E-02  4.6336E-02  4.3542E-01  4.4565E+02  1.9757E+00  3.0868E-02  1.1318E+00  1.9105E-02  6.9260E+00
             8.9477E+00
 PARAMETER: -3.9997E-01 -3.3458E+00 -2.9718E+00 -7.3145E-01  6.1995E+00  7.8090E-01 -3.3780E+00  2.2378E-01 -3.8578E+00  2.0353E+00
             2.2914E+00
 GRADIENT:  -2.2585E-02  1.1075E-01  7.5100E-02 -2.9744E-01 -2.8865E-04 -8.7257E-03  1.9586E-03 -4.1576E-02  7.6766E-03  3.7397E-04
            -1.1362E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -832.987776338035        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1146
 NPARAMETR:  6.0694E-01  3.1813E-02  4.6413E-02  4.3598E-01  2.6137E+02  1.9759E+00  2.5953E-02  1.1328E+00  1.0000E-02  6.9254E+00
             8.9473E+00
 PARAMETER: -3.9932E-01 -3.3479E+00 -2.9702E+00 -7.3016E-01  5.6659E+00  7.8103E-01 -3.5515E+00  2.2471E-01 -5.3801E+00  2.0352E+00
             2.2914E+00
 GRADIENT:  -2.5568E-03  2.0930E-03 -3.7469E-02  3.8595E-02 -5.4693E-04 -3.8424E-03  1.3624E-03 -1.7852E-03  0.0000E+00  1.0743E-03
            -3.6011E-02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -832.987777464007        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1203
 NPARAMETR:  6.0701E-01  3.1820E-02  4.6427E-02  4.3607E-01  2.5731E+02  1.9759E+00  2.5779E-02  1.1329E+00  1.0000E-02  6.9253E+00
             8.9477E+00
 PARAMETER: -3.9920E-01 -3.3477E+00 -2.9699E+00 -7.2995E-01  5.6503E+00  7.8104E-01 -3.5582E+00  2.2480E-01 -5.4232E+00  2.0352E+00
             2.2914E+00
 GRADIENT:  -1.2606E-03  1.3485E-03 -1.0206E-02  7.9924E-03 -5.7313E-04 -1.1163E-03  1.3445E-03 -2.9490E-04  0.0000E+00  1.1088E-03
            -5.3873E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1203
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.6665E-03 -2.0204E-04  1.2142E-02 -3.1998E-04  1.3626E-04
 SE:             2.9164E-02  5.3468E-05  2.1541E-02  2.7929E-04  6.4406E-05
 N:                     100         100         100         100         100

 P VAL.:         8.4594E-01  1.5775E-04  5.7299E-01  2.5192E-01  3.4374E-02

 ETASHRINKSD(%)  2.2961E+00  9.9821E+01  2.7833E+01  9.9064E+01  9.9784E+01
 ETASHRINKVR(%)  4.5395E+00  1.0000E+02  4.7920E+01  9.9991E+01  1.0000E+02
 EBVSHRINKSD(%)  2.0137E+00  9.9700E+01  2.7872E+01  9.9034E+01  9.9690E+01
 EBVSHRINKVR(%)  3.9868E+00  9.9999E+01  4.7976E+01  9.9991E+01  9.9999E+01
 RELATIVEINF(%)  5.2631E+00  8.8895E-06  8.2840E-01  7.9552E-05  5.8035E-06
 EPSSHRINKSD(%)  1.2557E+01
 EPSSHRINKVR(%)  2.3537E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -832.98777746400742     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -97.836950900269244     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.91
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.61
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -832.988       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.07E-01  3.18E-02  4.64E-02  4.36E-01  2.57E+02  1.98E+00  2.58E-02  1.13E+00  1.00E-02  6.93E+00  8.95E+00
 


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
+        7.35E+02
 
 TH 2
+       -2.98E+02  2.10E+04
 
 TH 3
+       -8.15E+02  4.27E+03  2.04E+05
 
 TH 4
+       -3.54E+02 -2.61E+03 -3.02E+04  5.20E+03
 
 TH 5
+        2.64E-04 -1.02E-02 -6.96E-03  1.92E-03  1.62E-08
 
 TH 6
+        4.68E+00  4.33E+01  1.46E+02 -3.87E+01  2.11E-05  4.57E+01
 
 TH 7
+        1.53E-01  6.90E+00 -1.58E+00 -5.46E-01 -7.53E-05 -6.99E-02  3.31E+00
 
 TH 8
+        2.34E-01 -5.44E+01 -4.28E+01 -4.63E+01  2.41E-05  1.54E+00 -1.45E+00  4.71E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -5.00E-03  1.75E-02  5.82E-03  2.15E-03 -2.72E-07  1.24E-03  4.44E-03 -5.84E-03  0.00E+00 -1.06E-04
 
 TH11
+       -9.47E+00 -1.68E+01  1.49E+02 -2.05E+01 -2.39E-06  1.06E+00  6.31E-03  2.94E+00  0.00E+00 -2.14E-06  4.79E+00
 
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
 #CPUT: Total CPU Time in Seconds,       27.025
Stop Time:
Sat Sep 25 14:36:56 CDT 2021
