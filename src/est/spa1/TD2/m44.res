Thu Sep 30 02:07:25 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat44.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m44.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2022.24982907087        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.6880E+02  6.6148E+01  1.1371E+01  1.2390E+02  2.6522E+01  2.4370E+01  5.2425E+00 -2.1053E+01  2.3803E+01 -2.9571E+00
            -3.8819E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2031.56899894263        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  9.3344E-01  9.8559E-01  9.5388E-01  9.8897E-01  9.6192E-01  1.0139E+00  9.9062E-01  1.0744E+00  9.3955E-01  1.0029E+00
             1.0755E+00
 PARAMETER:  3.1124E-02  8.5488E-02  5.2786E-02  8.8906E-02  6.1174E-02  1.1378E-01  9.0577E-02  1.7177E-01  3.7647E-02  1.0287E-01
             1.7275E-01
 GRADIENT:   1.8313E+00  1.3201E+01  3.2600E+00  1.7724E+01  2.4575E-01 -7.5491E+00  9.2145E-03 -1.1077E+01  7.7306E+00  6.8245E+00
             2.3517E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2032.74134204725        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      329
 NPARAMETR:  9.3607E-01  8.9360E-01  8.4157E-01  1.0290E+00  8.5753E-01  1.0289E+00  1.2656E+00  1.2944E+00  7.5646E-01  7.4887E-01
             1.0172E+00
 PARAMETER:  3.3937E-02 -1.2499E-02 -7.2489E-02  1.2862E-01 -5.3694E-02  1.2845E-01  3.3554E-01  3.5807E-01 -1.7910E-01 -1.8919E-01
             1.1709E-01
 GRADIENT:   7.9505E+00  3.4077E+00 -1.1582E+01  1.3050E+01  1.6029E+01 -2.1906E+00  9.2354E-01  6.2915E+00 -7.3053E+00 -5.6567E+00
            -2.0292E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2034.50136480285        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      506
 NPARAMETR:  9.2931E-01  7.5695E-01  1.1710E+00  1.1403E+00  9.4886E-01  1.0305E+00  1.3370E+00  1.3422E+00  7.7417E-01  9.5239E-01
             1.0472E+00
 PARAMETER:  2.6692E-02 -1.7846E-01  2.5782E-01  2.3128E-01  4.7506E-02  1.3009E-01  3.9042E-01  3.9432E-01 -1.5596E-01  5.1223E-02
             1.4608E-01
 GRADIENT:  -1.1084E+00  1.4208E+01  3.5434E+00  1.9165E+01 -7.8688E+00  1.2778E-01  6.9056E-01 -2.1487E+00 -4.9950E-01  1.2822E+00
             2.0080E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2034.93031629306        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      690             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2758E-01  6.1649E-01  1.3642E+00  1.2067E+00  9.8556E-01  1.0298E+00  1.4253E+00  1.4824E+00  7.6404E-01  9.9913E-01
             1.0460E+00
 PARAMETER:  2.4823E-02 -3.8372E-01  4.1054E-01  2.8786E-01  8.5452E-02  1.2935E-01  4.5438E-01  4.9364E-01 -1.6914E-01  9.9134E-02
             1.4500E-01
 GRADIENT:   3.7600E+02  4.2817E+01  6.2361E+00  2.7490E+02  1.2838E+01  6.7036E+01  1.0625E+01  2.9557E+00  9.3981E+00  1.0766E+00
             2.7634E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2035.17361920134        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      869
 NPARAMETR:  9.2871E-01  6.3872E-01  1.3562E+00  1.2120E+00  9.8650E-01  1.0291E+00  1.3763E+00  1.4938E+00  7.6034E-01  9.9436E-01
             1.0448E+00
 PARAMETER:  2.6045E-02 -3.4829E-01  4.0466E-01  2.9231E-01  8.6407E-02  1.2866E-01  4.1937E-01  5.0131E-01 -1.7399E-01  9.4340E-02
             1.4386E-01
 GRADIENT:   1.3017E+00  2.5142E+00 -2.3803E+00 -1.1225E+00  4.3091E+00  3.0427E-01 -2.2572E-01  3.6234E-01 -2.9128E+00  7.1532E-02
            -2.6780E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2037.24744076031        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1049
 NPARAMETR:  9.2249E-01  4.1242E-01  1.4601E+00  1.3582E+00  9.4713E-01  1.0235E+00  4.3019E-01  1.5214E+00  7.7403E-01  9.8234E-01
             1.0413E+00
 PARAMETER:  1.9320E-02 -7.8570E-01  4.7849E-01  4.0616E-01  4.5684E-02  1.2320E-01 -7.4353E-01  5.1963E-01 -1.5615E-01  8.2184E-02
             1.4044E-01
 GRADIENT:  -6.5048E+00  6.0159E+00 -4.3483E+00  1.8301E+01  1.9951E+00 -8.2543E-01  1.8559E-01 -2.0185E+00 -9.9517E+00 -2.0200E+00
            -2.7383E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2037.59244710545        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1225
 NPARAMETR:  9.3002E-01  3.3132E-01  1.5644E+00  1.4107E+00  9.5047E-01  1.0307E+00  2.5265E-01  1.6130E+00  7.6664E-01  9.9663E-01
             1.0442E+00
 PARAMETER:  2.7451E-02 -1.0047E+00  5.4749E-01  4.4405E-01  4.9206E-02  1.3020E-01 -1.2758E+00  5.7807E-01 -1.6574E-01  9.6629E-02
             1.4328E-01
 GRADIENT:   1.3157E+01  2.3339E+00 -1.0898E+00  1.1697E+01 -3.8624E+00  2.2587E+00  7.2353E-02 -4.6916E-01 -4.1414E-01 -2.6456E-01
            -1.6601E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2037.65841511039        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1403
 NPARAMETR:  9.2733E-01  3.0101E-01  1.7422E+00  1.4362E+00  9.9608E-01  1.0283E+00  2.1477E-01  1.7545E+00  7.5423E-01  1.0366E+00
             1.0448E+00
 PARAMETER:  2.4559E-02 -1.1006E+00  6.5515E-01  4.6198E-01  9.6076E-02  1.2788E-01 -1.4382E+00  6.6218E-01 -1.8205E-01  1.3591E-01
             1.4384E-01
 GRADIENT:   7.8132E+00  1.2957E+00 -6.4332E-01  7.7352E+00  1.4511E+00  1.3943E+00  4.6844E-02  1.1131E-01 -1.7886E-01 -5.1986E-01
            -2.0930E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2037.70095787363        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1583
 NPARAMETR:  9.2780E-01  3.0267E-01  1.7479E+00  1.4329E+00  9.9701E-01  1.0242E+00  5.2940E-02  1.7487E+00  7.5563E-01  1.0439E+00
             1.0441E+00
 PARAMETER:  2.5059E-02 -1.0951E+00  6.5842E-01  4.5967E-01  9.7006E-02  1.2394E-01 -2.8386E+00  6.5887E-01 -1.8020E-01  1.4294E-01
             1.4312E-01
 GRADIENT:   8.9114E+00  5.4643E-01  6.2968E-01  1.8437E+00 -8.8364E-02 -1.6493E-01  3.5779E-03 -4.3556E-01  6.5404E-02  1.1317E-01
            -7.1060E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2037.70283840984        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1763
 NPARAMETR:  9.2779E-01  3.0267E-01  1.7479E+00  1.4328E+00  9.9655E-01  1.0250E+00  1.0000E-02  1.7487E+00  7.5572E-01  1.0434E+00
             1.0441E+00
 PARAMETER:  2.5055E-02 -1.0951E+00  6.5842E-01  4.5963E-01  9.6542E-02  1.2466E-01 -7.3891E+00  6.5888E-01 -1.8009E-01  1.4252E-01
             1.4312E-01
 GRADIENT:   8.8887E+00  5.6718E-01  8.2001E-01  1.7756E+00 -5.1268E-01  1.1797E-01  0.0000E+00 -4.3579E-01  8.9793E-02  1.2711E-01
            -7.0367E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2037.71084757033        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1923
 NPARAMETR:  9.2372E-01  2.9751E-01  1.7449E+00  1.4345E+00  9.9488E-01  1.0248E+00  1.0000E-02  1.7521E+00  7.5397E-01  1.0409E+00
             1.0448E+00
 PARAMETER:  2.0657E-02 -1.1123E+00  6.5670E-01  4.6083E-01  9.4870E-02  1.2449E-01 -7.3891E+00  6.6079E-01 -1.8240E-01  1.4011E-01
             1.4382E-01
 GRADIENT:  -1.6624E-01  1.7889E-02  2.4887E-01 -1.9851E+00  3.0044E-01  1.0574E-01  0.0000E+00 -1.3354E-02  6.0732E-02  6.3088E-02
             2.1096E-02

0ITERATION NO.:   56    OBJECTIVE VALUE:  -2037.71084757033        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1945
 NPARAMETR:  9.2372E-01  2.9751E-01  1.7449E+00  1.4345E+00  9.9488E-01  1.0248E+00  1.0000E-02  1.7521E+00  7.5397E-01  1.0409E+00
             1.0448E+00
 PARAMETER:  2.0657E-02 -1.1123E+00  6.5670E-01  4.6083E-01  9.4870E-02  1.2449E-01 -7.3891E+00  6.6079E-01 -1.8240E-01  1.4011E-01
             1.4382E-01
 GRADIENT:  -1.6624E-01  1.7889E-02  2.4887E-01 -1.9851E+00  3.0044E-01  1.0574E-01  0.0000E+00 -1.3354E-02  6.0732E-02  6.3088E-02
             2.1096E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1945
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1091E-03 -1.2202E-04 -4.2916E-02 -8.4414E-03 -4.7721E-02
 SE:             2.9878E-02  5.6489E-05  1.9700E-02  2.9364E-02  2.0060E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7039E-01  3.0765E-02  2.9372E-02  7.7375E-01  1.7362E-02

 ETASHRINKSD(%)  1.0000E-10  9.9811E+01  3.4002E+01  1.6266E+00  3.2798E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  5.6443E+01  3.2267E+00  5.4839E+01
 EBVSHRINKSD(%)  3.5577E-01  9.9814E+01  3.7849E+01  2.1383E+00  2.8995E+01
 EBVSHRINKVR(%)  7.1028E-01  1.0000E+02  6.1372E+01  4.2309E+00  4.9584E+01
 RELATIVEINF(%)  9.7977E+01  2.0284E-05  1.2467E+01  6.2772E+00  1.3292E+01
 EPSSHRINKSD(%)  3.4513E+01
 EPSSHRINKVR(%)  5.7114E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2037.7108475703330     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1118.7723143656603     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.85
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2037.711       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.24E-01  2.98E-01  1.74E+00  1.43E+00  9.95E-01  1.02E+00  1.00E-02  1.75E+00  7.54E-01  1.04E+00  1.04E+00
 


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
+        1.23E+03
 
 TH 2
+       -2.55E+01  4.65E+02
 
 TH 3
+        2.27E-01 -9.98E+05  4.86E+01
 
 TH 4
+       -3.10E+06 -8.64E+05 -1.98E+01  9.11E+02
 
 TH 5
+        2.80E+00 -1.78E+02  6.25E+02 -5.33E+01  5.41E+02
 
 TH 6
+        2.56E+00 -3.53E+00  3.70E+00 -2.11E+00 -3.95E-01  1.87E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -2.76E+00  9.88E+05 -1.63E+01 -4.48E+00 -7.41E+02 -4.09E+00  0.00E+00  2.20E+01
 
 TH 9
+        3.16E+00 -1.10E+02 -2.91E+02 -4.77E+00  1.14E+00 -6.16E-01  0.00E+00  2.98E+02  3.24E+02
 
 TH10
+        1.91E+00  7.01E+00 -1.16E+03  1.71E+00 -7.33E+01  7.06E-01  0.00E+00  1.15E+03  7.51E-01  6.11E+01
 
 TH11
+       -9.12E+00 -1.29E+01  1.10E+06  1.91E+06 -1.04E+01  1.84E+00  0.00E+00 -1.10E+06  9.91E+00  9.53E+00  3.67E+02
 
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
 #CPUT: Total CPU Time in Seconds,       36.980
Stop Time:
Thu Sep 30 02:08:04 CDT 2021
