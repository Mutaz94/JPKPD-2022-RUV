Thu Sep 30 00:23:36 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat59.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:   387.472913462847        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3257E+02  1.1385E+02  3.4879E+02  2.3334E+01  2.7672E+02  3.0920E+01 -6.6105E+01 -5.0751E+02 -1.0006E+02 -2.2566E+02
            -4.0817E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1469.49346854078        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1122E+00  9.8383E-01  8.1967E-01  1.1796E+00  8.8510E-01  9.4617E-01  9.0471E-01  1.0408E+00  8.4326E-01  1.0399E+00
             5.2144E+00
 PARAMETER:  2.0636E-01  8.3694E-02 -9.8853E-02  2.6515E-01 -2.2049E-02  4.4664E-02 -1.3869E-04  1.3997E-01 -7.0484E-02  1.3914E-01
             1.7514E+00
 GRADIENT:   2.7836E+01  1.2827E+01 -2.5613E+01  5.7529E+01  1.9203E-01 -8.4260E+00  8.2365E+00  1.0111E+01  1.3402E+01  2.6299E+01
             2.7668E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1494.64710974466        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      168
 NPARAMETR:  1.0990E+00  6.9158E-01  3.0012E-01  1.2393E+00  3.7091E-01  1.0982E+00  2.4600E-01  2.5701E-01  1.1566E+00  4.6802E-01
             4.7043E+00
 PARAMETER:  1.9437E-01 -2.6878E-01 -1.1036E+00  3.1452E-01 -8.9178E-01  1.9365E-01 -1.3024E+00 -1.2586E+00  2.4544E-01 -6.5925E-01
             1.6485E+00
 GRADIENT:  -1.5328E-01  1.1305E+02  5.3858E+01  1.3922E+02 -1.2931E+02  2.3955E+01 -3.8553E-01  4.6646E-01  9.9442E+00 -1.9115E-01
             2.3367E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1537.20136182383        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      261
 NPARAMETR:  1.0493E+00  4.5866E-01  1.7139E-01  1.1893E+00  2.4539E-01  1.0276E+00  9.8972E-02  5.1546E-02  1.6218E+00  5.9645E-01
             3.3978E+00
 PARAMETER:  1.4814E-01 -6.7945E-01 -1.6638E+00  2.7335E-01 -1.3049E+00  1.2725E-01 -2.2129E+00 -2.8653E+00  5.8356E-01 -4.1676E-01
             1.3231E+00
 GRADIENT:  -6.6731E+01  6.3175E+01  6.7697E-01  1.0957E+02 -9.4190E+01 -1.4779E+01  9.4038E-03 -3.1709E-02  3.2431E+01 -4.5774E+00
             3.9575E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1558.91856437972        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      436
 NPARAMETR:  1.0805E+00  4.3198E-01  2.2249E-01  1.0976E+00  2.8010E-01  1.0432E+00  1.6951E-01  2.3411E-02  1.1715E+00  6.1670E-01
             3.2558E+00
 PARAMETER:  1.7742E-01 -7.3938E-01 -1.4029E+00  1.9309E-01 -1.1726E+00  1.4232E-01 -1.6749E+00 -3.6546E+00  2.5829E-01 -3.8338E-01
             1.2804E+00
 GRADIENT:  -2.1210E+00  3.9172E+00  3.6026E+00  2.9418E+00 -5.2719E+00 -1.0887E+00  1.1336E-01 -5.7587E-03 -2.6608E-01  5.5842E-01
            -5.8183E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1559.04301953535        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      615
 NPARAMETR:  1.0819E+00  4.1904E-01  2.1647E-01  1.0918E+00  2.7413E-01  1.0466E+00  5.3750E-02  2.0340E-02  1.1759E+00  6.1204E-01
             3.2659E+00
 PARAMETER:  1.7874E-01 -7.6980E-01 -1.4303E+00  1.8782E-01 -1.1941E+00  1.4555E-01 -2.8234E+00 -3.7952E+00  2.6207E-01 -3.9096E-01
             1.2835E+00
 GRADIENT:   4.4767E-01 -5.4434E-01 -1.6317E+00  1.1602E+00  2.0053E+00 -4.4220E-02  1.1275E-02 -4.6205E-03 -7.7999E-01 -2.7420E-01
            -1.5970E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1559.05220262010        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      792
 NPARAMETR:  1.0816E+00  4.1902E-01  2.1644E-01  1.0908E+00  2.7392E-01  1.0466E+00  1.4943E-02  2.1321E-02  1.1791E+00  6.1354E-01
             3.2704E+00
 PARAMETER:  1.7840E-01 -7.6983E-01 -1.4305E+00  1.8687E-01 -1.1949E+00  1.4559E-01 -4.1035E+00 -3.7481E+00  2.6479E-01 -3.8851E-01
             1.2849E+00
 GRADIENT:  -3.2646E-01 -1.7696E-01 -2.6958E-01 -3.3605E-03  2.9016E-01  6.1277E-03  9.0396E-04 -4.9149E-03 -5.6636E-02  5.7038E-02
             9.9783E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1559.79871505794        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      974
 NPARAMETR:  1.0780E+00  4.1463E-01  2.0922E-01  1.0759E+00  2.6898E-01  1.0490E+00  1.0000E-02  6.7632E-01  1.1918E+00  6.1474E-01
             3.2565E+00
 PARAMETER:  1.7515E-01 -7.8037E-01 -1.4644E+00  1.7318E-01 -1.2131E+00  1.4779E-01 -4.4645E+01 -2.9109E-01  2.7548E-01 -3.8656E-01
             1.2807E+00
 GRADIENT:  -8.0606E+00 -8.1723E+00 -2.1644E+00 -8.1470E+00  1.4940E+01 -2.5014E-02  0.0000E+00  1.3140E+00  2.1941E+00  1.1181E+01
             2.1027E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1560.93605645052        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1151
 NPARAMETR:  1.0831E+00  4.2531E-01  2.0550E-01  1.0761E+00  2.6800E-01  1.0528E+00  1.0000E-02  8.3053E-01  1.2010E+00  5.1451E-01
             3.1787E+00
 PARAMETER:  1.7980E-01 -7.5494E-01 -1.4823E+00  1.7334E-01 -1.2168E+00  1.5141E-01 -5.1795E+01 -8.5694E-02  2.8313E-01 -5.6455E-01
             1.2565E+00
 GRADIENT:   1.9782E+00  4.9847E-01  1.2786E+00  1.1261E+00  1.1340E+00  2.0819E-02  0.0000E+00 -1.2423E-01  2.2804E-01 -3.5697E-01
             3.4176E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1560.97585987850        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1328
 NPARAMETR:  1.0815E+00  4.1969E-01  2.0074E-01  1.0697E+00  2.6360E-01  1.0529E+00  1.0000E-02  8.4401E-01  1.2092E+00  5.1676E-01
             3.1684E+00
 PARAMETER:  1.7835E-01 -7.6823E-01 -1.5058E+00  1.6737E-01 -1.2333E+00  1.5154E-01 -8.1108E+01 -6.9596E-02  2.8996E-01 -5.6017E-01
             1.2532E+00
 GRADIENT:  -6.1075E-01  5.4293E-01  5.8046E-01  2.8879E-01 -1.0219E+00 -1.3098E-01  0.0000E+00 -9.9315E-04  3.3223E-02  1.4426E-02
            -1.5743E-02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1560.97684389748        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1385
 NPARAMETR:  1.0817E+00  4.1876E-01  2.0033E-01  1.0691E+00  2.6322E-01  1.0533E+00  1.0000E-02  8.4552E-01  1.2095E+00  5.1641E-01
             3.1677E+00
 PARAMETER:  1.7856E-01 -7.7046E-01 -1.5078E+00  1.6680E-01 -1.2347E+00  1.5189E-01 -8.3419E+01 -6.7804E-02  2.9020E-01 -5.6086E-01
             1.2530E+00
 GRADIENT:  -1.4750E-01  3.0688E-02  8.2112E-02 -5.4628E-02 -9.2930E-02 -9.1615E-03  0.0000E+00  1.0166E-02 -2.0496E-02  9.3516E-03
             3.9944E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1385
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3703E-03 -2.0040E-04  1.0074E-02 -7.5058E-03  9.5114E-03
 SE:             2.9093E-02  1.3247E-04  1.5475E-02  2.6430E-02  1.7629E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6243E-01  1.3034E-01  5.1507E-01  7.7642E-01  5.8953E-01

 ETASHRINKSD(%)  2.5348E+00  9.9556E+01  4.8155E+01  1.1455E+01  4.0940E+01
 ETASHRINKVR(%)  5.0053E+00  9.9998E+01  7.3121E+01  2.1598E+01  6.5119E+01
 EBVSHRINKSD(%)  2.5184E+00  9.9551E+01  4.8305E+01  9.8448E+00  4.0995E+01
 EBVSHRINKVR(%)  4.9733E+00  9.9998E+01  7.3276E+01  1.8720E+01  6.5184E+01
 RELATIVEINF(%)  9.4788E+01  2.8347E-04  2.6677E+00  4.5121E+01  1.5343E+00
 EPSSHRINKSD(%)  2.5322E+01
 EPSSHRINKVR(%)  4.4232E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1560.9768438974850     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -642.03831069281227     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.15
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1560.977       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.08E+00  4.19E-01  2.00E-01  1.07E+00  2.63E-01  1.05E+00  1.00E-02  8.46E-01  1.21E+00  5.16E-01  3.17E+00
 


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
+        8.14E+02
 
 TH 2
+       -4.75E+01  1.36E+03
 
 TH 3
+       -6.11E+01  2.30E+03  1.02E+04
 
 TH 4
+       -2.14E+01  1.98E+02 -3.89E+02  5.38E+02
 
 TH 5
+        1.75E+02 -4.66E+03 -1.26E+04 -5.35E+02  2.16E+04
 
 TH 6
+        8.92E-01 -1.32E+01  2.11E+01 -7.02E+00  2.30E+01  1.62E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -2.32E+00 -1.83E+01  1.51E+01  1.65E+00  3.46E+01 -1.80E-01  0.00E+00  1.73E+01
 
 TH 9
+        6.66E+00 -4.24E+01  1.22E+02 -8.71E+00  1.03E+02  7.78E-01  0.00E+00  5.37E-01  8.41E+01
 
 TH10
+       -1.52E+00 -5.28E+01 -3.70E+01  1.57E+00  1.60E+02  1.95E+00  0.00E+00  2.74E+01  4.90E+00  1.24E+02
 
 TH11
+       -1.36E+01 -1.70E+01 -4.66E+01 -6.86E+00  4.14E+01  1.78E+00  0.00E+00  9.72E+00  6.57E+00  1.15E+01  4.67E+01
 
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
 #CPUT: Total CPU Time in Seconds,       29.788
Stop Time:
Thu Sep 30 00:24:27 CDT 2021
