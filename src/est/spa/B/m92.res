Wed Sep 29 11:38:24 CDT 2021
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
$DATA ../../../../data/spa/B/dat92.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1653.64665841387        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4670E+02 -4.6047E+01 -6.3070E+01  2.7600E+01  8.3996E+01  5.4649E+01 -1.5972E+01  9.1992E+00 -8.1063E+00  1.7354E+01
             2.3190E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1662.52742467231        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.7862E-01  1.1098E+00  1.2164E+00  9.5436E-01  1.0491E+00  9.0869E-01  1.1304E+00  9.4313E-01  1.1194E+00  8.4748E-01
             9.9686E-01
 PARAMETER:  7.8391E-02  2.0419E-01  2.9593E-01  5.3290E-02  1.4795E-01  4.2442E-03  2.2256E-01  4.1454E-02  2.1281E-01 -6.5490E-02
             9.6854E-02
 GRADIENT:  -9.1202E+00 -1.1753E+01  2.0005E+01 -2.4015E+01 -9.1376E+00 -2.1452E+01 -2.6628E+00 -7.5515E+00  5.3616E+00 -1.0480E+01
            -8.6094E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1663.10087857835        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.7740E-01  1.2322E+00  1.1229E+00  8.8849E-01  1.0544E+00  9.4251E-01  1.1561E+00  1.0273E+00  1.1451E+00  7.8089E-01
             1.0164E+00
 PARAMETER:  7.7139E-02  3.0879E-01  2.1591E-01 -1.8232E-02  1.5297E-01  4.0793E-02  2.4508E-01  1.2694E-01  2.3546E-01 -1.4732E-01
             1.1622E-01
 GRADIENT:  -1.3932E+01  9.0364E+00  2.0334E+01 -1.1096E+01 -1.8424E+01 -6.6554E+00  4.5130E+00 -4.2599E+00  4.3303E+00 -1.5041E+01
            -1.3405E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.83951252895        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.8502E-01  1.2178E+00  1.0086E+00  9.0153E-01  1.0267E+00  9.6097E-01  1.1439E+00  8.7198E-01  1.0898E+00  8.6841E-01
             1.0072E+00
 PARAMETER:  8.4906E-02  2.9709E-01  1.0853E-01 -3.6600E-03  1.2637E-01  6.0187E-02  2.3443E-01 -3.6984E-02  1.8601E-01 -4.1097E-02
             1.0713E-01
 GRADIENT:   3.7680E+00  3.9230E+00  2.5214E+00  2.9199E+00 -5.3467E+00  1.0387E+00  3.8517E-01 -2.8874E-01  2.0634E-01  9.6153E-01
             1.9801E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1665.05240829919        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.8323E-01  1.4881E+00  8.2497E-01  7.3454E-01  1.0812E+00  9.5636E-01  9.9025E-01  8.5951E-01  1.2598E+00  8.6683E-01
             1.0073E+00
 PARAMETER:  8.3086E-02  4.9750E-01 -9.2405E-02 -2.0851E-01  1.7805E-01  5.5382E-02  9.0202E-02 -5.1389E-02  3.3094E-01 -4.2915E-02
             1.0727E-01
 GRADIENT:  -3.8381E+00  1.5255E+01 -1.5015E+00  1.4340E+01 -1.0686E+00 -1.4463E+00  4.0355E-01  9.4748E-01 -8.7250E-01 -8.5675E-01
            -5.1409E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1665.60079583232        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.8565E-01  1.7617E+00  6.0175E-01  5.5237E-01  1.1405E+00  9.5959E-01  8.7767E-01  5.0909E-01  1.5248E+00  8.9591E-01
             1.0103E+00
 PARAMETER:  8.5542E-02  6.6626E-01 -4.0791E-01 -4.9354E-01  2.3151E-01  5.8748E-02 -3.0487E-02 -5.7513E-01  5.2189E-01 -9.9161E-03
             1.1023E-01
 GRADIENT:  -5.2098E-02  1.6327E+01 -3.0852E+00  1.4362E+01  2.1740E+00 -4.0087E-01 -8.9435E-01  4.4986E-01 -1.0227E+00  5.4098E-01
             4.0668E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1666.04472597822        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  9.8621E-01  1.8726E+00  5.3697E-01  4.5425E-01  1.1863E+00  9.6073E-01  8.4164E-01  2.3520E-01  1.7243E+00  9.1956E-01
             1.0111E+00
 PARAMETER:  8.6110E-02  7.2731E-01 -5.2181E-01 -6.8910E-01  2.7086E-01  5.9943E-02 -7.2408E-02 -1.3473E+00  6.4483E-01  1.6138E-02
             1.1102E-01
 GRADIENT:   1.7699E+00 -1.2609E+01  2.6781E+00 -5.5725E+00 -5.8798E+00  1.8602E-01  9.2077E-02  7.3945E-02 -2.3988E+00  4.4144E-01
             3.0129E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1666.08603616856        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  9.8549E-01  1.8956E+00  5.2371E-01  4.4897E-01  1.1920E+00  9.6041E-01  8.3674E-01  1.2400E-01  1.7619E+00  9.1883E-01
             1.0105E+00
 PARAMETER:  8.5379E-02  7.3952E-01 -5.4681E-01 -7.0081E-01  2.7564E-01  5.9606E-02 -7.8241E-02 -1.9874E+00  6.6639E-01  1.5341E-02
             1.1040E-01
 GRADIENT:  -2.5431E-01 -1.9545E-01  1.0539E+00  1.9925E+00 -2.4428E+00 -2.0656E-02  4.8886E-02  2.3695E-02 -5.1835E-01 -1.9952E-01
            -1.9766E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1666.12754437909        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1424             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8623E-01  1.8907E+00  5.2138E-01  4.4536E-01  1.1921E+00  9.6073E-01  8.3510E-01  1.6245E-02  1.7693E+00  9.1831E-01
             1.0101E+00
 PARAMETER:  8.6133E-02  7.3693E-01 -5.5128E-01 -7.0887E-01  2.7571E-01  5.9942E-02 -8.0205E-02 -4.0200E+00  6.7057E-01  1.4781E-02
             1.1002E-01
 GRADIENT:   4.1279E+02  8.3615E+02  3.2979E+00  9.9737E+01  1.4573E+01  4.1686E+01  8.7901E+00  8.6041E-04  2.7312E+01  4.9824E-01
             7.1936E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1666.13018573198        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1608
 NPARAMETR:  9.8621E-01  1.8901E+00  5.2037E-01  4.4590E-01  1.1919E+00  9.6073E-01  8.3549E-01  1.0000E-02  1.7703E+00  9.1827E-01
             1.0102E+00
 PARAMETER:  8.6119E-02  7.3661E-01 -5.5322E-01 -7.0766E-01  2.7558E-01  5.9942E-02 -7.9736E-02 -4.7538E+00  6.7114E-01  1.4740E-02
             1.1011E-01
 GRADIENT:   1.6269E+00 -1.1632E+01  2.7775E-01 -7.5396E-01 -6.6651E-01  1.2237E-01 -2.7433E-02  0.0000E+00  1.9214E-02  2.7471E-02
            -4.8293E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1666.13018573198        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1630
 NPARAMETR:  9.8621E-01  1.8901E+00  5.2037E-01  4.4590E-01  1.1919E+00  9.6073E-01  8.3549E-01  1.0000E-02  1.7703E+00  9.1827E-01
             1.0102E+00
 PARAMETER:  8.6119E-02  7.3661E-01 -5.5322E-01 -7.0766E-01  2.7558E-01  5.9942E-02 -7.9736E-02 -4.7538E+00  6.7114E-01  1.4740E-02
             1.1011E-01
 GRADIENT:   1.6269E+00 -1.1632E+01  2.7775E-01 -7.5396E-01 -6.6651E-01  1.2237E-01 -2.7433E-02  0.0000E+00  1.9214E-02  2.7471E-02
            -4.8293E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1630
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1192E-04 -3.1913E-02 -2.3096E-04  3.2976E-02 -4.0475E-02
 SE:             2.9850E-02  2.4671E-02  8.4429E-05  2.2055E-02  2.1713E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9701E-01  1.9582E-01  6.2276E-03  1.3487E-01  6.2313E-02

 ETASHRINKSD(%)  1.0000E-10  1.7349E+01  9.9717E+01  2.6113E+01  2.7258E+01
 ETASHRINKVR(%)  1.0000E-10  3.1689E+01  9.9999E+01  4.5407E+01  4.7086E+01
 EBVSHRINKSD(%)  4.6942E-01  1.6280E+01  9.9759E+01  2.8944E+01  2.5752E+01
 EBVSHRINKVR(%)  9.3663E-01  2.9910E+01  9.9999E+01  4.9511E+01  4.4873E+01
 RELATIVEINF(%)  9.9022E+01  7.3246E+00  9.3772E-05  5.0026E+00  1.5783E+01
 EPSSHRINKSD(%)  4.4205E+01
 EPSSHRINKVR(%)  6.8869E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1666.1301857319752     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -930.97935916823701     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.43
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1666.130       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.89E+00  5.20E-01  4.46E-01  1.19E+00  9.61E-01  8.35E-01  1.00E-02  1.77E+00  9.18E-01  1.01E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.28E+03
 
 TH 2
+       -5.87E+01  3.42E+02
 
 TH 3
+        5.40E+01  1.37E+02  2.80E+02
 
 TH 4
+       -1.77E+02  2.91E+02 -2.85E+02  9.95E+02
 
 TH 5
+       -1.19E+01 -1.56E+02 -2.58E+02  2.78E+02  4.66E+02
 
 TH 6
+       -5.09E+01  2.88E-01 -1.64E+01  2.94E+01  9.52E+00  2.21E+02
 
 TH 7
+       -6.19E+00 -1.37E+01 -2.44E+01 -3.03E+01 -2.63E+01 -1.39E+00  1.52E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.26E+01 -5.89E+00 -4.05E+01  5.43E+01  1.29E+01  2.90E+00  1.73E+01  0.00E+00  9.45E+00
 
 TH10
+       -2.03E+01 -3.12E+01 -4.38E+01  9.95E-01 -6.08E+01  4.87E+00  2.63E+01  0.00E+00  1.62E+01  5.37E+01
 
 TH11
+       -7.36E+01  7.64E-01 -1.77E+01  3.17E+01 -1.94E+01  4.31E+00  4.12E+00  0.00E+00  5.95E+00  2.54E+01  3.15E+02
 
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
+       -5.25E+00  3.29E+02
 
 TH 3
+        7.39E+00  1.15E+02  2.77E+02
 
 TH 4
+       -1.45E+01  2.66E+02 -2.70E+02  9.48E+02
 
 TH 5
+       -5.63E+00 -1.43E+02 -2.36E+02  2.66E+02  5.15E+02
 
 TH 6
+       -2.57E-01 -6.76E-01  1.77E+00 -3.50E+00 -7.70E-01  2.12E+02
 
 TH 7
+        6.90E-01  6.73E+00 -1.02E+01 -1.50E+01 -1.34E+01 -5.92E-01  1.48E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.15E-01 -1.42E+01 -2.91E+01  6.00E+01  2.88E-01 -4.26E-01  1.39E+01  0.00E+00  2.42E+01
 
 TH10
+        2.98E-02 -1.40E+01 -2.51E+01 -3.54E+00 -6.63E+01  7.04E-01  9.39E+00  0.00E+00  5.15E+00  8.92E+01
 
 TH11
+       -7.51E+00 -1.56E+01 -1.82E+01  1.67E+00 -8.24E+00  2.61E+00  8.56E+00  0.00E+00  3.72E+00  1.83E+01  2.11E+02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.23E+03
 
 TH 2
+        4.53E+01  3.18E+02
 
 TH 3
+       -5.32E+01  1.04E+02  2.36E+02
 
 TH 4
+        1.27E+02  2.55E+02 -2.71E+02  9.27E+02
 
 TH 5
+       -1.94E+01 -1.50E+02 -2.33E+02  2.64E+02  5.31E+02
 
 TH 6
+        4.21E+01 -2.18E-01  1.30E+01 -2.90E+01 -1.12E+01  2.07E+02
 
 TH 7
+        9.14E+00  2.32E+01  5.16E+00  1.29E+00  1.08E+01  1.73E+00  1.49E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.26E-01 -1.55E+01 -3.02E+01  5.19E+01  7.32E+00 -3.37E+00  1.19E+01  0.00E+00  1.87E+01
 
 TH10
+        3.86E+00 -1.08E+00 -1.57E+01 -4.32E+00 -1.15E+02 -7.06E+00 -8.22E+00  0.00E+00  3.02E+00  1.37E+02
 
 TH11
+        3.53E+01 -2.64E+01 -1.12E+01 -9.93E+00  1.17E+00  3.39E+00  1.05E+01  0.00E+00  5.92E+00  1.51E+01  1.42E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       27.613
Stop Time:
Wed Sep 29 11:38:53 CDT 2021
