Wed Sep 29 13:53:00 CDT 2021
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
$DATA ../../../../data/spa/A3/dat88.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   513.195595635791        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5168E+02  1.2740E+02  9.0412E+01  5.6181E+01  2.9711E+02  4.6821E+01 -9.5228E+01 -4.5703E+01 -2.1330E+02 -2.4590E+02
            -3.7164E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -875.299218292519        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0647E+00  5.7661E-01  5.2284E-01  1.3235E+00  5.3559E-01  8.4202E-01  1.0085E+00  7.1873E-01  1.2766E+00  9.3359E-01
             1.7934E+00
 PARAMETER:  1.6271E-01 -4.5060E-01 -5.4848E-01  3.8027E-01 -5.2438E-01 -7.1955E-02  1.0849E-01 -2.3026E-01  3.4417E-01  3.1277E-02
             6.8411E-01
 GRADIENT:   3.7605E+02  5.2505E+01 -4.9867E+01  2.5432E+02  2.3423E+02 -6.6510E+01 -7.9849E+00  6.4598E+00 -3.4088E+01 -2.7301E+01
            -9.6920E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -931.359808578891        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      245
 NPARAMETR:  1.0884E+00  4.2869E-01  7.0251E-01  1.4541E+00  6.7112E-01  9.7811E-01  6.6717E-02  2.9518E-01  2.2031E+00  1.3991E+00
             1.9770E+00
 PARAMETER:  1.8469E-01 -7.4702E-01 -2.5310E-01  4.7437E-01 -2.9881E-01  7.7864E-02 -2.6073E+00 -1.1202E+00  8.8984E-01  4.3580E-01
             7.8159E-01
 GRADIENT:   1.9874E+02 -1.9561E+01 -6.7651E+01  5.0147E+01  1.5638E+02 -1.4712E+01  7.1497E-03  1.9749E-01  9.7974E+01 -7.2537E+00
            -7.6225E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1211.06782610823        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      426
 NPARAMETR:  1.0345E+00  2.9884E-01  4.9082E-01  1.4070E+00  4.1331E-01  8.9726E-01  1.1951E-01  4.9891E-02  1.1138E+00  8.0581E-01
             3.9149E+00
 PARAMETER:  1.3394E-01 -1.1079E+00 -6.1169E-01  4.4148E-01 -7.8355E-01 -8.4134E-03 -2.0243E+00 -2.8979E+00  2.0774E-01 -1.1590E-01
             1.4648E+00
 GRADIENT:  -5.5844E+00  7.0440E+00 -7.8951E+00  5.2990E+01  1.5045E+01 -1.7442E+01  6.1361E-03  4.7507E-02 -4.6486E+00  2.0439E+01
            -1.1327E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1219.35038894184        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      601
 NPARAMETR:  1.0307E+00  2.3757E-01  3.9677E-01  1.3559E+00  3.4735E-01  9.3868E-01  1.4220E+00  1.0000E-02  1.0232E+00  3.1629E-01
             4.1876E+00
 PARAMETER:  1.3028E-01 -1.3373E+00 -8.2440E-01  4.0449E-01 -9.5743E-01  3.6716E-02  4.5205E-01 -5.4443E+00  1.2289E-01 -1.0511E+00
             1.5321E+00
 GRADIENT:  -1.6504E+01  4.7468E+00  1.0886E+01  2.8665E+01 -1.4630E+01 -2.4945E+00  1.8071E-02  0.0000E+00 -1.2259E+01 -1.0133E+00
             6.0585E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1222.93144894226        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      778
 NPARAMETR:  1.0339E+00  1.5163E-01  2.8937E-01  1.2349E+00  2.7229E-01  9.2876E-01  1.3129E+00  1.0000E-02  1.2160E+00  3.9274E-01
             3.9682E+00
 PARAMETER:  1.3339E-01 -1.7863E+00 -1.1400E+00  3.1096E-01 -1.2009E+00  2.6100E-02  3.7227E-01 -1.0079E+01  2.9555E-01 -8.3460E-01
             1.4783E+00
 GRADIENT:   7.4449E+00 -2.3956E+00 -4.3377E+00 -1.9166E+01  1.1905E+01 -8.5274E+00  1.3922E-03  0.0000E+00 -1.7882E+00 -3.0587E+00
             5.5789E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1223.88943422418        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      953
 NPARAMETR:  1.0316E+00  1.5605E-01  2.9743E-01  1.2640E+00  2.7552E-01  9.5456E-01  7.6638E-01  1.0000E-02  1.2274E+00  5.0221E-01
             3.8540E+00
 PARAMETER:  1.3107E-01 -1.7576E+00 -1.1126E+00  3.3431E-01 -1.1891E+00  5.3493E-02 -1.6608E-01 -9.5495E+00  3.0489E-01 -5.8874E-01
             1.4491E+00
 GRADIENT:  -4.5207E-01 -9.7670E-01 -7.5760E-01  1.6916E+00  1.8371E+00 -3.7487E-01  2.2245E-02  0.0000E+00 -9.3614E-01 -1.5512E-01
            -4.4473E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1224.32404098185        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1131
 NPARAMETR:  1.0384E+00  2.2702E-01  2.8314E-01  1.2288E+00  2.7286E-01  9.5289E-01  1.1971E+00  1.0000E-02  1.2790E+00  5.1007E-01
             3.8358E+00
 PARAMETER:  1.3767E-01 -1.3827E+00 -1.1618E+00  3.0603E-01 -1.1988E+00  5.1746E-02  2.7993E-01 -9.5184E+00  3.4604E-01 -5.7321E-01
             1.4444E+00
 GRADIENT:   3.1942E+00 -6.6166E-01  1.7596E+00  4.7993E+00 -2.0536E+00 -2.6603E+00  5.5505E-01  0.0000E+00  1.5026E+00  1.2302E+00
             1.6934E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1225.28829251428        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1309
 NPARAMETR:  1.0388E+00  3.2269E-01  2.6284E-01  1.1782E+00  2.7210E-01  9.6207E-01  5.3009E-01  1.0000E-02  1.3438E+00  5.2945E-01
             3.8089E+00
 PARAMETER:  1.3811E-01 -1.0311E+00 -1.2362E+00  2.6395E-01 -1.2016E+00  6.1330E-02 -5.3470E-01 -9.7007E+00  3.9550E-01 -5.3592E-01
             1.4373E+00
 GRADIENT:  -3.7702E+00  2.8047E+00  1.4503E+01  8.6597E+00 -1.9571E+01 -2.5352E+00  2.2401E-01  0.0000E+00  1.3940E+00  2.9592E+00
             4.2635E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1229.60858971371        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1486
 NPARAMETR:  1.0286E+00  4.8767E-01  1.5640E-01  9.6994E-01  2.4634E-01  1.0490E+00  1.0000E-02  1.0000E-02  1.5966E+00  4.3616E-01
             3.6813E+00
 PARAMETER:  1.2822E-01 -6.1812E-01 -1.7553E+00  6.9482E-02 -1.3011E+00  1.4788E-01 -1.0071E+01 -1.3140E+01  5.6786E-01 -7.2974E-01
             1.4033E+00
 GRADIENT:   1.0433E+01  7.8791E+00  1.3479E+01 -4.1643E+00 -2.6206E+01  1.4167E+01  0.0000E+00  0.0000E+00 -2.4888E+00 -7.1108E+00
            -7.8927E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1231.05009560034        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1648
 NPARAMETR:  1.0164E+00  5.1305E-01  1.4167E-01  9.4378E-01  2.4402E-01  1.0067E+00  1.0000E-02  1.0000E-02  1.6834E+00  5.2983E-01
             3.6119E+00
 PARAMETER:  1.1628E-01 -5.6738E-01 -1.8543E+00  4.2138E-02 -1.3105E+00  1.0668E-01 -1.1697E+01 -1.3471E+01  6.2082E-01 -5.3521E-01
             1.3842E+00
 GRADIENT:   7.8372E-02  5.4094E-01  9.9206E-01 -1.3542E-01 -2.0724E+00  3.0174E-01  0.0000E+00  0.0000E+00  2.6568E-01 -1.1529E-02
            -2.3607E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1648
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.0015E-04 -1.2641E-04  1.7168E-04 -1.3625E-02  1.0926E-02
 SE:             2.8440E-02  2.0050E-04  1.6418E-04  2.5884E-02  1.8647E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7755E-01  5.2839E-01  2.9568E-01  5.9861E-01  5.5790E-01

 ETASHRINKSD(%)  4.7216E+00  9.9328E+01  9.9450E+01  1.3286E+01  3.7531E+01
 ETASHRINKVR(%)  9.2202E+00  9.9995E+01  9.9997E+01  2.4808E+01  6.0976E+01
 EBVSHRINKSD(%)  4.2269E+00  9.9336E+01  9.9517E+01  9.7538E+00  3.7681E+01
 EBVSHRINKVR(%)  8.2751E+00  9.9996E+01  9.9998E+01  1.8556E+01  6.1163E+01
 RELATIVEINF(%)  8.4681E+01  1.6481E-04  3.3552E-04  5.5268E+01  9.1477E-01
 EPSSHRINKSD(%)  2.7733E+01
 EPSSHRINKVR(%)  4.7775E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1231.0500956003407     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -495.89926903660250     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.06
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1231.050       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  5.13E-01  1.42E-01  9.44E-01  2.44E-01  1.01E+00  1.00E-02  1.00E-02  1.68E+00  5.30E-01  3.61E+00
 


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
+        9.78E+02
 
 TH 2
+       -5.84E+01  1.87E+03
 
 TH 3
+       -5.21E+02  3.66E+03  1.25E+04
 
 TH 4
+       -2.21E+01  1.03E+02 -5.76E+02  3.61E+02
 
 TH 5
+        4.89E+02 -6.57E+03 -1.50E+04 -9.97E+01  2.59E+04
 
 TH 6
+       -8.18E-01 -4.44E+00  8.13E+01 -1.17E+01  1.68E+01  1.64E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.12E+01 -2.67E+01  1.40E+02 -6.31E+00  7.73E+01  4.56E-01  0.00E+00  0.00E+00  4.19E+01
 
 TH10
+       -8.71E+00 -1.05E+02 -2.11E+02 -6.18E-01  4.78E+02  8.47E+00  0.00E+00  0.00E+00  4.46E-01  1.10E+02
 
 TH11
+       -1.75E+01 -1.26E+01 -2.63E+01 -3.91E+00  3.11E+01  2.07E+00  0.00E+00  0.00E+00  4.03E+00  2.07E+01  2.58E+01
 
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
 #CPUT: Total CPU Time in Seconds,       28.736
Stop Time:
Wed Sep 29 13:53:30 CDT 2021
