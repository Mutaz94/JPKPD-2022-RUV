Wed Sep 29 13:19:00 CDT 2021
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
$DATA ../../../../data/spa/A3/dat13.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m13.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   323.921743574233        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9282E+02  8.9471E+00  1.0465E+02 -8.7421E+01  6.3816E+01  2.9875E+01 -6.9830E+00 -6.0301E+01 -4.5055E+01 -3.8355E+01
            -3.7524E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1231.60930255962        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0913E+00  1.0656E+00  9.1270E-01  1.1695E+00  1.0993E+00  8.2398E-01  8.3389E-01  1.0273E+00  7.5656E-01  8.3253E-01
             5.4067E+00
 PARAMETER:  1.8738E-01  1.6356E-01  8.6529E-03  2.5660E-01  1.9469E-01 -9.3613E-02 -8.1651E-02  1.2695E-01 -1.7897E-01 -8.3289E-02
             1.7876E+00
 GRADIENT:   1.3441E+02 -1.6960E+01 -2.8299E+01  9.7552E+00  1.6783E+01 -2.1565E+01  1.2694E+01  6.7076E+00  2.0649E+01  1.0771E+01
             2.0719E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1255.30478837973        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0530E+00  1.2291E+00  5.1824E-01  9.9508E-01  7.1932E-01  9.5331E-01  2.1866E-01  8.3296E-01  7.1765E-01  4.0500E-01
             4.9855E+00
 PARAMETER:  1.5162E-01  3.0630E-01 -5.5732E-01  9.5068E-02 -2.2945E-01  5.2186E-02 -1.4202E+00 -8.2765E-02 -2.3177E-01 -8.0388E-01
             1.7065E+00
 GRADIENT:   3.5441E+01  9.4647E+01  3.7939E+01  3.6323E+01 -9.5531E+01  1.9084E+01  1.0212E-01  4.0880E+00  9.4177E+00  3.6215E+00
             1.4674E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1282.79573386107        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0123E+00  1.2355E+00  9.9183E-01  9.9836E-01  1.1790E+00  8.8257E-01  3.3288E-01  6.6682E-01  7.8619E-01  1.7681E-01
             4.0440E+00
 PARAMETER:  1.1226E-01  3.1149E-01  9.1796E-02  9.8363E-02  2.6464E-01 -2.4923E-02 -9.9997E-01 -3.0523E-01 -1.4056E-01 -1.6327E+00
             1.4972E+00
 GRADIENT:   1.2068E+01  2.0768E+01 -3.7962E+00  2.4992E+01  2.3779E+00 -4.2608E+00  5.4224E-01  1.4646E+00  2.2407E+00  2.8280E-01
             7.5328E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1283.91165385579        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0007E+00  1.0085E+00  1.2088E+00  1.1208E+00  1.1656E+00  8.8485E-01  2.1853E-01  1.9972E-01  6.8986E-01  1.3120E-01
             4.0290E+00
 PARAMETER:  1.0067E-01  1.0844E-01  2.8959E-01  2.1407E-01  2.5320E-01 -2.2335E-02 -1.4208E+00 -1.5108E+00 -2.7127E-01 -1.9310E+00
             1.4935E+00
 GRADIENT:  -2.8406E+00 -2.1698E+00 -9.7204E-01  1.4499E+00  2.6021E+00 -1.1568E+00  1.4266E-01  1.5992E-01  5.1946E-01  1.2927E-01
             4.2324E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1284.26676105616        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      480
 NPARAMETR:  1.0098E+00  9.0884E-01  1.3094E+00  1.1939E+00  1.1596E+00  8.8891E-01  1.3846E-01  7.7612E-02  6.4915E-01  1.0954E-01
             4.0765E+00
 PARAMETER:  1.0971E-01  4.4144E-03  3.6954E-01  2.7722E-01  2.4808E-01 -1.7761E-02 -1.8772E+00 -2.4560E+00 -3.3209E-01 -2.1114E+00
             1.5052E+00
 GRADIENT:  -4.2871E+00  1.7444E+00  1.2801E+00  1.3571E+00 -2.0526E+00 -9.1464E-01  8.2693E-03  2.2130E-02  2.0782E-01  7.6694E-02
             1.7835E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1284.39293817580        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      656
 NPARAMETR:  1.0103E+00  7.1561E-01  1.2311E+00  1.3109E+00  1.0460E+00  8.9382E-01  2.7789E-02  1.0000E-02  6.0116E-01  6.4178E-02
             4.0681E+00
 PARAMETER:  1.1022E-01 -2.3462E-01  3.0788E-01  3.7075E-01  1.4497E-01 -1.2256E-02 -3.4831E+00 -5.2507E+00 -4.0889E-01 -2.6461E+00
             1.5032E+00
 GRADIENT:  -1.4725E-01  1.9190E-01 -7.4428E-01  6.0864E-01  9.8411E-01 -1.1416E-02  8.5388E-05  0.0000E+00 -6.8940E-02  2.8510E-02
            -4.4576E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1284.43418535762        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  1.0091E+00  5.8067E-01  1.2598E+00  1.3930E+00  1.0037E+00  8.9392E-01  1.0000E-02  1.0000E-02  5.6524E-01  3.4794E-02
             4.0697E+00
 PARAMETER:  1.0902E-01 -4.4358E-01  3.3092E-01  4.3143E-01  1.0373E-01 -1.2141E-02 -5.0534E+00 -7.9504E+00 -4.7050E-01 -3.2583E+00
             1.5036E+00
 GRADIENT:   1.0630E+00  1.5725E-01  2.3096E-01 -8.5659E-02 -3.0000E-01  1.5308E-01  0.0000E+00  0.0000E+00  4.7228E-02  8.5555E-03
             5.9570E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1284.44415914559        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  1.0069E+00  4.8725E-01  1.2654E+00  1.4497E+00  9.7374E-01  8.9333E-01  1.0000E-02  1.0000E-02  5.4307E-01  1.9334E-02
             4.0625E+00
 PARAMETER:  1.0688E-01 -6.1897E-01  3.3537E-01  4.7136E-01  7.3386E-02 -1.2803E-02 -6.5668E+00 -1.0471E+01 -5.1051E-01 -3.8459E+00
             1.5018E+00
 GRADIENT:  -1.2073E+00  4.7424E-01  1.7139E-02  1.7489E+00 -1.5259E-01 -1.9908E-01  0.0000E+00  0.0000E+00 -9.9825E-02  2.6735E-03
            -5.8008E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1284.44631357160        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1198
 NPARAMETR:  1.0069E+00  4.6161E-01  1.2599E+00  1.4636E+00  9.6224E-01  8.9377E-01  1.0000E-02  1.0000E-02  5.3860E-01  1.0000E-02
             4.0631E+00
 PARAMETER:  1.0689E-01 -6.7304E-01  3.3102E-01  4.8093E-01  6.1514E-02 -1.2312E-02 -7.0391E+00 -1.1258E+01 -5.1878E-01 -4.5635E+00
             1.5020E+00
 GRADIENT:  -4.8598E-02 -7.0853E-02 -8.6654E-02 -4.6071E-01  2.0733E-01 -1.0299E-02  0.0000E+00  0.0000E+00  2.3332E-02  0.0000E+00
            -3.2360E-02

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1284.44633960807        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1290
 NPARAMETR:  1.0069E+00  4.6160E-01  1.2603E+00  1.4637E+00  9.6196E-01  8.9381E-01  1.0000E-02  1.0000E-02  5.3843E-01  1.0000E-02
             4.0631E+00
 PARAMETER:  1.0689E-01 -6.7305E-01  3.3135E-01  4.8095E-01  6.1217E-02 -1.2260E-02 -7.0391E+00 -1.1258E+01 -5.1910E-01 -4.5254E+00
             1.5020E+00
 GRADIENT:  -7.3678E-02  2.3214E-02  6.8293E-02 -2.2291E-01 -6.3105E-02  1.1190E-02  0.0000E+00  0.0000E+00 -3.7355E-03  2.2274E-04
            -7.9025E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1290
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3230E-03 -1.2074E-04  7.0607E-05 -1.2989E-02  4.6825E-06
 SE:             2.8383E-02  8.0246E-05  7.7558E-05  1.9055E-02  1.3355E-04
 N:                     100         100         100         100         100

 P VAL.:         9.3477E-01  1.3242E-01  3.6263E-01  4.9544E-01  9.7203E-01

 ETASHRINKSD(%)  4.9133E+00  9.9731E+01  9.9740E+01  3.6165E+01  9.9553E+01
 ETASHRINKVR(%)  9.5852E+00  9.9999E+01  9.9999E+01  5.9251E+01  9.9998E+01
 EBVSHRINKSD(%)  4.8399E+00  9.9738E+01  9.9693E+01  3.6301E+01  9.9515E+01
 EBVSHRINKVR(%)  9.4456E+00  9.9999E+01  9.9999E+01  5.9424E+01  9.9998E+01
 RELATIVEINF(%)  7.9917E+01  6.4098E-06  2.6031E-05  5.4422E-01  4.0241E-05
 EPSSHRINKSD(%)  1.6020E+01
 EPSSHRINKVR(%)  2.9473E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1284.4463396080680     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -549.29551304432982     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.11
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1284.446       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  4.62E-01  1.26E+00  1.46E+00  9.62E-01  8.94E-01  1.00E-02  1.00E-02  5.38E-01  1.00E-02  4.06E+00
 


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
+        1.25E+03
 
 TH 2
+       -1.34E+02  3.69E+02
 
 TH 3
+        3.88E+00  5.97E+01  5.91E+01
 
 TH 4
+       -1.82E+02  4.85E+02  2.90E+01  7.08E+02
 
 TH 5
+        2.91E+01 -1.99E+02 -1.23E+02 -1.63E+02  2.81E+02
 
 TH 6
+        4.09E+00 -2.17E+01  6.49E+00 -3.83E+01 -1.94E+00  2.02E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.99E+00 -4.98E+01 -2.33E+00 -3.41E+01  2.58E+01 -1.19E+00  0.00E+00  0.00E+00  9.50E+01
 
 TH10
+        5.50E-02 -2.32E-02 -1.55E-02 -7.63E-02  6.94E-02 -2.86E-02  0.00E+00  0.00E+00  1.67E-02  3.47E+01
 
 TH11
+       -1.63E+01 -1.73E+01 -1.11E+00 -1.67E+01  6.11E+00  4.99E+00  0.00E+00  0.00E+00  2.23E+01  2.09E-02  2.91E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.116
Stop Time:
Wed Sep 29 13:19:23 CDT 2021
